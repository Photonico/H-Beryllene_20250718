#!/usr/bin/env python3
"""Preflight or submit the failed-SCF recovery and independent topology analyses."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / '9.0_topology'
TOOLS = ROOT / 'tools'
PREFIX = ROOT / 'packages/soc_topology'
RECOVERY = TOOLS / 'recovery_20260923'
FAILED = {'beta': '41849.headnode', 'beta_1h': '41852.headnode', 'beta_2h': '41853.headnode'}
FINISHED = ('alpha', 'st', 'alpha_2h')


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def save(path, obj):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def history(jid):
    text = subprocess.check_output(['qstat', '-xf', jid], text=True)
    result = {'id': jid}
    for line in text.splitlines():
        m = re.match(r'\s*(job_state|Exit_status|resources_used.mem|resources_used.walltime|Resource_List.mem) = (.*)', line)
        if m:
            result[m[1]] = m[2]
    return result


def pbs(name, body):
    return ('#!/bin/bash\n#PBS -N ' + name + '\n#PBS -q cmt\n'
            '#PBS -l select=1:ncpus=1:mpiprocs=1:mem=12gb\n'
            '#PBS -l walltime=4:00:00\n#PBS -j oe\nset -euo pipefail\n'
            + (TOOLS / 'modules.sh').read_text() + '\n' + body + '\n')


def preflight():
    require(not RECOVERY.exists(), 'Recovery directory already exists; inspect its journal before any retry.')
    active = json.loads((TOOLS / 'ACTIVE_SUBMISSION.json').read_text())
    jobs = {j['label']: j for j in active['jobs']}
    require(all(jobs[s]['id'] == j for s, j in FAILED.items()), 'Previous campaign IDs changed.')
    statuses = {label: history(j['id']) for label, j in jobs.items()}
    require(all(x.get('job_state') == 'F' for x in statuses.values()), 'Previous campaign is still active.')
    for sid in FAILED:
        require(statuses[sid].get('Exit_status') == '255', sid + ': unexpected previous exit status')
    for sid in (*FINISHED, 'tools'):
        require(statuses[sid].get('Exit_status') == '0', sid + ': prerequisite did not finish normally')
    queue = subprocess.check_output(['qstat', '-u', 'lniu6305'], text=True)
    require('topo_' not in queue, 'Another topology job is in the queue.')
    entries = json.loads((ROOT / 'manifest.json').read_text())['structures']
    entries = {e['id']: e for e in entries}
    inputs = {}
    for sid, e in entries.items():
        d = ROOT / e['directory']
        src = REPO / e['reference']
        for stage in ('scf', 'trim', 'bands'):
            require(digest(d/stage/'POSCAR') == e['source_files']['CONTCAR']['sha256'], sid + ': geometry changed')
            require(digest(d/stage/'POTCAR') == e['source_files']['POTCAR']['sha256'], sid + ': potential changed')
        lines = (d/'scf/KPOINTS').read_text().splitlines()
        require(int(lines[1]) == 0 and lines[2].strip().lower().startswith('g')
                and [int(x) for x in lines[3].split()] == [105, 105, 1]
                and [float(x) for x in lines[4].split()] == [0, 0, 0], sid + ': unexpected SCF mesh')
        lines = (d/'trim/KPOINTS').read_text().splitlines()
        rows = [list(map(float, line.split()[:4])) for line in lines[3:] if line.strip()]
        require(int(lines[1]) == 4 and lines[2].lower().startswith('r') and len(rows) == 4
                and {tuple(x[:3]) for x in rows} == {(0,0,0),(.5,0,0),(0,.5,0),(.5,.5,0)}
                and all(x[3] == 1 for x in rows), sid + ': incorrect TRIM points')
        ref = 'c3-Beryllene_trilayer' if sid == 'st' else src.name
        require(digest(d/'bands/KPOINTS') == digest(REPO/'4.0_bandstructure'/ref/'KPOINTS'), sid + ': reference path changed')
        inputs[sid] = {stage: digest(d/stage/'KPOINTS') for stage in ('scf','trim','bands')}
        if sid in FINISHED:
            require((d/'PIPELINE_COMPLETE').is_file(), sid + ': completion marker absent')
            for stage in ('scf','trim','bands','gap_refine_1','gap_refine_2'):
                require(json.loads((d/stage/'topology_validation.json').read_text()).get('validation_passed') is True,
                        sid + '/' + stage + ': invalid completed stage')
        if sid in FAILED:
            require(not (d/'PIPELINE_COMPLETE').exists(), sid + ': already completed')
            require(not (d/'scf/topology_validation.json').exists(), sid + ': existing validated SCF')
            for stage in ('trim','bands','gap_refine_1','gap_refine_2'):
                require(not (d/stage/'OUTCAR').exists(), sid + ': later output needs manual resume review')
            require(digest(src/'CHGCAR') == e['source_files']['CHGCAR']['sha256'], sid + ': source density changed')
            require(digest(d/'scf/CHGCAR') == digest(src/'CHGCAR'), sid + ': initial density changed')
            incar = (d/'scf/INCAR').read_text()
            require(len(re.findall(r'^KPAR\s*=\s*7\s*$', incar, re.M)) == 1, sid + ': unexpected KPAR')
            require(re.search(r'^NCORE\s*=\s*6\s*$', incar, re.M), sid + ': unexpected NCORE')
            require('KILLED BY SIGNAL: 9' in (d/'scf/vasp.log').read_text(), sid + ': unexpected failure')
    require((PREFIX/'READY').is_file() and (PREFIX/'READY').stat().st_size > 0, 'Tools not ready')
    require(os.access(PREFIX/'vasp/bin/vasp_ncl', os.X_OK), 'Private SOC VASP missing')
    require(not (ROOT/entries['beta_1h']['directory']/'wcc_direct').exists(), 'WCC outputs already exist')
    return active, statuses, entries, inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--submit', action='store_true', help='Apply reviewed parallelism change and submit once')
    args = parser.parse_args()
    active, statuses, entries, inputs = preflight()
    print('PREFLIGHT PASSED: original jobs finished; six geometries, potentials and k-point inputs verified.', flush=True)
    print('PLAN: retry 3 beta SCFs with KPAR=1/NCORE=6, 42 ranks/200gb; 5 independent parity jobs; 1 dependent WCC job.', flush=True)
    if not args.submit:
        return
    RECOVERY.mkdir()
    save(RECOVERY/'audit.json', {'created_utc': stamp(), 'previous_submission': active, 'previous_jobs': statuses,
                              'kpoints_sha256': inputs, 'parallelism_change': 'Only failed beta SCF KPAR 7 -> 1; NCORE=6, 42 ranks, 200gb unchanged.',
                              'physics_changed': False})
    scripts = {}
    for sid, oldjob in FAILED.items():
        d = ROOT/entries[sid]['directory']
        (RECOVERY/(sid+'_failed_'+oldjob.split('.')[0]+'.log')).write_text((d/'scf/vasp.log').read_text())
        path = d/'scf/INCAR'
        before = path.read_text()
        after = re.sub(r'^KPAR\s*=\s*7\s*$', 'KPAR = 1', before, flags=re.M)
        path.write_text(after)
        save(RECOVERY/(sid+'_input_change.json'), {'before_incar': before, 'after_incar': after,
                                                'kpoints_sha256': digest(d/'scf/KPOINTS')})
        scripts[sid] = RECOVERY/(sid+'.pbs')
        scripts[sid].write_text((d/'run.pbs').read_text().replace('#PBS -N topo_'+sid, '#PBS -N r_'+sid))
    python = PREFIX/'venv/bin/python'
    for sid, e in entries.items():
        if not e['inversion_expected']:
            continue
        label = 'parity_'+sid
        scripts[label] = RECOVERY/(label+'.pbs')
        body = "cd '" + str(TOOLS) + "'\n'" + str(python) + "' -B parity_analysis.py '" + str(ROOT/e['directory']/'trim') + "' --expected-nelect " + str(e['expected_nelect'])
        scripts[label].write_text(pbs('p_'+sid, body))
    scripts['wcc'] = RECOVERY/'wcc.pbs'
    scripts['wcc'].write_text((TOOLS/'wcc.pbs').read_text().replace('#PBS -N topo_wcc', '#PBS -N r_wcc').replace('--gap-stage gap_refine_1', '--gap-stage bands --gap-stage gap_refine_1'))
    for p in scripts.values():
        subprocess.run(['bash','-n',str(p)], check=True)
    state = {'submitted_utc': stamp(), 'recovery_directory': str(RECOVERY.relative_to(REPO)),
             'previous_submission_record': str((RECOVERY/'audit.json').relative_to(REPO)), 'jobs': []}
    journal = (RECOVERY/'submission_journal.jsonl').open('x', buffering=1)
    def event(row):
        journal.write(json.dumps(dict(utc=stamp(), **row))+'\n'); os.fsync(journal.fileno())
    def submit(label, deps=()):
        command = ['qsub','-o',str(RECOVERY/(label+'.pbs.log'))]
        if deps:
            command += ['-W','depend=afterok:'+':'.join(deps)]
        command += [str(scripts[label])]
        event({'event':'attempt','label':label,'command':command})
        result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        event({'event':'response','label':label,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        jid = result.stdout.strip()
        require(result.returncode == 0 and re.fullmatch(r'\d+(?:\.[\w.-]+)?',jid), 'Submission uncertain/failed; inspect journal and queue before retry: '+result.stderr)
        state['jobs'].append({'label':label,'id':jid,'pbs':str(scripts[label].relative_to(REPO)),
                              'dependencies':list(deps),'submitted_utc':stamp()})
        save(RECOVERY/'submission.json', state)
        save(TOOLS/'ACTIVE_SUBMISSION.json', state)
        print('SUBMITTED',label,jid,flush=True)
        return jid
    for sid in FINISHED:
        submit('parity_'+sid)
    ids = {sid: submit(sid) for sid in FAILED}
    for sid in ('beta','beta_2h'):
        submit('parity_'+sid,[ids[sid]])
    submit('wcc',[ids['beta_1h']])
    state['submission_complete'] = True
    save(RECOVERY/'submission.json', state)
    save(TOOLS/'ACTIVE_SUBMISSION.json', state)
    journal.close()


if __name__ == '__main__':
    main()
