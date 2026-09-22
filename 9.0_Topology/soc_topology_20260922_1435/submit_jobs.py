#!/usr/bin/env python3
"""Submit this new campaign once, recording every accepted PBS identifier."""
import ast, json, subprocess, re, os
from pathlib import Path
from datetime import datetime, timezone
RUN=Path(__file__).resolve().parent

def main():
    manifest=json.loads((RUN/'manifest.json').read_text())
    tools=Path(manifest['tools_prefix'])
    if not (tools/'DOWNLOADS_READY').is_file() or not (tools/'DOWNLOADS_READY').stat().st_size:
        raise RuntimeError('Complete isolated tool downloads first')
    structures=manifest['structures']
    names=[s['id'] for s in structures]
    if set(names)!={'alpha','beta','st','alpha_2h','beta_1h','beta_2h'} or len(names)!=6:
        raise RuntimeError('Expected exactly the six unique campaign structure IDs')
    if (RUN/'SUBMISSION_COMPLETE').exists():
        raise RuntimeError('Campaign was already submitted')
    scripts=[RUN/'build_tools.pbs',RUN/'parity.pbs',RUN/'wcc.pbs']
    scripts += [RUN/name/'run.pbs' for name in names]
    scripts += [RUN/'build_tools.sh']
    python_scripts=[RUN/name for name in ('check_topology_outputs.py','parity_analysis.py',
                                         'wcc_driver.py','wcc_line_check.py')]
    # Validate the complete campaign before the first scheduler mutation.
    for path in scripts+python_scripts:
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError('Missing campaign script: '+str(path))
    for path in scripts:
        subprocess.run(['bash','-n',str(path)],check=True)
    for path in python_scripts:
        ast.parse(path.read_text(),filename=str(path))
    record=RUN/'submitted_jobs.jsonl'
    # Exclusive creation means an interrupted submission cannot silently duplicate jobs.
    with record.open('x',buffering=1) as log, (RUN/'submission_attempts.jsonl').open('x',buffering=1) as attempts:
        def journal(item):
            attempts.write(json.dumps(item)+'\n'); os.fsync(attempts.fileno())
        def submit(label,path,deps=()):
            command=['qsub','-o',str(RUN/(label+'.pbs.log'))]
            if deps: command+=['-W','depend=afterok:'+':'.join(deps)]
            command+=[str(path)]
            journal(dict(event='attempt',label=label,command=command,
                         utc=datetime.now(timezone.utc).isoformat()))
            response=subprocess.run(command,cwd=RUN,text=True,capture_output=True)
            journal(dict(event='response',label=label,returncode=response.returncode,
                         stdout=response.stdout,stderr=response.stderr,
                         utc=datetime.now(timezone.utc).isoformat()))
            if response.returncode:
                raise RuntimeError('PBS submission failed; inspect submitted_jobs.jsonl, '
                                   'submission_attempts.jsonl and the queue before retrying: '
                                   +response.stderr.strip())
            output=response.stdout.strip()
            if not re.fullmatch(r'\d+(?:\[\])?(?:\.[\w.-]+)?',output):
                raise RuntimeError('Unrecognized PBS response; inspect queue before any retry: '+output)
            item=dict(label=label,id=output,pbs=str(path.relative_to(RUN)),dependencies=list(deps),submitted_utc=datetime.now(timezone.utc).isoformat())
            log.write(json.dumps(item)+'\n'); os.fsync(log.fileno())
            print(label,output,flush=True)
            return output
        build=submit('tools',RUN/'build_tools.pbs')
        ids={s['id']:submit(s['id'],RUN/s['id']/'run.pbs') for s in structures}
        submit('parity',RUN/'parity.pbs',[build]+[ids[s['id']] for s in structures if s['inversion_expected']])
        submit('wcc',RUN/'wcc.pbs',[build,ids['beta_1h']])
    with (RUN/'SUBMISSION_COMPLETE').open('x') as complete:
        complete.write(datetime.now(timezone.utc).isoformat()+'\n')

if __name__=='__main__': main()
