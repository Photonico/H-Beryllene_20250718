import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path('/cmt2/lniu6305/H-Beryllene_20250718')
ROOT = REPO / '9.0_topology'
SCRIPTS = ROOT / 'tools'
PREFIX = Path('/cmt2/lniu6305/Packages/soc_topology_20260922_1435')
PYTHON = PREFIX / 'venv/bin/python'
VASP = Path('/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_ncl')
DIRECTORIES = dict(alpha='a-Beryllene', beta='b-Beryllene', st='c-Beryllene_trilayer',
                   alpha_2h='a-Beryllene_hh', beta_1h='b-Beryllene_b', beta_2h='b-Beryllene_bb')
STAMP = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
sys.path.insert(0, str(SCRIPTS))
from check_topology_outputs import read_poscar, inversion_check

def pbs(name, body, ranks=42, mem='200gb', hours=48):
    modules = (SCRIPTS / 'modules.sh').read_text()
    return f'''#!/bin/bash
#PBS -N {name}
#PBS -q cmt
#PBS -l select=1:ncpus={ranks}:mpiprocs={ranks}:mem={mem}
#PBS -l walltime={hours}:00:00
#PBS -j oe
set -euo pipefail
{modules}
{body}
'''

manifest = json.loads((SCRIPTS / 'manifest.json').read_text())
manifest.update(submission_revision_utc=STAMP, campaign_root=str(ROOT))
for entry in manifest['structures']:
    sid = entry['id']
    directory = ROOT / DIRECTORIES[sid]
    entry['directory'] = DIRECTORIES[sid]
    geometry = inversion_check(read_poscar(directory/'scf/POSCAR'), 1e-5, False)
    if geometry['has_inversion'] != entry['inversion_expected']:
        raise RuntimeError(f'{sid}: inversion differs from planned method')
    entry['inversion_checked'] = geometry
    for stage in ('scf','trim','bands'):
        if (directory/stage/'OUTCAR').exists():
            raise RuntimeError(f'{sid}/{stage}: existing output requires a resume decision')
    (directory/'manifest.json').write_text(json.dumps(entry, indent=2)+'\n')
    body = f'''cd '{directory}'
for stage in scf trim bands; do
  cd '{directory}'/"$stage"
  if [ "$stage" != scf ]; then cp ../scf/CHGCAR CHGCAR; fi
  mpirun -np 42 '{VASP}' > vasp.log 2>&1
  extra=()
  if [ "$stage" = trim ]; then extra+=(--require-wavecar); fi
  '{PYTHON}' -B '{SCRIPTS}/check_topology_outputs.py' . --manifest '{directory}/manifest.json' --structure-id '{sid}' "${{extra[@]}}"
done
cd '{directory}'
'{PYTHON}' -B '{SCRIPTS}/refine_gap.py' '{directory}' --vasp '{VASP}' --ranks 42
printf 'Completed {sid}\\n' > PIPELINE_COMPLETE
'''
    (directory/'run.pbs').write_text(pbs('topo_'+sid, body))
    print('PREPARED', sid, str(directory), 'P='+str(geometry['has_inversion']))

(ROOT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
(SCRIPTS/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
(SCRIPTS/'build_tools.pbs').write_text(pbs('topo_tools',
    f"cd '{SCRIPTS}'\nbash build_tools.sh '{PREFIX}' > '{PREFIX}/build.log' 2>&1",12,'48gb',12))
(SCRIPTS/'parity.pbs').write_text(pbs('topo_parity',
    f"cd '{SCRIPTS}'\n'{PYTHON}' -B parity_analysis.py '{ROOT}' > '{ROOT}/parity.log' 2>&1",1,'12gb',4))
(SCRIPTS/'wcc.pbs').write_text(pbs('topo_wcc',
    f"cd '{SCRIPTS}'\n'{PYTHON}' -B wcc_driver.py '{ROOT}/b-Beryllene_b' '{PREFIX}' --gap-stage gap_refine_1 --gap-stage gap_refine_2 > '{ROOT}/b-Beryllene_b/wcc_driver.log' 2>&1"))

record_path = SCRIPTS / f'submitted_jobs_{STAMP}.jsonl'
if (SCRIPTS/'ACTIVE_SUBMISSION.json').exists():
    raise RuntimeError('Current submission record exists; do not duplicate jobs')
jobs=[]
def submit(label, script, dependencies=()):
    command=['qsub']
    if dependencies:
        command += ['-W', 'depend=afterok:'+':'.join(dependencies)]
    command += [str(script)]
    job_id=subprocess.check_output(command,cwd=str(ROOT),text=True).strip()
    if not re.fullmatch(r'\d+(?:\.[A-Za-z0-9_.-]+)?',job_id):
        raise RuntimeError('Unexpected qsub response: '+job_id)
    row=dict(label=label,id=job_id,pbs=str(script.relative_to(REPO)),dependencies=list(dependencies),submitted_utc=datetime.now(timezone.utc).isoformat())
    jobs.append(row)
    with record_path.open('a') as stream: stream.write(json.dumps(row)+'\n')
    (SCRIPTS/'ACTIVE_SUBMISSION.json').write_text(json.dumps(dict(submitted_utc=STAMP,jobs=jobs),indent=2)+'\n')
    print('SUBMITTED',label,job_id,flush=True)
    return job_id

build=submit('tools',SCRIPTS/'build_tools.pbs')
ids={sid:submit(sid,ROOT/folder/'run.pbs') for sid,folder in DIRECTORIES.items()}
submit('parity',SCRIPTS/'parity.pbs',[ids[x] for x in DIRECTORIES if x!='beta_1h'])
submit('wcc',SCRIPTS/'wcc.pbs',[build,ids['beta_1h']])

rows='\n'.join(f"| {j['label']} | `{j['id']}` | {', '.join(j['dependencies']) or 'None'} | `{j['pbs']}` |" for j in jobs)
structures='\n'.join(f"| {v['id']} | `{v['directory']}` | {v['expected_nelect']} | {'Yes' if v['inversion_checked']['has_inversion'] else 'No'} | {'Fu–Kane parity' if v['inversion_expected'] else 'Direct SOC-DFT Wilson loop / WCC, lowest 4 and 6 bands'} |" for v in manifest['structures'])
summary=f'''# Beryllene SOC topology: submitted calculation campaign

Submitted UTC: **{STAMP}**. Repository: `/cmt2/lniu6305/H-Beryllene_20250718`.

All paths use the current `9.0_topology` directory. `20260922_draft` is preserved as collaborator-provided completed legacy calculations; its manuscript topology claims are not treated as validated for the current reference geometries.

## Structures and analysis methods

| Structure | Directory under 9.0_topology | Neutral spinor filling | Geometric inversion | Method |
|---|---|---:|---|---|
{structures}

Geometry is taken from the existing frozen `6.0_dielectric_selection` CONTCARs. The existing `4.0_bandstructure` k paths and PAW potentials are reused. Inversion was independently checked on each complete Be/H input structure at 1e-5 angstrom tolerance, allowing shifted inversion centers. Electronic time reversal still requires checking the converged state.

## Submitted PBS jobs

| Task | PBS job ID | afterok dependencies | Script |
|---|---|---|---|
{rows}

Submission means accepted by PBS, not calculation completion. The authoritative new records are `9.0_topology/tools/ACTIVE_SUBMISSION.json` and `{record_path.relative_to(REPO)}`. Legacy submission records are not used.

## Calculation sequence and parameters

Each material pipeline runs SOC SCF on a Gamma-centered 105 x 105 x 1 grid; explicit Gamma/X/Y/M TRIM; the original band path; then two local direct-gap refinement stages around sampled minima. Parameters: PBE, ENCUT 600 eV, Accurate precision, EDIFF 1e-8 eV, SIGMA 0.01 eV, 28 spinor bands, zero initial moments, existing slab/dipole setup, fixed structures. Site VASP: `{VASP}`. A scalar-relativistic CHGCAR supplies the initial density; its WAVECAR is not reused for SOC.

For the lowest N=2M spinor bands, every stage reports min_k(E_(N+1)-E_N), min_k E_(N+1)-max_k E_N, and the minimizing k point. A discrete grid and local refinement provide numerical evidence, not a rigorous all-k proof. Very small gaps or unstable minima remain unresolved pending convergence study. Negative indirect gap must not be described as an insulating Fermi-level gap.

Fu–Kane analysis uses IrRep 2.1.3 with SOC spinors and four distinct TRIM, retaining the inversion operation, degenerate blocks/Kramers-pair parities, delta_i and conditional Z2. The two square-lattice X-type TRIM are retained separately. The parity job does not depend on Wannier90 compilation.

1H-beta has five neutral electrons per primitive cell: lowest-4/lowest-6 manifolds are conditional isolated subspaces, not the ordinary insulating occupied manifold. WCC runs only when all sampled stages keep direct separation above 0.1 meV. It uses direct SOC-DFT overlaps through the VASP–Wannier90 interface and Z2Pack; it does not claim an already-fitted Wannier Hamiltonian. Loop-position, gap and movement convergence are retained.

## Software and next review

Existing Python environment: IrRep 2.1.3, Z2Pack 2.2.1, spglib 2.5.0, NumPy 1.26.4 and SciPy 1.13.1. The tools job builds downloaded Wannier90 3.1.0 and a private VASP 6.5.0 SOC executable linked to it in `{PREFIX}`. WCC waits for both that build and the 1H-beta pipeline. No edge-state jobs are included at this stage.

After the jobs finish, run the following read-only notebook cells to collect completion, gap, parity and WCC reports. Inspect failed prerequisites, magnetism, local-gap convergence and matching band indices before making final topology claims.
'''

def cell(kind, source):
    value={'cell_type':kind,'metadata':{},'source':source.splitlines(True)}
    if kind=='code': value.update(execution_count=None,outputs=[])
    return value

status_code='''from pathlib import Path
import json, subprocess
repo = Path.cwd()
if not (repo / "9.0_topology").is_dir():
    repo = Path("/cmt2/lniu6305/H-Beryllene_20250718")
topology = repo / "9.0_topology"
submission = json.loads((topology / "tools/ACTIVE_SUBMISSION.json").read_text())
print(subprocess.run(["qstat", "-x", *[j["id"] for j in submission["jobs"]]], capture_output=True, text=True).stdout)
for entry in json.loads((topology / "manifest.json").read_text())["structures"]:
    folder = topology / entry["directory"]
    print(entry["id"], "pipeline complete:", (folder / "PIPELINE_COMPLETE").exists())
'''
analysis_code='''for entry in json.loads((topology / "manifest.json").read_text())["structures"]:
    folder = topology / entry["directory"]
    print("\\n", entry["id"])
    for stage in ("scf", "trim", "bands", "gap_refine_1", "gap_refine_2"):
        report = folder / stage / "topology_validation.json"
        if report.exists():
            data = json.loads(report.read_text())
            print(stage, "valid:", data.get("validation_passed"))
            for gap in data.get("eigenvalues", {}).get("manifolds", []):
                print("  N=", gap["N"], "direct eV=", gap["direct_gap_ev"], "indirect eV=", gap["indirect_gap_ev"], "k=", gap["direct_gap_k"])
    for pattern in ("trim/parity-analysis-*/parity_summary.json", "trim/parity-analysis-*/PARITY_NOT_VALIDATED.json", "wcc_direct/status_*.json"):
        for report in sorted(folder.glob(pattern)):
            print(report.relative_to(repo))
            print(json.dumps(json.loads(report.read_text()), indent=2))
'''
nb_path=REPO/'9.0_topology.ipynb'
nb=json.loads(nb_path.read_text())
old_cells=nb.get('cells',[])
nb.update(nbformat=4,nbformat_minor=5)
nb['cells']=[cell('markdown',summary),cell('code',status_code),cell('code',analysis_code)]
if any(''.join(c.get('source',[])).strip() for c in old_cells):
    nb['cells'] += [cell('markdown','## Previous notebook notes (preserved)\n')] + old_cells
for index,c in enumerate(nb['cells']): c.setdefault('id',f'topology-{index}-{STAMP.lower()}')
nb.setdefault('metadata',{}).setdefault('kernelspec',dict(display_name='Python 3',language='python',name='python3'))
nb_path.write_text(json.dumps(nb,indent=1,ensure_ascii=False)+'\n')
print('NOTEBOOK_UPDATED',nb_path)
print(subprocess.check_output(['qstat','-x',*[j['id'] for j in jobs]],text=True))
