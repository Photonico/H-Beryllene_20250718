#!/usr/bin/env python3
"""Create a new immutable-input topology campaign without touching old runs."""
import hashlib, json, os, re, shutil, subprocess
from pathlib import Path

REPO = Path('/cmt2/lniu6305/H-Beryllene_20250718')
RUN = Path(__file__).resolve().parent
TOOLS = Path('/cmt2/lniu6305/Packages') / RUN.name
VASP = Path('/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_ncl')
SPECS = [('alpha','a-Beryllene',2,1,True),('beta','b-Beryllene',4,2,True),('st','c3-Beryllene',6,3,True),('alpha_2h','a-Beryllene_hh',4,3,True),('beta_1h','b-Beryllene_b',5,3,False),('beta_2h','b-Beryllene_bb',6,4,True)]
MODULES = '''set +u
source /etc/profile.d/modules.sh
module load pbspro oneapi-2024.2/compiler/2024.2.1 oneapi-2024.2/mkl/latest oneapi-2024.2/mpi/latest hdf/5/1.14.1-2_intel2021
set -u
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
'''

def write(path,text):
    with path.open('x') as f: f.write(text)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()

def pbs(name,body,ranks=42,mem='200gb',hours=48):
    return f'''#!/bin/bash
#PBS -N {name}
#PBS -q cmt
#PBS -l select=1:ncpus={ranks}:mpiprocs={ranks}:mem={mem}
#PBS -l walltime={hours}:00:00
#PBS -j oe
set -euo pipefail
{MODULES}
{body}
'''

def main():
    assert REPO.resolve() in RUN.parents
    assert not (RUN/'manifest.json').exists(), 'Campaign already prepared'
    assert not TOOLS.exists(), 'Refusing an existing tools prefix'
    manifests=[]
    # Preflight every reference before creating a calculation directory.
    for sid,ref,nelect,natom,inversion in SPECS:
        src=REPO/'6.0_dielectric_selection'/ref
        for name in ('CONTCAR','INCAR','POTCAR','KPOINTS','CHGCAR'):
            assert (src/name).is_file() and (src/name).stat().st_size>0,(src,name)
        band_ref = 'c3-Beryllene_trilayer' if sid == 'st' else ref
        assert (REPO/'4.0_bandstructure'/band_ref/'KPOINTS').is_file()
        p=src.joinpath('POTCAR').read_text()
        zvals=[float(v) for v in re.findall(r'ZVAL\s*=\s*([\d.]+)',p)]
        assert zvals==([2.] if natom==nelect//2 else [2.,1.]),(ref,zvals)
    for sid,ref,nelect,natom,inversion in SPECS:
        src=REPO/'6.0_dielectric_selection'/ref
        band_ref = 'c3-Beryllene_trilayer' if sid == 'st' else ref
        dest=RUN/sid; dest.mkdir()
        dipol='0.5 0.5 0.09644275954187936' if sid=='beta_1h' else '0.5 0.5 0.10'
        entry=dict(id=sid,reference=str(src.relative_to(REPO)),expected_nelect=nelect,atom_count=natom,inversion_expected=inversion,source_files={n:dict(bytes=(src/n).stat().st_size,sha256=sha(src/n)) for n in ('CONTCAR','INCAR','POTCAR','KPOINTS','CHGCAR')})
        manifests.append(entry); write(dest/'manifest.json',json.dumps(entry,indent=2)+'\n')
        common=f'''SYSTEM = {sid} SOC topology; frozen dielectric-selection geometry
PREC = Accurate
ENCUT = 600
GGA = PE
ALGO = Normal
NELM = 240
EDIFF = 1E-8
ISMEAR = 0
SIGMA = 0.01
NBANDS = 28
NSW = 0
IBRION = -1
LREAL = .FALSE.
LASPH = .TRUE.
GGA_COMPAT = .FALSE.
LSORBIT = .TRUE.
LNONCOLLINEAR = .TRUE.
SAXIS = 0 0 1
MAGMOM = {3*natom}*0
NELECT = {nelect}
NCORE = 6
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL = {dipol}
LOPTICS = .FALSE.
NWRITE = 3
'''
        for stage in ('scf','trim','bands'):
            d=dest/stage; d.mkdir()
            shutil.copyfile(src/'CONTCAR',d/'POSCAR')
            shutil.copyfile(src/'POTCAR',d/'POTCAR')
            if stage=='scf':
                inc=common+'ISTART = 0\nICHARG = 1\nISYM = 2\nKPAR = 1\nLWAVE = .FALSE.\nLCHARG = .TRUE.\n'
                shutil.copyfile(src/'CHGCAR',d/'CHGCAR')
                write(d/'KPOINTS','Dense SOC charge and sampled gap\n0\nGamma\n105 105 1\n0 0 0\n')
            else:
                inc=common+'ISTART = 0\nICHARG = 11\nISYM = -1\nKPAR = 1\nLCHARG = .FALSE.\nLORBIT = 11\nLWAVE = '+('.TRUE.' if stage=='trim' else '.FALSE.')+'\n'
                if stage=='trim': write(d/'KPOINTS','All four 2D TRIM, explicit order Gamma X Y M\n4\nReciprocal\n0 0 0 1\n0.5 0 0 1\n0 0.5 0 1\n0.5 0.5 0 1\n')
                else: shutil.copyfile(REPO/'4.0_bandstructure'/band_ref/'KPOINTS',d/'KPOINTS')
            write(d/'INCAR',inc)
        body=f'''cd '{dest}'
for stage in scf trim bands; do
  cd '{dest}'/"$stage"
  test ! -e OUTCAR || {{ echo 'Refusing to overwrite an existing VASP run'; exit 20; }}
  if [ "$stage" != scf ]; then cp ../scf/CHGCAR CHGCAR; fi
  mpirun -np 42 '{VASP}' > vasp.log 2>&1
  validation_args=(--manifest '{dest}/manifest.json' --structure-id '{sid}')
  if [ "$stage" = trim ]; then validation_args+=(--require-wavecar); fi
  python3 -B '{RUN}/check_topology_outputs.py' . "${{validation_args[@]}}"
done
printf 'Completed %s\\n' '{sid}' > '{dest}/PIPELINE_COMPLETE'
'''
        write(dest/'run.pbs',pbs('topo_'+sid,body))
    write(RUN/'manifest.json',json.dumps(dict(created_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),reference_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),tools_prefix=str(TOOLS),structures=manifests),indent=2)+'\n')
    TOOLS.mkdir(parents=True)
    write(RUN/'modules.sh',MODULES)
    write(RUN/'build_tools.pbs',pbs('topo_tools',f"cd '{RUN}'\nbash build_tools.sh '{TOOLS}' > '{TOOLS}/build.log' 2>&1",ranks=12,mem='48gb',hours=12))
    write(RUN/'wcc.pbs',pbs('topo_wcc',f"cd '{RUN}'\n'{TOOLS}/venv/bin/python' -B wcc_driver.py '{RUN}/beta_1h' '{TOOLS}' > '{RUN}/beta_1h/wcc_driver.log' 2>&1"))
    write(RUN/'parity.pbs',pbs('topo_parity',f"cd '{RUN}'\n'{TOOLS}/venv/bin/python' -B parity_analysis.py '{RUN}'",ranks=1,mem='12gb',hours=4))
    print('PREPARED',RUN, 'TOOLS',TOOLS,flush=True)

if __name__=='__main__': main()
