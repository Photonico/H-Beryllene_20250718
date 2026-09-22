#!/usr/bin/env python3
"""Commit only the submitted campaign and additive large-file ignore rules."""
import os, subprocess
from pathlib import Path
RUN=Path(__file__).resolve().parent
REPO=RUN.parent.parent
LIMIT=99_000_000

def git(*args):
    return subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()

def main():
    assert (RUN/'SUBMISSION_COMPLETE').is_file(), 'Submit every intended job first'
    assert len((RUN/'submitted_jobs.jsonl').read_text().splitlines())==9
    assert not git('diff','--cached','--name-only'), 'Existing staged changes require separate handling'
    assert not git('diff','--','.gitignore'), 'Preserve unrelated .gitignore edits; review manually'
    ignore=REPO/'.gitignore'
    large=[]
    for base,dirs,files in os.walk(REPO):
        dirs[:]=[d for d in dirs if d!='.git']
        for name in files:
            p=Path(base)/name
            if p.is_symlink(): continue
            if p.stat().st_size>LIMIT: large.append(p.relative_to(REPO).as_posix())
    def escape(path):
        return ''.join('\\'+c if c in '\\ !#[]*?' else c for c in path)
    rules='\n# Topology campaign: ignore files larger than 99,000,000 bytes.\n'
    rules+='\n'.join('/'+escape(p) for p in sorted(large))+'\n'
    with ignore.open('a') as f: f.write(rules)
    small=[RUN/'.gitignore',RUN/'README.md',RUN/'manifest.json',RUN/'submitted_jobs.jsonl',RUN/'SUBMISSION_COMPLETE']
    small+=list(RUN.glob('*.py'))+list(RUN.glob('*.sh'))+list(RUN.glob('*.pbs'))
    for d in RUN.iterdir():
        if d.is_dir() and (d/'manifest.json').is_file():
            small += [d/'manifest.json',d/'run.pbs']
            for stage in ('scf','trim','bands'):
                small += [d/stage/n for n in ('POSCAR','INCAR','KPOINTS')]
    small.append(ignore)
    assert all(p.is_file() and p.stat().st_size<=LIMIT for p in small)
    relative=sorted(set(str(p.relative_to(REPO)) for p in small))
    subprocess.run(['git','add','--',*relative],cwd=REPO,check=True)
    staged=git('diff','--cached','--name-only').splitlines()
    assert set(staged)<=set(relative), 'Unexpected staged path'
    for path in staged:
        size=int(git('cat-file','-s',':'+path))
        assert size<=LIMIT,(path,size)
        assert not path.endswith(('/POTCAR','/WAVECAR','/CHGCAR'))
    print('Ignored large paths:',len(large),'Staged campaign paths:',len(staged),flush=True)
    subprocess.run(['git','diff','--cached','--stat'],cwd=REPO,check=True)
    subprocess.run(['git','commit','-m','submit topology calculations'],cwd=REPO,check=True)
    print('COMMIT',git('rev-parse','HEAD'),flush=True)

if __name__=='__main__': main()
