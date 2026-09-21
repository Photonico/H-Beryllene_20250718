#!/usr/bin/env python3
"""Run the repository's optical notebooks and package checked thesis figures.

Run from any directory: python scripts/regenerate_optics.py
Use --baseline-source with a git-extracted vmatplot and 9.0 notebook to reproduce
legacy processing in an isolated output directory before applying corrections.
"""
from __future__ import annotations
import argparse
import base64
import ast
import gzip
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ['6.0_dielectric_function_convergence.ipynb', '6.1_dielectric_function.ipynb',
             '7.1_linear_optical_properties.ipynb', '7.2_liner_optical_properties_with_hydrogen.ipynb',
             '9.0_thesis_demo.ipynb', '9.1_thesis_demo_H.ipynb']
A_THESIS = ['3.1_dielectric.pdf', '4.1_absorption.pdf', '4.2_refractive.pdf',
            '4.3_extinction.pdf', '4.4_reflectivity.pdf', '4.5_energy-loss.pdf']
A_PUBLICATION = ['fig3.13_dielectric.pdf', 'fig3.14a_absorption.pdf', 'fig3.14b_refractive.pdf',
                 'fig3.15_extinction.pdf', 'S3.16_reflectivity.pdf', 'S3.17_energy-loss.pdf']
B_PUBLICATION = ['fig3.16_dielec.pdf', 'fig3.17a_abs.pdf', 'fig3.17b_energy-loss.pdf',
                 'S3.18a_refractive.pdf', 'S3.18b_reflectivity.pdf', 'S3.18c_extinction.pdf',
                 'fig3.18_dielec_H.pdf', 'fig3.19a_abs_H.pdf', 'fig3.19b_energy-loss_H.pdf',
                 'S3.19a_refractive_H.pdf', 'S3.19b_reflectivity_H.pdf', 'S3.19c_extinction_H.pdf']
C_PUBLICATION = [f'S3.{i}_dielec_{phase}.pdf' for i, phase in
                 zip(range(3,13), ['hcp','hcp','alpha','alpha','beta','beta','bcc','bcc','cubic','cubic'])]
UNCHANGED = {'figures_for_publication/S3.3_dielec_hcp.pdf',
             'figures/6_dielectric_test/6.10_dielec_a-Beryllium_nbands.pdf'}


def sha256(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def input_hashes():
    return {str(p.relative_to(ROOT)): sha256(p)
            for folder in ('6.0_dielectric_selection', '6.0_dielectric_function')
            for p in sorted((ROOT/folder).rglob('*'))
            if p.is_file() and p.name in {'CONTCAR','POSCAR','INCAR','KPOINTS','vaspout.h5','vasprun.xml'}}


def describe_figure(fig):
    import numpy as np
    from matplotlib.colors import to_rgb
    fig.canvas.draw()
    axes = []
    for ax in fig.axes:
        curves = []
        for line in ax.lines:
            if not line.get_visible() or line.get_label().startswith('_'):
                continue
            xy = np.asarray(line.get_xydata())
            page = ax.transData.transform(xy) * 72 / fig.dpi
            page[:,1] = fig.get_figheight()*72-page[:,1]
            curves.append({'label': line.get_label(), 'color': list(to_rgb(line.get_color())),
                           'xy': xy.tolist(), 'page_xy': page.tolist()})
        axes.append({'title': ax.get_title(), 'xlim': list(ax.get_xlim()),
                     'ylim': list(ax.get_ylim()), 'curves': curves})
    return {'size_inches': fig.get_size_inches().tolist(), 'axes': axes}


def verify_pdf(path, description, simplified=False):
    """Check exported vector vertices against plotted numerical arrays, in points."""
    import fitz
    import numpy as np
    with fitz.open(path) as doc:
        if len(doc) != 1:
            raise AssertionError(f'{path}: expected one page')
        drawings = doc[0].get_drawings()
        if simplified:
            from verify_optics_exports import numeric_spans, calibrate_axis
            spans=numeric_spans(doc[0])
            frames=[d for d in drawings if d['type']=='f' and len(d['items'])==1
                    and d['items'][0][0]=='re' and d['rect'].width>100
                    and d['rect'].height>100 and d['rect']!=doc[0].rect]
            if len(frames)!=len(description['axes']):
                raise AssertionError(f'{path}: unexpected preserved PDF axis count')
    checked = 0
    points = 0
    maximum = 0.0
    for axis_index,ax in enumerate(description['axes']):
        axis_drawings=drawings
        if simplified:
            frame=frames[axis_index]
            stop=frames[axis_index+1]['seqno'] if axis_index+1<len(frames) else float('inf')
            axis_drawings=[d for d in drawings if frame['seqno']<d['seqno']<stop]
            sx,tx,sy,ty,_=calibrate_axis(frame['rect'],axis_drawings,spans)
            ax['xlim']=sorted([(frame['rect'].x0-tx)/sx,(frame['rect'].x1-tx)/sx])
            ax['ylim']=sorted([(frame['rect'].y0-ty)/sy,(frame['rect'].y1-ty)/sy])
        for curve in ax['curves']:
            expected = np.asarray(curve['page_xy'])
            if simplified:
                xy=np.asarray(curve['xy'])
                expected=np.column_stack([sx*xy[:,0]+tx,sy*xy[:,1]+ty])
                curve['page_xy']=expected.tolist()
            candidates = []
            for drawing in axis_drawings:
                if drawing['color'] is None or not np.allclose(drawing['color'], curve['color'], atol=2e-6):
                    continue
                segments = [item for item in drawing['items'] if item[0] == 'l']
                if not segments or (not simplified and len(segments) != len(expected)-1):
                    continue
                vertices = np.array([tuple(segments[0][1])] + [tuple(item[2]) for item in segments])
                if simplified:
                    from scipy.spatial import cKDTree
                    if len(vertices)<3 or np.max(np.abs(vertices[[0,-1]]-expected[[0,-1]]))>0.01:
                        continue
                    residual=float(cKDTree(expected).query(vertices)[0].max())
                    # Check the omitted samples against the exported polyline too.
                    distances=[]
                    for point in expected:
                        delta=vertices[1:]-vertices[:-1]
                        length2=np.sum(delta*delta,axis=1)
                        fraction=np.clip(np.sum((point-vertices[:-1])*delta,axis=1)/np.maximum(length2,1e-30),0,1)
                        distances.append(np.linalg.norm(point-(vertices[:-1]+fraction[:,None]*delta),axis=1).min())
                    if max(distances)>0.25:
                        continue
                    candidates.append(residual)
                else:
                    candidates.append(float(np.max(np.abs(vertices-expected))))
            error = min(candidates, default=float('inf'))
            if error > 0.01:
                raise AssertionError(f'{path}: curve {curve["label"]!r} PDF residual {error} pt')
            maximum = max(maximum, error)
            checked += 1
            points += len(expected)
    return {'curves': checked, 'sample_vertices': points, 'max_residual_pt': maximum,
            'simplified_legacy_pdf': simplified,
            'omitted_sample_polyline_tolerance_pt': 0.25 if simplified else 0.0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-source', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    baseline = args.baseline_source is not None
    if baseline:
        args.baseline_source = args.baseline_source.resolve()
    args.output = (args.output.resolve() if args.output else
                   ROOT/'.agent/optics_audit/baseline' if baseline else ROOT/'exported_figures')
    if baseline and args.output == ROOT/'exported_figures':
        parser.error('Baseline output must be isolated from exported_figures')
    sys.path.insert(0, str(args.baseline_source if baseline else ROOT))
    os.chdir(ROOT)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import nbformat
    if not baseline:
        from vmatplot.plot_validation import assert_curve_visibility
    matplotlib.rcParams['path.simplify'] = False
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    validation = output/'validation'
    validation.mkdir(exist_ok=True)
    hashes = input_hashes()
    hash_file = validation/'input_hashes_before.json'
    if hash_file.exists() and json.loads(hash_file.read_text()) != hashes:
        raise RuntimeError('Inputs changed since the first run; inspect before regenerating')
    hash_file.write_text(json.dumps(hashes, indent=2)+'\n')
    originals = {str(p.relative_to(ROOT)): hashlib.sha256(subprocess.check_output(
        ['git','show','HEAD:'+str(p.relative_to(ROOT))])).hexdigest()
        for directory in ('figures_for_thesis','figures_for_publication')
        for p in (ROOT/directory).glob('*.pdf')}
    for rel in UNCHANGED:
        originals[rel] = hashlib.sha256(subprocess.check_output(['git','show','HEAD:'+rel])).hexdigest()
    records = {}
    current_notebook = None
    current_cell = None
    active_namespace = {}
    current_property = None
    original_savefig = plt.savefig

    def checked_savefig(filename, *save_args, **kwargs):
        rel = Path(filename).as_posix()
        fig = plt.gcf()
        if not baseline:
            assert_curve_visibility(fig)
        description = describe_figure(fig)
        destination = output/rel if baseline else ROOT/rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        unchanged = not baseline and rel in UNCHANGED
        if unchanged and sha256(destination) != originals[rel]:
            raise AssertionError(f'{rel}: expected unchanged Git PDF was modified')
        if not unchanged:
            original_savefig(destination, *save_args, **kwargs)
        check = None if baseline else verify_pdf(destination, description, simplified=unchanged)
        description['systems'] = active_namespace.get('systems', [])
        description['components'] = active_namespace.get('components', [])
        description['property'] = current_property
        records[rel] = {'notebook':current_notebook, 'cell':current_cell, 'sha256':sha256(destination),
                        'status':'verified_unchanged' if unchanged else 'regenerated',
                        'pdf_verification':check, **description}
        print(f'{records[rel]["status"]}: {rel}', file=sys.__stdout__, flush=True)

    plt.savefig = checked_savefig
    notebooks = ['9.0_thesis_demo.ipynb'] if baseline else NOTEBOOKS
    for name in notebooks:
        current_notebook = name
        source = (args.baseline_source if baseline else ROOT)/name
        nb = nbformat.read(source, as_version=4)
        namespace = {'__name__':'__main__'}
        active_namespace = namespace
        execution_count = 0
        for index, cell in enumerate(nb.cells):
            if cell.cell_type != 'code':
                continue
            current_cell = index
            current_property = None
            for node in ast.walk(ast.parse(cell.source)):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                    if node.func.id.startswith('plot_dielectric'):
                        current_property = 'dielectric'
                    elif node.func.id == 'plot_linear_optical_property':
                        current_property = ast.literal_eval(node.args[2])
            execution_count += 1
            plt.close('all')
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                exec(compile(cell.source, f'{name}:cell{index}', 'exec'), namespace)
            cell.execution_count = execution_count
            cell.outputs = []
            if stream.getvalue():
                cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=stream.getvalue()))
            for number in plt.get_fignums():
                fig = plt.figure(number)
                image = io.BytesIO()
                fig.savefig(image,format='png',dpi=90)
                cell.outputs.append(nbformat.v4.new_output('display_data',data={
                    'image/png':base64.b64encode(image.getvalue()).decode('ascii'),
                    'text/plain':repr(fig)},metadata={}))
        if not baseline:
            nbformat.write(nb, source)
    plt.close('all')
    after = input_hashes()
    if hashes != after:
        raise AssertionError('Raw calculation inputs changed during regeneration')
    with gzip.open(validation/'figure_data.json.gz', 'wt', encoding='utf-8') as stream:
        json.dump(records, stream, ensure_ascii=False, separators=(',',':'))
    (validation/'input_hashes_after.json').write_text(json.dumps(after, indent=2)+'\n')
    environment = {'python':sys.version, 'head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                   'path_simplify':False, 'raw_inputs_unchanged':True}
    import importlib.metadata
    environment['packages']={p:importlib.metadata.version(p) for p in ['numpy','scipy','matplotlib','h5py','nbformat','PyMuPDF']}
    source_root=args.baseline_source if baseline else ROOT
    environment['source_sha256']={str(p.relative_to(source_root)):sha256(p) for p in sorted((source_root/'vmatplot').glob('*.py'))}
    environment['notebook_sha256']={name:sha256(source_root/name) for name in notebooks}
    (validation/'environment.json').write_text(json.dumps(environment, indent=2)+'\n')
    if baseline:
        return
    manifest=[]
    for group, folder, names, target_folder in [
        ('A','figures_for_thesis',A_THESIS,'figures_ch3'),
        ('A','figures_for_publication',A_PUBLICATION,'figures_proj3'),
        ('B','figures_for_publication',B_PUBLICATION,'figures_proj3'),
        ('C','figures_for_publication',C_PUBLICATION,'figures_proj3')]:
        for name in names:
            rel=f'{folder}/{name}'
            target=output/folder/name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(ROOT/rel,target)
            record=records[rel]
            manifest.append({'group':group,'source':rel,'export':rel,'thesis_target':f'{target_folder}/{name}',
                             'status':record['status'],'notebook':record['notebook'],'cell':record['cell'],
                             'before_sha256':originals[rel],'after_sha256':sha256(target),
                             'pdf_verification':record['pdf_verification']})
    (output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(f'Packaged {len(manifest)} figures into {output}')


if __name__=='__main__':
    main()
