#!/usr/bin/env python3
"""Independently compare every exported curve with its raw density response."""
from __future__ import annotations
import dataclasses
import gzip
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import h5py
import numpy as np
from scipy.signal import find_peaks
from vmatplot.structure import slab_optical_geometry

H_EV = 4.135667662e-15
C_NM = 2.99792458e17
OLD_THICKNESS = {'a-Beryllene':3.96, 'b-Beryllene':.0473771641975619+3.96,
 'c3-Beryllene':.0751622445053417+3.96, 'a-Beryllene_hh':.04361280878130169+1.4,
 'b-Beryllene_b':.07108860259770651+2.68, 'b-Beryllene_bb':.0948207854203246+1.4}


def spectra(directory, component, factor, legacy=False):
    with h5py.File(ROOT/directory/'vaspout.h5','r') as f:
        e = f['results/linear_response/energies_dielectric_function'][()]
        index = {'xx':0,'yy':1,'zz':2}[component]
        raw = f['results/linear_response/density_density_dielectric_function'][index,index,:,:]
    ratio = factor[0]/factor[1]
    eps = 1+ratio*(raw[:,0]-1)+1j*ratio*raw[:,1]
    nk = np.sqrt(eps.astype(complex))
    n,k=nk.real,nk.imag
    # Independent vacuum-wavelength relation; no library frequency helpers.
    inverse_wavelength=e/(H_EV*C_NM)
    absorption=4*np.pi*k*inverse_wavelength
    if legacy:
        absorption/=2*np.pi
    values={'eps1':eps.real,'eps2':eps.imag,'refractive':n,'extinction':k,
            'absorption':absorption,'reflectivity':np.abs((nk-1)/(nk+1))**2,
            'energy-loss':(-1/eps).imag}
    np.testing.assert_allclose(n*n-k*k,eps.real,rtol=3e-13,atol=3e-12)
    np.testing.assert_allclose(2*n*k,eps.imag,rtol=3e-13,atol=3e-12)
    return e,values


def metadata(directory):
    path=ROOT/directory
    with h5py.File(path/'vaspout.h5','r') as f:
        version='.'.join(str(int(f['version/'+k][()])) for k in ('major','minor','patch'))
        bands=int(f['results/electron_eigenvalues/nb_tot'][()])
        requested={k:f['input/incar/'+k][()].item() for k in ('NBANDS','CSHIFT','SIGMA','ISMEAR','LOPTICS')}
    data={'vasp_version':version,'actual_total_bands':bands,'requested_incar':requested}
    xml=path/'vasprun.xml'
    compressed=path/'vasprun_.xml.gz'
    if xml.exists() or compressed.exists():
        stream=open(xml,'rb') if xml.exists() else gzip.open(compressed,'rb')
        data['xml_source']=str((xml if xml.exists() else compressed).relative_to(ROOT))
        data['xml_sha256']=hashlib.sha256((xml if xml.exists() else compressed).read_bytes()).hexdigest()
        for event,node in ET.iterparse(stream,events=['end']):
            if node.tag=='parameters':
                tags=['ISPIN','LSORBIT','ISMEAR','SIGMA','CSHIFT','NBANDS','LOPTICS','WPLASMAI','WPLASMA']
                data['actual_xml_parameters']={key:' '.join(el.text.split())
                    for key in tags if (el:=node.find(f'.//*[@name="{key}"]')) is not None}
                break
        stream.close()
    return data


def main():
    directory=ROOT/'exported_figures'
    with gzip.open(directory/'validation/figure_data.json.gz','rt') as stream:
        records=json.load(stream)
    from regenerate_optics import sha256, input_hashes
    expected_inputs=json.loads((directory/'validation/input_hashes_before.json').read_text())
    assert expected_inputs==input_hashes(), 'Raw inputs changed after generation'
    manifest=json.loads((directory/'manifest.json').read_text())
    assert len(manifest)==34
    for item in manifest:
        expected=item['after_sha256']
        assert expected==sha256(directory/item['export'])==sha256(ROOT/item['source'])
        assert expected==records[item['source']]['sha256']
    environment=json.loads((directory/'validation/environment.json').read_text())
    for path,expected in {**environment['source_sha256'],**environment['notebook_sha256']}.items():
        assert sha256(ROOT/path)==expected, f'Source changed after regeneration: {path}'
    for path,record in records.items():
        assert sha256(ROOT/path)==record['sha256'], f'Figure changed after verification: {path}'
    checked_curves=0
    checked_points=0
    maximum=0
    per_file={}
    cache={}
    for path,record in records.items():
        systems={s[0]:s for s in record['systems']}
        components=[next(iter(c)) if isinstance(c,dict) else c for c in record['components']]
        count=0
        max_error=0
        for axis_index,ax in enumerate(record['axes']):
            if record['property']=='dielectric':
                aliases={alias:comp for item in record['components'] for comp,alias in item.items()}
                part,alias=ax['title'].split(' part for ',1)
                component=aliases[alias]
                property_key='eps1' if part=='Real' else 'eps2'
            else:
                component=components[axis_index]
                property_key=record['property']
            for curve in ax['curves']:
                system=systems[curve['label']]
                material_name=Path(system[1]).name
                expected_factor=(1,1) if material_name.startswith(('a-Beryllium','c0-Beryllium')) else slab_optical_geometry(ROOT/system[1]).factor
                np.testing.assert_allclose(system[6],expected_factor,rtol=1e-12,atol=1e-12,
                    err_msg=f'{path}: thickness factor disagrees with actual structure {system[1]}')
                key=(system[1],component,tuple(system[6]))
                if key not in cache:
                    cache[key]=spectra(system[1],component,system[6])
                e,values=cache[key]
                actual=np.asarray(curve['xy'])
                indices=np.searchsorted(e,actual[:,0])
                np.testing.assert_allclose(e[indices],actual[:,0],rtol=0,atol=1e-12)
                target=values[property_key][indices]
                np.testing.assert_allclose(actual[:,1],target,rtol=2e-10,atol=2e-11,
                                           err_msg=f'{path}: {component} {curve["label"]}')
                low,high=ax['ylim']
                assert actual[:,1].min()>=low-1e-10 and actual[:,1].max()<=high+1e-10
                error=float(np.max(np.abs(target-actual[:,1])))
                max_error=max(max_error,error)
                checked_points+=len(actual)
                checked_curves+=1
                count+=1
        maximum=max(maximum,max_error)
        per_file[path]={'curves':count,'max_absolute_difference':max_error,'no_y_clipping':True,
                        'status':record['status']}
    # The two A layouts must use bit-identical data and physical limits.
    from regenerate_optics import A_THESIS,A_PUBLICATION
    paired=0
    for a,b in zip(A_THESIS,A_PUBLICATION):
        left=records['figures_for_thesis/'+a]['axes']
        right=records['figures_for_publication/'+b]['axes']
        for l,r in zip(left,right):
            assert l['xlim']==r['xlim'] and l['ylim']==r['ylim']
            assert [c['xy'] for c in l['curves']]==[c['xy'] for c in r['curves']]
        paired+=1
    materials={}
    npz={}
    for path in sorted((ROOT/'6.0_dielectric_selection').iterdir()):
        if not path.is_dir():continue
        rel=path.relative_to(ROOT)
        geometry=slab_optical_geometry(path) if path.name in OLD_THICKNESS else None
        factor=geometry.factor if geometry else (1,1)
        old_factor=(40,OLD_THICKNESS[path.name]) if geometry else (1,1)
        material={'geometry':dataclasses.asdict(geometry) if geometry else None,
                  'calculation':metadata(rel),'components':{}}
        for component in ('xx','yy','zz'):
            e,new=spectra(rel,component,factor)
            _,old=spectra(rel,component,old_factor,legacy=True)
            mask=(e>=0)&(e<=12)
            energies=e[mask]
            component_data={'last_sample_eV':float(energies[-1]),'sample_count':len(energies),'properties':{}}
            for prop in new:
                values=new[prop][mask]
                previous=old[prop][mask]
                imax=int(np.argmax(values))
                old_imax=int(np.argmax(previous))
                peaks,_=find_peaks(values)
                component_data['properties'][prop]={
                    'old_max':float(previous[old_imax]),'old_max_at_eV':float(energies[old_imax]),
                    'new_max':float(values[imax]),'new_max_at_eV':float(energies[imax]),
                    'last_value':float(values[-1]),'max_at_boundary':imax in (0,len(values)-1),
                    'local_peaks':[{'energy_eV':float(energies[i]),'value':float(values[i])} for i in peaks]}
                npz[f'{path.name}__{component}__{prop}']=values
            npz[f'{path.name}__{component}__energy_eV']=energies
            material['components'][component]=component_data
            if path.name=='a-Beryllene':
                for prop in new:
                    np.testing.assert_allclose(new[prop],old[prop]*(2*np.pi if prop=='absorption' else 1),rtol=1e-12)
        if path.name=='b-Beryllene':
            ex=material['components']
            e,vals_x=spectra(rel,'xx',factor)
            _,vals_z=spectra(rel,'zz',factor)
            for label,fac in [('old',old_factor),('new',factor)]:
                _,xx=spectra(rel,'xx',fac)
                _,zz=spectra(rel,'zz',fac)
                select=(e>=0)&(e<=3)&(xx['eps1']*zz['eps1']<0)
                material[label+'_low_energy_opposite_sign_samples_eV']=e[select].tolist()
            material['relative_eps_minus_identity_scale']=old_factor[1]/factor[1]
        materials[path.name]=material
    report={'checked_exported_curves':checked_curves,'checked_sample_values':checked_points,
            'max_absolute_array_difference':maximum,'identical_A_layout_pairs':paired,
            'all_plotted_samples_inside_axes':True,'files':per_file,'materials':materials,
            'scope':'Original independent-particle density-density response; original surface-centre thickness convention; no new DFT, Drude, BSE or local-field correction.'}
    (directory/'validation/numeric_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    np.savez_compressed(directory/'validation/corrected_spectra_0-12eV.npz',**npz)
    print(f'Independent audit passed: {checked_curves} curves, {checked_points} samples; max error {maximum:.3g}')
    print(f'A layouts: {paired} identical pairs')


if __name__=='__main__':main()
