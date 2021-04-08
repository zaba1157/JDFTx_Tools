#!/usr/bin/env python3
"""
Created on Tue Mar 30 11:33:52 2021

@author: NSing
"""

from pymatgen.core.surface import Structure
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import numpy as np
from pymatgen.core.structure import Molecule
import os


def save_structures(st_list, location = './', skip_existing = False, single_loc = False):
    # st_list = list of pymatgen structures
    if not os.path.exists(location):
        os.mkdir(location)
    folders = []
    for i, st in enumerate(st_list):
        if not single_loc:
            root = os.path.join(location, str(i+1).zfill(2))
        else:
            root = location
        if not os.path.exists(root):
            os.mkdir(root)
        if skip_existing and os.path.isfile(os.path.join(root, 'POSCAR')):
            continue
        st.to('POSCAR', os.path.join(root, 'POSCAR'))
        folders.append(root)
    return folders

def setup_neb(initial, final, nimages, save_loc, linear = False):
    si = Structure.from_file(initial)
    sf = Structure.from_file(final)
    # setup interpolated structures
    structures = si.interpolate(sf, nimages+1, autosort_tol=0)
    if not linear:
        from pymatgen.io.ase import AseAtomsAdaptor
        from ase.neb import NEB
        structures_ase = [ AseAtomsAdaptor.get_atoms(struc) for struc in structures ]
        neb = NEB(structures_ase)
        neb.interpolate('idpp') # type: NEB
        structures = [ AseAtomsAdaptor.get_structure(atoms) for atoms in neb.images ]
    # save folders with POSCARs
    for i, st in enumerate(structures):
        sub_folder = os.path.join(save_loc, str(i).zfill(2))
        if not os.path.exists(sub_folder):
            os.mkdir(sub_folder)
        st.to('POSCAR',os.path.join(sub_folder, 'POSCAR'))
        if i == 0 or i == nimages-1:
            st.to('POSCAR',os.path.join(sub_folder, 'CONTCAR'))
    return True

def place_ads(loc, ads_sts, surface_st, mol, sites_allowed, 
              ads_distance = 2.0, min_dist = 0.5, freeze_depth = 1.8, _z_dir = 2):
    # place adsorbate on all sites within height
    if loc == 'All':
        temp_list = []
        for st in ads_sts:
            height = 4
            asf = AdsorbateSiteFinder(st, height = height)
            sites = asf.find_adsorption_sites(distance = ads_distance, 
                                              symm_reduce=0.05, near_reduce=0.05,
                                              positions = sites_allowed)
            for site in sites['all']:
                new_st = asf.add_adsorbate(mol, site)
                if any([any([ np.sqrt(np.sum([(x1.coords[i] - x2.coords[i])**2 for i in range(3)]))
                             < min_dist for x2 in new_st.sites if x2 != x1]) for x1 in new_st.sites]):
                    print('Distance Error in Adsorbate Adding')
                    continue
                temp_list.append(assign_selective_dynamics(new_st, freeze_depth))
        return temp_list
    # place adsorbate on highest atom of type "loc"
    elif type(loc) == str and loc != 'center':
        el = loc
        max_site = None
        for site in surface_st.sites:
            if site.species_string != el:
                continue
            if max_site is None or site.coords[_z_dir] > max_site.coords[_z_dir]:
                max_site = site
        if max_site is None:
            print('Adsorbate location not found! ', el)
            return []
        temp_list = []
        ads_shift = np.array([0,0,0])
        ads_shift[_z_dir] = ads_distance
        for st in ads_sts:
            asf = AdsorbateSiteFinder(st)
            site = max_site.coords + ads_shift #np.array([0, 0, ads_distance])
            new_st = asf.add_adsorbate(mol.copy(), site, reorient = True)
            if any([any([ np.sqrt(np.sum([(x1.coords[i] - x2.coords[i])**2 for i in range(3)]))
                             < min_dist for x2 in new_st.sites if x2 != x1]) for x1 in new_st.sites]):
                    print('Distance Error in Adsorbate Adding')
                    continue
            temp_list.append(assign_selective_dynamics(new_st, freeze_depth))
        return temp_list
    # place adsorbate on specific site
    elif type(loc) == int:
        temp_list = []
        ads_shift = np.array([0,0,0])
        ads_shift[_z_dir] = ads_distance
        for st in ads_sts:
            asf = AdsorbateSiteFinder(st)
            site = surface_st.cart_coords[loc-1] + ads_shift # VESTA indexes to 1:, convert to 0:
            new_st = asf.add_adsorbate(mol.copy(), site, reorient = True)
            if any([any([ np.sqrt(np.sum([(x1.coords[i] - x2.coords[i])**2 for i in range(3)]))
                             < min_dist for x2 in new_st.sites if x2 != x1]) for x1 in new_st.sites]):
                    print('Distance Error in Adsorbate Adding')
                    continue
            temp_list.append(assign_selective_dynamics(new_st, freeze_depth))
        return temp_list
    # place adsorbate at center of lattice
    elif loc == 'center':
        max_site = surface_st.sites[0]
        for site in surface_st.sites:
            if site.coords[_z_dir] > max_site.coords[_z_dir]:
                max_site = site
        temp_list = []
        ads_shift = np.array([0,0,0])
        ads_shift[_z_dir] = ads_distance
        st = ads_sts[0]
        asf = AdsorbateSiteFinder(st)
        site = np.array([st.lattice.a/2, st.lattice.b/2, max_site.coords[2] + ads_distance]) 
        new_st = asf.add_adsorbate(mol.copy(), site, reorient = True)
        if any([any([ np.sqrt(np.sum([(x1.coords[i] - x2.coords[i])**2 for i in range(3)]))
                         < min_dist for x2 in new_st.sites if x2 != x1]) for x1 in new_st.sites]):
                print('Distance Error in Adsorbate Adding')
                return []
        return [assign_selective_dynamics(new_st, freeze_depth)]

def add_adsorbates(surface_st, adsorbates, ads_distance = 2.0, sites_allowed = ['ontop', 'bridge','hollow'],
                   min_dist = 0.5, freeze_depth = 1.8, molecules_loc = '', z_dir = 2):
    adsorbate_sts = {}
#    for adss, locs in adsorbates.items():
#        ads_st = surface_st.copy()
#        ads_sts = [ads_st]
#        if '_' in adss:
#            for ia, ads in enumerate(adss.split('_')):
#                mol = molecule_from_poscar(ads, location=molecules_loc)
#                ads_sts = place_ads(locs[ia], ads_sts, mol, sites_allowed, 
#                                    ads_distance, min_dist, freeze_depth, z_dir)
#        else:
#            ads_sts = []
#            mol = molecule_from_poscar(adss, location=molecules_loc)
#            for loc in locs:
#                ads_sts += place_ads(loc, [ads_st], surface_st, mol, sites_allowed, 
#                                     ads_distance, min_dist, freeze_depth, z_dir)
#        adsorbate_sts[adss] = ads_sts
#    return adsorbate_sts
    for adss, locs in adsorbates.items():
        # TODO: add back in multi molecule adsorption
#        print('Adding adsorbate: ', adss)
        ads_st = surface_st.copy()
        ads_sts = []
        mol = molecule_from_poscar(adss, location=molecules_loc)
        for loc in locs:
            ads_sts += place_ads(loc, [ads_st], surface_st, mol, sites_allowed, 
                                 ads_distance, min_dist, freeze_depth, z_dir)
        adsorbate_sts[adss] = ads_sts
    return adsorbate_sts

def assign_selective_dynamics(slab, depth):
    min_depth = min([x.coords[2] for x in slab.sites])
    sd_list = []
    sd_list = [[False, False, False] if site.coords[2] - min_depth < depth
               else [True, True, True] for site in slab.sites]
    new_sp = slab.site_properties
    new_sp['selective_dynamics'] = sd_list
    return slab.copy(site_properties=new_sp)

def molecule_from_poscar(adsorbate, location = ''):
    """
    Creates pymatgen molecule object from file.
    Args:
        filename: name of molecule file - can be path
    """
    if location == '':
        location = 'molecules/'+adsorbate+'/POSCAR'
    st = Structure.from_file(location)
    atoms, coords = [], []
    lattice = np.array([st.lattice.a, st.lattice.b, st.lattice.c,])
    mid_point = np.array([st.lattice.a/2, st.lattice.b/2, st.lattice.c/2,])
    for site in st.sites:
        atoms.append(site.species_string)
        nc = []
        for i,c in enumerate(site.coords):
            if c > mid_point[i]:
                nc.append(c - lattice[i])
            else:
                nc.append(c)
        coords.append(nc)
    return Molecule(atoms, coords)

def data_analysis(all_data, ref_mols):
    '''
    Helper function to analyze converged structures from jdft_manager
    Analysis includes:
        1) getting adsorption energies from adsorbed calcs + surfs + molecules
        2) get adsorption site and adsorbate binding distance (ensure it is bound)
    TODO
        3) ddec analysis of ads site
        4) DOS analysis of surf vs. adsorbed (plots)
        5) plots: ads energy vs. bias; ads energy of mol on many surfs (w/wo bias)
    '''
    #MOL: all_data[mol_name][bias_str] = data
    #SURF: all_data[surf_name]['surf'][bias_str] = data
    #ADSORBED: all_data[surf_name]['adsorbed'][mol_name][bias_str][mol_config] = data
    #DESORBED: all_data[surf_name]['desorbed'][mol_name][bias_str] = data
    #NEB: all_data[surf_name]['neb'][mol_name][bias_str][neb_path] = data
    analysis = {'ads_energies': {},
                'surf_data': {}}
    # get adsorption energies
    for k, v in all_data.items():
        if k in ref_mols or k in ['converged']:
            continue # k is a molecule or converged list
        if 'adsorbed' not in v:
            continue # no adsorbed data
        for mol, mv in v['adsorbed'].items():
            if mol not in ref_mols:
                print('Molecule '+mol+' should be added to reference molecules dictionary. Contact Nick.')
                continue
            refs = ref_mols[mol]
            for bias, configs in mv.items():
                if not v['surf'][bias]['converged']:
                    continue # surface not converged at relevant bias
                for ref in refs['refs']:
                    if ref not in all_data or bias not in all_data[ref] or not all_data[ref][bias]['converged']:
                        continue # molecule not converged at respective bias
                ads_list = [] # list of ads energies for a given molecule, surface and bias
                min_ads = None
                surf_energy = v['surf'][bias]['final_energy']
                ref_energy = np.sum([refs['coeffs'][i] * all_data[r][bias]['final_energy'] 
                                    for i,r in enumerate(refs['refs'])])
                for config, cv in configs.items():
                    ads_en = cv['final_energy'] - surf_energy - ref_energy
                    if np.abs(ref_energy) > 2.0:
                        print('Warning: Adsorption energy for '+mol+' on '+k+' at '+bias+' is > 2.0,'+
                              ' check convergence.')
                    ads_list.append(ads_en)
                    if min_ads is None or ads_en < min_ads['energy']:
                        min_ads = {'energy': ads_en, 'config': config}
                # TODO: finish ads energy analysis 
    
    return all_data

def write_parallel(roots, cwd, total_cores, cores_per_job, time, out, shell_folder,
                   qos = None, nodes=1):
    # get all necessary inputs
    script = os.path.join(os.environ['JDFTx_Tools_dir'], 'run_JDFTx.py')
    try:
        modules=' '.join(os.environ['JDFTx_mods'].split('_'))
    except:
        modules=('comp-intel/2020.1.217 intel-mpi/2020.1.217 cuda/10.2.89 vasp/6.1.1 mkl/2020.1.217'+
                 ' gsl/2.5/gcc openmpi/4.0.4/gcc-8.4.0 gcc/7.4.0')
    try:
        comp=os.environ['JDFTx_Computer']
    except:
        comp='Eagle'
    alloc = None
    if comp == 'Eagle':
        try:
            alloc = os.environ['JDFTx_allocation']
        except:
            alloc = 'electrobuffs'
    partition = None
    if comp == 'Bridges2':
        partition = 'RM-shared'
    
    # create shell file
    writelines = '#!/bin/bash'+'\n'
    writelines+='#SBATCH -J '+out+'\n'
    if comp == 'Bridges2':
        writelines+='#SBATCH -t '+str(time)+':00:00'+'\n'
    else:
        writelines+='#SBATCH --time='+str(time)+':00:00'+'\n'
    writelines+='#SBATCH -o '+out+'-%j.out'+'\n'
    writelines+='#SBATCH -e '+out+'-%j.err'+'\n'
    
    if partition is not None:
        writelines+='#SBATCH -p '+partition+'\n'
    if alloc is not None:
        writelines+='#SBATCH --account='+alloc+'\n'

    if comp == 'Eagle':
        writelines+='#SBATCH --tasks '+str(nodes * total_cores)+'\n'
    writelines+='#SBATCH --nodes '+str(nodes)+'\n'
    writelines+='#SBATCH --ntasks-per-node '+str(total_cores)+'\n'

    if qos=='high' and comp == 'Eagle':
        writelines+='#SBATCH --qos=high'+'\n'

    if time == 1 and comp == 'Eagle':
        writelines+='#SBATCH --partition=debug\n'
    
    writelines+='\nexport JDFTx_NUM_PROCS='+str(cores_per_job)+'\n'
    writelines+='module load '+modules+'\n\n'

    for i, root in enumerate(roots):
        if i+1 < len(roots):
            add = ' &'
        else:
            add = ' && fg'
        writelines+=('python ' + script +' -d '+ os.path.join(cwd, root) + ' > '
                     + os.path.join(cwd, root, 'out_file') + add + '\n')
        
    writelines+='exit 0'+'\n'

    with open(os.path.join(shell_folder, out+'.sh'),'w') as f:
        f.write(writelines)

