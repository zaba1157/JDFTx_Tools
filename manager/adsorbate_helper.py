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
import json
from itertools import permutations
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcdefaults()

hartree_to_ev = 27.2114
hartree_to_ev = 1

reference_molecules = {'H': {'refs': ['H2'], 'coeffs': [0.5]},
                       'H2': {'refs': ['H2'], 'coeffs': [1]},
#                       'H2_des': {'refs': ['H2'], 'coeffs': [1]},
                       'H2O': {'refs': ['H2O'], 'coeffs': [1]},
                       'H3O':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]},
#                       'H3O_des':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]},
                       'H_H2O':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]},
                       'H_H3O':{'refs': ['H2O', 'H2'], 'coeffs': [1, 1]},
                       'H2_H2O':{'refs': ['H2O', 'H2'], 'coeffs': [1, 1]},
                       'O': {'refs': ['H2O','H2'], 'coeffs': [1,-1]},
                       'CO2': {'refs': ['CO2'], 'coeffs': [1]},
                       'CO': {'refs': ['CO2','O'], 'coeffs': [1,-1]},
                       'CHO': {'refs': ['CO','H'], 'coeffs': [1,1]},
                       'COH': {'refs': ['CHO'], 'coeffs': [1]},
                       'OCH': {'refs': ['CHO'], 'coeffs': [1]},
                       'N': {'refs': ['N2'], 'coeffs': [0.5]},
                       'N2': {'refs': ['N2'], 'coeffs': [1]},
                       'N2H': {'refs': ['N2','H'], 'coeffs': [1,1]},
                       'NH': {'refs': ['N2','H'], 'coeffs': [0.5,1]},
                       'NH2': {'refs': ['N2','H'], 'coeffs': [0.5,2]},
                       'NH3': {'refs': ['N2','H'], 'coeffs': [0.5,3]},
                       'OC': {'refs': ['CO'], 'coeffs': [1,]},
                       'OCO': {'refs': ['CO2'], 'coeffs': [1]},
                       'OH': {'refs': ['O','H'], 'coeffs': [1,1]},
                       'OOH': {'refs': ['O','H'], 'coeffs': [2,1]},
                       'S2': {'refs': ['S8'], 'coeffs': [0.25]},
                       'S4': {'refs': ['S8'], 'coeffs': [0.5]},
                       'S6': {'refs': ['S8'], 'coeffs': [0.75]},
                       'S8': {'refs': ['S8'], 'coeffs': [1]},
                       }

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
#            print('ADS DIST:' + str(ads_distance))
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

def data_analysis(all_data, ref_mols, __ads_warning_dist = 2.5, verbose = True):
    '''
    Helper function to analyze converged structures from jdft_manager
    Analysis includes:
        1) getting adsorption energies from adsorbed calcs + surfs + molecules
        2) Surface charge density over bias
    TODO
        2) get adsorption site and adsorbate binding distance (ensure it is bound)
        2.1) get ads distance and # of bonds (within 0.5A of min length?)
        3) get charge density difference between ads site and surf site
        4) DOS analysis of surf vs. adsorbed (plots)            *** TODO: add DOS sp calc option to manager
        5) plots: ads energy vs. bias; ads energy of mol on many surfs (w/wo bias)
    '''
    #MOL: all_data[mol_name][bias_str] = data
    #SURF: all_data[surf_name]['surf'][bias_str] = data
    #ADSORBED: all_data[surf_name]['adsorbed'][mol_name][bias_str][mol_config] = data
    #DESORBED: all_data[surf_name]['desorbed'][mol_name][bias_str] = data
    #NEB: all_data[surf_name]['neb'][mol_name][bias_str][neb_path] = data
    
    #ANALYSIS: analysis[surf][DATA]
    #ANALYSIS: analysis[surf]['ads'][mol][bias][DATA]
    analysis = {}
    
    if verbose: print('\n----- ANALYSIS -----')
    skip_surfs= ['Mo6S8_100_h2o']
    # get adsorption data
    for surf, sv in all_data.items():
        if surf in skip_surfs:
            continue
        if surf in ref_mols or surf in ['converged']:
            continue # k is a molecule or converged list
        
        # initialize surface analysis data
        if surf not in analysis:
            analysis[surf] = {'surf': {}}
        if verbose: print('\nSURF: '+surf)
        
        analysis[surf]['surf']['net_charge'] = {}
        analysis[surf]['surf']['el_charge'] = {}
        for bias, bv in sv['surf'].items():
            analysis[surf]['surf']['net_charge'][bias] = bv['net_oxidation']
            analysis[surf]['surf']['el_charge'][bias] = {}
            for site, sitev in bv['site_data'].items():
                atom = sitev['atom']
                if atom in analysis[surf]['surf']['el_charge'][bias]:
                    analysis[surf]['surf']['el_charge'][bias][atom] += [sitev['oxi_state']]
                else:
                    analysis[surf]['surf']['el_charge'][bias][atom] = [sitev['oxi_state']]
            for atom in analysis[surf]['surf']['el_charge'][bias]:
                analysis[surf]['surf']['el_charge'][bias][atom] = np.average(
                        analysis[surf]['surf']['el_charge'][bias][atom])
        if verbose: print('Added surface charge data.')
        
        # no adsorbed data
        if 'adsorbed' not in sv:
            continue 
        
        # add adsorbate data
        for mol, mv in sv['adsorbed'].items():
            if mol not in ref_mols:
                print('Molecule '+mol+' should be added to reference molecules dictionary. Contact Nick.')
                continue
                        
#            refs = ref_mols[mol]
            refs = get_ref_mol_dic(mol, ref_mols) 
            no_ref = False
            for bias, configs in mv.items():
                if not sv['surf'][bias]['converged']:
                    continue # surface not converged at relevant bias
                for ref in refs:
                    if ref not in all_data or bias not in all_data[ref] or not all_data[ref][bias]['converged']:
                        no_ref = True # molecule not converged at respective bias
                if no_ref:
                    continue
                ads_list = {} # list of ads energies for a given molecule, surface and bias
                min_ads = None
                surf_energy = sv['surf'][bias]['final_energy']
                surf_data = sv['surf'][bias]
                ref_energy = np.sum([coef * all_data[ref][bias]['final_energy'] 
                                    for ref, coef in refs.items()])
                
                for config, cv in configs.items():
                    if not cv['converged']:
                        continue
                    
                    surf_comp = check_compatible(cv['inputs'], sv['surf'][bias]['inputs'])
                    if not surf_comp:
                        print("WARNING: Surface "+surf+" and adsorbed "+mol+"("+bias+
                              ") inputs not compatible")
                        continue
                    for r in refs:
                        mol_comp = check_compatible(cv['inputs'], all_data[r][bias]['inputs'])
                        if not mol_comp:
                            print("WARNING: Surface "+surf+" and adsorbed "+r+
                                  "("+bias+") inputs not compatible")
                            continue
                    
                    ads_en = (cv['final_energy'] - surf_energy - ref_energy) / hartree_to_ev
#                    print('ADS:',surf, mol, bias, ads_en)
                    ads_data = get_ads_data(surf_data, cv)
                    
                    if np.abs(ads_en) > 2.0:
                        print('Warning: Adsorption energy for '+mol+' on '+surf+' at '+bias+' is > 2.0 eV,'+
                              ' check convergence.')
                    ads_list[config] = ads_en
                    if min_ads is None or ads_en < min_ads['energy']:
                        min_ads = {'energy': ads_en, 'config': config, 'ads_data': ads_data}
                
                if min_ads is None:
                    continue
                if 'ads' not in analysis[surf]:
                    analysis[surf]['ads'] = {}
                if mol not in analysis[surf]['ads']:
                    analysis[surf]['ads'][mol] = {}
                if bias not in analysis[surf]['ads'][mol]:
                    analysis[surf]['ads'][mol][bias] = {}
                analysis[surf]['ads'][mol][bias]['all_energies'] = {'min_energy': min_ads, 'all_ens': ads_list,
                                                         '__Note__': 'ads energies (eV) vs ref mols',}
                analysis[surf]['ads'][mol][bias]['ads_energy'] = min_ads['energy']
                if min_ads['ads_data'] != False:
                    analysis[surf]['ads'][mol][bias]['ads_data'] = min_ads['ads_data']
                    analysis[surf]['ads'][mol][bias]['bond_len'] = min_ads['ads_data']['bond_len']
                    analysis[surf]['ads'][mol][bias]['is_bound'] = min_ads['ads_data']['is_bound']
                    analysis[surf]['ads'][mol][bias]['bond_charge_diff'] = (
                            min_ads['ads_data']['bond_charge_diff'])
                    analysis[surf]['ads'][mol][bias]['surf_site_charge'] = min_ads['ads_data']['surf_site_charge']
        if verbose: print('Added adsorbate data.')
        
        if 'desorbed' in sv:
            pass
    
    return analysis

def get_shared_site(main_site, all_sites, threshold = 0.75):
    for ind, site in all_sites.items():
        if site['atom'] != main_site['atom']:
            continue
        if all([np.abs(float(p) - float(site['positions'][i])) < threshold 
                for i,p in enumerate(main_site['positions'])]):
            return True, site
    return False, None   

def get_ads_data(surf, ads, lattice_diff = 0.01, max_bond_len = 2.5):
    # returns: surf_ads_site, bound_atom, bond_len, n_bonds, surf_jitter
    #          surf_site_charge, surf_charge_diff, bond_charge_diff
    surf_st = Structure.from_dict(surf['contcar'])
    ads_st = Structure.from_dict(ads['contcar'])
    surf_sites = surf['site_data']
    ads_sites = ads['site_data']
    
    if not (np.abs(surf_st.lattice.a - ads_st.lattice.a) < lattice_diff and
            np.abs(surf_st.lattice.b - ads_st.lattice.b) < lattice_diff and
            np.abs(surf_st.lattice.c - ads_st.lattice.c) < lattice_diff):
        print('ERROR: Surface and adsorbate have different lattice parameters.')
        return False
    
    surf_jitter = {}
    mol_sites = {ind: site for ind, site in ads_sites.items() if
                 not get_shared_site(site, surf_sites)[0]}
    ads_dist = None
    all_bonds = []
    for ind, site in ads_sites.items():
        shared = get_shared_site(site, surf_sites)
        if not shared[0]:
            continue
        shared_site = shared[1]
        surf_jitter[ind] = [float(p) - float(shared_site['positions'][i]) 
                            for i,p in enumerate(site['positions'])]
        
        # get all ads distances
        site_pos = np.array([float(p) for p in site['positions']])
        mol_dists = []
        mol_dists = {mi: np.linalg.norm(np.array([float(p) for p in mol['positions']]) - site_pos) 
                     for mi, mol in mol_sites.items()}
        min_mol_dist = min(list(mol_dists.values()))
        min_mol_site = mol_sites[[mi for mi, dist in mol_dists.items() if dist == min_mol_dist][0]]
        # if this is min_dist, set sites
        if ads_dist is None or min_mol_dist < ads_dist:
            ads_dist = min_mol_dist
            ads_surf_site = site
            init_surf_site = shared_site
            ads_mol_site = min_mol_site
        if min_mol_dist < max_bond_len:
            all_bonds += [site]
        
    if ads_dist is None:
        print('ERROR: ads_data is None')
        print(len(mol_sites))
        return False
    return {'bond_len': ads_dist, 'surf_site': ads_surf_site, 'mol_site': ads_mol_site,
            'n_bonds': len(all_bonds), 'all_bond_sites': all_bonds, 'is_bound': ads_dist < max_bond_len,
            'surf_site_rearrangement': surf_jitter, 'surf_site_charge': ads_surf_site['oxi_state'],
            'surf_charge_diff': ads_surf_site['oxi_state'] - init_surf_site['oxi_state'],
            'bond_charge_diff': np.abs(ads_surf_site['oxi_state'] - ads_mol_site['oxi_state'])}

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
        
def check_compatible(inputs1, inputs2, tags = ['fluid','elec-ex-corr','fluid-solvent',
                                               'pcm-variant','pseudos','elec-cutoff']):
    for t in tags:
        if (t in inputs1 and t not in inputs2) or (t not in inputs1 and t in inputs2):
            return False
    if not all([inputs1[t] == inputs2[t] for t in tags if t in inputs1]):
        return False
    return True

def set_rc_params():
    params = {'axes.linewidth' : 1.5,'axes.unicode_minus' : False,
              'figure.dpi' : 100,
              'font.size' : 18,'font.family': 'sans-serif','font.sans-serif': 'Verdana',
              'legend.frameon' : False,'legend.handletextpad' : 0.2,'legend.handlelength' : 0.9,
              'legend.fontsize' : 14,
              'mathtext.default' : 'regular','savefig.bbox' : 'tight',
              'xtick.labelsize' : 16,'ytick.labelsize' : 16,
              'xtick.major.size' : 6,'ytick.major.size' : 6,
              'xtick.major.width' : 1.5,'ytick.major.width' : 1.5,
              'xtick.top' : True,'xtick.bottom' : True,'ytick.right' : True,'ytick.left' : True,
              'xtick.direction': 'out','ytick.direction': 'out','axes.edgecolor' : 'black'}
    for p in params:
        mpl.rcParams[p] = params[p]
    return params

def get_ref_mol_dic(mol, reference_molecules):
    ref_dic = {}
    for ir, ref_mol in enumerate(reference_molecules[mol]['refs']):
        ref_c = reference_molecules[mol]['coeffs'][ir]
        if reference_molecules[ref_mol]['refs'] == [ref_mol]:
            if ref_mol in ref_dic:
                ref_dic[ref_mol] += ref_c
            else:
                ref_dic[ref_mol] = ref_c
        else:
            mrd = get_ref_mol_dic(ref_mol, reference_molecules)
            for k,v in mrd.items():
                if k in ref_dic:
                    ref_dic[k] += v * ref_c
                else:
                    ref_dic[k] = v * ref_c
    return ref_dic

def plot_scaling(analysis, mol1, mol2, width_height = (5,5), bias = 'all'):
    set_rc_params()
    plt.figure(figsize=width_height)
    colors = ['r','b','g','c','m','y','k']
    index = 0
    for surf, data in analysis.items():
        if mol1 not in data or mol2 not in data:
            continue
        c = colors[index % len(colors)]
        index += 1
        dt1 = data[mol1]
        dt2 = data[mol2]
        legend = True
        for bias1 in dt1:
            if bias1 not in dt2:
                continue
            if bias is not 'all':
                if type(bias) == str and bias != bias1:
                    continue
                elif type(bias) == list and bias1 not in bias:
                    continue
            if legend:
                plt.plot(dt1[bias1]['ads_energy'], dt2[bias1]['ads_energy'], label = surf,
                     color = c, marker = 'o')
                legend = False
            else:
                plt.plot(dt1[bias1]['ads_energy'], dt2[bias1]['ads_energy'], color = c, marker = 'o')
    plt.legend()
    plt.xlabel(mol1)
    plt.ylabel(mol2)
    plt.show()

def plot_tafel(analysis, mol, width_height = (5,5)):
    set_rc_params()
    mpl.rcParams['legend.frameon'] = True
    plt.figure(figsize=width_height)
    colors = ['r','b','g','c','m','y','k']
    bias_all = ['No_bias', '0.00V', '-0.25V', '-0.50V']
    allowed_surfs = ['Pt_111','Cu_111','Ni_111','Pd_111','Ag_111','Au_111']
    
    index = 0
    for surf, data in analysis.items():
        if surf not in allowed_surfs:
            continue
        if 'ads' not in data or mol not in data['ads']:
            continue
        c = colors[index % len(colors)]
        index += 1
        dt = data['ads'][mol]
        
        biases = [b for b in bias_all if b in dt and dt[b]['ads_energy'] < 2.0]
        energies = [dt[b]['ads_energy'] for b in bias_all if b in dt and dt[b]['ads_energy'] < 2.0]
        plt.plot(biases, energies, color = c, label = surf, ls = '-', marker = 'o')
        
#    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.07),
               ncol=4, fancybox=True, shadow=True)
    plt.xlabel('Bias (V)')
    plt.ylabel('Ads. Energy (eV)')
    
#    plt.ylim((-0.5, 1.1))
    plt.show()

def minimum_movement_strs(init, final, same_threshold = 0.4):
    isites = init.sites
    fsites = final.sites
    fnew = [0] * len(isites)
    used_sites = []
    incomplete = []
    for i,site in enumerate(isites):
        el = site.species_string
        added = False
        for ii,fs in enumerate(fsites):
            if fs.species_string != el or ii in used_sites:
                continue
            dist = fs.distance_from_point(site.coords)
            if dist < same_threshold:
                fnew[i] = fs
                used_sites.append(ii)
                added = True
                break
        if not added:
            incomplete.append(i)
            
    if len(incomplete) == 0:
        print('Inital and Final structures were successfully sorted.')
        return init, Structure.from_sites(fnew)
    if len(incomplete) > 7:
        print('WARNING: 8 or more sites are different between init and final ('+
              str(len(incomplete))+'). Check structures!')
        return init, final
    
    # 1 or more sites moved > same_threshold, go by nearest distance score!
    unused = [i for i,fs in enumerate(fsites) if i not in used_sites]
    if len(unused) != len(incomplete):
        print('METAERROR: variables unused and incomplete are different lengths.')
        return init, final
    # generate all possible ion pairs, find pair with minimum total distance
    all_sets = [list(zip(x,unused)) for x in permutations(incomplete,len(unused))]
    dists = []
    for setx in all_sets:
        dist = 0
        for pair in setx:
            si = isites[pair[0]]
            sf = fsites[pair[1]]
            if si.species_string != sf.species_string:
                dist += 1000
            dist += si.distance_from_point(sf.coords)
        dists.append(dist)
    min_dist_set = all_sets[dists.index(min(dists))]
    for pair in min_dist_set:
        if fnew[pair[0]] != 0:
            print('METAERROR: writing over existing site in fnew.')
        fnew[pair[0]] = fsites[pair[1]]
    if any([type(ff) == int for ff in fnew]):
        print('METAERROR: Not all sites in fnew were assigned from init.')
        return init, final
    print('Inital and Final structures were successfully sorted.')
    print('**** Inital and Final structures re-written! ****')
    return init, Structure.from_sites(fnew)


if __name__ == '__main__':
#    assert False
    jdft_data = {}
    data_files = [#'jdft_manager/backup_updates/Eagle/all_data.json',
                  #'jdft_manager/backup_updates/Bridges/all_data.json',
                  #'jdft_manager/backup_updates/Summit/all_data.json',
                  'all_data.json']
    for file in data_files:
        with open(file,'r') as f:
            add_data = json.load(f)
        for k,v in add_data.items():
            jdft_data[k] = v
    
    analysis = data_analysis(jdft_data, reference_molecules)
    
#    plot_scaling(analysis, 'H', 'CO', bias = 'all')
#    plot_tafel(analysis, 'H')
    