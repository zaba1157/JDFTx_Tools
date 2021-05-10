#!/usr/bin/env python3

import os
import argparse
import subprocess
import json
from pymatgen.core.structure import Structure
import numpy as np
from adsorbate_helper import save_structures, add_adsorbates, data_analysis, setup_neb, write_parallel

reference_molecules = {'H': {'refs': ['H2'], 'coeffs': [0.5]},
                       'H2': {'refs': ['H2'], 'coeffs': [1]},
                       'H2O': {'refs': ['H2O'], 'coeffs': [1]},
                       'H3O':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]}, # TODO: run calc of this above Pt(111)
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

class jdft_manager():
    '''
    Author: Nick Singstock
    Date: March 30th 2021
    
    Class jdft_manager allows for management of a large number of JDFTx surface calculations 
    in a properly setup folder using the manager_control.txt file. Please run the "-s True" 
    command in a new folder to setup a managed folder and obtain the important subdirectories
    and files. More information is available in the manager_control.txt file.
    
    Please email Nick with questions about proper use and to alert about bugs.
    '''
    
    def __init__(self):
        # initialize parameters
        self.cwd = os.getcwd()
        self.__get_user_inputs__()
        self.__get_run_cmd__()
        self.initialize_vars()
        self.proper_setup = self.check_setup()
        
        if self.args.setup == 'True':
            print('Setting up initial directory for management.')
            self.__setup__()
    
    def read_inputs(self, folder):
        # returns a dictionary of inputs
        with open(os.path.join(folder, 'inputs'), 'r') as f:
            input_text = f.read()
        tags = {}
        for line in input_text.split('\n'):
            if line == '':
                continue
            if line[0] == '#':
                continue
            words = line.split(' ')
            if words[0] in tags:
                if type(tags[words[0]]) != list:
                    tags[words[0]] = [tags[words[0]]]
                tags[words[0]].append(' '.join(words[1:]))
            else:
                tags[words[0]] = ' '.join(words[1:])
        return tags
    
    def write_inputs(self, inputs_dic, folder):
        # expects a dictionary of inputs with either strings or lists of strings as vals
        text = ''
        for k,v in inputs_dic.items():
            if type(v) == list:
                for vi in v:
                    text += k + ' ' + vi + '\n'
            else:
                text += k + ' ' + v + '\n'
        with open(os.path.join(folder, 'inputs'), 'w') as f:
            f.write(text)
    
    def read_optlog(self, folder):
        try:
            with open(os.path.join(folder, 'opt.log'), 'r') as f:
                opt_text = f.read()
        except:
            return False
        steps = []
        ionic_step = 0
        for i, line in enumerate(opt_text.split('\n')):
            if ' Step ' in line or '*Force-consistent' in line:
                continue
            if line == '':
                continue
            var = line.split()
            step_dic = {'opt_method': var[0][:-1],
                        'energy': float(var[3][:-1]),
                        'force': float(var[4]),
                        'step': ionic_step}
            ionic_step += 1
            steps.append(step_dic)
        return steps
    
    def read_outfile(self, folder):
        try:
            with open(os.path.join(folder, 'out'), 'r', errors='ignore') as f:
                out_text = f.read()
        except:
            print('Error reading out file.')
            return False
        site_data = {}
        record_forces, record_ions = False, False
        el_counter = 0
        net_oxidation = 0
        net_mag = 0
        for line in out_text.split('\n'):
            if 'Ionic positions in cartesian coordinates' in line:
                record_ions = True
                # reset net oxidation states and magnetization for new group
                net_oxidation = 0
                net_mag = 0
                continue
            if record_ions:
                if 'ion' in line:
                    atom, xpos, ypos, zpos = line.split()[1:5]
                    site_data[str(el_counter)] = {'atom': atom, 'positions': [xpos, ypos, zpos]}
                    el_counter += 1
                else: #if line == '' or 'ion' not in line:
                    record_ions = False
                    el_counter = 0
            if 'Forces in Cartesian coordinates' in line:
                record_forces = True
                continue
            if record_forces:
                if 'force' in line:
                    atom, xfor, yfor, zfor = line.split()[1:5]
                    if atom != site_data[str(el_counter)]['atom']:
                        print('Error reading out file.')
                        return False
                    site_data[str(el_counter)]['forces'] =  [xfor, yfor, zfor]
                    el_counter += 1
                else: #if line == '' or 'force' not in line:
                    record_forces = False
                    el_counter = 0
            if 'magnetic-moments' in line:
                el = line.split()[2]
                mags = line.split()[3:]
                mag_counter = 0
                for site,sv in site_data.items():
                    if sv['atom'] == el:
                        try:
                            site_data[site]['mag_mom'] = float(mags[mag_counter])
                            net_mag += float(mags[mag_counter])
                            mag_counter += 1
                        except:
                            print('Error reading magnetic moments')
                            return False
            if 'oxidation-state' in line:
                el = line.split()[2]
                oxis = line.split()[3:]
                oxi_counter = 0
                for site,sv in site_data.items():
                    if sv['atom'] == el:
                        try:
                            site_data[site]['oxi_state'] = float(oxis[oxi_counter])
                            net_oxidation += float(oxis[oxi_counter])
                            oxi_counter += 1
                        except:
                            print('Error reading oxidation states')
                            return False
#        print('Site data read')
#        site_data['net_oxidation'] = net_oxidation
#        site_data['net_magnetization'] = net_mag
        return site_data, net_oxidation, net_mag
    
    def read_contcar(self, folder):
        st = Structure.from_file(os.path.join(folder, 'CONTCAR'))
        return st.as_dict()
            
    def get_force(self, steps):
        return steps[-1]['force']
    
    def check_convergence(self, inputs, steps):
        force = self.get_force(steps)
        fmax = float(inputs['fmax'])
        if force <= fmax:
            return True
        return False
    
    def __setup__(self, verbose = True, overwrite = False):
        # make necessary subfolders
        for dir_to_make in [calc_folder, molecule_folder, inputs_folder, results_folder]:
            if not os.path.exists(dir_to_make):
                os.mkdir(dir_to_make)
            else:
                print('Directory: '+dir_to_make+' already exists. No changes made.')
        for cs in self.calc_subfolders:
            folder = os.path.join(calc_folder, cs)
            if not os.path.exists(folder):
                os.mkdir(folder)
            else:
                print('Directory: '+dir_to_make+' already exists. No changes made.')
        if verbose: print('Head directories created.')
        
        # copy default files to molecules and inputs
        def_inputs_folder = os.path.join(defaults_folder, 'inputs')
        for file in os.listdir(def_inputs_folder):
            file_loc = os.path.join(def_inputs_folder, file)
            self.run('cp ' + file_loc + ' ' + os.path.join(inputs_folder, file))
        
        def_mols_folder = os.path.join(defaults_folder, 'molecules')
        for folder in os.listdir(def_mols_folder):
            if folder in os.listdir(molecule_folder):
                print('Molecule: '+folder+' already exists. No changes made.')
                continue
            folder_loc = os.path.join(def_mols_folder, folder)
            self.run('cp -r ' + folder_loc + ' ' + os.path.join(molecule_folder, folder))
        if verbose: print('Default input files and molecules added.')
        
        if overwrite or not os.path.exists('./manager_control.txt'):
            self.run('cp ' + os.path.join(defaults_folder + 'manager_control.txt') + ' ./manager_control.txt')
        print('\nSuccessfully setup directory! Please see manager_control.txt for help.\n')
    
    def check_setup(self):
        if not all([os.path.exists(f) for f in [calc_folder, molecule_folder, 
                    inputs_folder, results_folder]]):
#            print('1')
            return False
#        for subdir in os.listdir(calc_folder):
        if any([x not in os.listdir(calc_folder) for x in self.calc_subfolders]):
#            print('2')
            return False
        return True
    
    def run(self, cmd):
        subprocess.call(cmd, shell=True)
        
    def initialize_vars(self):
        self.calc_subfolders = ['surfs', 'molecules', 'adsorbed', 'desorbed', 'neb']
        self.data_file = os.path.join(results_folder, 'all_data.json')
        
    
    def __get_user_inputs__(self):
        # get all user inputs from command line
        parser = argparse.ArgumentParser()
        
        parser.add_argument('-s', '--setup', help='Setup downstream folders for management. (-s True)',
                            type=str, default='False')
        parser.add_argument('-cc', '--check_calcs', help='Check convergence of all managed calculations. '+
                            'Default True.',type=str, default='True')
        parser.add_argument('-u', '--rerun_unconverged', help='Rerun all unconverged calculations being managed, '+
                            'requires "check_calcs". Default True.',type=str, default='True')
        parser.add_argument('-v', '--save', help='Save all newly processed data, requires "check_calcs".'+
                            ' Default True.',type=str, default='True')
        parser.add_argument('-a', '--analyze', help='Runs analysis on converged calcs, requires "save".'+
                            ' Default False.',type=str, default='False')
        parser.add_argument('-m', '--make_new', help='Make new calculations based on requested calcs.'+
                            ' Default True.',type=str, default='True')
        parser.add_argument('-rn', '--run_new', help='Run all newly setup calculations, requires "make_new".'+
                            ' Default True.',type=str, default='True')
        parser.add_argument('-ads', '--add_adsorbed', help='Add all requested adsorbates to converged surfs, '
                            +'requires "make_new". Default True.',type=str, default='True')
        parser.add_argument('-dist', '--adsorbate_distance', help='Standard distance from surface to adsorbate,'
                            +' requires "add_adsorbed". Default 2.0',type=float, default=2.0)
        parser.add_argument('-des', '--add_desorbed', help='Add all requested desorbed calcs to converged surfs, '
                            +'requires "make_new" and "add_adsorbed". Needed for NEB. Default True.',
                            type=str, default='True')
        parser.add_argument('-mol', '--add_molecules', help='Add all requested molecules, requires "make_new".'+
                            ' Default True.',type=str, default='True')
        parser.add_argument('-neb', '--make_neb', help='Makes NEB calculations from adsorbed+desorbed, '+
                            'requires "make_new". Default False.',type=str, default='True')
        parser.add_argument('-nebc', '--neb_climbing', help='If True, uses NEB climbing image. Requires'+
                            ' "make_neb". Default False.',type=str, default='False')
        parser.add_argument('-t', '--run_time', help='Time to run jobs. Default 48 (hours).',
                            type=int, default=48)
        parser.add_argument('-n', '--nodes', help='Nodes per job. Default 1.',
                            type=int, default=1)
        parser.add_argument('-c', '--cores', help='Cores per node.',
                            type=int, default=core_architecture)
        parser.add_argument('-r', '--short_recursive', help='Run jobs recursively on short queue until complete.'+
                            ' Very helpful when queue is busy. Default False.',type=str, default='False')
        parser.add_argument('-ra', '--read_all', help='Read all folders for new data. Does not use convergence '+
                            'file to speed up reading. Default False.', type=str, default='False')
        parser.add_argument('-rhe', '--rhe_zeroed', help='If True (default), converts all biases to be zeroed '+
                            'at 0V vs. RHE rather than 0V vs. SHE (if False).', type=str, default='True')
#        parser.add_argument('-ph', '--ph_rhe', help='pH for calculating SHE bias from RHE.'+
#                            ' Default 7.0. Be careful and consistent if changing! Requires rhe_zeroed = True',
#                            type=float, default=7.0)
        parser.add_argument('-q', '--qos', help='Whether qos should be high (True) or standard. Default False.',
                            type=str, default='False')
        parser.add_argument('-p', '--parallel', help='Runs multiple calcs together on a single node. Input'+
                            ' should be max number (int) of calcs to run together per node. Default 1.',
                            type=int, default=1)
        parser.add_argument('-b', '--backup', help='Whether to backup calcs folder. Default False.',
                            type=str, default='False')
        self.args = parser.parse_args()
    
    def __get_run_cmd__(self):
        self.run_cmd = run_command
        # add user inputs
        self.run_cmd += ' -t '+str(self.args.run_time)
        self.run_cmd += ' -n '+str(self.args.nodes)
        self.run_cmd += ' -c '+str(self.args.cores)
        if self.args.qos == 'True':
            self.run_cmd += ' -q high'
        if self.args.short_recursive == 'True':
            self.run_cmd += ' -r True'
    
    def get_bias(self, bias_str):
        if bias_str in ['No_bias']:
            return 'No_bias'
        return float(bias_str[:-1])
    
    def get_bias_str(self, bias):
        if bias == 'No_bias':
            return 'No_bias'
        return '%.2f'%bias + 'V'
    
    def scan_calcs(self, all_data, running_dirs, verbose = True):
        '''
        Main function for scanning through all previously-created sub-directories
        Functions independently from manager_control
        Functions:
            1) Scans through "calcs" directory 
            2) Saves data for all converged calculations
            3) Lists unconverged directories for rerunning
            4) Lists directories without inputs or CONTCARs to be started
        '''
        # look through all calc folders for converged calcs, unconverged calcs, and calcs to setup
        add_inputs = []
        run_new = []
        rerun = []
        running_parallel = self.get_parallel_running()
        if 'converged' not in all_data:
            all_data['converged'] = []
        if verbose: print('\n----- Scanning Through Calculations -----')
        for root, folders, files in os.walk(calc_folder):
            if 'POSCAR' not in files and 'inputs' not in files:
                continue
            if 'neb' in root and len(root.split(os.sep)) >= 7:
                # ignore neb subdirs
                continue
            if '__' in root:
                continue
            if verbose: print('\nPOSCAR/inputs found at:', root)
            if root in all_data['converged']:
                if verbose: print('Previously Converged.')
                continue
            full_root = os.path.join(self.cwd, root)
            if full_root in running_dirs:
                if verbose: print('Currently Running.')
                continue
            
            # get type of calculation
            calc_type = None
            for subf in self.calc_subfolders:
                tag = os.path.join(calc_folder, subf)
                if tag in root:
                    calc_type = subf
                    continue
            if calc_type is None:
                print('Error: No calc_type found for root: '+root+' Skipping.')
                continue
#            sub_dirs = root.split(tag)[-1].split(os.sep)
            sub_dirs = root.split(os.sep)
            
            if len(sub_dirs) < 4:
                if verbose: print('Not a calculation directory.')
                continue
            if 'inputs' not in files:
                add_inputs.append(root)
                if verbose: print('Adding inputs.')
                continue
            if 'CONTCAR' not in files and calc_type != 'neb':
                run_new.append(root)
                if verbose and self.args.run_new == 'True': print('Running unstarted job.')
                continue
            
            if calc_type == 'neb' and len(sub_dirs) > 6:
                print('BUG CHECKING: this is an neb subdir: ' + root)
                continue
            
            # read calc data at root
            if calc_type not in ['neb']:
                data = self.read_data(root)
            
            # save molecule data
            if calc_type == 'molecules':
                if verbose: print('Molecule calc read.')
                mol_name = sub_dirs[2]
                bias_str = sub_dirs[3]
                bias = self.get_bias(bias_str)
                if mol_name not in all_data:
                    all_data[mol_name] = {}
                data['bias'] = bias
                all_data[mol_name][bias_str] = data
                if data['converged']:
                    all_data['converged'].append(root)
                    if verbose: print('Molecule calc converged.')
                else:
                    rerun.append(root)
                    if verbose: print('Molecule calc not converged. Adding to rerun.')
                continue
            
            # save surface data
            elif calc_type == 'surfs':
                if verbose: print('Surface calc read.')
                surf_name = sub_dirs[2]
                bias_str = sub_dirs[3]
                bias = self.get_bias(bias_str)
                if surf_name not in all_data:
                    all_data[surf_name] = {}
                if 'surf' not in all_data[surf_name]:
                    all_data[surf_name]['surf'] = {}
                data['bias'] = bias
                all_data[surf_name]['surf'][bias_str] = data
                if data['converged']:
                    all_data['converged'].append(root)
                    if verbose: print('Surface calc converged.')
                else:
                    rerun.append(root)
                    if verbose: print('Surface calc not converged. Adding to rerun.')
                continue
            
            # save adsorbate and desorbed state calcs
            elif calc_type in ['adsorbed', 'desorbed']:
                # dic = surf: calc_type(s): mol: biases: configs: data. no configs for desorbed
                if verbose: print('Adsorbed/Desorbed calc read.')
                surf_name = sub_dirs[2]
                mol_name = sub_dirs[3]
                bias_str = sub_dirs[4]
                mol_config = None
                if calc_type == 'adsorbed': 
                    mol_config = sub_dirs[5]
                bias = self.get_bias(bias_str)
                if surf_name not in all_data or 'surf' not in all_data[surf_name]:
#                    print('Surface: '+surf_name+' must be created before adsorbed/desorbed calcs can be read!')
                    continue
#                print(surf_name)
#                if surf_name not in all_data:
#                    all_data[surf_name] = {}
                if calc_type not in all_data[surf_name]:
                    all_data[surf_name][calc_type] = {}
                if mol_name not in all_data[surf_name][calc_type]:
                    all_data[surf_name][calc_type][mol_name] = {}
                if bias_str not in all_data[surf_name][calc_type][mol_name]:
                    all_data[surf_name][calc_type][mol_name][bias_str] = {}
                    
                if bias_str in all_data[surf_name]['surf'] and all_data[surf_name]['surf'][bias_str]['converged']:
                    data['bias'] = bias
                    if calc_type == 'adsorbed':
                        all_data[surf_name][calc_type][mol_name][bias_str][mol_config] = data
                    elif calc_type == 'desorbed':
                        all_data[surf_name][calc_type][mol_name][bias_str] = data
                    if data['converged']:
                        all_data['converged'].append(root)
                        if verbose: print('Adsorbed/Desorbed calc converged.')
                    else:
                        rerun.append(root)
                        if verbose: print('Adsorbed/Desorbed calc not converged. Adding to rerun.')
                else:
                    print('Surface: '+surf_name+' at bias '+bias_str+
                          ' must be converged before adsorbed/desorbed can be saved.')
                continue
                    
            elif calc_type == 'neb':
                if verbose: print('NEB calc read.')
                # dic = surf: 'neb': mol: bias: path:
                surf_name = sub_dirs[2]
                mol_name = sub_dirs[3]
                bias_str = sub_dirs[4]
                neb_path = sub_dirs[5]
                bias = self.get_bias(bias_str)
                if surf_name not in all_data or 'surf' not in all_data[surf_name]:
                    continue
                if calc_type not in all_data[surf_name]:
                    all_data[surf_name][calc_type] = {}
                if mol_name not in all_data[surf_name][calc_type]:
                    all_data[surf_name][calc_type][mol_name] = {}
                if bias_str not in all_data[surf_name][calc_type][mol_name]:
                    all_data[surf_name][calc_type][mol_name][bias_str] = {}
                neb_data = self.get_neb_data(root, bias)
                all_data[surf_name][calc_type][mol_name][bias_str][neb_path] = neb_data
                if neb_data['converged']:
                    print('NEB calc for '+surf_name+' and '+mol_name+' at '+bias_str+' converged.')
                else:
                    print('NEB calc for '+surf_name+' and '+mol_name+' at '+bias_str+' not converged.'
                          +' Added to rerun.')
                    rerun.append(root)
                continue
        return all_data, add_inputs, rerun, run_new
    
    def read_data(self, folder):
        # currently reads inputs, opt_log for energies, and CONTCAR. Also checks convergence based on forces
        # reads out file for oxidation states and magentic moments
        inputs = self.read_inputs(folder)
        opt_steps = self.read_optlog(folder)
        if opt_steps == False:
            return {'opt': 'None', 'inputs': inputs, 'converged': False,
                    'final_energy': 'None', 'contcar': 'None'}
        # check if calc has high forces
        if opt_steps[-1]['force'] > 10:
            print('WARNING: High forces (> 10) in current step! May be divergent.')
        # check for convergence
        conv = self.check_convergence(inputs, opt_steps)
        contcar = 'None'
        if 'CONTCAR' in os.listdir(folder):
            contcar = self.read_contcar(folder)
        out_sites = self.read_outfile(folder)
        if out_sites == False:
            sites = {}
            net_oxi = 'None'
            net_mag = 'None'
        else:
            sites = out_sites[0]
            net_oxi = out_sites[1]
            net_mag = out_sites[2]
        return {'opt': opt_steps,
                'inputs': inputs,
                'converged': conv,
                'final_energy': 'None' if not conv else opt_steps[-1]['energy'],
                'contcar': contcar,
                'site_data': sites, 'net_oxidation': net_oxi, 'net_magmom': net_mag}
    
    def get_neb_data(self, folder, bias):
        # reads neb folder and returns data as a dictionary
        inputs = self.read_inputs(folder)
        opt_steps = self.read_optlog(folder)
        if opt_steps == False:
            return {'opt': 'None', 'inputs': inputs, 'converged': False,
                    'final_energy': 'None', 'images': {}}
        # check if calc has high forces
        if opt_steps[-1]['force'] > 10:
            print('WARNING: High forces (> 10) in current step! May be divergent.')
        # check for convergence
        conv = self.check_convergence(inputs, opt_steps)
        images = {}
        for root, folders, files in os.walk(folder):
            if 'CONTCAR' not in files or 'opt.log' not in files:
                continue
            image_num = root.split(os.sep)[-1]
            contcar = self.read_contcar(root)
            energy = self.read_optlog(root)[-1]['energy']
            images[image_num] = {'contcar': contcar, 'energy': energy}
        return {'opt': opt_steps,
                'inputs': inputs,
                'converged': conv,
                'final_energy': 'None' if not conv else opt_steps[-1]['energy'],
                'images': images}
    
    def analyze_data(self, data, ref_mols):
        '''
        Main function for analyzing converged data from scan_calcs function
        Functions:
            1) Creates analyzed.json file containing:
                - binding energies mapped over biases for each system
                - NEB barriers mapped over biases for each NEB system
        '''
        print('Data analysis not yet available. Please contact Nick to add.')
        return data_analysis(data)
    
    def rerun_calcs(self, rerun):
        print('\n----- Rerunning unconverged calcs -----\n')
        
        for root in rerun:
            os.chdir(root)
#            inputs = self.read_inputs('./')
#            inputs['restart'] = 'True'
#            self.write_inputs(inputs, './')
            self.run(self.run_cmd + ' -o '+self.get_job_name(root))
            os.chdir(self.cwd)
            
    def update_rerun(self, rerun):
        for root in rerun:
            os.chdir(root)
            self.failed_rerun_fixer(auto_delete = False)
            inputs = self.read_inputs('./')
            inputs['restart'] = 'True'
            self.write_inputs(inputs, './')
            os.chdir(self.cwd)
    
    def failed_rerun_fixer(self, auto_delete = False):
        '''
        Tries to fix errors that show up when rerunning a calc that previously failed.
        Current fixes:
            1) Length of fillings is incorrect
                Fix: delete fillings and ??? (just fillings doesn't seem to fix)
        '''
        with open('out','r') as f:
            outf = f.read()
        end_lines = outf.split('\n')[-10:]
        if 'Failed.' not in end_lines:
            return True
        for line in end_lines:
            # fillings is wrong size, delete
            if "Length of 'fillings' was" in line:
                if not auto_delete:
                    print('fillings is incorrect size, remove fillings and ???')
                else:
                    self.run('rm fillings')
        return True
                
    
    def get_job_name(self, root):
        return '-'.join(root.split(os.sep)[1:])
    
    def make_new_calcs(self, converged):
        '''
        Main management function for creating and upgrading calculations based on manager_control.txt
        Main functionalities:
            1) Upgrades surf calcs from No_bias -> 0V and 0V to other biases
            2) Adds requested molecules as adsorbates on converged surfaces at same bias
                2.1) Creates single point calculations of molecules above converged surface
                     at same bias. These are desorbed calculations and are needed for NEB.
                2.2) Runs molecules at same bias and solvent to allow for binding energy analysis
            3) Sets up NEB calculations from converged adsorbed+desorbed calculations at same bias
        '''
        # read manager_control.txt
        with open('manager_control.txt', 'r') as f:
            mc = f.read()
        print('\n----- Manager Control -----')
        # get dictionary of managed calcs
        managed_calcs = self.read_manager_control(mc)
        if managed_calcs == False:
            return False
        print('Manager control file successfully read.')
        # setup sub dirs based on converged and managed calcs
        new_folders = self.setup_managed_calcs(managed_calcs, converged)
        return new_folders
    
    def setup_managed_calcs(self, managed, converged, ads_distance = None, 
                            des_distance = 8, desorbed_single_point = True):
        '''
        Sets up all new calculations based on inputs from manager_control.txt
        '''
        if ads_distance is None:
            ads_distance = self.args.adsorbate_distance
        # creates new calc folders with POSCARs and inputs 
        # depends on args: add_adsorbed, add_desorbed, add_molecules, make_neb
        new_roots = []
        managed_mols = []
        # setup molecule folders in calc_folder
        if self.args.add_molecules == 'True':
            for mol in managed['molecules']:
                # ref_mols used for binding analysis, mol used for desorb SP calcs
                ref_mols = self.get_ref_mols(mol) + [mol]
                for ref_mol in ref_mols:
#                    if ref_mol in managed_mols:
#                        continue # It seems like this is taken care of by os.path.exists(bias_dir)
                    biases = list(set(managed['molecules'][mol]['biases']))
                    mol_location = self.get_mol_loc(ref_mol)
                    if mol_location == False:
                        continue
                    if not os.path.exists(os.path.join(calc_folder, 'molecules', ref_mol)):
                        os.mkdir(os.path.join(calc_folder, 'molecules', ref_mol))
                    for bias in biases:
                        bias_dir = os.path.join(calc_folder, 'molecules', ref_mol, self.get_bias_str(bias))
                        if bias_dir in converged:
                            continue
                        if os.path.exists(bias_dir):
                            # do not setup folders that exist, they are already setup
                            continue
                        os.mkdir(bias_dir)
                        # add POSCAR to bias folder
                        self.run('cp '+mol_location+' '+os.path.join(bias_dir, 'POSCAR'))
                        # add inputs to bias folder
                        self.run('cp '+os.path.join(inputs_folder, 'molecules_inputs')+' '
                                 +os.path.join(bias_dir, 'inputs'))
                        tags = ['target-mu '+ ('None' if bias in ['None','none','No_bias'] 
                                else '%.4f' % self.get_mu(bias, self.read_inputs(bias_dir)))]
                        self.add_tags(bias_dir, tags)
                        new_roots.append(bias_dir)
                        print('Reference molecule '+ref_mol+' at bias '+self.get_bias_str(bias)+
                              ' properly setup for mol '+mol)
                    managed_mols.append(ref_mol)
        
        # setup managed surfs, adsorbates, desorbed states, and NEB jobs
        for surf,v in managed.items():
            # 1) add new surfs and perform surf upgrading
            if surf in ['molecules'] or 'biases' not in v:
                continue
            print('\nManagement of surface: ' + surf)
            surf_roots = [os.path.join(calc_folder, 'surfs', surf, self.get_bias_str(bias)) 
                          for bias in v['biases']]
            for i, root in enumerate(surf_roots):
                # check if root is converged, if so it can be skipped
                if root in converged:
                    continue
                # check if root exists with POSCAR, if so it is being managed by rerun_calcs or add_inputs
                if os.path.exists(os.path.join(root, 'POSCAR')):
                    continue
                # root does not exist yet, check on bias dependency
                bias = v['biases'][i]
                if type(bias) == float and bias != 0.0:
                    # bias is a non-zero voltage, ensure 0 is converged
                    zero_root = os.path.join(calc_folder, 'surfs', surf, self.get_bias_str(0.0))
                    if zero_root not in converged:
                        # 0 V not yet converged! 
                        continue
                    # upgrade from 0 V (which exists!)
                    self.upgrade_calc(root, zero_root, bias, v['tags'] if 'tags' in v else [])
                    new_roots.append(root)
                elif bias == 0.0:
                    # bias is zero, ensure no-mu is converged
                    nm_root = os.path.join(calc_folder, 'surfs', surf, 'No_bias') 
                    if nm_root not in converged:
                        continue
                    self.upgrade_calc(root, nm_root, bias, v['tags'] if 'tags' in v else [])
                    new_roots.append(root)
                elif bias == 'No_bias':
                    # initial setup of no-mu surface
                    head_root = os.path.join(calc_folder, 'surfs', surf)
                    if not os.path.exists(os.path.join(head_root, 'POSCAR')):
                        print('POSCAR must be added to folder: '+head_root)
                        continue
                    else:
                        if not os.path.exists(root):
                            os.mkdir(root)
                        self.run('cp '+os.path.join(head_root, 'POSCAR')+' '+os.path.join(root, 'POSCAR'))
                    if not self.check_surface(os.path.join(root, 'POSCAR')):
                        continue
                    # copy inputs from inputs_folder and update based on tags
                    self.run('cp '+os.path.join(inputs_folder, 'surfs_inputs')+' '+os.path.join(root, 'inputs'))
                    if 'tags' in v:
                        self.add_tags(root, v['tags'])
                    new_roots.append(root)
                
            # 2) add adsorbates to converged surfaces at same bias
            if self.args.add_adsorbed == 'True':
                for mol, mv in v.items():
                    if mol in ['biases','tags']:
                        continue
                    if mol not in managed_mols:
                        print('Cannot add adsorbate, molecule '+mol+' not setup correctly.')
                        continue
                    # add each molecule requested at each bias on surface of corresponding bias
                    # mol = name of molecule to add, mv is dict with 'sites', 'biases' and 'tags' (optional)
                    
                    # create all adsorbate roots to make
                    if 'biases' not in mv:
                        print('Error: No biases found for mol '+mol+' for surf '+surf)
#                    ads_roots = [os.path.join(calc_folder, 'adsorbed', surf, mol, self.get_bias_str(bias)) 
#                          for bias in mv['biases']]
                    if not os.path.exists(os.path.join(calc_folder, 'adsorbed', surf)):
                        os.mkdir(os.path.join(calc_folder, 'adsorbed', surf))
                    if not os.path.exists(os.path.join(calc_folder, 'adsorbed', surf, mol)):
                        os.mkdir(os.path.join(calc_folder, 'adsorbed', surf, mol))
                    for bias in mv['biases']:
                        surf_root = os.path.join(calc_folder, 'surfs', surf, self.get_bias_str(bias))
                        if surf_root not in converged:
                            continue
                        # surf at same bias exists
                        st_surf = Structure.from_file(os.path.join(surf_root, 'CONTCAR'))
                        ads_dic = {mol: mv['sites']}
                        head_folder = os.path.join(calc_folder, 'adsorbed', surf, mol, self.get_bias_str(bias))
                        
                        # add molecule as adsorbate at all requested destinations
                        surf_ads = add_adsorbates(st_surf, ads_dic, ads_distance = ads_distance)
                        # save any new folders created, *** ANY EXISITING FOLDERS ARE IGNORED *** 
                        # does not use converged but has same functionality
                        # skipping earlier will miss new sites
                        save_locs = save_structures(surf_ads[mol], head_folder, skip_existing = True)
                        # add surfs to setup_new
                        if len(save_locs) == 0:
                            continue
                        for sl in save_locs:
                            print('Added adsorbate folder: '+sl)
                            self.run('cp '+os.path.join(inputs_folder, 'adsorbed_inputs')
                                     +' '+os.path.join(sl, 'inputs'))
                            # tags is a list
                            tags = mv['tags'] if 'tags' in mv else []
                            if 'tags' in v: tags += v['tags']
#                            tags += 'target-mu '+ ('None' if bias in ['None','none','No_bias'] else '%.2f'%bias)
                            tags += ['target-mu '+ ('None' if bias in ['None','none','No_bias'] 
                                     else '%.4f' % self.get_mu(bias, self.read_inputs(sl), tags))]
                            self.add_tags(sl, tags)
                            new_roots.append(sl)
                
                    # setup single point calcs of molecules (at bias) above surface (at bias)
                    # is single point a good estimate?
                    if self.args.add_desorbed == 'True' and 'desorb' in mv:
                        desorb_biases = mv['desorb']
                        if not os.path.exists(os.path.join(calc_folder, 'desorbed', surf)):
                            os.mkdir(os.path.join(calc_folder, 'desorbed', surf))
                        if not os.path.exists(os.path.join(calc_folder, 'desorbed', surf, mol)):
                            os.mkdir(os.path.join(calc_folder, 'desorbed', surf, mol))
                        for bias in desorb_biases:
                            surf_root = os.path.join(calc_folder, 'surfs', surf, self.get_bias_str(bias))
                            mol_root = os.path.join(calc_folder, 'molecules', mol, self.get_bias_str(bias))
                            if surf_root not in converged or mol_root not in converged:
#                                print('Desorbed waiting on converged surface and molecule.')
                                continue
                            # surf and mol at same bias exists
                            st_surf = Structure.from_file(os.path.join(surf_root, 'CONTCAR'))
                            ads_dic = {mol: ['center']}
                            head_folder = os.path.join(calc_folder, 'desorbed', surf, mol, self.get_bias_str(bias))
                            surf_ads = add_adsorbates(st_surf, ads_dic, ads_distance = des_distance,
                                                      molecules_loc = os.path.join(mol_root,'CONTCAR'))
                            save_locs = save_structures(surf_ads[mol], head_folder, skip_existing = True,
                                                        single_loc=True)
                            if len(save_locs) == 0:
                                continue
                            for sl in save_locs:
                                print('Added desorbed folder: '+sl)
                                self.run('cp '+os.path.join(inputs_folder, 'desorbed_inputs')
                                         +' '+os.path.join(sl, 'inputs'))
                                # tags is a list
                                tags = mv['tags'] if 'tags' in mv else []
                                if 'tags' in v: tags += v['tags']
                                tags += ['target-mu '+ ('None' if bias in ['None','none','No_bias'] 
                                         else '%.4f' % self.get_mu(bias, self.read_inputs(sl), tags))]
                                if desorbed_single_point:
                                    tags += ['max_steps 0']
                                self.add_tags(sl, tags)
                                new_roots.append(sl)
                        
                    # Create NEB calcs from converged ads. and des. calcs
                    if self.args.make_neb == 'True' and 'NEB' in mv:
                        # managed[surf][mol]['NEB']=[
                        #       {'site': site, 'bias': bias, 'images': images, 'fmax': fmax} ]
                        for neb in mv['NEB']:
                            # add NEB folders: calcs/neb/surf/mol/bias/number/
                            if not os.path.exists(os.path.join(calc_folder, 'neb', surf)):
                                os.mkdir(os.path.join(calc_folder, 'neb', surf))
                            if not os.path.exists(os.path.join(calc_folder, 'neb', surf, mol)):
                                os.mkdir(os.path.join(calc_folder, 'neb', surf, mol))
                            # get and check dependent dirs
                            bias = self.get_bias_str(neb['bias'])
                            ads_num = neb['site']#.zfill(2)
                            ads_root = os.path.join(calc_folder, 'adsorbed', surf, mol, bias, ads_num)
                            des_root = os.path.join(calc_folder, 'desorbed', surf, mol, bias)
                            if ads_root not in converged or des_root not in converged:
                                continue
                            # both ads and des calcs are converged at same bias, setup NEB folder
                            if not os.path.exists(os.path.join(calc_folder, 'neb', surf, mol, bias)):
                                os.mkdir(os.path.join(calc_folder, 'neb', surf, mol, bias))
                            neb_folder = os.path.join(calc_folder, 'neb', surf, mol, bias, ads_num)
                            if not os.path.exists(neb_folder):
                                os.mkdir(neb_folder)
                            if os.path.exists(os.path.join(neb_folder, 'inputs')):
                                # calculation is already set up
                                continue
                            nimages = neb['images']
                            
                            neb_set = setup_neb(os.path.join(ads_root, 'CONTCAR'), 
                                                os.path.join(des_root, 'CONTCAR'),
                                                nimages, neb_folder)
                            if not neb_set:
                                print('METAERROR: Failed to correctly setup NEB directory.')
                                continue
                            
#                            # setup neb calc with Neb_Make.py
#                            neb_make_cmd = ('Neb_Make.py '+ads_root+' '+des_root+' -i '+str(nimages)
#                                            +' -d '+neb_folder)
#                            if self.args.neb_climbing == 'True':
#                                neb_make_cmd += ' -c'
#                            self.run(neb_make_cmd)
                            # TODO: add check that folders and files were correctly setup
                            # all image sub_dirs should be setup! Still need to add inputs.
                            self.run('cp '+os.path.join(inputs_folder, 'neb_inputs')
                                     +' '+os.path.join(neb_folder, 'inputs'))
                            # tags is a list
                            tags = mv['tags'] if 'tags' in mv else []
                            if 'tags' in v: tags += v['tags']
                            tags += ['target-mu '+ ('None' if bias in ['None','none','No_bias'] 
                                     else '%.4f' % self.get_mu(neb['bias'], self.read_inputs(neb_folder), tags))]
                            tags += ['nimages '+str(nimages)]
                            self.add_tags(neb_folder, tags)
                            new_roots.append(neb_folder)
                            print('Setup NEB folder: '+neb_folder)
#                        print('NEB management not yet configured. Please contact Nick to add.')
        return new_roots
    
    def get_mol_loc(self, mol):
        if os.path.exists(os.path.join(molecule_folder, mol, 'POSCAR')):
            return os.path.join(molecule_folder, mol, 'POSCAR')
        elif os.path.exists(os.path.join(molecule_folder, mol, 'a', 'POSCAR')):
            return os.path.join(molecule_folder, mol, 'a', 'POSCAR')
        else:
            print('ERROR: Molecule not found: '+mol+'. Please add new folder with POSCAR to '+molecule_folder)
            return False
    
    def get_ref_mols(self, mol):
        refs = []
        for ref_mol in reference_molecules[mol]['refs']:
            if reference_molecules[ref_mol]['refs'] == [ref_mol]:
                refs += [ref_mol]
            else:
                refs += self.get_ref_mols(ref_mol)
        return refs
#        return reference_molecules[mol]['refs']
        
        
    def add_tags(self, root, tags):
        # inputs is a dictionary, tags is a list of lines to add
        inputs = self.read_inputs(root)
        inputs['restart'] = 'False'
        assert type(tags) == list, 'METAERROR: type for variable "tags" must be a list.'
        for tag in tags:
            tag_k = tag.split(' ')[0]
            tag_v = ' '.join(tag.split(' ')[1:])
            if tag_v in ['None']:
                if tag_k in inputs:
                    del inputs[tag_k]
            elif tag_v in ['pH','ph']:
                inputs['pH'] = tag_v
            else:
                inputs[tag_k] = tag_v
        self.write_inputs(inputs, root)
    
    def upgrade_calc(self, new_root, old_root, bias, tags, verbose = True):
        os.mkdir(new_root)
        self.run('cp ' + os.path.join(old_root, 'CONTCAR') + ' ' + os.path.join(new_root, 'POSCAR'))
        inputs = self.read_inputs(old_root)
        if bias in ['No_bias']:
            if 'target-mu' in inputs:
                del inputs['target-mu']
        else:
            mu = self.get_mu(bias, inputs, tags)
            inputs['target-mu'] = '%.4f'%(mu)
        inputs['restart'] = 'False'
        self.write_inputs(inputs, new_root)
        if verbose: print('Upgraded '+old_root+' to '+new_root)
    
    def get_mu(self, bias, inputs, tags = []):
        if bias == 'None':
            return 'None'
        assert 'fluid' in inputs, 'ERROR: fluid tag must be in inputs files to run biases!'
        if inputs['fluid'] == 'LinearPCM' and 'pcm-variant' in inputs and inputs['pcm-variant'] == 'CANDLE':
            Vref = 4.66
        elif inputs['fluid'] == 'LinearPCM' and 'pcm-variant' in inputs and inputs['pcm-variant'] == 'GLSSA13':
            Vref = 4.68
        elif inputs['fluid'] == 'NonlinearPCM' and ('pcm-variant' in inputs 
                   and inputs['pcm-variant'] == 'GLSSA13'):
            Vref = 4.62
        elif inputs['fluid'] == 'SaLSA':
            Vref = 4.54
        elif inputs['fluid'] == 'ClassicalDFT':
            Vref = 4.44
        else:
            assert False, ('ERROR: Fluid model must be in [CANDLE, SaLSA, ClassicalDFT]. '+
                           'Other models not yet configured.')
        rhe_shift = 0
        if self.args.rhe_zeroed:
            #JDFT uses SHE as zero point. 0V vs RHE === (-0.0591 * pH) V vs SHE
#            rhe_shift = -0.0591 * self.args.ph_rhe # input is RHE/V bias, output is SHE/JDFT/Hartree bias
            pH = 7.0
            for tag in tags:
                if 'pH' in tag or 'ph' in tag: 
                    pH = float(tag.split()[-1])
            rhe_shift = -0.0591 * pH #7 if 'pH' not in inputs else -0.0591 * float(inputs['pH'])
        return -(Vref + bias + rhe_shift)/27.2114 
        
    def read_manager_control(self, mc_text):
        '''
        Reads manager_control.txt file so that it can be used to setup and upgrade calculations
        '''
        managed = {'molecules': {}}
        ignore = True
        surf = None
        surf_bias = None
        mol = None
        error = False

        for line in mc_text.split('\n'):
            if '----- CALCULATIONS BELOW -----' in line:
                ignore = False
                continue
            if ignore or line == '':
                continue
            if line[0] == '#':
                continue
            if line[0] == '=':
                # new surface 
                surf = line[1:]
                managed[surf] = {}
                mol = None
                surf_bias = None
                mol_bias = None
                continue
            if line[0] == '-':
                # molecule for surface
                mol = line.split(':')[0][1:]
                # DONE: fix so this can read other types of inputs, even double inputs
#                sites = [int(x) for x in self.read_bias(line)]
                sites = self.get_sites(line)
                if surf is None:
                    print('Error in manager_control: surface must be listed before molecule '+mol)
                    error = True
                    continue
                managed[surf][mol] = {'sites': sites}
                managed['molecules'][mol] = {'biases': []}
                continue
            if line[0] == '+':
                # tag for inputs file
                if surf is None:
                    print('Error in manager_control: surface must be listed before adding "+" tags')
                    error = True
                    continue
                if mol is None:
                    if 'tags' not in managed[surf]:
                        managed[surf]['tags'] = []
                    managed[surf]['tags'].append(line[1:])
                else:
                    if 'tags' not in managed[surf][mol]:
                        managed[surf][mol]['tags'] = []
                    managed[surf][mol]['tags'].append(line[1:])
                continue
            if 'Biases:' in line:
                if surf is None:
                    print('Error in manager_control: surface must be listed before biases using "=".')
                    error = True
                    continue
                if mol is None:
                    # biases for surface
                    surf_bias = self.read_bias(line)
                    managed[surf]['biases'] = surf_bias
                else:
                    # biases for molecule
                    mol_bias = self.read_bias(line)
                    if surf_bias is None:
                        print("Error in manager_control: surface "+surf+" has no biases listed before mol "+mol)
                        error = True
                        continue
                    if any([b not in surf_bias for b in mol_bias]):
                        print("Error in manager_control: biases for mol "+mol+" not all in surface "+surf+" biases")
                        error = True
                        continue
                    managed[surf][mol]['biases'] = mol_bias
                    managed['molecules'][mol]['biases'] += mol_bias
                continue
            if 'Desorb:' in line:
                if mol is None:
                    print("Error in manager_control: molecule must be listed before Desorb")
                    error = True
                    continue
                if mol_bias is None:
                    print("Error in manager_control: biases for mol "+mol+" must be listed before desorb")
                    error = True
                    continue
                desorb_bias = self.read_bias(line)
                if any([b not in mol_bias for b in desorb_bias]):
                    print("Error in manager_control: desorb biases for mol "+mol+" must be in mol biases")
                    error = True
                    continue
                managed[surf][mol]['desorb'] = desorb_bias
                continue
            if 'NEB:' in line:
                if surf is None or surf_bias is None or mol is None or mol_bias is None or desorb_bias is None:
                    print("Error in manager_control: NEB tag must be nested under a molecule "+
                          "with Biases and Desorb")
                    error = True
                    continue
                try:
                    site = line.split()[1]
                    bias = 'No_bias' if line.split()[2] in ['None','No_bias','none'] else float(line.split()[2])
                    images = int(line.split()[3])
                    fmax = float(line.split()[4])
                except:
                    print("Error in manager_control: NEB line for mol "+mol+" entered incorrectly: "+line)
                    error = True
                    continue
                if 'NEB' not in managed[surf][mol]:
                    managed[surf][mol]['NEB'] = []
                managed[surf][mol]['NEB'].append({'site': site, 'bias': bias, 'images': images, 'fmax': fmax})
                continue
            
        if error:
            print('\nPlease fix errors in manager_control.txt and rerun\n')
            return False
        return managed
    
    def get_sites(self, line):
        sites = line.split('[')[-1].replace(']','').split(', ')
        for i,s in enumerate(sites):
            try:
                sites[i] = int(s)
            except:
                pass
        return sites
        
    def check_surface(self, file, dist = 8):
        # check that surface seems reasonable for manager to handle
        st = Structure.from_file(file)
        # check that z-direction is long distance
        if not (st.lattice.c > st.lattice.a and st.lattice.c > st.lattice.b):
            print('Bad Surface: z-direction is not longest side, skipping: '+file)
            return False
        # check that 10A (dist) of space is available in the vacuum above surface
        coords = st.cart_coords
        max_lens = [np.max(coords[:,0]), np.max(coords[:,1]), np.max(coords[:,2])]
#        if max_lens[2] <= max_lens[0] or max_lens[2] <= max_lens[1]:
#            print('Bad Surface: surface is wider than it is tall, skipping: '+file)
#            return False
        if st.lattice.c - max_lens[2] < dist:
            print('Bad Surface: Needs at least '+str(dist)+'A of vacuum space above surface, skipping: '+file)
            return False
        return True
    
    def add_calc_inputs(self, folders):
        for root in folders:
            # get type of calculation
            calc_type = None
            for subf in self.calc_subfolders:
                tag = os.path.join(calc_folder, subf)
                if tag in root:
                    calc_type = subf
                    continue
            if calc_type is None:
                print('Error: No calc_type found for root: '+root+' Skipping.')
                continue
            inputs_file = calc_type + '_inputs'
            self.run('cp '+os.path.join(inputs_folder, inputs_file)+' '+os.path.join(root, 'inputs'))
    
    def get_running_jobs_dirs(self):
        p = subprocess.Popen(['squeue' ,'-o', '"%Z %T"'],   #can be user specific, add -u username 
                         stdout=subprocess.PIPE)
        jobs_running = []
        for i,line in enumerate(p.stdout):
            if i == 0:
                continue
            jobs_running.append(str(line, 'utf-8').replace('"', '').split()[0])
        return jobs_running
        
        
    def read_bias(self, line):
        bias_str = line.split('[')[-1].replace(']','').split(', ')
        return ['No_bias' if x in ['None','none','No_bias'] else float(x) for x in bias_str]
    
    def update_run_new(self, run_new):
        for root in run_new:
            inputs = self.read_inputs(root)
            inputs['restart'] = 'False'
            self.write_inputs(inputs, root)
        
    def run_new_calcs(self, new_calcs):
        print('\n----- Running new calcs -----\n')
        for root in new_calcs:
            os.chdir(root)
            self.run(self.run_cmd + ' -o '+self.get_job_name(root))
            os.chdir(self.cwd)
            print('Calculation run: '+self.get_job_name(root))
    
    def get_parallel_running(self): # TODO: finish this so that -p jobs don't run over eachother
        shell_folder = 'tmp_parallel'
        if not os.path.exists(shell_folder):
            return []
        if 'running.txt' not in os.listdir(shell_folder):
            return []
    
    def run_all_parallel(self, roots):
        '''
        This function handles submitting clusters of calculations together on single nodes
        '''
        print('\n----- Running All Calcs in Parallel -----\n')
        
        max_per_node = self.args.parallel
        total_calcs = len(roots) 
        if total_calcs < max_per_node:
            max_per_node = total_calcs
        cores_per_node = core_architecture
        total_nodes = int(np.ceil(total_calcs / max_per_node))
        cores_per_calc = int(np.floor(cores_per_node / max_per_node))
        
        shells = []
        shell_folder = 'tmp_parallel'
        if os.path.exists(shell_folder):
            for shell_root, dirs, files in os.walk(shell_folder):
                for name in files:
                    os.remove(os.path.join(shell_root, name))
        else:
            os.mkdir(shell_folder)
        
        # write calcs for each node to tmp_parallel folder
        for i in range(total_nodes):
            # TODO: set cores per calc to update based on number of calcs in sub_roots
            sub_roots = roots[i*max_per_node:(i+1)*max_per_node]
            out_file = 'submit_'+str(i)
            write_parallel(sub_roots, self.cwd, cores_per_node, cores_per_calc, self.args.run_time, out_file, 
                           shell_folder, self.args.qos)
            shells.append(out_file + '.sh')
        
        # submit shell scripts
        os.chdir(shell_folder)
        for shell in shells:
            os.system('sbatch '+shell)
        os.chdir(self.cwd)
    
    def backup_calcs(self):
        if not os.path.exists(backup_folder):
            os.mkdir(backup_folder)
        # scan all sub_dirs in calcs and copy to backup_folder   
        for root, folders, files in os.walk(calc_folder):
            # copy folder over
            for i, sub in enumerate(root.split(os.sep)):
                if i == 0: continue
                sub_folder = os.sep.join(root.split(os.sep)[1:i+1])
                backup_f = os.path.join(backup_folder, sub_folder)
                if not os.path.exists(backup_f):
                    os.mkdir(backup_f)
            # copy over files
            for file in files_to_backup:
                if file in files:
                    self.run('cp '+os.path.join(root, file)+' '+os.path.join(backup_f, file))
        print('\nCalculation files backed up successfully.')
        
    def manager(self):
        '''
        Main class function for jdft_manager. Linearly runs the main sub-functions of jdft_manager 
        based on user inputs from command line. Main sub-functions also have descriptions.
        
        '''
        # ensure subfolders are correctly setup
        if self.args.setup == 'True':
            return
        assert self.proper_setup, ('ERROR: jdft_manager not yet setup! Please run with '+
                                   '-h to check input parameters or -s to setup folder.')
        
        # read through current results data 
        all_data = {}
        if os.path.isfile(self.data_file) and self.args.read_all != 'True':
            with open(self.data_file, 'r') as f:
                all_data = json.load(f)
                
        # get parallel tag
        parallel = self.args.parallel
        
        # scan through all subfolders to check for converged structures 
        add_inputs = []
        if self.args.check_calcs == 'True':
            # scan through folders
            running_jobs_dirs = self.get_running_jobs_dirs()
            all_data, add_inputs, rerun, run_new  = self.scan_calcs(all_data, running_jobs_dirs)
            # save all data 
            if self.args.save == 'True':
                # run analysis of converged calcs
                if self.args.analyze == 'True':
                    print('\n----- Running Calculation Analysis -----')
                    if len(all_data.keys()) < 2:
                        print('No data available yet!')
                    else:
                        all_data = self.analyze_data(all_data, reference_molecules)
                # save 
#                if len(all_data.keys()) < 2:
                print('\nSaving converged calculations.')
                with open(self.data_file, 'w') as f:
                    json.dump(all_data, f)
        
            # save and rerun unconverged (if requested)
            with open(os.path.join(results_folder, 'unconverged.txt'), 'w') as f:
                f.write('\n'.join(rerun))
            if self.args.rerun_unconverged == 'True' and len(rerun) > 0: #self.args.check_calcs == 'True' and 
#                print('\nRerunning unconverged calculations')
                self.update_rerun(rerun)
                if parallel == 1:
                    self.rerun_calcs(rerun)
        
        # make new surfaces based on manager_control.txt file, add new calcs to add_inputs
        if self.args.make_new == 'True':
            new_folders = self.make_new_calcs(all_data['converged'])
            if new_folders == False:
                print('Exiting.\n\n')
                return
            
            # add inputs to any created file with a POSCAR and no 'inputs'. 
            # Does not need to be in manager_control.txt
            if len(add_inputs) > 0:
                self.add_calc_inputs(add_inputs)
        
        # run jobs with added 'inputs'
        # molecules should only be run if requested or in manager_control
        # neb calcs should be handled differently ? submission is the same, inputs is all that's changed.
        #   sub_folders should start with single point calcs after setup! 
        # bias calcs should pass structure forward from nomu -> 0V -> other biases
        if self.args.run_new == 'True' and len(new_folders + add_inputs + run_new) > 0:
            if len(run_new) > 0:
                self.update_run_new(run_new)
            if parallel == 1:
                self.run_new_calcs(new_folders + add_inputs + run_new)
        
        all_roots = []
        if self.args.run_new == 'True':
            all_roots += new_folders + add_inputs + run_new
        if self.args.rerun_unconverged == 'True':
            all_roots += rerun
        if parallel > 1:
            self.run_all_parallel(all_roots)
            
        if self.args.backup == 'True':
            self.backup_calcs()
            
        print('----- Done -----\n\n')


# Eagle defaults
calc_folder = 'calcs/'
molecule_folder = 'molecules/'
results_folder = 'results/'
inputs_folder = 'inputs/'
backup_folder = 'backup/'

try:
    home_dir = os.environ['JDFTx_home']
except:
    home_dir = '/home/nicksingstock'

try:
    core_architecture = int(os.environ['CORES_PER_NODE'])
except:
    core_architecture = 36
defaults_folder = os.path.join(home_dir, 'bin/JDFTx_Tools/manager/defaults/')
run_command = 'python '+ os.path.join(home_dir, 'bin/JDFTx_Tools/sub_JDFTx.py')

files_to_backup = ['POSCAR','CONTCAR','out','inputs','opt.log','Ecomponents']

'''
TODOs:
DONE    1. make sure everything runs so far
DONE    2. add/create default inputs files and molecules 
DONE        2.1. add template file for adding surfs, sites and mol binding to manage
DONE    3. add scan_calcs funtion and data saving
DONE        3.1. add warnings for divergent calcs
DONE        3.2. add ability to check whether calc is already running in scan_calcs, do NOT rerun!
DONE        3.3. add surface checking
DONE    4. add compatibility with 'surface_adsorbate_setup_CP_HER_v2' and ability to 
DONE       4.1. create desorbed calcs (single point?)
DONE    5. add structure creation
DONE    6. add molecule running
DONE    7. add dos options to inputs file
    8. make interpreter that analyzes all output data
DONE    9. setup NEB calcs and manage
    10. Add calculation_info.json (or.txt) to each folder with useful info about
        the calculation setup and initial structure. Update each run. 
    11. setup function that checks whether managed has changed that may cause errors
    12. setup function to check that different inputs files are compatible
DONE    13. Add ability to combine multiple calcs into one node for job submission efficiency
*    14. Run DOS SP calcs on converged structures (if user requested)
    15. Make ref_mols a separate file so user can edit
    16. Optional pre-NEB image SP calcs for wfns.
*    17. Add NEB building from any two directories (with warnings/errors about structure diffs)
        to enable non-adsorption barriers to be studied
    18. Add back in ability to run two+ adsorbates as a single calculation with "_" separator
*    19. Add charge density analysis using Aziz script on Summit. Oxi-states available in out file.
    
DONE   FIX ERROR: -r jobs are running with resart False so they keep starting from POSCAR

surf -> run calc nomu -> run calc 0 V -> run other biases 
                      -> add adsorbates above converged surf at same bias
                      -> add desorbed state single point calcs at same bias with mol at same bias
mols -> converge if requested, runs all biases from same starting POSCAR. Use these for binding energy analysis.
neb -> setup from converged adsorbate+desorbed state at same bias
'''


if __name__ == '__main__':
    jm = jdft_manager()
    jm.manager()