#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Musgrave Group

@author: zaba1157, nisi6161
"""

from JDFTx import JDFTx
import os
from ase.io import read
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, GPMin, MDMin, FIRE
from ase.io.trajectory import Trajectory
from ase.neb import NEB
import argparse


opj = os.path.join
ope = os.path.exists

def force_checker(max_force = 500):
    if 'opt.log' in os.listdir():
        try:
            with open('opt.log', 'r') as f:
                opt_text = f.read()
        except:
            print('Unable to read opt.log for force checker.')
    elif 'neb.log' in os.listdir():
        try:
            with open('neb.log', 'r') as f:
                opt_text = f.read()
        except:
            print('Unable to read neb.log for force checker.')
    else:
        return
    if len(opt_text) == 0: return 
    
    opt = [line.split() for line in opt_text.split('\n') if line != '' and '*Force-consistent' not in line]
    force = float(opt[-1][4])
    assert force < max_force, 'ERROR: Calculation ended due to high forces. Edit Structure.'

def insert_el(filename):
    """
    Inserts elements line in correct position for Vasp 5? Good for
    nebmovie.pl script in VTST-tools package
    Args:
        filename: name of file to add elements line


    """
#    fp = open(filename,"r")
#    contents = fp.readlines()
#    line = 0
#    for i in contents:
#        if line == 0:
#            ele_line = i
#        else:
#            break
#        line += 1
#    fp.close()
    with open(filename, 'r') as f:
        file = f.read()
    contents = file.split('\n')
    ele_line = contents[0]
    if contents[5].split() != ele_line.split():
        contents.insert(5, ele_line)
    with open(filename, 'w') as f:
        f.write('\n'.join(contents))
#    contents.insert(5, ele_line)
#    fp = open(filename, "w")
#    fp.writelines(contents)
#    fp.close()

def add_dos(cmds):
    from pymatgen.core.structure import Structure
    if not ope('./inputs_dos'):
        return cmds
    new_cmds = []
    for cmd in cmds:
        if 'density-of-states' in cmd:
            print('WARNING: command density-of-states in inputs file is being overwritten by inputs_dos!')
            continue
        new_cmds += [cmd]
    with open('./inputs_dos','r') as f:
        dos = f.read()
    st = Structure.from_file('./POSCAR')
    
    dos_line = ''
    for line in dos.split('\n'):
        if 'Orbital' in line and len(line.split()) == 3:
            # Format: Orbital Atom_type orbital_type(s)
            # Adds DOS for ALL atoms of this type and all requested orbitals
            atom_type = line.split()[1]
            indices = [i for i,el in enumerate(st.species) if str(el) == atom_type]
            orbitals = [orb for orb in line.split()[2]]
            assert all([orb in ['s','p','d','f'] for orb in orbitals]), 'ERROR: Not all orbital types allowed! ('+','.join(orbitals)+')'
            for i in range(len(indices)):
                for orbital in orbitals:
                    dos_line += ' \\\nOrbital ' +atom_type + ' ' + str(i+1) + ' ' + orbital 
        else:
            dos_line += ' \\\n' + line
    new_cmds += [('density-of-states', dos_line)]
    return new_cmds

def initialize_calc(command_file, jdftx_exe):

    notinclude = ['ion-species',
                  #'include',
                  #'lcao-params',
                  'ionic-minimize',
                  #'dump-name',
                  #'converge-empty-states',
                  'latt-scale',
                  'latt-move-scale',
                  'coulomb-interaction',
                  'coords-type',
                  'ion',
                  'climbing','pH','ph',
                  'logfile','pseudos','nimages','max_steps','fmax','optimizer','restart','parallel']

    def read_commands(command_file,notinclude):
        cmds = []
        script_cmds = {}
        with open(command_file,'r') as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0 or line[0] == '#': continue
                linelist = line.split()
                if linelist[0] in notinclude:
                    if len(linelist) > 1:
                        script_cmds[linelist[0]] = ' '.join(linelist[1:])
                else:
                    if len(linelist) > 1:
                        tpl = (linelist[0], ' '.join(linelist[1:]))
                        cmds.append(tpl)

        """ Defaults """

        #cmds['core-overlap-check']  = 'none'
        cmds += [('core-overlap-check', 'none')]

        return cmds, script_cmds



    cmds, script_cmds = read_commands(command_file,notinclude)
    cmds = add_dos(cmds)

    def calc_type(script_cmds):
        if 'nimages' in script_cmds.keys():
            calc = 'neb'
        else:
            calc = 'opt'

        return calc




    ctype = calc_type(script_cmds)
    psuedos = script_cmds['pseudos']
    max_steps = int(script_cmds['max_steps'])
    fmax = float(script_cmds['fmax'])
    restart = True if ('restart' in script_cmds and script_cmds['restart'] == 'True') else False
    parallel_neb = True if ('parallel' in script_cmds and script_cmds['parallel'] == 'True') else False
    climbing_neb = True if ('climbing' in script_cmds and script_cmds['climbing'] == 'True') else False

    if ctype == 'opt':
        if restart:
            atoms = read('CONTCAR',format='vasp')
        else:
            atoms = read('POSCAR',format='vasp')
        jdftx_num_procs = os.environ['JDFTx_NUM_PROCS']
        exe_cmd = 'mpirun -np '+str(jdftx_num_procs)+' '+jdftx_exe

        #Set up JDFTx calculator
        calculator = JDFTx(
            executable = exe_cmd,
            pseudoSet=psuedos,
            commands=cmds,
            outfile = os.getcwd()
            )

        #Set calculator
        atoms.set_calculator(calculator)


        def optimizer(opt='BFGS',imag_atoms=atoms,logfile='opt.log'):

            """
            ASE Optimizers:
                BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, GPMin, MDMin and FIRE.
            """

            opt_dict = {'BFGS':BFGS, 'BFGSLineSearch':BFGSLineSearch,
                        'LBFGS':LBFGS, 'LBFGSLineSearch':LBFGSLineSearch,
                        'GPMin':GPMin, 'MDMin':MDMin, 'FIRE':FIRE}

            dyn = opt_dict[opt](imag_atoms,logfile=logfile,restart='hessian.pckl')
            return dyn


        #Structure optimization
        #dyn = BFGS(atoms,logfile='opt.log')
        if 'optimizer' in script_cmds.keys():
            dyn = optimizer(opt=script_cmds['optimizer'], imag_atoms=atoms)
        else:
            dyn = optimizer(imag_atoms=atoms)



        def write_contcar(a=atoms):
            a.write('CONTCAR',format="vasp", direct=True)
            insert_el('CONTCAR')


        '''
        Only energy and forces seem to be implemented in JDFTx.py. A Trajectory
        object must be attached so JDFTx does not error out trying to get stress
        '''

        traj = Trajectory('opt.traj', 'w', atoms, properties=['energy', 'forces'])
        dyn.attach(traj.write, interval=1)
        dyn.attach(write_contcar,interval=1)
        dyn.attach(force_checker,interval=1)

        dyn.run(fmax=fmax,steps=max_steps)

        traj.close()

    elif ctype == 'neb':

        #jdftx_num_procs = 36 #1
        jdftx_num_procs = os.environ['JDFTx_NUM_PROCS']
        exe_cmd = 'mpirun -np '+str(jdftx_num_procs)+' '+jdftx_exe

        nimages = int(script_cmds['nimages'])
        image_dirs = [str(i).zfill(2) for i in range(0,nimages+2)]
        
        try:
            initial = read('00/POSCAR')
            final = read(opj(str(nimages+1).zfill(2),'POSCAR'))
        except:
            initial = read('00/CONTCAR')
            final = read(opj(str(nimages+1).zfill(2),'CONTCAR'))
        
        images = [initial]
        if restart:
            try:
                images += [read(opj(i,'CONTCAR')) for i in image_dirs[1:-1]]
            except:
                images += [read(opj(i,'POSCAR')) for i in image_dirs[1:-1]]
                print('WARNING: CONTCAR files not found, starting from POSCARs')
        else:
            images += [read(opj(i,'POSCAR')) for i in image_dirs[1:-1]]
        images += [final]
        #images = [read(opj(i,'POSCAR')) for i in image_dirs ]

        neb = NEB(images, parallel=parallel_neb, climb=climbing_neb) 
        for i, image in enumerate(images[1:-1]):
            #if i == j:
            image.calc = JDFTx(
                            executable = exe_cmd,
                            pseudoSet=psuedos,
                            commands=cmds,
                            outfile = opj(os.getcwd(),image_dirs[i+1])
                            )

        def optimizer(opt='BFGS',imag_atoms=images,logfile='neb.log'):

            """
            ASE Optimizers:
                BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, GPMin, MDMin and FIRE.
            """

            opt_dict = {'BFGS':BFGS, 'BFGSLineSearch':BFGSLineSearch,
                        'LBFGS':LBFGS, 'LBFGSLineSearch':LBFGSLineSearch,
                        'GPMin':GPMin, 'MDMin':MDMin, 'FIRE':FIRE}

            dyn = opt_dict[opt](imag_atoms,logfile=logfile)
            return dyn



        if 'optimizer' in script_cmds.keys():
            dyn = optimizer(opt=script_cmds['optimizer'], imag_atoms=neb)
        else:
            dyn = optimizer(imag_atoms=neb)


        '''
        Only energy and forces seem to be implemented in JDFTx.py. A Trajectory
        object must be attached so JDFTx does not error out trying to get stress
        '''

        traj = Trajectory('neb.traj', 'w', neb, properties=['energy', 'forces'])


        def write_contcar(img_dir, image):
            image.write(opj(img_dir,'CONTCAR'),format="vasp", direct=True)
            insert_el(opj(img_dir,'CONTCAR'))



        for i,image in enumerate(images[1:-1]):
            img_dir = image_dirs[i+1]
            dyn.attach(Trajectory(opj(os.getcwd(),img_dir,'opt-'+img_dir+'.traj'), 'w', image,
                                  properties=['energy', 'forces']))

            dyn.attach(write_contcar, interval=1, img_dir=img_dir, image=image)

        dyn.attach(force_checker,interval=1)

        dyn.run(fmax=fmax,steps=max_steps)
        traj.close()



if __name__ == '__main__':

    jdftx_exe = os.environ['JDFTx']
    
    # optional, change to another directory (for parallel runs)
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='Directory to run in and save files to.',
                        type=str, default='./')
    parser.add_argument('-g', '--gpu', help='If True, runs GPU install of JDFTx.',
                        type=str, default='False')


    args = parser.parse_args()
    if args.dir != './':
        os.chdir(args.dir)
    if args.gpu == 'True':
        try:
            jdftx_exe = os.environ['JDFTx_GPU']
        except:
            print('Environment variable "JDFTx_GPU" not found, running standard JDFTx.')

    command_file = 'inputs'
    initialize_calc(command_file, jdftx_exe)
