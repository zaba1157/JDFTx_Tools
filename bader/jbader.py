#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:46:13 2020

@author:Jacob and edited by Yousef
"""


import os, sys
import argparse
import shutil
from shutil import copy2

import gen_bgw


def readacf(acffile,Nat):
    #get atom electrons from ACF.dat
    tempfile = gen_bgw.readfile(acffile)
    tempfile = [x.split() for x in tempfile]
    chg = [float(tempfile[x][4]) for x in range(2,2+Nat)]
    line = 2+Nat+1
    vacchg = float(tempfile[line][2])
    vacvol = float(tempfile[line+1][2])
    Nelec = float(tempfile[line+2][3])
    return chg,[vacchg,vacvol,Nelec]

def parse_args():
    #get arguments from command line
    parser = argparse.ArgumentParser(description="Script to do Bader analysis from JDFTx calculation")
    parser.add_argument("-f", type=str,default='out', help="output file to read (default: out)")
    parser.add_argument('-o', type=str,default='CHGCAR', help='name of CHGCAR file to write (default: CHGCAR)')
    parser.add_argument('-p', type=str,default='jdft', help='prefix of .n files (default: jdft)')
    parser.add_argument('-i', type=int,default=2,help='set Ninterp in createVASP (default: 2)')
    args = parser.parse_args()
    return args

def main():

    #store arguments
    args = parse_args()
    outfile = args.f
    chgcar = args.o
    prefix = args.p
    Ninterp = args.i

    tempfile = gen_bgw.readfile(outfile)
    #get number atoms
    line = gen_bgw.find_key('total atoms',tempfile)
    Nat = int(float(tempfile[line].split()[4]))
    #determine number of spins so whether need to read .n or .n_up and .n_dn files
    line = gen_bgw.find_key('spintype',tempfile)
    spintype = tempfile[line].split()[1]
    #get neutral atom valence electrons
    endline = gen_bgw.find_key('Setting up symmetries',tempfile)
    ionslines = gen_bgw.find_first_range_key('ion ',tempfile,startline=0,endline=endline)
    ions = [tempfile[x].split()[1] for x in ionslines]
    atomlines = gen_bgw.find_all_key(' pseudopotential, ',tempfile)
    atoms = [tempfile[x].split()[0].strip('\'') for x in atomlines]
    vallines = gen_bgw.find_all_key('valence electrons',tempfile)
    valelec = [float(tempfile[x].split()[0]) for x in vallines]
    atomtypes = dict(zip(atoms,valelec))

    if spintype == 'z-spin' or spintype == 'vector-spin':
        Nspin = 2
    elif spintype == 'no-spin' or spintype == 'spin-orbit':
        Nspin = 1
    else:
        print('Could not determine spin, assuming no spin')
        Nspin = 1

    if Nspin == 1:
        #make VASP CHGCAR
        print('Making CHGCAR')
        com = ['createVASP',outfile,chgcar,prefix + '.n',str(Ninterp)]
        gen_bgw.run_command(com)
        #run bader
        print('Run Bader on CHGCAR')
        com = ['bader',chgcar]
        gen_bgw.run_command(com)

    else:
        #do same process for both spins
        print('Making CHGCAR_up')
        com = ['createVASP',outfile,chgcar+'_up',prefix + '.n_up',str(Ninterp)]
        gen_bgw.run_command(com)
        print('Making CHGCAR_dn')
        com = ['createVASP',outfile,chgcar+'_dn',prefix + '.n_dn',str(Ninterp)]
        gen_bgw.run_command(com)
        print('Run Bader on CHGCAR_up')
        com = ['bader',chgcar+'_up']
        gen_bgw.run_command(com)
        shutil.move('ACF.dat','ACF.dat_up')
        print('Run Bader on CHGCAR_dn')
        com = ['bader',chgcar+'_dn']
        gen_bgw.run_command(com)
        shutil.move('ACF.dat','ACF.dat_dn')

        atomchgsup,totalchgsup = readacf('ACF.dat_up',Nat)
        atomchgsdn,totalchgsdn = readacf('ACF.dat_dn',Nat)
        totalchg = [atomchgsup[i]+atomchgsdn[i] for i in range(len(atomchgsup))]

        #write summed ACF.dat file
        print('Making final summed ACF.dat file')
        copy2('ACF.dat_up','ACF.dat')
        tempfile = gen_bgw.readfile('ACF.dat')
        tempfile[0] = '    #         X           Y           Z        CHARGE     OXIDATION STATE\n'
        for i in range(Nat):
            line = tempfile[i+2].split()[:4]
            line.append(str(totalchg[i]))
            line.append(str(atomtypes[ions[i]]-totalchg[i]))
            line = [float(x) for x in line]
            line[0] = int(line[0])
            tempfile[i+2] = '{:5d} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:15.4f}\n'.format(line[0],line[1],line[2],line[3],line[4],line[5])
        i = 2+Nat+1
        tempfile[i] = '    VACUUM CHARGE: {:20.4f}\n'.format(totalchgsup[0]+totalchgsdn[0])
        tempfile[i+1] = '    VACUUM VOLUME: {:20.4f}\n'.format(totalchgsup[1]+totalchgsdn[1])
        tempfile[i+2] = '    NUMBER OF ELECTRONS: {:14.4f}\n'.format(totalchgsup[2]+totalchgsdn[2])
        gen_bgw.writefile('ACF.dat',tempfile)


if __name__ == '__main__':
    main()


