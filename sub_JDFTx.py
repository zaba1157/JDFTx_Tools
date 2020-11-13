#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 20:43:05 2020

@author: zaba1157
"""
import argparse
import os


def write(nodes,cores,time,out,alloc,script):
    np=nodes*cores
    writelines = '#!/bin/bash'+'\n'
    writelines+='#SBATCH -J '+out+'\n'
    writelines+='#SBATCH --time='+str(time)+':00:00'+'\n'
    writelines+='#SBATCH -N '+str(nodes)+'\n'
    writelines+='#SBATCH --tasks '+str(np)+'\n'
    writelines+='#SBATCH -o '+out+'-%j.out'+'\n'
    writelines+='#SBATCH -e '+out+'-%j.err'+'\n'
    writelines+='#SBATCH --account='+alloc+'\n'
    if time == 1:
        writelines+='#SBATCH --partition=debug\n'
    writelines+='export JDFTx_NUM_PROCS='+str(np)+'\n'
    writelines+='python '+script+' > '+out+'\n'
    writelines+='exit 0'+'\n'

    with open('submit.sh','w') as f:
        f.write(writelines)


if __name__ == '__main__':
    
    script = '/home/zaba1157/jdftx_scripts/run_JDFTx.py'

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nodes', help='Number of nodes',
                        type=int, default=1)
    parser.add_argument('-c', '--cores', help='Number of cores',
                        type=int, default=36)
    parser.add_argument('-t', '--time', help='Time limit',
                        type=int, default=48)
    parser.add_argument('-o', '--outfile', help='Outfile name',
                        type=str, required=True)
    parser.add_argument('-a', '--allocation', help='Allocation',
                        type=str,default='custws')
    #parser.add_argument('-s', '--script', help='Python script to submit.',
    #                    type=str, required=True)
    
    args = parser.parse_args()

    write(args.nodes,args.cores,args.time,args.outfile,args.allocation,
          script)
    os.system('sbatch submit.sh')
    
