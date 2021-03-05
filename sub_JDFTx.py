#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 20:43:05 2020

@author: zaba1157
"""
import argparse
import os


def write(nodes,cores,time,out,alloc,qos,script):
    np=nodes*cores
    writelines = '#!/bin/bash'+'\n'
    writelines+='#SBATCH -J '+out+'\n'
    writelines+='#SBATCH --time='+str(time)+':00:00'+'\n'
    writelines+='#SBATCH -o '+out+'-%j.out'+'\n'
    writelines+='#SBATCH -e '+out+'-%j.err'+'\n'

    writelines+='#SBATCH --tasks '+str(np)+'\n'
    writelines+='#SBATCH --nodes '+str(nodes)+'\n'
    writelines+='#SBATCH --ntasks-per-node '+str(cores)+'\n'

    if alloc=='environ':
        writelines+='#SBATCH --account='+os.environ['JDFTx_allocation']+'\n'
    else:
        writelines+='#SBATCH --account='+alloc+'\n'
    if qos=='high':
        writelines+='#SBATCH --qos=high'+'\n'

    if time == 1:
        writelines+='#SBATCH --partition=debug\n'
    writelines+='\nexport JDFTx_NUM_PROCS='+str(np)+'\n'
    writelines+='module load comp-intel/2020.1.217 intel-mpi/2020.1.217 cuda/10.2.89 vasp/6.1.1 mkl/2020.1.217 gsl/2.5/gcc openmpi/4.0.4/gcc-8.4.0 gcc/7.4.0'+'\n\n'
    writelines+='python '+script+' > '+out+'\n'
    writelines+='exit 0'+'\n'

    with open('submit.sh','w') as f:
        f.write(writelines)


if __name__ == '__main__':
    
    script = '/home/nicksingstock/bin/jdft/run_JDFTx.py'

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
                        type=str,default='environ')
    parser.add_argument('-q', '--qos', help='Priority / QOS (e.g. high)',
                        type=str,default='standard')

    #parser.add_argument('-s', '--script', help='Python script to submit.',
    #                    type=str, required=True)
    
    args = parser.parse_args()

    write(args.nodes,args.cores,args.time,args.outfile,args.allocation,args.qos,
          script)
    os.system('sbatch submit.sh')
    
