#!/usr/bin/env python3

#Created on Wed Oct  2 10:14:24 2019
#@author: Nick

from pymatgen.core.surface import SlabGenerator, Structure
import os
import argparse

def make_surface(file, index, slab_height, vac_space):
    st = Structure.from_file(file)

    slabgen = SlabGenerator(st, index, slab_height, vac_space, center_slab=True)
    all_slabs = slabgen.get_slabs(symmetrize = False)
    print("The slab has %s termination." %(len(all_slabs)))

    if not os.path.exists('surf_subfolder/'):
        os.mkdir('surf_subfolder/')

    for i, slab in enumerate(all_slabs):
        slab.to('POSCAR', 'surf_subfolder/POSCAR_'+str(i).zfill(2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='File to read',
                        type=str, default="POSCAR")
    parser.add_argument('-i', '--index', help='Surface Miller Index',
                        type=str, default="100")
    parser.add_argument('-sh', '--slab_height', help='Height of the slab',
                        type=int, default=15)
    parser.add_argument('-vh', '--vac_space', help='Height of vacuum space',
                        type=int, default=20)

    args = parser.parse_args()
	
    mindex = tuple([int(x) for x in args.index])
	
    make_surface(args.file, mindex, args.slab_height, args.vac_space)
