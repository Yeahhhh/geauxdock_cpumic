#!/usr/bin/env python

from __future__ import print_function

import os
import sys


def write_file(filename, str):
    #    print(filename, ':')
    #    print(str)
    with open(filename, 'w') as f:
        print(str, file=f)


def generate_lig_txt(sdf, ff):
    sdf_path = os.path.abspath(sdf)
    ff_path = os.path.abspath(ff)
    lig_id = os.path.basename(sdf_path).split('.')[0]
    str = lig_id + ' ' + sdf_path + " " + ff_path
    write_file('lig.txt', str)


def generate_prt_txt(myfile):
    mypath = os.path.abspath(myfile)
    myid = os.path.basename(mypath).split('.')[0]
    str = myid + ' ' + mypath
    write_file('prt.txt', str)


def run(bin, prt_file, lig_ff_file):
    cmd = bin
    cmd = cmd + ' -nt ' + "1"
    cmd = cmd + ' -lp ' + prt_file
    cmd = cmd + ' -ll ' + lig_ff_file
    print(cmd)
    os.system(cmd)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit('Usage: %s dock.exe prt.pdb lig.sdf ff.ff' % sys.argv[0])

    bin = sys.argv[1]
    prt = sys.argv[2]
    lig = sys.argv[3]
    ff = sys.argv[4]

    generate_prt_txt(prt)
    generate_lig_txt(lig, ff)
    run(bin, 'prt.txt', 'lig.txt')

# example command:
# ./rundock.py ./dock 1a07C-1.pdb 1a07C1.sdf 1a07C1-0.8.ff
