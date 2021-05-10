#!/usr/bin/env python 

# print out 10x fofn with offset data

# Usage: python 10x_fofn.py 10x_dir 

import os
import sys
import glob

dirname = sys.argv[1]

assert os.path.isdir(dirname)

fns = glob.glob(dirname + '/' + '*R*gz')

for fn in fns:
    if 'R1' in fn:
        print(os.path.abspath(fn),23,sep='\t')
    elif 'R2' in fn:
        print(os.path.abspath(fn),0,sep='\t')