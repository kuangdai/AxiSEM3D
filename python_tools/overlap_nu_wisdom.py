#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
max_nu_wisdom.py

Overlap several Nu-wisdoms (wavefield learning results of Fourier expansion 
order fields in NetCDF format) by summation or maximum.

To see usage, type
python max_nu_wisdom.py -h
'''

################### PARSER ###################

aim = '''Overlap several Nu-wisdoms (wavefield learning results of Fourier expansion 
order fields in NetCDF format) by summation or maximum.'''

notes = '''Nu-wisdoms to be overlapped must be learnt with the same mesh 
(same background model and period).
  
'''

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=aim, epilog=notes, 
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-i', '--input', dest='in_wisdom_files', action='store', 
                    type=str, required=True, nargs='+',
                    help='Nu-wisdoms to be overlapped <required>') 
parser.add_argument('-o', '--output', dest='out_wisdom_file', action='store', 
                    type=str, required=True,
                    help='resultant Nu-wisdom <required>\n') 
parser.add_argument('-s', '--summation', dest='summation', action='store_true', 
                    help='overlap by summation; default = False (by maximum)\n')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset

# check size
assert len(args.in_wisdom_files) > 1, 'At least two input Nu-wisdoms are needed.'

# read
nus = []
for infile in args.in_wisdom_files:
    if args.verbose:
        print('Reading ' + infile)
    ncdf = Dataset(infile, 'r')
    data = ncdf.variables['axisem3d_wisdom'][:, :]
    ind = np.lexsort((data[:, 0], data[:, 1]))
    nus.append(data[ind])
    ncdf.close()
    if len(nus) > 1:
        assert data.shape == nus[0].shape, 'Incompatible mesh.'

# overlap
nu_res = nus[0]
for nu in nus[1:]:
    if args.summation:
        nu_res[:, 2] += nu[:, 2]
    else:
        nu_res[:, 2] = np.maximum(nu_res[:, 2], nu[:, 2])

# set starting field to -1 as an indication of overlap
nu_res[:, 3].fill(-1) 

# write
length = len(nu_res)
ncdf = Dataset(args.out_wisdom_file, 'w')
ncdf.createDimension('ncdim_4', size=4)
ncdf.createDimension('ncdim_' + str(length), size=length)
var = ncdf.createVariable('axisem3d_wisdom', float, 
    ('ncdim_' + str(length), 'ncdim_4'))
var[:, :] = nu_res[:, :]
ncdf.close()
if args.verbose:
    print('Result written to ' + args.out_wisdom_file)

