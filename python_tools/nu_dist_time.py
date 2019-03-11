#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
nu_dist_time.py

Compute Nu(distance, time) from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver). 

To see usage, type
python nu_dist_time.py -h
'''
    
################### PARSER ###################
aim = '''Compute Nu(distance, time) from a NetCDF database of surface wavefield
created by AxiSEM3D (named axisem3d_surface.nc by the solver).'''

notes = '''Parallelise data processing using --nporc option.
 
'''

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=aim, epilog=notes, 
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-i', '--input', dest='in_surface_nc', 
                    action='store', type=str, required=True,
                    help='NetCDF database of surface wavefield\n' + 
                         'created by AxiSEM3D <required>')
parser.add_argument('--multi_file', dest='multi_file', action='store_true', 
                    help='Does the NetCDF database consist of\n' +
                         'multiple files; default = False')                                                    
parser.add_argument('-o', '--output', dest='out_nu', 
                    action='store', type=str, required=True,
                    help='output NetCDF database <required>') 
parser.add_argument('-m', '--min_dist', dest='min_dist', 
                    action='store', type=float, default=0.,
                    help='minimum distance (deg); default = 0')
parser.add_argument('-M', '--max_dist', dest='max_dist', 
                    action='store', type=float, default=180.,
                    help='maximum distance (deg); default = 180')
parser.add_argument('-D', '--dist_interval', dest='dist_interval', 
                    action='store', type=float, default=1.,
                    help='distance interval (deg); default = 1')  
parser.add_argument('-t', '--tstart', dest='tstart', 
                    action='store', type=float, required=True,
                    help='start time (sec) <required>') 
parser.add_argument('-d', '--time_interval', dest='time_interval', 
                    action='store', type=float, required=True,
                    help='time interval (sec) <required>') 
parser.add_argument('-n', '--ntimesteps', dest='ntimesteps',
                    action='store', type=int, required=True,
                    help='number of timesteps <required>')
parser.add_argument('-c', '--component', dest='component', action='store', 
                    type=str, default='Z', choices=['R', 'T', 'Z'],
                    help='seismogram component; default = Z')
parser.add_argument('-e', '--error', dest='error', 
                    action='store', type=float, default=1e-3,
                    help='truncation error; default = 1e-3')                             
parser.add_argument('-p', '--nproc', dest='nproc', action='store', 
                    type=int, default=1, 
                    help='number of processors; default = 1')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')        
args = parser.parse_args()

################### PARSER ###################


import numpy as np
from netCDF4 import Dataset
import os, shutil
from multiprocessing import Pool
import time

###### read surface database
if args.verbose:
    clock0 = time.clock()
    print('Reading global parameters...')
    
if args.multi_file:
    # create first
    nc_surfs = []
    for irank in np.arange(0, 99999):
        fname = args.in_surface_nc + str(irank)
        if os.path.isfile(fname):
            nc = Dataset(fname, 'r')
            print(str(nc.variables['time_points'][0]))
            if nc.variables['time_points'][-1] < -1. or '--' in str(nc.variables['time_points'][0]):
                print('Skip opening nc file %s' % (fname))
                continue 
            nc_surfs.append(nc)
            if args.verbose:
                print('Done opening nc file %s' % (fname))
    nc_surf = nc_surfs[0]
else:
    nc_surf = Dataset(args.in_surface_nc, 'r')
    if args.verbose:
        print('Done opening nc file %s' % (args.in_surface_nc))
    nc_surfs = [nc_surf]

# global attribute
srclat = nc_surf.source_latitude
srclon = nc_surf.source_longitude
srcdep = nc_surf.source_depth
srcflat = nc_surf.source_flattening
surfflat = nc_surf.surface_flattening
r_outer = nc_surf.radius
# time
var_time = nc_surf.variables['time_points'][:]
nstep = len(var_time)
assert nstep > 0, 'Zero time steps'
t0 = var_time[0]
# theta
var_theta = nc_surf.variables['theta'][:]
nele = len(var_theta)
# GLL and GLJ
var_GLL = nc_surf.variables['GLL'][:]
var_GLJ = nc_surf.variables['GLJ'][:]
nPntEdge = len(var_GLL)

# element on rank
irank_ele = np.zeros(nele, dtype='int')
irank_ele.fill(-1)
for inc, nc in enumerate(nc_surfs):
    keys = nc.variables.keys()
    for eleTag in np.arange(nele):
        key = 'edge_' + str(eleTag) + 'r'
        if key in keys:
            irank_ele[eleTag] = inc


if args.verbose:
    elapsed = time.clock() - clock0
    print('Reading global parameters done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### surface sampling
if args.verbose:
    clock0 = time.clock()
    print('Sampling distance...')
ndists = int((args.max_dist - args.min_dist) / args.dist_interval) + 1
dists = np.radians(np.linspace(args.min_dist, args.max_dist, ndists))    
if args.verbose:
    elapsed = time.clock() - clock0
    print('    Number of distances: %d' % (ndists))
    print('Sampling distance done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### prepare theta
if args.verbose:
    clock0 = time.clock()
    print('Locating points in distance...')
# locate element
max_theta = np.amax(var_theta, axis=1)
eleTags = np.searchsorted(max_theta, dists)
# locate GLL-point
pntTags = np.zeros(ndists, dtype=int)
real_theta = np.zeros(ndists)
for idist, dist in enumerate(dists):
    if eleTags[idist] == 0 or eleTags[idist] == nele - 1:
        lbases = var_GLJ[:]
    else:
        lbases = var_GLL[:]
    theta_bounds = var_theta[eleTags[idist]]
    eta = (dist - theta_bounds[0]) / (theta_bounds[1] - theta_bounds[0]) * 2. - 1.
    pntTags[idist] = np.searchsorted(lbases, eta)
    if pntTags[idist] == nPntEdge:
        pntTags[idist] = nPntEdge - 1
    real_eta = lbases[pntTags[idist]]    
    real_theta[idist] = (real_eta + 1.) * .5 * (theta_bounds[1] - theta_bounds[0]) + theta_bounds[0]
if args.verbose:
    elapsed = time.clock() - clock0
    print('Locating points in distance done, ' + 
          '%f sec elapsed.\n' % (elapsed))   

###### prepare time steps
if args.verbose:
    clock0 = time.clock()
    print('Preparing timesteps...')
if nstep == 1:
    steps = np.array([0])
    dt = 0.
else:
    tend = args.tstart + (args.ntimesteps - 1) * args.time_interval
    times = np.arange(args.tstart, tend + args.time_interval * .1, args.time_interval)
    steps = np.searchsorted(var_time[:], times)
    # remove values out of range
    steps = steps[steps<nstep]
    if steps[0] == 0:
        steps = steps[steps>0]
        steps = np.insert(steps, 0, 0)
    else:
        steps = steps[steps>0]
    dt = var_time[1] - t0    
    nsteps = len(steps)
if args.verbose:
    elapsed = time.clock() - clock0
    print('    Number of snapshots: %d' % (nsteps))
    print('Preparing timesteps done, ' + 
          '%f sec elapsed.\n' % (elapsed))

# close serial input    
nu = np.zeros((nsteps, ndists), dtype=int)

def compute_nu(iproc):
    # if args.nproc == 1:
    #     nc_surf_local = Dataset(args.in_surface_nc, 'r')
    #     iproc = 0
    # else:
    #     # copy netcdf file for parallel access
    #     tempnc = args.out_vtk + '/surface_temp.nc' + str(iproc)
    #     shutil.copy(args.in_surface_nc, tempnc)
    #     nc_surf_local = Dataset(tempnc, 'r')

    # compute nu
    if args.verbose and iproc == 0:
        clock0 = time.clock()
        print('Computing nu...')
    for it, istep in enumerate(steps):
        if it % args.nproc != iproc: 
            continue
        eleTag_last = -1
        fourier_last = None
        for idist, dist in enumerate(dists):
            if eleTags[idist] == eleTag_last:
                fourier = fourier_last
            else:
                nce = nc_surfs[irank_ele[eleTags[idist]]]
                fourier_r = nce.variables['edge_' + str(eleTags[idist]) + 'r'][istep, :]
                fourier_i = nce.variables['edge_' + str(eleTags[idist]) + 'i'][istep, :]
                fourier = fourier_r[:] + fourier_i[:] * 1j
                fourier_last = fourier
                eleTag_last = eleTags[idist]
            nu_p_1 = int(len(fourier) / nPntEdge / 3)
            fmat = np.zeros((3, nu_p_1), dtype=fourier.dtype)
            for idim in np.arange(0, 3):
                for alpha in np.arange(nu_p_1):
                    fmat[idim, alpha] = fourier[idim * nPntEdge * nu_p_1 + pntTags[idist] * nu_p_1 + alpha]
            if args.component == 'R':
                fcoef = fmat[0, :] * np.cos(dist) - fmat[2, :] * np.sin(dist)
            elif args.component == 'T':
                fcoef = fmat[1, :]
            else:
                fcoef = fmat[0, :] * np.sin(dist) + fmat[2, :] * np.cos(dist)
            for ic in np.arange(nu_p_1 - 1, 1, -1):
                if np.abs(fcoef[ic]) >= np.abs(fcoef[0]) * args.error:
                    nu[it, idist] = ic
                    break
        if args.verbose:
            print('    Done with timestep t = %f s; tstep = %d / %d; iproc = %d' \
                % (var_time[istep], it + 1, nsteps, iproc))
    # close
    for nc_s in nc_surfs:
        nc_s.close() 
    
    # remove temp nc
    if args.nproc > 1:
        os.remove(tempnc)
    
    if args.verbose and iproc == 0:
        elapsed = time.clock() - clock0
        print('Computing nu done, ' + 
              '%f sec elapsed.' % (elapsed))

# compute nu in parallel
args.nproc = max(args.nproc, 1)
if args.nproc == 1:
    compute_nu(0)
else:
    with Pool(args.nproc) as p:
        p.map(compute_nu, range(0, args.nproc))
        
###### prepare output
nc_nu = Dataset(args.out_nu, 'w')
# time
ncdim_nstep = 'ncdim_' + str(nsteps)
nc_nu.createDimension(ncdim_nstep, size=nsteps)
var_time_out = nc_nu.createVariable('time_points', 
    float, (ncdim_nstep,))
var_time_out[:] = var_time[steps]
# dist
ncdim_ndist = 'ncdim_' + str(ndists)
nc_nu.createDimension(ncdim_ndist, size=ndists)
var_dist_out = nc_nu.createVariable('dist_points', 
    float, (ncdim_ndist,))
var_dist_out[:] = real_theta
# nu
var_nu_out = nc_nu.createVariable('nu', 
    int, (ncdim_nstep, ncdim_ndist))
var_nu_out[:, :] = nu
nc_nu.close()

