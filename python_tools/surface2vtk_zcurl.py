#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface2vtk_zcurl.py

Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver). Data
are presented on an unstructured mesh.
-------------------------------------------------------------------
This script visualises displacement curl in the vertical direction;
gradient calculation is FD-based and is NOT accurate.
-------------------------------------------------------------------

To see usage, type
python surface2vtk_zcurl.py -h
'''
    
################### PARSER ###################
aim = '''Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver). Data
are presented on an unstructured mesh.
-------------------------------------------------------------------
This script visualises displacement curl in the vertical direction;
gradient calculation is FD-based and is NOT accurate.
-------------------------------------------------------------------'''

notes = '''Parallelise data processing using --nporc option.
Animate the VKT files with Paraview.
 
'''

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=aim, epilog=notes, 
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-i', '--input', dest='in_surface_nc', 
                    action='store', type=str, required=True,
                    help='NetCDF database of surface wavefield\n' + 
                         'created by AxiSEM3D <required>')
parser.add_argument('-o', '--output', dest='out_vtk', 
                    action='store', type=str, required=True,
                    help='directory to store the vtk files\n' +
                         '<required>') 
parser.add_argument('-s', '--spatial_sampling', dest='spatial_sampling', 
                    action='store', type=float, required=True,
                    help='spatial sampling on surface (km)\n' +
                         '<required>') 
parser.add_argument('-m', '--min_dist', dest='min_dist', 
                    action='store', type=float, default=0.,
                    help='minimum distance (deg); default = 0')
parser.add_argument('-M', '--max_dist', dest='max_dist', 
                    action='store', type=float, default=180.,
                    help='maximum distance (deg); default = 180')                          
parser.add_argument('-t', '--tstart', dest='tstart', 
                    action='store', type=float, required=True,
                    help='start time of animation (sec)\n' +
                         '<required>') 
parser.add_argument('-d', '--time_interval', dest='time_interval', 
                    action='store', type=float, required=True,
                    help='time interval between snapshots (sec)\n' +
                         '<required>') 
parser.add_argument('-n', '--nsnapshots', dest='nsnapshots',
                    action='store', type=int, required=True,
                    help='number of snapshots <required>')
parser.add_argument('-p', '--nproc', dest='nproc', action='store', 
                    type=int, default=1, 
                    help='number of processors; default = 1')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')        
# hidden options
# timestep range, used to take a subset of the timesteps
# determined by tstart, time_interval and nsnapshots
parser.add_argument('--min_step', dest='min_step',
                    action='store', type=int, default=None,
                    help=argparse.SUPPRESS)
parser.add_argument('--max_step', dest='max_step',
                    action='store', type=int, default=None,
                    help=argparse.SUPPRESS)
args = parser.parse_args()

################### PARSER ###################


import numpy as np
from netCDF4 import Dataset
import pyvtk, os, shutil
from multiprocessing import Pool
import time

# slightly increase the radius for plot
r_plot = 1.0001

################### MESH TOOLS ###################
# https://github.com/caosdoar/spheres/blob/master/src/spheres.cpp

CubeToSphere_origins = np.array([
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, -1.0, 1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0]])

CubeToSphere_rights = np.array([
    [2.0, 0.0, 0.0],
    [0.0, 0.0, 2.0],
    [-2.0, 0.0, 0.0],
    [0.0, 0.0, -2.0],
    [2.0, 0.0, 0.0],
    [2.0, 0.0, 0.0]])
    
CubeToSphere_ups = np.array([
    [0.0, 2.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 0.0, 2.0],
    [0.0, 0.0, -2.0]])
    
def SpherifiedCube(divisions, zmin, zmax, zsort=True):
    step = 1.0 / divisions
    step3 = np.array([step, step, step])
    # vertices
    k = divisions + 1
    n_vertices_full = 6 * k ** 2
    xyz_full = np.zeros((n_vertices_full, 3))
    mask_full = np.zeros(n_vertices_full, dtype=int)
    ipnt = 0
    for face in np.arange(6):
        origin = CubeToSphere_origins[face]
        right = CubeToSphere_rights[face]
        up = CubeToSphere_ups[face]
        for j in np.arange(k):
            j3 = np.array([j, j, j])
            for i in np.arange(k):
                i3 = np.array([i, i, i])
                p = origin + step3 * (i3 * right + j3 * up)
                p2 = p * p
                z = p[2] * np.sqrt(1.0 - 0.5 * (p2[0] + p2[1]) + p2[0] * p2[1] / 3.0)
                if z >= zmin and z <= zmax:
                    xyz_full[ipnt, 0] = p[0] * np.sqrt(1.0 - 0.5 * (p2[1] + p2[2]) + p2[1] * p2[2] / 3.0)
                    xyz_full[ipnt, 1] = p[1] * np.sqrt(1.0 - 0.5 * (p2[2] + p2[0]) + p2[2] * p2[0] / 3.0)
                    xyz_full[ipnt, 2] = z
                    mask_full[ipnt] = 1
                else:
                    mask_full[ipnt] = 0
                ipnt += 1
    
    # form index map from full to need and xyz_need
    n_vertices_need = np.sum(mask_full)
    xyz_need = np.zeros((n_vertices_need, 3))
    idx_need = np.zeros(n_vertices_need, dtype=int)
    f2n = np.zeros(n_vertices_full, dtype=int)
    idx_current = 0
    for imask, needed in enumerate(mask_full):
        if needed:
            xyz_need[idx_current] = xyz_full[imask]
            idx_need[idx_current] = imask
            f2n[imask] = idx_current
            idx_current += 1
        else:
            f2n[imask] = -1
            
    # sort by z (theta) to reduce IO access to netcdf database
    if zsort:
        idx_sort = np.argsort(xyz_need[:, 2])
        xyz_need = xyz_need[idx_sort]
        # no idea how this works, but it does...
        idx_revs = np.argsort(idx_sort)
        f2n[idx_need] = f2n[idx_need][idx_revs]
                
	# cells
    connectivity = []
    for face in np.arange(6):
        for j in np.arange(divisions):
            for i in np.arange(divisions):
                a = f2n[(face * k + j) * k + i]
                b = f2n[(face * k + j) * k + i + 1]
                c = f2n[(face * k + j + 1) * k + i]
                d = f2n[(face * k + j + 1) * k + i + 1]
                # skip if one of the vertices are out of z-range
                if a < 0 or b < 0 or c < 0 or d < 0:
                    continue
                connectivity.append([a, c, d, b])
    aconnectivity = np.array(connectivity)
    
    xyz_need *= r_plot
    return xyz_need, aconnectivity

################### MESH TOOLS ###################


###### read surface database
if args.verbose:
    clock0 = time.clock()
    print('Reading global parameters...')
nc_surf = Dataset(args.in_surface_nc, 'r')
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
if args.verbose:
    elapsed = time.clock() - clock0
    print('Reading global parameters done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### create output directory
try:
    os.makedirs(args.out_vtk)
except OSError:
    pass
          
###### surface sampling
if args.verbose:
    clock0 = time.clock()
    print('Sampling surface...')
divisions = int(0.5 * np.pi * r_outer / (args.spatial_sampling * 1e3)) + 1
zmin = np.cos(np.radians(args.max_dist))
zmax = np.cos(np.radians(args.min_dist))
xyz, connect = SpherifiedCube(divisions, zmin, zmax)
nstation = len(xyz)
ncell = len(connect)
if args.verbose:
    elapsed = time.clock() - clock0
    print('    Number of sampling points: %d' % (nstation))
    print('    Number of quad cells: %d' % (ncell))
    print('Sampling surface done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### generate mesh vtk
if args.verbose:
    clock0 = time.clock()
    print('Generating vtk mesh...')
vtk_points = pyvtk.UnstructuredGrid(list(zip(xyz[:,0], xyz[:,1], xyz[:,2])), quad=connect)
if args.verbose:
    elapsed = time.clock() - clock0
    print('Generating vtk mesh done, ' + 
          '%f sec elapsed.\n' % (elapsed))
    
###### dist, azim
if args.verbose:
    clock0 = time.clock()
    print('Computing (distances, azimuths) of points...')    
# dists
dists = np.arccos(xyz[:, 2] / r_plot)
azims = np.arctan2(xyz[:, 1], xyz[:, 0])    
if args.verbose:
    elapsed = time.clock() - clock0
    print('Computing (distances, azimuths) of points done, ' + 
          '%f sec elapsed.\n' % (elapsed))

###### prepare theta
def interpLagrange(target, lbases):
    nrow, ncol = lbases.shape
    results = np.zeros((nrow, ncol))
    target_dgr = np.tile(np.array([target]).T, (1, ncol - 1))
    for dgr in np.arange(0, ncol):
        lbases_dgr = np.tile(lbases[:, [dgr]], (1, ncol - 1))
        lbases_sub = np.delete(lbases, dgr, axis=1)
        results[:, dgr] = np.prod(target_dgr - lbases_sub, axis=1) / \
                          np.prod(lbases_dgr - lbases_sub, axis=1)
    return results

if args.verbose:
    clock0 = time.clock()
    print('Locating points in distance...')
# locate element
max_theta = np.amax(var_theta, axis=1)
eleTags = np.searchsorted(max_theta, dists)
# compute weights
lbases = np.tile(var_GLL, (nstation, 1))
for ist in np.arange(nstation):
    if eleTags[ist] == 0 or eleTags[ist] == nele - 1:
        lbases[ist, :] = var_GLJ[:]
theta_bounds = var_theta[eleTags, :]
etas = (dists - theta_bounds[:, 0]) / (theta_bounds[:, 1] - theta_bounds[:, 0]) * 2. - 1.
weights = interpLagrange(etas, lbases)
######################################################
# another distance nearby to compute gradient
delta = 1e-6
# try plus first
dists1 = dists + delta
eleTags1 = np.searchsorted(max_theta, dists1)
# use minus if plus locates in another element
idx_minus = np.nonzero(eleTags1 - eleTags)
dists1[idx_minus] = dists[idx_minus] - delta
# check
eleTags1 = np.searchsorted(max_theta, dists1)
assert len(eleTags[np.nonzero(eleTags1 - eleTags)]) == 0
# weights
etas1 = (dists1 - theta_bounds[:, 0]) / (theta_bounds[:, 1] - theta_bounds[:, 0]) * 2. - 1.
weights1 = interpLagrange(etas1, lbases)
######################################################
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
    tend = args.tstart + (args.nsnapshots - 1) * args.time_interval
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
if args.verbose:
    elapsed = time.clock() - clock0
    print('    Number of snapshots: %d' % (len(steps)))
    print('Preparing timesteps done, ' + 
          '%f sec elapsed.\n' % (elapsed))

# close serial input    
nc_surf.close()
def write_vtk(iproc):
    if args.nproc == 1:
        nc_surf_local = Dataset(args.in_surface_nc, 'r')
        iproc = 0
    else:
        # copy netcdf file for parallel access
        tempnc = args.out_vtk + '/surface_temp.nc' + str(iproc)
        shutil.copy(args.in_surface_nc, tempnc)
        nc_surf_local = Dataset(tempnc, 'r')

    # write vtk
    if args.verbose and iproc == 0:
        clock0 = time.clock()
        print('Generating snapshot...')
    for it, istep in enumerate(steps):
        if it % args.nproc != iproc: 
            continue
        if args.min_step is not None:
            if it < args.min_step:
                continue
        if args.max_step is not None:
            if it > args.max_step:
                continue  
        disp_curl = np.zeros(nstation)
        eleTag_last = -1
        fourier_last = None
        for istation, dist in enumerate(dists):
            # poles cannot be computed correctly
            if (dist < delta or dist > np.pi - delta):
                disp_curl[istation] = 0.
                continue
            if eleTags[istation] == eleTag_last:
                fourier = fourier_last
            else:
                fourier_r = nc_surf_local.variables['edge_' + str(eleTags[istation]) + 'r'][istep, :]
                fourier_i = nc_surf_local.variables['edge_' + str(eleTags[istation]) + 'i'][istep, :]
                fourier = fourier_r[:] + fourier_i[:] * 1j
                fourier_last = fourier
                eleTag_last = eleTags[istation]
            nu_p_1 = int(len(fourier) / nPntEdge / 3)
            wdotf = np.zeros((3, nu_p_1), dtype=fourier.dtype)
            wdotf1 = np.zeros((3, nu_p_1), dtype=fourier.dtype)
            for idim in np.arange(0, 3):
                start = idim * nPntEdge * nu_p_1
                end = idim * nPntEdge * nu_p_1 + nPntEdge * nu_p_1
                fmat = fourier[start:end].reshape(nPntEdge, nu_p_1)
                wdotf[idim] = weights[istation].dot(fmat)
                wdotf1[idim] = weights1[istation].dot(fmat)
            exparray = 2. * np.exp(np.arange(0, nu_p_1) * 1j * azims[istation])
            exparray[0] = 1.
            exparray1 = 2. * np.exp(np.arange(0, nu_p_1) * 1j * (azims[istation] + delta))
            exparray1[0] = 1.
            spz = wdotf.dot(exparray).real
            spz_dist1 = wdotf1.dot(exparray).real
            spz_azim1 = wdotf.dot(exparray1).real
            uR = spz[0] * np.cos(dist) - spz[2] * np.sin(dist)
            uR_azim1 = spz_azim1[0] * np.cos(dist) - spz_azim1[2] * np.sin(dist)
            duR = (uR_azim1 - uR) / delta / np.sin(dist)
            uT = spz[1]
            uT_dist1 = spz_dist1[1]
            duT = (uT_dist1 - uT) / (dists1[istation] - dist)
            disp_curl[istation] = duR - duT
        vtk = pyvtk.VtkData(vtk_points,
            pyvtk.PointData(pyvtk.Scalars(disp_curl, name='disp_curl')), 
            'surface animation')
        vtk.tofile(args.out_vtk + '/surface_vtk_zcurl.' + str(it) + '.vtk', 'binary')
        if args.verbose:
            print('    Done with snapshot t = %f s; tstep = %d / %d; iproc = %d' \
                % (var_time[istep], it + 1, len(steps), iproc))
    # close
    nc_surf_local.close()
    
    # remove temp nc
    if args.nproc > 1:
        os.remove(tempnc)
    
    if args.verbose and iproc == 0:
        elapsed = time.clock() - clock0
        print('Generating snapshots done, ' + 
              '%f sec elapsed.' % (elapsed))

# write_vtk in parallel
args.nproc = max(args.nproc, 1)
if args.nproc == 1:
    write_vtk(0)
else:
    with Pool(args.nproc) as p:
        p.map(write_vtk, range(0, args.nproc))
        

