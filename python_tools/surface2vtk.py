#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface2vtk.py

Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver). Data
are presented on an unstructured mesh.

To see usage, type
python surface2vtk.py -h
'''
    
################### PARSER ###################
aim = '''Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver). Data
are presented on an unstructured mesh.'''

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
parser.add_argument('-m', '--multi_file', dest='multi_file', action='store_true', 
                    help='NetCDF database consists of multiple\n' +
                         'files (using netcdf_no_assemble);\n' +
                         'default = False')
parser.add_argument('-o', '--output', dest='out_vtk', 
                    action='store', type=str, required=True,
                    help='directory to store the vtk files\n' +
                         '<required>') 
parser.add_argument('-s', '--spatial_sampling', dest='spatial_sampling', 
                    action='store', type=float, required=True,
                    help='spatial sampling on surface (km)\n' +
                         '<required>') 
parser.add_argument('-md', '--min_dist', dest='min_dist', 
                    action='store', type=float, default=0.,
                    help='minimum distance (deg); default = 0')
parser.add_argument('-Md', '--max_dist', dest='max_dist', 
                    action='store', type=float, default=180.,
                    help='maximum distance (deg); default = 180')                          
parser.add_argument('-t0', '--tstart', dest='tstart', 
                    action='store', type=float, required=True,
                    help='start time of animation (sec)\n' +
                         '<required>') 
parser.add_argument('-dt', '--time_interval', dest='time_interval', 
                    action='store', type=float, required=True,
                    help='time interval between snapshots (sec)\n' +
                         '<required>') 
parser.add_argument('-nt', '--nsnapshots', dest='nsnapshots',
                    action='store', type=int, required=True,
                    help='number of snapshots <required>')
parser.add_argument('-N', '--norm', dest='norm', action='store_true', 
                    help='only dump displacement norm;\n' +
                         'default = False (dump 3D vector)')
parser.add_argument('-p', '--using_mpi', dest='using_mpi', action='store_true', 
                    help='parallel mode with MPI')                           
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
import pyvtk
import os
import time

# mpi
if args.using_mpi:
    from mpi4py import MPI
    mpi_size = MPI.COMM_WORLD.Get_size()
    mpi_rank = MPI.COMM_WORLD.Get_rank()
else:
    mpi_size = 1
    mpi_rank = 0

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
if args.multi_file:
    # read index file
    rank_edge = np.genfromtxt(args.in_surface_nc + '/rank_edge.txt', skip_header=1, dtype=str)
    ranks = rank_edge[:, 0].astype(int)
    edges = np.core.defchararray.replace(rank_edge[:, 1], 'edge_', '').astype(int)
    ranks_unique = np.unique(ranks)
    # open
    nc_surfs = []
    for irank in ranks_unique:
        fname = args.in_surface_nc + '/axisem3d_surface.nc.rank' + str(irank)
        nc = Dataset(fname, 'r')
        nc_surfs.append(nc)
        if args.verbose and mpi_rank == 0:
            print('Done opening nc file %s' % (fname))
    nc_surf = nc_surfs[0]
    # map
    edge_nc = {}
    for i in np.arange(len(edges)):
        edge_nc[edges[i]] = ranks_unique.tolist().index(ranks[i])
else:
    nc_surf = Dataset(args.in_surface_nc, 'r')
    if args.verbose and mpi_rank == 0:
        print('Done opening nc file %s' % (args.in_surface_nc))
    nc_surfs = [nc_surf]
    edge_nc = None


###### read surface database
if args.verbose and mpi_rank == 0:
    print()
    clock0 = time.clock()
    print('Reading global parameters...')
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
if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('Reading global parameters done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### create output directory
try:
    os.makedirs(args.out_vtk)
except OSError:
    pass
          
###### surface sampling
if args.verbose and mpi_rank == 0:
    clock0 = time.clock()
    print('Sampling surface...')
divisions = int(0.5 * np.pi * r_outer / (args.spatial_sampling * 1e3)) + 1
zmin = np.cos(np.radians(args.max_dist))
zmax = np.cos(np.radians(args.min_dist))
xyz, connect = SpherifiedCube(divisions, zmin, zmax)
nstation = len(xyz)
ncell = len(connect)
if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('    Number of sampling points: %d' % (nstation))
    print('    Number of quad cells: %d' % (ncell))
    print('Sampling surface done, ' + 
          '%f sec elapsed.\n' % (elapsed))
          
###### generate mesh vtk
if args.verbose and mpi_rank == 0:
    clock0 = time.clock()
    print('Generating vtk mesh...')
vtk_points = pyvtk.UnstructuredGrid(list(zip(xyz[:,0], xyz[:,1], xyz[:,2])), quad=connect)
if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('Generating vtk mesh done, ' + 
          '%f sec elapsed.\n' % (elapsed))
    
###### dist, azim
if args.verbose and mpi_rank == 0:
    clock0 = time.clock()
    print('Computing (distances, azimuths) of points...')    
# dists
dists = np.arccos(xyz[:, 2] / r_plot)
azims = np.arctan2(xyz[:, 1], xyz[:, 0])    
if args.verbose and mpi_rank == 0:
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

if args.verbose and mpi_rank == 0:
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
if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('Locating points in distance done, ' + 
          '%f sec elapsed.\n' % (elapsed))    

###### prepare time steps
if args.verbose and mpi_rank == 0:
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
if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('    Number of snapshots: %d' % (len(steps)))
    print('Preparing timesteps done, ' + 
          '%f sec elapsed.\n' % (elapsed))


# write vtk
if args.verbose and mpi_rank == 0:
    clock0 = time.clock()
    print('Generating snapshot...')
    
for it, istep in enumerate(steps):
    # min/max step
    if args.min_step is not None:
        if it < args.min_step:
            continue
    if args.max_step is not None:
        if it > args.max_step:
            continue  
    if it % mpi_size != mpi_rank:
        continue         
    
    if args.verbose:
        clock0s = time.clock()
        
    # output        
    if args.norm:
        disp_norm = np.zeros(nstation)
    else:
        disp = np.zeros((nstation, 3))
        
    eleTag_last = -1
    fmat_last = None
    for istation, dist in enumerate(dists):
        etag = eleTags[istation]
        if etag == eleTag_last:
            fmat = fmat_last
            nu_p_1 = fmat.shape[-1]
        else:
            if (edge_nc is None):
                nc = nc_surf
            else:
                nc = nc_surfs[edge_nc[etag]]
            fourier_r = nc.variables['edge_' + str(etag) + 'r'][istep, :]
            fourier_i = nc.variables['edge_' + str(etag) + 'i'][istep, :]
            fourier = fourier_r[:] + fourier_i[:] * 1j
            nu_p_1 = int(len(fourier) / nPntEdge / 3)
            fmat = fourier.reshape(3, nPntEdge, nu_p_1)
            fmat_last = fmat
            eleTag_last = etag
        wdotf = np.tensordot(weights[istation], fmat, ([0], [1]))
        exparray = 2. * np.exp(np.arange(0, nu_p_1) * 1j * azims[istation])
        exparray[0] = 1.
        spz = wdotf.dot(exparray).real
        if args.norm:
            disp_norm[istation] = np.linalg.norm(spz)
        else:
            disp[istation, 0] = spz[0] * np.cos(dist) - spz[2] * np.sin(dist)
            disp[istation, 1] = spz[1]
            disp[istation, 2] = spz[0] * np.sin(dist) + spz[2] * np.cos(dist)
    if args.norm:
        vtk = pyvtk.VtkData(vtk_points,
            pyvtk.PointData(pyvtk.Scalars(disp_norm, name='disp_norm')), 
            'surface animation')
    else:
        vtk = pyvtk.VtkData(vtk_points,
            pyvtk.PointData(pyvtk.Vectors(disp, name='disp_RTZ')),
            'surface animation')
    vtk.tofile(args.out_vtk + '/surface_vtk.' + str(it) + '.vtk', 'binary')
    if args.verbose:
        elapsed = time.clock() - clock0s
        print('    Done with snapshot t = %f s; tstep = %d / %d, rank = %d, elapsed = %f' \
            % (var_time[istep], it + 1, len(steps), mpi_rank, elapsed))

if args.verbose and mpi_rank == 0:
    elapsed = time.clock() - clock0
    print('Generating snapshots done, ' + 
          '%f sec elapsed.' % (elapsed))

