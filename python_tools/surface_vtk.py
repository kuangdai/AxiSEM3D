#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface_vtk.py

Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver).

To see usage, type
python surface_vtk.py -h
'''
	
################### PARSER ###################
aim = '''Generate VTK animations from a NetCDF database of surface wavefield 
created by AxiSEM3D (named axisem3d_surface.nc by the solver).'''

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
parser.add_argument('-c', '--channels', dest='channels', action='store', 
					type=str, default='RTZ', choices=['RTZ', 'ENZ', 'SPZ'],
					help='channels to be animated; default = RTZ')
parser.add_argument('-p', '--nproc', dest='nproc', action='store', 
					type=int, default=1, 
					help='number of processors; default = 1')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
					help='verbose mode')		
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset
from surface_utils import interpLagrange
from surface_utils import SurfaceStation, latlon2thetaphi, thetaphi2xyz
import pyvtk, os, shutil
from multiprocessing import Pool

###### read surface database
nc_surf = Dataset(args.in_surface_nc, 'r', format='NETCDF4')
# global attribute
srclat = nc_surf.source_latitude
srclon = nc_surf.source_longitude
srcdep = nc_surf.source_depth
srcflat = nc_surf.source_flattening
surfflat = nc_surf.surface_flattening
r_outer = nc_surf.radius
# time
var_time = nc_surf.variables['time_points']
nstep = len(var_time)
assert nstep > 0, 'Zero time steps'
t0 = var_time[0]
solver_dtype = var_time.datatype
# theta
var_theta = nc_surf.variables['theta']
nele = len(var_theta)
# GLL and GLJ
var_GLL = nc_surf.variables['GLL']
var_GLJ = nc_surf.variables['GLJ']
nPntEdge = len(var_GLL)

###### create stations
ndist = int(np.pi * r_outer / (args.spatial_sampling * 1e3)) + 1
dists = np.linspace(0, np.pi, num=ndist, endpoint=True)
stations = []
for dist in dists:
	r = r_outer * np.sin(dist)
	nazi = int(2. * np.pi * r / (args.spatial_sampling * 1e3))
	azis = np.linspace(0, 2. * np.pi, num=nazi, endpoint=False)
	for azi in azis:
		st = SurfaceStation('', '')
		st.setloc_source_centered(dist, azi, surfflat, 
			srclat, srclon, srcflat)
		stations.append(st)

if args.verbose:
	print('Number of sampling points on surface: %d' % (len(stations)))
	
###### prepare theta
nstation = len(stations)
eleTags = np.zeros(nstation, dtype=int)
weights = np.zeros((nstation, nPntEdge))
distLast = -1.
x = np.zeros(nstation)
y = np.zeros(nstation)
z = np.zeros(nstation)
for ist, station in enumerate(stations):
	# coordinates
	theta, phi = latlon2thetaphi(station.lat, station.lon, surfflat)
	xyz = thetaphi2xyz(theta, phi)
	x[ist] = xyz[0]
	y[ist] = xyz[1]
	z[ist] = xyz[2]
	# ele and weights
	if np.isclose(station.dist, distLast):
		weights[ist, :] = weights[ist - 1, :]
		eleTags[ist] = eleTags[ist - 1]
		continue
	# locate station
	eleTag = -1
	for iele in np.arange(0, nele):
		if station.dist <= max(var_theta[iele, 0:1]):
			eleTag = iele
			break
	assert eleTag >= 0, 'Fail to locate point, dist = %f' \
		% (station.dist)
	theta0 = var_theta[eleTag, 0]
	theta1 = var_theta[eleTag, 1]
	eta = (station.dist - theta0) / (theta1 - theta0) * 2. - 1.
	# weights considering axial condition
	if eleTag == 0 or eleTag == nele - 1:
		weights[ist, :] = interpLagrange(eta, var_GLJ)
	else:
		weights[ist, :] = interpLagrange(eta, var_GLL)
	eleTags[ist] = eleTag
	distLast = station.dist
vtk_points = pyvtk.UnstructuredGrid(list(zip(x, y, z)), range(len(stations)))
	
###### prepare time steps
if nstep == 1:
	steps = np.array([0])
	dt = 0.
else:
	dt = var_time[1] - t0
	istart = max(int(round((args.tstart - t0) / dt)), 0)
	dtsteps = max(int(round(args.time_interval / dt)), 1)
	iend = min(istart + dtsteps * (args.nsnapshots - 1) + 1, nstep)
	steps = np.arange(istart, iend, dtsteps)
		
if args.verbose:
	print('Number of snapshots: %d' % (len(steps)))
nc_surf.close()
	
###### IO
# create output directory
try:
	os.makedirs(args.out_vtk)
except OSError:
	pass

def write_vtk(iproc):
	if args.nproc == 1:
		nc_surf_local = Dataset(args.in_surface_nc, 'r', format='NETCDF4')
	else:
		# copy netcdf file for parallel access
		tempnc = args.out_vtk + '/surface_temp.nc' + str(iproc)
		shutil.copy(args.in_surface_nc, tempnc)
		nc_surf_local = Dataset(tempnc, 'r', format='NETCDF4')

	# write vtk
	for it, istep in enumerate(steps):
		if (it % args.nproc != iproc): 
			continue
		disp = np.zeros((nstation, 3))
		eleTagLast = -1
		for ist, station in enumerate(stations):
			# Fourier
			if eleTags[ist] == eleTagLast:
				fourier = fourierLast
			else: 
				fourier_r = nc_surf_local.variables['edge_' + str(eleTags[ist]) + 'r'][istep, :]
				fourier_i = nc_surf_local.variables['edge_' + str(eleTags[ist]) + 'i'][istep, :]
				fourier = fourier_r[:] + fourier_i[:] * 1j
				fourierLast = fourier
				eleTagLast = eleTags[ist]
			nu_p_1 = int(len(fourier) / nPntEdge / 3)
			exparray = 2. * np.exp(np.arange(0, nu_p_1) * 1j * station.azimuth)
			exparray[0] = 1.
			# compute disp
			spz = np.zeros(3)
			for idim in np.arange(0, 3):
				start = idim * nPntEdge * nu_p_1
				end = idim * nPntEdge * nu_p_1 + nPntEdge * nu_p_1
				fmat = fourier[start:end].reshape(nPntEdge, nu_p_1)
				spz[idim] = weights[ist].dot(fmat.dot(exparray).real)
			if args.channels == 'SPZ':
				disp[ist, :] = spz
			else:	
				ur = spz[0] * np.sin(station.dist) + spz[2] * np.cos(station.dist)
				ut = spz[0] * np.cos(station.dist) - spz[2] * np.sin(station.dist)	
				if args.channels == 'ENZ':
					disp[ist, 0] = -ut * np.sin(self.baz) + spz[1] * np.cos(self.baz)
					disp[ist, 1] = -ut * np.cos(self.baz) - spz[1] * np.sin(self.baz)
					disp[ist, 2] = ur	
				else:
					disp[ist, 0] = ut
					disp[ist, 1] = spz[1]
					disp[ist, 2] = ur
		vtk = pyvtk.VtkData(vtk_points,
			pyvtk.PointData(pyvtk.Vectors(disp, name='disp_'+args.channels)),
			'surface animation')
		vtk.tofile(args.out_vtk + '/surface_vtk.' + str(it) + '.vtk', 'binary')
		if args.verbose:
			print('Done with snapshot t = %f s; tstep = %d / %d; iproc = %d' \
				% (istep * dt + t0, it + 1, len(steps), iproc))
	# close
	nc_surf_local.close()
	
	# remove temp nc
	if args.nproc > 1:
		os.remove(tempnc)

# write_vtk in parallel
args.nproc = max(args.nproc, 1)
if args.nproc == 1:
	write_vtk(0)
else:
	with Pool(args.nproc) as p:
		p.map(write_vtk, range(0, args.nproc))
		

