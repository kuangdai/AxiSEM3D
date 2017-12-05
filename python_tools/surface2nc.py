#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface2nc.py

Extract synthetics at given stations from a NetCDF database of surface
wavefield created by AxiSEM3D (named axisem3d_surface.nc by the solver)
and save them into a NetCDF waveform database (same as the built-in
NetCDF output axisem3d_synthetics.nc).

To see usage, type
python surface2nc.py -h
'''

################### PARSER ###################
aim = '''Extract synthetics at given stations from a NetCDF database of surface
wavefield created by AxiSEM3D (named axisem3d_surface.nc by the solver)
and save them into a NetCDF waveform database (same as the built-in
NetCDF output axisem3d_synthetics.nc).'''

notes = '''One may further use nc2ascii.py to convert the output into ascii format.
'''

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=aim, epilog=notes, 
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-i', '--input', dest='in_surface_nc', 
					action='store', type=str, required=True,
                    help='NetCDF database of surface wavefield\n' + 
					     'created by AxiSEM3D <required>')
parser.add_argument('-o', '--output', dest='out_waveform_nc', 
					action='store', type=str, required=True,
                    help='NetCDF waveform database to store the\n' + 
					     'extracted synthetics <required>') 
parser.add_argument('-s', '--stations', dest='stations', 
					action='store', type=str, required=True,
                    help='list of stations, see OUT_STATIONS_FILE\n' + 
					     'in inparam.time_src_recv <required>') 
parser.add_argument('-r', '--crdsys', dest='crdsys', action='store', 
					type=str, default='geographic', 
					choices=['geographic', 'source-centered'],
					help='coordinate system used in the station list,\n' + 
					     'see OUT_STATIONS_SYSTEM in inparam.time_src_recv;\n' +
						 'default = geographic')
parser.add_argument('-d', '--duplicated', dest='duplicated', action='store', 
                    type=str, default='rename', 
					choices=['ignore', 'rename', 'error'],
                    help='how to haddle stations with the same network\n' + 
                         'and name, see OUT_STATIONS_DUPLICATED in\n' + 
                         'inparam.time_src_recv; default = rename') 
parser.add_argument('-c', '--components', dest='components', action='store', 
					type=str, default='RTZ', choices=['RTZ', 'ENZ', 'SPZ'],
					help='seismogram components, see OUT_STATIONS_COMPONENTS\n' + 
					     'in inparam.time_src_recv; default = RTZ')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')		
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset
from surface_utils import interpLagrange
from surface_utils import SurfaceStation
from obspy.geodetics.base import gps2dist_azimuth

###### read surface database
nc_surf = Dataset(args.in_surface_nc, 'r', format='NETCDF4')
# global attribute
srclat = nc_surf.source_latitude
srclon = nc_surf.source_longitude
srcdep = nc_surf.source_depth
srcflat = nc_surf.source_flattening
surfflat = nc_surf.surface_flattening
# time
var_time = nc_surf.variables['time_points']
nstep = len(var_time)
solver_dtype = var_time.datatype
# theta
var_theta = nc_surf.variables['theta']
nele = len(var_theta)
# GLL and GLJ
var_GLL = nc_surf.variables['GLL']
var_GLJ = nc_surf.variables['GLJ']
nPntEdge = len(var_GLL)

###### read station info
station_info = np.loadtxt(args.stations, dtype=str)
stations = {}
buried = 0
largest_depth = -1.
largest_depth_station = ''
for ist in np.arange(0, len(station_info)):
	name = station_info[ist, 0]
	network = station_info[ist, 1]
	lat_theta = float(station_info[ist, 2])
	lon_phi = float(station_info[ist, 3])
	depth = float(station_info[ist, len(station_info[ist]) - 1])
	key = network + '.' + name
	# ignore buried depth
	if depth > 0.: 
		buried += 1
		if (largest_depth < depth):
			largest_depth = depth
			largest_depth_station = key
	# duplicated stations	
	if key in stations:
		if args.duplicated == 'error':
			assert False, 'Duplicated station %s found in %s' \
				% (key, args.stations)
		if args.duplicated == 'rename':
			append = 0
			nameOriginal = name
			while key in stations:
				append += 1
				name = nameOriginal + '__DUPLICATED' + str(append)
				key = network + '.' + name
		if args.duplicated == 'ignore':
			continue
	st = SurfaceStation(network, name)
	# coordinate system
	if args.crdsys == 'geographic':
		st.setloc_geographic(lat_theta, lon_phi, surfflat, 
			srclat, srclon, srcflat)
	else:
		lat_theta = np.radians(lat_theta)
		lon_phi = np.radians(lon_phi)
		st.setloc_source_centered(lat_theta, lon_phi, surfflat, 
			srclat, srclon, srcflat)
	stations[key] = st
if buried > 0:
	print('Warning: Ignoring buried depth of %d stations;' % (buried))	
	print('         the deepest station %s is buried at %.0f m' \
		% (largest_depth_station, largest_depth))	
	
###### prepare output
nc_wave = Dataset(args.out_waveform_nc, 'w', format='NETCDF4')
# global attribute
nc_wave.source_latitude = srclat
nc_wave.source_longitude = srclon
nc_wave.source_depth = srcdep
# time
ncdim_nstep = 'ncdim_' + str(nstep)
nc_wave.createDimension(ncdim_nstep, size=nstep)
var_time_out = nc_wave.createVariable('time_points', 
	solver_dtype, (ncdim_nstep,))
var_time_out[:] = var_time[:]

# waveforms
nc_wave.createDimension('ncdim_3', size=3)
for ist, station in enumerate(stations.values()):
	# waveform
	key = station.network + '.' + station.name
	var_wave = nc_wave.createVariable(key + '.' + args.components, 
		solver_dtype, (ncdim_nstep, 'ncdim_3'))
	# station info
	var_wave.latitude = station.lat
	var_wave.longitude = station.lon
	var_wave.depth = 0.
	# locate station
	eleTag = -1
	for iele in np.arange(0, nele):
		if station.dist <= var_theta[iele, 1]:
			eleTag = iele
			break
	assert eleTag >= 0, 'Fail to locate station %s, dist = %f' \
		% (key, station.dist)
	theta0 = var_theta[eleTag, 0]
	theta1 = var_theta[eleTag, 1]
	eta = (station.dist - theta0) / (theta1 - theta0) * 2. - 1.
	# weights considering axial condition
	if eleTag == 0 or eleTag == nele - 1:
		weights = interpLagrange(eta, var_GLJ)
	else:
		weights = interpLagrange(eta, var_GLL)
		
	# Fourier
	# NOTE: change to stepwise if memory issue occurs
	fourier_r = nc_surf.variables['edge_' + str(eleTag) + 'r']
	fourier_i = nc_surf.variables['edge_' + str(eleTag) + 'i']
	fourier = fourier_r[:, :] + fourier_i[:, :] * 1j
	nu_p_1 = int(fourier_r[:, :].shape[1] / nPntEdge / 3)
	exparray = 2. * np.exp(np.arange(0, nu_p_1) * 1j * station.azimuth)
	exparray[0] = 1.
	# back azimuth
	if args.components == 'ENZ':
		d, az, baz = gps2dist_azimuth(srclat, srclon, 
			station.lat, station.lon, a=1., f=surfflat)
		baz = np.radians(baz)
	# compute disp
	disp = np.zeros((nstep, 3), dtype=solver_dtype)
	for istep in np.arange(0, nstep):
		spz = np.zeros(3)
		for idim in np.arange(0, 3):
			start = idim * nPntEdge * nu_p_1
			end = idim * nPntEdge * nu_p_1 + nPntEdge * nu_p_1
			fmat = fourier[istep, start:end].reshape(nPntEdge, nu_p_1)
			spz[idim] = weights.dot(fmat.dot(exparray).real)
		if args.components == 'SPZ':
			disp[istep, :] = spz
		else:	
			ur = spz[0] * np.sin(station.dist) + spz[2] * np.cos(station.dist)
			ut = spz[0] * np.cos(station.dist) - spz[2] * np.sin(station.dist)	
			if args.components == 'ENZ':
				disp[istep, 0] = -ut * np.sin(baz) + spz[1] * np.cos(baz)
				disp[istep, 1] = -ut * np.cos(baz) - spz[1] * np.sin(baz)
				disp[istep, 2] = ur	
			else:
				disp[istep, 0] = ut
				disp[istep, 1] = spz[1]
				disp[istep, 2] = ur
	var_wave[:, :] = disp[:, :]
	if args.verbose:
		print('Done with station %s, %d / %d' % (key, ist + 1, len(stations)))
		
nc_surf.close()
nc_wave.close()


