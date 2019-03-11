#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface2stations.py

Extract synthetics at given stations from a NetCDF database of surface
wavefield created by AxiSEM3D (named axisem3d_surface.nc by the solver)
and save them into a NetCDF waveform database (same as the built-in
NetCDF output axisem3d_synthetics.nc).

To see usage, type
python surface2stations.py -h
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
parser.add_argument('-m', '--multi_file', dest='multi_file', action='store_true', 
                    help='Does the NetCDF database consist of\n' +
                         'multiple files; default = False')                         
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
parser.add_argument('-F', '--order_Fourier', dest='order_Fourier', 
                    action='store', type=int, default=-1, 
                    help='only compute the specified order;\n' + 
                         'default = -1 (sum up all the orders)') 
parser.add_argument('-f', '--factor_Fourier', dest='factor_Fourier', 
                    action='store', type=complex, default=1.0, 
                    help='factor to scale the Fourier coefficients,\n' + 
                         'can be a complex number; default = 1.0')                         
parser.add_argument('-l', '--source_lat_lon', dest='source_lat_lon', 
                    action='store', nargs=2, type=float, default=None, 
                    help='specify source latitude and longitude;\n' + 
                         'default = None (use those in the solver)')                         
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')        
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset
from obspy.geodetics import gps2dist_azimuth
import os

################### TOOLS ###################
def rotation_matrix(theta, phi):
    return np.array([[np.cos(theta) * np.cos(phi), -np.sin(phi), np.sin(theta) * np.cos(phi)],
                     [np.cos(theta) * np.sin(phi), np.cos(phi), np.sin(theta) * np.sin(phi)],
                     [-np.sin(theta), 0., np.cos(theta)]])
                     
def latlon2thetaphi(lat, lon, flattening):
    temp = (1. - flattening) * (1. - flattening)
    return np.pi / 2. - np.arctan(temp * np.tan(np.radians(lat))), np.radians(lon)

def thetaphi2latlon(theta, phi, flattening):
    temp = 1. / (1. - flattening) / (1. - flattening)
    return np.degrees(np.arctan(temp * np.tan(np.pi / 2. - theta))), np.degrees(phi)
           
def thetaphi2xyz(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])

def xyz2thetaphi(xyz):
    theta = np.arccos(xyz[2])
    phi = np.arctan2(xyz[1], xyz[0])
    return theta, phi 

class SurfaceStation:
    # class variables and methods
    src_lat = None
    src_lon = None
    src_rmat = None
    def setSource(srclat, srclon, srcflattening):
        SurfaceStation.src_lat, SurfaceStation.src_lon = srclat, srclon
        srctheta, srcphi = latlon2thetaphi(srclat, srclon, srcflattening)
        SurfaceStation.src_rmat = rotation_matrix(srctheta, srcphi)
    
    def __init__(self, network, name):
        self.network = network
        self.name = name

    def setloc_geographic(self, lat, lon, flattening):
        self.lat = lat
        self.lon = lon
        theta, phi = latlon2thetaphi(lat, lon, flattening)
        xglb = thetaphi2xyz(theta, phi)
        xsrc = SurfaceStation.src_rmat.T.dot(xglb)
        self.dist, self.azimuth = xyz2thetaphi(xsrc)
        d, az, baz = gps2dist_azimuth(SurfaceStation.src_lat, SurfaceStation.src_lon, 
            self.lat, self.lon, a=1., f=flattening)
        self.baz = np.radians(baz)
    
    def setloc_source_centered(self, dist, azimuth, flattening):
        self.dist = dist
        self.azimuth = azimuth
        xsrc = thetaphi2xyz(dist, azimuth)
        xglb = SurfaceStation.src_rmat.dot(xsrc)
        theta, phi = xyz2thetaphi(xglb)
        self.lat, self.lon = thetaphi2latlon(theta, phi, flattening)
        d, az, baz = gps2dist_azimuth(SurfaceStation.src_lat, SurfaceStation.src_lon, 
            self.lat, self.lon, a=1., f=flattening)
        self.baz = np.radians(baz)

def interpLagrange(target, bases):
    nbases = len(bases)
    results = np.zeros(nbases)
    for dgr in np.arange(0, nbases):
        prod1 = 1.
        prod2 = 1.
        for i in np.arange(0, nbases):
            if i != dgr:
                prod1 *= target - bases[i]
                prod2 *= bases[dgr] - bases[i]
        results[dgr] = prod1 / prod2
    return results    
################### TOOLS ###################


###### read surface database
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
if args.source_lat_lon is not None:
    srclat = args.source_lat_lon[0]
    srclon = args.source_lat_lon[1]
else:    
    srclat = nc_surf.source_latitude
    srclon = nc_surf.source_longitude
srcdep = nc_surf.source_depth
srcflat = nc_surf.source_flattening
surfflat = nc_surf.surface_flattening
# time
var_time = nc_surf.variables['time_points'][:]
nstep = len(var_time)
solver_dtype = nc_surf.variables['time_points'].datatype
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

# set source
SurfaceStation.setSource(srclat, srclon, srcflat)

###### read station info
station_info = np.loadtxt(args.stations, dtype=str, ndmin=2)
stations = {}
buried = 0
largest_depth = -1.
largest_depth_station = ''
for ist in np.arange(0, len(station_info)):
    name = station_info[ist, 0]
    network = station_info[ist, 1]
    lat_theta = float(station_info[ist, 2])
    lon_phi = float(station_info[ist, 3])
    depth = float(station_info[ist, 5])
    key = network + '.' + name
    # ignore buried depth
    if depth > 0.: 
        buried += 1
        if largest_depth < depth:
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
        st.setloc_geographic(lat_theta, lon_phi, surfflat)
    else:
        lat_theta = np.radians(lat_theta)
        lon_phi = np.radians(lon_phi)
        st.setloc_source_centered(lat_theta, lon_phi, surfflat)
    stations[key] = st
    
# sort stations by distance    
station_list = list(stations.values())
station_list.sort(key=lambda st: st.dist)    
    
if buried > 0:
    print('Warning: Ignoring buried depth of %d stations;' % (buried))    
    print('         the deepest station %s is buried at %.0f m' \
        % (largest_depth_station, largest_depth))    
    
###### prepare output
nc_wave = Dataset(args.out_waveform_nc, 'w')
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
max_theta = np.amax(var_theta, axis=1)
eleTag_last = -1
for ist, station in enumerate(station_list):
    # waveform
    key = station.network + '.' + station.name
    var_wave = nc_wave.createVariable(key + '.' + args.components, 
        solver_dtype, (ncdim_nstep, 'ncdim_3'))
    # station info
    var_wave.latitude = station.lat
    var_wave.longitude = station.lon
    var_wave.depth = 0.
    # locate station
    eleTag = np.searchsorted(max_theta, station.dist)
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
    if eleTag != eleTag_last:
        nce = nc_surfs[irank_ele[eleTag]]
        fourier_r = nce.variables['edge_' + str(eleTag) + 'r'][:, :]
        fourier_i = nce.variables['edge_' + str(eleTag) + 'i'][:, :]
        fourier = fourier_r + fourier_i * 1j
        nu_p_1 = int(fourier_r.shape[1] / nPntEdge / 3)
        eleTag_last = eleTag
        
    # exp array    
    exparray = 2. * np.exp(np.arange(0, nu_p_1) * 1j * station.azimuth)
    exparray[0] = 1.
    if args.order_Fourier >= 0:
        assert args.order_Fourier < nu_p_1, 'Specified Fourier order %d greater than maximum %d' \
            % (args.order_Fourier, nu_p_1 - 1)
        exparray = exparray[args.order_Fourier]
    exparray *= args.factor_Fourier    
        
    # compute disp
    spz = np.zeros((nstep, 3), dtype=solver_dtype)
    fmat = fourier[:, :].reshape(nstep, 3, nPntEdge, nu_p_1)
    if args.order_Fourier >= 0:
        frow = fmat[:, :, :, args.order_Fourier]
        spz[:, :] = (frow * exparray).real.dot(weights)
    else:
        spz[:, :] = fmat.dot(exparray).real.dot(weights)
    # rotate
    disp = np.zeros((nstep, 3), dtype=solver_dtype)
    if args.components == 'SPZ':
        disp = spz
    else:
        ur = spz[:, 0] * np.sin(station.dist) + spz[:, 2] * np.cos(station.dist)
        ut = spz[:, 0] * np.cos(station.dist) - spz[:, 2] * np.sin(station.dist)    
        if args.components == 'ENZ':
            disp[:, 0] = -ut * np.sin(self.baz) + spz[:, 1] * np.cos(self.baz)
            disp[:, 1] = -ut * np.cos(self.baz) - spz[:, 1] * np.sin(self.baz)
            disp[:, 2] = ur    
        else:
            disp[:, 0] = ut
            disp[:, 1] = spz[:, 1]
            disp[:, 2] = ur
    var_wave[:, :] = disp[:, :]
    if args.verbose:
        print('Done with station %s, %d / %d' % (key, ist + 1, len(station_list)))

for nc_s in nc_surfs:
    nc_s.close() 
nc_wave.close()


