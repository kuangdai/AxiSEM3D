#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
nc2ascii.py

Extract synthetics from a NetCDF waveform database created by AxiSEM3D
(named axisem3d_synthetics.nc by the solver) and save them in ascii format.

To see usage, type
python nc2ascii.py -h
'''

################### PARSER ###################

aim = '''Extract synthetics from a NetCDF waveform database created by AxiSEM3D
(named axisem3d_synthetics.nc by the solver) and save them in ascii format.'''

notes = '''String replacement rules for filename, header and footer:
  @NW@ -> network name
  @ST@ -> station name
  @CH@ -> channel
  @NS@ -> number of steps
  @DT@ -> time step
  @SR@ -> sampling rate
  @T0@ -> start time w.r.t. source origin
  @T1@ -> end time w.r.t. source origin
  
'''

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=aim, epilog=notes, 
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-i', '--input', dest='in_nc_file', action='store', 
                    type=str, required=True,
                    help='NetCDF waveform database created by AxiSEM3D\n' +
                         '<required>') 
parser.add_argument('-o', '--output', dest='out_ascii_dir', action='store', 
                    type=str, required=True,
                    help='directory to store the ascii output files;\n' +
                         '<required>') 
parser.add_argument('-s', '--stations', dest='stations', action='store', 
                    nargs='+', type=str, default=['*.*'],
                    help='stations to be extracted, given as a\n' +  
                         'list of regexes of "Network.Station";\n' + 
                         'default = *.* (all stations)')
parser.add_argument('-c', '--channels', dest='channels', action='store', 
                    type=str, default='ENZRTSP',
                    help='channels to be extracted, depending on\n' + 
                         'OUT_STATIONS_COMPONENTS in inparam.time_src_recv;\n' +  
                         'default = ENZRTSP (all channels)')
parser.add_argument('-e', '--numfmt', dest='number_format', action='store', 
                    type=str, default='%.6e',
                    help='output number format;\n' +
                         'default = %%.6e') 
parser.add_argument('-n', '--fnamefmt', dest='fname_format', action='store', 
                    type=str, default='@NW@.@ST@..BH@CH@.ascii',
                    help='output filename format;\n' +  
                         'default = @NW@.@ST@..BH@CH@.ascii')
parser.add_argument('-H', '--headerfmt', dest='header_format', action='store', 
                    type=str, default='',
                    help='format of the header to be placed\n' + 
                         'at the beginning of each ascii file;\n' + 
                         'default = ""')
parser.add_argument('-F', '--footerfmt', dest='footer_format', action='store', 
                    type=str, default='',
                    help='format of the footer to be placed\n' +  
                         'at the end of each ascii file;\n' + 
                         'default = ""')
parser.add_argument('-p', '--nproc', dest='nproc', action='store', 
                    type=int, default=1, 
                    help='number of processors; default = 1')                         
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                    help='verbose mode')                         
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset
import fnmatch, os, shutil
from multiprocessing import Pool

# open file
ncdf = Dataset(args.in_nc_file, 'r')

# time info
vartime = ncdf.variables['time_points'][:]
nsteps = len(vartime)
strnsteps = str(nsteps)
assert nsteps > 0, 'Zero time steps'
if nsteps == 1:
    dt = 'n.a.'
    sampling_rate = 'n.a.'
else:
    dt = str(vartime[1] - vartime[0])
    sampling_rate = str(1. / (vartime[1] - vartime[0]))
tstart = str(vartime[0])
tend = str(vartime[-1])

if args.verbose:
    print('--- Time info ---')
    print('Number of steps: ' + strnsteps)
    print('Time step: ' + dt)
    print('Sampling rate: ' + sampling_rate)
    print('Start time: ' + tstart)
    print('End time: ' + tend)
    print()

# create output directory
try:
    os.makedirs(args.out_ascii_dir)
except OSError:
    pass

# names
fname = args.fname_format.replace('@NS@', strnsteps)
fname = fname.replace('@DT@', dt)
fname = fname.replace('@SR@', sampling_rate)
fname = fname.replace('@T0@', tstart)
fname = fname.replace('@T1@', tend)

header = args.header_format.replace('@NS@', strnsteps)
header = header.replace('@DT@', dt)
header = header.replace('@SR@', sampling_rate)
header = header.replace('@T0@', tstart)
header = header.replace('@T1@', tend)

footer = args.footer_format.replace('@NS@', strnsteps)
footer = footer.replace('@DT@', dt)
footer = footer.replace('@SR@', sampling_rate)
footer = footer.replace('@T0@', tstart)
footer = footer.replace('@T1@', tend)

# stations
stations = []
for var in ncdf.variables:
    if var == 'time_points':
        continue
    
    # station info
    varlist = var.split('.')
    assert len(varlist) == 3
    network = var.split('.')[0]
    station = var.split('.')[1]
    synchannels = var.split('.')[2]
    
    # match station regexes
    key = network + '.' + station
    for regex in args.stations:
        if fnmatch.fnmatch(key, regex):
            stations.append([network, station, var])
            break
    
# close
ncdf.close()    

# extract waveforms
nstations = len(stations)
if args.verbose:
    print('--- Extracting Waveforms at %d Stations ---' % (nstations))
    
def write_ascii(iproc):
    if args.nproc == 1:
        ncdf = Dataset(args.in_nc_file, 'r')
        iproc = 0
    else:
        # copy netcdf file for parallel access
        tempnc = args.out_ascii_dir + '/nc_temp.nc' + str(iproc)
        shutil.copy(args.in_nc_file, tempnc)
        ncdf = Dataset(tempnc, 'r')

    for ist, st in enumerate(stations):
        if ist % args.nproc != iproc: 
            continue
        
        # station info
        network = st[0]
        station = st[1]
        var = st[2]
        
        # names
        fname_st = fname.replace('@NW@', network)
        fname_st = fname_st.replace('@ST@', station)
        header_st = header.replace('@NW@', network)
        header_st = header_st.replace('@ST@', station)
        footer_st = footer.replace('@NW@', network)
        footer_st = footer_st.replace('@ST@', station)

        # data
        wave = ncdf.variables[var][:, :]
        
        # waveform
        for channel in args.channels:
            dim = synchannels.find(channel)
            if dim < 0:
                continue
            # names
            fname_stc = fname_st.replace('@CH@', channel)
            header_stc = header_st.replace('@CH@', channel)
            footer_stc = footer_st.replace('@CH@', channel)

            # save to ascii
            fname_stc = args.out_ascii_dir + '/' + fname_stc
            np.savetxt(fname_stc, wave[:, dim], 
                fmt=args.number_format, 
                header=header_stc, footer=footer_stc, comments='')
            if args.verbose:
                print('Done with ascii file = %s; no. station = %d / %d; iproc = %d' \
                    % (fname_stc, ist + 1, nstations, iproc))

    # close
    ncdf.close()
    # remove temp nc
    if args.nproc > 1:
        os.remove(tempnc)
    
# write_ascii in parallel
args.nproc = max(args.nproc, 1)
if args.nproc == 1:
    write_ascii(0)
else:
    with Pool(args.nproc) as p:
        p.map(write_ascii, range(0, args.nproc))
        


