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

notes = '''String replacement rules for FILENAME_FORMAT, HEADER_FORMAT and FOOTER_FORMAT:
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
                         'default = ENZRT (all channels)')
parser.add_argument('-e', '--numfmt', dest='number_format', action='store', 
                    type=str, default='%.6e',
                    help='output number format;\n' +
                         'default = %%.6e') 
parser.add_argument('-n', '--fnamefmt', dest='filename_format', action='store', 
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
args = parser.parse_args()

################### PARSER ###################

import numpy as np
from netCDF4 import Dataset
import fnmatch, os

# open file
ncdf = Dataset(args.in_nc_file, 'r', format='NETCDF4')

# time info
vartime = ncdf.variables['time_points']
assert len(vartime) >= 2
nsteps = len(vartime)
strnsteps=str(nsteps)
dt = str(vartime[1] - vartime[0])
sampling_rate = str(1. / (vartime[1] - vartime[0]))
tstart = str(vartime[0])
tend = str(vartime[-1])

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

# extract waveforms
print('--- Extracting Waveforms ---')
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
    found = False
    for regex in args.stations:
        if fnmatch.fnmatch(key, regex):
            found = True
            break
    if not found:
        continue
    
    # waveform
    for channel in args.channels:
        dim = synchannels.find(channel)
        if dim < 0:
            continue
        # data
        wave = ncdf[var][0:nsteps, dim]
        # filename
        fname = args.filename_format.replace('@NW@', network)
        fname = fname.replace('@ST@', station)
        fname = fname.replace('@CH@', channel)
        fname = fname.replace('@NS@', strnsteps)
        fname = fname.replace('@DT@', dt)
        fname = fname.replace('@SR@', sampling_rate)
        fname = fname.replace('@T0@', tstart)
        fname = fname.replace('@T1@', tend)
        fname = args.out_ascii_dir + '/' + fname
        # header
        header = args.header_format.replace('@NW@', network)
        header = header.replace('@ST@', station)
        header = header.replace('@CH@', channel)
        header = header.replace('@NS@', strnsteps)
        header = header.replace('@DT@', dt)
        header = header.replace('@SR@', sampling_rate)
        header = header.replace('@T0@', tstart)
        header = header.replace('@T1@', tend)
        # footer
        footer = args.footer_format.replace('@NW@', network)
        footer = footer.replace('@ST@', station)
        footer = footer.replace('@CH@', channel)
        footer = footer.replace('@NS@', strnsteps)
        footer = footer.replace('@DT@', dt)
        footer = footer.replace('@SR@', sampling_rate)
        footer = footer.replace('@T0@', tstart)
        footer = footer.replace('@T1@', tend)
        # save to ascii    
        np.savetxt(fname, wave, fmt=args.number_format, header=header, footer=footer, comments='')
        print('Done with ascii output: ' + fname)

# close
ncdf.close()


