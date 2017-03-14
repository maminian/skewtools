#/usr/bin/python
# 
# Convert a text history into an hdf file.
# Need to modify this on a case by case basis,
# since my output format has changed.
#
# This is expecting a name without the ".txt"
# Output will be the same name but with ".h5"

from pylab import *
from numpy import *
import matplotlib.pyplot as pyplot
import h5py

# --------------


# Read data from txt
fname = sys.argv[1]+'.txt'
fo = open(fname,'r');
raw = fo.readlines();
fo.close();

# Open h5 file for writing 
h5file = sys.argv[1]+'.h5'
fo = h5py.File(h5file,'w')

fdata = {}

fdata["Nterms"] = int(raw[0]);
fdata["Peclet"] = float(raw[1]);
fdata["aratio"] = float(raw[2]);
fdata["dt"] = float(raw[3]);
fdata["tsteps"] = int(raw[4]);



fdata["t"] = zeros( (fdata["tsteps"],1) );
if (size(raw[5].split()[1:]) > 1):
     print 'nope'
     # This (most likely) is a position history file
else:          # Skewness.
     fdata["Sk"] = zeros( (fdata["tsteps"],1) );
     i = 0;
     for line in raw[5:]:
          derp = line.split();
          fdata["t"][i] = float(derp[0]);
          fdata["Sk"][i] = float(derp[1]);
          i += 1;
     #
     #
     
     dsetSk = fo.create_dataset('Skew', shape(transpose(fdata['Sk'])), dtype='d')
     dsett = fo.create_dataset('Time', shape(transpose(fdata['t'])), dtype='d')
     
     dsetSk[:] = transpose(fdata['Sk'])
     dsett[:] = transpose(fdata['t'])
#     dsetSk.attrs["Peclet"] = data['Peclet']
#     dsetSk.attrs['Number_of_walks'] = data['Nwalks']
#     dsetSk.attrs['Aspect_ratio'] = data['aratio']
#     dsetSk.attrs['Timestep'] = data['dt']
     
     fo.close()
#
#
#
