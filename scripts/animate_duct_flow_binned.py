#!/usr/bin/python

import numpy as np
import scipy.stats
from pylab import sys,transpose,size,cla,hold

import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

#
# For simulations with position histories, 
# bin the 2d data.
#

walks = h5py.File(sys.argv[1])

nwalkers = int(walks['nTrials'].value)
#tsteps = int(walks['timesteps'].value)

X = transpose(walks['X'].value)[:,:nwalkers]
Y = transpose(walks['Y'].value)[:,:nwalkers]
Z = transpose(walks['Z'].value)[:,:nwalkers]

t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value
aratio = walks['aratio'].value

walks.close()

# Make this odd to surround (y,z) = (0,0),
# make it even to split it.

nb = 21
ya = np.linspace(-1.,1.,nb+1)
za = np.linspace(-1.,1.,nb+1)

tsteps = np.size(t)

mean_sl = np.zeros((tsteps,nb,nb))
var_sl = np.zeros((tsteps,nb,nb))
sk_sl = np.zeros((tsteps,nb,nb))
kurt_sl = np.zeros((tsteps,nb,nb))


# Main loop. Get the statistics in y,z bins for all time values.
print "Binning data..."
for kt in range(tsteps):
     for i in range(nb):
          for j in range(nb):
               # Bin the data. Minor annoyance of handling the 
               # last bin.
               if (i == nb-1):
                    masky = ((Y[kt,:] >= ya[i]) & (Y[kt,:] <= ya[i+1]))
               else:
                    masky = ((Y[kt,:] >= ya[i]) & (Y[kt,:] < ya[i+1]))
               # end if
               
               if (j == nb-1):
                    maskz = ((Z[kt,:] >= za[j]) & (Z[kt,:] <= za[j+1]))
               else:
                    maskz = ((Z[kt,:] >= za[j]) & (Z[kt,:] < za[j+1]))
               # end if
               
               # Create the combined mask, apply it to X, and 
               # calculate the statistics.
               mask = (masky & maskz)
               xsamp = X[kt,mask]

#               _,_,mmean,mvar,mskew,mkurt = scipy.stats.describe(xsamp)
               _,_,_,_,mskew,_ = scipy.stats.describe(xsamp)
               
#               mean_sl[kt,i,j] = mmean
#               var_sl[kt,i,j] = mvar
               sk_sl[kt,i,j] = mskew
#               kurt_sl[kt,i,j] = mkurt
               
          # end for
     # end for
     print float(kt)/tsteps*100
# end for

#
# Look at the statistics at a few sample bins.
# These are indices for the stats.
#

# Center
ci = nb/2
cj = nb/2

# Short-end wall
wi = nb-1
wj = nb/2

# Corner
coi = nb-1
coj = nb-1

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')

ax.hold(True)
# Skewness at the center
ax.plot(t,sk_sl[:,ci,cj],label="Center")
# Skewness at short wall
ax.plot(t,sk_sl[:,wi,wj],label="Wall")
# Skewness at corner
ax.plot(t,sk_sl[:,coi,coj],label="Corner")

ax.hold(False)
pyplot.legend(loc='best')
pyplot.xlabel(r'Time $\tau$')
pyplot.ylabel('Skewness')

#pyplot.show(block=False)
pyplot.savefig('skewness_sl.png',dpi=150)
