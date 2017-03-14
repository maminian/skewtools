#!/usr/bin/python

import numpy as np

from numpy import transpose,size,sqrt
import sys

import matplotlib.pyplot as pyplot
from matplotlib.pyplot import cla,hold

import matplotlib.animation as anim
import h5py

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

import string

# ------------------------------------------------

def update(i,axXY,axXZ,axYZ,axAvSk,x,y,z,AvSk,Pe,t,mmap):

     # x-bounds for the axes.
     maxx = max(abs(x[i,:]))
     
     # Change the current axis, clear it, and update it with the new
     # scatter info.
     pyplot.sca(axXY)
     cla()

     axXY.hist2d(x[i,:], y[i,:],bins=201,cmap=mmap,vmin=1,vmax=100)
#     pylab.xlim([-maxx,maxx])  # X axis
     axXY.set_ylim([-1.,1.])       # Y axis
     pyplot.title('XY projection')

     # Repeat for the other two projections.
     pyplot.sca(axXZ)
     cla()

     axXZ.hist2d(x[i,:], z[i,:],bins=201,cmap=mmap,vmin=1,vmax=100)
     axXZ.set_xlim([-maxx,maxx])  # X axis
     axXZ.set_ylim([-1.,1.])       # Z axis
     pyplot.title('XZ projection')
     
     pyplot.sca(axYZ)
     cla()

     axYZ.hist2d(y[i,:], z[i,:],bins=201,cmap=mmap,vmin=1,vmax=100)
     axYZ.set_xlim([-1.,1.])  # Y axis
     axYZ.set_ylim([-1.,1.])  # Z axis
     pyplot.title('YZ projection')
     
     # For the skewness, the plot stays the same, but we need to 
     # re-draw the plot with a new position for the red line.
     pyplot.sca(axAvSk)
     cla()
     
     axAvSk.plot(t,AvSk)
     axAvSk.set_xscale('log')
     axAvSk.grid(True)
     hold(True)     
     axAvSk.plot([t[i],t[i]],[min(AvSk),max(AvSk)],color='red')
     hold(False)
     pyplot.title('Skewness over time')
     
     return axXY,axXZ,axYZ,axAvSk

#def animate_mc(x,y,z,Sk,t):
def animate_mc(x,y,z,AvSk,Mean,Pe,t):

     fig = pyplot.figure()
     # Generate a subplot with the particle animation on top
     # and an animated tracker of the skewness over time
     # on the bottom.

     gs = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[1,1])
     axXY = fig.add_subplot(gs[0],axisbg='white')
     axXZ = fig.add_subplot(gs[1],axisbg='white')
     axYZ = fig.add_subplot(gs[2],axisbg='white')
     axAvSk = fig.add_subplot(gs[3])
     
     mmap = pyplot.cm.Greens
     mmap.set_under('white')
     
#     axXY.hist2d(x[0,:], y[0,:], bins=201, norm=LogNorm(),alpha=0.2)
#     axXZ.hist2d(x[0,:], z[0,:], bins=201, norm=LogNorm(),alpha=0.2)
#     axYZ.hist2d(y[0,:], z[0,:], bins=201, norm=LogNorm(),alpha=0.2)
     axXY.hist2d(x[0,:], y[0,:], bins=201, cmap=mmap, vmin=1, vmax=100)
     axXZ.hist2d(x[0,:], z[0,:], bins=201, cmap=mmap, vmin=1, vmax=100)
     axYZ.hist2d(y[0,:], z[0,:], bins=201, cmap=mmap, vmin=1, vmax=100)

     axAvSk.plot(t,AvSk)
     
#     delay = 1000./(fps)
     ani = anim.FuncAnimation( fig, update, frames=size(t), 
                                   fargs=(axXY,axXZ,axYZ,axAvSk,x,y,z,AvSk,Pe,t,mmap),repeat=False)
     
     return ani
     
# --------------------------------------------

# Hooray hdf
walks = h5py.File(sys.argv[1])

print "Constructing animation..."
nwalkers = int(walks['nTrials'].value)
tsteps = int(walks['timesteps'].value)

X = transpose(walks['X'].value)[:,:nwalkers]
Y = transpose(walks['Y'].value)[:,:nwalkers]
Z = transpose(walks['Z'].value)[:,:nwalkers]

AvSk = transpose(walks['Avgd_Skewness'].value)[:nwalkers]
AvVar = transpose(walks['Avgd_Variance'].value)[:nwalkers]

XMean = transpose(walks['Mean'].value)

t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value

walks.close()

ani = animate_mc(X,Y,Z,AvSk,XMean,Pe,t)

# Show the animation, or write it to a movie file.
if False:
     print "Saving animation to disk (may take quite a while)..."
     #
     # dpi=200 gives decent quality for the filesize. Takes about 5-10 minutes to make.
     
     fname = sys.argv[2]
     ani.save(fname+'.mp4',dpi=200,fps=30,writer=anim.FFMpegWriter(fps=24))
else:
     pyplot.show(block=False)

print "Done."


