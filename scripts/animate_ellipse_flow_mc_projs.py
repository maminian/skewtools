#!/usr/bin/python

import numpy as np
import pylab
from pylab import *
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

# ------------------------------------------------

def update(i,axXY,axXZ,axYZ,axSk,x,y,z,Sk,t,aratio):

     # x-bounds for the axes.
     maxx = max(abs(x[i,:]))
     
     # Change the current axis, clear it, and update it with the new
     # scatter info.
     pyplot.sca(axXY)
     cla()
     axXY.scatter(x[i,:], y[i,:], s=5, marker='.', alpha=0.2)
     pylab.xlim([-maxx,maxx])  # X axis
     pylab.ylim([-1.,1.])       # Y axis
     pyplot.title('XY projection')

     # Repeat for the other two projections.
     pyplot.sca(axXZ)
     cla()
     axXZ.scatter(x[i,:], z[i,:], s=5, marker='.', alpha=0.2)
     pylab.xlim([-maxx,maxx])  # X axis
     pylab.ylim([-1./aratio,1./aratio])       # Z axis
     pyplot.title('XZ projection')
     
     # For the YZ projection, draw the unit circle as well.
     pyplot.sca(axYZ)
     cla()
     axYZ.scatter(y[i,:], z[i,:], s=5, marker='.', alpha=0.2)
     hold(True)
     th = linspace(0,2*pi,361)
     axYZ.plot(cos(th),1./aratio*sin(th),color='red')
     hold(False)
     pylab.xlim([-1.,1.])  # Y axis
     pylab.ylim([-1./aratio,1./aratio])  # Z axis
     pyplot.title('YZ projection')
     
     # For the skewness, the plot stays the same, but we need to 
     # re-draw the plot with a new position for the red line.
     pyplot.sca(axSk)
     cla()
     
     axSk.plot(t,Sk)
     axSk.set_xscale('log')
     axSk.grid(True)
     hold(True)     
     axSk.plot([t[i],t[i]],[min(Sk),max(Sk)],color='red')
     hold(False)
     pyplot.title('Skewness over time')
     
     return axXY,axXZ,axYZ,axSk

# end update

def update_slow(i,axXY,axXZ,axYZ,axSk,walks,Sk,t,aratio):

     print i
     x = walks['X'][:,i]
     y = walks['Y'][:,i]
     z = walks['Z'][:,i]
     
     # x-bounds for the axes.
#     maxx = max(abs(x).max(),maxx)
     
     # Change the current axis, clear it, and update it with the new
     # scatter info.
     pyplot.sca(axXY)
     cla()
     axXY.hist2d(x, y,bins=201)
#     pylab.xlim([-maxx,maxx])  # X axis
     pylab.ylim([-1.,1.])       # Y axis
     pyplot.title('XY projection')

     # Repeat for the other two projections.
     pyplot.sca(axXZ)
     cla()
     axXZ.hist2d(x, z,bins=201)
#     pylab.xlim([-maxx,maxx])  # X axis
     pylab.ylim([-1./aratio,1./aratio])       # Z axis
     pyplot.title('XZ projection')
     
     # For the YZ projection, draw the unit circle as well.
     pyplot.sca(axYZ)
     cla()
     axYZ.hist2d(y, z,bins=201)
     hold(True)
     th = linspace(0,2*pi,361)
     axYZ.plot(cos(th),1./aratio*sin(th),color='red')
     hold(False)
     pylab.xlim([-1.,1.])  # Y axis
     pylab.ylim([-1./aratio,1./aratio])  # Z axis
     pyplot.title('YZ projection')
     
     # For the skewness, the plot stays the same, but we need to 
     # re-draw the plot with a new position for the red line.
     pyplot.sca(axSk)
     cla()
     
     axSk.plot(t,Sk)
     axSk.set_xscale('log')
     axSk.grid(True)
     hold(True)     
     axSk.plot([t[i],t[i]],[min(Sk),max(Sk)],color='red')
     hold(False)
     pyplot.title('Skewness over time')
     
     return axXY,axXZ,axYZ,axSk

# end update_slow

def animate_mc(x,y,z,Sk,t,aratio):

     fig = pyplot.figure()
     # Generate a subplot with the particle animation on top
     # and an animated tracker of the skewness over time
     # on the bottom.

     gs = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[1,1])
     axXY = fig.add_subplot(gs[0])
     axXZ = fig.add_subplot(gs[1])
     axYZ = fig.add_subplot(gs[2])
     axSk = fig.add_subplot(gs[3])
     
     axXY.scatter(x[0,:], y[0,:], s=5, marker='.',alpha=0.2)
     axXZ.scatter(x[0,:], z[0,:], s=5, marker='.',alpha=0.2)
     axYZ.scatter(y[0,:], z[0,:], s=5, marker='.',alpha=0.2)
     th = linspace(0,2*pi,361)
     axYZ.plot(cos(th),1./aratio*sin(th),color='red')

     axSk.plot(t,Sk)
     
     
#     delay = 1000./(fps)
     ani = anim.FuncAnimation( fig, update, frames=size(t), 
                                   fargs=(axXY,axXZ,axYZ,axSk,x,y,z,Sk,t,aratio),repeat=False)
     
     return ani
     
# --------------------------------------------

def animate_mc_slow(walks):

     fig = pyplot.figure()
     # Generate a subplot with the particle animation on top
     # and an animated tracker of the skewness over time
     # on the bottom.

     t = walks['Time'].value
     Sk = walks['Avgd_Skewness'].value

     aratio = walks['aratio'].value

     gs = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[1,1])
     axXY = fig.add_subplot(gs[0])
     axXZ = fig.add_subplot(gs[1])
     axYZ = fig.add_subplot(gs[2])
     axSk = fig.add_subplot(gs[3])
     
     axXY.scatter(walks['X'][:,0], walks['Y'][:,0], s=5, marker='.',alpha=0.2)
     axXZ.scatter(walks['X'][:,0], walks['Z'][:,0], s=5, marker='.',alpha=0.2)
     axYZ.scatter(walks['Y'][:,0], walks['Z'][:,0], s=5, marker='.',alpha=0.2)
     th = linspace(0,2*pi,361)
     axYZ.plot(cos(th),1./aratio*sin(th),color='red')

     axSk.plot(t,Sk)
     
#     delay = 1000./(fps)
#     ani = anim.FuncAnimation( fig, update, frames=size(t), 
#                                   fargs=(axXY,axXZ,axYZ,axSk,x,y,z,Sk,t,aratio),repeat=False)


     ani = anim.FuncAnimation( fig, update_slow, frames=size(t), 
                                   fargs=(axXY,axXZ,axYZ,axSk,walks,Sk,t,aratio),repeat=False)

     return ani
     
# --------------------------------------------

walks = h5py.File(sys.argv[1])

print "Constructing animation..."

# Need to split this into two versions, depending on whether I have the memory capacity 
# to store the entire position histories.
if False:
     X = transpose(walks['X'].value)
     Y = transpose(walks['Y'].value)
     Z = transpose(walks['Z'].value)
     Sk = transpose(walks['Avgd_Skewness'].value)
     t = transpose(walks['Time'].value)

     aratio = walks['aratio'].value

     ani = animate_mc(X,Y,Z,Sk,t,aratio)
else:
     ani = animate_mc_slow(walks)
# end if


# Show the animation, or write it to a movie file.
if True:
     print "Saving animation to disk (may take quite a while)..."
     #dpi=200 gives decent quality for the filesize. Takes about 5-10 minutes to make.
     ani.save(sys.argv[2],dpi=120,writer=anim.FFMpegWriter(fps=24))

else:
     pyplot.show(block=False)
#

#walks.close()

print "Done."


