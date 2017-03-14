#!/usr/bin/python

import numpy as np
import pylab

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

def snapshot(i,axXY,axXZ,axPDF,axAvSk,xi,yi,zi,AvVar,AvSk,Pe,aratio,t,mmap):

     # x-bounds for the axes.
#     maxx = max(abs(x[i,:]))
     xlims = [-0.8,0.6]
     
     # Get the x variable scaling for this timestep.
     # Includes both advective and diffusive components.
     scaling = 4*sqrt(AvVar[i])
     if (i == 0):
          scaling = 1.
     # end if


     if (scaling==0.):
          scaling = 4*sqrt(np.var(xi))
     # end if
     
     print sqrt(np.var(xi))

     # Change the current axis, clear it, and update it with the new
     # scatter info.
     pyplot.sca(axXY)
     cla()
     

     
#     axXY.hist2d(xi/scaling,yi,bins=201,cmap=mmap,vmin=1,vmax=100)
     axXY.hist2d(xi,yi,bins=201,cmap=mmap,vmin=1,vmax=100)

#     pylab.xlim([-maxx,maxx])  # X axis
     axXY.set_ylim([-1.,1.])       # Y axis
     pyplot.title('XY projection')

     # Repeat for the other two projections.
     pyplot.sca(axXZ)
     cla()
#     axXZ.hist2d(xi/scaling,zi/aratio,bins=201,cmap=mmap,vmin=1,vmax=100)
     axXZ.hist2d(xi,zi/aratio,bins=201,cmap=mmap,vmin=1,vmax=100)
          
#     pylab.xlim([-maxx,maxx])  # X axis
     axXZ.set_ylim([-1./aratio,1./aratio])       # Z axis
     pyplot.title('XZ projection')

     pyplot.sca(axPDF)
     cla()

#     axPDF.hist2d(yi, zi/aratio,bins=201,cmap=mmap,vmin=1,vmax=100)
#     pylab.xlim([-1.,1.])  # Y axis
#     pylab.ylim([-1./aratio,1./aratio])  # Z axis
#     pyplot.title('YZ projection')

     # Instead of showing the cross-section distribution, which 
     # will be trivial for the majority of cases, look at 
     # the averaged distribution.
     nbins = 101

     mybins = np.linspace(scaling*xlims[0],scaling*xlims[1],nbins)
#     axPDF.hist(xi/scaling,bins=mybins,normed=True,histtype='step',align='mid')
     axPDF.hist(xi,bins=mybins,normed=True,histtype='step',align='mid')
     pyplot.title(r'Probability density $\langle C \rangle$')
     
     # For the skewness, the plot stays the same, but we need to 
     # re-draw the plot with a new position for the red line.
     pyplot.sca(axAvSk)
     cla()
     
     axAvSk.plot(t[1:],AvSk[1:])
     axAvSk.set_xscale('log')
     axAvSk.grid(True)
     hold(True)     
     axAvSk.plot([t[i],t[i]],[min(AvSk),max(AvSk)],color='red',linewidth=3)
     hold(False)
     pyplot.title('Full skewness')
     
     
     # Plot decorations

     xticks_real = np.array(xlims)*scaling
     horiticks = np.array([-0.8,-0.4,0.,0.4])*scaling
     horiticklabels = [r'$-0.8$',r'$-0.4$',r'$0$',r'$0.4$']

     axXY.set_xlabel(r'$X / 4 \sigma(\tau)$')
     axXY.set_xlim(xticks_real)
     axXY.set_xticks(horiticks)
     axXY.set_xticklabels(horiticklabels)
          
     axXZ.set_xlabel(r'$X / 4 \sigma(\tau)$')
     axXZ.set_xlim(xticks_real)
     axXZ.set_xticks(horiticks)
     axXZ.set_xticklabels(horiticklabels)
     
     axPDF.set_xlabel(r'$X / 4 \, \sigma(\tau)$')
     axPDF.set_xlim(xticks_real)
     axPDF.set_xticks(horiticks)
     axPDF.set_xticklabels(horiticklabels)
     
     axAvSk.set_xlabel(r'$\tau$')
     axAvSk.set_xlim([t[1],t[-1]])
     return axXY,axXZ,axPDF,axAvSk
     
# end snapshot

# Hooray hdf
print "Reading preliminaries..."
walks = h5py.File(sys.argv[1])


tsteps = int(walks['timesteps'].value)

AvSk = transpose(walks['Avgd_Skewness'].value)
AvVar = transpose(walks['Avgd_Variance'].value)

# XMean = transpose(walks['Mean'].value)

t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value
aratio = walks['aratio'].value

# Generate a subplot with the particle animation on top
# and an animated tracker of the skewness over time
# on the bottom.


mmap = pyplot.cm.inferno
undercolor = pyplot.cm.inferno([0.,1.])[0]
mmap.set_under(color=undercolor)
#mmap.set_under('white')

fig = pyplot.figure()
gs = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[1,1])
axXY = fig.add_subplot(gs[0],axisbg=undercolor)
axXZ = fig.add_subplot(gs[1],axisbg=undercolor)
axPDF = fig.add_subplot(gs[2])
axAvSk = fig.add_subplot(gs[3])

# Do the first timestep to set up the sizing of the axes for tight_layout.
# It will be repeated in the loop.
xi = walks['X'][0]
yi = walks['Y'][0]
zi = walks['Z'][0]
snapshot(0,axXY,axXZ,axPDF,axAvSk,xi,yi,zi,AvVar,AvSk,Pe,aratio,t,mmap)

pyplot.tight_layout()

print "Constructing animation..."
for i in range(tsteps):

     axAvSk.plot(t,AvSk)

     xi = walks['X'][i]
     yi = walks['Y'][i]
     zi = walks['Z'][i]
     
     snapshot(i,axXY,axXZ,axPDF,axAvSk,xi,yi,zi,AvVar,AvSk,Pe,aratio,t,mmap)

     outfile = 'movie_%s.png'%string.zfill(i,4)     
     print outfile
     pyplot.savefig(outfile,dpi=150)

     
# end for

walks.close()

