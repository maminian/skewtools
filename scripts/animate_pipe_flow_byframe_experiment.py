#!/usr/bin/python

import numpy as np
import sys
from numpy import *

import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

import string

# ------------------------------------------------

def generate_radial_pdf(rbins):
     # Given a set of radial bin boundaries, produce the 
     # probability of a uniformly distributed particle landing 
     # in each of the annuli. 
     
     nbins = len(rbins)-1
     rpdf = zeros(nbins)
     
     for i in range(nbins):
          rpdf[i] = rbins[i+1]**2 - rbins[i]**2
     # end for
     rpdf = rpdf/sum(rpdf)

     return rpdf
# end def

def apply_radial_weights(ri,rpdf,rbins):
     # Given radial locations ri, 
     # generate weights for them inversely proportional 
     # to the area of the bin that they land in.
     
     weights = zeros(len(ri)*2)
     nbins = len(rbins)-1
     
     for i in range(len(ri)):
          for k in range(nbins):
               # Can take a small shortcut here since the bins are in increasing order.
               if ( (ri[i] <= rbins[k+1]) & (ri[i] >= rbins[k]) ):
                    weights[i] = 1./rpdf[k]
                    weights[i+len(ri)] = 1./rpdf[k]
                    break
               # end if
          # end for

     # end for
     
     return weights
# end def


def rescale_coord(x,Pe,t):
     # Transform x coordinates from local back to lab coordinates.
     umean = 0.25    # Nondimensionalized lab frame average flow speed.
     return (x + Pe*umean*t)
# end def

def snapshot(i,axI,axISl,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer):

     # x-bounds for the axes.

     xlims = [0.,Lscale]
     
     timer.set_text(r'$\tau = %3.2e$'%t[i])
     
     # Keep a smaller number of bins for the transverse direction.
     binref=21
     mybins=[linspace(xlims[0],xlims[1],int(binref*Lscale)),linspace(-1.,1.,binref)]
     slsidx = binref/2 + np.floor(sls*binref/2)

     # Change the current axis, clear it, and update it with the new
     # scatter info.
     pyplot.sca(axI)
     pyplot.cla()
     
     # Plot "intensity" and save the bins for the following plot.
     Islices,xloc,yloc,_ = axI.hist2d(xi,yi,bins=mybins,cmap=mmap,normed=True)
     axI.hold(True)
     
     # Plot the locations of the slice lines (as well as computing them.
     for sl in sls:
          axI.plot(xlims,[sl,sl],color='white',linewidth=0.5)
     # end for
     axI.hold(False)

     # Set the window, if it hasn't been already.
     
     axI.set_xlim(xlims)           # X axis
     axI.set_ylim([-1.,1.])       # Y axis
     
     realtime = tscale*t[i]

     pyplot.title('Intensity')

     # 
     pyplot.sca(axISl)
     pyplot.cla()
     axISl.hold(True)
     # Plot the previous hist2d cut along the specified slices.

     # COLORS TO MATCH SARAH'S PLOT
     slcolors=["#7AFF04","#0C04FF","#790776"]

     for i in range(len(sls)):
          idx = int(slsidx[i])
          axISl.plot(xloc[:-1],Islices[:,idx], color=slcolors[i], label='y = %3.2f a' % (sls[i]))
     # end for
     axISl.hold(False)
     axISl.legend(loc='upper right')
     pyplot.title('Intensity on vertical slices')
     axISl.grid(True)

     # Now we need to work out the concentration density.     
     pyplot.sca(axC)
     pyplot.cla()
     axC.hold(True)
     
     ri = sqrt(yi**2 + zi**2)
     myweights = apply_radial_weights(ri,rpdf,linspace(0.,1.,binref/2))
     
     Cslices,xloc,yloc,_ = axC.hist2d(np.concatenate([xi,xi]),np.concatenate([ri,-ri]),bins=mybins,weights=myweights,cmap=mmap,normed=True)

     # Plot the locations of the slice lines.
     for sl in sls:
          axC.plot(xlims,[sl,sl],color='white',linewidth=0.5)
     # end for
     axC.hold(False)
     pyplot.title('Concentration')
     
     # 
     pyplot.sca(axCSl)
     pyplot.cla()
     axCSl.hold(True)
     
     # Plot the previous hist2d cut along the specified slices.
     for i in range(len(sls)):
          idx = int(slsidx[i])
          axCSl.plot(xloc[:-1],Cslices[:,idx], color=slcolors[i], label='r = %3.2f a' % (sls[i]))
     # end for
     axCSl.hold(False)
     axCSl.legend(loc='upper right')
     pyplot.title('Concentration on radial slices')
     pyplot.xlabel(r'$x/a$')
     axCSl.grid(True)

     return axI,axISl,axC,axCSl
     
# end snapshot

# Hooray hdf
print "Reading preliminaries..."
walks = h5py.File(sys.argv[1])

nwalkers = int(walks['nTrials'].value)
tsteps = int(walks['timesteps'].value)

t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value
aratio = walks['aratio'].value

mmap = pyplot.cm.jet
undercolor = pyplot.cm.jet([0.,1.])[0] # Dark blue
mmap.set_under(color=undercolor)
#mmap.set_under('white')

fig = pyplot.figure(figsize=(16,9))

timer=fig.text(0.5,0.975,'Time = %3.2e'%t[0],horizontalalignment='center',verticalalignment='top')
timer.set_fontsize(24)

# Goal is to match Sarah's videos.
gs = gridspec.GridSpec(5,1,height_ratios=[0.2,1,1,1,1])
axI = fig.add_subplot(gs[1],axisbg=undercolor)
axISl = fig.add_subplot(gs[2])
axC = fig.add_subplot(gs[3],axisbg=undercolor)
axCSl = fig.add_subplot(gs[4])

# The pipe in the experiment is 5*10**-4m in radius. She 
# takes slices at y=r=0.0001159, y=r=0.0002507, and y=r=0.0003585, 
# relative to origin at the center of the pipe.
# Pipe length is 0.3m, but the video only looks at the first ~0.1m of pipe.
R = 5.*10**-4
sls = np.array([0.0001159,0.0002507,0.0003585])/R
L = 0.1
Lscale = L/R
kappa = 4.9*10**-6*(10**-4)   # Units m^2/s

tscale = R**2/kappa

# Do the first timestep to set up the sizing of the axes for tight_layout.
# It will be repeated in the loop.
xi = walks['X'][:nwalkers,0]
yi = walks['Y'][:nwalkers,0]
zi = walks['Z'][:nwalkers,0]

binref = 21
rpdf = generate_radial_pdf(linspace(0.,1.,binref/2))


snapshot(0,axI,axISl,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer)

pyplot.tight_layout()

print "Constructing animation..."
for i in range(tsteps):
     xi = walks['X'][:nwalkers,i]
     for j in range(nwalkers):
          xi[j] = rescale_coord(xi[j],Pe,t[i])
     # end for
     yi = walks['Y'][:nwalkers,i]
     zi = walks['Z'][:nwalkers,i]
     
     snapshot(i,axI,axISl,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer)

     outfile = 'movie_%s.png'%string.zfill(i,4)     
     print outfile
     pyplot.savefig(outfile,dpi=150)

     
# end for

walks.close()

