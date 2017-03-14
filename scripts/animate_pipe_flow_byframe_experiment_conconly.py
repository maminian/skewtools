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

def snapshot(i,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer):

     # x-bounds for the axes.

     xlims = [0.,Lscale]
     
     timer.set_text(r'$\tau = %3.2e$'%t[i])
     
     # Keep a smaller number of bins for the transverse direction.
     binref=41
     mybins=[linspace(xlims[0],xlims[1],int(binref*Lscale)),linspace(-1.,1.,binref)]
     slsidx = binref/2 + np.floor(sls*binref/2)


     # COLORS TO MATCH SARAH'S PLOT
     slcolors=["#7AFF04","#0C04FF","#790776"]

     # Now we need to work out the concentration density.     
     pyplot.sca(axC)
     pyplot.cla()
     axC.hold(True)
     
     ri = sqrt(yi**2 + zi**2)
     myweights = apply_radial_weights(ri,rpdf,linspace(0.,1.,binref/2))
     
     Cslices,xloc,yloc,_ = axC.hist2d(np.concatenate([xi,xi]),np.concatenate([ri,-ri]), bins=mybins, weights=myweights, cmap=mmap, normed=True)

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

     # Slice locations (in cm). Copied from the main section of code.
     slicesdim = 100.*np.array([0.0001159,0.0002507,0.0003585])

     # Plot the previous hist2d cut along the specified slices.
     for i in range(len(sls)):
          idx = int(slsidx[i])
#          axCSl.plot(xloc[:-1],Cslices[:,idx], color=slcolors[i], label='r = %3.2f a' % (sls[i]))
          axCSl.plot(xloc[:-1],Cslices[:,idx], color=slcolors[i], label='r = %6.4f cm' % (slicesdim[i]))
     # end for
     axCSl.hold(False)
     axCSl.legend(loc='upper right')
     pyplot.title('Concentration on radial slices')
     pyplot.xlabel(r'$x/a$')

     # Hard-coded, simulation specific!
     axCSl.set_ylim([0.,0.5])
     axCSl.set_ylabel(r'$C(x,r,t)$')

     axCSl.grid(True)

     return axC,axCSl
     
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
gs = gridspec.GridSpec(3,1,height_ratios=[0.2,1,1])
axC = fig.add_subplot(gs[1],axisbg=undercolor)
axCSl = fig.add_subplot(gs[2])

# The pipe in the experiment is 5*10**-4m in radius. She 
# takes slices at y=r=0.0001159, y=r=0.0002507, and y=r=0.0003585, 
# relative to origin at the center of the pipe.
# Pipe length is 0.3m, but the video only looks at the first ~0.1m of pipe.
R = 5.*10**-4
sls = np.array([0.0001159,0.0002507,0.0003585])/R
#sls = np.array([0.0001159,0.0002507,0.0003585])
L = 0.1
#L = 0.05
Lscale = L/R
kappa = 4.9*10**-6*(10**-4)   # Units m^2/s

tscale = R**2/kappa

# Do the first timestep to set up the sizing of the axes for tight_layout.
# It will be repeated in the loop.
xi = walks['X'][:nwalkers,0]
yi = walks['Y'][:nwalkers,0]
zi = walks['Z'][:nwalkers,0]

binref = 41
rpdf = generate_radial_pdf(linspace(0.,1.,binref/2))

snapshot(0,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer)

pyplot.tight_layout()

print "Constructing animation..."
for i in range(tsteps):
     xi = walks['X'][:nwalkers,i]
     for j in range(nwalkers):
          xi[j] = rescale_coord(xi[j],Pe,t[i])
     # end for
     yi = walks['Y'][:nwalkers,i]
     zi = walks['Z'][:nwalkers,i]
     
     snapshot(i,axC,axCSl,xi,yi,zi,Pe,aratio,t,sls,Lscale,tscale,rpdf,mmap,timer)

     outfile = 'movie_%s.png'%string.zfill(i,4)     
     print outfile
     pyplot.savefig(outfile,dpi=150)

     
# end for

walks.close()

