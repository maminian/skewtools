#!/usr/bin/python

import numpy as np
from numpy import *
import matplotlib.pyplot as pyplot
from matplotlib.pyplot import cla,hold
import matplotlib.animation as anim
import h5py
import time

import sys

from matplotlib.colors import LogNorm
from matplotlib import gridspec

# ------------------------------------------------

def uc(y):
# Channel flow for which laplace(uc)=-1 and
# the mean of uc over the domain is zero.
     return (1./3 - y**2)
# end def

def plotPiecewiseConstBins(i,ax,nbins,a,b,array):
# Plot piecwise constant functions specified by array.
     
     db = float(b-a)/nbins
    
     ax.hold(True)
     for j in range(nbins):
          horiz = [array[i,j],array[i,j]]
          vert = [a+j*db,a+(j+1)*db]
          
          ax.plot(horiz,vert,color='red')
     # end for
     ax.hold(False)

     return 0
# end def

def update(i,axPos,axSk,x,y,AvSk,Mean,nbins,x0width,Pe,t):

     # Be clever about the axes.
     
     pyplot.sca(axPos)
     cla()

     # Operating under the assumption that u is zero-average.
     idx=i

     xmin = x[0:(idx+1),:].min()*1.1
     xmax = x[0:(idx+1),:].max()*1.1

     ymin = -1.
     ymax = 1.

     if (xmin==xmax):
          xmin = x[1,:].min()*1.1
          xmax = x[1,:].max()*1.1
     # end if

     # Predictively change the axes;
     # lower and upper bounds are set 
     # by the combination of flow and Wiener motion.
     u_min=Pe*uc(1.)
     u_max=Pe*uc(0.)
     stds=6.
     
     #axPos.set_xlim([(u_min*t[idx]-stds*sqrt(t[idx])),u_max*t[idx]+stds*sqrt(t[idx])])
     axPos.set_xlim([xmin,xmax])
     axPos.set_ylim([ymin,ymax])
     
     mmap = pyplot.cm.hot;
     mmap.set_under(color='black')
     
     
#     axPos.scatter(x[i,:], y[i,:], s=5, marker='o', alpha=1.0)
     axPos.hist2d(x[idx,:],y[idx,:],range=[[xmin,xmax],[ymin,ymax]],bins=201,cmap=mmap,normed=True)
     
     
     plotPiecewiseConstBins(i,axPos,nbins,-1.,1.,XMean)
     
     # I think we have to repeat the axis positioning.
#     axPos.set_xlim([(u_min*t[idx]-stds*sqrt(t[idx])),u_max*t[idx]+stds*sqrt(t[idx])])
     axPos.set_xlim([xmin,xmax])
     axPos.set_ylim([ymin,ymax])
     
     axPos.set_xlabel('X')
     axPos.set_ylabel('Y')
     
     xmean = mean(x[idx,:])
     ymean = mean(y[idx,:])
     

     # Plot the flow profile.
#     n_arrows = 11

#     xrange = xmean*ones(n_arrows)
#     yrange = linspace(yminbnd,ymaxbnd,n_arrows)
#     
#     u = uc(yrange)
#     u -= mean(u)
#     v = zeros(n_arrows)
#     
#     winy = linspace(-1.,1.,201)
#     winx = Pe*uc(winy)*t[idx]
#     
#     hold(True)
#     
#     axPos.quiver(xrange,yrange,Pe*u*t[idx],v)
#     axPos.plot([Pe*uc(ymean)*t[idx],Pe*uc(ymean)*t[idx]],[yminbnd,ymaxbnd],color='red')
#     axPos.plot(winx,winy,color='red')

     hold(False)
     
     axPos.set_title('Time = ' + str(t[idx]))

     pyplot.sca(axSk)
     cla()
     
     axSk.plot(t,AvSk)
     axSk.set_xscale('log')
     axSk.grid(True)
     hold(True)     
     axSk.plot([t[idx],t[idx]],[min(AvSk),max(AvSk)],color='red')
     hold(False)
     
#     time.sleep(max(0.,t[idx]-t[idx-1]))
     
     return axPos,axSk

def animate_mc(x,y,AvSk,Mean,nbins,x0width,Pe,t):

     fig = pyplot.figure()
     # Generate a subplot with the particle animation on top
     # and an animated tracker of the skewness over time
     # on the bottom.

     gs = gridspec.GridSpec(2,1, height_ratios=[2,1])
     axPos = fig.add_subplot(gs[0],axisbg='black')
#     axPos.scatter(x[0,:], y[0,:], s=5, marker='.',alpha=0.2)
     mmap = pyplot.cm.hot
     mmap.set_under(color='white')
     axPos.hist2d(x[0,:],y[0,:],bins=201,cmap=mmap,vmin=0.,vmax=100.)
     
     axSk = fig.add_subplot(gs[1])
     
#     hist_stuff = axPos.hist2d(x[1,:],y[1,:],bins=201,norm=LogNorm())
#     im = imshow(hist_stuff[0],cmap=pyplot.cm.Blues,norm=LogNorm())
     
#     cbar=fig.colorbar(im,ax=axPos)
     
     axSk.plot(t,AvSk)
     
     ani = anim.FuncAnimation( fig, update, frames=size(t), 
                                   fargs=(axPos,axSk,x,y,AvSk,Mean,nbins,x0width,Pe,t),repeat=False)

     return ani
     
# --------------------------------------------

print "Reading data..."

walks = h5py.File(sys.argv[1])

print "Constructing animation..."
tsteps = int(walks['timesteps'].value)
nbins = int(walks['nBins'].value)
x0width = walks['x0width'].value

X = transpose(walks['X'].value)
Y = transpose(walks['Y'].value)

AvSk = transpose(walks['Avgd_Skewness'].value)
AvVar = transpose(walks['Avgd_Variance'].value)

XMean = transpose(walks['Mean'].value)
#Sk = transpose(walks['Skewness'].value)

t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value

walks.close()


ani = animate_mc(X,Y,AvSk,XMean,nbins,x0width,Pe,t)

if (False):
     print "Saving animation to disk (may take quite a while)..."
     #
     # dpi=200 gives decent quality for the filesize. Takes about 5-10 minutes to make.
     # assume a file input '***.h5', chop off the suffix and add the '.mp4'.
     fname = sys.argv[2]
     mywriter = anim.ImageMagickFileWriter()
     ani.save(fname,dpi=200,fps=24,writer=mywriter)
else:

     pyplot.show(block=False)
# end if

print "Done."


