#!/usr/bin/python

import numpy as np
from pylab import *
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py
import time

from matplotlib.colors import LogNorm
from matplotlib import gridspec

# ------------------------------------------------


def update(i,fig,axMe,axVar,axSk,axKu,Me,Var,Sk,Ku,mtext,nbins,Pe,t):

#     fig.clear()
#     gs = gridspec.GridSpec(2,2)
#     axMe = fig.add_subplot(gs[0])
#     axVar = fig.add_subplot(gs[1])
#     axSk = fig.add_subplot(gs[2])
#     axKu = fig.add_subplot(gs[3])
#     mtext.remove()
     mtext.set_text(r'Channel stats on slices; $\tau = %.2e$'%t[i])

#     fig.canvas.draw()

     ya = np.linspace(-1.,1.,nbins)

     pyplot.sca(axMe)
     cla()
     pyplot.step(ya,Me[i,:],where='mid')
#     axMe.set_xticks([])
#     pyplot.xlabel(r'$y$')
     memin = Me[:(i+1),:].min()
     memax = Me[:(i+1),:].max()
     axMe.set_ylim([1.1*memin,1.1*memax])
     
     pyplot.title(r'Mean')
     
#     axMe.set_ylim([-2./3*Pe*t[i],1./3*Pe*t[i]])

     pyplot.sca(axVar)
     cla()
     pyplot.step(ya,Var[i,:]/(2.*t[i]),where='mid')
#     axVar.set_xticks([])
#     pyplot.xlabel(r'$y$')
     merp=max(1.,Var[:(i+1),:].max()/(2.*t[i]))
     axVar.set_ylim([0.,1.1*merp])
          
     pyplot.title(r'Variance/$2 \tau$')
     
#     axVar.set_ylim([2.*t[i]])
     
     pyplot.sca(axSk)
     cla()
     pyplot.step(ya,Sk[i,:],where='mid')
#     axSk.set_xticks([])
     pyplot.xlabel(r'$y$')
     pyplot.title(r'Skewness')
     
     axSk.set_ylim([-3.,2.])
     
     pyplot.sca(axKu)
     cla()
     pyplot.step(ya,Ku[i,:],where='mid')
     pyplot.title(r'Kurtosis')
     
     pyplot.xlabel(r'$y$')
     
     axKu.set_ylim([-2.,5.])
     
     return axMe,axVar,axSk,axKu

def animate_stats(Me,Var,Sk,Ku,nby,nbz,Pe,t):

     # We're going to have four plots, so need to create space 
     # accordingly.
     fig = pyplot.figure(figsize=(12,8))

     gs = gridspec.GridSpec(2,2)
     axMe = fig.add_subplot(gs[0])
     axVar = fig.add_subplot(gs[1])
     axSk = fig.add_subplot(gs[2])
     axKu = fig.add_subplot(gs[3])
     
	 nt = length(t)
     # Collect the data for the four sample points:
	 Mes = np.zeros((4,nt))
	 Vas = np.zeros((4,nt))
	 Sks = np.zeros((4,nt))
	 Kus = np.zeros((4,nt))
	 
	 # Center of duct
	 ci = nby/2
	 cj = nbz/2
	 Mes[0,:] = Me[ci,cj,:]
	 Vas[0,:] = Var[ci,cj,:]
	 Sks[0,:] = Sk[ci,cj,:]
	 Kus[0,:] = Ku[ci,cj,:]
	 
	 # Short-end wall
	 si = 0
	 sj = nbz/2
	 Mes[1,:] = Me[si,sj,:]
	 Vas[1,:] = Var[si,sj,:]
	 Sks[1,:] = Sk[si,sj,:]
	 Kus[1,:] = Ku[si,sj,:]
	 
	 # Far-end wall
	 fi = nby/2
	 fj = 0
	 Mes[2,:] = Me[fi,fj,:]
	 Vas[2,:] = Var[fi,fj,:]
	 Sks[2,:] = Sk[fi,fj,:]
	 Kus[2,:] = Ku[fi,fj,:]
	 
	 # Corner
	 ki = 0
	 kj = 0
	 Mes[3,:] = Me[ki,kj,:]
	 Vas[3,:] = Var[ki,kj,:]
	 Sks[3,:] = Sk[ki,kj,:]
	 Kus[3,:] = Ku[ki,kj,:]
	 
     
     mtext = fig.suptitle(r'Duct pointwise stats, ; $\tau = %.2e$'%t[0],fontsize=20)
       
     
#     mmap = pyplot.cm.hot
#     mmap.set_under(color='white')
     
#     pyplot.tight_layout()
     
     ani = anim.FuncAnimation( fig, update, frames=size(t), 
                                   fargs=(fig,axMe,axVar,axSk,axKu,Me,Var,Sk,Ku,mtext,nbins,Pe,t),repeat=False)

     return ani
     
# --------------------------------------------

print "Reading data..."

walks = h5py.File(sys.argv[1])

print "Constructing animation..."
nby = int(walks['nBinsY'].value)
nbz = int(walks['nBinsZ'].value)

Me = transpose(walks['Mean'].value)
Var = transpose(walks['Variance'].value)
Sk = transpose(walks['Skewness'].value)
Ku = transpose(walks['Kurtosis'].value)


t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value

walks.close()


ani = animate_stats(Me,Var,Sk,Ku,nby,nbz,Pe,t)

if (True):
     print "Saving animation to disk (may take quite a while)..."
     #
     # dpi=200 gives decent quality for the filesize. Takes about 5-10 minutes to make.
     # assume a file input '***.h5', chop off the suffix and add the '.mp4'.
     fname = sys.argv[2]
     mywriter = anim.ImageMagickFileWriter()
#     mywriter = anim.FFMpegWriter()
     ani.save(fname,dpi=120,fps=12,writer=mywriter)
else:

     pyplot.show(block=False)
# end if

print "Done."


