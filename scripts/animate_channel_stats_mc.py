#!/usr/bin/python

import sys
import numpy as np
from numpy import *
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
     pyplot.cla()
     pyplot.step(ya,Me[i,:],where='mid')
#     axMe.set_xticks([])
#     pyplot.xlabel(r'$y$')
     memin = Me[:(i+1),:].min()
     memax = Me[:(i+1),:].max()
     axMe.set_ylim([1.1*memin,1.1*memax])
     
     pyplot.title(r'Mean')
     
#     axMe.set_ylim([-2./3*Pe*t[i],1./3*Pe*t[i]])

     pyplot.sca(axVar)
     pyplot.cla()
     pyplot.step(ya,Var[i,:]/(2.*t[i]),where='mid')
#     axVar.set_xticks([])
#     pyplot.xlabel(r'$y$')
     merp=max(1.,Var[:(i+1),:].max()/(2.*t[i]))
     axVar.set_ylim([0.,1.1*merp])
          
     pyplot.title(r'Variance/$2 \tau$')
     
#     axVar.set_ylim([2.*t[i]])
     
     pyplot.sca(axSk)
     pyplot.cla()
     pyplot.step(ya,Sk[i,:],where='mid')
#     axSk.set_xticks([])
     pyplot.xlabel(r'$y$')
     pyplot.title(r'Skewness')
     
     axSk.set_ylim([-3.,2.])
     
     pyplot.sca(axKu)
     pyplot.cla()
     pyplot.step(ya,Ku[i,:],where='mid')
     pyplot.title(r'Kurtosis')
     
     pyplot.xlabel(r'$y$')
     
     axKu.set_ylim([-2.,5.])
     
     return axMe,axVar,axSk,axKu

def animate_stats(Me,Var,Sk,Ku,nbins,Pe,t):

     # We're going to have four plots, so need to create space 
     # accordingly.
     fig = pyplot.figure(figsize=(12,8))

     gs = gridspec.GridSpec(2,2)
     axMe = fig.add_subplot(gs[0])
     axVar = fig.add_subplot(gs[1])
     axSk = fig.add_subplot(gs[2])
     axKu = fig.add_subplot(gs[3])
     
     
     
     mtext = fig.suptitle(r'Channel stats on slices, ; $\tau = %.2e$'%t[0],fontsize=20)
       
     
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
nbins = int(walks['nBins'].value)

Me = transpose(walks['Mean'].value)
Var = transpose(walks['Variance'].value)
Sk = transpose(walks['Skewness'].value)
Ku = transpose(walks['Kurtosis'].value)


t = transpose(walks['Time'].value)
Pe = walks['Peclet'].value

walks.close()


ani = animate_stats(Me,Var,Sk,Ku,nbins,Pe,t)

if (True):
     print "Saving animation to disk (may take quite a while)..."
     #
     # dpi=200 gives decent quality for the filesize. Takes about 5-10 minutes to make.
     # assume a file input '***.h5', chop off the suffix and add the '.mp4'.
     fname = sys.argv[2]
#     mywriter = anim.ImageMagickFileWriter()
     mywriter = anim.FFMpegWriter()
     ani.save(fname,dpi=120,fps=12,writer=mywriter)
else:

     pyplot.show(block=False)
# end if

print "Done."


