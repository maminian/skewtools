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
from matplotlib import gridspec,text

import scripts.skewtools as st

import string

def remove_decorations(axin):
     axin.set_xticks([])
     axin.set_yticks([])

     axin.spines['top'].set_visible(False)
     axin.spines['bottom'].set_visible(False)
     axin.spines['left'].set_visible(False)
     axin.spines['right'].set_visible(False)
#end def


# ------------------------
# Hooray hdf

walks = h5py.File(sys.argv[1])

tsteps = int(walks['timesteps'].value)

AvSk = transpose(walks['Avgd_Skewness'].value)
AvVar = transpose(walks['Avgd_Variance'].value)

# XMean = transpose(walks['Mean'].value)

t = transpose(walks['Time'].value)
Pe = walks['Peclet'][0]
aratio = walks['aratio'][0]

# Generate a subplot with the particle animation on top
# and an animated tracker of the skewness over time
# on the bottom.


mmap = pyplot.cm.inferno
#mmap = pyplot.cm.magma
#mmap = pyplot.cm.plasma
#mmap = pyplot.cm.viridis
bgcolor = mmap([0.,1.])[0]
mmap.set_under(color=bgcolor)

# Set plotting parameters; number of rows/column panels, bins in each direction.
m,n,nb = 3,3,180

X = walks['X']
Y = walks['Y']
Z = walks['Z']

tstep = int(float(len(t))/(m*n-1))
tidxs=range(0,len(t),tstep)
tidxs = tidxs[:m*n]
#tidxs = range(0,len(t),20)
#tidxs.pop(0) # Remove the value at zero.

fig,ax = st.panels_hist2d(t,tidxs,m,n,X,Y,ysc=1.,bins=nb,cmap=mmap,normed=True)
fig2,ax2 = st.panels_hist2d(t,tidxs,m,n,X,Z,ysc=1./aratio,bins=nb,cmap=mmap,normed=True)
fig3,ax3 = st.panels_hist2d(t,tidxs,m,n,Z,Y,xsc=1./aratio,ysc=1.,bins=nb,cmap=mmap,normed=True)

# Modify the plots; reduce spacing and put the titles inside.

fig.subplots_adjust(wspace = 0.05, hspace = 0.05, left=0.05, right=0.95, bottom=0.05, top=0.95)
fig2.subplots_adjust(wspace = 0.05, hspace = 0.05, left=0.05, right=0.95, bottom=0.05, top=0.95)
fig3.subplots_adjust(wspace = 0.05, hspace = 0.05, left=0.05, right=0.95, bottom=0.05, top=0.95)

xpos = 0.1
ypos = 0.7
props = dict(boxstyle='round', facecolor='black', alpha=0.8)

po=0
for k in range(len(tidxs)):

     j = np.mod(k,n)
     i = (k - j)/n
     idx = tidxs[k]

#     print i,j,k
     
     # Re-center the axes about x=0, and draw a dashed black vertical line 
     # at x=0.
     ax[i,j].plot([0,0],[-1,1],color='white',linestyle='--',linewidth=1)
     ax2[i,j].plot([0,0],[-1./aratio,1./aratio],color='white',linestyle='--',linewidth=1)

     xl,xu = ax[i,j].get_xlim()
     xsc = max(abs(xl),abs(xu))
     ax[i,j].set_xlim([-xsc,xsc])
     ax[i,j].set_ylim([-1,1])

     xl,xu = ax2[i,j].get_xlim()
     print xl,xu,1./aratio
     xsc = max(abs(xl),abs(xu))
     ax2[i,j].set_xlim([-xsc,xsc])
#     ax2[i,j].set_ylim([-1./aratio,1./aratio])

     ax3[i,j].set_xlim([-1./aratio,1./aratio])
     ax3[i,j].set_ylim([-1.,1.])

     ax[i,j].set_axis_bgcolor(bgcolor)
     ax2[i,j].set_axis_bgcolor(bgcolor)
     ax3[i,j].set_axis_bgcolor(bgcolor)

     remove_decorations(ax[i,j])
     remove_decorations(ax2[i,j])
     remove_decorations(ax3[i,j])


     ax[i,j].annotate(s=r'$\tau=%.2g$'%t[idx], xy=(xpos,ypos), xycoords='axes fraction', color='white',ha='left',va='center',bbox=props)
     ax2[i,j].annotate(s=r'$\tau=%.2g$'%t[idx], xy=(xpos,ypos), xycoords='axes fraction', color='white',ha='left',va='center',bbox=props)
     ax3[i,j].annotate(s=r'$\tau=%.2g$'%t[idx], xy=(xpos,ypos), xycoords='axes fraction', color='white',ha='left',va='center',bbox=props)
# end for


# Show this stuff.
pyplot.show(block=False)

pyplot.sca(ax[0,0])
pyplot.savefig('hist2d_panels_xy.png',dpi=120,bbox_inches='tight')
pyplot.sca(ax2[0,0])
pyplot.savefig('hist2d_panels_xz.png',dpi=120,bbox_inches='tight')
walks.close()

