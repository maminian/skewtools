# Generate a 4x4 plot of the partial variance and skewness (as a 
# function of y and z), imported from a simulation file.
#
# Usage:
# 
# python moments_2d_panels.py simulationfile.h5
#
# The simulation file is assumed to have the proper binned data.

import h5py
import sys
import matplotlib

import matplotlib.pyplot as pyplot
from numpy import *


mf = h5py.File(sys.argv[1])

mycmap = pyplot.get_cmap('seismic')

# Get the problem geometry;
# Channel=0, Duct=1, Ellipse=2
#geometry = int(mf['geometry'][0])

# Hard-coded for now until the simulations contain the geometry.
geometry = 2

aratio = mf['aratio'][0]
Pe = mf['Peclet'][0]

nby = int(mf['nBinsY'].value)
nbz = int(mf['nBinsZ'].value)
va = mf['Variance'].value
sk = mf['Skewness'].value
t = mf['Time'].value

mf.close()

mxsk = abs(sk).max() # Largest absolute skewness
nt = size(t) # Number of time samples.

npan = 4         # Number of panels in one direction.
npt = min(npan**2,len(t))

fheight = 8    # Figure height in "inches"
fwidth = fheight/aratio + 0.05*npan + 1     # Corresponding figure width to keep physical dimensions.
fig,ax = pyplot.subplots(npan,npan,sharex=True,sharey=True,figsize=(fwidth,fheight))

# This might not cooperate nicely with the duct.
# Also, need to add an extra "bin" for pcolor to cooperate.
Y,Z = meshgrid(linspace(-1.,1.,nby+1),linspace(-1./aratio,1./aratio,nbz+1))

# Need a *second* meshgrid just for the contouring!
Y2,Z2 = meshgrid(linspace(-1.,1.,nby),linspace(-1./aratio,1./aratio,nbz))
contours = linspace(-ceil(mxsk),ceil(mxsk),80)

# Define the mask for the ellipses to remove the unphysical domain.
# Custom colormap via:
# http://stackoverflow.com/questions/9707676/defining-a-discrete-colormap-for-imshow-in-matplotlib
Yref,Zref = meshgrid(linspace(-1.,1.,200),linspace(-1./aratio,1./aratio,200))

if (geometry == 2):
     mmask = (Yref**2 + (Zref*aratio)**2 <= 1.)
     #maskcmap = matplotlib.colors.ListedColormap([(1.,1.,1.,1.),(0.,0.,0.,0.)])
     maskcmap = matplotlib.colors.ListedColormap([(0.,0.,0.,1.),(1.,1.,1.,0.)])
# end if

tsamps = linspace(0,len(t)-1,npt).astype(int)

#for i in range(np):
#     for j in range(np):
for k in range(npt):
     j = mod(k,npan)
     i = (k-j)/npan

     idx = tsamps[k]

     # Apparently "nearest" interpolation actually means "none".
     ax[i,j].hold(True)
     im = ax[i,j].pcolor(Z,Y,sk[:,:,idx],vmin=-mxsk,vmax=mxsk,cmap=mycmap)
#          im = ax[i,j].contourf(Z2,Y2,sk[:,:,idx],levels=contours,vmin=-mxsk,vmax=mxsk,cmap=mycmap)

     if False:
          # Mask out the exterior of the ellipse.
          ax[i,j].pcolor(Zref,Yref,mmask,cmap=maskcmap)
     # end if

#          im = ax[i,j].imshow(sk[:,:,idx],vmin=-mxsk,vmax=mxsk,cmap=mycmap,interpolation='nearest')

     # Mask out the exterior for the pipes.
#          if (geometry==2):          
#               ax[i,j].pcolor(Z2,Y2,mmask,cmap=maskcmap)
     # end if

     ax[i,j].hold(False)
	  
     prettytime=r'$%.2g$'%t[idx]

     ax[i,j].set_title( r'$t=$'+prettytime )

     # Pretty the plots.
     padding = 0.1
     limy = 1 + padding
     limx = 1./aratio + padding
#          limx=limy
#          ax[i,j].set_aspect('equal')
#          ax[i,j].set_xlim([-limx,limx])
#          ax[i,j].set_ylim([-limy,limy])

     ax[i,j].set_xticks([int(-1./aratio),0.,int(1./aratio)])
     ax[i,j].set_yticks([-1.,0.,1.])
     ax[i,j].set_xticklabels([str(int(-1./aratio)),str(0), str(int(1./aratio))])
     ax[i,j].set_yticklabels([str(-1),str(0), str(1)])

     ax[i,j].grid(True)

# end for

# Clear the remaining panels.
for k in range(npt,npan**2):
     j = mod(k,npan)
     i = (k-j)/npan
     fig.delaxes(ax[i,j])
# end for

pyplot.tight_layout()

# Add an overall title.
pyplot.subplots_adjust(top=0.88)
prettytime=r'$%.2g$'%t[idx]
fig.suptitle(r'Partial skewness evolution, $\lambda=%.2g$, $Pe=%.2g$' % (aratio,Pe) )

# Add a single colorbar for the panels.
# Has to be AFTER the pyplot.tight_layout() to play nicely.
# (all the plots share the same "value" span, so the fact that I'm using the 
# ax[npan-1,npan-1] axis doesn't matter.)
cax,kw = matplotlib.colorbar.make_axes([axs for axs in ax.flat])
pyplot.colorbar(im, cax=cax, **kw)

for i in range(len(tsamps)):
     print tsamps[i],t[tsamps[i]]
# end

pyplot.savefig('slices_2d.png',dpi=160,bbox_inches='tight')
pyplot.show(block=False)
