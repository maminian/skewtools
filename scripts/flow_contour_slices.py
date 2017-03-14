#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as pyplot
import h5py
import sys

# hdf version!
# This assumes a particular format of the input file.
u = h5py.File(sys.argv[1])['Flow_values'].value
y = h5py.File(sys.argv[1])['y_mesh'].value
z = h5py.File(sys.argv[1])['z_mesh'].value

aratio = h5py.File(sys.argv[1])['aratio'].value

ulap = h5py.File(sys.argv[1])['Laplacian_values'].value 
#

Y,Z = meshgrid(y,z)

fig,ax=pyplot.subplots(1,2,figsize=(3./aratio+3.,3.))

maxval = abs(u).max()

con = ax[0].pcolor(1./aratio*Z,Y,u,cmap=pyplot.cm.bwr,vmin=-maxval,vmax=maxval)
#ax.set_aspect('equal', adjustable='box')
ax[0].hold(True)
#cbar = pyplot.colorbar(con)

# Plot the slice lines on the left
# y=0, y=0.5, y=0.8 slice lines
yl = [len(y)/2, len(y)*15/20, len(y)*9/10]

slcolors = [[0.,0.,0.8],[1.,0.5,0.],[0.,0.5,0.]]

for i in range(3):
     ax[0].plot([-1./aratio,1./aratio],[y[yl[i]],y[yl[i]]],linewidth=2,color=slcolors[i])
# end for

# Plot the corresponding flow values in the second window.
for i in range(3):
     ax[1].plot(1./aratio*z, u[:,yl[i]],color=slcolors[i],label=r'$y=%.1f$'%y[yl[i]])
# end for


#ax[1].set_aspect('equal', adjustable='box')

ax[0].set_xlim([-1./aratio,1./aratio])
ax[0].set_ylim([-1.,1.])

#ax[1].set_xlim([-1./aratio,1./aratio])
ax[1].set_xlim([-4.,4.]) # Temporary
ax[1].set_ylim([-1.5*maxval,1.5*maxval])

ax[0].set_title(r'Flow profile for rectangle aspect ratio $\lambda = %.3g$'%aratio)
ax[0].set_aspect('equal', adjustable='box')

ax[1].set_title(r'Flow profiles along $y=const.$')
ax[1].legend(loc='upper right',fancybox=True,framealpha=0.5)
ax[1].grid(True)

pyplot.tight_layout()
pyplot.show(block=False)

pyplot.savefig('u_slices_%.2f.png'%aratio,bbox_inches='tight')
