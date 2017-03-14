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

fig = pyplot.figure()
ax = pyplot.gca()

maxval = u.max()
minval = u.min()

vmaxval = maxval
if (minval >0):
     vminval = 0.
     vmaxval = u.max()
else:
     vminval = -abs(u).max()
     vmaxval = -vminval
# end if

# HACK FOR TRIANGLE
vminval = -0.3
vmaxval = u.max()

ax.contour(Z,Y,u-vminval,61,colors='white',vmin=vminval,vmax=vmaxval)
ax.contour(Z,Y,u-vminval,[0.],colors='red',vmin=vminval,vmax=vmaxval)
con=ax.contourf(Z,Y,u,61,cmap=pyplot.cm.viridis,vmin=vminval,vmax=vmaxval)
#ax.set_aspect('equal', adjustable='box')
ax.hold(True)

uell = 0.5-Y**2-(Z/aratio)**2
#ax.contour(Z,Y,uell,colors='black')

cbar = pyplot.colorbar(con)

ulapdiff = log10(abs(ulap + 2.0))

fig2,ax2 = pyplot.subplots(1,1)
con2=ax2.contourf(Z[1:-1,1:-1],Y[1:-1,1:-1],ulapdiff,13,cmap=pyplot.cm.Purples)
cbar2 = pyplot.colorbar(con2)

#ax.set_aspect('equal', adjustable='box')
#ax2.set_aspect('equal', adjustable='box')

ax.set_xlim([z.min(),z.max()])
ax.set_ylim([y.min(),y.max()])

ax2.set_xlim([z.min(),z.max()])
ax2.set_ylim([y.min(),y.max()])

#ax.set_title(r'Flow profile for rectangle aspect ratio $\lambda = %.3g$'%aratio)
ax2.set_title('Log-error of the numerical Laplacian from -2')

pyplot.show(block=False)
