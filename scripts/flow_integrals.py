#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as pyplot
import h5py
import sys

# Define my simple 2D integrator over the rectangular region 
# defined through arrays x and y, based on sampled values
# f(x[i],y[j]) = f[i,j].
#
def my2dint(f,x,y):
     dx = diff(x)
     dy = diff(y)
     
     s = 0.
     for i in range(len(dx)):
          for j in range(len(dy)):
               s += dx[i]*dy[j]*f[i,j]
          # end for
     # end for
     
     return s
# end def

# hdf version!
# This assumes a particular format of the input file.
u = h5py.File(sys.argv[1])['Flow_values'].value
y = h5py.File(sys.argv[1])['y_mesh'].value
z = h5py.File(sys.argv[1])['z_mesh'].value

aratio = h5py.File(sys.argv[1])['aratio'].value

print ""
print "Duct aspect ratio: %.6f." % aratio
print "Mesh size: %i." % len(y)
print "Calculating averages of powers of the flow..."
print "" 

fpowint = zeros(5)

for i in range(len(fpowint)):
     fpowint[i] = my2dint(u**(i+1),y,z)
     area = (y[-1]-y[0])*(z[-1]-z[0])
     print "< u**%i > = %.3e" % (i+1,fpowint[i]/area)
# end for
