# Construct an arbitrary initial condition that can be 
# initially skewed. This is done by 
# distributing points according to the corresponding histogram.

from scipy.special import erf
from numpy import *
import h5py
from matplotlib import pyplot

aratio = 0.1
#nyz = 1000000
nTot = 100
fname = 'icond.h5'

def skewnormaldist(x,alpha=0.,sigma=1.):
     # Via https://en.wikipedia.org/wiki/Skew_normal_distribution

     delta = alpha/sqrt(1+alpha**2)

     # Choose a shift and scale to make this thing mean zero, variance one.
     omega = sigma/sqrt(1-2.*delta**2/pi)

     xi = -omega*delta*sqrt(2./pi)
     xtil = (x-xi)/omega

     return (2./omega)*(1./sqrt(2.*pi)*exp(-xtil**2/2.))*(0.5*(1. + erf(alpha*xtil/sqrt(2.))))
# end def

fx = linspace(-20.,40.,1001)
fy = skewnormaldist(fx,5.,6.)

fig,ax = pyplot.subplots(1,1)

ax.plot(fx,fy)

mos = array([0.,0.,0.,0.])

mos[0] = trapz(fy,fx)
mos[1] = trapz(fy*fx,fx)
mos[2] = trapz(fy*fx**2,fx) - mos[1]**2
mos[3] = (trapz(fy*fx**3,fx) - 3.*trapz(fy*fx**2,fx)*mos[1] + 2.*mos[1]**3)/mos[2]**1.5

print "Area: " + str(mos[0])
print "Mean: " + str(mos[1])
print "Variance: " + str(mos[2])
print "Skewness: " + str(mos[3])



# Create the data. Start by coarsening the 
# data to a small number of bins.
nbinsx = 101
xcoarse = linspace(fx.min(),fx.max(),nbinsx+1)
ycoarse = skewnormaldist(xcoarse,5.,6.)

#for i in range(len(ycoarse)-1):
#     il = argmin(abs(fx - xcoarse[i]))
#     ir = argmin(abs(fx - xcoarse[i+1]))
     
#     print il,ir
     
#     ycoarse[i] = trapz(fy[il:ir],fx[il:ir])
#end for

fig2,ax2 = pyplot.subplots(1,1)

ax2.step(xcoarse + diff(xcoarse[:2])/2.,ycoarse,where='mid')

# Heuristically round and re-examine the moments.
yidx = 15 #Heuristic as fuck
yround = around(ycoarse/ycoarse[yidx]).astype(int)

nx = sum(yround)
ny = 50
nz = int(ny/aratio)

nTot = nx*ny*nz

xa = zeros(nTot)
ya = zeros(nTot)
za = zeros(nTot)

yv = linspace(-1.,1.,ny)
zv = linspace(-1./aratio,1./aratio,nz)

idx = -1
for j in range(len(yround)-1):
#     for k in range(yround[j]):
     blergh = linspace(xcoarse[j],xcoarse[j+1],yround[j])
     for xval in blergh:
          for m in range(ny):
               for n in range(nz):
                    idx += 1
                    xa[idx] = xval
                    ya[idx] = yv[m]
                    za[idx] = zv[n]
               #
          # end for
     # end for
     #
# end for

fig3,ax3 = pyplot.subplots(3,1)
ax3[0].hist(xa,nbinsx)
ax3[1].hist(ya,ny)
ax3[2].hist(za,nz)

if True:
     # Save to file.
     mf = h5py.File(fname,'w')

     mf.create_dataset('nTot',data=nTot)
     mf.create_dataset('X',data=xa)
     mf.create_dataset('Y',data=ya)
     mf.create_dataset('Z',data=za)

     mf.close()
     print 0
# end if

ax2.grid(True)
pyplot.show(block=False)

