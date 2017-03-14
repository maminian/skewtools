# Compare Francesca's formulas to a channel simulation.

from numpy import *
import h5py
from matplotlib import pyplot
import sys

def c1(y,t,nmax,Pe):

     p1 = y**4/12. - y**2/6. + 7./180.

     c1sum = Pe*p1

     for n in range(1,nmax+1):
          p2 = 4.*( (-1)**n/(n*pi)**4 )*exp(-(n*pi)**2*t)
          c1sum += Pe*p2*cos(n*pi*y)
     # end for

     return c1sum
# end def

def c2(y,t,nmax,Pe):
     q1 = 1./226800.*(-413. + 3840.*t - 1020.*y**2 + 3570.*y**4 - 2940.*y**6 + 675.*y**8)

     c2sum = Pe**2*q1
     c2sum += 2.*t

     for n in range(1,nmax+1):
          q2 = 2.*(-1)**n/(n*pi)**6*exp(-(n*pi)**2*t)*(17./3. - 64./(n*pi)**2 - 2.*t + y**2)
          q3 = 4.*y*(-1)**n/(n*pi)**5*exp(-(n*pi)**2*t)*( -1./(n*pi)**2 + 1./3.*(y**2-1.))

          c2sum += Pe**2*q2*cos(n*pi*y) + Pe**2*q3*sin(n*pi*y)
     # end for

     return c2sum
# end def

def c3(y,t,nmax,Pe):
     r1 = 1./30.*t*(7. - 30.*y**2 + 15.*y**4)
     r5 = -4076777./13621608000. + 8447./4989600.*y**2 - 713./907200.*y**4 - 1./1200.*y**6 
     r5 += 13./12096.*y**8 - 211./453600.*y**10 + 1./14784.*y**12 + 244./155925.*t - 8./945.*t*y**2 + 4./945.*t*y**4
     
     c3sum = Pe*r1 + Pe**3*r5
     
     for n in range(1,nmax+1):

          expo = exp(-(n*pi)**2*t)

          # ------------------
          r2 = 24.*(-1)**n/(n*pi)**4*expo*t
          # ------------------


          # ------------------
          r3 = 14640./(n*pi)**12 - 3705./(2.*(n*pi)**10) + 691./(20.*(n*pi)**8) + 3.*t**2/(n*pi)**8 
          r3 += y**2*( -231./(2.*(n*pi)**10) + 9./(2.*(n*pi)**8) - 1./(3.*(n*pi)**6) ) 
          r3 += y**4*( 23./(4.*(n*pi)**8) + 2./(3.*(n*pi)**6)) - y**6/(3.*(n*pi)**6)
          r3 += t*(231./(n*pi)**10 - 31./(n*pi)**8 - 8./(15.*(n*pi)**6) - 3.*y**2/(n*pi)**8) 

          r3 *= (-1)**n*expo
          # ------------------


          # ------------------
          r4 = 231./(n*pi)**8 + 44./(n*pi)**6 - 28./(5.*(n*pi)**4)
          r4 += y**2*( -76./(n*pi)**6 + 4./(n*pi)**4 ) 
          r4 += 8.*y**4/(5.*(n*pi)**4) + t*( 6./(n*pi)**6 + 2./(n*pi)**4 - 2.*y**2/(n*pi)**4 )


          r4 *= (-1)**n*expo/(n*pi)**3*y 
          # ------------------


          c3sum += Pe*r2*cos(n*pi*y)
          c3sum += Pe**3*(r3*cos(n*pi*y) + r4*sin(n*pi*y))
     # end for

     return c3sum
# end def

def sk_exact(y,t,nmax,Pe):
     c1e = c1(y,t,nmax,Pe)
     c2e = c2(y,t,nmax,Pe)
     c3e = c3(y,t,nmax,Pe)

     numer = c3e - 3.*c1e*c2e + 2.*c1e**3
     denom = (c2e - c1e**2)**1.5

     return numer/denom
# end def

# Specify parameters
nmax = 100

# Read simulation data.

sim = h5py.File(sys.argv[1])

tsteps = int(sim['timesteps'].value)
nbins = int(sim['nBins'].value)

sk = transpose(sim['Skewness'].value)

mean = transpose(sim['Mean'].value)

t = transpose(sim['Time'].value)
Pe = sim['Peclet'].value

sim.close()

# Define y values to evaluate the exact moments, 
# and the bin centers for the simulation moments.
binshift = 1./nbins
bincenters = linspace(-1.,1.,nbins+1)[:-1] + binshift
yv = linspace(-1.,1.,10*nbins+1)

# Set up a 3x3 grid of plots
ng = 4
fig,ax = pyplot.subplots(ng,ng)

if True:
     for i in range(ng):
          for j in range(ng):
               
               # Evaluate the exact skewness.
               idx = int( (j+ng*i)/float(ng**2)*(size(t)) )
               idx = max(idx,1)
               ske = sk_exact(yv,t[idx],nmax,Pe)


               ax[i,j].hold(True)
               ax[i,j].plot(yv,ske,'r-')
               ax[i,j].plot(bincenters,sk[idx,:],'b+')
               ax[i,j].set_title( 't='+str(t[idx])[:5] )
          # end for
     # end for
# end if

pyplot.tight_layout()
pyplot.show(block=False)

