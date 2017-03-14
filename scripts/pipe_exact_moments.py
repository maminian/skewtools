import h5py
from numpy import *
import sys
import scripts.skewtools as st
from matplotlib import pyplot
import scipy.special as sf	# Special functions

#

nmax=1000

bzs = -sf.jn_zeros(1,nmax) # Tabulate the first nmax zeroes of Bessel -J_1 for the sums.

def m2(t,Pe=10.**4):
     out = zeros(shape(t))
     for i in range(nmax):
          out += bzs[i]**-8 * exp(-bzs[i]**2*t)
     # end for
     
     out *= 128*Pe**2
     out += -1./360*Pe**2 + 2*(1. + 1./48*Pe**2)*t

     return out
# end def

def m3(t,Pe=10.**4):
     out = zeros(shape(t))
     for i in range(nmax):
          out += (t*bzs[i]**-8 + 18*bzs[i]**-10 - 240*bzs[i]**-12)*exp(-bzs[i]**2*t)
     # end for
     
     out *= 128*Pe**3
     out += 1./480*Pe**3*(t-17./112)

     return out
# end def

# --------------------------------------------- 
#
# Keep in mind these are for the flow u=2(1-r**2).
#

t,var,sk = st.importDatasets(sys.argv[1],'Time','Avgd_Variance','Avgd_Skewness')

fig,ax = pyplot.subplots(2,1)

tv = t[1:]
skv = sk[1:]

exact = m3(tv)/m2(tv)**1.5

ax[0].plot(tv,exact,label='Exact')
ax[0].scatter(tv,skv,facecolor=[0,0,0,0],edgecolor=[0,0,0,1],s=40,label='Monte Carlo')

ax[1].plot(tv,abs(exact-skv),color='red',label='Absolute error')
ax[1].plot(tv,abs(exact-skv)/abs(exact),color=[0.4,0,0.4],label='Relative error')

ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[1].set_yscale('log')

ax[0].set_xlim([tv[0],tv[-1]])
ax[1].set_xlim([tv[0],tv[-1]])

ax[0].set_ylim([-0.2,0.5])

ax[0].grid(True)
ax[1].grid(True)

ax[0].legend(loc='best')
ax[1].legend(loc='best')

pyplot.tight_layout()
pyplot.show(block=False)
