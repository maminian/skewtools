import matplotlib.pyplot as pyplot
from matplotlib import cm
from numpy import linspace,sqrt

import numpy as np

# h5py must be either installed on the system using a generic 
# 'pip ...' command, or locally, making sure to specify the correct 
# settings with pip so that the 'built' package is installed locally.
# Then make sure that the installed directory is in the PYTHONPATH.
import sys
import h5py
import glob

import skewtools

def slant_grid(ax,**kwargs):

     ax.hold(True)
     
     # Construct gridlines with slopes 1 and 2.
     xstart = ax.xaxis.get_majorticklocs()
     ystart = ax.yaxis.get_majorticklocs()

     mxlim = ax.get_xlim()
     mylim = ax.get_ylim()

     print mxlim

     xv = np.array(mxlim)
     yv = np.array(mylim)
     
     xv2 = xv.copy()
     xv2[0] *= yv[0]/yv[1]

     dx = xstart[2]/xstart[0]
     val = xv2[0]
     
     while (val < xv2[1]):
          x = xv2/xv2[0]*val
          y1 = x/val*xv[0]
          y2 = (x/val)**2*xv[0]
          ax.plot(x,y1,**kwargs)
          ax.plot(x,y2,**kwargs)
          val *= dx

     # end while

#     print ystart
#     for i in range(1,len(ystart),2):
#          print ystart[i]
#          ax.plot(xv,xv/xv[0]*ystart[i],color='black',linewidth=0.4,linestyle=':')
#          ax.plot(xv,xv**2/xv[0]**2*ystart[i],color='black',linewidth=0.4,linestyle=':')
     # end for

     ax.hold(False)
     
     ax.set_xlim(mxlim)
     ax.set_ylim(mylim)
     
     return ax
# end def

# Read the input arguments as a list of files to pull data from.
# Then pull the information from each of those files.

if len(sys.argv) == 2:
     files = glob.glob(sys.argv[1])
else:
     files = sys.argv[1:]
# end if

files.sort()

Pe = skewtools.gatherDataFromFiles(files,'Peclet')[0]

variances = skewtools.gatherDataFromFiles(files,'Avgd_Variance')

# Assume the same times for all files to save space.
times = skewtools.gatherDataFromFiles(files,'Time')

# Set up a colormap.
cnum=0
cmaps = cm.jet(linspace(0,1,len(variances)))
idx = 0

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

ax.set_xlabel('Time')
ax.set_ylabel( r'$\sigma^2$')
ax.set_xscale('log')
ax.set_yscale('log')


ax.hold(True)

ch_eff = 1. + Pe**2*8./945.
pipe_eff = 1. + Pe**2*1./768.

for i in range(len(variances)):

#     Pe = skewtools.gatherDataFromFiles(files,'Peclet')[cnum]
     
#     time = skewtools.gatherDataFromFiles([files[cnum]],'Time')[0]
#     ax.plot(time, variance/(2.*time*Pe**2), markersize=5, color=(0.,0.,1.,min(1.,2./len(variances))))
     ax.plot(times[i][1:], (variances[i][1:]), markersize=5, color=cmaps[cnum],label=files[cnum])

     print ax.get_xlim()

     cnum += 1
# end for


# Displaying mean and +/- 3 stdev when we have
# multiple input files.

limlist = np.zeros( (np.shape(times)[0],2) )
for i in range(np.shape(times)[0]):
     limlist[i,0] = times[i][1]
     limlist[i,1] = times[i][-1]
# end for

horizx = [min(limlist[:,0]),max(limlist[:,1])]


# Plot decorations


#ax.plot([time[1],time[-1]],[ch_eff,ch_eff],color='red')
#ax.plot([time[1],time[-1]],[pipe_eff,pipe_eff],color='green')

#ax.hold(False)


#ax.grid(True)

ax.set_xlim(horizx)
#ax.set_ylim([])

slant_grid(ax,color=[0.1,0.1,0.1],linewidth=0.2,linestyle='--')

ax.set_xscale('log')
ax.set_yscale('log')

#pyplot.legend(loc='best')

# Wrap-up
#title('Skewness for point-source initial data '+r'$\delta(x)\delta(y-a)$, $a=0,0.1,...,0.5, Pe=10^9.$')
pyplot.show(block=False)
pyplot.savefig('mc_variance.png',dpi=150)
