import matplotlib.pyplot as pyplot
from matplotlib import cm
from numpy import linspace,sqrt

# h5py must be either installed on the system using a generic 
# 'pip ...' command, or locally, making sure to specify the correct 
# settings with pip so that the 'built' package is installed locally.
# Then make sure that the installed directory is in the PYTHONPATH.
import sys
import h5py
import glob

import skewtools

# Read the input arguments as a list of files to pull data from.
# Then pull the information from each of those files.

if len(sys.argv) == 2:
     files = glob.glob(sys.argv[1])
else:
     files = sys.argv[1:]
# end if

files.sort()

kurts = skewtools.gatherDataFromFiles(files,'Avgd_Kurtosis')

# Assume the same times for all files to save space.
time = skewtools.gatherDataFromFiles([files[0]],'Time')[0]

# Set up a colormap.
cnum=0
cmaps = cm.jet(linspace(0,1,len(kurts)))
idx = 0

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.hold(True)

for kurt in kurts:
     time = skewtools.gatherDataFromFiles([files[cnum]],'Time')[0]
     ax.plot(time, kurt, markersize=5, color=(0.,0.,1.,min(1.,2./len(kurts))))
     ax.plot(time, kurt, markersize=5, color=cmaps[cnum])
     cnum += 1
# end for


# Displaying mean and +/- 3 stdev when we have
# multiple input files.

horizx = [time[1],time[-1]]

if ( False & (len(kurts)>2) ):

     samp_mean = skewtools.meanLine(kurts)
     stdl,stdu = skewtools.stdLines(kurts,3)

     ax.plot(time,samp_mean,linewidth=3,linestyle='-',color='black')
     ax.plot(time,stdl,linewidth=1,linestyle='-',color='black')
     ax.plot(time,stdu,linewidth=1,linestyle='-',color='black')

     samp_min,idx = skewtools.constrainedMinimum(time,samp_mean,time[0],time[500])


#     ax.plot(horizx,[samp_min,samp_min],linewidth=1,linestyle='-',color='red')
     
# end if

# Plot decorations

ax.hold(False)

pyplot.xlabel('Time')
pyplot.ylabel('Kurtosis')
pyplot.grid(True)
pyplot.xscale('log')
#legend(loc='best')

# Wrap-up
#title('Skewness for point-source initial data '+r'$\delta(x)\delta(y-a)$, $a=0,0.1,...,0.5, Pe=10^9.$')
pyplot.show(block=False)
pyplot.savefig('mc_kurtosis.png',dpi=150)
