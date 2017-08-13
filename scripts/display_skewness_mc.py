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

print sys.argv

if len(sys.argv) == 2:
     if (type(sys.argv[1]) == str):
          files = glob.glob(sys.argv[1])
     elif (type(sys.argv[1]) == list):
          files = sys.argv[1]
     else:
          print "Unrecognized input."
     # end if
else:
     files = sys.argv[1:]
# end if

files.sort()
print files

skews = skewtools.gatherDataFromFiles(files,'Avgd_Skewness')

# Assume the same times for all files to save space.
time = skewtools.gatherDataFromFiles([files[0]],'Time')[0]

# Set up a colormap.
cnum=0
cmaps = cm.jet(linspace(0,1,len(skews)))
idx = 0

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.hold(True)



for skew in skews:
     time = skewtools.gatherDataFromFiles([files[cnum]],'Time')[0]
#     ax.plot(time[1:], skew[1:], markersize=5, color=(0.,0.,1.,min(1.,2./len(skews))))
#     ax.plot(time[1:], skew[1:], markersize=5, color=cmaps[cnum], label=files[cnum])
     ax.plot(time[1:], skew[1:], linewidth=0.5, color=cmaps[cnum])
     ax.scatter(time[1:],skew[1:],s=4,color=cmaps[cnum],label=files[cnum])
     cnum += 1
# end for


# Displaying mean and +/- 3 stdev when we have
# multiple input files.

horizx = [time[1],time[-1]]

if ( False & (len(skews)>2) ):

     samp_mean = skewtools.meanLine(skews)
     stdl,stdu = skewtools.stdLines(skews,3)

     ax.plot(time,samp_mean,linewidth=3,linestyle='-',color='black')
     ax.plot(time,stdl,linewidth=1,linestyle='-',color='black')
     ax.plot(time,stdu,linewidth=1,linestyle='-',color='black')

     samp_min,idx = skewtools.constrainedMinimum(time,samp_mean,time[0],time[500])


#     ax.plot(horizx,[samp_min,samp_min],linewidth=1,linestyle='-',color='red')
     
# end if

# Plot decorations

ax.hold(False)

ax.set_xlabel('Time')
ax.set_ylabel('Skewness')
ax.grid(True)

ax.set_xscale('log')

#if len(files)<5:
#     pyplot.legend(loc='best')
# end if

# Wrap-up
#title('Skewness for point-source initial data '+r'$\delta(x)\delta(y-a)$, $a=0,0.1,...,0.5, Pe=10^9.$')
pyplot.show(block=False)
pyplot.savefig('mc_skewness.png',dpi=150)
