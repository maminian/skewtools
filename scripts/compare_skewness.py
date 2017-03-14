import matplotlib.pyplot as pyplot
from matplotlib import cm
from numpy import linspace,sqrt,log10,ones,mean,shape

# add a path for h5py for Kure.
import sys

# = %pwd
#h5pydir = workdir+'/lib/python-2.7/site-packages/h5py/'
#sys.path.append(str(h5pydir))

#import highlevel as h5py

import h5py

import skewtools

# Read the input files. The first file will be used as a reference,
# and the rest will be compared against it.
reference = [sys.argv[1]]
files = sys.argv[2:]

# Generate a list
tref = skewtools.gatherDataFromFiles(reference,'Time')[0]

# May need to change the name of the dataset here.
skref = skewtools.gatherDataFromFiles(reference,'Avgd_Skewness_Hybrid')[0]
kuref = skewtools.gatherDataFromFiles(reference,'Avgd_Kurtosis_Hybrid')[0]-3.

times = skewtools.gatherDataFromFiles(files,'Time')
skews = skewtools.gatherDataFromFiles(files,'Avgd_Skewness')
kurts = skewtools.gatherDataFromFiles(files,'Avgd_Kurtosis')


# Get a time interval which all the simulations share.
tmin = tref.min()
tmax = tref.max()

for ta in times:
     if (min(ta[1:]) > tmin):
          tmin = min(ta[1:])
     # end if
     if (max(ta[1:]) < tmax):
          tmax = max(ta[1:])
     # end if
# end for

# Generate a common set of times to compare against, 
# then loop through and compare them all to the reference.
tcomp = 10.**linspace(log10(tmin),log10(tmax),1000)


# Set up the plotting.
fig = pyplot.figure()
ax1 = fig.add_subplot(2,1,1)


# Set up a colormap.
cnum=0
cmaps = cm.rainbow(linspace(0,1,len(skews)+1))

ax1.plot(tref,skref,label='Reference',color=cmaps[0])
pyplot.ylabel('Skewness')
pyplot.grid(True)
pyplot.xscale('log')

ax2 = fig.add_subplot(2,1,2)
pyplot.grid(True)
pyplot.xscale('log')
pyplot.yscale('log')

pyplot.xlabel('Time')
pyplot.ylabel('Absolute Error')

# Finally, loop through, and plot the results.
ax1.hold(True)
ax2.hold(True)

cnum+=1
for i,_ in enumerate(times):
     errs = skewtools.compare_functional_diff(tcomp,tref,skref,times[i],skews[i])
     ax1.plot(times[i],skews[i],label=files[i],color=cmaps[cnum])
     
     ax2.plot(tcomp,abs(errs),linewidth=1,label=files[i],color=cmaps[cnum])
     
     avgline = mean(abs(errs))*ones(shape(tcomp))
#     ax2.plot(tcomp,avgline,linestyle='--',linewidth=4,color=cmaps[cnum])
     
     cnum+=1
# end for

pyplot.sca(ax1)
pyplot.legend(loc='best')

ax1.hold(False)
ax2.hold(False)

ax1.set_xlim([tmin,tmax])
ax1.set_ylim([1.1*min(skews[0]),1.1*max(skews[0])])

ax2.set_xlim([tmin,tmax])

pyplot.savefig('skewness_comp.png',dpi=150)

# -----------------------
# Repeat everything for the kurtosis.
#

# Set up the plotting.
fig2 = pyplot.figure()
ax3 = fig2.add_subplot(2,1,1)


# Set up a colormap.
cnum=0
cmaps = cm.rainbow(linspace(0,1,len(kurts)+1))

ax3.plot(tref,kuref,label='Reference',color=cmaps[0])
pyplot.ylabel('Kurtosis')
pyplot.grid(True)
pyplot.xscale('log')

ax4 = fig2.add_subplot(2,1,2)
pyplot.grid(True)
pyplot.xscale('log')
pyplot.yscale('log')

pyplot.xlabel('Time')
pyplot.ylabel('Absolute Error')

# Finally, loop through, and plot the results.
ax3.hold(True)
ax4.hold(True)


cnum+=1
for i,_ in enumerate(times):
     errs = skewtools.compare_functional_diff(tcomp,tref,kuref,times[i],kurts[i])
     ax3.plot(times[i],kurts[i],label=files[i],color=cmaps[cnum])
     
     ax4.plot(tcomp,abs(errs),linewidth=1,label=files[i],color=cmaps[cnum])
     
#     avgline = mean(abs(errs))*ones(shape(tcomp))
#     ax2.plot(tcomp,avgline,linestyle='--',linewidth=4,color=cmaps[cnum])
     
     cnum+=1
# end for

pyplot.sca(ax3)
pyplot.legend(loc='best')

ax3.hold(False)
ax4.hold(False)

ax3.set_xlim([tmin,tmax])
ax3.set_ylim([1.1*min(kurts[0]),1.1*max(kurts[0])])

ax4.set_xlim([tmin,tmax])

# Wrap-up
#title('Skewness for point-source initial data '+r'$\delta(x)\delta(y-a)$, $a=0,0.1,...,0.5, Pe=10^9.$')
pyplot.show(block=False)

pyplot.savefig('kurtosis_comp.png',dpi=150)
