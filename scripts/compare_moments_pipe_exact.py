import scripts.skewtools as st
from numpy import *
from matplotlib import pyplot,cm
import sys
import glob

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

variances = st.gatherDataFromFiles(files,'Avgd_Variance')
skews = st.gatherDataFromFiles(files,'Avgd_Skewness')
times = st.gatherDataFromFiles(files,'Time')
Pes = st.gatherDataFromFiles(files,'Peclet')

cmaps = cm.jet(linspace(0,1,len(files)))

fig,ax = pyplot.subplots(1,2,figsize=(12,9))

for i in range(len(files)):
     t = times[i]

     mvar = variances[i]
     msk = skews[i]
     mm3 = msk*mvar**1.5
     
     exactm2 = st.m2_pipe(t,Pe=Pes[i],nmax=1000,m20=mvar[0])
     exactm3 = st.m3_pipe(t,Pe=Pes[i],nmax=1000,m30=mm3[0])

     ax[0].plot(t[1:],mvar[1:],color=cmaps[i],marker='.')
     if (i==0):
          ax[0].scatter(t[1:],exactm2[1:],marker='o',s=100,edgecolor='black',facecolor=[0,0,0,0],label='Exact')
     else:
          ax[0].scatter(t[1:],exactm2[1:],marker='o',s=100,edgecolor='black',facecolor=[0,0,0,0])
     # end if

     ax[1].plot(t[1:],mm3[1:],color=cmaps[i],marker='.')
     if (i==0):
          ax[1].scatter(t[1:],exactm3[1:],marker='o',s=100,edgecolor='black',facecolor=[0,0,0,0],label='Exact')
     else:
          ax[1].scatter(t[1:],exactm3[1:],marker='o',s=100,edgecolor='black',facecolor=[0,0,0,0])
     # end if

# end for

tmin=100.
tmax=0.
for time in times:
     tmin = min(tmin,time[1])
     tmax = max(tmax,time[-1])
# end for

ax[0].set_xlim([tmin,tmax])
ax[1].set_xlim([tmin,tmax])
#ax[1].set_ylim([10.**-6,10.**1])

ax[0].set_xscale('log')
ax[1].set_xscale('log')

ax[0].set_yscale('log')
ax[1].set_yscale('log')

ax[0].grid(True)
ax[1].grid(True)

ax[0].legend(loc='upper left')
ax[1].legend(loc='upper left')

fs=36
ax[0].text(0.5,0.8,r'$M_2$',transform=ax[0].transAxes,fontsize=fs)
ax[1].text(0.5,0.8,r'$M_3$',transform=ax[1].transAxes,fontsize=fs)

pyplot.tight_layout()

pyplot.show(block=False)
pyplot.savefig('pipe_moments_compare.png',bbox_inches='tight',dpi=120)
