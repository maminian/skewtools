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

skews = st.gatherDataFromFiles(files,'Avgd_Skewness')
times = st.gatherDataFromFiles(files,'Time')
Pes = st.gatherDataFromFiles(files,'Peclet')

cmaps = cm.jet(linspace(0,1,len(files)))

fig,ax = pyplot.subplots(2,1,figsize=(12,9))

for i in range(len(files)):
     t = times[i][1:]
     sk = skews[i][1:]

     exact = st.m3_channel(t,Pes[i])/((st.m2_channel(t,Pes[i]))**1.5)
     ax[0].plot(t,sk,color=cmaps[i],marker='.')
     ax[0].plot(t,exact,color='black',linewidth=1)
     if (i==0):
          ax[1].plot(t,abs(exact-sk),color=cmaps[i],marker='.',label='Abs. error')
          ax[1].plot(t,abs(exact-sk)/abs(exact),color=cmaps[i],marker='*',label='Rel. error')
     else:
          ax[1].plot(t,abs(exact-sk),color=cmaps[i],marker='.')
          ax[1].plot(t,abs(exact-sk)/abs(exact),color=cmaps[i],marker='*')
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
ax[1].set_ylim([10.**-6,10.**1])

ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[1].set_yscale('log')

ax[0].grid(True)
ax[1].grid(True)
ax[1].legend(loc='upper left')

pyplot.tight_layout()

pyplot.show(block=False)
pyplot.savefig('mc_skewness_compex.png',bbox_inches='tight',dpi=120)
