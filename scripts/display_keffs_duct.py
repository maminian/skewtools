import sys
from numpy import *
from matplotlib import pyplot
import glob

import scripts.skewtools as st

if len(sys.argv) == 2:
     files = glob.glob(sys.argv[1])
else:
     files = sys.argv[1:]
# end if

aratios = st.gatherDataFromFiles(files,'aratios')
keffs = st.gatherDataFromFiles(files,'keffs')
keffs_asymp = st.gatherDataFromFiles(files,'keffs_asymp')


fig,ax = pyplot.subplots(1,1)

chkeff = 8./945

for i in range(len(aratios)):
     ax.plot(aratios[i],keffs[i]/chkeff,label='Exact')
     ax.plot(aratios[i],keffs_asymp[i]/chkeff,label='Asymptotic')
# end for

ax.set_xticks(arange(0.,1.,0.1))
ax.set_yticks(arange(0.,9.,1.))

ax.grid(True)
ax.legend()
pyplot.show(block=False)
