from numpy import *
from matplotlib import pyplot
import scripts.skewtools as st
import sys


Y,Z,t,l,q = st.importDatasets(sys.argv[1],'Y','Z','Time','aratio','q')


i=0

figscale = 5.
fig,ax = pyplot.subplots(1,1,figsize=(figscale/l,figscale))

ax.cla()


ax.scatter(Y[:,i],Z[:,i],edgecolor=None)

winscale = 1.5
ax.set_xlim([-winscale/l,winscale/l])
ax.set_ylim([-winscale,winscale])

pyplot.savefig('cloudframe_'+str(i).zfill(4)+'.png',dpi=80)
