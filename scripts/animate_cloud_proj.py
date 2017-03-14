# Visualize side view cloud animation 
# via internally calculated hist2d instead of saving 
# particle histories and doing the calculation with 
# python's hist2d.
#

from numpy import *
from matplotlib import pyplot
import scripts.skewtools as st
import sys

mf = sys.argv[1]

t,hcx,hcy,hh = st.importDatasets(mf,'Time','hist2dcx','hist2dcy','hist2d')

nhb = shape(hcx)[0]
nby = shape(hcy)[0]

fig,ax = pyplot.subplots(1,1)

#for i in [61]:
for i in range(len(t)):
     ax.cla()

     Xc,Yc = meshgrid(hcx[:,i],hcy[:,i])
     ax.pcolor(Xc,Yc,hh[:,:,i],cmap=pyplot.cm.inferno)

     ax.set_xlim(Xc.min(),Xc.max())
     ax.set_ylim(Yc.min(),Yc.max())

     pyplot.savefig('cloud_proj_'+str(i-7).zfill(4)+'.png')
     print '%i of %i'%(i,len(t)-1)
#

pyplot.show(block=False)
