from numpy import *
from matplotlib import pyplot
import scripts.skewtools as st
import sys


X,Y,t,Pe = st.importDatasets(sys.argv[1],'X','Y','Time','Peclet')


figscale = 5.
fig,ax = pyplot.subplots(1,1,figsize=(4*figscale,figscale))

uwall = 2./3.
xmax = X.max() + uwall*Pe*t[-1]

for i in range(len(t)):
#for i in [52]:

     ax.cla()
#     ax.hold(True)
     ax.plot([0,xmax],[1,1],linewidth=0.5,color='k')
     ax.plot([0,xmax],[-1,-1],linewidth=0.5,color='k')
     subset1 = ((Y[:,i]<1.)*(Y[:,i] > -1.))
     subset2 = ~subset1


     ax.scatter(X[subset1,i],Y[subset1,i],facecolor=[0,0,0.9],edgecolor=[0,0,0,0],s=1,alpha=0.2)
     ax.scatter(X[subset2,i],Y[subset2,i],facecolor=[0.9,0,0],edgecolor=[0,0,0,0],s=1,alpha=0.2)
#     ax.hist2d(X[subset,i] + uwall*Pe*t[i],Y[subset,i],cmap=pyplot.cm.inferno,)
#     ax.hold(False)
     
     ax.set_xlim([0.,xmax])
     ax.set_ylim([-1.05,1.05])
     print '%i active particles, %i of %i frames'%(sum(subset1),i,len(t)-1)
     
     pyplot.savefig('cloudframe_'+str(i).zfill(4)+'.png',dpi=80,bbox_inches='tight')
# end for


#pyplot.tight_layout()

