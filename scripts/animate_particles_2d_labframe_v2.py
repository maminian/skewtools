from numpy import *
from matplotlib import pyplot
import scripts.skewtools as st
import sys


X,Y,t,Pe = st.importDatasets(sys.argv[1],'X','Y','Time','Peclet')


figscale = 5.
fig,ax = pyplot.subplots(2,1,figsize=(4*figscale,1.5*figscale))

orig0 = ax[0].get_position()
orig1 = ax[1].get_position()
ax[0].set_position([orig0.x0,orig0.y0-0.2,orig0.width,orig0.height+0.2])
ax[1].set_position([orig1.x0,orig1.y0,orig1.width,orig1.height-0.2])

np = shape(X)[0]
uwall = 2./3.
xmax = X.max()

nabsorbed = zeros(shape(t))

for i in range(len(t)):
#for i in [12]:

     ax[0].cla()
#     ax.hold(True)
     ax[0].plot([0,xmax],[1,1],linewidth=0.5,color='k')
     ax[0].plot([0,xmax],[-1,-1],linewidth=0.5,color='k')
     subset1 = ((Y[:,i]<1.)*(Y[:,i] > -1.))

     subset2 = ~subset1
     subset2a = (Y[:,i]>1.)
     subset2b = (Y[:,i]<-1.)

     nbinsx = 401
     nbinsy = 101
     vmaxval = np/nbinsx
     

#     ax[0].scatter(X[subset2,i],Y[subset2,i],facecolor=[1,0,0],edgecolor=[0,0,0,0],s=1,alpha=0.1)
     ax[0].hist2d(X[subset2a,i],Y[subset2a,i],cmap=pyplot.cm.Reds,bins=[linspace(0,xmax,nbinsx),linspace(1,1.05,2)],vmin=0,vmax=vmaxval)
     ax[0].hist2d(X[subset2b,i],Y[subset2b,i],cmap=pyplot.cm.Reds,bins=[linspace(0,xmax,nbinsx),linspace(-1.05,-1,2)],vmin=0,vmax=np/nbinsx)

     ax[0].scatter(X[subset1,i],Y[subset1,i],facecolor=[1,0,0,0.2],edgecolor='none',s=1)
#     ax[0].hist2d(X[subset1,i],Y[subset1,i],cmap=pyplot.cm.inferno,bins=[linspace(0,xmax,nbinsx),linspace(-1,1,nbinsy)])
#     ax.hold(False)
     
     ax[0].set_xlim([0.,xmax])
     ax[0].set_ylim([-1.05,1.05])

     ax[1].cla()

     ax[1].hist(X[subset2,i],bins=linspace(0.,xmax,nbinsx),facecolor=[0,0,0,0],edgecolor='k')
     ax[1].set_xlim([0.,xmax])
     ax[1].set_ylim([0.,4*np/nbinsx])


     
     pyplot.savefig('particleframe_'+str(i).zfill(4)+'.png',dpi=80,bbox_inches='tight')
     print '%i active particles, %i of %i frames'%(sum(subset1),i,len(t)-1)
# end for


#pyplot.tight_layout()

