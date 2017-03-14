from numpy import *
from matplotlib import pyplot
import scripts.skewtools as st
import sys


#X,Y,t,avvar,avmean = st.importDatasets(sys.argv[1],'X','Y','Time','Avgd_Variance','Avgd_Mean')
histcx,histcy,hists,t,avvar,avmean = st.importDatasets(sys.argv[1],'hist2dcx','hist2dcy','hist2d','Time','Avgd_Variance','Avgd_Mean')

figscale = 5.
fig,ax = pyplot.subplots(1,1,figsize=(5*figscale,figscale))

#xmin = X.min()
#xmax = X.max()

mycm = pyplot.cm.viridis

for i in range(len(t)):

     ax.cla()

#     ax.scatter(X[:,i],Y[:,i],edgecolor=None,s=1.,alpha=0.1)
     Xg,Yg = meshgrid(histcx[:,i],histcy[:,i])
#     ax.pcolor(Xg,Yg,hists[:,:,i],cmap=mycm)
     if (i < 5):
          mmax = 0.
          for j in range(i-4,i+1):
               mmax += 1./5.*hists[:,:,j].max()
          #
     else:
          mmax = hists[:,:,(i-5):i].max()
     #

     ax.contourf(Xg,Yg,hists[:,:,i],cmap=mycm,levels=linspace(0,1,11)*mmax)
#     ax.plot([avmean[i],avmean[i]],[-1.,1.],linewidth=1,color='red')
#     xmin = min(xmin,min(X[:,i]))
#     xmax = max(xmin,max(X[:,i]))
     ax.set_facecolor(mycm.colors[0])
     ax.set_xlim([histcx.min(),histcx.max()])
     ax.set_ylim([-1.,1.])

     pyplot.tight_layout()
     pyplot.savefig('cloudframe_'+str(i).zfill(4)+'.png',dpi=80,bbox_inches='tight')
     print '%i of %i'%(i,len(t)-1)
#
