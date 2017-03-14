import scripts.skewtools as st
from numpy import linspace,unique,shape,size,abs,argmin,mod
import glob
import matplotlib.pyplot as pyplot
import sys


fl = glob.glob(sys.argv[1])

fl.sort()

Pes = st.gatherDataFromFiles(fl,'Peclet')
#aratios = st.gatherDataFromFiles(fl,'aratio')
hcs = st.gatherDataFromFiles(fl,'Hist_centers')
hhs = st.gatherDataFromFiles(fl,'Hist_heights')

# Assume the arrays share a common time stepping.
t = st.gatherDataFromFiles([fl[0]],'Time')[0]


m,n = 3,3

tidxs = linspace(size(t)/2.,size(t)-1,m*n).astype('int')
for i in range(len(tidxs)):
	print 'tau=%.2g'%t[tidxs[i]]
# end for


fig,ax = pyplot.subplots(m,n,figsize=(3*m,3*n))

#alist = unique(aratios)
#alist.sort()
#mycm=pyplot.cm.gnuplot2(linspace(0,0.7,len(alist)))

# for outlining the time in the figures..
props = dict(boxstyle='round',facecolor='white',alpha=0.8)


for idx in range(m*n):

	j = idx%n
	i = (idx-j)/n

	for k in range(len(fl)):

#		cidx = argmin(abs(alist-aratios[k]))
#		ax[i,j].step(hcs[k][:,tidxs[idx]],hhs[k][:,tidxs[idx]],where='mid',lw=0.5,color=mycm[cidx],alpha=0.4)
		ax[i,j].step(hcs[k][:,tidxs[idx]],hhs[k][:,tidxs[idx]],where='mid',lw=0.5,alpha=0.4)
	# end for
	ax[i,j].set_yticklabels([])
	ax[i,j].set_xticklabels([])
	ax[i,j].grid()
	
	ax[i,j].text(0.05,0.8,r'$\tau=%.2g$'%t[tidxs[idx]],transform=ax[i,j].transAxes,
					color='black',fontsize=16,ha='left',va='center',bbox=props)
	# end for
	
#	if (idx==0):
#		for p in range(len(alist)):
#			ax[i,j].plot(0,0,lw=1,color=mycm[p],label=r'$\lambda=%.4f$'%alist[p])
		# end for
	# end if
# end for
ax[0,0].legend(loc='lower right')

fig.suptitle(r'$Pe = %i$'%Pes[0],ha='center',va='center',fontsize=36)

pyplot.tight_layout()
pyplot.subplots_adjust(top=0.95)
pyplot.savefig('hist_panels.png',dpi=120,bbox_inches='tight')
pyplot.show(block=False)
