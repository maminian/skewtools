#!/usr/bin/python

import numpy as np
#import scipy.stats
import sys

import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec


# ------------------------------------------------

def update(i,ax,hist_centers,hist_heights,t,avsk):

     pyplot.cla()
     
#     samp_mean = mean(x[i,:])
     if (i<2):
          xmax = max( [abs(hist_centers[1,:].max()), abs(hist_centers[1,:].min())] )
     else:
          xmax = 0.
          for k in range(i):
               xmax = max( [xmax,abs(hist_centers[i,:].max()), abs(hist_centers[i,:].min())] )
          #
     # end if
     ymax = hist_heights[i,:].max()

     ax.step(hist_centers[i,:],hist_heights[i,:],where='mid',color='blue',label=r'$C(x,\tau)$')
     # TEMPORARY PLOTTING C(-x)
     ax.step(-hist_centers[i,:],hist_heights[i,:],where='mid',color='red',linewidth=0.5,label=r'$C(-x,\tau)$')

     # Find the mean, median, and plot.
     xmed = find_median(hist_heights[i,:],hist_centers[i,:])
     ax.plot([xmed,xmed],[0,ymax],color='black',label='median')

     xmean = sum(hist_heights[i,:]*hist_centers[i,:])/sum(hist_heights[i,:])
     ax.plot([xmean,xmean],[0,ymax],color='green',label='mean')
     
     ax.set_title(r'$\tau = %.3g, \, \, Sk = %.3g$'%(t[i],avsk[i]) )

     ax.set_xlim([-1.1*xmax,1.1*xmax])
     ax.set_ylim([0,1.1*ymax])

     ax.legend(loc='upper left',fancybox=True,framealpha=0.5)
#

def animate_hist(hist_centers,hist_heights,t,avsk):

     # Generate a subplot with the particle animation on top
     # and an animated tracker of the skewness over time
     # on the bottom.
#     gs = gridspec.GridSpec(2,1, height_ratios=[3,1])

     fig,ax = pyplot.subplots(1,1)

     ani = anim.FuncAnimation( fig, update, frames=np.size(t), 
                                   fargs=(ax,hist_centers,hist_heights,t,avsk),repeat=False)
     
     return ani

def find_median(hght,centers):
     # Finds the approximate x location centers[i] of a histogram
     # with i approximately satisfying sum(hght[:i])==sum(hght[:i]).
     s=0
     shot = sum(hght)/2.
     loc = -1.
     for i in range(len(centers)):
          s += hght[i] 
          if (s > shot):
               loc=(centers[i-1]+centers[i])/2.
               break
          elif (s==shot):
               loc=centers[i]
               break
          # end if
     # end for

     return loc
# end def

# --------------------------------------------


walks = {};

print "Reading data..."

#
# For simulations with position histories, 
# bin the 2d data.
#

walks = h5py.File(sys.argv[1])


t = np.transpose(walks['Time'].value)
hist_centers = np.transpose(walks['Hist_centers'].value)
hist_heights = np.transpose(walks['Hist_heights'].value)
avsk = np.transpose(walks['Avgd_Skewness'].value)

Pe = walks['Peclet'].value
#aratio = walks['aratio'].value

walks.close()
     
print "Constructing animation..."
fps = 24
ani = animate_hist(hist_centers,hist_heights,t,avsk)


#print "Saving animation to disk (may take quite a while, may hang, watch the process)..."
okska=1.
ani.save('conc_pdf_%.2f.mp4'%okska,dpi=240,fps=12)

pyplot.show(block=False)

print "Done."

