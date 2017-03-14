#!/usr/bin/python

import numpy as np
import sys

import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
import h5py

from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec


# ------------------------------------------------
# Given the cross-sectionally averaged distribution, 
# construct the symmetric and antisymmetric parts, 
# and see what the mean and median look like.
#
# Functionally related to animate_pdf_mc.py.
# --------------------------------------------


walks = h5py.File(sys.argv[1])


t = np.transpose(walks['Time'].value)
hist_centers = np.transpose(walks['Hist_centers'].value)
hist_heights = np.transpose(walks['Hist_heights'].value)
avsk = np.transpose(walks['Avgd_Skewness'].value)
avmean = np.transpose(walks['Avgd_Mean'].value)

nwalkers = int(walks['nTrials'].value)
Pe = walks['Peclet'].value
aratio = walks['aratio'].value

walks.close()

# Construct symmetric and antisymmetric parts. This has to be done on 
# a symmetrized grid (via padding the histograms with zeros)

hc2 = {}
hh2 = {}

if False:
     for i in range(len(t)):
          # Determine whether the histogram centers lean to left or right.
          hc = hist_centers[i,:]
          ht = hist_heights[i,:]

          dx = hc[1]-hc[0]
          if (hc[-1]+hc[0])<0:
               hct = np.concatenate( (hc,np.arange(hc[-1],-hc[0]-dx,dx)), 0)
               hct = np.concatenate( (hct,[-hc[0]]), 0)
               
#               hht = np.concatenate( (ht,np.zeros( len(hct

               hc2.append( [hct,hht] ) 
          # end if
          
     # end for
