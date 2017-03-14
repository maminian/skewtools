#!/usr/bin/python

from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Script to read the output file from cchoose_test, which
# contains the index set maximizing the l_1 norm of
# over partial sums of a fixed number of terms.
#
# Only supports reading a single output file.


idxset = {};

fo = open(sys.argv[1],'r');
raw = fo.readlines();
fo.close();


idxset["aratio"] = float(raw[0]);
idxset["Nterms"] = int(raw[1]);

idxset["i"] = zeros( (idxset["Nterms"],1) );
idxset["j"] = zeros( (idxset["Nterms"],1) );
idxset["u_ij"] = zeros( (idxset["Nterms"],1) );

k = 0;
for line in raw[2:]:
     triple = line.split();
     idxset["i"][k] = int(triple[0]);
     idxset["j"][k] = int(triple[1]);
     idxset["u_ij"][k] = float(triple[2]);
     k += 1;

# Scatterplot the indices and produce a contour map of
# the Fourier coefficient function.
imax = max(idxset['i']);
jmax = max(idxset['j']);

X,Y = meshgrid(linspace(1.,imax+1,101),linspace(1.,jmax+1,101));
Z = (-8./pi**4)/((X-0.5)*(Y-0.5)*((X-0.5)**2 + (idxset['aratio']*(Y-0.5))**2));

fig = plt.figure()
ax = plt.gca()

con=ax.contourf(X,Y,abs(Z),10.**arange(-10.,1.,1.),cmap=plt.cm.Accent,norm=LogNorm())

# Setting the axes to be equal may become unreasonable for small aspect ratio;
# be careful with this setting.
#ax.set_aspect('equal', adjustable='box')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(con,cax=cax)

ax.scatter(idxset['i'],idxset['j'],marker='.',color='black')

sca(ax)
xlabel('Summation index i')
ylabel('Summation index j')
title(str(idxset['Nterms'])+' summation terms on the lattice, with contour map\n of Fourier coefficient with aspect ratio '+str(idxset["aratio"]))


show(block=False)

savefig('terms'+str(idxset['Nterms'])+'.png',dpi=150,bbox_inches='tight')
