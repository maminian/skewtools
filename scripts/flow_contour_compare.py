from numpy import *
import scripts.skewtools as st
from matplotlib import pyplot
import sys

aratio,uv,y,z = st.importDatasets(sys.argv[1],'aratio','Flow_values','y_mesh','z_mesh')

print aratio

z = z/aratio

fig,ax = pyplot.subplots(1,1, figsize=(6/aratio,6))

yv,zv = meshgrid(y,z)

mycm = pyplot.cm.viridis

uell = 1.-yv**2-aratio**2*zv**2

ax.contourf(zv,yv,uv,cmap=mycm,levels=linspace(uv.min(),uv.max(),9))
ax.contour(zv,yv,uell,colors='white',levels=linspace(0,uell.max(),9),lw=1.5)

ax.set_xticks([])
ax.set_yticks([])

ax.set_aspect('equal')

pyplot.tight_layout()

pyplot.show(block=False)
pyplot.savefig('compare_contours.png',dpi=120,bbox_inches='tight')
