#!/usr/bin/python

from pylab import *
import matplotlib.pyplot as pyplot
from matplotlib import cm
import h5py

temp = h5py.File(sys.argv[1])
aratio = temp['Aratio'].value
Sk = temp['Skew'].value
temp.close()

pyplot.plot(aratio,Sk, markersize=4,marker='.');

hx = [0,0.3]
hy = [-2*sqrt(5)/7,-2*sqrt(5)/7]
pyplot.plot(hx,hy,linestyle='--',linewidth=2,color='red',label=r'$\frac{-2\sqrt{5}}{7}$')
hx = [0.533352,0.533352]
hy = [-0.2,0.08]
pyplot.plot(hx,hy,linestyle='--',linewidth=2,color='purple',label=r'$0.533352...$')

xlabel('Aspect ratio '+r'$\epsilon = \frac{H}{L}$')
ylabel('Geometric Skewness')
title('Geometric Skewness as a function of aspect ratio')
grid(True)
legend(loc='best')


# Wrap-up
show(block=False)
savefig('skew_vs_aratio.png',dpi=150)
