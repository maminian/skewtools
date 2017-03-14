import numpy as np
from numpy import zeros,size,shape,sqrt,std,log,exp
from os import listdir
from os.path import isfile,join

import h5py
from matplotlib import pyplot

def meanLine(data):
# Takes a list of 1d arrays in and returns their pointwise average.
     
     meanline = zeros(size(data[0]))
     for line in data:
          meanline += line
     # end for
     
     meanline /= max(len(data),1)
     
     return meanline
     
# end def
     
def stdLines(data,nstdev):
# Takes a list of 1d arrays in and returns two arrays of +/- nstdev pointwise
# standard deviations.
     
     sampmean = meanLine(data)
     
     nt = size(data[0])
     nsamp = len(data)
     
     upper = zeros(nt)
     lower = zeros(nt)
     
     samps = zeros(nsamp)
     
     for i in range(nt):

          for j in range(nsamp):
#               stdev += (data[j][i]-sampmean[i])**2
               samps[j] = data[j][i]
          # end for
          stdev = std(samps)
          
          upper[i] = sampmean[i] + nstdev*stdev
          lower[i] = sampmean[i] - nstdev*stdev
     # end for
     
     return lower,upper
     
# end def

def averageDataFromFiles(file_list,fname,key_list=['Mean','Variance','Skewness','Kurtosis','Avgd_Mean','Avgd_Variance','Avgd_Skewness','Avgd_Kurtosis','Hist_centers','Hist_heights']):
# Takes a list of h5 files (for example, from generateH5FileList)
# and a list of keys, and averages across the keys specified.
# All other fields in the h5 file are assumed to be identical across 
# the entire file list, 
# and are copied over only from the first file in file_list.
#
# The result is saved in a new h5 file given by fname.
# (include the .h5).
#
# Unfortunately this can only save in the current directory for now.
     
     out_file = h5py.File(fname,'w')
     
     # Save relevant datasets and initialize the array to be averaged.
     fref = h5py.File(file_list[0],'r')
     to_avg = {}
          
     other_keys = fref.keys()
     
     for key in key_list:
          to_avg[key] = zeros(shape(fref[key].value))
          other_keys.remove(key)
     # end for
     

     
     for okey in other_keys:
          # This will probably break if there are arrays of more than one dimension!
          okey_container = out_file.create_dataset(fref[okey].name, fref[okey].shape, fref[okey].dtype)
          
          okey_container.write_direct(fref[okey].value)
     # end for
     
     idx = 0
     for files in file_list:
          tempfile = h5py.File(files,'r')
          
          for key in key_list:
               to_avg[key] += tempfile[key].value
          # end for
          
          tempfile.close()
          
          idx += 1
     # end for
     
     for key in key_list:
          to_avg[key] /= idx
          avg_container = out_file.create_dataset(fref[key].name, fref[key].shape, fref[key].dtype)
          avg_container.write_direct(to_avg[key])
     # end for
     
     fref.close()
     out_file.close()
     
# end def

def gatherDataFromFiles(file_list,key):
# Takes a list of our h5 files, spits out a list of arrays 
# for the specified key ('Skew', for instance, or 'X' or 'Z', 
# if positions have been saved).

     data = []
     
     # FILE_LIST SHOULD BE A LIST, NOT A STRING.
     for hfile in file_list:

          h5obj = h5py.File(hfile)
          data.append(h5obj[key].value)
          h5obj.close()
     # end for

     return data

# end def

def generateH5FileList(directory):
# Takes in a (relative) directory, spits out a list of (relative)
# h5 files located in that directory.

     elements = listdir(directory)
     
     outlist = []
     
     for element in elements:
          thing = join(directory,element)
          if (isfile(thing) & (thing[-3:]=='.h5')):
               outlist.append(thing)
          # end if
     # end for

     outlist.sort()

     return outlist
# end def

def constrainedMinimum(x,y,xmin,xmax):
# Takes a pair of arrays x,y, and finds the smallest value
# of y in the interval [xmin,xmax], along with the index.

     ixmin = locateIndex(x,xmin)
     ixmax = locateIndex(x,xmax)
     
     
     yval = min(y[ixmin:ixmax])
     iyval = np.where(y==yval)[0][0]
     
     return yval,iyval

# end def

def locateIndex(x,x0):
# Returns the index j of array x which satisfies 
# x[j] <= x0 < x[j+1]. Assumes x is monotone. Adapted from
# Numerical Recipes. Also known as binary search. Similar
# in spirit to bisection.
#

     n = len(x)
     jl=0
     ju=n
     x_is_incr = (x[-1] >= x[0])
     
     while (ju-jl > 1):
          jm = (ju+jl)/2
          
          # Cute equivalence to account for both monotone
          # increasing and monotone decreasing cases.
          if (x_is_incr == (x0 > x[jm])):
               jl = jm
          else:
               ju = jm
          # end if
     # end while
     
     return max(0,jl-1)
     
# end def

def compare_functional_diff(t,t1,y1,t2,y2):
# Given two pairs of arrays t1,y1 and t2,y2, evaluate 
# the difference of them on a common array t, and 
# return the array of differences.
#

     nt = len(t)
     errs = zeros(nt)
     
     for i in range(nt):
          y1v = interpolate_1d_fn(t[i],t1,y1)
          y2v = interpolate_1d_fn(t[i],t2,y2)
          
          errs[i] = y2v-y1v
     # end for
     
     return errs
# end def

def interpolate_1d_fn(tval,t,y):
# Given arrays t and y, return an interpolated yval
# at requested tval. Linear interpolation.

     idx1 = locateIndex(t,tval)
     idx2 = idx1+1
     
     m = (y[idx2] - y[idx1])/(t[idx2] - t[idx1])
     yval = y[idx1] + m*(tval - t[idx1])
     
     return yval
     
# end def

def power_law(x,y,x1,x2):
# From two arrays x and y assumed to have the same length,
# and two values x1, x2 assumed to be contained in x, 
# get the slope of the line in log-log space,
# which would indicate a power law in normal space.
#
# For instance, we would expect an output of 2 if x and y 
# were related as y=x^2.
#

     i1 = locateIndex(x,x1)
     i2 = locateIndex(x,x2)
     
     m = (log(y[i2])-log(y[i1]))/(log(x[i2])-log(x[i1]))
     b = exp(y[i1] + m*(0.-x[i1]))
     
     return b,m

# end def

def panels_1d(m,n,plotcmd,t,x,A,ax=None,fig=None,**kwargs):
#
# Produces a "pretty" m-by-n paneled result of a plot-like command 
# of the input arrays of the form (x,A[i,:]).
# Should work for plot, scatter, hist, possibly others.
# Output is the figure and axes, which can be edited before 
# pyplot.show()ing. Intended to display things over time.
#
# t and A are expected to have the same first dimension, 
# and the arguments of A[i,...] are commensurate for plotcmd.
#
# It is also assumed that m*n < size(t).
#
     import numpy
     if (ax==None & fig==None):
          fig,ax = pyplot.subplots(m,n)
     # end if
	 
     nt = len(t)

     idxs = numpy.arange(0,nt,(float(nt)/(m*n)))

     for j in range(n):
          for i in range(m):

               idx = int( idxs[j + n*i] )

               pyplot.sca(ax[i,j])
               if (plotcmd == pyplot.hist):
                    plotcmd(A[idx],x,**kwargs)
               else:
                    plotcmd(x,A[idx],**kwargs)
               # end if
               
               ax[i,j].set_xticklabels([])
               ax[i,j].set_yticklabels([])
               ax[i,j].set_title(r'$t=%.3g$'%t[idx])
          # end for
     # end for
     pyplot.subplots_adjust(wspace = 0.3, hspace = 0.3)

     return fig,ax

# end def

def panels_hist2d(t,tidxs,m,n,X,Y,xsc=None,ysc=1.,**kwargs):
     """Produces an m-by-n paneled result of a hist2d command
     
     of the input arrays of the form (t,X,Y,A), 
     where panels show (roughly)
     
     hist2d(X,Y,A[i,:,:],**kwargs)
     
     for i=0,...,len(tidxs), where A[i,:,:] is the evaluated at t[i].
     Output is the figure and axes, which can be further 
     modified before pyplot.show()ing.
     
     No input checking is done; the arrays should be commensurate 
     for the process described above.
     
     Must have len(tidxs)<=m*n. If len(tidxs)<m*n, 
     it will fill up the panels by row until all 
     the time indices are used up.
     """
     
     import numpy
     fig,ax = pyplot.subplots(m,n)

     for i in range(m):
          for j in range(n):
               if ( j+n*i <= (len(tidxs)-1) ):
                    idx = tidxs[j + n*i]

                    pyplot.sca(ax[i,j])

                    if True:
                         xl = min(X[:,idx])
                         xu = max(X[:,idx])
                         xs = max(abs(xl),abs(xu))
                         
                         ax[i,j].hist2d(X[:,idx],Y[:,idx],range=[[-xs,xs],[-ysc,ysc]],**kwargs)
                    else:
                         ax[i,j].hist2d(X[:,idx],Y[:,idx],**kwargs)
                    # end if                    
                    ax[i,j].set_xticklabels([])
                    ax[i,j].set_yticklabels([])
               else:
                    ax[i,j].set_frame_on(False)
                    ax[i,j].set_xticks([])
                    ax[i,j].set_yticks([])
               # end if
          # end for
     # end for
     
     pyplot.subplots_adjust(wspace = 0.3, hspace = 0.3)
     return fig,ax

# end def

def add_prefix(mlist,prefix):
	# For an input list of strings, adds the specified prefix to each.
	for i in range(len(mlist)):
		mlist[i] = prefix+mlist[i]
	# end for
	return mlist
# end def

def add_suffix(mlist,suffix):
	# For an input list of strings, adds the specified suffix to each.
	for i in range(len(mlist)):
		mlist[i] = mlist[i]+suffix
	# end for
	return mlist
# end def

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

def plotDataFromH5List(flist,xax,yax,fig=None,ax=None,**kwargs):
     # Plots everything on flist on a common axis; passes **kwargs in 
     # to the plot command (e.g., linewidth, color, etc). 
     # Returns the figure and axis in case you want to edit further, 
     # or do another call to plot another set of data.
     
     import matplotlib.pyplot as pyplot
     
     xvals = gatherDataFromFiles(flist,xax)
     yvals = gatherDataFromFiles(flist,yax)
     
     if ((fig == None) & (ax == None)):
          fig,ax = pyplot.subplots(1,1)
     # end if
     
     ax.hold(True)
     for i in range(len(flist)):
          ax.plot(xvals[i],yvals[i],**kwargs)
     # end for
     
     return fig,ax
# end def

def createTstepFile(filename,tvals,dsetname='target_times'):
     # Creates an hdf file containing by default a dataset
     # named 'target_times' containing the values tvals.
     #
     # Will write over anything previously in filename.
     
     mf = h5py.File(filename,'w')
     mf.create_dataset(dsetname,data=tvals)
     mf.close()
     
     return 0
# end def

def importDatasets(fname,*args):
     # Takes in the location of an h5 file,
     # and an arbitrary number of names of datasets in,
     # reads them all, passes them to output.
     # 
     # e.g., 
     #
     # t,mmean,mvar,msk = importDatasets('myfile.h5','Time','Mean','Variance','Skewness')
     #
     # looks in 'myfile.h5' for the appropriate datasets and puts them 
     # in the corresponding output variables.
     #
     
     mf = h5py.File(fname,'r')
     out = []
     
     for dset in args:
          out.append(mf[dset].value)
     # end for

     return out
# end def

def m2_channel(t,Pe=10.**4,nmax=1000,cutoff=10.**-3):
     # Calculates values for m2, the second moment of the averaged distribution
     # in the channel. Switches from short time asymptotics to series solution
     # at a hard-coded value.
     from numpy import argmin,zeros,shape,pi,sqrt
     
     cidx = argmin(abs(t-cutoff))
     
     m2 = zeros(shape(t))
     te = t[:cidx]
     tl = t[cidx:]
     
     # Early time component, Poisson resummed
     m2[:cidx] = 2.*te + Pe**2*(4./45.*te**2 - 4./9.*te**3 + 128./(105.*sqrt(pi))*te**3.5 - 1./3.*te**4)
     
     # Intermediate/late time component.
     for i in range(1,nmax):
          m2[cidx:] += (i*pi)**-8 * exp(-(i*pi)**2*tl)
     # end for
     
     m2[cidx:] *= 16.*Pe**2
     m2[cidx:] += -8./4725.*Pe**2 + 2*(1. + 8./945.*Pe**2)*tl 
     
     return m2
# end def

def m3_channel(t,Pe=10.**4,nmax=1000,cutoff=10.**-3):
     # Calculates values for m3, the third moment of the averaged distribution
     # in the channel. Switches from short time asymptotics to series solution
     # at a hard-coded value.
     from numpy import argmin,zeros,shape,pi,sqrt
          
     cidx = argmin(abs(t-cutoff))
     
     m3 = zeros(shape(t))
     te = t[:cidx]
     tl = t[cidx:]
     
     # Early time component, Poisson resummed
     m3[:cidx] = Pe**3*(-16./945.*te**3 + 16./45.*te**4 - 256./(105.*sqrt(pi))*te**4.5 + 5./2.*te**5 - 14464./(3465.*sqrt(pi))*te**5.5 + 14./15.*te**6)
     
     # Intermediate/late time component.
     for i in range(1,nmax):
          factor = -1488./(i*pi)**12 + 144./(i*pi)**10 - 24.*tl/(i*pi)**10
          m3[cidx:] += factor * Pe**3 * exp(-(i*pi)**2*tl)
     # end for
     
     m3[cidx:] += Pe**3*(-64./155925.*tl + 1376./19348875.)     

     return m3
# end def

def hk(tv,y):
     '''
     hk(tv,y). Assumes a scalar input tv and array input y. Returns 
     the one dimensional heat kernel, G(y,t) = exp(-y**2/(4t))/sqrt(4 pi t), 
     of the same size as y. The function G satisfies the free space problem
     
     G_t - G_yy = 0, G(y,t) = delta(y), delta being the Dirac delta function.

     '''
     from numpy import exp,pi,sqrt

     return exp(-y**2/(4.*tv))/sqrt(4.*pi*tv)
#

def p41(tv,d):
     '''
     p41(tv,d). Assumes scalar tv, array y. Returns a particular 
     polynomial in 
     tv and d>0 which is part of the Poisson resummation of the 
     channel first moment solution.
     
     d is understood to be abs(y-y0), for some scalar y0.
     '''
     
     return -1./3.*(d**2*tv + 4*tv**2)
#

def p42(tv,d):
     '''
     p42(tv,d). Assumes scalar tv, array y. Returns a particular 
     polynomial in 
     tv and d>0 which is part of the Poisson resummation of the 
     channel first moment solution.
     
     d is understood to be abs(y-y0), for some scalar y0.
     '''
     
     return 1./12.*(d**3 + 6*d*tv)
#

def c1_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2):
     ''' 
      c1_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2)
      Calculates values for c1, the pointwise first moment 
      in the channel. Switches from short time asymptotics to series solution
      at the cutoff point. 

      t and y are assumed to be one dimensional arrays. The 
      output c1 is a two dimensional array of shape (len(y),len(t)).
      
     '''
     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange
     from scipy.special import erfc
     
     cidx = argmin(abs(t-cutoff))

     if (len(shape(t))==0):

          if (t > cutoff):
               c1 = zeros( len(y) )
               
               for k in range(1,nmax):
#                    print k,(-1)**k/((k*pi)**4)*exp(-(k*pi)**2*t)
                    c1 += (-1)**k/((k*pi)**4)*exp(-(k*pi)**2*t)*cos(k*pi*y)
               # end for
               c1 *= 4.
               c1 += (y**4/2.-y**2+7./30.)/6.
          else:
               c1 = t*(1./3.-y**2) - t**2
               for m in arange(-3,4,2):
                    c1 += -4*p41(t,abs(y-m))*hk(t,y-m)
                    c1 += -4*p42(t,abs(y-m))*erfc(abs(y-m)/sqrt(4.*t))
               #
          #
     else:
          
          c1 = zeros( (len(y),len(t)) )
          
          te = t[:cidx]
          tl = t[cidx:]
          
          # Early time component, Poisson resummed
     #     m2[:cidx] = 2.*te + Pe**2*(4./45.*te**2 - 4./9.*te**3 + 128./(105.*sqrt(pi))*te**3.5 - 1./3.*te**4)
          
          # Intermediate/late time component.
          for j in range(cidx,len(t)):
               for i in range(1,nmax):
                    c1[:,j] += (-1)**i/(i**4)*exp(-(i*pi)**2*t[j])*cos(i*pi*y)
               # end for
               c1[:,j] *= 4./(pi**4)
               c1[:,j] += (y**4/2.-y**2+7./30.)/6.
          # end for

          # Early time component. Can get away with only using a few images!
          for j in range(cidx):
               c1[:,j] = t[j]*(1./3.-y**2) - t[j]**2 - 1./6.*(y**4/2. - y**2 + 7./30.)
               for m in arange(-3,4,2):
                    c1[:,j] -= p41(t[j],abs(y-m))*hk(t[j],abs(y-m))
                    c1[:,j] -= p42(t[j],abs(y-m))*erfc(abs(y-m)/sqrt(4.*t[j]))
               #
          #

     #


     c1 *= Pe
     return c1
# end def

def c2_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2):
     ''' 
      c2_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2)
      Calculates values for c2, the pointwise central second moment 
      in the channel. Switches from short time asymptotics to series solution
      at the cutoff point. 

      t and y are assumed to be one dimensional arrays. The 
      output c2 is a two dimensional array of shape (len(y),len(t)).
      
     '''
     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange
     from scipy.special import erfc
     
     cidx = argmin(abs(t-cutoff))

     if (len(shape(t))==0):
          c2 = zeros( len(y) )
          if (t > cutoff):
               for k in range(1,nmax):
                    c2 += chanQ2(t,y,k)*cos(k*pi*y) + chanQ3(t,y,k)*sin(k*pi*y)
               #
               c2 += chanQ1(t,y)
          else:
               print 'what'
               print 
               c2 += 0.
          #
          c2 *= Pe**2
          c2 += 2*t
     else:
          print 'what'
          c2 = zeros( (len(y),len(t)) )
          
          for j in range(cidx,len(t)):
               for k in range(1,nmax):
                    c2[:,j] += chanQ2(t[j],y,k)*cos(k*pi*y) + chanQ3(t[j],y,k)*sin(k*pi*y)
               #
               c2[:,j] += chanQ1(t[j],y)
               c2[:,j] *= Pe**2
               c2[:,j] += 2*t[j]
          #
          for j in range(cidx):
               c2[:,j] += 0.

               c2 *= Pe**2
               c2 += 2*t
          #
     #


     return c2
# end def

def chanQ1(t,y):
     '''
     Outputs the Q1 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''
     
     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange
     
     Q1 = (-413 + 3840*t - 1020*y**2 + 3570*y**4 -2940*y**6 + 675*y**8)/226800.
     
     return Q1
#

def chanQ2(t,y,n):
     '''
     Outputs the Q2 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     prefactor = 2*(-1)**n/((n*pi)**6) * exp(-(n*pi)**2*t)
     
     Q2 = prefactor*(17./3 - 64./((n*pi)**2) -2*t + y**2)
     
     return Q2
#

def chanQ3(t,y,n):
     '''
     Outputs the Q3 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     prefactor = 4*y*((-1)**n)/((n*pi)**5) * exp(-(n*pi)**2*t)
     Q3 = prefactor*(-1./((n*pi)**2) + 1./3*(y**2-1.))
     
     return Q3
#

def c3_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2):
     ''' 
      c3_channel(t,y,Pe=10.**4,nmax=1000,cutoff=10.**-2)
      Calculates values for c3, the pointwise central second moment 
      in the channel. Switches from short time asymptotics to series solution
      at the cutoff point. 

      t and y are assumed to be one dimensional arrays. The 
      output c3 is a two dimensional array of shape (len(y),len(t)).
      
     '''
     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange
     from scipy.special import erfc
     
     cidx = argmin(abs(t-cutoff))

     if (len(shape(t))==0):
          c3 = zeros( len(y) )
          if (t > cutoff):
               for k in range(1,nmax):
                    c3 += Pe*chanR2(t,y,k)*cos(k*pi*y)
                    c3 += Pe**3*chanR3(t,y,k)*cos(k*pi*y)
                    c3 += Pe**3*chanR4(t,y,k)*sin(k*pi*y)
               #
               c3 += Pe*chanR1(t,y) + Pe**3*chanR5(t,y)
          else:
               c3 += 0.
          #

     else:
          c3 = zeros( (len(y),len(t)) )
          
          for j in range(cidx,len(t)):
               for k in range(1,nmax):
                    c3 += 0.
               #
          #
          else:
               c3 += 0.
          #
          for j in range(cidx):
               c3[:,j] += 0.
          #
     #


     return c3
# end def

def chanR1(t,y):
     '''
     Outputs the R1 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''
     
     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange
     
     R1=1./30*t*(7. - 30.*y**2 + 15.*y**4)
     
     return R1
#

def chanR2(t,y,n):
     '''
     Outputs the R2 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     R2=24.*(-1)**n/((n*pi)**4)*t*exp(-(n*pi)**2*t)
     
     return R2
#
def chanR3(t,y,n):
     '''
     Outputs the R3 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     R3 = 14640./((n*pi)**12) - 3705./(2.*(n*pi)**10) + 691./(20*(n*pi)**8) + 3.*t**2/((n*pi)**8)
     R3+= y**2*(-231./(2.*(n*pi)**10)+9./(2.*(n*pi)**8)-1./(3.*(n*pi)**6))
     R3+= y**4*(23./(4.*(n*pi)**8) + 2./(3.*(n*pi)**6))
     R3+= -y**6/(3.*(n*pi)**6)
     R3+= t*(231./((n*pi)**10) - 31./((n*pi)**8) - 8./(15.*(n*pi)**6) - 3.*y**2/((n*pi)**8))
     
     
     R3*=(-1)**n*exp(-(n*pi)**2*t)
     
     return R3
#
def chanR4(t,y,n):
     '''
     Outputs the R4 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     R4 = 231./((n*pi)**8)+44./((n*pi)**6) -28./(5.*(n*pi)**4)
     R4+= y**2*(-76./((n*pi)**6) + 4./((n*pi)**4))
     R4+= 8*y**4/(5*(n*pi)**4)
     R4+= t*(6./((n*pi)**6) + 2/((n*pi)**4) -2*y**2/((n*pi)**4))
     
     R4*= (-1)**n*exp(-(n*pi)**2*t)*y/(n*pi)**3
     return R4
#
def chanR5(t,y):
     '''
     Outputs the R5 coefficient needed in the channel second moment formula, 
     detailed in the PRL paper.
     
     ASSUMES t IS A SCALAR!
     '''

     from numpy import argmin,zeros,shape,pi,sqrt,sin,cos,abs,arange,exp
     
     R5 = -4076777./13621608000. + 8447./4989600.*y**2 - 713./907200.*y**4
     R5 += -1./1200.*y**6 + 13./12096.*y**8 - 211./453600.*y**10 + 1./14784.*y**12
     R5 += 244./155925.*t - 8./945.*t*y**2 + 4./945.*t*y**4
     
     return R5
#

def centered_to_noncentered(mu,var,skew):
     # Given the centered statistics mean, var, skew, 
     # return the corresponding triple
     # (m1,m2,m3).

     m1 = mu
     m2 = var + m1**2
     m3 = skew*var**1.5 + 3.*m1*m2 - 2.*m1**3

     return m1,m2,m3

# end def

def plotSummarizedDataFromH5List(flist,xax,yax,fig=None,ax=None,myc=[1,0,0,1]):
     # Functionally similar to plotDataFromH5List, 
     # but assumes the files all share the same X data.
     # Plots the mean +/- 2 standard deviations of the y data.
     #
     # Returns the figure and axis in case you want to edit further, 
     # or do another call to plot another set of data.
     
     import matplotlib.pyplot as pyplot
     
     xvals = gatherDataFromFiles([flist[0]],xax)[0]
     yvals = gatherDataFromFiles(flist,yax)
     
     ml = meanLine(yvals)
     psd,msd = stdLines(yvals,2)

     if ((fig == None) & (ax == None)):
          fig,ax = pyplot.subplots(1,1)
     # end if
     
     ax.plot(xvals,ml,color=myc,lw=2)

     myc2 = [myc[0],myc[1],myc[2],0.1]
     ax.fill_between(xvals,msd,psd,color=myc2)
     
     return fig,ax
# end def

def m2_pipe(t,Pe=10.**4,nmax=1000,m20=0.):
     # Calculates exact M2 in the pipe. Optional 
     # initial M2(t=0), m20 included.
     from scipy.special import jnp_zeros
     
     mu = jnp_zeros(0,nmax)
     Pe_fran = Pe/2.

     m2 = m20 + 2.*t + Pe_fran**2*(-1./1440. + 1./96.*t)

     for n in range(nmax):
          m2 += Pe_fran**2 * 32.*(exp(-(mu[n]**2)*t)/(mu[n]**8))
     # end for

     return m2
# end def

def m3_pipe(t,Pe=10.**4,nmax=1000,m30=0.):
     # Calculates exact M3 in the pipe. Optional 
     # initial M3(t=0), m30 included.

     from scipy.special import jnp_zeros
     
     mu = jnp_zeros(0,nmax)
     Pe_fran = Pe/2.
     
     m3 = 1./480.*(-17./896.+112./896.*t)
     for n in range(nmax):
          m3 += 16.*exp(-(mu[n]**2)*t)*(-240./(mu[n]**12) + 18./(mu[n]**10) + t/(mu[n]**8))
     # end for
     
     m3 *= Pe_fran**3

     return m3
# end def

def c1_pipe(r,t,Pe=10.**4,nmax=1000):
     # Calculates exact C1(r,t) in the pipe. Assumes t scalar.
     
     from scipy.special import jnp_zeros,j0
     mu = jnp_zeros(0,nmax)
     Pe_fran = Pe/2.
     
     c1 = zeros(shape(r))
     
     for n in range(nmax):
          c1 += 1./(j0(mu[n])*mu[n]**4)*(1.-exp(-mu[n]**2*t)*j0(mu[n]*r))
     #
     c1 *= -4*Pe_fran
     
     return c1
#

def c2_pipe(r,t,Pe=10.**4,nmax=1000):
     # Calculates exact C2(r,t) in the pipe. Assumes t scalar.
     
     from scipy.special import jnp_zeros,j0,j1
     mu = jnp_zeros(0,nmax)
     Pe_fran = Pe/2.
     
     c2 = zeros(shape(r))
     
     for n in range(nmax):
          c2 += (-24.+6*mu[n]**2+72*mu[n]**2*t-9*mu[n]**4*t+6*mu[n]**4*r**2*t)/(3*mu[n]**8*j0(mu[n]))*j0(mu[n]*r)
          c2 += (24-24*mu[n]**2*t)/(3*mu[n]**7*j0(mu[n]))*r*j1(mu[n]*r)
          
          expterm = (-120+9*mu[n]**2+mu[n]**2*r**2 + mu[n]**4*t)/(3*mu[n]**8*j0(mu[n]))*j0(mu[n]*r)
          expterm += (-2. -mu[n]**2 + mu[n]**2*r**2)/(3*mu[n]**7*j0(mu[n]))*r*j1(mu[n]*r)
          c2 += exp(-mu[n]**2*t)*expterm
     #
     c2 *= 4
     c2 += 1./5760.*(29.-240.*t -60*r**2 + 15*r**4 + 240*r**2*t)
     c2 *= Pe_fran**2
     #
     c2 += 2*t 
     
     return c2
#

def duct_meanspeed(aratio,kmax=1000):
     # Calculates the mean speed 
     # for the flow with laplacian -2.
     val = 2./3.
     for k in range(1,kmax):
          val += -4.*aratio/(((k-0.5)*pi)**5)*tanh((k-0.5)*pi/aratio)
     # end for

     return val
# end def

def ellipse_meanspeed(aratio):
     # Calculates the mean speed
     # for the ellipse flow with laplacian -2.
     return 0.25
# end def
