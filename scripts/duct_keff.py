from numpy import *
import matplotlib.pyplot as pyplot


def gameval(m,n,lam,kmax):
     gam = 0.
     if (m==0 and n==0):
          gam=0.
     elif (m==0 or n==0):
          for kmh in arange(0.5,kmax+1,1.):
               gam += tanh(kmh*pi/lam)/(kmh*(kmh**2 - m**2)*(kmh**2 + lam**2*n**2))
          # end for
          gam *= -16.*lam*(-1)**(m+n)/(2*pi**5)
     else:
          for kmh in arange(0.5,kmax+1,1.):
               gam += tanh(kmh*pi/lam)/(kmh*(kmh**2 - m**2)*(kmh**2 + lam**2*n**2))
          # end for
          gam *= -16.*lam*(-1)**(m+n)/(pi**5)
     # end if
     
     return gam
# end def

def coeff(lam,kmax):
     cval = 0.

     return 0.
# end def

###################################################

aratio = 1.0
nmax = 100
gamma = zeros( (nmax,nmax) )

# Calculate gamma values.
for m in range(nmax):
     for n in range(nmax):
          gamma[m,n] = gameval(m,n,aratio,2*nmax)
     # end for
# end for

term1 = 8./945.
term2 = 0.
term3 = 0.

for m in range(1,nmax):
     term2 += 2.*(-2.*(-1)**m*gamma[m,0]/(m*pi)**4)
# end for

for m in range(nmax):
     for n in range(nmax):
          if (m==0 and n==0):
               term3 += 0.
          elif (m==0 or n==0):
               term3 += 2.*gamma[m,n]**2/(4.*pi**2*(m**2+aratio**2*n**2))
          else:
               term3 += gamma[m,n]**2/(4.*pi**2*(m**2+aratio**2*n**2))
          # end if
     # end for
# end for

# Correct for fixing mean speed and cross sectional area

astc = sqrt(pi/4.)
pxstc = 
correction = astc**6*pxstc**2

keff = term1 + term2 + term3
keff *= correction

print "Aspect ratio: "+str(aratio)
print "Effective diffusivity boost: "+str(keff)
print "Circular pipe diffusiviy: "+str(1./768.)
