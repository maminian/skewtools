! Functions using Francesca's exact formulas for the moments 
! of the diffusive scalar in the channel.

double precision function exact_channel_m1(Nterms,Pe,t)
implicit none

     integer             :: Nterms     
     double precision    :: Pe,t
     
     exact_channel_m1 = 0.0d0
     
end function exact_channel_m1
!
! --------------------------------
!
double precision function exact_channel_m2(Nterms,Pe,t)
implicit none

     integer                       :: Nterms     
     double precision              :: Pe,t,npi
     
     integer                       :: n
     double precision, parameter   :: pi = 4.0d0*datan(1.0d0)
     
     exact_channel_m2 = 0.0d0
     
     do n=Nterms,1,-1
          npi = n*pi
          exact_channel_m2 = exact_channel_m2 + &
               (16.0d0)*dexp(-(npi**2)*t)/(npi**8)
     end do

     exact_channel_m2 = Pe**2*(exact_channel_m2) + 2.0d0*t + Pe**2*(-8.0d0/4725.0d0 + 16.0d0/945.0d0*t)
     
end function exact_channel_m2
!
! --------------------------------
!
double precision function exact_channel_m3(Nterms,Pe,t)
implicit none

     integer                       :: Nterms     
     double precision              :: Pe,t,npi
     
     integer                       :: n
     double precision, parameter   :: pi = 4.0d0*datan(1.0d0)

     exact_channel_m3 = 0.0d0
     
     do n=Nterms,1,-1
          npi = n*pi
          exact_channel_m3 = exact_channel_m3 + &
               ((-1488.0d0/(npi**2)+144.0d0-24.0d0*t)*dexp(-(npi**2)*t)/(npi**10))
     end do
     
     exact_channel_m3 = Pe**3*(exact_channel_m3 +(-64.0d0/155925.0d0)*t+(1376.0d0/19348875.0d0))
     
end function exact_channel_m3
!
! --------------------------------
!
double precision function exact_channel_m4(Nterms,Pe,t)
implicit none
! Be careful; we've done a pre-step where the Peclet number 
! cancels on numerator and denominator when calculating skewness.
!
! The result here will NOT be correct if we only seek m4.

     integer                       :: Nterms     
     double precision              :: Pe,t,npi
     
     integer                       :: n
     double precision, parameter   :: pi = 4.0d0*datan(1.0d0)

     exact_channel_m4 = 0.0d0
     
     do n=Nterms,1,-1
          npi = n*pi
          exact_channel_m4 = exact_channel_m4 + dexp(-npi**2*t) * ( &
               Pe**2*t*192.0d0/(npi**8) + Pe**4 * ( &
               276096.0d0/(npi**16) - 38560.0d0/(npi**14) + 15488.0d0/(15.0d0*npi**12) + &
               3288.0d0*t/(npi**14) - 400.0d0*t/(npi**12) - 64.0d0*t/(15.0d0*npi**10) + &
               24.0d0*t**2/(npi**12) &
               ) &
               )
          
     end do
     
     exact_channel_m4 = exact_channel_m4 + 12*t**2 + Pe**4*(&
               (7038848.0d0/162820783125.0d0) - (4352.0d0/14189175.0d0*t) + (256.0d0/297675.0d0*t**2)) + &
               Pe**2*(-(32.0d0/1575.0d0*t) + (64.0d0/315.0d0*t**2))
               
     
end function exact_channel_m4
!
! --------------------------------
!
double precision function exact_channel_c1(Nterms,Pe,t,y)
implicit none

     integer             :: Nterms     
     double precision    :: Pe,t,y
     
     exact_channel_c1 = 0.0d0
end function exact_channel_c1
!
! --------------------------------
!
double precision function exact_channel_c2(Nterms,Pe,t,y)
implicit none

integer             :: Nterms     
     double precision    :: Pe,t,y
     
     exact_channel_c2 = 0.0d0
end function exact_channel_c2
!
! --------------------------------
!
double precision function exact_channel_c3(Nterms,Pe,t,y)
implicit none

     integer             :: Nterms     
     double precision    :: Pe,t,y
     
     exact_channel_c3 = 0.0d0
end function exact_channel_c3
!
! --------------------------------
!
subroutine channel_full_moments_asymp(Pe,t,mean,vari,skew,kurt)
implicit none

     double precision, intent(in)  :: Pe,t
     double precision, intent(out) :: mean,vari,skew,kurt
     
     mean = 0.0d0

     vari = (2.0d0)*t + Pe**2*((4.0d0/45.0d0)*t**2 - (4.0d0/9.0d0)*t**3)

     skew = Pe**3*(-(16.0d0/945.0d0)*t**3 + (16.0d0/45.d00)*t**4)
     skew = skew/vari**1.5d0
     
     kurt = 12.0d0*t**2 + (16.0d0/15.0d0)*Pe**2*t**3 + (16.0d0/945.0d0)*Pe**4*t**4 - &
               (352.0d0/945.0d0)*Pe**4*t**5
               
     kurt = kurt/(vari**2)
     
end subroutine channel_full_moments_asymp
!
! ----------------------------------------------------------------------------------------
!
subroutine channel_full_moments_exact(Nterms,Pe,t,mean,vari,skew,kurt)
implicit none

     integer, intent(in)           :: Nterms
     double precision, intent(in)  :: Pe,t
     
     double precision, intent(out) :: mean,vari,skew,kurt
     
     double precision              :: m1,m2,m3,m4
     
     double precision exact_channel_m1,exact_channel_m2,exact_channel_m3,exact_channel_m4
     
     m1 = exact_channel_m1(Nterms,Pe,t)
     m2 = exact_channel_m2(Nterms,Pe,t)
     m3 = exact_channel_m3(Nterms,Pe,t)
     m4 = exact_channel_m4(Nterms,Pe,t)

     mean = m1
     vari = m2-mean**2
     skew = (m3-3.0d0*m1*m2+2.0d0*m1**3)/(vari)**1.5d0
     kurt = (m4-4.0d0*m1*m3+6.0d0*m1**2*m2-3.0d0*m1**4)/vari**2

     
end subroutine channel_full_moments_exact
