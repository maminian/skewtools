! Functions for the short-time asymptotics of the moments in the channel. 

double precision function asymp_st_channel_m2(Pe,t)
implicit none
! Be careful; we've done a pre-division where the Peclet number 
! cancels on numerator and denominator when calculating skewness.
!
! The result here will NOT be correct if we only seek m2.
   
     double precision              :: Pe,t
     
     double precision, parameter   :: rtpi = dsqrt(4.0d0*datan(1.0d0))
     
     
     asymp_st_channel_m2 = 2.0d0*t + Pe**2*((4.0d0/45.0d0)*t**2 - (4.0d0/9.0d0)*t**3 &
                    +(128.0d0/(105.0d0*rtpi))*t**3.5d0 -(1.0d0/3.0d0)*t**4)
     
     
end function asymp_st_channel_m2

double precision function asymp_st_channel_m3(Pe,t)
implicit none
! Be careful; we've done a pre-step where the Peclet number 
! cancels on numerator and denominator when calculating skewness.
!
! The result here will NOT be correct if we only seek m3.
   
     double precision              :: Pe,t
     
     double precision, parameter   :: rtpi = dsqrt(4.0d0*datan(1.0d0))
     
     
     asymp_st_channel_m3 = Pe**3*((-16.0d0/945.0d0)*t**3 + (16.0d0/45.0d0)*t**4 &
                    -(256.0d0/(105.0d0*rtpi))*t**4.5d0 + (5.0d0/2.0d0)*t**5 &
                    -(14464.0d0/(3465.0d0*rtpi))*t**5.5d0 + (14.0d0/15.0d0)*t**6)
     
end function asymp_st_channel_m3

