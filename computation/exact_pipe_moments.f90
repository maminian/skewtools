! Functions using Francesca's exact formulas for the moments 
! of the diffusive scalar in the pipe.

double precision function exact_pipe_m1(nTerms,Pe,t)
implicit none

     integer             :: nTerms     
     double precision    :: Pe,t
     
     exact_pipe_m1 = 0.0d0
     
end function exact_pipe_m1
!
! --------------------------------
!
double precision function exact_pipe_m2(nTerms,Pe,t)
implicit none
! Be careful; we've done a pre-step where the Peclet number 
! cancels on numerator and denominator when calculating skewness.
!
! The result here will NOT be correct if we only seek m2.

     integer                       :: nTerms     
     double precision              :: Pe,t
     
     double precision, parameter   :: pi = 4.0d0*datan(1.0d0)
     
     exact_pipe_m2 = 0.0d0

end function exact_pipe_m2
!
! --------------------------------
!
double precision function exact_pipe_m3(nTerms,Pe,t)
implicit none
! Be careful; we've done a pre-step where the Peclet number 
! cancels on numerator and denominator when calculating skewness.
!
! The result here will NOT be correct if we only seek m3.

     integer                       :: nTerms     
     double precision              :: Pe,t
     
     double precision, parameter   :: pi = 4.0d0*datan(1.0d0)

     exact_pipe_m3 = 0.0d0
     
end function exact_pipe_m3
!
! --------------------------------
!
double precision function exact_pipe_c1(nTerms,Pe,t,y)
implicit none

     integer             :: nTerms
     double precision    :: Pe,t,y
     
     exact_pipe_c1 = 0.0d0
     
end function exact_pipe_c1
!
! --------------------------------
!
double precision function exact_pipe_c2(nTerms,Pe,t,y)
implicit none

     integer             :: nTerms     
     double precision    :: Pe,t,y
     
     exact_pipe_c2 = 0.0d0
     
end function exact_pipe_c2
!
! --------------------------------
!
double precision function exact_pipe_c3(nTerms,Pe,t,y)
implicit none

     integer             :: nTerms     
     double precision    :: Pe,t,y
     
     exact_pipe_c3 = 0.0d0
     
end function exact_pipe_c3
!
! --------------------------------
!
double precision function asymp_pipe_m1(Pe,t)
implicit none

     double precision    :: Pe,t
     
     asymp_pipe_m1 = 0.0d0
end function asymp_pipe_m1
!
! --------------------------------
!
double precision function asymp_pipe_m2(Pe,t)
implicit none

     double precision    :: Pe,t
     
     asymp_pipe_m2 = 2.0d0*t + Pe**2*( 1.0d0/12.0d0*t**2 -2.0d0/3.0d0*t**3 )
end function asymp_pipe_m2
!
! --------------------------------
!
double precision function asymp_pipe_m3(Pe,t)
implicit none

     double precision    :: Pe,t
     
     asymp_pipe_m3 = Pe**3*( 1.0d0/3.0d0*t**4 )
end function asymp_pipe_m3
