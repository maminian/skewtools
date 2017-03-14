double precision function linear_interp_1d(un,u,y,y0)
! Linear interpolation in 1d, with the accompanying index-locator
! function.
!
! Given an array size n of "exact" u values, defined on an interval of 
! y values, and interpolation point
! y0, returns the approximate u(y0) value using a simple linear interpolation.


implicit none
     ! Input vars
     integer                            :: un
     double precision, dimension(un)    :: u
     double precision, dimension(un)    :: y
     double precision                   :: y0
     
     ! Internal vars
     double precision                   :: u1,u2,p,q
     integer                            :: j,k
     double precision, parameter        :: one = 1.0d0
     
     ! Functions
     integer locate3
     
     ! If on a uniform grid, should use locate2 instead to reduce
     ! the lookup cost.
     k = locate3(un,y(1),y(un),y0)

     u1 = u(k)
     u2 = u(k+1)

     linear_interp_1d = u1 + (u2-u1)/(y(k+1)-y(k))*(y0 - y(k))

end function linear_interp_1d
!
! ----------------------------------------------------
!
integer function locate3(n,xl,xr,x0)
! I THINK WE CAN DO BETTEH!
! https://www.youtube.com/watch?v=NTpptLoUEk8
!
! Assuming the array x is a uniform mesh, we can get the 
! precise index with some modular arithmetic. xl and xr are 
! the lower and upper bounds of the array.
!
! Also assumes arrays start counting at 1.
implicit none

     ! Input
     double precision    :: xl,xr,x0
     integer             :: n
     ! Internal
     double precision    :: h

     h = (xr-xl)/(n-2)

     ! Add one because indexing starts at 1.
     locate3 = floor( ((x0-xl) - dmod(x0-xl,h)) / h ) + 1

end function locate3
