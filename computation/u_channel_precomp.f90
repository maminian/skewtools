double precision function u_channel_precomp(y,a,b)
! Calculate the approximate value of the flow u(y), 
! having already precomputed the flow on a fine grid.
! Essentially this is a wrapper function for linear_interp_1d.
! All the arrays are held in the module mod_ductflow.
!
! a and b are dummy variables here for the purpose of interfacing.

use mod_oscillatory_channel

implicit none

     double precision, intent(in)  :: y,a,b
     double precision linear_interp_1d

     u_channel_precomp = linear_interp_1d(ngo, mocflow, ygrid, y)

end function u_channel_precomp
