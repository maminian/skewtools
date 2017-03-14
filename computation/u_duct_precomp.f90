double precision function u_duct_precomp(y,z)
! Calculate the approximate value of the flow u(y,z), 
! having already precomputed the flow on a fine grid.
! Essentially this is a wrapper function for linear_interp_2d.
! All the arrays are held in the module mod_ductflow.
!

use mod_ductflow

implicit none

     double precision, intent(in)  :: y,z
     double precision linear_interp_2d

     u_duct_precomp = linear_interp_2d(ui, uj, u_precomp, ya, za, y, z)

end function u_duct_precomp
