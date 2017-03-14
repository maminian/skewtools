subroutine precompute_uvals_ss(a,b,aratio)
! Precomputing array u, ya, and za, for use in bilinear interpolation
! to optimize function calls in the Monte Carlo iteration.
! Assumes we are working on the square [-1,1]x[-1,1].

use mod_ductflow
implicit none
     ! Input
     double precision, intent(in)                      :: a,b,aratio
        
     
     ! Internal
     integer                                           :: i,j
     double precision, dimension(ui)                   :: yatta
     double precision, dimension(uj)                   :: zatta
     double precision, dimension(ui,uj)                :: udumb
     
     ! Functions
     double precision u_duct_ss
     
     ! -------------------------------------
     
     ! First construct the y,z arrays.
     if (.true.) then
          call uniform_nodes(ui,yatta,-a,a)
          call uniform_nodes(uj,zatta,-b,b)
     else
          ! Chebyshev nodes. Be aware you need to use the general
          ! index locator for linear interpolation if you use this,
          ! which will slow down function evaluations.
          call padded_cheb_nodes(ui,ya,-a,a)
          call padded_cheb_nodes(uj,za,-b,b)
     end if
     
     
     ya = yatta
     za = zatta

     ! Now compute the corresponding u values.
     do i=1,ui
          do j=1,uj
               udumb(i,j) = u_duct_ss(ya(i),za(j),nTerms,aratio)
!               write(*,*) ya(i),za(j),u_duct_ss(ya(i),za(j),nTerms,aratio)
          end do
     end do

     u_precomp = udumb

end subroutine precompute_uvals_ss
