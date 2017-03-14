subroutine precompute_uvals(ui,uj,u,ya,za,nTerms,idxlist,uij_vals)
! Precomputing array u, ya, and za, for use in bilinear interpolation
! to optimize function calls in the Monte Carlo iteration.
! Assumes we are working on the square [-1,1]x[-1,1].

implicit none
     ! Input
     integer                                 :: ui,uj,nTerms
     integer, dimension(nTerms)              :: idxlist
     double precision, dimension(nTerms)     :: uij_vals
     
     ! Input/Output
     double precision, dimension(ui,uj)      :: u
     double precision, dimension(ui)         :: ya
     double precision, dimension(uj)         :: za     
     
     ! Internal
     integer                                 :: i,j
     
     ! Functions
     double precision u_duct
     
     ! -------------------------------------
     
     ! First construct the y,z arrays.
     if (.false.) then
          call uniform_nodes(ui,ya,-1.0d0,1.0d0)
          call uniform_nodes(uj,za,-1.0d0,1.0d0)
     else
          ! Chebyshev nodes. Be aware you need to use the general
          ! index locator for linear interpolation if you use this,
          ! which will slow down function evaluations.
          call padded_cheb_nodes(ui,ya,-1.0d0,1.0d0)
          call padded_cheb_nodes(uj,za,-1.0d0,1.0d0)
     end if
     
     
     ! Now compute the corresponding u values.
     do i=1,ui
          do j=1,uj
               u(i,j) = u_duct(ya(i),za(j),nTerms,idxlist,uij_vals)
          end do
     end do

end subroutine precompute_uvals
