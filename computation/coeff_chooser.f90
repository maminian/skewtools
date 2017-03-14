subroutine coeff_chooser(nTerms,aratio,idx_list,uij_vals)
! Returns the (i,j) indices of the
! largest nTerms u_ij coefficients in decreasing order, with the purpose
! of speeding up convergence of relevant series for duct flow.
!
! Eventual goal is to input an arbitrary function as
! one of the arguments, but hard-code the coefficients 
! for now.
!

implicit none
     
     ! Inputs
     integer                                      :: nTerms
     double precision                             :: aratio
     integer, dimension(1:nTerms,1:2)             :: idx_list
     double precision, dimension(1:nTerms)        :: uij_vals
     
     ! Internal
     integer                            :: i,idx,next,nmax,itmp
     integer, dimension(:), allocatable :: bins
     double precision                   :: aval
     
     
     ! Bound on the number of bins necessary, assuming
     ! the indices favor m for aratio in [0,1], and
     ! the coefficient function is symmetric in m,n
     ! at aratio=1.
!     nmax = floor( dsqrt(dble(2*nTerms)) ) + 1

     ! Have reason to believe this might not be working
     nmax = nTerms
     
     allocate(bins(1:nmax))
     
     do i=1,nmax
          bins(i) = 1
     end do

     do idx=1,nTerms

          call get_largest_idx(bins,nmax,aratio,next,aval)
          itmp = next
          
          idx_list(idx,1) = bins(itmp)
          idx_list(idx,2) = next
          uij_vals(idx) = aval
          
          ! Push the boundary of the index set at next.
          
          bins(itmp) = bins(itmp) + 1
          
     end do
     
     deallocate(bins)
     
end subroutine coeff_chooser
!
! ------------------------------------------------------
!
subroutine get_largest_idx(bins,nmax,aratio,next,aval)
implicit none
     
     double precision              :: aratio, aval, acompare
     integer                       :: next,nmax,i
     integer, dimension(1:nmax)    :: bins
     
     double precision uij_exact

     next = 1
     aval = uij_exact(bins(next),next,aratio)
     do i=1,nmax
          acompare = uij_exact(bins(i),i,aratio)
          if (dabs(acompare) .gt. dabs(aval)) then
               next = i
               aval = acompare
          end if
     end do
     
end subroutine get_largest_idx
!
! ------------------------------------------------------
!
double precision function uij_exact(i,j,aratio)

implicit none
     
     ! Inputs
     integer             :: i,j
     double precision    :: aratio, iprod1, iprod2, ew
     double precision    :: twooverpi,half
     
     ! Functions
     double precision laplace_evalue
     
     parameter( twooverpi = 0.5d0/datan(1.0d0) )
     parameter( half = 0.5d0 )
     
     ! (phi_ij,-1) = -4/pi**2*(-1)**(i+j)/((i-1/2)(j-1/2))
     iprod1 = -twooverpi**2*(-1)**(i+j)/((i-half)*(j-half))

     ! (phi_ij,phi_ij) = 1
     iprod2 = 1.0d0    
     
     ! Eigenvalue of the nondimensionalized Laplacian.
     ew = laplace_evalue(i,j,aratio)
     
     uij_exact = iprod1/(ew*iprod2)

end function uij_exact
!
! -------------------------------
!
double precision function laplace_evalue(i,j,aratio)
! Eigenvalue for the nondimensionalized Laplace operator.
implicit none
     integer             :: i,j
     double precision    :: aratio,pi,negpisq,half
     parameter(pi = 4.0d0*datan(1.0d0))
     parameter(negpisq = -pi**2)
     parameter(half = 0.5d0)
     
     laplace_evalue = negpisq*( (i-half)**2 + aratio**2*(j-half)**2 )
     
end function laplace_evalue
