subroutine make_histogram(n,X,nhb,centers,heights)
! Bins the array X with specified bin centers. This subroutine uses 
! equally spaced bins.
!
! In:
!
!    n
!    X(n)     - array of particle positions to bin across the second 
!                          dimension, then averaged across the first dimension.
!    nhb                 - number of histogram bins
!
! In/out:
!
!    centers(nhb), heights(nhb)    - centers and heights of the bins. Normalized
!                                    to be a probability density.
! 


implicit none

     integer, intent(in)                                    :: n,nhb
     double precision, dimension(n), intent(in)             :: X
     double precision, dimension(nhb), intent(out)          :: centers,heights


     integer, dimension(n)                                  :: Xidx
     double precision                                       :: xmin,xmax,db
     integer                                                :: j

     ! Set up the binning.
     xmin = minval(X)
     xmax = maxval(X)
!     write(*,*) xmin,xmax
     if (xmin .eq. xmax) then
          xmin = xmin - 1.0d0
          xmax = xmax + 1.0d0
     end if


     db = (xmax-xmin)/nhb

     do j=1,nhb
          centers(j) = xmin + (j-0.5d0)*db
          heights(j) = 0.0d0
     end do


     call uniform_bins_idx(n,X,xmin,xmax,nhb,Xidx)


     do j=1,n
!          write(*,*) j,X(j),nhb,Xidx(j)
          heights(Xidx(j)) = heights(Xidx(j)) + 1.0d0
     end do


!     heights = heights/sum(heights)

end subroutine make_histogram
