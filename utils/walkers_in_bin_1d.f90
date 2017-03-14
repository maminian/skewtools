subroutine walkers_in_bin_1d(nGates,X,Y,bin_lo,bin_hi,X_bin,bin_count)
! Takes arrays X,Z, collects all indices bin_lo < Z(i) < bin_hi
! and saves them sequentially in X_bin.
! 
! The actual number of relevant values in X_bin 
! is unknown a priori, but at most nGates. Hence we keep track of 
! the actual number of X positions in a bin with bin_count.
!

implicit none

     integer, intent(in)                               :: nGates
     double precision, dimension(nGates), intent(in)   :: X,Y
     double precision, intent(in)                      :: bin_lo,bin_hi
     
     double precision, dimension(nGates), intent(out)  :: X_bin
     integer, intent(inout)                            :: bin_count

     integer                                           :: i
     
     bin_count=0
     do i=1,nGates
          if ((bin_lo < Y(i)) .and. (Y(i) < bin_hi)) then
               bin_count = bin_count + 1
               X_bin(bin_count) = X(i)
          end if     
     end do

end subroutine walkers_in_bin_1d
