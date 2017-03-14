subroutine walkers_in_bin_2d(nTot,X,Y,Z,yl,yh,zl,zh,X_bin,bin_count)

! Takes arrays X,Y,Z, collects all indices satisfying 
! yl <= Y(i) < yh and
! zl <= Z(i) < zh, 
! and saves them sequentially in X_bin(1:bin_count).
! 
! The actual number of relevant values in X_bin 
! is unknown a priori, but at most nGates. Hence we keep track of 
! the actual number of X positions in a bin with bin_count.
!

implicit none

     integer, intent(in)                               :: nTot
     double precision, dimension(nTot), intent(in)     :: X,Y,Z
     double precision, intent(in)                      :: yl,yh,zl,zh
     
     double precision, dimension(nTot), intent(out)    :: X_bin
     integer, intent(out)                              :: bin_count

     ! Internal
     integer                                           :: i
     
     bin_count=0

     do i=1,nTot
          if ((yl <= Y(i)) .and. (Y(i) < yh) .and. (zl <= Z(i)) .and. (Z(i) < zh)) then
               bin_count = bin_count + 1
               X_bin(bin_count) = X(i)
          end if     
     end do

end subroutine walkers_in_bin_2d
