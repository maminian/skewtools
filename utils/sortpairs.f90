subroutine sortpairs(nTot,X,Xdup,bin_idxs,nbins)
!
! Given a pair of arrays X, Xdup (double), bin_idxs (integer), 
! sorts (buckets) the bin_idxs and carries along the associated 
! Xdup value. bin_idxs is known to have integer values from 1 to nbins.
!
! This is a "partial" sorting technically, since we don't need 
! the Xdup values to be sorted in any way.
! 
! Accomplished by doing a first pass of bin_idxs to see 
! how many of each integer there are, then a second pass of 
! copying over the contents of X into Xdup in an appropriate order.
!

implicit none
     integer, intent(in)                               :: nTot,nbins
     double precision, dimension(nTot), intent(in)     :: X

     double precision, dimension(nTot), intent(out)    :: Xdup
     integer, dimension(nTot), intent(out)             :: bin_idxs

     integer, dimension(nTot)                          :: bin_idxs2
     integer, dimension(nbins)                         :: binTots,binCurr
     integer, dimension(nbins+1)                       :: binCum
     integer                                           :: i,idx,ptr

     do i=1,nbins
          binTots(i) = 0
     end do

     !
     ! Do a count of the number of things in each bin.
     !

     do i=1,nTot
          ptr = bin_idxs(i)
          binTots(ptr) = binTots(ptr) + 1
     end do

     !
     ! Generate a cumulative count as pointers.
     ! Makes it easy to reference the start and end elements
     ! of the subset.
     !
     binCum(1) = 1
     do i=2,nbins+1
          binCum(i) = binCum(i-1) + binTots(i-1)
          binCurr(i-1) = binCum(i-1)
     end do
     
     bin_idxs2 = bin_idxs
     !
     ! Now do a second loop, placing values of X into Xdup
     ! and updating the pointer binCurr along the way.
     !

     do i=1,nTot
          idx = bin_idxs2(i)

          ptr = binCurr(idx)
          binCurr(idx) = binCurr(idx) + 1

          Xdup(ptr) = X(i)
          bin_idxs(ptr) = bin_idxs2(i)
     end do
     

end subroutine sortpairs
