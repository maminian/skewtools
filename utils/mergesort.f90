recursive subroutine mergesort(n,p,r,x)
!
! An implementation of mergesort following pseudocode 
! at http://www.personal.kent.edu/~rmuhamma/Algorithms/MyAlgorithms/Sorting/mergeSort.htm
! 
! Input:
!    integer n
!    double precision, dimension(n), intent(in)   x
!

implicit none
     integer, intent(in)                               :: n,p,r
     double precision, dimension(n), intent(inout)     :: x
     
     integer q

     if (p .lt. r) then
          q = floor((p+r)/2.0d0)

          call mergesort(n,p,q,x)
          call mergesort(n,q+1,r,x)
          call mmerge(n,p,q,r,x)
     end if



end subroutine mergesort
!
! ---------------------------------------
!
subroutine mmerge(n,p,q,r,y)
implicit none
     integer, intent(in)                               :: n,p,q,r
     double precision, dimension(n), intent(inout)     :: y
     double precision, dimension(:), allocatable       :: l1,l2

     integer i,j,k
     integer n1,n2

!     integer, parameter                                :: n1 = q-p+1
!     integer, parameter                                :: n2 = r-q

!     double precision, dimension(n1+1)                 :: l1
!     double precision, dimension(n2+1)                 :: l2



     n1 = q-p+1
     n2 = r-q
     allocate(l1(n1+1), l2(n2+1))

     l1(1:n1) = y(p:p+n1-1)
     l2(1:n2) = y(q+1:q+n2)

     l1(n1+1) = maxval(l1(1:n1)) + 1.0d0
     l2(n2+1) = maxval(l2(1:n2)) + 1.0d0
     
     i=1
     j=1
     do k=p,r

          if (l1(i) .le. l2(j)) then
               y(k) = l1(i)
               i = min(i,n1) +1    ! need to have an in-bounds index for comparison.
          else
               y(k) = l2(j)
               j = min(j,n2) +1
          end if
     end do

     deallocate(l1,l2)

end subroutine mmerge
