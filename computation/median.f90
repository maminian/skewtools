subroutine median(n,x,med)
! Combination function to compute the median 
! of an array.
!
! Done by using a generic sort on input array x
! and returning the n/2 element (or average for n odd)
!
implicit none
     ! In/out
     integer, intent(in)                          :: n
     double precision, dimension(1:n), intent(in)   :: x
     double precision, intent(out)                :: med
     
     ! Internal
     integer                            :: i
     double precision, dimension(1:n)     :: xtemp

     ! Handle degenerate cases.
     if (n .eq. 0) then
          med = 0.0d0
     else if (n .eq. 1) then
          med = x(1)
     end if

     ! Sort the array, store in xtemp
     xtemp = x
     call mergesort(n,1,n,xtemp)
     
     ! Take the n/2 element or the average of the two.
     i = int(n/2.0d0)
     

     if (mod(n,2) .eq. 0) then
          med = (xtemp(i) + xtemp(i+1))/2.0d0
     else
          med = xtemp(i+1)
     end if
          
end subroutine median
