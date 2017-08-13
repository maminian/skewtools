subroutine median(n,x,median)
! Combination function to compute the median 
! of an array.
!
! Sort followed 
implicit none
     ! In/out
     integer, intent(in)                          :: n
     double precision, dimension(1:n), intent(in) :: x
     double precision, intent(out)                :: median
     
     ! Internal
     integer                            :: i
     double precision                   :: temp
     double precision, dimension(1:n)   :: xtemp

     ! Handle degenerate cases.
     if (n .eq. 0) then
          median = 0.0d0
     else if (n .eq. 1) then
          median = x(1)
     else

     ! Sort the array, store in xtemp
     call mergesort(n,x,xtemp)
     
     ! Take the n/2 element or the average of the two.
     i = int(n/2.0d0)
     if (mod(n,2) .eq. 0) then
          median = (xtemp(i) + xtemp(i+1))/2.0d0
     else
          median = xtemp(i)
     end if
          
end subroutine moments
