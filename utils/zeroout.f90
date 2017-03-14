subroutine zeroout(n,x)
! Zeroes the elements of x.
implicit none
     integer, intent(in)                          :: n
     double precision, dimension(n), intent(out)  :: x
     integer                                      :: i

     do i=1,n
          x(i) = 0.0d0
     end do
end subroutine zeroout
