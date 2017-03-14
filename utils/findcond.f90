subroutine findcond(n,b,q,aptr)
!
! Given a boolean array b, dimension(n), 
! outputs an array aptr where the 
! first q elements are integer pointers 
! to the elements of b which are .true.
!

implicit none
     integer, intent(in)                :: n
     logical, dimension(n), intent(in)  :: b

     integer, intent(out)               :: q
     integer, dimension(n), intent(out) :: aptr

     integer                            :: i
     
     q = 0
     
     do i=1,n
          if (b(i)) then
               q = q+1
               aptr(q) = i
          end if
     end do

end subroutine findcond
