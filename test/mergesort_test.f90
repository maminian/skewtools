! Test of the mergesort.
! Requires an input n.
!
program mergesorttest

use mtmod

implicit none

     double precision, dimension(:), allocatable  :: x,y
     integer                                      :: i,m
     character(len=128)                           :: buffer
     logical, parameter                           :: printarrays = .false.

     integer(8), parameter                        :: mt_seed = 26827445

     logical is_ascending

     ! Set a seed for the rng.
     call sgrnd(mt_seed)

     ! Get the size of the array
     call get_command_argument(1,buffer)
     read(buffer,*) m
     write(*,*) "Array size: ",m

     allocate(x(m), y(m))

     ! Make a random list of numbers     
     call my_normal_rng(m,x,0.0d0,1.0d0)

     write(*,*) "Original array is sorted?", is_ascending(m,x)

     if (printarrays) then
          write(*,*) ""
          write(*,*) "Unsorted list:"
          write(*,*) ""
          do i=1,m
               write(*,*) x(i)
          end do
     end if

     y = x

!mergesort(n,p,r,x)
     call mergesort(m,1,m,y)

     if (printarrays) then
          write(*,*) ""
          write(*,*) "Sorted list:"
          write(*,*) ""
          do i=1,m
               write(*,*) y(i)
          end do
     end if

     
     write(*,*) "New array is sorted?", is_ascending(m,y)

     deallocate(x,y)

end program mergesorttest
!
! -------------------------------------------------------------------
!
logical function is_ascending(m,y)
implicit none
     integer, intent(in)                          :: m
     double precision, dimension(m), intent(in)   :: y

     integer i

     is_ascending = .true.
     do i=1,m-1
          if (y(i+1) - y(i) .lt. 0.0d0) then
               is_ascending = .false.
               go to 666
          end if
     end do

666  continue
     
end function is_ascending
