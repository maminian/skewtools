! Test of the median subroutine.
!
program median_test

use mtmod

implicit none

     integer, parameter                           :: m = 8
     integer, parameter                           :: p = 7
     double precision, dimension(m)               :: x
     double precision, dimension(p)               :: y

     integer                                      :: i,j

     integer(8), parameter                        :: mt_seed = 26827445
     double precision                             :: med

     ! Set a seed for the rng.
     call sgrnd(mt_seed)

     write(*,*) "--------------------------"
     write(*,*) ""
     write(*,*) "Test 1: array all the same value."
     write(*,*) ""
     write(*,*) "--------------------------"
     x = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)
     y = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)

     write(*,*) "Array size:",m
     write(*,*) "Elements:"
     write(*,*) ""
     call pev(m,x)
     write(*,*) ""

     call median(m,x,med)
     write(*,*) "Median:",med

     write(*,*) "Array size:",p
     write(*,*) "Elements:"
     call pev(p,y)
     write(*,*) ""

     call median(p,y,med)
     write(*,*) "Median:",med


     write(*,*) "--------------------------"
     write(*,*) ""
     write(*,*) "Test 2: ascending elements centered around zero."
     write(*,*) ""
     write(*,*) "--------------------------"

     x = (/ -4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0 /)
     y = (/ -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0 /)

     write(*,*) "Array size:",m
     write(*,*) "Elements:"
     write(*,*) ""
     call pev(m,x)
     write(*,*) ""

     call median(m,x,med)
     write(*,*) "Median:",med

     write(*,*) "Array size:",p
     write(*,*) "Elements:"
     call pev(p,y)
     write(*,*) ""

     call median(p,y,med)
     write(*,*) "Median:",med


     write(*,*) "--------------------------"
     write(*,*) ""
     write(*,*) "Test 3: assorted elements."
     write(*,*) ""
     write(*,*) "--------------------------"

     call my_normal_rng(m,x,0.0d0,1.0d0)
     call my_normal_rng(p,y,0.0d0,1.0d0)

     write(*,*) "Array size:",m
     write(*,*) "Elements:"
     write(*,*) ""
     call pev(m,x)
     write(*,*) ""

     call median(m,x,med)
     write(*,*) "Median:",med

     write(*,*) "Array size:",p
     write(*,*) "Elements:"
     call pev(p,y)
     write(*,*) ""

     call median(p,y,med)
     write(*,*) "Median:",med


end program median_test

subroutine pev(m,x)
! Print elements vertically
implicit none
     integer, intent(in)                          :: m
     double precision, dimension(m), intent(in)   :: x

     integer i

     do i=1,m
          write(*,*) x(i)
     end do
end subroutine pev
