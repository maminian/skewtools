program sort_tester
! Tests the program sortpairs on three examples of varying
! degeneracy.
!
! See the file sortpairs.f90 for description
! of what it does. Relative to this file it 
! is located at ../utils/sortpairs.f90.
implicit none

     integer, dimension(7)              :: bin_idxs1
     integer, dimension(3)              :: bin_idxs2
     integer, dimension(2)              :: bin_idxs3
     integer, dimension(50)             :: bin_idxs4

     double precision, dimension(7)     :: X1,X1dup
     double precision, dimension(3)     :: X2,X2dup
     double precision, dimension(2)     :: X3,X3dup
     double precision, dimension(50)    :: X4,X4dup

     integer                            :: i,nbins,nTot
     
     character(len=14)                  :: st

     st = "(I2,A3,ES10.3)"

     ! -------------------------------------------
     ! Case 1:
     !
     ! There are 4 bins, and each is represented.
     !

     bin_idxs1 = (/ 2, 1, 2, 4, 3, 4, 4 /)
     X1 = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 /)

     nTot = 7
     nbins = 4
     
     write(*,*) "Before:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     
     do i=1,nTot
          write(*,st) bin_idxs1(i),"   ",X1(i)
     end do

     write(*,*) ""
     write(*,*) "After:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     write(*,*) ""

     call sortpairs(nTot,X1,X1dup,bin_idxs1,nbins)

     do i=1,nTot
          write(*,st) bin_idxs1(i),"   ",X1dup(i)
     end do

     write(*,*) ""
     write(*,*) "Press Enter to continue."
     write(*,*) ""

     read(*,*) 

     ! -------------------------------------------
     ! Case 2:
     !
     ! There are 12 bins, only two isolated bins are represented.
     !
     bin_idxs2 = (/ 10, 7, 10 /)
     X2 = X1(1:3)     

     nTot = 3
     nbins = 12
     
     write(*,*) "Before:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     
     do i=1,nTot
          write(*,st) bin_idxs2(i),"   ",X2(i)
     end do

     write(*,*) ""
     write(*,*) "After:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     write(*,*) ""

     call sortpairs(nTot,X2,X2dup,bin_idxs2,nbins)

     do i=1,nTot
          write(*,st) bin_idxs2(i),"   ",X2dup(i)
     end do

     write(*,*) ""
     write(*,*) "Press Enter to continue."
     write(*,*) ""

     read(*,*)

     ! -------------------------------------------
     ! Case 3:
     !
     ! There are 2 bins, and each is represented.
     !
     bin_idxs3 = (/ 2,1 /)
     X3 = X1(1:2)
     
     nTot = 2
     nbins = 2
     
     write(*,*) "Before:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     
     do i=1,nTot
          write(*,st) bin_idxs3(i),"   ",X3(i)
     end do

     write(*,*) ""
     write(*,*) "After:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     write(*,*) ""

     call sortpairs(nTot,X3,X3dup,bin_idxs3,nbins)

     do i=1,nTot
          write(*,st) bin_idxs3(i),"   ",X3dup(i)
     end do

     write(*,*) ""
     write(*,*) "End of test."
     write(*,*) ""

     ! -------------------------------------------
     ! Case 4:
     !
     ! Five bins, many entries.
     !

     nTot = 50
     nbins = 5

     do i=1,nTot
          X4(i) = dsin(dble(i))
          bin_idxs4(i) = floor( 5*abs(dcos(dble(i))) ) + 1
     end do

     
     write(*,*) "Before:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     
     do i=1,nTot
          write(*,st) bin_idxs4(i),"   ",X4(i)
     end do

     write(*,*) ""
     write(*,*) "After:"
     write(*,*) ""
     write(*,"(A5,A5)") "idx","X"
     write(*,*) ""

     call sortpairs(nTot,X4,X4dup,bin_idxs4,nbins)

     do i=1,nTot
          write(*,st) bin_idxs4(i),"   ",X4dup(i)
     end do

     write(*,*) ""
     write(*,*) "Press Enter to continue."
     write(*,*) ""

end program sort_tester
