program hdf_1d_rw_test
!
! Tests the subroutines 
!
! hdf_create_file.f90,
! hdf_add_1d_darray_to_file.f90,
! hdf_read_1d_darray.f90
!
! in the process of reading and writing 
! to an h5 file.
!
! The array to be read in is handled in the module almod.
!

use mod_readbuff
use hdf5
implicit none

     integer                                      :: i,j,m,n
     double precision, dimension(:), allocatable  :: arA
!     double precision, dimension(:), allocatable  :: arB
     character(len=1024)                :: fname,dsetname,description
     character(len=1024)                :: fname2,dsname2 ! temporary

     parameter(m=5)

     ! Make a fake dataset.
     allocate( arA(1:m) )

     arA = (/ 1.0d0, 4.0d0, 2.0d0, 8.0d0, 6.0d0 /)

     ! Specify filenames, create the file, and write to it.
     fname = "test_1d.h5"
     dsetname = "sample_array"
     description = "...sister? Your feelings have now betrayed her too."
     call hdf_create_file(fname)
     call hdf_add_1d_darray_to_file(m,arA,fname,dsetname,description)

     ! The way these functions are designed, the file is closed after the write.
     ! Now open it, read the size of the array inside, allocate 
     ! space for B, fill it, and compare the new array to the old one.
     
     fname2 = "c_final.h5"
     dsname2 = "Avgd_Skewness"
     call hdf_read_1d_darray(n,fname,dsetname)
!     call hdf_read_1d_darray(n,fname2,dsname2)

     
     write(*,*) "Array A: ", arA
     write(*,*) "Array B: ", arB

     deallocate(arA,arB)

end program hdf_1d_rw_test
