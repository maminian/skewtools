program hdf_test1
! Testing writing an array by going line by line 
! using hdf_write_to_open_2d_darray.

  USE HDF5 ! This module contains all necessary modules

  IMPLICIT NONE
! moo
     integer   :: m,n,i,j
     double precision, dimension(:,:), allocatable     :: array
     character(len=10)   :: temp
! hdf
  CHARACTER(LEN=12) :: filename  ! File name
  CHARACTER(LEN=9) :: dsetname  ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


  INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
  INTEGER     ::   rank = 2                        ! Dataset rank

  INTEGER     ::   error ! Error flag

! --------------------

     filename = "hdf_test1.h5"
     dsetname = "hdf_test1"

     ! Read inputs
     call get_command_argument(1,temp)
     read(temp,*) m
     call get_command_argument(2,temp)
     read(temp,*) n

     allocate(array(m,n))
     
     ! Make an example array.
     do i=1,m
          do j=1,n
               array(i,j) = dble(i)/dble(j)
          end do
     end do
     
     dims = (/m,n/)
! -----------------------------------

  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(error)

  !
  ! Create a new file using default properties.
  !
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

  !
  ! Create the dataspace.
  !
  CALL h5screate_simple_f(rank, dims, dspace_id, error)

  !
  ! Create the dataset with default properties.
  !
  CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
       
! Modify line by line.
     do i=1,m
          call hdf_write_to_open_2d_darray(m,n,i,array(i,:),dset_id)
     end do
  !
  ! End access to the dataset and release resources used by it.
  !
  CALL h5dclose_f(dset_id, error)

  !
  ! Terminate access to the data space.
  !
  CALL h5sclose_f(dspace_id, error)

  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id, error)

  !
  ! Close FORTRAN interface.
  !
  CALL h5close_f(error)
     
     ! Deallocate.
!     deallocate(array)

end program hdf_test1
