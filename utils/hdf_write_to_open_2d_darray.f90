subroutine hdf_write_to_open_2d_darray(m,n,i,nrow,array,dset_id)
! As the name suggests; given an opened dataset id,
! modifies the i-th through (i+nrow-1) row of it.
!
! The h5 dataset is assumed to have total dimension (m,n).
!
! Does no error checking whatsoever. Don't be dumb!

use hdf5
implicit none
     ! Inputs
     integer                                      :: m,n,i,nrow
     double precision, dimension(1:nrow,1:n)      :: array
     INTEGER(HID_T)                               :: dset_id          ! Dataset identifier 
     
     
     ! HDF things
     INTEGER(HID_T)                               :: dataspace        ! Dataspace identifier 
     INTEGER(HID_T)                               :: memspace         ! memspace identifier
     integer                                      :: error,rank
     
     integer(hsize_t), dimension(2)               :: offset,stride,block,steps,dimsm
     
     rank=2
     offset = (/i-1,0/)                 ! Which element to start at. HDF counts from zero.
     stride = (/1,1/)                   ! Write sequential elements
     block = (/1,1/)                    ! No blocks.
     steps = (/nrow,n/)                 ! How many times to 'repeat the pattern'
                                        ! in each direction.
                                        ! In this case, the size of the array.

     dimsm=(/nrow,n/)                   ! Dimensions of subset to write to dataset.


     !
     ! Get dataset's dataspace identifier and select subset.
     !
     if (i+nrow-1 .gt. m) then
          write(*,*) "Warning: possibly writing past the end of a HDF dataset."
     end if

     CALL h5dget_space_f(dset_id, dataspace, error)

     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       offset, steps, error, stride, block) 

     !
     ! Create memory dataspace.
     !
     CALL h5screate_simple_f(rank, dimsm, memspace, error)

     !
     ! Write subset to dataset  
     !
     
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, error, &
       memspace, dataspace)


end subroutine hdf_write_to_open_2d_darray
