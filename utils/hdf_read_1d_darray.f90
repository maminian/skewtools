subroutine hdf_read_1d_darray(m,filename,dsetname)
! Read a 1d double array from an h5 file under the corresponding 
! dsetname into a temporary buffer from the module mod_readbuff.
use mod_readbuff
use hdf5
implicit none

     
     integer(HSIZE_T), intent(out)                 :: m
     character(len=1024), intent(in)               :: filename,dsetname
     
!     double precision, dimension(m), intent(out)  :: A

     ! HDF variables.
     integer(hid_t)                               :: file_id,dset_id,dspace_id
     integer                                      :: hdferror,rank

     integer(HSIZE_T), dimension(1)               :: dims,maxdims

     
     ! ------------------------------
     ! Initialize the hdf interface.
     call h5open_f(hdferror)
     
     ! Open the file.
     call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferror)
     call h5dopen_f(file_id, dsetname, dset_id, hdferror)

     ! Read the file, figuring out the dimensions, allocating, 
     ! then copying over the array.

     call h5dget_space_f(dset_id, dspace_id, hdferror)                  ! Getting the dataspace ID
     call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
     call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferror)  ! Getting dims from dataspace


     m = dims(1)
     allocate(readbuff_double(m))

     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, readbuff_double, dims, hdferror, h5S_ALL_F, dspace_id)  ! Reading array of size dims.
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, hdferror, h5S_ALL_F, dspace_id)  ! Reading array of size dims.

     ! Close the dataset, file, and hdf interface.
     call h5sclose_f(dspace_id, hdferror)
     call h5dclose_f(dset_id, hdferror)
     call h5fclose_f(file_id, hdferror)
     call h5close_f(hdferror)
     
!     A = readbuff_double

!     deallocate(readbuff_double)
     
end subroutine hdf_read_1d_darray
