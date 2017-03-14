subroutine hdf_add_2d_darray_to_file(m,n,A,filename,arrayname,description)
! Given an hdf file already created, takes an array with 
! dimensions m,n and writes it to the hdf file.
!
! Baby steps. Test it with the dump_uduct_flow code.


use hdf5
implicit none

     ! Inputs
     integer                            :: m,n
     double precision, dimension(m,n)   :: A
     character(len=1024)                :: filename,arrayname,attrname
     character(len=1024)                :: description
     
     ! HDF variables.
     integer(hid_t)                     :: file_id,dset_id,dspace_id, &
                                             attr_id,aspace_id,atype_id
     integer(hsize_t), dimension(1)     :: adims
     integer                            :: rank,error
     integer                            :: arank
     integer(HSIZE_T), dimension(2)     :: data_dims
     integer(size_t)                    :: attrlen
     
     ! Misc.
     
     
     parameter(rank=2)   ! Dimension of array.
     parameter(arank=1)  ! Rank of attribute (size of attribute array?)
     parameter(adims=(/1/))  ! Size of array of hdf attributes. For our purposes,
                         ! only using 1.
     parameter(attrname="Description")
     
     ! ------------------------------
     
     ! Do some preliminary work
     data_dims(1) = m
     data_dims(2) = n
     
     
     ! Initialize interface, open the file.
     call h5open_f(error)
     call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
     
     ! Create the dataset and dataspace and all that.
     call h5screate_simple_f(rank, data_dims, dspace_id, error)
     call h5dcreate_f(file_id, arrayname, H5T_NATIVE_DOUBLE, dspace_id, &
                         dset_id, error)
     
     ! Write array.
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, data_dims, error)

     ! Write the text description for the array.
     ! Turns out this requires making the datatype and whatnot.
     description = trim(description)
     attrlen = len_trim(description)
     
     call h5screate_simple_f(arank, adims, aspace_id, error)
     call h5screate_simple_f(arank, adims, aspace_id, error)
     call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     call h5tset_size_f(atype_id, attrlen, error)
     call h5acreate_f(dset_id, attrname, atype_id, aspace_id, attr_id, error)
     ! The write happens here.
     call h5awrite_f(attr_id, atype_id, description, adims, error)
     call h5aclose_f(attr_id, error)
     
     ! Close the dataset, file, and hdf interface.
     call h5dclose_f(dset_id, error)
     call h5fclose_f(file_id, error)
     call h5close_f(error)
     
     ! EXIT
     
end subroutine hdf_add_2d_darray_to_file
