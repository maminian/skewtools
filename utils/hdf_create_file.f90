subroutine hdf_create_file(filename)
! Creates a blank h5 file with the given filename.

     USE HDF5 ! This module contains all necessary modules 
        
     IMPLICIT NONE

     CHARACTER(LEN=1024)      :: filename
     INTEGER(HID_T)           :: file_id     ! File identifier
 
     INTEGER                  :: error  ! Error flag
     
     !
     !    Initialize FORTRAN interface.
     !
     CALL h5open_f(error)
     
     !
     ! Create a new file using default properties.
     ! 
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

     !
     ! Terminate access to the file.
     !
     CALL h5fclose_f(file_id, error)
     
     !
     !    Close FORTRAN interface.
     !
     CALL h5close_f(error)
     
end subroutine hdf_create_file
