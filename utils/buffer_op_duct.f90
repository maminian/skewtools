subroutine buffer_op_duct(bk,nTot,buffer_len,Xbuffer,Ybuffer,Zbuffer,&
                              X,Y,Z,tsteps,inext,dset_id_X,dset_id_Y,dset_id_Z)
! The basic buffered write operation.
! Saves the most recent buffer_len timesteps in an array 
! before writing to the .h5 file (otherwise file i/o dominates computation time).

use HDF5
implicit none

     integer, intent(in)                                                   :: nTot,tsteps,buffer_len
     integer, intent(inout)                                                :: bk,inext
     integer(hid_t), intent(inout)                                         :: dset_id_X,dset_id_Y,dset_id_Z

     double precision, dimension(1:nTot), intent(in)                       :: X,Y,Z
     double precision, dimension(1:buffer_len,1:nTot), intent(inout)       :: Xbuffer,Ybuffer,Zbuffer

     

     bk = bk + 1

     Xbuffer(bk,:) = X
     Ybuffer(bk,:) = Y
     Zbuffer(bk,:) = Z


     ! If we've hit the end of the buffer, write it to the appropriate
     ! location in the h5 file, and "reset" the buffer (by setting bk=0).
     if (bk .eq. buffer_len) then
          call hdf_write_to_open_2d_darray(tsteps,nTot,inext,buffer_len,Xbuffer,dset_id_X)
          call hdf_write_to_open_2d_darray(tsteps,nTot,inext,buffer_len,Ybuffer,dset_id_Y)
          call hdf_write_to_open_2d_darray(tsteps,nTot,inext,buffer_len,Zbuffer,dset_id_Z)

          bk = 0
          inext = inext + buffer_len
     end if

     
end subroutine buffer_op_duct
