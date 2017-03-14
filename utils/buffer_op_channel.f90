subroutine buffer_op_channel(bk,nTot,buffer_len,Xbuffer,Ybuffer,&
                              X,Y,tsteps,inext,dset_id_X,dset_id_Y)
! The basic buffered write operation.

use HDF5
implicit none

     integer, intent(in)                                                   :: nTot,tsteps,buffer_len
     integer, intent(inout)                                                :: bk,inext
     integer(hid_t), intent(inout)                                         :: dset_id_X,dset_id_Y

     double precision, dimension(1:nTot), intent(in)                       :: X,Y
     double precision, dimension(1:buffer_len,1:nTot), intent(inout)       :: Xbuffer,Ybuffer
     

     ! We need to reshape X and Y to put them in the buffer; it doesn't 
     ! care about nRounds.
     !
     ! If we really do, each round will be packed in blocks size ng;
     ! the first from indices 1,...,nGates, the second round nGates+1,..,2*nGates, etc.

     bk = bk + 1

     Xbuffer(bk,:) = X
     Ybuffer(bk,:) = Y
     
     ! If we've hit the end of the buffer, write it to the appropriate
     ! location in the h5 file, and "reset" the buffer (by setting bk=0).
     if (bk .eq. buffer_len) then

          call hdf_write_to_open_2d_darray(tsteps,nTot,inext,buffer_len,Xbuffer,dset_id_X)

          call hdf_write_to_open_2d_darray(tsteps,nTot,inext,buffer_len,Ybuffer,dset_id_Y)

          bk = 0
          inext = inext + buffer_len
     end if

     
end subroutine buffer_op_channel
