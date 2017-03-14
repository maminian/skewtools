subroutine generate_target_times(tmin,tstep_type,Tfinal,other_file)
! From the input options, fill in an array 
! target_times which will be the times on which 
! output data is saved.

use mod_readbuff ! For passing an unallocated array between here and the subroutine
             ! hdf_read_1d_darray().

use mod_time     ! For time-related variables and arrays.

implicit none

!     integer, intent(inout)                            :: ntt
!     double precision, dimension(1:ntt), intent(inout) :: target_times
     double precision, intent(inout)                   :: Tfinal
     double precision, intent(in)                      :: tmin
     character(len=1024), intent(in)                   :: tstep_type,other_file

     
     character(len=1024)                               :: expo,unif,supplied,stt
     integer                                           :: tt_idx
     double precision                                  :: kscale,dt

     logical                                           :: flag
          
     parameter(expo='expo',unif='unif',supplied="supplied",stt="target_times")
     

     if (tstep_type .eq. unif) then
          ! Uniform timestepping using the specified tmin and tfinal with 
          ! ntt timesteps.
          allocate(target_times(ntt))
          target_times(1) = 0.0d0
     
          dt = (Tfinal-tmin)/dble(ntt-1)
          ! Generate the target times
          do tt_idx=2,ntt
               target_times(tt_idx) = tmin + dt*(tt_idx-1)
          end do

     else if (tstep_type .eq. expo) then
          ! Exponential timestepping; these are uniformly spaced in a log scale.
          ! This is the usual setting for our purposes.

          allocate(target_times(ntt))
          target_times(1) = 0.0d0
     
          target_times(2) = tmin
          kscale = dble(Tfinal/tmin)**(1.0d0/dble(ntt-2))
          
          do tt_idx=3,ntt-1
               target_times(tt_idx) = kscale*target_times(tt_idx-1)
          end do
          
          target_times(ntt) = Tfinal

     else if (tstep_type .eq. supplied) then
          ! User has supplied their own timesteps in the specified h5 file.
          ! Override all other settings and use them.


          call hdf_read_1d_darray(ntt,other_file,stt)

          flag = (.not. (readbuff_double(1) .eq. 0.0d0))
          if (flag) then
               ntt = ntt + 1
          end if

          allocate(target_times(ntt))

          if (flag) then
               target_times(1) = 0.0d0
               target_times(2:ntt) = readbuff_double
          else
               target_times = readbuff_double
          end if

          Tfinal = target_times(ntt)

          deallocate(readbuff_double)

     else 
          write(*,*) "Unrecognized tstep_type. Use either ''unif'', ''expo'', or ''supplied''."
     end if
     
     Tfinal = target_times(ntt)

end subroutine generate_target_times
