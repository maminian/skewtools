program channel_mc
! Program to do Monte Carlo in a channel.

use HDF5
use mtmod
use mod_time
use mod_duration_estimator
use mod_readbuff
use mod_parameters
use mod_active_particles
use mod_oscillatory_channel

implicit none


     ! -------------------------------
     ! Things moved to time_module.f90
     !
     !     integer                                           :: ntt
     !     double precision, dimension(:), allocatable       :: target_times
     !
     ! ----------------------------------------------

!     integer, parameter                           :: i64 = selected_int_kind(18)
     integer(kind=i64)                            :: mc_n
     integer                                      :: nt,kt,ny,i,nTot
     double precision                             :: t

     
     ! Internal variables for checking if we're at a target time during the timestepping.
     double precision                             :: next_tt
     integer                                      :: tt_idx
     
     character(len=1024)                               :: param_file, filename, out_msg
     
     double precision, dimension(:), allocatable       :: X,Y,Xtemp,Ytemp    ! Position arrays
     double precision, dimension(:,:), allocatable     :: Xbuffer,Ybuffer

     integer(kind=i64)                                 :: bk,inext,rem
     double precision, dimension(:), allocatable       :: means,vars,skews,kurts,t_hist
     double precision, dimension(:,:), allocatable     :: hist_centers,hist_heights

     
     ! Moments through slices and distributions
     double precision, dimension(:,:), allocatable     :: means_sl,vars_sl,skews_sl,kurts_sl
     integer                                           :: bin_count,nbx,nby
     double precision                                  :: dby     ! Bin width

     ! Stuff for 2d histogram looking into the short direction.
     double precision, dimension(:,:,:), allocatable   :: hist2d
     double precision, dimension(:,:), allocatable     :: hist2dcx, hist2dcy ! bin centers.

     ! HDF
     integer                                           :: rank,h5error
     character(len=1024)                               :: fname2,descr,arrayname
     
     character(len=1024)                :: dsetname
     
     ! HDF variables.
     integer(hid_t)                     :: file_id
     integer(hid_t)                     :: dset_id_X,dset_id_Y,dset_id_Z,&
                                             dset_id_Sk,dset_id_t
     integer(hid_t)                     :: dspace_id_X,dspace_id_Y,dspace_id_Z,&
                                             dspace_id_Sk,dspace_id_t
     
     integer(hsize_t), dimension(2)     :: data_dims


     ! Flags to save position histories and read IC from a file.
     logical                            :: check_ic_channel,robin

     double precision u_channel_precomp
     

     ! Advection/diffusion functions!
     external  :: impose_reflective_BC_rect_robin, impose_reflective_BC_rect, u_channel, u_channel_precomp

     ! TEMPORARY MANUAL PARAMETER
     !
     robin = .false.

     ! Read all parameters from file.
     call get_command_argument(1,param_file)
     call get_command_argument(2,filename)
     
     call read_inputs_mc2(param_file) ! Parameters are passed using mod_parameters.          
     if (use_oscillatory) then
          call read_oscillatory_params(param_file)
     end if

     if (filename=="") then
          out_msg = 'missing_args'
          call channel_mc_messages(out_msg)
          go to 1234
     end if
     
	! Assign the number of bins in each direction.
     if (nbins .eq. 0) then
          nby = 0
          dby = 0.0d0
     else
          nby = nbins
          dby = 2.0d0*a/nby 
     end if
     
     ! ----------
     ! Other hard coded stuff
     nbx = nby

     ! --------------------------------
     !
     ! Based on the input, generate the array target_times, (times to save output)
     ! and get information for the internal array t_hist.
     !
     
     call generate_target_times(dt,tstep_type,Tfinal,other_file)
     
     ! Get the value of nt before allocating arrays.
     call correct_tstep_info(ntt,nt,target_times,dtmax)

     ! Initialize the Mersenne Twister RNG with seed read in from the 
     ! input files.
     
     call sgrnd(mt_seed)

     ! Get the total number of simulations
     ny = nGates
     nTot = nGates*x0n
     
     ! Correct nTot if we are using an external file.
     if (use_external_ic) then
          dsetname = "nTot"
          call hdf_read_1d_darray(i,ic_file,dsetname)
          nTot = int(readbuff_double(1))
          deallocate(readbuff_double)
     else
     ! Recalculate nGates and nTot based on the values for ny and nz.
          nGates = ny
          nTot = nGates*x0n
     end if
     
     ! ------------------
     ! Allocate arrays.
     !
     
     ! Time
     allocate(t_hist(nt))
     
     ! Positions
     allocate(X(nTot),Y(nTot),Xtemp(nTot),Ytemp(nTot))
     if (robin) then
          ! Allocate the array for tracking the active set.
          allocate(aset(nTot),asettemp(nTot))
          nactive = nTot
          do i=1,nTot
               asettemp(i) = i
          end do
          aset = asettemp
     end if

     ! Channel-averaged stats
     allocate(means(ntt),vars(ntt),skews(ntt),kurts(ntt))

     ! Stats on slices
     if (.not. (nbins .eq. 0)) then
          allocate(means_sl(ntt,nbins),vars_sl(ntt,nbins),skews_sl(ntt,nbins),kurts_sl(ntt,nbins))
     end if

     ! If the 2d histogram (looking in the short direction) is desired, 
     ! allocate.
     if (save_hist2d) then
          allocate(hist2d(ntt,nbx,nby))
          allocate(hist2dcx(ntt,nbx))
          allocate(hist2dcy(ntt,nby))
     end if
     
     ! Cross-sectionally averaged distribution
     allocate(hist_centers(ntt,nhb),hist_heights(ntt,nhb))
     


     ! --------------------------
     ! Generate initial conditions and internal timestepping.
     !
     if (use_oscillatory) then
          call precalculate_efuncs()
     end if

     call set_initial_conds_channel_mc(ny,nGates,x0n,nTot,X,Y,y0,a,x0width,t_warmup,use_external_ic,ic_file)
     
     if (.not. check_ic_channel(nTot,Y,a)) then
          write(*,*) "Part of the initial condition lies outside the domain. Exiting."
          go to 1234
     end if
     
     call generate_internal_timestepping(ntt,nt,target_times,t_hist,dtmax)
     
     call print_parameters(aratio,q,Pe,nGates,nTot,y0,z0,save_hist,&
                              t_hist,dtmax,nt,ntt,mt_seed,geometry,use_external_ic)

     ! --------------------------------------------------------
     ! Initialize HDF with appropriate dataset, etc.

     fname2 = trim(filename)

     call hdf_create_file(fname2)
     
     ! We need to open the h5 file after hdf_create_file
     ! because the interface is "global" amongst all files 
     ! containing the hdf5 module.
     
     call h5open_f(h5error)
     
     
     if (save_hist) then
          ! Set up dataspaces in the hdf file for:
          ! X, Y.
          ! 
          ! Allocate memory for the memory buffers here, too.
          
          allocate(Xbuffer(buffer_len,nTot))
          allocate(Ybuffer(buffer_len,nTot))
          
          call h5fopen_f(fname2, H5F_ACC_RDWR_F, file_id, h5error)

          data_dims(1) = ntt
          data_dims(2) = nTot
          rank = 2
          
          dsetname = "X"
          call h5screate_simple_f(rank, data_dims, dspace_id_X, h5error)
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id_X, &
                         dset_id_X, h5error)

          dsetname = "Y"
          call h5screate_simple_f(rank, data_dims, dspace_id_Y, h5error)
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id_Y, &
                         dset_id_Y, h5error)
          
     end if

     ! Save a history of time per iteration for predicting time to completion.
     mde_ntt = nt-1
     allocate(mde_dts(mde_ntt))
     

     inext = 1

     tt_idx = 1
!     call accumulate_moments_1d(tt_idx,ntt,nTot,X,Y,a,means,vars,skews,&
!                    kurts,nby,means_sl,vars_sl,skews_sl,kurts_sl)
     call accumulate_moments_1d(tt_idx,ntt,nTot,X,Y,a,means,vars,skews,&
                    kurts,nby,means_sl,vars_sl,skews_sl,kurts_sl)
     call make_histogram(nTot,X,nhb,hist_centers(tt_idx,1:nhb),hist_heights(tt_idx,1:nhb))

     if (save_hist2d) then
          call make_histogram2d(nTot,X,Y,nbx,nby,hist2dcx(tt_idx,1:nbx),&
                    hist2dcy(tt_idx,1:nby), hist2d(tt_idx,1:nbx,1:nby))
     end if


	tt_idx = 2
	next_tt = target_times(tt_idx)
     
     ! -----------
     ! Prepare the buffer to save position histories if requested.

     bk = 0
     if (save_hist) then
          call buffer_op_channel(bk,nTot,buffer_len,Xbuffer,Ybuffer,&
                                   X,Y,ntt,inext,dset_id_X,dset_id_Y)
     end if


     ! --------------
     ! Start the timestepping.
     
     out_msg = 'simul_start'
     call channel_mc_messages(out_msg)

     do kt=2,nt

          call system_clock(mde_t1,count_rate)    ! Time for progress.

          ! Push forward time.
          t = t_hist(kt)
          dt = t_hist(kt) - t_hist(kt-1)

          if (use_oscillatory) then
               call update_flow(t)
               call apply_advdiff1_chan(nTot,X,Y,Pe,dt,a, &
                              u_channel_precomp,impose_reflective_BC_rect)
          else
               !
               ! Primary timestep
               !
               if (.not. robin) then
                    call apply_advdiff1_chan(nTot,X,Y,Pe,dt,a, &
                                   u_channel,impose_reflective_BC_rect)
               else
                    if (nactive .gt. 0) then
                         call apply_advdiff_chan_robin(nTot,X,Y,Pe,dt,a, &
                                        u_channel,impose_reflective_BC_rect_robin)
                    end if
               end if
          end if

          
          !
          ! If we're at a target time,  
          ! calculate and save moments (and positions, if requested),
          ! then increment tt_idx and update next_tt.
          !
          if ( t .eq. next_tt ) then
               ! Output centerline flow velocity.
!               write(1,*) next_tt,u_channel_precomp(0.0d0,-1.0d0,1.0d0)

               do i=1,nactive

                    Xtemp(i) = X(aset(i))
                    Ytemp(i) = Y(aset(i))
               end do
               
               if (.not. robin) then 
                    ! Update the statistics.
                    call accumulate_moments_1d(tt_idx,ntt,nTot,X,Y,a,means,vars,skews,&
                                   kurts,nby,means_sl,vars_sl,skews_sl,kurts_sl)
                    ! Update the histogram centers and heights.
                    call make_histogram(nTot,X,nhb,hist_centers(tt_idx,1:nhb),hist_heights(tt_idx,1:nhb))

                    if (save_hist2d) then
                         call make_histogram2d(nTot,X,Y,nbx,nby,hist2dcx(tt_idx,1:nbx),&
                                   hist2dcy(tt_idx,1:nby), hist2d(tt_idx,1:nbx,1:nby))
                    end if
               else
                    if (nactive .gt. 0) then
                         call accumulate_moments_1d(tt_idx,ntt,nactive,&
                                   Xtemp(1:nactive),Ytemp(1:nactive),a,means,vars,skews,&
                                   kurts,nby,means_sl,vars_sl,skews_sl,kurts_sl)
                         call make_histogram(nactive,Xtemp(1:nactive),nhb,&
                                        hist_centers(tt_idx,1:nhb),hist_heights(tt_idx,1:nhb))

                         if (save_hist2d) then
                              call make_histogram2d(nactive,X(1:nactive),Y(1:nactive),&
                                        nbx,nby,hist2dcx(tt_idx,1:nbx),&
                                        hist2dcy(tt_idx,1:nby), hist2d(tt_idx,1:nbx,1:nby))
                         end if
                    end if
               end if

               !
               ! Write history if requested.
               !
               
               if (save_hist) then
                    call buffer_op_channel(bk,nTot,buffer_len,&
                                        Xbuffer,Ybuffer,X,Y,ntt,inext,dset_id_X,dset_id_Y)
               end if

               !
               ! Update the target time and array index.
               !
               if (next_tt .lt. Tfinal) then
                    tt_idx = tt_idx + 1
                    next_tt = target_times(tt_idx)
               end if
          end if
          
          call system_clock(mde_t2,count_rate)    ! Time in milliseconds

          mde_ntc = kt-1          
          mde_dts(mde_ntc) = (mde_t2-mde_t1)/dble(count_rate)  ! Time in seconds

          ! Display percentage progress. The last argument as .true. should be used with gfortran
          ! (or any other compiler that supports "\b"), or .false. with ifort.
          call progress_meter(kt,nt,.true.)
               
     end do
     
     out_msg = 'simul_done'
     call channel_mc_messages(out_msg)
     
     if (save_hist) then
     ! Write the remainder of the buffer, then close the file.
     
          rem = ntt-inext+1
          if (rem .gt. 0) then
               call hdf_write_to_open_2d_darray(ntt,nTot,inext,rem,Xbuffer(1:rem,1:nTot),dset_id_X)
               call hdf_write_to_open_2d_darray(ntt,nTot,inext,rem,Ybuffer(1:rem,1:nTot),dset_id_Y)
          end if

          call h5dclose_f(dset_id_X,h5error)
          call h5dclose_f(dset_id_Y,h5error)
          call h5fclose_f(file_id,h5error)
     end if
     
     ! Because of the nature of hdf5 mod for fortran,
     ! we close the interface here, since it gets
     ! re-opened in the calls below.
     call h5close_f(h5error)


     if (save_hist2d) then
          arrayname = "hist2dcx"
          descr = "Array tracking bin centers in the x direction for hist2d"
          call hdf_add_2d_darray_to_file(ntt,nbx,hist2dcx,fname2,arrayname,descr)

          arrayname = "hist2dcy"
          descr = "Array tracking bin centers in the y direction for hist2d"
          call hdf_add_2d_darray_to_file(ntt,nby,hist2dcy,fname2,arrayname,descr)

          arrayname = "hist2d"
          descr = "Array tracking the density for hist2d"
          call hdf_add_3d_darray_to_file(ntt,nbx,nby,hist2d,fname2,arrayname,descr)
     end if

     ! Save all the remaining arrays. It's a lot of fluff so it's been 
     ! given its own subroutine.

     call save_the_rest_channel(fname2,geometry,ntt,target_times,means,vars,skews,kurts,&
                    nby,means_sl,vars_sl,skews_sl,kurts_sl,nhb,hist_centers,hist_heights,&
                    Pe,nTot,mt_seed,dtmax,t_warmup)

     ! -------
     
     deallocate(X,Y)

     deallocate(means,vars,skews,kurts,target_times)

     if (.not. (nbins .eq. 0)) then
          deallocate(means_sl,vars_sl,skews_sl,kurts_sl)
     end if

     deallocate(hist_centers,hist_heights)
     deallocate(t_hist)
     
     if (save_hist) then
          deallocate(Xbuffer,Ybuffer)
     end if

     out_msg = 'done'
     call channel_mc_messages(out_msg)
     
     write(*,*) ""
     
1234 continue

end program channel_mc
