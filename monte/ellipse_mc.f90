program ellipse_mc
! Program to do Monte Carlo in an ellipse.

use HDF5
use mtmod
use mod_time
use mod_duration_estimator
use mod_parameters

implicit none

     ! Array sizes, parameters, local vars
     
     integer                                           :: nTot,nt,kt,ny,nz,tt_idx
     double precision                                  :: next_tt
     double precision                                  :: t
     
                                                            
     integer                                           :: maxrefl     

     ! Positions, position/statistic histories
     integer                                           :: nbx,nby,nbz
     double precision                                  :: dby,dbz

	 ! Stuff for 2d histogram looking into the short direction.
     double precision, dimension(:,:,:), allocatable   :: hist2d
     double precision, dimension(:,:), allocatable     :: hist2dcx, hist2dcy ! bin centers.
	 
	 
     double precision, dimension(:), allocatable       :: X,Y,Z
     double precision, dimension(:,:), allocatable     :: Xbuffer,Ybuffer,Zbuffer
     integer                                           :: bk,inext,rem
     double precision, dimension(:), allocatable       :: means,vars,skews,kurts,medians,t_hist
     
     double precision, dimension(:,:,:), allocatable   :: means_sl,vars_sl,skews_sl,kurts_sl
     double precision, dimension(:,:), allocatable     :: hist_centers,hist_heights
     
     ! i/o
     character(len=1024)                               :: param_file,filename
     character(len=1024)                               :: out_msg

	 character(len=1024)                               :: arrayname,descr
		  
     ! HDF
     integer                                           :: rank,h5error
     character(len=1024)                               :: fname2
     
     character(len=1024)                               :: dsetname
     
     ! HDF variables.
     integer(hid_t)                                    :: file_id
     integer(hid_t)                                    :: dset_id_X,dset_id_Y,dset_id_Z
     integer(hid_t)                                    :: dspace_id_X,dspace_id_Y,dspace_id_Z
     
     integer(hsize_t), dimension(2)                    :: data_dims

     ! References to functions that go in arguments.
     external  :: impose_reflective_BC_ellipse, u_ellipse

          
     ! Parameters.
     !

     ! Maximum reflections applied before giving up and stopping the point on the boundary.
     parameter(maxrefl = 10)

     geometry = "ellipse"

	 ! -------------------------------------------------------

     ! Read all parameters from file.
     call get_command_argument(1,param_file)
     call get_command_argument(2,filename)
     
     call read_inputs_mc(param_file,aratio,q,Pe,nGates,x0n,x0width,y0,z0,save_hist,nbins,save_hist2d,t_warmup,&
                              use_external_ic,ic_file,tstep_type,dt,dtmax,Tfinal,ntt,other_file,mt_seed)
     
     if (filename=="") then
          out_msg = 'missing_args'
          call duct_mc_messages(out_msg,nz)
          go to 1234
     end if

     ! Set the dimensions of the thing.
     b = a/aratio
     
	! Assign the number of bins in each direction for ptwise stats.
     if (nbins .eq. 0) then
          nby = 0
          nbz = 0

          dby = 0.0d0
          dbz = 0.0d0
     else
          nby = nbins
          nbz = ceiling(nbins/aratio) ! Could also make this the same as nby.
          
          dby = (2.0d0*a)/nby
          dbz = (2.0d0*b)/nbz
     end if
	 
	 nbx = nby
     write(*,*) nbins,nbx,nby,nbz
	 
	 
     !
     ! Generate the target times; times at which output is saved.
     ! Internal timestepping is created after.
     ! 
     
     call generate_target_times(dt,tstep_type,Tfinal,other_file)

     ! Get the value of nt before allocating arrays.
     call correct_tstep_info(ntt,nt,target_times,dtmax)
     
     ! Initialize the Mersenne Twister RNG with seed read in from the 
     ! input files.
     
     call sgrnd(mt_seed)
     

     if (nGates .gt. 1) then
          ! Need to adjust this to make the spacing
          ! uniform taking into account the aspect ratio.
          
          ny = floor(dsqrt(dble(nGates)*aratio))+1
          nz = floor(dsqrt(dble(nGates)/aratio))+1
!          ny = floor(dsqrt(dble(nGates)))+1
!          nz = floor(dsqrt(dble(nGates)))+1

          !
          ! Because of the rejection method used to generate uniform points,
          ! each of ny,nz needs to be scaled up appropriately so that 
          ! the number of points that are actually simulated is genuinely 
          ! ~nGates, as input from the user's file.
          ! 
          ! The factor is the ratio of the circle to its circumscribing square.
          ny = floor(ny/dsqrt(datan(1.0d0)))+1
          nz = floor(nz/dsqrt(datan(1.0d0)))+1
     else
          ny = 1
          nz = 1
     end if
     
     ! If resolution is an issue, exit.
     if ( (ny .lt. 7) .and. (nGates .gt. 1) ) then
          out_msg = 'resolution'
          call duct_mc_messages(out_msg,nz)
          go to 1234
     end if


     ! Calculate nGates and nTot (total number of particles).
     ! In the ellipse this is done with two "sweeps." 
     ! The number of points should be approximately ~(pi/4)*ny*nz. 
     ! Sweep the grid and update the values for nGates and nTot,
     ! then allocate memory for X,Y,Z, then fill in the values.

     call get_pts_in_ellipse(ny,nz,x0n,a,b,nGates,nTot)

     ! Internal time
     allocate(t_hist(nt))

     ! Positions
     allocate(X(nTot), Y(nTot), Z(nTot))

     ! Channel-averaged stats
     allocate(means(ntt), vars(ntt), skews(ntt),kurts(ntt),medians(ntt))
	 
     ! Stats on Y slices (integrated across Z)
     if (.not. (nbins .eq. 0)) then
          allocate(means_sl(ntt,nby,nbz),vars_sl(ntt,nby,nbz),&
                    skews_sl(ntt,nby,nbz),kurts_sl(ntt,nby,nbz))
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
          ! X, Y, Z.
          ! 
          ! Allocate memory for the memory buffers here, too.
          
          allocate(Xbuffer(buffer_len,nTot))
          allocate(Ybuffer(buffer_len,nTot))
          allocate(Zbuffer(buffer_len,nTot))
          
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
          dsetname = "Z"
          call h5screate_simple_f(rank, data_dims, dspace_id_Z, h5error)
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id_Z, &
                         dset_id_Z, h5error)
          
     end if
     
     !
     ! -----------------------------------------------
     !
     ! Set initial conditions and internal timestepping.
     !
     
     call set_initial_conds_ellipse_mc(ny,nz,x0n,a,b,nGates,nTot,X,Y,Z, &
                              y0,z0,x0width,t_warmup)

     call generate_internal_timestepping(ntt,nt,target_times,t_hist,dtmax)

     call print_parameters(aratio,q,Pe,nGates,nTot,y0,z0,save_hist,&
                              t_hist,dtmax,nt,ntt,mt_seed,geometry,use_external_ic)

     ! Save a history of time per iteration for predicting time to completion.
     mde_ntt = nt-1
     allocate(mde_dts(mde_ntt))

     ! --------------------------------
     ! 
     ! Calculate the statistics of the initial condition.
     !
     inext = 1
     tt_idx = 1
     
     call accumulate_moments_2d(tt_idx,ntt,nTot,X,Y,Z,-a,a,-b,b,means,vars,skews,&
               kurts,nby,nbz,means_sl,vars_sl,skews_sl,kurts_sl)
     call median(nTot,X,medians(tt_idx))

     call make_histogram(nTot,X,nhb,hist_centers(tt_idx,1:nhb),hist_heights(tt_idx,1:nhb))

	 if (save_hist2d) then
          call make_histogram2d(nTot,X,Y,nbx,nby,hist2dcx(tt_idx,1:nbx),&
                    hist2dcy(tt_idx,1:nby), hist2d(tt_idx,1:nbx,1:nby))
     end if
	 
     tt_idx = 2
     next_tt = target_times(tt_idx)

     ! -----------
     ! Prepare the buffer to save position histories if requested.
     ! The subroutine is geometry independent once the initial conditions are set, 
     ! so there's no need to make a "buffer_op_ellipse" subroutine.
     
     bk = 0
     
     if (save_hist) then
          call buffer_op_duct(bk,nTot,buffer_len,Xbuffer,Ybuffer,Zbuffer,&
                              X,Y,Z,ntt,inext,dset_id_X,dset_id_Y,dset_id_Z)
     end if

     ! --------------
     ! Start the timestepping.
     
     out_msg = 'simul_start'
     call duct_mc_messages(out_msg,nz)
     
     do kt=2,nt

          call system_clock(mde_t1,count_rate)    ! Time for progress.     

          ! Push forward time.
          t = t_hist(kt)
          dt = t_hist(kt) - t_hist(kt-1)

          call apply_advdiff1_ellipse(nTot,X,Y,Z,Pe,dt,a,b, &
                              u_ellipse,impose_reflective_BC_ellipse,maxrefl)

          
          !
          ! Check if we're at a target time. If we are,  
          ! and calculate and save moments (and positions, if requested),
          ! then increment tt_idx and update next_tt.
          !

          if (t .eq. next_tt) then

               call accumulate_moments_2d(tt_idx,ntt,nTot,X,Y,Z,-a,a,-b,b,means,vars,skews,&
                         kurts,nby,nbz,means_sl,vars_sl,skews_sl,kurts_sl)
               call median(nTot,X,medians(tt_idx))

               
               call make_histogram(nTot,X,nhb,hist_centers(tt_idx,1:nhb),hist_heights(tt_idx,1:nhb))
               
               if (save_hist2d) then
                    call make_histogram2d(nTot,X,Y,nbx,nby,hist2dcx(tt_idx,1:nbx),&
                              hist2dcy(tt_idx,1:nby), hist2d(tt_idx,1:nbx,1:nby))
               end if
			   
               !
               ! Write history if requested.
               !
               ! Need to buffer writes to the hard drive so that we don't lock up the 
               ! computation with file opens/closes. Ideally the buffer should be as 
               ! large as possible while fitting into RAM; modify the relevant parameter 
               ! at the end of the variable definitions.
               !
               
               if (save_hist) then
                    call buffer_op_duct(bk,nTot,buffer_len,Xbuffer,Ybuffer,Zbuffer,&
                                   X,Y,Z,ntt,inext,dset_id_X,dset_id_Y,dset_id_Z)
               end if
               
               !
               ! Update the target time and array index.
               !
               if (next_tt .lt. Tfinal) then
                    tt_idx = tt_idx + 1
                    next_tt = target_times(tt_idx)
               end if
               
          end if

          ! Display percentage progress. The last argument as .true. should be used with gfortran
          ! (or any other compiler that supports "\b"), or .false. with ifort.

          call system_clock(mde_t2,count_rate)    ! Time in milliseconds

          mde_ntc = kt-1
          mde_dts(mde_ntc) = (mde_t2-mde_t1)/dble(count_rate)  ! Time in seconds

          call progress_meter(kt,nt,.true.)

     end do

     out_msg = 'simul_done'
     call duct_mc_messages(out_msg,nz)
     

     if (save_hist) then
     ! Write the remainder of the buffer, then close the file.

          rem = ntt-inext+1
          
          if (rem .gt. 0) then
               call hdf_write_to_open_2d_darray(ntt,nTot,inext,rem,Xbuffer(1:rem,1:nTot),dset_id_X)
               call hdf_write_to_open_2d_darray(ntt,nTot,inext,rem,Ybuffer(1:rem,1:nTot),dset_id_Y)
               call hdf_write_to_open_2d_darray(ntt,nTot,inext,rem,Zbuffer(1:rem,1:nTot),dset_id_Z)
          end if

          call h5dclose_f(dset_id_X,h5error)
          call h5dclose_f(dset_id_Y,h5error)
          call h5dclose_f(dset_id_Z,h5error)
          call h5fclose_f(file_id,h5error)
     end if
     
     ! Because of the nature of hdf5 mod for fortran,
     ! we close the interface here, since it gets
     ! re-opened in the calls below.

     call h5close_f(h5error)

     if (save_hist2d) then
	      write(*,*) nbx,nby,shape(hist2dcx),shape(hist2dcy),shape(hist2d)
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
     !
     ! The ellipse and duct implementations are identical again here, so no use making another subroutine.
     
     call save_the_rest_duct(fname2,geometry,ntt,target_times,means,vars,skews,kurts,medians,nby,nbz,&
                              means_sl,vars_sl,skews_sl,kurts_sl,nhb,hist_centers,hist_heights,&
                              Pe,nTot,mt_seed,aratio,q,dtmax,t_warmup)

     
     ! ------
     deallocate(X,Y,Z)

     deallocate(means,vars,skews,kurts,medians,target_times)
     
     if (.not. (nbins .eq. 0)) then
          deallocate(means_sl,vars_sl,skews_sl,kurts_sl)
     end if
     
     deallocate(hist_centers,hist_heights)
     deallocate(t_hist)
     
     if (save_hist) then
          deallocate(Xbuffer,Ybuffer,Zbuffer)
     end if

     if (save_hist2d) then
          deallocate(hist2d,hist2dcx,hist2dcy)
     end if
	 
     out_msg = 'done'
     call duct_mc_messages(out_msg,nz)
     
1234 continue

end program ellipse_mc
