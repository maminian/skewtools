subroutine save_the_rest_duct(fname,geometry,ntt,t_hist,means,vars,skews,kurts,nby,nbz,&
                              means_sl,vars_sl,skews_sl,kurts_sl,nhb,hist_centers,hist_heights,&
                              Pe,nTot,mt_seed,aratio,q,dtmax,t_warmup)
! Saves the remainder of the calculated data (moments, problem parameters, solver settings, 
! etc) in the h5 file.

implicit none
     integer, parameter                 :: i64 = selected_int_kind(18)
     
     
     character(len=1024), intent(in)                                       :: fname,geometry
     character(len=1024)                                                   :: dsetname, descr
     integer, intent(in)                                                   :: ntt, nby, nbz, nhb, nTot
     double precision, dimension(1:ntt), intent(in)                        :: t_hist,means,vars,skews,kurts
     double precision, dimension(1:ntt,1:nby,1:nbz), intent(in)            :: means_sl,vars_sl,skews_sl,kurts_sl
     double precision, dimension(1:ntt,1:nhb), intent(in)                  :: hist_centers,hist_heights
     double precision, intent(in)                                          :: Pe, aratio, q, dtmax, t_warmup
     
     integer(i64), intent(in)                                              :: mt_seed

     character(len=1024)                                                   :: channel,duct,ellipse
     channel = "channel"
     duct = "duct"
     ellipse = "ellipse"

     ! ---------------------------------------------

     dsetname = "geometry"
     descr = "Problem geometry: Channel=0, Duct=1, Ellipse=2, Other=3"
     if (geometry .eq. channel) then
          call hdf_add_1d_darray_to_file(1,0,fname,dsetname,descr)
     else if (geometry .eq. duct) then
          call hdf_add_1d_darray_to_file(1,1,fname,dsetname,descr)
     else if (geometry .eq. ellipse) then
          call hdf_add_1d_darray_to_file(1,2,fname,dsetname,descr)
     end if


     dsetname = "Time"
     descr = "Nondimensionalized time, t = a^2/kappa t'"
     call hdf_add_1d_darray_to_file(ntt,t_hist,fname,dsetname,descr)

     dsetname = "Avgd_Mean"
     descr = "Cross-section averaged mean for nTot particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,means,fname,dsetname,descr)

     dsetname = "Avgd_Variance"
     descr = "Cross-section averaged variance for nTot particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,vars,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness"
     descr = "Cross-section averaged skewness for nTot particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,skews,fname,dsetname,descr)
	 
     dsetname = "Avgd_Kurtosis"
     descr = "Cross-section averaged skewness for nTot particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,kurts,fname,dsetname,descr)
     
     if (.not. (nby .eq. 0)) then
     
          dsetname = "nBinsY"
          descr = "Number of bins used in Y direction to calculate pointwise statistics."
          call hdf_add_1d_darray_to_file(1,dble(nby),fname,dsetname,descr)
	      
          dsetname = "nBinsZ"
          descr = "Number of bins used in Z direction to calculate pointwise statistics."
          call hdf_add_1d_darray_to_file(1,dble(nbz),fname,dsetname,descr)
          
          dsetname = "Mean"
          descr = "Pointwise mean in the X direction."
          call hdf_add_3d_darray_to_file(ntt,nby,nbz,means_sl,fname,dsetname,descr)
          
          dsetname = "Variance"
          descr = "Pointwise variance in the X direction."
          call hdf_add_3d_darray_to_file(ntt,nby,nbz,vars_sl,fname,dsetname,descr)

          dsetname = "Skewness"
          descr = "Pointwise skewness in the X direction."
          call hdf_add_3d_darray_to_file(ntt,nby,nbz,skews_sl,fname,dsetname,descr)
	      
          dsetname = "Kurtosis"
          descr = "Pointwise kurtosis in the X direction."
          call hdf_add_3d_darray_to_file(ntt,nby,nbz,kurts_sl,fname,dsetname,descr)
          
     end if
     
     dsetname = "Hist_centers"
     descr = "Bin centers for the cross-sectionally averaged distribution."
     call hdf_add_2d_darray_to_file(ntt,nhb,hist_centers,fname,dsetname,descr)

     dsetname = "Hist_heights"
     descr = "Bin heights for the cross-sectionally averaged distribution."
     call hdf_add_2d_darray_to_file(ntt,nhb,hist_heights,fname,dsetname,descr)

     dsetname = "Peclet"
     descr = "Peclet number."
     call hdf_add_1d_darray_to_file(1,Pe,fname,dsetname,descr)
     
     dsetname = "aratio"
     descr = "Aspect ratio of the domain."
     call hdf_add_1d_darray_to_file(1,aratio,fname,dsetname,descr)

     dsetname = "q"
     descr = "Shape parameter (only relevant for racetrack)."
     call hdf_add_1d_darray_to_file(1,q,fname,dsetname,descr)

     dsetname = "nTot"
     descr = "Number of particles used."
     call hdf_add_1d_darray_to_file(1,dble(nTot),fname,dsetname,descr)
     
     dsetname = "mt_seed"
     descr = "Integer seed used in the Mersenne Twister (RNG)."
     call hdf_add_1d_darray_to_file(1,dble(mt_seed),fname,dsetname,descr)
     
     dsetname = "timesteps"
     descr = "Number of large timesteps."
     call hdf_add_1d_darray_to_file(1,dble(ntt),fname,dsetname,descr)

     dsetname = "dtmax"
     descr = "Maximum internal timestep."
     call hdf_add_1d_darray_to_file(1,dtmax,fname,dsetname,descr)

     dsetname = "t_warmup"
     descr = "Duration initial condition was let sit before turning on the flow."
     call hdf_add_1d_darray_to_file(1,t_warmup,fname,dsetname,descr)

     dsetname = "nhb"
     descr = "Number of bins used for the cross-sectionally averaged distribution."
     call hdf_add_1d_darray_to_file(1,dble(nhb),fname,dsetname,descr)

end subroutine save_the_rest_duct
