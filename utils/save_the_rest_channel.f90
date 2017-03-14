subroutine save_the_rest_channel(fname,geometry,ntt,target_times,means,vars,skews,kurts,n_bins,&
                              means_sl,vars_sl,skews_sl,kurts_sl,nhb,hist_centers,hist_heights,&
                              Pe,nGates,x0n,x0width,mt_seed)
! Saves the remainder of the calculated data (moments, problem parameters, solver settings, 
! etc) in the h5 file.

implicit none
     integer, parameter                 :: i64 = selected_int_kind(18)
     
     
     character(len=1024), intent(in)                             :: fname,geometry
     character(len=1024)                                         :: dsetname, descr
     integer, intent(in)                                         :: ntt, n_bins, nGates, x0n, nhb
     double precision, dimension(1:ntt), intent(in)              :: target_times,means,vars,skews,kurts
     double precision, dimension(1:ntt,1:n_bins), intent(in)     :: means_sl,vars_sl,skews_sl,kurts_sl
     double precision, dimension(1:ntt,1:nhb), intent(in)        :: hist_centers,hist_heights
     double precision, intent(in)                                :: Pe, x0width
     
     integer(i64), intent(in)                                    :: mt_seed

     character(len=1024)                                         :: channel,duct,ellipse
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
     descr = "Nondimensionalized time, t = \\frac{H^2}{\\kappa} t'"
     call hdf_add_1d_darray_to_file(ntt,target_times,fname,dsetname,descr)

     
     dsetname = "Avgd_Mean"
     descr = "Cross-section averaged mean for Nwalkers particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,means,fname,dsetname,descr)

     dsetname = "Avgd_Variance"
     descr = "Cross-section averaged variance for Nwalkers particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,vars,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness"
     descr = "Cross-section averaged skewness for Nwalkers particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,skews,fname,dsetname,descr)
	 
     dsetname = "Avgd_Kurtosis"
     descr = "Cross-section averaged kurtosis for Nwalkers particles, in X direction."
     call hdf_add_1d_darray_to_file(ntt,kurts,fname,dsetname,descr)
     
     if (.not. (n_bins .eq. 0)) then
          
          dsetname = "nBins"
          descr = "Number of bins used to calculate statistics across slices."
          call hdf_add_1d_darray_to_file(1,dble(n_bins),fname,dsetname,descr)
          
          dsetname = "Mean"
          descr = "Mean on slices for Nwalkers particles, in X direction."
          call hdf_add_2d_darray_to_file(ntt,n_bins,means_sl,fname,dsetname,descr)
          
          dsetname = "Variance"
          descr = "Variance on slices for Nwalkers particles, in X direction."
          call hdf_add_2d_darray_to_file(ntt,n_bins,vars_sl,fname,dsetname,descr)

          dsetname = "Skewness"
          descr = "Skewness on slices for Nwalkers particles, in X direction."
          call hdf_add_2d_darray_to_file(ntt,n_bins,skews_sl,fname,dsetname,descr)

          dsetname = "Kurtosis"
          descr = "Kurtosis on slices for Nwalkers particles, in X direction."
          call hdf_add_2d_darray_to_file(ntt,n_bins,kurts_sl,fname,dsetname,descr)
     
     end if
     
     dsetname = "Hist_centers"
     descr = "Bin centers for the cross-sectionally averaged distribution."
     call hdf_add_2d_darray_to_file(ntt,nhb,hist_centers,fname,dsetname,descr)

     dsetname = "Hist_heights"
     descr = "Bin heights for the cross-sectionally averaged distribution (normalized to PDF)."
     call hdf_add_2d_darray_to_file(ntt,nhb,hist_heights,fname,dsetname,descr)

     dsetname = "Peclet"
     descr = "Peclet number."
     call hdf_add_1d_darray_to_file(1,Pe,fname,dsetname,descr)
     
     dsetname = "x0width"
     descr = "Initial conditionn longitudinal width."
     call hdf_add_1d_darray_to_file(1,x0width,fname,dsetname,descr)

     dsetname = "nGates"
     descr = "Number of random walkers used in this trial."
     call hdf_add_1d_darray_to_file(1,dble(nGates),fname,dsetname,descr)
     
     dsetname = "mt_seed"
     descr = "Integer seed used in the Mersenne Twister (RNG)."
     call hdf_add_1d_darray_to_file(1,dble(mt_seed),fname,dsetname,descr)
     
     dsetname = "timesteps"
     descr = "Number of timesteps."
     call hdf_add_1d_darray_to_file(1,dble(ntt),fname,dsetname,descr)
     
     dsetname = "x0n"
     descr = "Number of discretization points discretizing the longitudinal IC."
     call hdf_add_1d_darray_to_file(1,dble(x0n),fname,dsetname,descr)

     dsetname = "nhb"
     descr = "Number of bins used for the cross-sectionally averaged distribution."
     call hdf_add_1d_darray_to_file(1,dble(nhb),fname,dsetname,descr)

end subroutine save_the_rest_channel
