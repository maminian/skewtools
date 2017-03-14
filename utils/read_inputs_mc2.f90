subroutine read_inputs_mc2(param_file)
! Reads the parameter input file specified on the calling of the program, which 
! sets all of the problem and timestepping parameters.

use mod_parameters
use mod_time

implicit none
!     integer, parameter                 :: i64 = selected_int_kind(18)
     
     
     character(len=1024), intent(in)    :: param_file
     
!     character(len=1024), intent(out)   :: tstep_type, other_file, ic_file
     
!     double precision, intent(out)      :: aratio, q, Pe, y0, z0, dt, dtmax, Tfinal, x0width
!     logical, intent(out)               :: save_hist, save_hist2d, use_external_ic
!     integer, intent(out)               :: nGates, tsteps, x0n, n_bins
!     integer(i64), intent(out)          :: mt_seed
!     double precision, intent(out)      :: t_warmup
     
     ! Internal
     integer                            :: funit,i
     character(len=1)                   :: dummy
     
     funit=55
     ! Initial condition/problem parameters.

     open(funit,file=param_file)
          ! Six initial lines to skip

          do i=1,6
               read(funit,*) dummy
          end do

          read(funit,*) aratio     ! Aspect ratio (ignored in channel)
          read(funit,*) q          ! Shape parameter (for racetrack)
          read(funit,*) Pe         ! Peclet number
          read(funit,*) nGates     ! Approximate number of discretization points to use in the transverse direction.
          read(funit,*) x0n        ! Number of discretization points to use in the longitudinal direction.
                                   ! Total number of points will be (roughly) nGates*x0n.
                                   ! Memory requirement is then to leading order 3*nGates*x0n*8 bytes.
          read(funit,*) x0width    ! IC characteristic width (relative to short side length 2)
          read(funit,*) y0         ! If nGates=1, specify initial y position for point source release.
          read(funit,*) z0         ! If nGates=1, specify initial z position for point source release.
          read(funit,*) save_hist  ! Flag to save full particle position histories.
          read(funit,*) nbins      ! Number of bins to use in the short direction for pointwise statistics. 
                                   ! If zero, not saved. Computationally expensive.
          read(funit,*) save_hist2d  ! Flag to save 2d histogram data.

          read(funit,*) t_warmup   ! After setting the initial condition with the given 
                                   ! parameters above, allow it to diffuse with no 
                                   ! flow for nondimensional time t_warmup. 
                                   ! Used for setting a plug IC which is Gaussian in x, or 
                                   ! diffusing a point source injection.
                                   ! No statistics are collected during this time.

          read(funit,*) use_external_ic      ! Whether to look at an external h5 file for the initial condition.
          read(funit,*) ic_file              ! Name of the h5 file to look at if use_external_ic is true.

          ! Five lines to skip between the initial condition and timestepping
          do i=1,5
               read(funit,*) dummy
          end do

     ! Timestepping parameters.
     ! Basically, you specify the values of time you want data to be saved, 
     ! and the internal timesteps will be chosen as min(t_{n}-t,dtmax).

          read(funit,*) tstep_type ! 'unif', 'expo', or 'supplied'
          read(funit,*) dt         ! Timestep to save statistics for 'unif', initial time for 'expo'
          read(funit,*) dtmax      ! Maximum internal timestep
          read(funit,*) Tfinal     ! Final time
          read(funit,*) ntt        ! Total number of timesteps to save.
          read(funit,*) other_file ! If tstep_type=='supplied', the name of the h5 file containing the specified times.
          read(funit,*) mt_seed    ! Seed for the RNG. Usually filled in with a random 
                                   ! seed by batch_submit.py.

          ! Skip 5 lines
          do i=1,5
               read(funit,*) dummy     
          end do

          ! Parameters for other setups -- oscillatory flow, absorbing boundaries in the future.
          read(funit,*) use_oscillatory
     
     close(funit)
     
     
end subroutine read_inputs_mc2
