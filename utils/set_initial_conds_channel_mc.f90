subroutine set_initial_conds_channel_mc(ny,nGates,x0n,nTot,X,Y,y0,a,x0width,t_warmup,use_external_ic,ic_file)
! The purpose of the subroutine is in the name.
!
! y0 refers to point-source initial condition at location y0.
use mtmod
use mod_readbuff

implicit none
     
     integer, intent(in)                                         :: ny,nGates,x0n,nTot
     double precision, dimension(nTot), intent(inout)            :: X,Y
     
     double precision, intent(in)                                :: y0,a,x0width,t_warmup
     logical, intent(in)                                         :: use_external_ic
     character(len=1024), intent(in)                             :: ic_file

     ! Internal
     integer                                                     :: idx,i,j,k,nsteps
     double precision                                            :: xl,yl,dx,dy,dtw
     character(len=1024)                                         :: dsetname

     parameter(nsteps=10)
     ! Advection/diffusion functions!     
     external  :: impose_reflective_BC_rect, u_dummy

     ! ------------------------------------------------------------

     if (.not. use_external_ic) then
          ! Construct the initial condition from the parameters specified.
          if (x0n .gt. 1) then 
               dx = x0width/(x0n-1)
               xl = -x0width/2
          else
               dx = 0.0d0
               xl = 0.0d0
          end if

          if (nGates .gt. 1) then
               dy = 2.0d0*a/(ny-1)
               yl = -a
          else
               dy = 0.0d0
               yl = y0
          end if
          
          idx = 0

          do j=0,ny-1
               do k=0,x0n-1
                    idx = idx + 1
                    X(idx) = xl + k*dx
                    Y(idx) = yl + j*dy

               end do
          end do

          
     else
          ! Skip all this and read the x,y,z initial data in from the file.
          dsetname = "X"
          call hdf_read_1d_darray(nTot,ic_file,dsetname)
          X = readbuff_double
          deallocate(readbuff_double)

          dsetname = "Y"
          call hdf_read_1d_darray(nTot,ic_file,dsetname)
          Y = readbuff_double
          deallocate(readbuff_double)

     end if


     ! Diffuse the initial condition by calling the advection diffusion operator 
     ! with Pe = 0.

     if (t_warmup .gt. 0.0d0) then
          
          dtw = t_warmup/nsteps

          do i=1,nsteps
               call apply_advdiff1_chan(nTot,X,Y,0.0d0,dtw,a, &
                              u_dummy,impose_reflective_BC_rect)
          end do

     end if
     
end subroutine set_initial_conds_channel_mc
