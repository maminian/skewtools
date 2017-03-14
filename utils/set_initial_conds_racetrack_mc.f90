subroutine set_initial_conds_racetrack_mc(x0n,aratio,q,nGates,nTot,X,Y,Z,y0,z0,x0width,t_warmup)
! Given all the data, specify the initial conditions in the arrays X,Y,Z.

use mtmod

implicit none
     
     double precision, intent(in)                                :: y0,z0,aratio,q,x0width,t_warmup
     integer, intent(in)                                         :: x0n,nGates,nTot
     double precision, dimension(nTot), intent(out)              :: X,Y,Z
     
     ! Internal
     integer                                                     :: idx,ix,iy,iz
     double precision                                            :: dist,dx,dy,dz,xl,yl,zl,yw,zw

     integer                                                     :: nsteps,i,j
     double precision                                            :: dtw

     ! Advection/diffusion functions!     
     external  :: impose_reflective_BC_racetrack, u_racetrack
     
     double precision bdistfun_rt

     parameter(nsteps=10)
     
     ! For the moment, do things differently: Just do random placings with a 
     ! rejection method in the transverse coordinates.
     if (x0width .eq. 0.0d0) then
          xl = 0.0d0
          dx = 0.0d0
     else
          xl = -x0width/2.0d0
          dx = x0width/(x0n-1)
     end if

     if (nGates .le. 1) then
          idx = 0
          do i=0,x0n-1
               do j=1,nGates
                    idx = idx + 1

                    X(idx) = xl + dx*i

                    Y(idx) = y0
                    Z(idx) = z0

               end do
          end do
     else
          yl = -1.3d0
          yw = -2*yl
          zl = -1.3d0/aratio
          zw = -2*zl

          idx = 0
          do i=0,x0n-1
               do j=1,nGates
                    idx = idx + 1

                    X(idx) = xl + dx*i

                    Y(idx) = yl + yw*grnd()
                    Z(idx) = zl + zw*grnd()
                    do while (bdistfun_rt(Y(idx),Z(idx),aratio,q) .lt. 0.0d0)
                         Y(idx) = yl + yw*grnd()
                         Z(idx) = zl + zw*grnd()
                    end do
               end do
          end do

     end if
     
     ! Diffuse the initial condition by calling the advection diffusion operator 
     ! with Pe = 0.
     if (t_warmup .gt. 0.0d0) then
          
          dtw = t_warmup/nsteps

          do i=1,nsteps
               call apply_advdiff1_ellipse(nTot,X,Y,Z,0.0d0,dtw,aratio,q, &
                              u_racetrack,impose_reflective_BC_racetrack,floor(10*dtw))
          end do

     end if

end subroutine set_initial_conds_racetrack_mc
