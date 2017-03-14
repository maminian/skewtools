subroutine set_initial_conds_ellipse_mc(ny,nz,x0n,a,b,nGates,nTot,X,Y,Z,y0,z0,x0width,t_warmup)
! Given all the data, specify the initial conditions in the arrays X,Y,Z.

implicit none
     
     double precision, intent(in)                                :: y0,z0,a,b,x0width,t_warmup
     integer, intent(in)                                         :: ny,nz,x0n,nGates,nTot
     double precision, dimension(nTot), intent(out)              :: X,Y,Z
     
     ! Internal
     integer                                                     :: idx,ix,iy,iz
     double precision                                            :: dist,dx,dy,dz,xl,yl,zl

     integer                                                     :: nsteps,i
     double precision                                            :: dtw

     parameter(nsteps=10)
     
     ! Advection/diffusion functions!     
     external  :: impose_reflective_BC_ellipse, u_ellipse

     if (x0n .gt. 1) then
          dx = x0width/dble(x0n-1)
          xl = -x0width/2
     else
          dx = 0.0d0
          xl = 0.0d0
     end if

     if (nGates .gt. 1) then
          yl = -a
          dy = 2.0d0*a/(ny-1)
          zl = -b
          dz = 2.0d0*b/(nz-1)
     else
          yl = y0
          dy = 0.0d0
          zl = z0
          dz = 0.0d0
     end if


     idx = 0
     do iz=0,nz-1
          do iy=0,ny-1
               dist = ((yl + iy*dy)**2)/(a**2) + ((zl + iz*dz)**2)/(b**2)
               if (dist .le. 1.0d0) then
                    do ix=0,x0n-1
                         
                         idx = idx + 1
                         X(idx) = xl + ix*dx
                         Y(idx) = yl + iy*dy
                         Z(idx) = zl + iz*dz
                         
                    end do
               end if
          end do
     end do

     ! Diffuse the initial condition by calling the advection diffusion operator 
     ! with Pe = 0.
     if (.true.) then
          if (t_warmup .gt. 0.0d0) then
               
               dtw = t_warmup/nsteps

               do i=1,nsteps
                    call apply_advdiff1_ellipse(nTot,X,Y,Z,0.0d0,dtw,a,b, &
                                   u_ellipse,impose_reflective_BC_ellipse,floor(10*dtw))
               end do

          end if
     else
          call apply_advdiff1_ellipse(nTot,X,Y,Z,0.0d0,t_warmup,a,b, &
                              u_ellipse,impose_reflective_BC_ellipse,floor(10*dtw))
     end if

end subroutine set_initial_conds_ellipse_mc
