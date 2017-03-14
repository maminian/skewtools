subroutine set_initial_conds_triangle_mc(ny,nz,x0n,a,nGates,nTot,X,Y,Z,y0,z0,x0width,t_warmup,&
                    use_external_ic,ic_file,nl,lls)
! Given all the data, specify the initial conditions in the arrays X,Y,Z.

!use mtmod
use mod_readbuff

implicit none
     
     double precision, intent(in)                                :: y0,z0,a,x0width,t_warmup
     integer, intent(in)                                         :: ny,nz,x0n,nGates,nTot,nl
     double precision, dimension(nl,3), intent(in)               :: lls
     logical, intent(in)                                         :: use_external_ic
     character(len=1024), intent(in)                             :: ic_file
     
     double precision, dimension(nTot), intent(out)              :: X,Y,Z

     
     ! Internal
     integer                                                     :: idx,ix,iy,iz
     double precision                                            :: dist,dx,dy,dz,xl,yl,zl

     integer                                                     :: nsteps,i
     double precision                                            :: dtw
     double precision, dimension(3)                              :: tempv
     double precision, dimension(nl)                             :: bvals

     character(len=1024)                                         :: dsetname

     logical cond
     double precision    :: zr,yr

     parameter(nsteps=10)
     
     ! Advection/diffusion functions!     
     external  :: impose_reflective_BC_polygon, u_triangle


     if (.not. use_external_ic) then
          if (x0n .gt. 1) then
               dx = x0width/dble(x0n-1)
               xl = -x0width/2
          else
               dx = 0.0d0
               xl = 0.0d0
          end if
              
          if (nGates .gt. 1) then
               yl = -1.0d0
               dy = a*2.0d0*dsqrt(3.0d0)/(ny-1)
               zl = -a*dsqrt(3.0d0)
               dz = a*2.0d0*dsqrt(3.0d0)/(nz-1)

               yr = -1.0d0 + a*2.0d0*dsqrt(3.0d0)
               zr = a*dsqrt(3.0d0)
          else
               yl = y0
               dy = 0.0d0
               zl = z0
               dz = 0.0d0
          end if
          
          
          idx = 0
          do iz=0,nz-1
               do iy=0,ny-1
                    
                    tempv(1) = 1.0d0
                    tempv(2) = yl + iy*dy
                    tempv(3) = zl + iz*dz

                    call matvec(nl,3,lls,tempv,bvals)
                    
                    if (all(bvals .ge. 0.0d0)) then

                         do ix=0,x0n-1
                              
                              idx = idx + 1
                              X(idx) = xl + ix*dx
                              Y(idx) = tempv(2)
                              Z(idx) = tempv(3)

                         end do
                    end if
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

          dsetname = "Z"
          call hdf_read_1d_darray(nTot,ic_file,dsetname)
          Z = readbuff_double
          deallocate(readbuff_double)
     end if

!     write(*,*) minval(Y),maxval(Y)
!     write(*,*) minval(Z),maxval(Z)

     ! Diffuse the initial condition by calling the advection diffusion operator 
     ! with Pe = 0.

!     write(*,*) t_warmup

     if (t_warmup .gt. 0.0d0) then
          
          dtw = t_warmup/nsteps

          do i=1,nsteps
               call apply_advdiff1_triangle(nTot,X,Y,Z,0.0d0,dtw,a, &
                              u_triangle,impose_reflective_BC_polygon,nl,lls)
          end do

     end if

!     write(*,*) minval(Y),maxval(Y)
!     write(*,*) minval(Z),maxval(Z)

end subroutine set_initial_conds_triangle_mc
