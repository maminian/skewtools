subroutine apply_advdiff1_triangle(nTot,xv,yv,zv,Pe,dt,a, & 
                    flow,reflector,nl,lls)

! Does the basic advection diffusion operation in the triangle.
!
! This is what has been used up to this point (11 May 2016); 
! essentially what's been done is an operator splitting where 
! advection operator is done first, then diffusion operator.
!
! Arrays X,Y,Z (dimension n) 
! scalars Pe, dt, a
!
! double precision function flow, 
! subroutine reflector.
!

use mtmod

implicit none

     ! Inputs/outputs
     integer, intent(in)                                         :: nTot
     double precision, dimension(nTot), intent(inout)            :: xv,yv,zv
     double precision, intent(in)                                :: Pe,dt,a
     integer, intent(in)                                         :: nl
     double precision, dimension(nl,3), intent(in)               :: lls

     ! Internal
     double precision, dimension(3,nTot)                         :: W
     double precision                                            :: mcvar,yprev,zprev
     integer                                                     :: i

     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number of arguments and their type.
     interface
          double precision function flow(a1,a2,a3)
               implicit none
               double precision    :: a1,a2,a3
          end function flow

          subroutine reflector(a1,a2,a3,a4,a5,a6)
               implicit none
               double precision, intent(inout)                :: a1,a2
               double precision, intent(in)                   :: a3,a4
               integer, intent(in)                            :: a5
               double precision, dimension(a5,3), intent(in)  :: a6
          end subroutine reflector
     end interface


     ! Generate the proper white noise in advance.
     ! Note in the ellipse that the variance of the white noise is the 
     ! same in all directions because the nondimensionalization is 'isotropic'.
     mcvar = 2.0d0*dt

!     write(*,*) dt,mcvar

     call my_normal_rng(nTot,W(1,1:nTot),0.0d0,mcvar)
     call my_normal_rng(nTot,W(2,1:nTot),0.0d0,mcvar)
     call my_normal_rng(nTot,W(3,1:nTot),0.0d0,mcvar)
     
     do i=1,nTot
          ! Advection, then diffusion.
          xv(i) = xv(i) + Pe*flow(yv(i),zv(i),a)*dt + W(1,i)

          yprev = yv(i)
          zprev = zv(i)
          yv(i) = yv(i) + W(2,i)
          zv(i) = zv(i) + W(3,i)
          
!          write(*,*) yv(i),zv(i),yprev,zprev
          call reflector(yv(i),zv(i),yprev,zprev,nl,lls)
!          write(*,*) yv(i),zv(i),yprev,zprev
     end do

end subroutine apply_advdiff1_triangle
