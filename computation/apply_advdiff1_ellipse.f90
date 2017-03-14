subroutine apply_advdiff1_ellipse(nTot,xv,yv,zv,Pe,dt,a,b, & 
                    flow,reflector,maxrefl)

! Does the basic advection diffusion operation in the ellipse.
!
! This is what has been used up to this point (11 May 2016); 
! essentially what's been done is an operator splitting where 
! advection operator is done first, then diffusion operator.
!
! Arrays X,Y,Z (dimension n) 
! scalars Pe, dt, a, b
!
! double precision function flow, 
! subroutine reflector.
!

use mtmod

implicit none

     ! Inputs/outputs
     integer, intent(in)                                         :: nTot,maxrefl
     double precision, dimension(nTot), intent(inout)            :: xv,yv,zv
     double precision, intent(in)                                :: Pe,dt,a,b

     ! Internal
     double precision, dimension(3,nTot)                         :: W
     double precision                                            :: mcvar,yprev,zprev
     integer                                                     :: i

     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number and type of arguments.
     interface
          double precision function flow(a1,a2,a3,a4)
               implicit none
               double precision :: a1,a2,a3,a4
          end function flow

          subroutine reflector(a1,a2,a3,a4,a5,a6,a7)
               implicit none
               double precision, intent(out) :: a1,a2
               double precision, intent(in)  :: a3,a4,a5,a6
               integer, intent(in)           :: a7
          end subroutine reflector
     end interface


     ! Generate the proper white noise in advance.
     ! Note in the ellipse that the variance of the white noise is the 
     ! same in all directions because the nondimensionalization is 'isotropic'.
     mcvar = 2.0d0*dt

     call my_normal_rng(nTot,W(1,1:nTot),0.0d0,mcvar)
     call my_normal_rng(nTot,W(2,1:nTot),0.0d0,mcvar)
     call my_normal_rng(nTot,W(3,1:nTot),0.0d0,mcvar)
     
     do i=1,nTot
          ! Advection, then diffusion.
          xv(i) = xv(i) + Pe*flow(yv(i),zv(i),a,b)*dt + W(1,i)

          yprev = yv(i)
          zprev = zv(i)
          yv(i) = yv(i) + W(2,i)
          zv(i) = zv(i) + W(3,i)

          call reflector(yv(i),zv(i),yprev,zprev,a,b,maxrefl)
     end do

end subroutine apply_advdiff1_ellipse
