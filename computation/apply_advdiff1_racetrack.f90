subroutine apply_advdiff1_racetrack(nTot,xv,yv,zv,Pe,dt,aratio,q, & 
                    flow,reflector,maxrefl)

! Does the basic advection diffusion operation in the racetrack.
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
     double precision, intent(in)                                :: Pe,dt,aratio,q

     ! Internal
     double precision, dimension(3,nTot)                         :: W
     double precision                                            :: mcvar,yprev,zprev,y1,z1
     integer                                                     :: i

     double precision bdistfun_rt
     
     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number and type of arguments.
     interface
          double precision function flow(a1,a2,a3,a4)
               implicit none
               double precision :: a1,a2,a3,a4
          end function flow

          subroutine reflector(a1,a2,a3,a4,a5,a6,a7,a8,a9)
               implicit none
               double precision, intent(out) :: a1,a2
               double precision, intent(in)  :: a3,a4,a5,a6,a7,a8
               integer, intent(in)           :: a9
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
          xv(i) = xv(i) + Pe*flow(yv(i),zv(i),aratio,q)*dt + W(1,i)

          yprev = yv(i)
          zprev = zv(i)
          y1 = yv(i) + W(2,i)
          z1 = zv(i) + W(3,i)

          call reflector(yv(i),zv(i),y1,z1,yprev,zprev,aratio,q,maxrefl)
!          if (bdistfun_rt(yv(i),zv(i),aratio,q) .lt. 0.0d0) then
!               write(*,*) yv(i),zv(i),bdistfun_rt(yv(i),zv(i),aratio,q)
!          end if
     end do

end subroutine apply_advdiff1_racetrack
