subroutine apply_advdiff1_duct(n,xv,yv,zv,Pe,dt,a,b, & 
                    flow,reflector)
! Does the basic advection diffusion operation in the duct.
!
! This is what has been used up to this point (11 May 2016); 
! essentially what's been done is an operator splitting where 
! advection operator is done first, then diffusion operator.
!
! Arrays X,Y,Z (dimension n) 
! scalars Pe, dt, a, b
! double precision function flow, 
! subroutine reflector.
!

use mod_ductflow

implicit none

     ! Inputs/outputs
     integer, intent(in)                                         :: n
     double precision, dimension(n), intent(inout)               :: xv,yv,zv
     double precision, intent(in)                                :: Pe,dt,a,b

     ! Internal
     double precision, dimension(3,n)                            :: W
     double precision                                            :: mcvar
     integer                                                     :: i

     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number and type of arguments.
     interface
          double precision function flow(p,q)
               double precision :: p,q
          end function flow

          subroutine reflector(p,q,r)
               double precision :: p,q,r
          end subroutine reflector
     end interface


     ! Generate the proper white noise in advance.
     mcvar = 2.0d0*dt
     call my_normal_rng(n,W(1,1:n),0.0d0,mcvar)
     call my_normal_rng(n,W(2,1:n),0.0d0,mcvar)
     call my_normal_rng(n,W(3,1:n),0.0d0,mcvar)
     
     do i=1,n
          ! Advection, then diffusion.

          xv(i) = xv(i) + Pe*flow(yv(i),zv(i))*dt + W(1,i)
          yv(i) = yv(i) + W(2,i)
          zv(i) = zv(i) + W(3,i)

          ! Call the generic reflection subroutine passed in.
          call reflector(yv(i),-a,a)
          call reflector(zv(i),-b,b)
     end do

end subroutine apply_advdiff1_duct
