subroutine apply_advdiff1_chan(n,xv,yv,Pe,dt,a,flow,reflector)
! Does the basic advection diffusion operation in the channel.
!
! This is what has been used up to this point (11 May 2016); 
! essentially what's been done isan operator splitting where 
! advection operator is done first, then diffusion operator.
!
! Arrays X,Y (dimension n) 
! scalars Pe, dt, a, 
! double precision function flow, 
! subroutine reflector.
!

implicit none

     ! Inputs/outputs
     integer, intent(in)                                         :: n
     double precision, dimension(n), intent(inout)               :: xv,yv
     double precision, intent(in)                                :: Pe,dt,a

     ! Internal
     double precision, dimension(2,n)                            :: W
     double precision                                            :: mcvar
     integer                                                     :: i

     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number and type of arguments.
     interface
          double precision function flow(p,q,r)
               double precision :: p,q,r
          end function flow

          subroutine reflector(p,q,r)
               double precision :: p,q,r
          end subroutine reflector
     end interface


     ! Generate the proper white noise in advance.
     mcvar = 2.0d0*dt
     call my_normal_rng(n,W(1,1:n),0.0d0,mcvar)
     call my_normal_rng(n,W(2,1:n),0.0d0,mcvar)
     
     do i=1,n
          ! Advection, then diffusion.
          xv(i) = xv(i) + Pe*flow(yv(i),-a,a)*dt + W(1,i)
          yv(i) = yv(i) + W(2,i)
          call reflector(yv(i),-a,a)
     end do

end subroutine apply_advdiff1_chan
