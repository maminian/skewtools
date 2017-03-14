subroutine apply_advdiff_chan_robin(n,xv,yv,Pe,dt,a,flow,reflector)
! Advection diffusion operator, now with Robin boundary conditions.
!
! Arrays X,Y (dimension n) 
! scalars Pe, dt, a, 
! double precision function flow, 
! subroutine reflector.
!
! Robin boundary conditions are done by:
! 
!    (1) Introducing a parameter pexit indicating the threshold
!        required for a particle to exit. For instance, 
!         if pexit=0.6, when a particle leaves the domain, 
!         it has a probability 0.6 to permanently leave the 
!         domain when it exits (it becomes inactive);
!    (2) Using an integer array indicating the location of the 
!         next active particle. Initially this is the array
!         (2,3,...,nTot,0), but changes when particles become 
!         inactive (this is a linked list, I think). 
!         It may also be possible to do this by 
!         shuffling array entries (moving inactives to the 
!         bottom), but then tracking individual paths becomes 
!         impossible, and may bottleneck things.

use mod_active_particles

implicit none

     ! Inputs/outputs
     integer, intent(in)                                         :: n
     double precision, dimension(n), intent(inout)               :: xv,yv
     double precision, intent(in)                                :: Pe,dt,a

     ! Internal
     double precision, dimension(:,:), allocatable               :: W
     double precision                                            :: mcvar
     integer                                                     :: i,idx,k
     logical                                                     :: absorbed
     
     double precision, parameter                                 :: uwall = 2.0d0/3.0d0

     ! Interface necessary for the passed function and subroutine;
     ! only specifies the number and type of arguments.
     interface
          double precision function flow(p,q,r)
               double precision :: p,q,r
          end function flow

          subroutine reflector(p,q,r,s)
               double precision    :: p,q,r
               logical             :: s
          end subroutine reflector
     end interface

     ! -----------------------
     !
     allocate(W(2,nactive))

     ! Generate the proper white noise in advance.
     mcvar = 2.0d0*dt
     call my_normal_rng(nactive,W(1,1:nactive),0.0d0,mcvar)
     call my_normal_rng(nactive,W(2,1:nactive),0.0d0,mcvar)
     
     k = 1
     do i=1,nactive
          
          idx = aset(i)
          
          xv(idx) = xv(idx) + Pe*flow(yv(idx),-a,a)*dt + W(1,i) + Pe*uwall*dt
          yv(idx) = yv(idx) + W(2,i)

          call reflector(yv(idx),-a,a,absorbed)
          
          ! Remove from the active set if the particle has been absorbed.
          if (absorbed) then
               nactive = nactive - 1
!               write(*,*) xv(idx), yv(idx)
          else
               asettemp(k) = idx
               k = k+1
          end if
     end do
     aset = asettemp


     deallocate(W)

end subroutine apply_advdiff_chan_robin
