subroutine get_pts_in_triangle(ny,nz,x0n,a,nGates,nTot,nl,lls)
! In the triangle, the technique to set the initial condition is to proceed with 
! uniform spacing as if we were in the rectangle, then exclude points that lie 
! outside the triangle. This leaves an issue of not knowing exactly how many 
! points will be remaining. This function is a trimmed down version of the 
! set_initial_conditions_triangle where only the *number* of points in the domain 
! is counted.
!
! Could probably be done in a single call if I was clever. But this isn't a bottleneck 
! in computations.

implicit none

     integer, intent(in)                                    :: ny,nz,x0n,nl
     double precision, intent(in)                           :: a
     double precision, dimension(nl,3), intent(in)          :: lls
     integer, intent(out)                                   :: nGates,nTot


     integer                                                :: idx,iy,iz
     double precision                                       :: hy,hz,rl,rb
     double precision, dimension(3)                         :: tempv
     double precision, dimension(nl)                        :: bvals


     ! Sample points in the circumscribing square; throw out any that 
     ! lie outside the triangle.

     rl = -a*dsqrt(3.0d0)
     rb = -1.0d0

     hz = 2.0d0*a*dsqrt(3.0d0)/(nz-1)
     hy = 2.0d0*a*dsqrt(3.0d0)/(ny-1)
     idx = 0

     if (nGates .gt. 1) then
          do iz=0,nz-1
               do iy=0,ny-1
                    
                    tempv(1) = 1.0d0
                    tempv(2) = rb + iy*hy
                    tempv(3) = rl + iz*hz

                    call matvec(nl,3,lls,tempv,bvals)

                    if (all(bvals .ge. 0.0d0)) then
                         idx = idx + 1
                    end if
               end do
          end do
          
          nGates = idx
     else
          ! Nothing gets changed.
          nGates = 1
     end if

     nTot = nGates*x0n

!     write(*,*) nTot
end subroutine get_pts_in_triangle
