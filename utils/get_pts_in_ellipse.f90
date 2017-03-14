subroutine get_pts_in_ellipse(ny,nz,x0n,a,b,nGates,nTot)
! In the ellipse, the technique to set the initial condition is to proceed with 
! uniform spacing as if we were in the rectangle, then exclude points that lie 
! outside the circle. This leaves an issue of not knowing exactly how many 
! points will be remaining. This function is a trimmed down version of the 
! set_initial_conditions_ellipse where only the *number* of points in the domain 
! is counted.
!
! Could probably be done in a single call if I was clever. But this isn't a bottleneck 
! in computations.

implicit none

     integer, intent(in)                                    :: ny,nz,x0n
     double precision, intent(in)                           :: a,b
     integer, intent(out)                                   :: nGates,nTot


     integer                                                :: idx,iy,iz
     double precision                                       :: hy,hz,dist


     ! Sample points in the circumscribing square; throw out any that 
     ! lie outside the circle.
     
     hy = 2.0d0*a/dble(ny-1)
     hz = 2.0d0*b/dble(nz-1)
     idx = 0

     if (nGates .gt. 1) then
          do iz=0,nz-1
               do iy=0,ny-1
                    dist = ((-a + iy*hy)**2)/(a**2) + ((-b + iz*hz)**2)/(b**2)

                    if (dist .le. 1.0d0) then
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

end subroutine get_pts_in_ellipse
