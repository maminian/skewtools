subroutine impose_reflective_BC_polygon(p1y,p1z,p0y,p0z,nl,lls)
! Imposes reflective boundary
! conditions for the Monte Carlo simulation 
! for a general (convex) polygonal geometry.
! 
! The boundaries are specified as a set of 
! linear equations; the interior is described 
! as when all of them are positive.
!
! Boundary points are where any of them are zero.
!
! Inputs:
! 
! double precision,                     :: p0y,p0z
! integer                               :: nl
! double precision, dimension(nl,3)     :: lls
!
! Input/output:
!
! double precision, dimension(2)        :: p1
!

!use mod_triangle_bdry

implicit none
     double precision, intent(in)                      :: p0y,p0z
     double precision, intent(inout)                   :: p1y,p1z
     integer, intent(in)                               :: nl
     double precision, dimension(nl,3), intent(in)     :: lls

     double precision, dimension(2)                    :: pb,mp,p0,p1
     double precision, dimension(3)                    :: tempv,line
     double precision, dimension(nl)                   :: bvals
     logical, dimension(nl)                            :: bcond
     integer, dimension(nl)                            :: bidxs
!     double precision, dimension(nl)                   :: ctimes
     integer                                           :: nbi,bidx,i,minidx
     double precision                                  :: minct,ct

     double precision crosstime

     p0(1) = p0y
     p0(2) = p0z
     p1(1) = p1y
     p1(2) = p1z

     mp(1) = p1y
     mp(2) = p1z

     tempv(1) = 1.0d0
     tempv(2) = mp(1)
     tempv(3) = mp(2)

     call matvec(nl,3,lls,tempv,bvals)
     bcond = (bvals .lt. 0.0d0)


     ! Loop until there are no boundary crossings.
     do while (any(bcond))
!          write(*,*) "impose" 
!          write(*,*) lls(1,:)
!          write(*,*) lls(2,:)
!          write(*,*) lls(3,:)
!          write(*,*) " "
!          write(*,*) bcond
!          write(*,*) tempv
!          write(*,*) bvals
!          write(*,*) "moo"
!          read(*,*)
          ! Find which boundaries have been crossed.
          call findcond(nl,bcond,nbi,bidxs)

          ! Find the time of crossing on each crossed 
          ! boundary. Take the boundary crossed first.

          line = lls(bidxs(1),:)
          minct = crosstime(pb,mp,line)

          bidx = bidxs(1)
          do i=2,nbi
               line = lls(bidxs(i),:)
               ct = crosstime(pb,mp,line)
               if (ct .lt. minct) then
                    minct = ct
                    bidx = bidxs(i)
               end if
          end do
          
          ! Reflect across this boundary.
          line = lls(bidx,:)
          call reflector(pb,mp,line)

          tempv(2) = mp(1)
          tempv(3) = mp(2)

          ! Re-evaluate the new position.
          call matvec(nl,3,lls,tempv,bvals)
          bcond = (bvals .lt. 0.0d0)

     end do
     
     p1y = tempv(2)
     p1z = tempv(3)

end subroutine impose_reflective_BC_polygon
!
! ------------------------------------
!
double precision function crosstime(p,q,l)
! Calculates the time of intersection through line l 
! in a parameterized path going from point p to q.
implicit none
     double precision, dimension(2), intent(in)   :: p,q
     double precision, dimension(3), intent(in)   :: l
     
     crosstime = -(l(1) + l(2)*p(1) + l(3)*p(2))/(l(2)*(q(1)-p(1))+l(3)*(q(2)-p(2)))

end function crosstime
!
! ------------------------------------
!
subroutine reflector(p0,p1,l)
!
! Reflects the particle that would have gone from p0 to p1 
! across the line l. Should not end up in this subroutine 
! unless this actually happens.
!
! On output, the points are changed to
!
! p0: point of intersection with l
! p1: position after reflection; p1 = p0+v for a vector v.
!
implicit none
     double precision, dimension(2), intent(inout)     :: p0,p1
     double precision, dimension(3), intent(in)        :: l

     double precision, dimension(2)                    :: gradl,pb,v
     double precision                                  :: s

     double precision crosstime

     gradl(1) = l(2)
     gradl(2) = l(3)

     ! Find the time and location of intersection, take the component 
     ! of the vector that's outside the domain
     s = crosstime(p0,p1,l)
     pb = p0 + s*(p1-p0)
     v = p1 - pb

     ! Reflect this vector component across the plane
     call reflect(2,v,gradl)

     ! Update p1 based on this reflection.
     p1 = pb + v
     p0 = pb
     
end subroutine reflector

