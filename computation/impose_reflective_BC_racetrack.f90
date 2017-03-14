subroutine impose_reflective_BC_racetrack(yout,zout,y1,z1,y0,z0,aratio,q,maxrefl)
! Impose reflective boundaries for the racetrack.
! (yf,zf) is the position after taking a step,
! (y0,z0) is the previous position.

implicit none
     double precision, intent(in)       :: y0,z0,aratio,q
     double precision, intent(in)       :: y1,z1
     double precision, intent(out)      :: yout,zout
     integer, intent(in)                :: maxrefl
     double precision                   :: bdf,l
     double precision, dimension(2)     :: v1,nhat
     
     integer                            :: nrefl,mmr
     logical                            :: flag
     double precision                   :: s,yc,zc,yold,zold,yf,zf


     double precision bdistfun_rt
     
     l = aratio

     ! Idiot-proofing
     mmr = max(maxrefl,4)


     bdf = bdistfun_rt(y1,z1,aratio,q)
     nrefl = 0

     yold = y0
     zold = z0
     yf = y1
     zf = z1

     flag = ((bdf .lt. 0.0d0) .and. (nrefl .lt. mmr))

!     if (bdf .lt. 0.0d0) then
!          write(*,*) "moo",nrefl
!          write(*,*) y0,z0,bdistfun_rt(y0,z0,aratio,q)
!          write(*,*) y1,z1,bdistfun_rt(y1,z1,aratio,q)
!     end if

     do while ( flag )

          ! Find the point of intersection.
          ! Parameterize the line connecting (y0,z0) to (yf,zf), 
          ! calculate a double s indicating point of intersection.
          ! Essentially a 1D calculation, should be able to do 
          ! a couple iterations of Newton's method to capture.
          call calc_intersection_pt_rt(yold,zold,yf,zf,aratio,q,s,.false.)
!          if ( ( s .lt. 0.0d0) .or. (s .gt. 1.0d0) ) then
!               call calc_intersection_pt_rt(yold,zold,yf,zf,aratio,q,s,.true.)
!               write(*,*) "moo",nrefl
!               write(*,*) s
!               write(*,*) yold,zold,bdistfun_rt(yold,zold,aratio,q)
!               write(*,*) yf,zf,bdistfun_rt(yf,zf,aratio,q)
!               read(*,*)
!          end if
          ! Find the coordinates of intersection.
          yc = yold + s*(yf-yold)
          zc = zold + s*(zf-zold)

          ! Construct the component of the vector going out of the domain.
!          write(*,*) "moo"
!          write(*,*) yc,zc,bdistfun_rt(yc,zc,aratio,q)
          v1(1) = yf-yc
          v1(2) = zf-zc
          
          ! Construct the outward normal vector
          call bdistfun_rt_grad(yc,zc,aratio,q,nhat)


          call reflect(2,v1,nhat)

          ! Get new point an update the old point if it's needed.
          yold = yc
          zold = zc

          yf = yc + v1(1)
          zf = zc + v1(2)

          ! Update the distance of the new point and the number of reflections.
          bdf = bdistfun_rt(yf,zf,aratio,q)

          nrefl = nrefl + 1


          flag = ((bdf .lt. 0.0d0) .and. (nrefl .lt. mmr))
!          write(*,*) (bdf .lt. 0.0d0),nrefl,mmr,flag
     end do


     ! Last line of defense.
     if (bdf .lt. 0.0d0) then
          write(*,*) bdf
     end if

     if ((nrefl .eq. mmr) .and. (bdf .lt. 0.0d0)) then
!          write(*,*) "moo"
          yout = y0
          zout = z0
     end if
     
     yout = yf
     zout = zf
!     write(*,*) y0,z0,bdistfun_rt(y0,z0,aratio,q)
!     write(*,*) y1,z1,bdistfun_rt(y1,z1,aratio,q)
!     write(*,*) yout,zout,bdistfun_rt(yout,zout,aratio,q)

     
end subroutine impose_reflective_BC_racetrack

!
! -----------------------------
!
subroutine calc_intersection_pt_rt(y0,z0,yf,zf,aratio,q,s,diagnose)
! Calculates the intersection with the boundary 
! assuming the starting and ending points are in the interior 
! and exterior, appropriately.
!
! Essentially does a few iterations of bisection followed by 
! a few iterations of Newton's method.
!

implicit none
     double precision, intent(in)  :: y0,z0,yf,zf,aratio,q
     double precision, intent(out) :: s
     logical, intent(in)           :: diagnose

     integer, parameter            :: nbi=10 ! Number of bisection iterations
     integer, parameter            :: mnni=5 ! Max number of Newton iterations

     double precision, parameter   :: reltol = 1.0d-4  ! Relative tolerance for cv.
     double precision, parameter   :: abstol = 1.0d-8  ! Absolute tolerance for cv.

     double precision              :: tol
     double precision              :: yl,yr,zl,zr,yc,zc,sl,sr,sc
     double precision              :: vl,vr,vc
     integer                       :: nni,i
     logical                       :: flag

     double precision bdistfun_rt,dderiv_rt
     
     nni = 0

     sl = 0.0d0
     sr = 1.0d0
     sc = 0.5d0

     yl = y0
     zl = z0
     yr = yf
     zr = zf

     call lininterp_rt(yl,yr,sl,sr,sc,yc)
     call lininterp_rt(zl,zr,sl,sr,sc,zc)


     vl = bdistfun_rt(yl,zl,aratio,q)
     vr = bdistfun_rt(yr,zr,aratio,q)
     vc = bdistfun_rt(yc,zc,aratio,q)

     tol = min(reltol*(dabs(vl)+dabs(vr)),abstol)

     if (diagnose) then
          write(*,*) ""
     end if

     do i=1,nbi
          ! Bisection iteration
          if (vc*vl .lt. 0.0d0) then
               if (diagnose) then
                    write(*,*) "replace right end point"
                    write(*,*) sl,sc,sr
                    write(*,*) vl,vc,vr
               end if
               vr = vc
               sr = sc
               yr = yc
               zr = zc
               sc = (sl + sr)/2.0d0


          else

               if (diagnose) then
                    write(*,*) "replace left end point"
                    write(*,*) sl,sc,sr
                    write(*,*) vl,vc,vr
               end if    
               vl = vc
               sl = sc
               yl = yc
               zl = zc

               sc = (sl + sr)/2.0d0

          end if

          call lininterp_rt(yl,yr,sl,sr,sc,yc)
          call lininterp_rt(zl,zr,sl,sr,sc,zc)
          vc = bdistfun_rt(yc,zc,aratio,q)
          
!          write(*,*) dabs(vc),tol,yc,zc,sc
          if (diagnose) then
               write(*,*) i,sc,vc,tol
          end if
     end do

     if (diagnose) then     
          write(*,*) "end bisection, begin newton"
     end if

     flag = ((dabs(vc) .gt. tol) .and. (nni .le. mnni))
     do while (flag)
          ! Newton iteration.
          sc = sc - bdistfun_rt(yc,zc,aratio,q)/dderiv_rt(y0,yf,z0,zf,aratio,q,sc)

          call lininterp_rt(yl,yr,sl,sr,sc,yc)
          call lininterp_rt(zl,zr,sl,sr,sc,zc)
          vc = bdistfun_rt(yc,zc,aratio,q)

          nni = nni + 1
          if (diagnose) then
               write(*,*) nni,sc,vc,tol
          end if
!          write(*,*) dabs(vc),tol,yc,zc,sc
          flag = ((dabs(vc) .gt. tol) .and. (nni .le. mnni))
     end do

     if (diagnose) then
          write(*,*) "end newton"
     end if

     ! The rootfinder is having some difficulty with the non-convex domains. 
     ! If we get a solution outside of [0,1], hard limit it.
     s = max(sc,0.0d0)
     s = min(s,1.0d0)
     
end subroutine calc_intersection_pt_rt
!
! -------------------------------------------------
!
subroutine lininterp_rt(zl,zr,sl,sr,sc,zc)
! Linear interpolation.
implicit none
     double precision, intent(in)  :: zl,zr,sl,sr,sc
     double precision, intent(out) :: zc

     double precision    :: m

     m = (zr-zl)/(sr-sl)
     zc = zl + m*(sc-sl)
end subroutine lininterp_rt
!
!
!
double precision function dderiv_rt(y0,yf,z0,zf,aratio,q,s)
! Directional derivative for use in Newton's method above.
! Keep in mind y0,yf,z0,zf are parameters here; the 
! derivative is essentially in the direction of the 
! vector from (y0,z0) to (yf,zf), evaluated at the 
! point s.
!
! The formula was calculated in mathematica.
!

implicit none
     double precision, intent(in)  :: y0,yf,z0,zf,aratio,q,s
     double precision              :: x,l

     l = aratio
     
     x = 0.0d0

     x = x + (2*(y0*(-y0 + yf) + q**2*z0*(-z0 + zf) + 2*l**2*q**2* & 
     &       (y0**4 - y0**3*yf + 3*y0*yf*z0**2 + z0**3*(z0 - zf) + 3*y0**2*z0*(-2*z0 + zf)) + &
     &      l**4*(-2*y0**4 + 2*y0**3*yf - y0*(yf + 6*yf*z0**2) + z0*(q**2 - 2*z0**2)*(z0 - zf) + &
     &         y0**2*(1 + 12*z0**2 - 6*z0*zf))))/(-1 + l**2*q**2)

     x = x + (2*s*((y0 - yf)**2 + q**2*(z0 - zf)**2 - &
     &    6*l**2*q**2*(y0**2 - z0*(yf + z0 - zf) - y0*(yf - 2*z0 + zf))* &
     &       (y0**2 + y0*(-yf - 2*z0 + zf) + z0*(yf - z0 + zf)) + &
     &      l**4*(6*y0**4 - 12*y0**3*yf - yf**2*(1 + 6*z0**2) + 2*y0*yf*(1 + 6*z0*(3*z0 - 2*zf)) - &
     &         (q**2 - 6*z0**2)*(z0 - zf)**2 + y0**2*(-1 + 6*yf**2 - 6*(6*z0**2 - 6*z0*zf + zf**2)))))/ &
     & (-1 + l**2*q**2)

     x = x + (12*l**2*(-l + q)*(l + q)*s**2*(y0**4 - 3*y0**3*yf + z0*(-3*yf**2 + (z0 - zf)**2)*(z0 - zf) + &
     &      3*y0**2*(yf**2 - 2*z0**2 + 3*z0*zf - zf**2) - y0*yf*(yf**2 - 3*(3*z0**2 - 4*z0*zf + zf**2))))/ & 
     &      (-1 + l**2*q**2)

     x = x + (4*l**2*(l**2 - q**2)*s**3*(y0**4 - 4*y0**3*yf + yf**4 - 4*y0*yf*(yf**2 - 3*(z0 - zf)**2) - & 
     &      6*yf**2*(z0 - zf)**2 + &
     &      (z0 - zf)**4 + 6*y0**2*(yf + z0 - zf)*(yf - z0 + zf)))/(-1 + l**2*q**2)

     dderiv_rt = x
end function dderiv_rt
