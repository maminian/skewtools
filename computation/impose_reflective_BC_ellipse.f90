subroutine impose_reflective_BC_ellipse(yf,zf,y0,z0,a,b,maxrefl)
! Impose reflective boundaries for the ellipse.
! (yf,zf) is the position after taking a step (corrected on output),
! (y0,z0) is the previous position.

implicit none
     double precision, intent(in)       :: y0,z0,a,b
     double precision, intent(inout)    :: yf,zf
     integer, intent(in)                :: maxrefl
     double precision                   :: distsq
     double precision, dimension(2)     :: v1,nhat
     
     integer                            :: nrefl,mmr
     logical                            :: goodsol
     double precision                   :: s,yc,zc,yold,zold

     
     ! Idiot-proofing
     mmr = max(maxrefl,4)


     distsq = (yf/a)**2 + (zf/b)**2
     nrefl = 0

     yold = y0
     zold = z0

     do while ( (distsq .gt. 1.0d0) .and. (nrefl .lt. mmr) )
          ! Find the point 
          call ell_refl_ssols(yold,yf,zold,zf,a,b,s,goodsol)

          if (goodsol) then
!               write(*,*) "oooooop"
               ! Find the coordinates of intersection.
               yc = yold + s*(yf-yold)
               zc = zold + s*(zf-zold)

               ! Construct the vector going out of the domain.
               v1(1) = yf-yc
               v1(2) = zf-zc
               
               ! Construct the outward normal vector.
               nhat(1) = 2.0d0*yc/(a**2)
               nhat(2) = 2.0d0*zc/(b**2)
               call normalize(2,nhat)

               call reflect(2,v1,nhat)

               ! Get new point an update the old point if it's needed.
               yold = yc
               zold = zc

               yf = yc + v1(1)
               zf = zc + v1(2)

          else
               ! Send the thing back to where it started. Ugly, but it works.
!               write(*,*) "moo"
               yf = y0
               zf = z0
          end if

          ! Update the distance of the new point and the number of reflections.
          distsq = (yf/a)**2 + (zf/b)**2
          nrefl = nrefl + 1
     end do

     ! Last line of defense.
     if ((nrefl .eq. mmr) .and. (distsq .gt. 1.0d0)) then
          yf = yc
          zf = zc
     end if

     
end subroutine impose_reflective_BC_ellipse

!
! -----------------------------
!

subroutine ell_refl_ssols(y0,yf,z0,zf,a,b,sol,flag)
implicit none
     double precision, intent(in)  :: y0,yf,z0,zf,a,b
     double precision, intent(out) :: sol
     logical, intent(out)          :: flag

     double precision              :: discrim,denom,sol1,sol2,numterm1,numterm2

     double precision ell_refl_discrim

     discrim = ell_refl_discrim(y0,yf,z0,zf,a,b)

     if (discrim .lt. 0) then
!          write(*,*) "case 0!"
          sol1 = 0.0d0
          sol2 = 0.0d0
          flag = .false.
     else
          denom = b**2*(yf-y0)**2 + a**2*(zf-z0)**2
          numterm1 = b**2*y0*(y0-yf) + a**2*z0*(z0-zf)
          numterm2 = a*b*dsqrt(discrim)

          sol1 = (numterm1 + numterm2)/denom
          sol2 = (numterm1 - numterm2)/denom
          
          flag = .false.

          if ((sol1 .ge. 0.0d0) .and. (sol1 .le. 1.0d0)) then
!               write(*,*) "case 1!"
               sol = sol1
               flag = .true.
          end if

          if ((sol2 .ge. 0.0d0) .and. (sol2 .le. 1.0d0)) then
!               write(*,*) "case 2!"
               if (flag) then
                    ! Both solutions are in [0,1]! 
                    ! This seems to happen when the point is on the 
                    ! boundary. Take the larger time.
!                    write(*,*) "case 2a!"
                    sol = max(sol1,sol2)
                    flag = .true.
               else
                    ! All is well, only one solution.
                    sol = sol2
                    flag = .true.
!                    write(*,*) "case 2b!"
               end if
          end if

     end if

!     write(*,*) sol1,sol2,sol,flag

end subroutine ell_refl_ssols

!
! -----------------------------
!

double precision function ell_refl_discrim(y0,yf,z0,zf,a,b)

implicit none
     double precision, intent(in)  :: y0,yf,z0,zf,a,b
     
     ell_refl_discrim= b**2*(yf-y0)**2 + a**2*(zf-z0)**2 - (yf*z0 - y0*zf)**2

end function ell_refl_discrim


