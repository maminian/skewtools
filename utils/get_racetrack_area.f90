subroutine get_racetrack_area(q,aratio,area)
! Calculates the area in the racetrack with the given parameters.
! Simple Monte Carlo method with rejection.

use mtmod

implicit none
     double precision, intent(in)       :: q,aratio
     double precision, intent(out)      :: area

     integer, parameter                 :: ntot = 10**9
     integer                            :: nin,i
     double precision                   :: area_rec,yl,yr,zl,zr,my,mz,y,z


     double precision bdistfun_rt
     ! ---------------------
     
     yl = -1.3d0
     yr = 1.3d0

     zl = yl/aratio
     zr = yr/aratio

     my = yr-yl
     mz = zr-zl

     area_rec = my*mz

     do i=1,ntot
          ! Uniform random
          y = yl + my*grnd()
          z = zl + mz*grnd()
          
          if (bdistfun_rt(y,z,aratio,q) .ge. 0.0d0) then
               nin = nin + 1
          end if
          
     end do

     area = area_rec*dble(nin)/ntot

end subroutine get_racetrack_area
