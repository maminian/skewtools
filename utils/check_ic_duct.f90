logical function check_ic_duct(nTot,Y,Z,a,b)
! Checks that the initial data is contained within the cross section.
implicit none

     integer, intent(in)                               :: nTot
     double precision, dimension(nTot), intent(in)     :: Y,Z
     double precision, intent(in)                      :: a,b

     integer                                           :: i
     double precision                                  :: ymin,ymax,zmin,zmax
     

     ymin = minval(Y)
     ymax = maxval(Y)
     zmin = minval(Z)
     zmax = maxval(Z)


     if ((ymin .lt. -a) .or. (ymax .gt. a) .or. (zmin .lt. -b) .or. (zmax .gt. b)) then
          check_ic_duct = .false.
          write(*,*) ymin,ymax,a
          write(*,*) zmin,zmax,b
     else
          check_ic_duct = .true.
     end if
          
end function check_ic_duct
