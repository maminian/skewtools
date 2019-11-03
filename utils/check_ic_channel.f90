logical function check_ic_channel(nTot,Y,a)
! Checks that the initial data is contained within the cross section.
implicit none

     integer, intent(in)                               :: nTot
     double precision, dimension(nTot), intent(in)     :: Y
     double precision, intent(in)                      :: a

     integer                                           :: i
     double precision                                  :: ymin,ymax
     

     ymin = minval(Y)
     ymax = maxval(Y)

     if ((ymin .lt. -a) .or. (ymax .gt. a)) then
          check_ic_channel = .false.
     else
          check_ic_channel = .true.
     end if
          
end function check_ic_channel
