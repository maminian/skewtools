double precision function u_dummy(y,z)
! Dummy flow to be used when the flow is not important.

implicit none
     double precision    :: y,z

     u_dummy = 0.0d0

end function u_dummy
