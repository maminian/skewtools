double precision function u_ellipse(y,z,a,b)
! Calculate the pipe flow velocity.

implicit none
     double precision    :: y,z,a,b
     double precision    :: c,aratio
     
     aratio = a/b
     
     c = 0.5d0/(1.0d0 + aratio**2)
     
     u_ellipse = c*(0.5d0 - (y/a)**2 - (z/b)**2)
     
     ! Factor of two to make the Laplacian -2.
     u_ellipse = u_ellipse*2.0d0
     
end function u_ellipse
