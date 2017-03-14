module mod_triangle_bdry
     ! Defining the triangle boundary.

     integer, parameter  :: nl = 3
     double precision, dimension(nl,3), parameter :: lls=reshape( (/ 1.0d0, 1.0d0, 1.0d0, &
                  -0.5d0, -0.5d0, 1.0d0, &
                  dsqrt(3.0d0)/2.0d0, -dsqrt(3.0d0)/2.0d0, 0.0d0  /) , (/nl,3/) )

end module mod_triangle_bdry
