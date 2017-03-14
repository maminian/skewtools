double precision function linear_interp_2d(um,un,u,x,y,x0,y0)
! Bilinear interpolation in 2d, with the accompanying index-locator
! function.
!
! Given an m-by-n array of "exact" u values, rectangular domain
! defined on a grid with x and y values, and interpolation point
! x0, y0, returns the approximate u(x0,y0) value using a 4-point bilinear 
! interpolation.

implicit none
     ! Input vars
     integer                            :: um,un
     double precision, dimension(um,un) :: u
     double precision, dimension(um)    :: x
     double precision, dimension(un)    :: y
     double precision                   :: x0,y0
     
     ! Internal vars
     double precision                   :: u1,u2,u3,u4,p,q
     integer                            :: j,k
     double precision                   :: one
     parameter(one = 1.0d0)
     
     ! Functions
     integer locate,locate2
     

     ! If on a uniform grid, should use locate2 instead to reduce
     ! the lookup cost.
     if (.true.) then
          j = locate2(um,x(1),x(um),x0)
          k = locate2(un,y(1),y(un),y0)
     else
          j = locate(um,x,x0)
          k = locate(un,y,y0)
     end if

     u1 = u(j,k)
     u2 = u(j+1,k)
     u3 = u(j+1,k+1)
     u4 = u(j,k+1)
     

     p = ( x0 - x(j) )/( x(j+1) - x(j) )
     q = ( y0 - y(k) )/( y(k+1) - y(k) )
     
     linear_interp_2d = (one-p)*(one-q)*u1 + p*(one-q)*u2 + &
                    p*q*u3 + (one-p)*q*u4
     

end function linear_interp_2d
!
! ----------------------------------------------------
!
integer function locate(n,x,x0)
! Given 1d double precision array x dimension n,
! monotonically increasing or decreasing,
! and double precision x0, locate the index
! j which satisfies x(j) <= x0 <= x(j+1).
!
! This is not idiot-proof; if x0 lies outside
! the domain of x then you're SOL.
!
! Adapted from "Numerical Recipes." The basic idea
! behaves like bisection; the endpoints of the array
! serve as the positive/negative bounds of the
! "function" x-x0, then the bounds are iteratively
! refined until we reach a bound of (/j,j+1/).
! 
! Should converge in log_2(n) steps.

implicit none
     ! Input vars
     integer                            :: n
     double precision, dimension(n)     :: x
     double precision                   :: x0

     ! Internal
     integer                            :: jl,ju,jm
     logical                            :: xisincr
     
     ! Upper/lower limits on the containing index
     jl = 0
     ju = n+1
     
     xisincr = ( x(n) .gt. x(1) )
     
     do while (ju-jl .gt. 1)
          
          jm = (ju+jl)/2      ! Compute a midpoint of idx bound
          
          ! Choose the next bound depending on whether
          ! the array x is monotone increasing or decreasing.
          if ( xisincr .eqv. ( x0 .gt. x(jm) ) ) then
               jl = jm
          else
               ju = jm
          end if
          
     end do
          
     ! Return lower index bound.
     ! This will only appear if there is a double reflection
     ! that would be necessary in the code with MC, where
     ! jl=0 since x0 < min(x).
     !
     ! Should not occur with proper simulation, though; 
     ! dt should be chosen small enough relative to Peclet
     ! that probability of a Brownian motion that large is
     ! vanishingly small.
     locate = max(jl,1)
     
end function locate
!
! -----------------------
!
integer function locate2(n,xl,xr,x0)
! I THINK WE CAN DO BETTEH!
! https://www.youtube.com/watch?v=NTpptLoUEk8
!
! Assuming the array x is a uniform mesh, we can get the 
! precise index with some modular arithmetic. xl and xr are 
! the lower and upper bounds of the array.
!
! Also assumes arrays start counting at 1.
implicit none

     ! Input
     double precision    :: xl,xr,x0
     integer             :: n
     ! Internal
     double precision    :: h

     h = (xr-xl)/(n-2)

     ! Add one because indexing starts at 1.
     locate2 = floor( ((x0-xl) - dmod(x0-xl,h)) / h ) + 1

end function locate2
