double precision function u_channel(y,a,b)
! Calculate the channel flow velocity; the channel is in the interval
! a,b. This flow is guaranteed to be integral zero on [a,b] and 
! second derivative (Laplacian) -1.

implicit none
     double precision y,a,b
     
     u_channel = 0.5d0*( (y - a)*(b - y) - 1.0d0/6.0d0*(b-a)**2 )
     
     ! Multiply by a factor of two to match Francesca's.
     
     u_channel = 2.0d0*u_channel
     
end function u_channel
