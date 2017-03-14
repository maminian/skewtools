double precision function u_racetrack(y,z,aratio,q)
! Calculate the racetrack flow velocity.

implicit none
     double precision    :: y,z,aratio,q,l

     l = aratio
     
     u_racetrack = (y**4 - 6*y**2*z**2 + z**4)*l**2*(l**2-q**2)

     ! The q**2*z**2 here is not a typo
     u_racetrack = u_racetrack + (y**2 + q**2*z**2)*(1.0d0-l**4) 

     u_racetrack = u_racetrack*(-1.0d0)/(1.0d0-q**2*l**2)

     ! This is actually the lab-frame flow, but the central 
     ! statistics are calculated, the sample mean is subtracted 
     ! off anyways.
     u_racetrack = u_racetrack + 1.0d0
     
     ! Appropriate factor to make the Laplacian -2.
     u_racetrack = u_racetrack*(1.0d0+q**2)*(1.0d0-l**4)/(1.0d0-q**2*l**2)
     
end function u_racetrack
