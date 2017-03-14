double precision function bdistfun_rt(y,z,aratio,q)
! Boundary distance function for the racetrack.
! If positive, in the interior, if negative, in 
! the exterior, zero on the boundary.
!
! Essentially the Dirichlet flow solution.

implicit none
     double precision    :: y,z,aratio,q,l

     l = aratio
     
     bdistfun_rt = (y**4 - 6*y**2*z**2 + z**4)*l**2*(l**2-q**2)

     ! The q**2*z**2 here is not a typo
     bdistfun_rt = bdistfun_rt + (y**2 + q**2*z**2)*(1.0d0-l**4) 

     bdistfun_rt = bdistfun_rt*(-1.0d0)/(1.0d0-q**2*l**2)

     bdistfun_rt = bdistfun_rt + 1.0d0
     
end function bdistfun_rt
!
!-------------------------------------------------------
!
subroutine bdistfun_rt_grad(y,z,aratio,q,vec)
! Outward unit normal gradient for the racetrack.
! Partial in y, then partial in z.
!

implicit none
     double precision, intent(in)                 :: y,z,aratio,q
     double precision                             :: l

     double precision, dimension(2), intent(out)  :: vec

     l = aratio
     
     vec(1) = (-2*y*(-1 + 2*l**2*q**2*(y**2 - 3*z**2) + l**4*(1 - 2*y**2 + 6*z**2))) & 
                              &/(-1 + l**2*q**2)

     vec(2) = (2*z*(2*l**4*(-3*y**2 + z**2) + q**2*(1 - l**4 + 6*l**2*y**2 - 2*l**2*z**2))) & 
                              &/(-1 + l**2*q**2)

     call normalize(2,vec)

end subroutine bdistfun_rt_grad

