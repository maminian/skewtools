double precision function u_duct_ss(y,z,nTerms,aratio)
! Calculate the approximate value of the flow u(y,z), 
! using a single-series solution
! whose Laplacian is guaranteed -2 for any number of terms, 
! but boundary conditions at the far walls are only met approximately.
!
! However, the degree of this is relatively insignificant, even 
! for a reasonable number of terms, and there is 
! great added benefit when later moving to calculate the 
! moments of the flow.
!

implicit none
     double precision, intent(in)            :: y,z,aratio
     integer, intent(in)                     :: nTerms

     integer                                 :: k
     double precision                        :: q,yterm,zterm,betak
     double precision                        :: pi,absz
     
     parameter(pi = 4.0d0*datan(1.0d0))
     
     u_duct_ss = 0.0d0

     absz = dabs(z)

     do k = 1,nTerms

          q = (k-0.5d0)*pi

          ! The original equation needs to be shuffled around for it 
          ! to be numerically stable; cosh(qz)/cosh(q/aratio) is apparently ill-behaved.
          yterm = 4.0d0*(-1)**k/(q**3)*dcos(q*y)
          zterm = (dexp(-q*(1.0d0/aratio + z)) + dexp(q*(-1.0d0/aratio + z)) )/(1.0d0 + dexp(-2.0d0*q/aratio))
          
          betak = -4.0d0*aratio*dtanh(q/aratio)/(q**5)
          
          u_duct_ss = u_duct_ss + yterm*zterm - betak

     end do

     
     u_duct_ss = u_duct_ss + ( 1.0d0/3.0d0 - y**2 )

!(dabs(y) .eq. 1.0d0) .or.
!     write(*,*) y,z,u_duct_ss
!     if ( (dabs(z) .eq. 1.0d0/aratio)) then
!          write(*,*) y,z,u_duct_ss
!     end if

end function u_duct_ss
