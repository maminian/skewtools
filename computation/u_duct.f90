double precision function u_duct(y,z,nTerms,idxlist,uij_vals)
! Calculate the approximate value of the flow u(y,z), with precalculated
! Fourier coefficients Amn_vals on indices idxlist.
!
! This is the zero-average flow.

implicit none
     double precision                        :: y,z,pi,pisq,half
     integer                                 :: k,i,j,nTerms
     integer, dimension(1:nTerms,1:2)        :: idxlist
     double precision, dimension(1:nTerms)   :: uij_vals
     
     parameter(pi = 4.0d0*datan(1.0d0))
     parameter(pisq = pi**2)
     parameter(half = 0.5d0)
     
     ! Sum over index set.
     
     u_duct = 0.0d0
     do k=1,nTerms
          i = idxlist(k,1)
          j = idxlist(k,2)
          
          u_duct = u_duct + uij_vals(k)*( dcos((i-half)*pi*y)*dcos((j-half)*pi*z) &
                              - (-1)**(i+j)/(pisq*(i-half)*(j-half)) )
     end do
     
     ! Extra factor of two to match with Francesca's work;
     ! Her flow has Laplacian -2 while mine has Laplacian -1.
     !
     
     u_duct = 2.0d0*u_duct
     
end function u_duct
