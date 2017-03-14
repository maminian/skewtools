! This is a set of functions needed to calculate the exact statistics
! in the duct. The functions take the formulas for the noncentered
! moments \int_{-\infty}^\infty x^j (\int_-1^1 \int_-1^1 T_0(x-ut,y,z) dy dz) dx
!
! via eigenfunction expansion of u, derived in the paper, and returns
! the values as a result of truncating the relevant series.
!

!
! ------------------
!

double precision function avg_phij(i,j)

     implicit none
     integer                       :: i,j
     double precision, parameter   :: pisq = (4.0d0*datan(1.0d0))**2
     double precision, parameter   :: half = 0.5d0
     
     avg_phij = dble((-1)**(i+j)) / (pisq*(i-half)*(j-half))
     
end function avg_phij
!
! ------------------
!
double precision function avg_phij_sq(i,j)

     implicit none
     integer   :: i,j
     
     avg_phij_sq = 0.25d0
     
end function avg_phij_sq
!
! ------------------
!
double precision function compute_u1_noncenter_exact(nTerms,idxlist,uij_vals)
     implicit none
     ! Input
     integer                                 :: nTerms
     integer, dimension(1:nTerms,1:2)        :: idxlist
     double precision, dimension(1:nTerms)   :: uij_vals
     
     ! Internal
     integer   :: k,i,j
     double precision    :: twooverpi
     parameter(twooverpi = 0.5d0/datan(1.0d0))
     
     ! Functions
     double precision avg_phij
     
     compute_u1_noncenter_exact = 0.0d0
     do k=1,nTerms
          i = idxlist(k,1)
          j = idxlist(k,2)
          
          compute_u1_noncenter_exact = compute_u1_noncenter_exact + uij_vals(k)*avg_phij(i,j)
     end do
     

end function compute_u1_noncenter_exact

! ------------------

double precision function compute_u2_noncenter_exact(nTerms,idxlist,uij_vals)
     implicit none
     ! Input
     integer                                 :: nTerms
     integer, dimension(1:nTerms,1:2)        :: idxlist
     double precision, dimension(1:nTerms)   :: uij_vals
     
     ! Internal
     integer   :: k,i,j
     
     ! Functions
     double precision avg_phij_sq
     
     compute_u2_noncenter_exact = 0.0d0
     do k=1,nTerms
          i = idxlist(k,1)
          j = idxlist(k,2)
          
          compute_u2_noncenter_exact = compute_u2_noncenter_exact + (uij_vals(k)**2)*(avg_phij_sq(i,j))
     end do
     

end function compute_u2_noncenter_exact

! ------------------

double precision function compute_u3_noncenter_exact(nTerms,idxlist,uij_vals)
     implicit none
     ! Input
     integer                                                          :: nTerms
     integer, dimension(1:nTerms,1:2)                                 :: idxlist
     double precision, dimension(1:nTerms)                            :: uij_vals
     
     ! Internal
     integer   :: a,b,c,i,j,m,n,p,q
     
     ! Functions
     double precision Alpha_eval
     
     
     compute_u3_noncenter_exact = 0.0d0
     do a=1,nTerms
          i = idxlist(a,1)
          j = idxlist(a,2)
          do b=1,nTerms
               m = idxlist(b,1)
               n = idxlist(b,2)
               do c=1,nTerms
                    p = idxlist(c,1)
                    q = idxlist(c,2)
                    
                    compute_u3_noncenter_exact = compute_u3_noncenter_exact + &
                         uij_vals(a)*uij_vals(b)*uij_vals(c)*Alpha_eval(i,m,p)*Alpha_eval(j,n,q)
                    
               end do
          end do
     end do
     
     
end function compute_u3_noncenter_exact

! ----------------------------------------------------------------------
! BELOW ARE THE SUBROUTINES FOR THE ASYMPTOTIC FLOW.
! --------------

subroutine precompute_A_asymp(A_asymp,nTerms,aratio)

     implicit none
     
     ! Input
     integer             :: nTerms
     double precision    :: aratio
     
     ! Output
     double precision, dimension(1:nTerms)   :: A_asymp
     
     ! Internal
     integer             :: k
     double precision    :: negfouroverpicubed
     
     parameter(negfouroverpicubed = -0.12900613773279796d0)     
     
     do k=1,nTerms
          A_asymp(k) = 0.0d0
     end do
     do k=1,nTerms,2
          A_asymp(k) = negfouroverpicubed / dble(k**3)
     end do

end subroutine precompute_A_asymp

! ------------------

subroutine compute_ubar_asymp(A_asymp,nTerms,ubar_asymp,aratio)
     implicit none
     ! Input
     integer                                 :: nTerms
     double precision, dimension(1:nTerms)     :: A_asymp
     double precision                        :: aratio
     
     ! Output
     double precision    :: ubar_asymp
     
     ! Internal
     integer             :: k
     double precision    :: fouroverpisq, q, exp_q
     
     parameter(fouroverpisq = 0.405284734569351086d0)
     q = -3.141592653589793238d0 / aratio
     exp_q = dexp(q)
     
     ubar_asymp = 0.0d0
     do k=1,nTerms,2
          ubar_asymp = ubar_asymp + aratio * (fouroverpisq) * A_asymp(k)/dble(k**2) * (1.0d0 - exp_q**k)
     end do
     
     ubar_asymp = ubar_asymp + 1.0d0/12.0d0

end subroutine compute_ubar_asymp

! ------------------

subroutine compute_usqu_asymp(A_asymp,nTerms,usqu_asymp,aratio)
     implicit none
     ! Input
     integer   :: nTerms
     double precision, dimension(1:nTerms)     :: A_asymp
     double precision                        :: aratio
     
     ! Output
     double precision    :: usqu_asymp
     
     ! Internal
     integer   :: k
     double precision    :: twooverpi, q, exp_q
     
     parameter(twooverpi = 0.63661977236758134d0)
     
     q = -3.141592653589793238d0 / aratio
     exp_q = dexp(q) 
     
     usqu_asymp = 0.0d0
     do k=1,nTerms,2
          ! First part of sum
          usqu_asymp = usqu_asymp + 2.0d0/q * ( A_asymp(k)**2 / dble(k) * (1.0d0 - exp_q**k) )
          
          ! Second part of sum
          usqu_asymp = usqu_asymp + A_asymp(k)**2 * ( -1.0d0/(2.0d0*dble(k)*q)*(1.0d0 - exp_q**(2*k)) + exp_q**k )
     end do
     
     usqu_asymp = usqu_asymp + 1.0d0/120.0d0

end subroutine compute_usqu_asymp

! ------------------

subroutine compute_ucub_asymp(A_asymp,nTerms,Alpha,alphaMax,ucub_asymp,aratio)
     implicit none
     ! Input
     integer   :: nTerms,alphaMax
     double precision, dimension(1:nTerms)                              :: A_asymp
     double precision, dimension(1:alphaMax,1:alphaMax,1:alphaMax)    :: Alpha
     double precision                                                 :: aratio
     
     ! Output
     double precision    :: ucub_asymp
     
     ! Internal
     integer             :: k,l,m
     double precision    :: twooverpi, q, exp_q, pi
     
     double precision Beta_tilde,Alpha_eval
     
     parameter(twooverpi = 0.63661977236758134d0)
     parameter(pi = 3.141592653589793238d0)
     
     q = - pi / aratio
     exp_q = dexp(q)
     
     
     ucub_asymp = 0.0d0

     ! Breaking full expression into its parts.
     ! You could theoretically could wrap these into each other, but it probably isn't worth it,
     ! and the compiler will probably do some wizardry with the summations regardless.
     
     ! Un
     do k=1,nTerms,2
          ucub_asymp = ucub_asymp + 3.0d0 * A_asymp(k) * &
                         ( -((pi*dble(k))**2 - 12.0d0)) / ((pi*dble(k))**5) * &
                         ( -2.0d0/q*(1.0d0 - exp_q**k)/dble(k) )
     end do

     ! Deux
     do k=1,nTerms,2
          do l=1,nTerms,2
          
               if (k .eq. l) then
                    ucub_asymp = ucub_asymp - 3.0d0 * A_asymp(k)*A_asymp(l) * &
                    ( -(1.0d0)/(8.0d0*(pi*dble(k))**2) - 1.0d0/24.0d0 ) * &
                    ( 2.0d0*exp_q**k - 1.0d0/q * (1.0d0 - exp_q**(2*k)/dble(k)) )
               else
                    ! Differs slightly from expression; keep in mind k,l both odd.
                    ucub_asymp = ucub_asymp - 3.0d0 * A_asymp(k)*A_asymp(l) * &
                              ( dble(k*l)*2.0d0 / (pi*(k-l)*(k+l))**2 ) * &
                              ( -2.0d0/q * ( (exp_q**l - exp_q**k)/dble(k-l) + (1.0d0-exp_q**(k+l))/(k+l) ) )
               end if
          end do
     end do

     ! Trois
     do k=1,nTerms,2
          do l=1,nTerms,2
               do m=1,nTerms,2
                    ucub_asymp = ucub_asymp + A_asymp(k)*A_asymp(l)*A_asymp(m) * &
                              Alpha_eval(k,l,m) * Beta_tilde(k,l,m,aratio)
               end do
          end do
     end do
     
     ucub_asymp = ucub_asymp + 1.0d0/1120.0d0
     
end subroutine compute_ucub_asymp
