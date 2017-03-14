subroutine my_normal_rng(na,array,mean,variance)
! Populates array size n with normally distributed
! random variables with given mean and variance.
!
! Done by generating pairs of uniform random [0,1]
! with Fortran's built-in function, then doing
! the Box-Muller transform to get iid normal vars.
!
! The Mersenne Twister is assumed already initialized/seeded
! in a parent function.

use mtmod ! Mersenne Twister module

implicit none

     integer, intent(in)                               :: na
     double precision, dimension(1:na), intent(inout)  :: array
     double precision, intent(in)                      :: mean,variance
     
     double precision                                  :: stdev,modulus,modsq,mult,twopi
     double precision, dimension(1:2)                  :: pair
     logical                                           :: nisodd,good
     integer                                           :: krng,nac
     
     parameter( twopi = 6.283185307179586d0 )
     
     stdev = dsqrt(variance)

     nac=na
     
     nisodd = (mod(nac,2) .eq. 1)
     if (nisodd) then
          nac=nac-1
     end if


     ! Testing "Numerical Recipes" version, avoiding calculating
     ! cosine and sine with a variation on Box-Muller.


     do krng=1,nac,2
     
          ! Grab pairs of points until you get a pair 
          ! that lies in the unit ball.
          good = .false.
          do while (.not. good)
!               call random_number(pair)
               pair(1) = grnd()
               pair(2) = grnd()
               pair = 2.0d0*pair - 1.0d0
               modsq = pair(1)**2 + pair(2)**2
               
               ! Logical evaluation
               good = (modsq .lt. 1.0d0)
          end do
          
          ! Then, apply the formula.
          mult = dsqrt(-2.0d0*dlog(modsq)/modsq)
          array(krng) = mean + stdev*mult*pair(1)
          array(krng+1) = mean + stdev*mult*pair(2)
          
     end do
          
          
     ! Handle the last point with a normal Box-Muller if necessary.
     if (nisodd) then
          nac = nac+1
!          call random_number(pair)
          pair(1) = grnd()
          pair(2) = grnd()
          modulus = dsqrt(-2.0d0 * dlog(pair(1)))
          array(nac) = mean + stdev*modulus*dcos(twopi*pair(2))
     end if
     
end subroutine my_normal_rng
