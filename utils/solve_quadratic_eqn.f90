subroutine solve_quadratic_eqn(a,b,c,t1,t2)
! Solves the quadratic equation
! a*t**2 + b*t + c .eq. 0.
!
! The two solutions get saved in t1 and t2.
!
! The solutions are assumed real.

implicit none
     double precision, intent(in)  :: a,b,c
     double precision, intent(out) :: t1,t2

     double precision              :: discr
     
     discr = b**2 - 4.0d0*a*c
     
     t1 = (-b - dsqrt(discr))/(2.0d0*a)
     t2 = (-b + dsqrt(discr))/(2.0d0*a)

end subroutine solve_quadratic_eqn
