! Scripts for basic operations with vectors and pairs of vectors; 
! dot products, lengths,
! applying orthogonal projections, projections orthogonal,
! and reflections. These all MODIFY the first input, so be careful.

double precision function dot_prod(n,u,v)
! Why am I not using BLAS or similar for this? Laziness.
     implicit none
     integer                            :: n,i
     double precision, dimension(n)     :: u,v
     
     dot_prod=0.0d0
     do i=1,n
          dot_prod = dot_prod + u(i)*v(i)
     end do
     
end function dot_prod
!
!
!
double precision function norm(n,u)
     implicit none
     integer                            :: n
     double precision, dimension(n)     :: u
     
     double precision dot_prod

     norm = dsqrt(dot_prod(n,u,u))
end function norm
!
! ----------------------------------------
!
subroutine normalize(n,u)
! Normalizes 2d vector u.
     implicit none
     integer                            :: n
     double precision, dimension(n)     :: u
     double precision                   :: s
     
     double precision norm
     
     s = norm(n,u)
     u = u/s
end subroutine normalize
!
! ----------------------------------------
!
subroutine orth_proj(n,v,u)
! Projects v onto u (u not necessarily unit).
! v is changed on output.
     implicit none
     integer, intent(in)                               :: n
     double precision, dimension(n), intent(in)        :: u
     double precision, dimension(n), intent(inout)     :: v
     
     double precision dot_prod
     
     v = dot_prod(n,u,v)/dot_prod(n,u,u)*u

end subroutine orth_proj
!
! ----------------------------------------
!
subroutine proj_orth(n,v,u)
! The projection of v orthogonal to u,
! u not necessarily unit.
     implicit none
     integer, intent(in)                :: n
     double precision, dimension(n), intent(in)        :: u
     double precision, dimension(n), intent(inout)     :: v
     double precision, dimension(n)                    :: temp
     
     temp = v
     call orth_proj(n,temp,u)
     v = v - temp
end subroutine proj_orth
!
! ----------------------------------------
!
subroutine reflect(n,v,u)
! Reflects v across hyperplane defined by vector u.
! If u is the normal to a surface, it should be the
! _OUTWARD NORMAL_. u does not need to be unit.
     implicit none
     integer, intent(in)                               :: n
     double precision, dimension(n), intent(in)        :: u
     double precision, dimension(n), intent(inout)     :: v
     double precision, dimension(n)                    :: temp
     
     temp = v
     call orth_proj(n,temp,u)
     v = v - 2.0d0*temp
end subroutine reflect            

