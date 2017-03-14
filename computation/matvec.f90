subroutine matvec(m,n,A,u,v)
!
! A simple matrix-vector multiplication subroutine.
!
! Inputs:
!
!    m,n: integers
!    A: double precision array, dimension(m,n)
!    u: double precision array, dimension(n)
!
! Outputs: 
!
!    v: double precision array, dimension(m)
!

implicit none
     integer, intent(in)                          :: m,n
     double precision, dimension(m,n), intent(in) :: A
     double precision, dimension(n), intent(in)   :: u
     
     double precision, dimension(m), intent(out)  :: v

     integer                                      :: i,j

     do i=1,m
          v(i) = 0.0d0
          do j=1,n
               v(i) = v(i) + A(i,j)*u(j)
          end do
     end do

end subroutine matvec
