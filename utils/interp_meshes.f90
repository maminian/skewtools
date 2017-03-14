subroutine padded_cheb_nodes(n,x,x0,xf)
! Fills array x, length n, with the n-2 Chebyshev nodes
! on the interval (x0,xf), with the first and
! last nodes as x0,xf padded on.
implicit none

     ! In
     integer                            :: n,i
     double precision                   :: x0,xf
     
     ! In/out
     double precision, dimension(n)     :: x
     
     ! Internal
     double precision                   :: pi
     
     parameter(pi = 4.0d0*datan(1.0d0))
     
     x(1) = x0
     do i=2,n-1
          x(i) = 0.5d0*(x0+xf) + 0.5d0*(xf-x0)*dcos(dble(2*(n-i)-1)/dble(2*(n-2))*pi)
     end do
     x(n) = xf


end subroutine padded_cheb_nodes
!
! -----------------------------------------
!
subroutine uniform_nodes(n,x,x0,xf)
! Fills array x, length n, with the uniformly distributed nodes
! on the interval (x0,xf).
implicit none

     ! In
     integer                            :: n,i
     double precision                   :: x0,xf
     
     ! In/out
     double precision, dimension(n)     :: x
     
     ! Internal
     double precision                   :: dx
     
     ! Uniform nodes.
     dx = (xf-x0)/(n-1)
     
     do i=1,n
          x(i) = x0 + (i-1)*dx
     end do
     
end subroutine uniform_nodes
