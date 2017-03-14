subroutine precalculate_Alpha(Alpha,alphaMax)

     implicit none
     
     ! Input variables
     integer :: alphaMax
     
     ! Input/output variables
     double precision, dimension(1:alphaMax,1:alphaMax,1:alphaMax) :: Alpha
     
     ! Internal variables
     integer             :: i,m,p
     
     ! Functions
     double precision Alpha_eval
     
     ! This loop could possibly be optimized further, but it will probably never
     ! be a bottleneck, so it's a low priority.
     do i=1,alphaMax
          do m=1,alphaMax
               do p=1,alphaMax
                    Alpha(i,m,p) = Alpha_eval(i,m,p)
               end do
          end do
     end do
     
end subroutine precalculate_Alpha
