double precision function Alpha_eval(i,m,p)

     implicit none
     
     ! Input variables
     integer :: i,m,p
     
     ! Internal variables
     double precision, parameter   :: half = 0.5d0
     double precision, parameter   :: twooverpi = half/datan(1.0d0)

     Alpha_eval = twooverpi*(i-half)*(m-half)*(p-half)*(-1)**(i+m+p) &
                                        /dble( (i-m-p+half)*(i+m-p-half)*(i-m+p-half)*(i+m+p-3*half) )
end function Alpha_eval
