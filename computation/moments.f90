subroutine moments(n,x,mean,var,skew,kurt)
! Combination function to compute the mean, variance,
! skewness, kurtosis of an array x.
!
! Adapted from Numerical Recipes.
implicit none
     ! In/out
     integer, intent(in)                          :: n
     double precision, dimension(1:n), intent(in) :: x
     double precision, intent(out)                :: mean,var,skew,kurt
     
     ! Internal
     integer                            :: i
     double precision                   :: temp

     ! Handle degenerate cases.
     if (n .eq. 0) then
          mean = 0.0d0
          var = 0.0d0
          skew = 0.0d0
          kurt = 0.0d0
     else if (n .eq. 1) then
          mean = x(1)
          var = 0.0d0
          skew = 0.0d0
          kurt = 0.0d0
     else

     ! Usual algorithm, build the mean, then 
     ! build the centralized statistics based off of that.
          mean = 0.0d0
          do i=1,n
               mean = mean + x(i)
          end do
          mean = mean/n
          
          var = 0.0d0
          skew = 0.0d0
          kurt = 0.0d0
          
          do i=1,n

               temp = x(i) - mean
               var = var + temp**2
               skew = skew + temp**3
               kurt = kurt + temp**4
          end do
          
          var = var/dble(n-1)
          skew = skew/(n*(var**1.5d0))
          kurt = kurt/(n*(var**2)) - 3.0d0
     
     end if
     
     if (var .eq. 0.0d0) then
          skew = 0.0d0
          kurt = 0.0d0
     end if
          
end subroutine moments
