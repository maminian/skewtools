program rng_tester
! Takes n, mean, variance as arguments, samples a normal
! distribution with those quantities, and spits back out
! the sample mean, variance, skewness.

implicit none
     
     integer                                      :: n
     double precision                             :: inmean,invar
     double precision, dimension(:), allocatable  :: rng
     character(len=1024)                          :: inarg
     double precision                             :: mean,var,skew,kurt

     

     call get_command_argument(1,inarg)
     read(inarg,*) n
     call get_command_argument(2,inarg)
     read(inarg,*) inmean
     call get_command_argument(3,inarg)
     read(inarg,*) invar

!     n = 10**7
!     inmean = -1.2d0
!     invar = 4.3d0

     write(*,*) ""
     write(*,*) "Input parameters:"
     write(*,*) (/ dble(n), inmean,invar /)
     write(*,*) ""

     allocate( rng(1:n) )
     
     ! Generate numbers from my rng and take a look!
     call my_normal_rng(n,rng,inmean,invar)
     
     ! Get statistics
     call moments(n,rng,mean,var,skew,kurt)
     
     write(*,*) "Sample statistics (mean, variance, skew, kurtosis):"
     write(*,*) mean
     write(*,*) var
     write(*,*) skew
     write(*,*) kurt

     deallocate(rng)
     
end program rng_tester
