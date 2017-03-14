module mod_active_particles
     ! Contains definitions for all the parameters, 
     ! in addition to some hard-coded values.
     
     integer, dimension(:), allocatable      :: aset
     integer, dimension(:), allocatable      :: asettemp
     integer                                 :: nactive
     double precision, parameter             :: pexit = 0.0d0 ! May get moved later.
     
     contains 
	 
     subroutine get_subset_2d(nTot,X,Y,aset,nactive,Xtemp,Ytemp)
          ! Given full array X,Y, and the active set, generate the list of active particles 
          ! for the purposes of collecting statistics.
          integer, intent(in)                               :: nTot,nactive
          double precision, dimension(nTot), intent(in)     :: X,Y
          double precision, dimension(nTot), intent(out)    :: Xtemp,Ytemp
          integer, dimension(nTot), intent(in)              :: aset
          
          integer                                           :: i
          
          do i=1,nactive
               Xtemp(i) = X(aset(i))
               Ytemp(i) = Y(aset(i))
          end do
     end subroutine get_subset_2d
end module mod_active_particles
