program largepe_root
! Script to calculate the time of the first extremum in 
! the duct, if it exists, in the limit of large Pe 
! (though in principle any Pe can be input to this code.)
! 
! Input arguments are necessary; after compiling, run as
!
! ./largepe_root nTerms amin amax na Pe
!
! Where na aspect ratios are sampled from amin to amax, 
! outputting the time to the extremum, for a given Pe.

implicit none
     
     character(len=1024)           :: inarg, cchoose_mode, ofmt
     integer                       :: nTerms,na,i
     double precision              :: aratio,time,amin,amax,da,Pe

     double precision              :: u1nc,u2nc,u3nc,u2c,u3c,c0,c2,lapu
     
     ! Arrays to be allocated.
     double precision, dimension(:), allocatable  :: aratios, uij_vals, times
     integer, dimension(:,:), allocatable         :: idx_list

     ! functions
     double precision compute_u1_noncenter_exact, compute_u2_noncenter_exact, compute_u3_noncenter_exact

     ! parameters
     parameter(cchoose_mode="regular")
     parameter(lapu = -2.0d0)

     ! Read the input -- number of terms to use in sums
     call get_command_argument(1,inarg)
     read(inarg,*) nTerms

     call get_command_argument(2,inarg)
     read(inarg,*) amin

     call get_command_argument(3,inarg)
     read(inarg,*) amax

     call get_command_argument(4,inarg)
     read(inarg,*) na

     call get_command_argument(5,inarg)
     read(inarg,*) Pe

     if (na .gt. 1) then
          da = (amax-amin)/(na-1)
     else
          da = 0.0d0
     end if
     
     ! Allocate memory to arrays.
     allocate( aratios(na), times(na), uij_vals(nTerms) )
     allocate( idx_list(nTerms,2) )


     ! Write a header. 
     ofmt = "(A5,A12,A5,A13)"

     write(*,"(A5)") "     "
     write(*,ofmt) "     ", "Aspect ratio","     ","Extremum time"
     write(*,ofmt) "     ", "------------","     ","------------------"
     write(*,*) ""

     ofmt = "(A5,ES10.3,A7,ES10.3)"

     ! Calculate the moments over a span of aspect ratios.
     do i=1,na
          aratios(i) = amin + (i-1)*da
          
          ! Get the coefficients uij and index set {(i_k,j_k)}.
          call coeff_chooser(nTerms,aratios(i),idx_list,uij_vals,cchoose_mode)
     

          u1nc = compute_u1_noncenter_exact(nTerms,idx_list,uij_vals)
          u2nc = compute_u2_noncenter_exact(nTerms,idx_list,uij_vals)
          u3nc = compute_u3_noncenter_exact(nTerms,idx_list,uij_vals)

          u2c = u2nc - u1nc**2
          u3c = u3nc - 3.0d0*u1nc*u2nc + 2.0d0*u1nc**3

          c0 = 36.0d0*u3c
          c2 = 6.0d0*lapu * ( u3c*u1nc - u2c*(u2nc - 2.0d0*u1nc**2) )
     
          times(i) = dsqrt(c0/c2)/Pe 

          write(*,ofmt) "     ", aratios(i), "          ", times(i)
     end do

     write(*,*) ""

     deallocate( aratios, times, uij_vals, idx_list )
     
end program largepe_root
