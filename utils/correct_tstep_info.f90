subroutine correct_tstep_info(ntt,nt,target_times,dtmax)
! This is a 'dry run' version of generate_internal_timestepping, which 
! calculates the internal time array size, nt.
!
! ntt is the number of target time points which will be saved to file. 
! Essentially, the target times will grow exponentially, and the internal 
! timestepping will prevent timesteps from exceeding dtmax.
! 
! The resulting nt is an _upper bound_ for the number of 
! internal timesteps needed.

implicit none
     
     integer, intent(in)                               :: ntt
     integer, intent(inout)                            :: nt
     double precision, intent(in)                      :: dtmax
     double precision, dimension(1:ntt), intent(in)    :: target_times

     double precision                   :: dist,tnext,tcurr
     integer                            :: k,tt_idx,idx,i
     

     idx = 0
     do i=2,ntt
          tcurr = target_times(i-1)
          tnext = target_times(i)
          do while (tcurr .lt. tnext)
               idx = idx + 1
               tcurr = tcurr + dtmax
          end do
     end do

     idx = idx + 1

     nt = idx
     
end subroutine correct_tstep_info

