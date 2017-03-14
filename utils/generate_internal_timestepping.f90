subroutine generate_internal_timestepping(ntt,nt,target_times,t_hist,dtmax)
!
! Generate the array of internal time values based on dtmax and target_times.
! 

implicit none
     
     integer, intent(in)                               :: ntt
     integer, intent(inout)                            :: nt
     double precision, dimension(1:ntt), intent(in)    :: target_times
     double precision, dimension(1:nt), intent(inout)  :: t_hist
     double precision, intent(in)                      :: dtmax
     
     integer                                           :: idx,i
     double precision                                  :: tcurr,tnext

     idx = 0
     do i=2,ntt
          tcurr = target_times(i-1)
          tnext = target_times(i)
          do while (tcurr .lt. tnext)
               idx = idx + 1
               t_hist(idx) = tcurr
               tcurr = tcurr + dtmax
          end do
     end do

     idx = idx + 1
     t_hist(idx) = target_times(ntt)
     
end subroutine generate_internal_timestepping
