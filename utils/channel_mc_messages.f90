subroutine channel_mc_messages(msg)
! Is this bad practice?
implicit none

     character(len=1024), intent(in)    :: msg
     
     character(len=1024)                :: missing_args,simul_start,simul_done,done
     
     parameter(missing_args = 'missing_args')
     parameter(simul_start = 'simul_start')
     parameter(simul_done = 'simul_done')
     parameter(done = 'done')
     
     if (msg .eq. missing_args) then
          write(*,*) "You must specify the name of an output file in the third argument. Exiting."
     else if (msg .eq. simul_start) then
          write(*,"(A14)",advance="no") "Simulating... "
          write(*,"(I3,A1)",advance="no") 0,"%"
     else if (msg .eq. simul_done) then
          write(*,"(A8)") "\b\b\b done."
          write(*,"(A19)",advance="no") "Writing to file... "
     else if (msg .eq. done) then
          write(*,*) " done."
     end if
     
end subroutine channel_mc_messages
