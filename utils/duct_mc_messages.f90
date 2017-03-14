subroutine duct_mc_messages(msg,ny)
! Is this bad practice?
implicit none

     character(len=1024), intent(in)    :: msg
     integer, intent(in)                :: ny
     
     character(len=1024)                :: missing_args,resolution,&
                                             simul_start,simul_done,done
     
     parameter(missing_args = 'missing_args')
     parameter(resolution = 'resolution')
     parameter(simul_start = 'simul_start')
     parameter(simul_done = 'simul_done')
     parameter(done = 'done')
     
     if (msg .eq. missing_args) then
          write(*,*) "You must specify the name of an output file in the third argument. Exiting."
     else if (msg .eq. resolution) then
          write(*,*) ""
          write(*,*) "You need to specify a larger number of points to accurately"
          write(*,*) "resolve the y direction for your aspect ratio. Current rule"
          write(*,*) "of thumb is that you need nGates > 7/sqrt(aratio)."
          write(*,*)
          write(*,*) "Currently, ny=",ny,"."
          write(*,*)
          write(*,*) "Exiting."
          write(*,*) ""
     else if (msg .eq. simul_start) then
          write(*,"(A14)",advance="no") "Simulating... "
          write(*,"(I3,A1)",advance="no") 0,"%"
     else if (msg .eq. simul_done) then
          write(*,"(A8)") "\b\b\bdone."
          write(*,"(A19)",advance="no") "Writing to file... "
     else if (msg .eq. done) then
          write(*,"(A5)") " done."
     end if
     
end subroutine duct_mc_messages


