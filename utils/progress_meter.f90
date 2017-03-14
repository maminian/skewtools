subroutine progress_meter(k,n,in_place)
! Displays a simple percentage progress meter, 
! assuming this is placed inside a do loop, k=1,...,n.

use mod_duration_estimator

implicit none
     integer, intent(in) :: k,n
     logical, intent(in) :: in_place

     mde_pttc = predict_completion(mde_ntt,mde_ntc,mde_dts)     

     call mde_pretty_print_time(mde_pttc,mde_pttc_pretty,mde_time_unit)

     if (in_place) then
          ! Use the gfortran functionality to shift the cursor 
          ! to the left, to have a "dynamic" percentage.
          ! Does not work with ifort, which is what the case below handles.
          if (k .eq. 2) then
               write(*,"(A27)",advance="no") "   Time remaining:         "
          end if

          if (mod(k*100,n) .lt. 100) then
               write(*,"(A27)",advance="no") "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
               write(*,"(A4)",advance="no") "\b\b\b\b"
               write(*,"(I3,A1,A3,A16,F5.1,A3)",advance="no") floor(dble(k)/dble(n)*100.0d0),"%",&
                                   "   ","Time remaining: ",mde_pttc_pretty,mde_time_unit
          else
               write(*,"(A8)",advance="no") "\b\b\b\b\b\b\b\b"
               write(*,"(F5.1,A3)",advance="no") mde_pttc_pretty,mde_time_unit
          end if
          

     else
          if (mod(k*100,n) .lt. 100) then
               write(*,"(I3,A1,A3,A16,F5.1,A3)") floor(dble(k)/dble(n)*100.0d0),"%",&
                                   "   ","Time remaining: ",mde_pttc_pretty,mde_time_unit
          end if
     end if

end subroutine progress_meter
