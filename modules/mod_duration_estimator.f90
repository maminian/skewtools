module mod_duration_estimator
     ! For estimating how much longer the full program will 
     ! take to complete based on history of timesteps.

     integer                                      :: count_rate  ! Count used in system_clock
     integer                                      :: mde_ntt,mde_ntc
     double precision, allocatable, dimension(:)  :: mde_dts
     double precision                             :: mde_pttc,mde_pttc_pretty
     character(len=3)                             :: mde_time_unit

     integer                                      :: mde_t1,mde_t2

contains
     double precision function predict_completion(mde_ntt,mde_ntc,mde_dts)
     ! Based on the total number of timesteps and number of timesteps 
     ! completed, and a filled history mde_dts(1:mde_ntc), 
     ! predict how much time is remaining assuming time per timestep is 
     ! relatively consistent.

          integer                                 :: mde_ntt,mde_ntc
          double precision, dimension(1:mde_ntt)  :: mde_dts
          
          predict_completion = sum(mde_dts(1:mde_ntc)) * (mde_ntt-mde_ntc)/mde_ntc
          
     end function predict_completion

     subroutine mde_pretty_print_time(pttc,pretty_val,time_unit)
     ! Converts pttc into the smallest "standard" time units so that it is bounded by 100. The 
     ! input pttc is assumed in units of seconds. The time units going out 
     ! are also output ("sec","min","hrs","dys","yrs")

          double precision, intent(in)  :: pttc
          double precision, intent(out) :: pretty_val
          character(len=3)              :: time_unit
          
          pretty_val = pttc

          ! Seconds
          if (pretty_val .lt. 60.0d0) then
               time_unit = "sec"
               go to 10
          end if
          
          ! Minutes
          pretty_val = pretty_val/60.0d0
          if (pretty_val .lt. 60.0d0) then
               time_unit = "min"
               go to 10
          end if

          ! Hours
          pretty_val = pretty_val/60.0d0
          if (pretty_val .lt. 24.0d0) then
               time_unit = "hrs"
               go to 10
          end if

          ! Days
          pretty_val = pretty_val/24.0d0
          if (pretty_val .lt. 365.25d0) then
               time_unit = "day"
               go to 10
          end if

          ! Years
          pretty_val = pretty_val/365.25d0
          time_unit = "yrs"
          
10        continue
     end subroutine mde_pretty_print_time

end module mod_duration_estimator
