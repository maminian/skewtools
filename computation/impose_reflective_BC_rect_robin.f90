subroutine impose_reflective_BC_rect_robin(z,lower,upper,absorbed)
! Imposes reflective boundary
! conditions for the Monte Carlo simulation on channel
! geometry (or duct, if applied in each direction) 
! 
! For now, this assumes there is no double reflecting.
! This relies on a small enough dt that the likelihood
! is outlandishly small.

use mtmod
use mod_active_particles

implicit none
     double precision z,lower,upper,residue,v
     logical absorbed

     absorbed = .false.
     do while ( ((z .gt. upper) .or. (z .lt. lower)) .and. (.not. absorbed))
          ! Check if we need to actually impose the reflection,
          ! or just exit. For reasons, we just return a flag
          ! if this is the case.
          v = grnd()
          absorbed = (v .lt. pexit)
          
          if (.not. absorbed) then
               if (z .gt. upper) then
                    
                    residue = z-upper
                    z = upper - residue
                    
               else if (z .lt. lower) then
                    
                    residue = lower-z
                    z = lower + residue
                    
               end if
          end if
     end do

end subroutine impose_reflective_BC_rect_robin

