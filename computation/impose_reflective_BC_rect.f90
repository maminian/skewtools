subroutine impose_reflective_BC_rect(z,lower,upper)
! Imposes reflective boundary
! conditions for the Monte Carlo simulation on channel
! geometry (or duct, if applied in each direction) 
! 
! For now, this assumes there is no double reflecting.
! This relies on a small enough dt that the likelihood
! is outlandishly small.

implicit none
     double precision z,lower,upper,residue
     
     do while ((z .gt. upper) .or. (z .lt. lower))
          if (z .gt. upper) then
          
               residue = z-upper
               z = upper - residue
               
          else if (z .lt. lower) then
               
               residue = lower-z
               z = lower + residue
               
          end if
     end do

end subroutine impose_reflective_BC_rect

