module mod_ductflow
     ! For arrays relevant to calculate the duct flow.


     integer                                           :: ui,uj,nTerms
     double precision, dimension(:), allocatable       :: ya,za
     double precision, dimension(:,:), allocatable     :: u_precomp


     ! Gridsize for the precalculation of the flow.
     !
     ! 256 in both directions corresponds to spatial step ~0.004,
     ! for uniform grid size.
     ! "Need" 10 points in the boundary layer to be fair,
     ! so this can resolve boundary layer fairly up to aratio = 0.04.
     !
     ! With non-uniform mesh (currently being used) this isn't as 
     ! much of a problem, but I haven't attempted any analysis.
     !

     parameter(ui = 256)
     parameter(uj = 256)

     !
     ! Total number of terms to use in the series when 
     ! precalculating u. Don't need much using the single series formulation.
     !

     parameter(nTerms = 256)
     
end module mod_ductflow
