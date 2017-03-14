double precision function skew_Calc(aratio,nTerms,idxlist,uij_vals,Alpha,alphaMax,flow_type)
! Parent function for calculating the formula for the geometric skewness,
! with a given rectangular aspect ratio aratio, and precomputed Alpha.
!
! NOTE that we do not do a "rectangular sum" over m,n;
! instead, we do a sum over an index set which maximizes the magnitudes of
! the Fourier coefficients for a fixed number of terms.
!
! The hope is to improve the rate of convergence for small aratio.

     implicit none
     
     ! Input variables
     double precision                                                 :: aratio
     integer                                                          :: nTerms,alphaMax
     double precision, dimension(1:alphaMax,1:alphaMax,1:alphaMax)    :: Alpha
     integer, dimension(1:nTerms,1:2)                                 :: idxlist
     double precision, dimension(1:nTerms)                            :: uij_vals
     character(len=1024)                                              :: flow_type
     

     ! Internal variables
     double precision                   :: u1_nc_exact,u2_nc_exact,u3_nc_exact
     double precision                   :: u2_c_exact, u3_c_exact
     
     double precision                   :: ubar_asymp,usqu_asymp,ucub_asymp
     double precision                   :: usqu_zavg_asymp, ucub_zavg_asymp
     
     double precision, dimension(1:nTerms)                            :: A_asymp
     
     ! Aspect ratio cutoff when using the hybrid method.
     double precision, parameter   :: hybrid_cutoff = 0.2d0
     
     ! Functions
     double precision compute_u1_noncenter_exact,&
                         compute_u2_noncenter_exact,&
                         compute_u3_noncenter_exact
     
     if (trim(flow_type) .eq. "exact") then

          u1_nc_exact = compute_u1_noncenter_exact(nTerms,idxlist,uij_vals)
          u2_nc_exact = compute_u2_noncenter_exact(nTerms,idxlist,uij_vals)
          u3_nc_exact = compute_u3_noncenter_exact(nTerms,idxlist,uij_vals)
          
          ! From the three terms above calculate the centered moments.
          u3_c_exact = u3_nc_exact - 3*u1_nc_exact*u2_nc_exact + 2*u1_nc_exact**3
          u2_c_exact = u2_nc_exact - u1_nc_exact**2
          
          skew_Calc = u3_c_exact/(u2_c_exact**(1.5d0))
     
     else if (trim(flow_type) .eq. "asymp") then

          ! For the asymptotic formula, numerical performance is a non-issue,
          ! so we keep the old implementation.
          !
          ! In any case, the summmand is monotone decreasing, so there is
          ! nothing to optimize.
          call precompute_A_asymp(A_asymp,nTerms,aratio)
          
          call compute_ubar_asymp(A_asymp,nTerms,ubar_asymp,aratio)
          call compute_usqu_asymp(A_asymp,nTerms,usqu_asymp,aratio)
          call compute_ucub_asymp(A_asymp,nTerms,Alpha,alphaMax,ucub_asymp,aratio)
          
          ucub_zavg_asymp = ucub_asymp - 3.0d0*ubar_asymp*usqu_asymp + 2.0d0*ubar_asymp**3
          usqu_zavg_asymp = usqu_asymp - ubar_asymp**2
  
          skew_Calc = ucub_zavg_asymp/(usqu_zavg_asymp**(1.5d0))
     
     else if (trim(flow_type) .eq. "hybrid") then
          ! Combine the two. Set an (arbitrary) threshold of aratio=0.2;
          ! greater than that, use the exact, and less, use the asymptotic result.

          if (aratio .lt. hybrid_cutoff) then

               call precompute_A_asymp(A_asymp,nTerms,aratio)
          
               call compute_ubar_asymp(A_asymp,nTerms,ubar_asymp,aratio)
               call compute_usqu_asymp(A_asymp,nTerms,usqu_asymp,aratio)
               call compute_ucub_asymp(A_asymp,nTerms,Alpha,alphaMax,ucub_asymp,aratio)

               
               ucub_zavg_asymp = ucub_asymp - 3.0d0*ubar_asymp*usqu_asymp + 2.0d0*ubar_asymp**3
               usqu_zavg_asymp = usqu_asymp - ubar_asymp**2

               skew_Calc = ucub_zavg_asymp/(usqu_zavg_asymp**(1.5d0))
               
          else


               u1_nc_exact = compute_u1_noncenter_exact(nTerms,idxlist,uij_vals)
               u2_nc_exact = compute_u2_noncenter_exact(nTerms,idxlist,uij_vals)
               u3_nc_exact = compute_u3_noncenter_exact(nTerms,idxlist,uij_vals)
               
               ! From the three terms above calculate the centered moments.
               u3_c_exact = u3_nc_exact - 3*u1_nc_exact*u2_nc_exact + 2*u1_nc_exact**3
               u2_c_exact = u2_nc_exact - u1_nc_exact**2
               
               skew_Calc = u3_c_exact/(u2_c_exact**(1.5d0))
     
          end if
     
     else
          skew_Calc = 0.0d0
          write(*,*) "Unrecognized flow calculation type - use one of ''exact'', ''asymp'', or ''hybrid''."
     end if

end function skew_Calc
