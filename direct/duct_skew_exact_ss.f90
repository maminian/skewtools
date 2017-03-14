program duct_skew_exact_ss
! Spits out an h5 file containing a list of time values and 
! skewness values, whose particular form is specified by 
! an input file.
!
! For now, specified manually.

implicit none

     integer                                      :: nTerms,nt,i
     double precision                             :: aratio,Pe,t,tmin,tmax,tscale
     double precision, dimension(:), allocatable  :: tvals,asympSk,uij_vals
     integer, dimension(:,:), allocatable         :: idxlist
     character(len=1024)                          :: fomat,inarg,fname,dsetname,descr
     
     double precision, dimension(4)               :: moment_vals
     
     ! Chosen value of the Laplacian. The calculations have it set as -2, but 
     ! that's with an ad-hoc multiplier
     double precision, parameter                  :: laplacianu = -2.0d0
     
     double precision duct_skew_asymp_skew_calc
     double precision compute_u1_noncenter_exact,compute_u2_noncenter_exact,&
                         compute_u3_noncenter_exact
     
     ! --------------------
     ! Set parameters via command-line arguments.
     
!     call get_command_argument(1,inarg)
!     read(inarg,*) Pe
!     call get_command_argument(2,inarg)
!     read(inarg,*) nTerms

     write(*,"(A14)",advance="no") "Aspect ratio? "
     read(*,*) aratio
     write(*,"(A15)",advance="no") "Peclet number? "
     read(*,*) Pe
     write(*,"(A27)",advance="no") "Number of terms in series? "
     read(*,*) nTerms
     write(*,"(A7)",advance="no") "t_min? "
     read(*,*) tmin
     write(*,"(A7)",advance="no") "t_max? "
     read(*,*) tmax
     write(*,"(A8)",advance="no") "tsteps? "
     read(*,*) nt
     write(*,"(A17)",advance="no") "Output filename? "
     read(*,*) fname
     
     ! ------------------
     
     
     allocate(tvals(nt),asympSk(nt))
     allocate(idxlist(nTerms,2),uij_vals(nTerms))
     
     ! tmin*tscale**(nt-1) = tmax
     tscale = (tmax/tmin)**(1.0d0/dble(nt-1))

     write(*,*) ""
     write(*,"(A14,ES10.3)") "Aspect ratio: ",aratio
     write(*,"(A8,ES10.3)") "Peclet: ",Pe
     write(*,"(A23,I10)") "Terms used for Series: ",nTerms

     write(*,"(A57)") "Computing u values and averages needed for the moments..."
     call coeff_chooser(nTerms,aratio,idxlist,uij_vals,"special")
     
     uij_vals = uij_vals
     
     moment_vals(1) = compute_u1_noncenter_exact(nTerms,idxlist,uij_vals)
     moment_vals(2) = compute_u2_noncenter_exact(nTerms,idxlist,uij_vals)
     moment_vals(3) = compute_u3_noncenter_exact(nTerms,idxlist,uij_vals)
     moment_vals(4) = laplacianu
     
     write(*,"(A5)") "done."

     fomat = "(A10,A3,A10,A3,A10,A3,A10)"
     write(*,*) ""
     write(*,fomat) "t_value","Asymp"
     write(*,fomat) "-------","   ","-------"

     i=1
     fomat = "(ES10.3,A3,ES10.3)"
     tvals(i) = tmin
     asympSk(i) = duct_skew_asymp_skew_calc(moment_vals,Pe,tvals(i))
     
     
     write(*,fomat) tvals(i),"   ",asympSk(i)
     
     do i=2,nt
     
          tvals(i) = tvals(i-1)*tscale
          asympSk(i) = duct_skew_asymp_skew_calc(moment_vals,Pe,tvals(i))
          
          
          write(*,fomat) tvals(i),"   ",asympSk(i)
          
     end do
     
     write(*,*) ""
     
     ! Write stuff to file to plot later.
     
!     fname = "direct_duct_skew.h5"
     call hdf_create_file(fname)
     
     dsetname = "Time"
     descr = "Time values evaluated"
     call hdf_add_1d_darray_to_file(nt,tvals,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness_Asymp"
     descr = "Average skewness using series formula"
     call hdf_add_1d_darray_to_file(nt,asympSk,fname,dsetname,descr)
     
     dsetname = "nt"
     descr = "Number of time values"
     call hdf_add_1d_darray_to_file(1,nt,fname,dsetname,descr)
     
     dsetname = "Peclet"
     descr = "Peclet number"
     call hdf_add_1d_darray_to_file(1,Pe,fname,dsetname,descr)
     
     dsetname = "aratio"
     descr = "Aspect ratio"
     call hdf_add_1d_darray_to_file(1,aratio,fname,dsetname,descr)
     
     dsetname = "nTerms"
     descr = "Number of terms used in the series solution"
     call hdf_add_1d_darray_to_file(1,nTerms,fname,dsetname,descr)
     
     deallocate(tvals,asympSk,idxlist,uij_vals)
     
end program duct_skew_exact_ss
!
! ----------------------------------------------------------------------------------------
!
double precision function duct_skew_asymp_skew_calc(moment_vals,Pe,t)
implicit none

     double precision, dimension(4)     :: moment_vals
     integer             :: nTerms
     double precision    :: m1,m2,m3,Pe,t
     
     double precision    :: u1avg,u2avg,u3avg,laplu
     double precision    :: u2zeroavg,u3zeroavg
     
     u1avg = moment_vals(1)
     u2avg = moment_vals(2)
     u3avg = moment_vals(3)
     laplu = moment_vals(4)
     
     u2zeroavg = u2avg - u1avg**2
     u3zeroavg = u3avg - 3.0d0*u1avg*u2avg + 2.0d0*u1avg**3
     
     m1 = 0.0d0
!     m2 = 2.0d0*t + Pe**2*( t**2*u2zeroavg + t**3*(laplu**3/3.0d0)*u1avg )
!     m3 = Pe**3*( t**3*(u3zeroavg) + t**4*(laplu/2.0d0)*(u2avg-2.0d0*u1avg**2) )

     m2 = (2.0d0*t)/(Pe**2) + ( t**2*u2zeroavg + t**3*(laplu**3/3.0d0)*u1avg )
     m3 = ( t**3*(u3zeroavg) + t**4*(laplu/2.0d0)*(u2avg-2.0d0*u1avg**2) )
     
     duct_skew_asymp_skew_calc = m3/m2**1.5d0
     
end function duct_skew_asymp_skew_calc
