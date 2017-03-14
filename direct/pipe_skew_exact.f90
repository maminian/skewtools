program pipe_skew_exact
! Spits out an h5 file containing a list of time values and 
! skewness values, whose particular form is specified by 
! an input file.
!
! For now, specified manually.

implicit none

     integer                                      :: nTerms,nt,i
     double precision                             :: Pe,t,tmin,tmax,tscale,rel_err,m2,m3
     double precision, dimension(:), allocatable  :: AvSk,tvals,asympSk,hybridSk
     character(len=1024)                          :: fomat,inarg,fname,dsetname,descr
     
     ! Cutoff between using the short time asymptotics and exact formula.
     double precision, parameter                  :: cutoff = 1.0d6
     
     double precision pipe_skew_exact_skew_calc, pipe_skew_asymp_skew_calc
     
     ! --------------------
     ! Set parameters via command-line arguments.
     
!     call get_command_argument(1,inarg)
!     read(inarg,*) Pe
!     call get_command_argument(2,inarg)
!     read(inarg,*) nTerms

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
     
     ! ------------------
     
     
     allocate(tvals(nt),AvSk(nt),asympSk(nt),hybridSk(nt))
     
     ! tmin*tscale**(nt-1) = tmax
     tscale = (tmax/tmin)**(1.0d0/dble(nt-1))

     write(*,*) ""
     write(*,"(A8,ES10.3)") "Peclet: ",Pe
     write(*,"(A23,I10)") "Terms used for Series: ",nTerms

     fomat = "(A10,A3,A10,A3,A10,A3,A10)"
     write(*,*) ""
     write(*,fomat) "t_value","   ","Series","   ","Asymp","   ","Rel_Error"
     write(*,fomat) "-------","   ","-------","   ","-------","   ","-------"

     i=1
     fomat = "(ES10.3,A3,ES10.3,A3,ES10.3,A3,ES10.3)"
     tvals(i) = tmin
     AvSk(i) = pipe_skew_exact_skew_calc(nTerms,Pe,tvals(i))
     asympSk(i) = pipe_skew_asymp_skew_calc(Pe,tvals(i))
     
     if (tvals(i) .lt. cutoff) then
          hybridSk(i) = asympSk(i)
     else
          hybridSk(i) = AvSk(i)
     end if
     
     rel_err = abs((asympSk(i)-AvSk(i))/AvSk(i))
     
     write(*,fomat) tvals(i),"   ",AvSk(i),"   ",asympSk(i),"   ",rel_err
     
     do i=2,nt
     
          tvals(i) = tvals(i-1)*tscale
          AvSk(i) = pipe_skew_exact_skew_calc(nTerms,Pe,tvals(i))
          asympSk(i) = pipe_skew_asymp_skew_calc(Pe,tvals(i))
          
          if (tvals(i) .lt. cutoff) then
               hybridSk(i) = asympSk(i)
          else
               hybridSk(i) = AvSk(i)
          end if
          
          rel_err = abs((asympSk(i)-AvSk(i))/AvSk(i))
          
          write(*,fomat) tvals(i),"   ",AvSk(i),"   ",asympSk(i),"   ",rel_err
          
     end do
     
     write(*,*) ""
     
     ! Write stuff to file to plot later.
     
     fname = "direct_pipe_skew.h5"
     call hdf_create_file(fname)
     
     dsetname = "Time"
     descr = "Time values evaluated"
     call hdf_add_1d_darray_to_file(nt,tvals,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness"
     descr = "Average skewness using series formula"
     call hdf_add_1d_darray_to_file(nt,AvSk,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness_Asymp"
     descr = "Average skewness using series formula"
     call hdf_add_1d_darray_to_file(nt,asympSk,fname,dsetname,descr)

     dsetname = "Avgd_Skewness_Hybrid"
     descr = "Average skewness, using the asymptotic below a cutoff time"
     call hdf_add_1d_darray_to_file(nt,hybridSk,fname,dsetname,descr)

     dsetname = "nt"
     descr = "Number of time values"
     call hdf_add_1d_darray_to_file(1,nt,fname,dsetname,descr)
     
     dsetname = "nTerms"
     descr = "Number of terms used in the series solution"
     call hdf_add_1d_darray_to_file(1,nTerms,fname,dsetname,descr)
     
     deallocate(tvals,AvSk,asympSk)
     
end program pipe_skew_exact
!
! ----------------------------------------------------------------------------------------
!
double precision function pipe_skew_asymp_skew_calc(Pe,t)
implicit none
     double precision    :: m1,m2,m3,Pe,t
     
     double precision asymp_pipe_m1,asymp_pipe_m2,asymp_pipe_m3
     
     m1 = asymp_pipe_m1(Pe,t)
     m2 = asymp_pipe_m2(Pe,t)
     m3 = asymp_pipe_m3(Pe,t)

     
     pipe_skew_asymp_skew_calc = (m3-3.0d0*m1*m2+2.0d0*m1**3)/(m2-m1**2)**1.5d0
     
end function pipe_skew_asymp_skew_calc
!
! ----------------------------------------------------------------------------------------
!
double precision function pipe_skew_exact_skew_calc(nTerms,Pe,t)
implicit none
     integer             :: nTerms
     double precision    :: m1,m2,m3,Pe,t
     
     double precision exact_pipe_m1,exact_pipe_m2,exact_pipe_m3
     
     m1 = exact_pipe_m1(nTerms,Pe,t)
     m2 = exact_pipe_m2(nTerms,Pe,t)
     m3 = exact_pipe_m3(nTerms,Pe,t)

     
     pipe_skew_exact_skew_calc = (m3-3.0d0*m1*m2+2.0d0*m1**3)/(m2-m1**2)**1.5d0
     
end function pipe_skew_exact_skew_calc

