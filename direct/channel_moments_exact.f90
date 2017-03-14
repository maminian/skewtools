program channel_moments_exact
! Spits out an h5 file containing a list of time values and 
! central moments values for the channel problem, where 
! we have explicit solutions.
!
! Additionally contains a manually entered cutoff to switch 
! from the short time asymptotics to the full solution.
! (There are convergence as well as numerical issues with the exact 
!  formulae at short time.)
!

implicit none

     integer                                      :: Nterms,nt,i
     double precision                             :: Pe,t,tmin,tmax,tscale,rel_err,m2,m3

     double precision, dimension(:), allocatable  :: tvals
     
     double precision, dimension(:), allocatable  :: AvMe,AvVa,AvSk,AvKu
     double precision, dimension(:), allocatable  :: asympMe,asympVa,asympSk,asympKu
     double precision, dimension(:), allocatable  :: hybridMe,hybridVa,hybridSk,hybridKu

     character(len=1024)                          :: fomat,inarg,fname,dsetname,descr
     
     ! Cutoff between using the short time asymptotics and exact formulae.
     double precision, parameter                  :: cutoff = 1.0d-3
     
     
     ! --------------------
     ! Set parameters via command-line arguments.
     
!     call get_command_argument(1,inarg)
!     read(inarg,*) Pe
!     call get_command_argument(2,inarg)
!     read(inarg,*) Nterms

     write(*,"(A15)",advance="no") "Peclet number? "
     read(*,*) Pe
     write(*,"(A27)",advance="no") "Number of terms in series? "
     read(*,*) Nterms
     write(*,"(A7)",advance="no") "t_min? "
     read(*,*) tmin
     write(*,"(A7)",advance="no") "t_max? "
     read(*,*) tmax
     write(*,"(A8)",advance="no") "tsteps? "
     read(*,*) nt
     
     ! ------------------
     
     
     allocate(tvals(nt))
     allocate(AvMe(nt),AvVa(nt),AvSk(nt),AvKu(nt))
     allocate(asympMe(nt),asympVa(nt),asympSk(nt),asympKu(nt))
     allocate(hybridMe(nt),hybridVa(nt),hybridSk(nt),hybridKu(nt))
     
     ! tmin*tscale**(nt-1) = tmax
     tscale = (tmax/tmin)**(1.0d0/dble(nt-1))

     write(*,*) ""
     write(*,"(A8,ES10.3)") "Peclet: ",Pe
     write(*,"(A23,I10)") "Terms used for Series: ",Nterms

     ! Unfortunately, the way I've set this up, I need to do this step manually.
     ! (too lazy to do it "cleanly.")
     
     i=1
     tvals(i) = tmin
     call channel_full_moments_exact(Nterms,Pe,tvals(i),AvMe(i),AvVa(i),AvSk(i),AvKu(i))
     call channel_full_moments_asymp(Pe,tvals(i),asympMe(i),asympVa(i),asympSk(i),asympKu(i))
     
     if (tvals(i) .lt. cutoff) then
          hybridMe(i) = asympMe(i)
          hybridVa(i) = asympVa(i)
          hybridSk(i) = asympSk(i)
          hybridKu(i) = asympKu(i)
     else
          hybridMe(i) = AvMe(i)
          hybridVa(i) = AvVa(i)
          hybridSk(i) = AvSk(i)
          hybridKu(i) = AvKu(i)
     end if
     
     ! Now start the proper loop.
     do i=2,nt
     
          tvals(i) = tvals(i-1)*tscale
          
          call channel_full_moments_exact(Nterms,Pe,tvals(i),AvMe(i),AvVa(i),AvSk(i),AvKu(i))
          call channel_full_moments_asymp(Pe,tvals(i),asympMe(i),asympVa(i),asympSk(i),asympKu(i))
          
          if (tvals(i) .lt. cutoff) then
               hybridMe(i) = asympMe(i)
               hybridVa(i) = asympVa(i)
               hybridSk(i) = asympSk(i)
               hybridKu(i) = asympKu(i)
          else
               hybridMe(i) = AvMe(i)
               hybridVa(i) = AvVa(i)
               hybridSk(i) = AvSk(i)
               hybridKu(i) = AvKu(i)
          end if
          
     end do
     
     ! Write stuff to file to plot later.
     
     fname = "direct_channel_moments.h5"
     call hdf_create_file(fname)
     
     dsetname = "Time"
     descr = "Time values evaluated"
     call hdf_add_1d_darray_to_file(nt,tvals,fname,dsetname,descr)

     dsetname = "Avgd_Mean"
     descr = "Average mean using series formula"
     call hdf_add_1d_darray_to_file(nt,AvMe,fname,dsetname,descr)

     dsetname = "Avgd_Variance"
     descr = "Average variance using series formula"
     call hdf_add_1d_darray_to_file(nt,AvVa,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness"
     descr = "Average skewness using series formula"
     call hdf_add_1d_darray_to_file(nt,AvSk,fname,dsetname,descr)
     
     dsetname = "Avgd_Kurtosis"
     descr = "Average kurtosis using series formula"
     call hdf_add_1d_darray_to_file(nt,AvKu,fname,dsetname,descr)

     !
     ! -------------------------------
     !
     
     dsetname = "Avgd_Mean_Asymp"
     descr = "Average mean using short time asymptotics"
     call hdf_add_1d_darray_to_file(nt,asympMe,fname,dsetname,descr)

     dsetname = "Avgd_Variance_Asymp"
     descr = "Average variance using short time asymptotics"
     call hdf_add_1d_darray_to_file(nt,asympVa,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness_Asymp"
     descr = "Average skewness using short time asymptotics"
     call hdf_add_1d_darray_to_file(nt,asympSk,fname,dsetname,descr)
     
     dsetname = "Avgd_Kurtosis_Asymp"
     descr = "Average kurtosis using short time asymptotics"
     call hdf_add_1d_darray_to_file(nt,asympKu,fname,dsetname,descr)

     !
     ! -------------------------------
     !

     dsetname = "Avgd_Mean_Hybrid"
     descr = "Average mean, using the asymptotic below a cutoff time"
     call hdf_add_1d_darray_to_file(nt,hybridMe,fname,dsetname,descr)
     
     dsetname = "Avgd_Variance_Hybrid"
     descr = "Average variance, using the asymptotic below a cutoff time"
     call hdf_add_1d_darray_to_file(nt,hybridVa,fname,dsetname,descr)
     
     dsetname = "Avgd_Skewness_Hybrid"
     descr = "Average skewness, using the asymptotic below a cutoff time"
     call hdf_add_1d_darray_to_file(nt,hybridSk,fname,dsetname,descr)
     
     dsetname = "Avgd_Kurtosis_Hybrid"
     descr = "Average kurtosis, using the asymptotic below a cutoff time"
     call hdf_add_1d_darray_to_file(nt,hybridKu,fname,dsetname,descr)

     !
     ! -------------------------------
     !

     dsetname = "nt"
     descr = "Number of time values"
     call hdf_add_1d_darray_to_file(1,nt,fname,dsetname,descr)
     
     dsetname = "Nterms"
     descr = "Number of terms used in the series solution"
     call hdf_add_1d_darray_to_file(1,Nterms,fname,dsetname,descr)
     
     dsetname = "Pe"
     descr = "Peclet value"
     call hdf_add_1d_darray_to_file(1,Pe,fname,dsetname,descr)
     
     deallocate(tvals,AvMe,AvVa,AvSk,AvKu,asympMe,asympVa,asympSk,&
                    asympKu,hybridMe,hybridVa,hybridSk,hybridKu)
     
end program channel_moments_exact

