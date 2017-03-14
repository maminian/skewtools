subroutine write_outputs_mc(Ntrials,Pe,aratio,dt,nArray,nVars,array,filename)
! Not to be confused with the similarly named function when
! doing the exact calculation in inviscid setting.

implicit none
     integer                                           :: Ntrials,nArray,nVars
     double precision                                  :: Pe,aratio,dt
     double precision, dimension(1:nArray,1:nVars)     :: array
     character(len=1024)                               :: filename

     integer                                 :: funit,k
     
     funit = 53
     
     open(funit,file=filename)
     
          write(funit,*) Ntrials
          write(funit,*) Pe
          write(funit,*) aratio
          write(funit,*) dt
          write(funit,*) nArray
          
          do k=0,nArray-1
               write(funit,*) k*dt, array(k+1,1:nvars)
          end do
     
     close(funit)

end subroutine write_outputs_mc
