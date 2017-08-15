subroutine accumulate_moments_2d(tt_idx,ntt,nTot,X,Y,Z,yl,yr,zl,zr,means,vars,skews,&
               kurts,nby,nbz,means_sl,vars_sl,skews_sl,kurts_sl)

! A subroutine to be used in the main loop of the Monte Carlo code.
! Calculates the moments for a specific time.
!
! This version is for 2d (duct/ellipse) geometry.
!


implicit none
     

     ! Input arguments
     integer, intent(in)                                              :: tt_idx,ntt,nTot,nby,nbz
     double precision, dimension(nTot), intent(in)                    :: X,Y,Z
     double precision, intent(in)                                     :: yl,yr,zl,zr
     double precision, dimension(ntt), intent(inout)                  :: means,vars,skews,kurts
     double precision, dimension(ntt,nby,nbz), intent(inout)          :: means_sl,vars_sl,skews_sl,kurts_sl

     ! Internal
     integer                                                          :: kb,jb,fi,li,idx,nbins,i
     double precision, dimension(nTot)                                :: Xsorted
     integer, dimension(nTot)                                         :: Yidx,Zidx,bin_idxs

     double precision bdistfun_rt
     logical check_ic_duct

     nbins = nby*nbz
     
     ! Set values to zero just in case they're not already.
     ! They will also default to in the situation where 
     ! the bins are empty and they aren't handled in the loops below.
     means(tt_idx) = 0.0d0
     vars(tt_idx) = 0.0d0
     skews(tt_idx) = 0.0d0
     kurts(tt_idx) = 0.0d0
     do kb=1,nby
          do jb=1,nbz
               means_sl(tt_idx,kb,jb) = 0.0d0
               vars_sl(tt_idx,kb,jb) = 0.0d0
               skews_sl(tt_idx,kb,jb) = 0.0d0
               kurts_sl(tt_idx,kb,jb) = 0.0d0
          end do
     end do

     !
     ! Calculate all the moments. 
     !


     call moments(nTot,X,means(tt_idx),vars(tt_idx),skews(tt_idx),kurts(tt_idx))
     
     !
     ! Next, moments across (y,z) bins.
     !

     if (nby .gt. 1) then
          ! Start by binning independently in Y and Z.
          ! Enumerate the array of bins from 1,...,nby*nbz in the usual way.
          call uniform_bins_idx(nTot,Y,yl,yr,nby,Yidx)
          call uniform_bins_idx(nTot,Z,zl,zr,nbz,Zidx)

     

          
          bin_idxs = Yidx + (nby-1)*(Zidx-1) + 1
          write(*,*) minval(bin_idxs),maxval(bin_idxs),1,nby*nbz
          do i=1,nTot
               if ((bin_idxs(i) .lt. 1) .or. (bin_idxs(i) .gt. nby*nbz)) then
                    write(*,*) 
                    write(*,*) "DANGER DANGER WILL ROBINSON"
                    write(*,*) i,Y(i),Z(i),bdistfun_rt(Y(i),Z(i),0.5d0,0.2d0)
                    write(*,*) bin_idxs(i),nby*nbz
                    write(*,*) minval(bin_idxs),maxval(bin_idxs)
                    read(*,*) 
               end if
          end do
          ! Now sort the list of X positions based on their bin index.
!          write(*,*) minval(Y),maxval(Y),yl,yr
!          write(*,*) minval(Z),maxval(Z),zl,zr
!          write(*,*) minval(bin_idxs),maxval(bin_idxs)
          call sortpairs(nTot,X,Xsorted,bin_idxs,nbins)


          ! Finally, loop over the y and z bin indexes, calculating the moments for
          ! the corresponding subset of X values (which now are grouped together 
          ! in the array Xsorted).
     
          fi = 0    ! First index of active subset
          li = 0    ! Last index of active subset
          
          do jb=1,nbz
               do kb=1,nby
          
                    idx = (kb-1) + (nby-1)*(jb-1) + 1
 
                    ! Get the bounds on the active subset, 
                    ! assuming bin_idxs has been sorted already.
                    fi = li+1
                    li = fi

                    do while ( (bin_idxs(li) .eq. idx) )
                         li = li + 1

                         if (li .eq. (nTot+1)) then
                              ! If the lower index has passed the size of 
                              ! the array, break out. We also should be at the last 
                              ! index if we land in here. 
                              ! Otherwise something went badly wrong.           
                              go to 10
                         end if
                    end do

10                  continue

                    li = li-1 ! Necessary to decrement because of the bookkeeping.

                    ! Calculate the moments of this subset of particles!
                    call moments(li-fi+1,Xsorted(fi:li),means_sl(tt_idx,kb,jb),vars_sl(tt_idx,kb,jb),&
                                   skews_sl(tt_idx,kb,jb),kurts_sl(tt_idx,kb,jb))
                    
                    ! If we've exhausted the entries of bin_idxs (i.e., there are more
                    ! bins but they're empty), break out of the double loop.
                    if (li .eq. nTot) then
                         go to 20
                    end if
               end do
          end do
          
     end if
     
20   continue




end subroutine accumulate_moments_2d
