subroutine accumulate_moments_1d(tt_idx,ntt,nTot,X,Y,a,means,vars,skews,&
               kurts,nby,means_sl,vars_sl,skews_sl,kurts_sl)
! A subroutine to be used in the main loop of the Monte Carlo code.
! Calculates the moments for a specific time.
!
! This version is for 1d (channel) geometry.
!

implicit none

     ! Input arguments
     integer, intent(in)                                         :: tt_idx,ntt,nTot,nby
     double precision, dimension(nTot), intent(in)               :: X,Y
     double precision, intent(in)                                :: a
     double precision, dimension(ntt), intent(inout)             :: means,vars,skews,kurts
     double precision, dimension(ntt,nby), intent(inout)         :: means_sl,vars_sl,skews_sl,kurts_sl


     ! Internal
     integer                                                     :: kb,bin_count,fi,li,idx,nbins 
     double precision                                            :: bin_lo,bin_hi

     !
     ! Array to hold X values located in a certain bin, for a fixed round.
     ! The array size here is really just an upper bound.
     double precision, dimension(nTot)                                :: Xsorted
     integer, dimension(nTot)                                         :: Yidx,bin_idxs
     
     nbins = nby

     ! Set values to zero just in case they're not already.
     means(tt_idx) = 0.0d0
     vars(tt_idx) = 0.0d0
     skews(tt_idx) = 0.0d0
     kurts(tt_idx) = 0.0d0
     
     do kb=1,nby
          means_sl(tt_idx,kb) = 0.0d0
          vars_sl(tt_idx,kb) = 0.0d0
          skews_sl(tt_idx,kb) = 0.0d0
          kurts_sl(tt_idx,kb) = 0.0d0
     end do
     
     !
     ! Calculate all the moments. 
     !

     call moments(nTot,X,means(tt_idx),vars(tt_idx),skews(tt_idx),kurts(tt_idx))

     if (nby .gt. 0) then
          ! Start by binning in Y.
          ! Enumerate the array of bins from 1,...,nby.
          call uniform_bins_idx(nTot,Y,-a,a,nby,Yidx)
          
          bin_idxs = Yidx
          
          ! Now sort the list of X positions based on their bin index.

          call sortpairs(nTot,X,Xsorted,bin_idxs,nbins)

          ! Finally, loop over the y and z bin indexes, calculating the moments for
          ! the corresponding subset of X values (which now are grouped together 
          ! in the array Xsorted).
     
          fi = 0    ! First index of active subset
          li = 0    ! Last index of active subset
          

          do kb=1,nby
     
               idx = kb

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

10             continue

               li = li-1 ! Necessary to decrement because of the bookkeeping.

               ! Calculate the moments of this subset of particles!
               call moments(li-fi+1,Xsorted(fi:li),means_sl(tt_idx,kb),vars_sl(tt_idx,kb),&
                              skews_sl(tt_idx,kb),kurts_sl(tt_idx,kb))
               
               ! If we've exhausted the entries of bin_idxs (i.e., there are more
               ! bins but they're empty), break out of the double loop.
               if (li .eq. nTot) then
                    go to 20
               end if
          end do
                    
     end if


20   continue

end subroutine accumulate_moments_1d
