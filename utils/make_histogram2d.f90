subroutine make_histogram2d(n,X,Y,nhb,nby,hcx,hcy,heights)
! Bins the pair of arrays X,Y in two dimensions, with nhb bins in the 
! x direction and nby bins in the y direction. hcx and hcy track 
! the locations of bin *centers*, not boundaries. heights is the count, 
! non-normalized.
!

implicit none

     integer, intent(in)                                    :: n,nhb,nby
     double precision, dimension(n), intent(in)             :: X,Y
     double precision, dimension(nhb), intent(out)          :: hcx
     double precision, dimension(nby), intent(out)          :: hcy
     double precision, dimension(nhb,nby), intent(out)      :: heights
     
     integer, dimension(n)                                  :: Xidx,Yidx,Bidx
     double precision                                       :: xmin,xmax,dbx,ymin,ymax,dby
     integer                                                :: j,i
     
     ! Set up the binning.
     xmin = minval(X)
     xmax = maxval(X)

     if (xmin .eq. xmax) then
          xmin = xmin - 1.0d0
          xmax = xmax + 1.0d0
     end if

     ymin = -1.0d0
     ymax = 1.0d0

     dbx = (xmax-xmin)/nhb
     dby = (ymax-ymin)/nby

     do j=1,nhb
          hcx(j) = xmin + (j-0.5d0)*dbx
     end do

     do j=1,nby
          hcy(j) = ymin + (j-0.5d0)*dby
     end do

     do i=1,nhb
          do j=1,nby
               heights(i,j) = 0.0d0
          end do
     end do


     call uniform_bins_idx(n,X,xmin,xmax,nhb,Xidx)
     call uniform_bins_idx(n,Y,ymin,ymax,nby,Yidx)

     do j=1,n
          heights(Xidx(j),Yidx(j)) = heights(Xidx(j),Yidx(j)) + 1.0d0
     end do


!     heights = heights/sum(heights)

end subroutine make_histogram2d
