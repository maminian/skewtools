subroutine uniform_bins_idx(nx,X,xmin,xmax,nb,Xidx)
!
! Given a double array X of size nx, and lower and upper 
! bounds xmin and xmax and number of bins nb, 
! output an integer array Xidx corresponding to 
! the bin number assignment, assuming uniformly spaced bins.
!
! In other words, maps the elements 'linearly' to the integers 
! 0,1,...,nb-1.

implicit none
     integer, intent(in)                          :: nx,nb
     double precision, intent(in), dimension(nx)  :: X
     double precision, intent(in)                 :: xmin,xmax
     
     integer, intent(out), dimension(nx)          :: Xidx
     
     double precision                             :: width
     
     width = xmax-xmin

     Xidx = floor(((nb-1)/width)*(X - xmin))+1

end subroutine uniform_bins_idx
