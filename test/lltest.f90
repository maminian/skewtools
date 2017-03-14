program lltest
! To prototype the 'linked list' idea.
use mtmod
implicit none
     integer(8), parameter              :: mtseed = 765874272
     integer, parameter                 :: na = 2**20
     double precision, parameter        :: cutoff = 0.1d0
     
     integer, dimension(na)             :: ll,lltemp
     double precision, dimension(na)    :: a
     integer                            :: i,j,k,nactive,moo,nstart,idx,ilast
     character(len=50)                  :: myfmt1,myfmt2
     
     
     myfmt1 = '(8F6.2)'
     myfmt2 = '(8I6)'
     
     call sgrnd(mtseed)
     
     do i=1,na
          lltemp(i) = i
     end do

     nactive = na
     

     j = 0
     do while (nactive .gt. 0)
          ll = lltemp
          call sample(na,a,nstart,nactive,ll)

          moo = nactive
          k = 1
          do i=1,moo
               idx = ll(i)

               if (a(idx) .lt. cutoff) then
                    nactive = nactive -1
               else
                    lltemp(k) = idx
                    k = k +1
               end if

          end do

          ll = lltemp
          j = j + 1

     end do

     write(*,'(A7,I5,A12)') "Done in",j,"iterations."
!     write(*,'(200F5.2)') a
     

     
end program lltest
!
!
!
subroutine sample(na,a,nstart,nactive,ll)
use mtmod
implicit none
     integer, intent(in)                          :: na,nstart,nactive
     integer, dimension(na), intent(in)           :: ll
     double precision, dimension(na), intent(inout) :: a
     integer i,idx
     
     do i=1,nactive
          idx = ll(i)
          a(idx) = grnd()

     end do

     
end subroutine sample
