program cchoose_test
! Quick testerino of coeff_chooser.

implicit none
     
     
     integer                                           :: nTerms,k,t1,t2
     double precision                                  :: aratio
     integer, dimension(:,:), allocatable              :: idx_list
     double precision, dimension(:), allocatable       :: uij_vals
     character(len=1024)                               :: tmp,cchoose_mode
     
     parameter(cchoose_mode="special")  ! Alternatively, "regular"
     
     call get_command_argument(1,tmp)
     read(tmp,*) nTerms
     call get_command_argument(2,tmp)
     read(tmp,*) aratio
     
     allocate( idx_list(1:nTerms,1:2) )
     allocate( uij_vals(1:nTerms) )

     call system_clock(t1)
     
     call coeff_chooser(nTerms,aratio,idx_list,uij_vals,cchoose_mode)
     
     call system_clock(t2)
     write(*,*) ""
     write(*,"(A25,D8.2,A9)") "Index set constructed in ",dble(t2-t1)*1.0d-3," seconds."
     
     call cchoose_print(nTerms,aratio,idx_list,uij_vals)


     deallocate(idx_list,uij_vals)
     
end program cchoose_test
!
! -----------------------
!
subroutine cchoose_print(nTerms,aratio,idx_list,uij_vals)
! Writing to file.

implicit none
     
     integer                                 :: nTerms
     double precision                        :: aratio
     integer, dimension(1:nTerms,1:2)        :: idx_list
     double precision, dimension(1:nTerms)   :: uij_vals

     integer                                 :: funit,i,k
     character(len=1024)                     :: out_file,fomat,temp
     
     ! Hacky file name, won't work for larger numbers.
     
     write(temp,*) "(A5,I",(floor(log10(dble(nTerms)))+1),",A4)"
     
     
     write(out_file,temp) trim("terms"),nTerms,trim(".txt")
     funit=53
     
     open(funit,file=out_file)
     
          ! Write solution parameters.
          write(funit,*) aratio
          write(funit,*) nTerms
          
          do i=1,nTerms
               write(funit,"(I5,I5,E26.17)") idx_list(i,1),idx_list(i,2),uij_vals(i)
          end do
          
     close(funit)
     
     write(*,*) "Results saved to ",trim(out_file)
     write(*,*) ""

end subroutine cchoose_print
