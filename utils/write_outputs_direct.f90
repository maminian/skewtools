subroutine write_outputs(output_combo,nAratios,nTerms,flow_type,output_file)

     implicit none
     
     ! Internal vars
     integer                                           :: funit,i
     character(len=1024)                               :: output_file
     
     ! Inputs
     integer                                           :: nAratios,nTerms
     double precision, dimension(1:nAratios,1:2)       :: output_combo
     character(len=1024)                               :: flow_type


     ! Make a filename. Messier than I'd like. Eh.

     funit = 53 ! Arbitrary
     
     open(funit,file=output_file)
     
          ! Write solution parameters.
          write(funit,*) nTerms
          write(funit,*) trim(flow_type)
          write(funit,*) nAratios       ! Or, more generally, the number of lines below this one.
          
          do i=1,nAratios
               write(funit,"(ES26.17,ES26.17)") output_combo(i,1),output_combo(i,2)
          end do
          
     close(funit)

end subroutine write_outputs
