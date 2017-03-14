subroutine read_inputs_direct(input_file,nAratios,aratios,mMax,nMax,flow_type)

implicit none
     
     ! Internal vars
     integer                                           :: funit
     
     ! Inputs
     character(len=1024)                               :: input_file
     
     ! Inputs/outputs
     integer                                           :: nAratios,mMax,nMax
     double precision, dimension(1:1024)               :: aratios
     character(len=1024)                               :: flow_type

     funit = 53 ! Arbitrary
     
     open(funit,file=input_file)
     
          read(funit,*) nAratios 
          read(funit,*) aratios(1:nAratios)
          read(funit,*) mMax
          read(funit,*) nMax
          read(funit,*) flow_type
          
     close(funit)

end subroutine read_inputs_direct
