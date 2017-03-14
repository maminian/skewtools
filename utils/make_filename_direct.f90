subroutine make_filename(mMax, nMax, flow_type, suffix, output_file)

! Makes a succinct filename describing the max summation indices,
! problem type, flow type.
!
! PROBLEM TYPE MUST BE 4 CHARACTERS LONG FOR NOW ("root" or "eval")

     implicit none

     ! Inputs
     integer                  :: mMax, nMax, sm, sn
     character(len=1024)      :: flow_type, suffix, temp1, output_file


     ! This is a "general" implementation, unless you want to sum to more than 10**10
     ! on every single index. (The "I1" in the fortran format descriptor assumes
     ! a 1-digit integer, which is sm and sn, the number of digits in the 
     ! truncation indices M and N.
     
     sm = floor(log10(dble(mMax)))+1               ! Count the number of digits in nMax.
     sn = floor(log10(dble(nMax)))+1               ! Count the number of digits in nMax.
     
     
     write(temp1,"(A5,I1,A5,I1,A13)") "(A1,I",sm,",A1,I",sn,",A5,A1,A4,A4)"   ! Create the formatting string with this number of digits.
     
     
     ! Write the max index on the first line, then the array in columns
     ! on the following lines.
     write(output_file,temp1) trim("m"),mMax,trim("n"),nMax,trim(flow_type),"_",trim(suffix),trim(".txt")   



end subroutine make_filename
