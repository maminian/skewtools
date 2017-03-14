program rect_duct
! Main program for calculating geometric skewness for the duct
! as a function of the aspect ratio of the rectangular cross section.

     implicit none
     
     ! Input variables; problem parameters to be read from file.
     integer                                           :: nAratios ! Number of aspect ratios
     double precision, dimension(1:1024)               :: aratios ! Aspect ratios to loop over. Note the hard cap.
     character(len=1024)                               :: input_file
     character(len=1024)                               :: flow_type ! Which expr. to use for the flow u.
     
     
     ! Internal variables.
     integer :: i, j, k, mMax, nMax, alphaMax, t1, t2, nTerms
     double precision, dimension(1:1024)               :: skewnesses
     character(len=1024)                               :: output_file, problem_type,arrayname,description
     
     character(len=1024)                               :: cchoose_mode
     integer, dimension(:,:), allocatable              :: idxlist
     double precision, dimension(:), allocatable       :: uij_vals
     
     ! Functions
     double precision skew_Calc

     parameter (problem_type="eval")
     parameter (cchoose_mode="special") ! Alternatively, "regular" for a rectangular sum.
     
     double precision, dimension(:,:,:), allocatable    :: Alpha
     
     ! Read input.txt to gather problem parameters.
     call get_command_argument(1,input_file)
     call read_inputs_direct(input_file, nAratios, aratios, mMax, nMax, flow_type)
     
     ! Allocate index set array and amn array.
     nTerms = mMax*nMax

     allocate(idxlist(1:nTerms,1:2))
     allocate(uij_vals(1:nTerms))

     ! Do a dry run of coeff_chooser in the minimum aspect ratio asked for, to 
     ! know the largest Alpha we need to precompute.
     call coeff_chooser(nTerms,0.2d0,idxlist,uij_vals)
     
     if (minval(aratios) .lt. 0.0d0) then
          write(*,*) "WARNING - NEGATIVE ASPECT RATIOS INPUT. CODE WILL LIKELY BREAK."
     end if     
     
     call system_clock(t1)
     
     call writeRectHeader(flow_type)
     
     ! Execute main loop, looping over aspect ratios.
     do i=1,nAratios
     
          ! Find the proper index set and calculate the relevant uij_vals to optimize convergence.
          call coeff_chooser(nTerms,aratios(i),idxlist,uij_vals,cchoose_mode)

          ! Then calculate skewness relative to this index set.
          ! The Alpha and 1 (=alphaMax) are not being used right now.

          skewnesses(i) = skew_Calc(aratios(i),nTerms,idxlist,uij_vals,Alpha,1,flow_type)

          
          write(*,"(ES11.3,A4,ES11.3)") aratios(i),"    ",skewnesses(i)

     end do
     
     ! Write the results to file.
     
     output_file = "direct.h5"
     call hdf_create_file(output_file)
     
     arrayname = "Aratio"
     description = "Aspect ratio of the duct"
     call hdf_add_2d_darray_to_file(1,nAratios,aratios,output_file,arrayname,description)
     
     arrayname = "Skew"
     description = "Geometric skewness"
     call hdf_add_2d_darray_to_file(1,nAratios,skewnesses,output_file,arrayname,description)
     
     call system_clock(t2)
     write(*,*) ""
     write(*,"(A16,ES8.3,A9)") "Done in roughly ",dble(t2-t1)*1.0d-3," seconds."
     
     
end program rect_duct
!
! ------------------------------------
!
subroutine writeRectHeader(flow_type)

     implicit none
     character(len=1024) flow_type

     if (trim(flow_type) .eq. "exact") then
          write(*,*) "Using the Fourier series expression for the flow."
     else if (trim(flow_type) .eq. "asymp") then
          write(*,*) "Using the asymptotic series expression for the flow."
     else if (trim(flow_type) .eq. "hybrid") then
          write(*,*) "Using a hybrid expression, with asymptotic below a threshold aratio."
     end if
     write(*,*) ""
     
     write(*,*) "Calculating skewness values... "
     write(*,*) ""
     write(*,"(A12,A5,A8)") "Aspect ratio","     ","Skewness"
     write(*,"(A13,A4,A8)") "-------------","    ","--------"
     write(*,*)
     
end subroutine writeRectHeader
