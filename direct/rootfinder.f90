program rootfinder
! Implementation of a simple secant method to get a higher precision
! estimate of the zero for the geometric skewness in the duct as a function 
! of aspect ratio.

     implicit none
     
     ! Inputs
     integer             :: mMax,nMax,alphaMax
     double precision    :: r1, r2, tol
     character(len=1024) :: input_file, flow_type
     
     ! Internal
     double precision, dimension(:,:,:), allocatable   :: Alpha
     integer                                           :: funit,i,nIter,maxIter,nTerms
     character(len=1024)                               :: temp1,output_file,problem_type
     integer                                           :: t1, t2
     integer, dimension(:,:), allocatable              :: idxlist
     double precision, dimension(:), allocatable       :: uij_vals
     
     double precision                                  :: fudge
     
     
     double precision skew_Calc,skew_Calc2
     
     parameter(maxIter = 20)
     parameter(problem_type = "root")
     parameter(fudge = 0.4d0)
     
     ! Iterative variables. r3, skew3 will be the outputs.
     double precision                                  :: r3, skew1, skew2, skew3
     
     
     
     call system_clock(t1)
     
     ! Read inputs to get left/right bounds, tolerance
     call get_command_argument(1,input_file)

     funit=53
     open(funit,file=input_file)
          read(funit,*) r1
          read(funit,*) r2
          read(funit,*) tol
          read(funit,*) nTerms
          read(funit,*) flow_type
     close(funit)

     ! Allocate index set array and uij array.

     allocate(idxlist(1:nTerms,1:2))
     allocate(uij_vals(1:nTerms))

     ! alphaMax kept as legacy; we aren't precomputing alpha values here.
     alphaMax = 1
     
     allocate(Alpha(1:alphaMax,1:alphaMax,1:alphaMax))
     
     write(*,*) ""
!     write(*,"(A21)") "Precomputing alpha..."
!     call precalculate_Alpha(Alpha,alphaMax)



     if (trim(flow_type) .eq. "exact") then
          write(*,"(A49)") "Using the Fourier series expression for the flow."
     else if (trim(flow_type) .eq. "asymp") then
          write(*,*) "Using the asymptotic series expression for the flow."
     else if (trim(flow_type) .eq. "hybrid") then
          write(*,*) "Using a hybrid expression, with asymptotic below a threshold aratio."
     end if
     write(*,*) ""

     write(*,"(A29)") "Calculating initial values..."
     call coeff_chooser(nTerms,r1,idxlist,uij_vals)
     skew1 = skew_Calc(r1, nTerms, idxlist, uij_vals, Alpha, alphaMax, flow_type)
     call coeff_chooser(nTerms,r2,idxlist,uij_vals)
     skew2 = skew_Calc(r2, nTerms, idxlist, uij_vals, Alpha, alphaMax, flow_type)

     skew3 = skew2
     ! Set initial values.
     if (abs(skew2)<abs(skew1)) then
          r3 = r2
          skew3 = skew2
     else
          r3 = r1
          skew3 = skew1
     end if 
     
     
     nIter = 0
     
     write(*,"(A27,I4,A17)") "Rootfinding skewness using ",nTerms," summation terms."
     write(*,*) ""
     write(*,"(A4,A6,A23,A9,A23)") "Iter","      ","Aspect ratio","         ","Skewness"
     write(*,"(A4,A6,A23,A9,A23)") "----","      ","------------","         ","--------"
     write(*,"(I3,A7,ES23.16,A9,ES23.16)") nIter,"      ",r3,"    ",skew3
               
     ! Loop secant method until within tolerance. Function is nice enough that there shouldn't
     ! be any quirks in the iterations.
     do while((abs(skew3) > tol) .and. (nIter < maxIter))
          
          nIter = nIter + 1

          r3 = r1 - skew1*(r1-r2)/(skew1-skew2)


          call coeff_chooser(nTerms,r3,idxlist,uij_vals)
          skew3 = skew_Calc(r3, nTerms, idxlist, uij_vals, Alpha, alphaMax, flow_type)

          ! Ouput.
          write(*,"(I3,A7,ES23.16,A9,ES23.16)") nIter,"      ",r3,"    ",skew3
          
          ! Hand-me-downs.
          r2 = r1
          skew2 = skew1
          
          r1 = r3
          skew1 = skew3
          
     enddo


     ! Write solution to file.
!     call make_filename(nTerms, 1, flow_type, problem_type,output_file)
!     call write_outputs((/r3,skew3/), 1, mMax, nMax, flow_type, output_file)

     call system_clock(t2)
     
     write(*,*) ""
     write(*,"(A16,ES8.2,A9)") "Done in roughly ",dble(t2-t1)*1.0d-3," seconds."
     
     ! Cleanup.
     deallocate(Alpha,idxlist,uij_vals)
     
end program rootfinder
