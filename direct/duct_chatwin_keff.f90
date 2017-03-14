program duct_chatwin_keff
! 
! Applying Chatwin's (1982) formulae for effective diffusivity 
! in the rectangular duct, to reconcile the difference between 
! all my results (consistent between summation, MC, and FEM) 
! and his summation.
!

implicit none

     character(len=1024)                          :: input_file,output_file,arrayname,description
     integer                                      :: funit,naratios,nmax,i
     double precision, dimension(:), allocatable  :: aratios,keffs

     double precision calculate_duct_keff_chat


     ! Read inputs to get left/right bounds, tolerance
     call get_command_argument(1,input_file)

     funit=53
     open(funit,file=input_file)
          read(funit,*) naratios 
          allocate(aratios(naratios),keffs(naratios))
          
          read(funit,*) aratios(1:naratios)
          read(funit,*) nmax
     close(funit)

     ! Loop over the aspect ratios.
     do i=1,naratios
          keffs(i) = calculate_duct_keff_chat(aratios(i),nmax)
          write(*,*) i,aratios(i),keffs(i)
     end do

     ! Write the output.
     output_file = "duct_chatwin_keff.h5"
     call hdf_create_file(output_file)
     
     arrayname = "aratios"
     description = "Aspect ratio of the duct"
     call hdf_add_1d_darray_to_file(naratios,aratios,output_file,arrayname,description)
     
     arrayname = "keffs"
     description = "Dimensionless effective diffusivity"
     call hdf_add_1d_darray_to_file(naratios,keffs,output_file,arrayname,description)

     deallocate(aratios,keffs)


end program duct_chatwin_keff
!
! ----------------------------------------
!
double precision function calculate_duct_keff_chat(aratio,nmax)
implicit none
     double precision         :: aratio,a,b,term1,term2,term3,d0
     integer                  :: nmax

     double precision duct_keff_chat_sum1,duct_keff_chat_sum2,duct_keff_chat_sum3

     a = 1.0d0/aratio
     b = 1.0d0

     
     term1 = duct_keff_chat_sum1(a,b,nmax)
     term2 = duct_keff_chat_sum2(a,b,nmax)
     term3 = duct_keff_chat_sum3(a,b,nmax)

     d0 = b**6/(30240.0d0)

     calculate_duct_keff_chat = (term1 + term2 + term3) / d0
     
end function calculate_duct_keff_chat
!
! ----------------------------------------
!
double precision function duct_keff_chat_sum1(a,b,nmax)
implicit none
     double precision    :: a,b,pi
     integer             :: nmax,p
     double precision chat_wp0

     pi = 4.0d0*datan(1.0d0)

     duct_keff_chat_sum1 = 0.0d0
     do p=2,nmax,2
          duct_keff_chat_sum1 = duct_keff_chat_sum1 + ( a/(p*pi)*chat_wp0(a,b,p,nmax) )**2
     end do

     duct_keff_chat_sum1 = 0.5d0*duct_keff_chat_sum1

end function duct_keff_chat_sum1
!
! ----------------------------------------
!
double precision function duct_keff_chat_sum2(a,b,nmax)
implicit none
     double precision    :: a,b,pi
     integer             :: nmax,q
     double precision chat_w0q

     pi = 4.0d0*datan(1.0d0)

     duct_keff_chat_sum2 = 0.0d0
     do q=2,nmax,2
          duct_keff_chat_sum2 = duct_keff_chat_sum2 + ( b/(q*pi)*chat_w0q(a,b,q,nmax) )**2
     end do

     duct_keff_chat_sum2 = 0.5d0*duct_keff_chat_sum2

end function duct_keff_chat_sum2
!
! ----------------------------------------
!
double precision function duct_keff_chat_sum3(a,b,nmax)
implicit none
     double precision    :: a,b,pi
     integer             :: nmax,p,q
     double precision chat_wpq

     pi = 4.0d0*datan(1.0d0)

     duct_keff_chat_sum3 = 0.0d0
     do p=2,nmax,2
          do q=2,nmax,2
               duct_keff_chat_sum3 = duct_keff_chat_sum3 + chat_wpq(a,b,p,q,nmax)**2 / &
                                   ( (p*pi/a)**2 + (q*pi/b)**2 )
          end do
     end do

     duct_keff_chat_sum3 = 0.25d0*duct_keff_chat_sum3

end function duct_keff_chat_sum3
!
! ----------------------------------------
!
double precision function chat_wp0(a,b,p,nmax)
implicit none
     double precision         :: a,b,pi
     integer                  :: p,nmax,n

     pi = 4.0d0*dtan(1.0d0)

     chat_wp0 = 0.0d0
     do n=1,4*nmax,2
          chat_wp0 = chat_wp0 + dtanh(n*pi*a/(2*b)) / (n**3*((n*a)**2 + (p*b)**2))
     end do

     chat_wp0 = chat_wp0*(-32.0d0*a*b**3/pi**5)
end function chat_wp0
!
! ----------------------------------------
!
double precision function chat_w0q(a,b,q,nmax)
implicit none
     double precision         :: a,b,pi
     integer                  :: q,nmax,n

     pi = 4.0d0*dtan(1.0d0)

     chat_w0q = 0.0d0
     do n=1,4*nmax,2
          chat_w0q = chat_w0q + ( q**2/dble(n**2*(n**2-q**2)) ) * &
                                   ( dtanh(n*pi*a/(2*b)) / (n*pi*a/(2*b)) )
     end do

     chat_w0q = chat_w0q*(2.0d0/pi**2)
     chat_w0q = chat_w0q + 1.0d0
     chat_w0q = chat_w0q*(-2.0d0*b**2/(pi*q)**2)

end function chat_w0q
!
! ----------------------------------------
!
double precision function chat_wpq(a,b,p,q,nmax)
implicit none
     double precision         :: a,b,pi
     integer                  :: p,q,nmax,n

     pi = 4.0d0*dtan(1.0d0)

     chat_wpq = 0.0d0
     do n=1,4*nmax,2
          chat_wpq = chat_wpq + dtanh(n*pi*a/(2*b)) / &
                         dble( n*(n**2-q**2) * ((n*a)**2 + (p*b)**2) )
     end do

     chat_wpq = -64*a*b**3/(pi**5)*chat_wpq

end function chat_wpq
