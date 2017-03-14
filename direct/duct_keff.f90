program duct_keff
!
! Calculates the effective diffusivity coefficient for a range of aspect ratios in the duct.
! Main purpose is to compare with Chatwin (1982).
!

implicit none

     character(len=1024)                          :: input_file,output_file,arrayname,description
     integer                                      :: funit,naratios,nmax,i
     double precision, dimension(:), allocatable  :: aratios,keffs,keffs_asymp

     double precision calculate_duct_keff,calculate_duct_keff_asymp


     ! Read inputs to get left/right bounds, tolerance
     call get_command_argument(1,input_file)

     funit=53
     open(funit,file=input_file)
          read(funit,*) naratios 
          allocate(aratios(naratios),keffs(naratios),keffs_asymp(naratios))
          
          read(funit,*) aratios(1:naratios)
          read(funit,*) nmax
     close(funit)

     ! Loop over the aspect ratios.
     do i=1,naratios
          write(*,*) i,aratios(i)
          keffs(i) = calculate_duct_keff(aratios(i),nmax)
          keffs_asymp(i) = calculate_duct_keff_asymp(aratios(i),nmax)
     end do

     ! Write the output.
     output_file = "duct_keff.h5"
     call hdf_create_file(output_file)
     
     arrayname = "aratios"
     description = "Aspect ratio of the duct"
     call hdf_add_1d_darray_to_file(naratios,aratios,output_file,arrayname,description)
     
     arrayname = "keffs"
     description = "Dimensionless effective diffusivity"
     call hdf_add_1d_darray_to_file(naratios,keffs,output_file,arrayname,description)

     arrayname = "keffs_asymp"
     description = "Dimensionless effective diffusivity asymptotic for small aspect ratio"
     call hdf_add_1d_darray_to_file(naratios,keffs_asymp,output_file,arrayname,description)

     deallocate(aratios,keffs)

end program duct_keff
!
! ----------------------------------------------------------
!
double precision function calculate_duct_keff(aratio,nmax)
implicit none
     
     double precision                        :: term0,term1,term2,term3,term4,aratio
     integer                                 :: nmax
     double precision, dimension(nmax,nmax)  :: umn

     double precision duct_keff_sum1,duct_keff_sum2,duct_keff_sum3,duct_keff_sum4


     term0 = 8.0d0/945.0d0


     term1 = duct_keff_sum1(aratio,nmax)
     term2 = duct_keff_sum2(aratio,nmax)
     term3 = duct_keff_sum3(aratio,nmax)
     term4 = duct_keff_sum4(aratio,nmax)

     calculate_duct_keff = term0 + term1 + term2 + term3 + term4
     
end function calculate_duct_keff
!
! ----------------------------------------------------------
!
double precision function calculate_duct_keff_asymp(aratio,nmax)
! Exactly the same as above, with a few terms removed.
implicit none
     
     double precision                        :: term0,term1,term2,term3,term4,aratio
     integer                                 :: nmax
     double precision, dimension(nmax,nmax)  :: umn

     double precision duct_keff_sum1,duct_keff_sum2,duct_keff_sum3,duct_keff_sum4


     term0 = 8.0d0/945.0d0


     term1 = duct_keff_sum1(aratio,nmax)
     term2 = duct_keff_sum2(aratio,nmax)
     term3 = duct_keff_sum3(aratio,nmax)
     term4 = duct_keff_sum4(aratio,nmax)

     calculate_duct_keff_asymp = term0 + term2 + term3 + term4
     
end function calculate_duct_keff_asymp
!
! ----------------------------------------------------------
!
double precision function duct_keff_sum1(aratio,nmax)
implicit none
     double precision    :: aratio,pi
     integer             :: nmax,m
     double precision dke_um0

     pi = 4.0d0*datan(1.0d0)

     duct_keff_sum1 = 0.0d0
     do m=1,nmax
          duct_keff_sum1 = duct_keff_sum1 + ( dke_um0(aratio,m,nmax)/(pi*m) )**2
     end do

     duct_keff_sum1 = 0.5d0*duct_keff_sum1

end function duct_keff_sum1
!
! ----------------------------------------------------------
!
double precision function duct_keff_sum2(aratio,nmax)
implicit none
     double precision    :: aratio,pi
     integer             :: nmax,n
     double precision dke_u0n

     pi = 4.0d0*datan(1.0d0)

     duct_keff_sum2 = 0.0d0
     do n=1,nmax
          duct_keff_sum2 = duct_keff_sum2 + ( dke_u0n(aratio,n,nmax)/(pi*n) )**2
     end do

     duct_keff_sum2 = 0.5d0/(aratio**2)*duct_keff_sum2

end function duct_keff_sum2
!
! ----------------------------------------------------------
!
double precision function duct_keff_sum3(aratio,nmax)
implicit none
     double precision    :: aratio,pi
     integer             :: nmax,m,n
     double precision dke_umn

     pi = 4.0d0*datan(1.0d0)

     duct_keff_sum3 = 0.0d0
     do m=1,nmax
          do n=1,nmax
               duct_keff_sum3 = duct_keff_sum3 + dke_umn(aratio,m,n,nmax)**2/(pi**2*(m**2 + (aratio*n)**2))
          end do
     end do

     duct_keff_sum3 = 0.25d0*duct_keff_sum3

end function duct_keff_sum3
!
! ----------------------------------------------------------
!
double precision function duct_keff_sum4(aratio,nmax)
implicit none
     double precision    :: aratio,pi
     integer             :: nmax,m,n
     double precision dke_um0

     pi = 4.0d0*datan(1.0d0)

     duct_keff_sum4 = 0.0d0
     do m=1,nmax
          duct_keff_sum4 = duct_keff_sum4 + dke_um0(aratio,m,nmax)*(-4)*(-1)**m/(pi**4*m**4)
     end do

     duct_keff_sum4 = 1.0d0*duct_keff_sum4

end function duct_keff_sum4
!
! ----------------------------------------------------------
!
double precision function dke_um0(aratio,m,nmax)
implicit none
     double precision    :: aratio,pi,numer,denom
     integer             :: m,n,k,nmax

     pi = 4.0d0*datan(1.0d0)
     
     dke_um0 = 0.0d0
     ! As a precaution, need to go well beyond m for the denominator (k-1/2)**2 - m**2 to force convergence.
     do k=1,2*nmax
          numer = (-8*(-1)**m*aratio*dtanh(((-0.5 + k)*pi)/aratio))
          denom = ((-0.5 + k)**3*((-0.5 + k)**2 - m**2)*pi**5)

          dke_um0 = dke_um0 + numer/denom
     end do

end function dke_um0
!
! ----------------------------------------------------------
!
double precision function dke_u0n(aratio,n,nmax)
implicit none
     double precision    :: aratio,pi,numer,denom
     integer             :: m,n,k,nmax

     pi = 4.0d0*datan(1.0d0)

     dke_u0n = 0.0d0
     ! As a precaution, need to go well beyond m for the denominator (k-1/2)**2 - m**2 to force convergence.
     do k=1,2*nmax
          numer = (-8*(-1)**n*aratio*dtanh(((-0.5 + k)*pi)/aratio))
          denom = ((-0.5 + k)**3*pi**5*((-0.5 + k)**2 + n**2*aratio**2))

          dke_u0n = dke_u0n + numer/denom
     end do

     dke_u0n = dke_u0n

end function dke_u0n
!
! ----------------------------------------------------------
!
double precision function dke_umn(aratio,m,n,nmax)
implicit none
     double precision    :: aratio,pi,numer,denom
     integer             :: m,n,k,nmax

     pi = 4.0d0*datan(1.0d0)

     dke_umn = 0.0d0
     ! As a precaution, need to go well beyond m for the denominator (k-1/2)**2 - m**2 to force convergence.
     do k=1,2*nmax
          numer = (-16*(-1)**(m + n)*aratio*dtanh(((-0.5 + k)*pi)/aratio))
          denom = ((-0.5 + k)*((-0.5 + k)**2 - m**2)*pi**5*((-0.5 + k)**2 + n**2*aratio**2))

          dke_umn = dke_umn + numer/denom

     end do

     dke_umn = dke_umn

end function dke_umn
!
! ----------------------------------------------------------
!
