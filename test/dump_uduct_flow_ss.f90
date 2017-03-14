program dump_uduct_flow_ss
! Takes input parameters of the mesh size, number of terms
! in Fourier series, and aspect ratio, and dumps the arrays 
! of (y,z,u) into an h5 file
! to be read/visualized in Python (or otherwise).

implicit none

     ! Function arguments
     integer             :: nMesh,nTerms
     double precision    :: aratio

     ! Fourier coefficients and index pairs
     integer, dimension(:,:), allocatable              :: idx_list
     double precision, dimension(:), allocatable       :: uij_vals
     
     ! Internal
     integer                                           :: m,n
     double precision                                  :: hy,hz
     character(len=1024)                               :: temp
     
     ! Output
     double precision, dimension(:,:), allocatable     :: u_out,u_lap_out
     double precision, dimension(:), allocatable       :: y,z
     integer                                           :: funit
     character(len=1024)                               :: fname
     
     
     ! hdf version
     character(len=1024)                               :: fname2,arrayname
     character(len=1024)                               :: description

     ! Functions
     double precision u_duct_ss


     ! -------------------------------------------
     
     call get_command_argument(1,temp)
     read(temp,*) nMesh
     call get_command_argument(2,temp)
     read(temp,*) nTerms
     call get_command_argument(3,temp)
     read(temp,*) aratio
     
     allocate(u_out(nMesh,nMesh))
     allocate(y(nMesh),z(nMesh))
     
     allocate(u_lap_out(nMesh-2,nMesh-2))

     ! Calculate flow.

     ! call coeff_chooser(nTerms,aratio,idx_list,uij_vals)

     hy = 2.0d0/dble(nMesh-1)
     hz = (2.0d0/aratio)/dble(nMesh-1)
     do m=0,nMesh-1
          y(m+1) = -1.0d0+hy*m
          z(m+1) = -1.0d0/aratio+hz*m
     end do
     
     do m=1,nMesh
          do n=1,nMesh
               u_out(m,n) = u_duct_ss(y(m),z(n),nTerms,aratio)
          end do
     end do
     
     ! Calculate the Laplacian on the interior.
     do m=2,nMesh-1
          do n=2,nMesh-1
               u_lap_out(m-1,n-1) = 1.0d0/(hy**2)*(u_out(m-1,n) -2.0d0*u_out(m,n) + u_out(m+1,n)) + &
                                   1.0d0/(hz**2)*(u_out(m,n-1) - 2.0d0*u_out(m,n) + u_out(m,n+1))

!               u_lap_out(m-1,n-1) = 1.0d0/(h**2) * u_lap_out(m-1,n-1)

          end do
     end do

     ! ------------------
     ! Write output

     fname = trim("u_duct.h5")
     
     arrayname = "Flow_values"
     description = "These are flow values for the duct. Soon to be added, including params."
     call hdf_create_file(fname)
     call hdf_add_2d_darray_to_file(nMesh,nMesh,u_out,fname,arrayname,description)
     
     arrayname = "Laplacian_values"
     description = "These are values of the 5-point (scaled) discrete Laplacian on the interior."
     call hdf_add_2d_darray_to_file(nMesh-2,nMesh-2,u_lap_out,fname,arrayname,description)
     
     arrayname = "y_mesh"
     description = "The y mesh used in calculating the flow."
     call hdf_add_1d_darray_to_file(nMesh,y,fname,arrayname,description)
     
     arrayname = "z_mesh"
     description = "The z mesh used in calculating the flow."
     call hdf_add_1d_darray_to_file(nMesh,z,fname,arrayname,description)

     arrayname = "aratio"
     description = "Aspect ratio of the domain."
     call hdf_add_1d_darray_to_file(1,aratio,fname,arrayname,description)

     deallocate(u_out,y,z,u_lap_out)
     
end program dump_uduct_flow_ss
