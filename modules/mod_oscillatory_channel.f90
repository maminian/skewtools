module mod_oscillatory_channel
     ! Contains subroutines to precalculate the 
     ! eigenfunctions cos((j-1/2)*pi*y) on a grid,
     ! and efficiently evaluate the flow at a given time
     ! represented as a matrix-vector product. The resulting
     ! flow at any point can then be found using 
     ! some interpolation scheme.
     !
     
     integer, parameter  :: nefuncs = 100   ! Number of eigenfunctions to use
     integer, parameter  :: ngo = 1001      ! Number of grid points for interpolation

     double precision    :: Pr,omega,amp
     logical             :: skip_transient

     double precision, dimension(ngo,nefuncs)     :: efuncs
     double precision, dimension(ngo)             :: ygrid,mocflow
     double precision, dimension(ngo)             :: term0,term1,term2

     double precision, parameter                  :: pi = 4.0d0*datan(1.0d0)
     double precision, parameter                  :: zero = 0.0d0
     double precision, parameter                  :: one = 1.0d0
     character(len=1), parameter                  :: no = 'n'

     contains
     
     subroutine read_oscillatory_params(param_file)
          ! Reads the oscillatory parameters from the parameter file.
          ! We assume we already know we want to use oscillatory flow.
          !
          implicit none
          character(len=1024), intent(in)    :: param_file
          
          ! Internal
          integer                            :: funit,i
          integer, parameter                 :: nls = 38 ! Number of lines to skip over to get to the good stuff
          character(len=1)                   :: dummy
          
          funit=55
          open(funit,file=param_file)
               do i=1,nls
                    read(funit,*) dummy
               end do
               
               read(funit,*) Pr
               read(funit,*) omega
               read(funit,*) amp
               read(funit,*) skip_transient
          close(funit)
          
     end subroutine read_oscillatory_params

     subroutine precalculate_efuncs()
          implicit none
          ! Does what it says, precalculates the eigenfunctions cos((j-1/2)*pi*y).
          
          integer                                                :: i,j

          double precision, dimension(nefuncs)                   :: fj1,fj2
          double precision lf,jmhp


          do j=1,ngo
               ygrid(j) = -1.0d0 + (j-1)*2.0d0/(ngo-1)
          end do

          do j=1,nefuncs
               efuncs(:,j) = dcos((j-0.5d0)*pi*ygrid)
          end do

          do j=1,nefuncs
               jmhp = pi*(j-0.5d0)
               lf = 4*(-1)**j

               fj1(j) = lf * (-amp*Pr*omega)/(jmhp*(jmhp**4*Pr**2 + omega**2))
               fj2(j) = lf * (jmhp*amp*Pr**2)/(jmhp**4*Pr**2 + omega**2)
          end do

          ! Call the BLAS routine DGEMV for the matrix-vector operations.

          call dgemv(no,ngo,nefuncs,one,efuncs,ngo,fj1,1,zero,term1,1)
          call dgemv(no,ngo,nefuncs,one,efuncs,ngo,fj2,1,zero,term2,1)

     end subroutine precalculate_efuncs
     !
     !----------------------
     !
     subroutine update_flow(tv)
          implicit none

          double precision, intent(in)            :: tv

          integer                                 :: i,j

          double precision                        :: jmhp,lf

          double precision, dimension(nefuncs)    :: fj0
          
          if (.not. skip_transient) then
               do j=1,nefuncs
                    jmhp = pi*(j-0.5d0)
                    lf = 4*(-1)**j
                    
                    fj0(j) = dexp(-jmhp**2*Pr*tv)*lf*(-1.0d0)* &
                              &(jmhp**4*Pr**2 - amp*jmhp**2*Pr*omega + omega**2)/ & 
                              &(jmhp**3*(jmhp**4*Pr**2 + omega**2))
               end do

               call dgemv(no,ngo,nefuncs,one,efuncs,ngo,fj0,1,zero,term0,1)
          else
               call dgemv(no,ngo,nefuncs,zero,efuncs,ngo,fj0,1,zero,term0,1)
          end if

          mocflow = 1.0d0 - ygrid**2 + term0 + term1*dcos(omega*tv) + term2*dsin(omega*tv)

     end subroutine

end module mod_oscillatory_channel
