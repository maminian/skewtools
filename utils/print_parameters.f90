subroutine print_parameters(aratio,q,Pe,nGates,nTot,y0,z0,save_hist,&
                              t_hist,dtmax,nt,ntt,mt_seed,geometry,use_external_ic)

! Spit out a nice header with all the simulation settings.

implicit none
     integer, parameter                                :: i64 = selected_int_kind(18)
     
     double precision, intent(in)                      :: aratio,q,Pe,y0,z0,dtmax
     integer, intent(in)                               :: nGates,nTot,nt,ntt
     logical, intent(in)                               :: save_hist,use_external_ic
     double precision, dimension(1:nt), intent(in)     :: t_hist
     integer(kind=i64), intent(in)                     :: mt_seed
     character(len=1024), intent(in)                   :: geometry
     character(len=15)                                 :: num2str

     character(len=1024), parameter                    :: racetrack = "racetrack"
     character(len=1024), parameter                    :: channel = "channel"

     write(*,*) ""
     write(*,"(80A)") REPEAT("=",80)
     write(*,*) ""
     
     write(*,"(A10,A7)") "Geometry: ",trim(geometry)
     write(*,*) ""

     ! Display differently for point source and uniform.
     if (use_external_ic) then
          write(*,"(A28)") "Initial data read from file."
          write(*,*) ""
     else     
          if (nGates .eq. 1) then
               write(*,"(A50)") "Point source initial data."
               write(*,*) ""
               
               ! Display point source different for channel and pipe/duct.
               if (geometry .eq. "channel") then
                    write(*,"(A50,ES10.3)") "Starting coordinate: ",y0
               else
                    write(*,"(A50,A1,ES10.3,A3,ES10.3,A2)") "Starting coordinates: ","(",y0," , ",z0," )"
               end if

          else
               write(*,"(A50)") "Uniform initial data."
          end if
     end if

     write(*,*) ""

     write(*,"(A50,ES10.3)") "Peclet: ", Pe
     if (.not. (geometry .eq. channel)) then
          write(*,"(A50,ES10.3)") "Aspect ratio: ",aratio
          write(*,*) ""
     end if
     if (geometry .eq. racetrack) then
          write(*,"(A50,ES10.3)") "Shape parameter: ",q
     end if
     write(*,*) ""
     
     write(num2str,"(I15)") nTot
     write(*,"(A50,A15)") "Number of particles: ", adjustl(num2str)
     
     write(*,"(A50,A1,ES10.3,A3,ES10.3,A2)") "Time interval: ","(", t_hist(1), " , ", t_hist(nt)," )"
     
     write(num2str,"(I15)") ntt
     write(*,"(A50,A15)") "Number of requested timesteps: ",adjustl(num2str)
     
     write(num2str,"(I15)") nt
     write(*,"(A50,A15)") "Number of internal timesteps: ",adjustl(num2str)
     write(*,"(A50,ES10.3)") "Largest internal timestep: ",dtmax
     
     write(*,*) ""

     if (save_hist) then
          write(*,"(A40)") "Saving position histories to file."
          write(*,*) ""
     end if
     
     write(num2str,"(I15)") mt_seed
     write(*,"(A40,A15)") "Mersenne Twister seed: ", adjustl(num2str)

1234 continue
     
     write(*,*) ""
     write(*,"(80A)") REPEAT("=",80)
     write(*,*) ""
     

     
end subroutine print_parameters
