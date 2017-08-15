module mod_parameters
     ! Contains definitions for all the parameters, 
     ! in addition to some hard-coded values.
     
     integer, parameter                 :: i64 = selected_int_kind(18)

     ! Length of the buffer before writing to disk.
     ! Only relevant if saving the entire position history.
     ! Run time is bottlenecked by read/writes to file to a severe degree, so in general 
     ! this should be made as large as possible while still fitting in memory.
     !
     ! Simulation memory with no buffer is
     ! approximately the size of a double (8 bytes) times 2 
     ! or 3 X,Y(,Z) times the number of particles. For example, with 1e6 particles, 
     ! approximate memory usage for channel simulation is 16e6 bytes, or 15.26MB.
     ! With a buffer length of 100, this becomes 1526MB, or 1.49GB.
     integer, parameter                 :: buffer_len = 100

     ! Number of bins when looking at the cross-sectionally averaged distribution.
     integer, parameter                 :: nhb = 400

     ! Scale of the channel: [-a,a]. Not much reason to use anything other than 1.0d0.
     double precision, parameter        :: a = 1.0d0
     double precision                   :: b
     character(len=1024)                :: geometry
     
     ! These will be read from file.
     double precision                   :: aratio, q, Pe, x0width, y0, z0, t_warmup, &
                                             dt, dtmax, Tfinal
     integer                            :: nGates, x0n, nbins
     logical                            :: save_hist, save_hist2d, use_external_ic, use_oscillatory
     character(len=1024)                :: ic_file, tstep_type, other_file
     
     integer(i64)                       :: mt_seed

end module mod_parameters
