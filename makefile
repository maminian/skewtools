# Makefile for skewtools
# 
# Last updated: March 2017
# 

# -----
# Subdirectories. Make sure HOME is correct, in particular.
# --------
HOME = $(shell pwd)
UTILS = $(HOME)/utils
COMPS = $(HOME)/computation
MONTE = $(HOME)/monte
TEST = $(HOME)/test
MODS = $(HOME)/modules

# --------------
# Fortran compiler, you _need_ to use h5fc (the hdf5 fortran compiler wrapper), as 
# all the i/o is done with hdf data containers.
# 
# ----------------------

#FC        = /home/katrina/a/aminian/local_installs/hdf5-1.10.1-build-serial/bin/h5fc
FC = /home/katrina/a/aminian/local_installs/hdf5-1.10.1-racecar/install/bin/h5fc

# ---------
# Link your BLAS and LAPACK libraries.
# Currently only BLAS is used, and for pulsatile (oscillatory) flow.
# ----------------

#BLAS = -L/usr/lib64/ -lblas
BLAS = -L/usr/lib64/ /usr/lib64/libblas.so.3
#LAPACK = -llapack

# -------
# Miscellaneous flags. fbackslash allows for 'backspacing' 
# in writes (gfortran only), letting you have a dynamic 
# progress bar (for instance).
# ---------------

OTHER2 = -fbackslash -fbounds-check -funroll-loops -O3

# Uncomment to enable profiling
OTHER = $(OTHER2)
#OTHER = $(OTHER2) -pg


CLEANUP = rm -f ./*.o ./*.mod

# ---------
# END OF CONFIGURABLES
# -----------------------

LINKS   = $(BLAS)

# exe names, object names.
EXES_MC = channel_mc duct_mc ellipse_mc triangle_mc racetrack_mc
EXES_DIRECT = rect_duct rootfinder duct_keff duct_chatwin_keff
EXES_COMPS = dump_uduct_flow dump_utriangle_flow dump_uracetrack_flow channel_moments_exact pipe_skew_exact duct_skew_exact largepe_root dump_uduct_flow_ss
EXES_TEST = cchoose_test rng_tester sortpairs_tester hdf_test1 h5_1d_rw_test lltest mergesort_test median_test

EXES = $(EXES_MC) $(EXES_DIRECT) $(EXES_COMPS) $(EXES_TEST)

rect_duct: $(DIRECT)/rect_duct.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

rootfinder: $(DIRECT)/rootfinder.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
channel_mc: $(MONTE)/channel_mc.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
duct_mc: $(MONTE)/duct_mc.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
ellipse_mc: $(MONTE)/ellipse_mc.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

triangle_mc: $(MONTE)/triangle_mc.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

racetrack_mc: $(MONTE)/racetrack_mc.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

cchoose_test: $(TEST)/cchoose_test.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
rng_tester: $(TEST)/rng_tester.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

mt_tester: $(TEST)/mt_tester.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90  -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

sortpairs_tester: $(TEST)/sortpairs_tester.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_uduct_flow: $(TEST)/dump_uduct_flow.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_uduct_flow_ss: $(TEST)/dump_uduct_flow_ss.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_utriangle_flow: $(TEST)/dump_utriangle_flow.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_uracetrack_flow: $(TEST)/dump_uracetrack_flow.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

hdf_test1: $(TEST)/hdf_test1.f90
	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

#triangle_mc: $(MONTE)/triangle_mc.f90
#	$(FC) $(LINKS) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
#	$(CLEANUP)

	
channel_moments_exact: $(DIRECT)/channel_moments_exact.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

pipe_skew_exact: $(DIRECT)/pipe_skew_exact.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_skew_exact: $(DIRECT)/duct_skew_exact.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_mc_alt: $(MONTE)/duct_mc_alt.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

largepe_root: $(TEST)/largepe_root.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 $(COMPS)/*.f90  -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_keff: $(DIRECT)/duct_keff.f90
	$(FC) $(LINKS) $(DIRECT)/*.f90 $(COMPS)/*.f90 $(MODS)/*.f90 $(UTILS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_chatwin_keff: $(DIRECT)/duct_chatwin_keff.f90
	$(FC) $(LINKS) $(UTILS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

h5_1d_rw_test: $(TEST)/h5_1d_rw_test.f90
	$(FC) $(LINKS) $(UTILS)/hdf_create_file.f90 $(UTILS)/hdf_add_1d_darray_to_file.f90 $(UTILS)/almod.f90 $(UTILS)/hdf_read_1d_darray.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

lltest: $(TEST)/lltest.f90
	$(FC) $(LINKS) $(MODS)/mtfort90.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

mergesort_test: $(TEST)/mergesort_test.f90
	$(FC) $(LINKS) $(MODS)/mtfort90.f90 $(UTILS)/my_normal_rng.f90 $(UTILS)/mergesort.f90 $(TEST)/$@.f90 -o $@ $(LINKS) $(OTHER)
	$(CLEANUP)

median_test: $(TEST)/median_test.f90
	$(FC) $(LINKS) $(MODS)/mtfort90.f90 $(UTILS)/my_normal_rng.f90 $(UTILS)/mergesort.f90 $(COMPS)/median.f90 $(TEST)/$@.f90 -o $@ $(LINKS) $(OTHER)
	$(CLEANUP)

clean:
	rm -f $(EXES) ./*.o ./*.mod


