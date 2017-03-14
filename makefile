# Makefile for geometric skewness code.
#
# Manuchehr Aminian, September 2014


# -----
# Subdirectories. Make sure HOME is correct, in particular.
# --------
HOME = $(shell pwd)
UTILS = $(HOME)/utils
COMPS = $(HOME)/computation
DIRECT = $(HOME)/direct
MONTE = $(HOME)/monte
TEST = $(HOME)/test
MODS = $(HOME)/modules

# --------------
# Fortran compiler, you _need_ to use h5fc (the hdf5 gfortran wrapper) if you want to use hdf5.
# (You'll need to edit looper.f90 and remove the call to dSaveArrayToHDF otherwise.)
#
# If you modify the code appropriately (remove calls to dSaveArrayToHDF), gfortran should work fine.
# ----------------------

FC        = h5fc
#FC        = gfortran

# ---------
# Link your BLAS library
# ----------------
BLAS = -lblas

# --------
# Link your LAPACK library (necessary if you plan to use it for the eigenvalue problem)
# ----------------
LAPACK = -llapack

# -------
# Miscellaneous flags. fbackslash allows for 'backspacing' in writes,
# letting you have a dynamic progress bar (for instance).
# ---------------

OTHER2 = -fbackslash -fbounds-check -funroll-loops -O3


# -fbounds-check -O3 -funroll-loops -Wall

# Uncomment to enable profiling
#OTHER = $(OTHER2) -pg
OTHER = $(OTHER2)

CLEANUP = rm -f ./*.o ./*.mod

# ---------
# END OF CONFIGURABLES
# -----------------------

LINKS   = $(BLAS)

# exe names, object names.
EXES = rect_duct rootfinder cchoose_test rng_tester sortpairs_tester duct_mc channel_mc dump_uduct_flow dump_utriangle_flow hdf_test1 ellipse_mc triangle_mc racetrack_mc channel_moments_exact pipe_skew_exact duct_skew_exact duct_mc_alt largepe_root dump_uduct_flow_ss duct_keff triangle_mc duct_chatwin_keff h5_1d_rw_test lltest

rect_duct: $(DIRECT)/rect_duct.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

rootfinder: $(DIRECT)/rootfinder.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
channel_mc: $(MONTE)/channel_mc.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
duct_mc: $(MONTE)/duct_mc.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
ellipse_mc: $(MONTE)/ellipse_mc.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

triangle_mc: $(MONTE)/triangle_mc.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

racetrack_mc: $(MONTE)/racetrack_mc.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

cchoose_test: $(TEST)/cchoose_test.f90
	$(FC) $(UTILS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
rng_tester: $(TEST)/rng_tester.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

mt_tester: $(TEST)/mt_tester.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90  -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

sortpairs_tester: $(TEST)/sortpairs_tester.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_uduct_flow: $(TEST)/dump_uduct_flow.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_uduct_flow_ss: $(TEST)/dump_uduct_flow_ss.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

dump_utriangle_flow: $(TEST)/dump_utriangle_flow.f90
	$(FC) $(MODS)/*.f90 $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

hdf_test1: $(TEST)/hdf_test1.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)
	
channel_moments_exact: $(DIRECT)/channel_moments_exact.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

pipe_skew_exact: $(DIRECT)/pipe_skew_exact.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_skew_exact: $(DIRECT)/duct_skew_exact.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_mc_alt: $(MONTE)/duct_mc_alt.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90 -o $@ $(MONTE)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

largepe_root: $(TEST)/largepe_root.f90
	$(FC) $(UTILS)/*.f90 $(COMPS)/*.f90  -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_keff: $(DIRECT)/duct_keff.f90
	$(FC) $(DIRECT)/*.f90 $(COMPS)/*.f90 $(MODS)/*.f90 $(UTILS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

duct_chatwin_keff: $(DIRECT)/duct_chatwin_keff.f90
	$(FC) $(UTILS)/*.f90 -o $@ $(DIRECT)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

h5_1d_rw_test: $(TEST)/h5_1d_rw_test.f90
	$(FC) $(UTILS)/hdf_create_file.f90 $(UTILS)/hdf_add_1d_darray_to_file.f90 $(UTILS)/almod.f90 $(UTILS)/hdf_read_1d_darray.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

lltest: $(TEST)/lltest.f90
	$(FC) $(MODS)/mtfort90.f90 -o $@ $(TEST)/$@.f90 $(LINKS) $(OTHER)
	$(CLEANUP)

clean:
	rm -f $(EXES) ./*.o ./*.mod


