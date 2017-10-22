FF = gfortran
FFLAGS = -ffree-form -ffree-line-length-none -fimplicit-none -fbounds-check -g3 -fbacktrace
LIBS = -llapack -lblas

FILES = det_n2___.f _rndm_mumbers.f lattice.f mf_cmd_phi_gf_oct21_2013_2D.f

pbc: $(FILES)
	$(FF) $(FFLAGS) $(FILES) -o hub_gf_2d_run_1 $(LIBS)

clean:
	#rm *.o
	rm *.mod
	ls -l hub_gf_2d_run
