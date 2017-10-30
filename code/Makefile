FF = gfortran
FFLAGS = -Wunused -ffree-form -ffree-line-length-none -fimplicit-none -fbounds-check -g3 -fbacktrace
LIBS = -llapack -lblas

FILES = det_n2___.f _rndm_mumbers.f
MAINFILE = mf_cmd_phi_gf_oct21_2013_2D.f
EXEFNAME = hub_gf_2d_run_2

obc: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice_obc.f greenfun_pbc_ED.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)

ed: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice.f greenfun_pbc_ED.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)

pbc: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice.f greenfun_pbc.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)

clean:
	#rm *.o
	rm *.mod
	ls -l $(EXEFNAME)
