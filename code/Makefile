FF = gfortran
FFLAGS = -Wunused -ffree-form -ffree-line-length-none -fimplicit-none -fbounds-check -g3 -fbacktrace
LIBS = -llapack -lblas

FILES = det_n2___.f _rndm_mumbers.f
MAINFILE = mf_cmd_phi_gf_oct21_2013_2D.f

circle: EXEFNAME = hub_gf_2d_run_circle
circle: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice_circle.f greenfun_pbc_ED.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)
	ls -l $(EXEFNAME)

obc: EXEFNAME = hub_gf_2d_run_obc
obc: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice_obc.f greenfun_pbc_ED.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)
	ls -l $(EXEFNAME)

ed: EXEFNAME = hub_gf_2d_run_pbc_ed
ed: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice.f greenfun_pbc_ED.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)
	ls -l $(EXEFNAME)

pbc: EXEFNAME = hub_gf_2d_run_pbc 
pbc: $(FILES)
	$(FF) $(FFLAGS) $(FILES) lattice.f greenfun_pbc.f $(MAINFILE) -o $(EXEFNAME) $(LIBS)
	ls -l $(EXEFNAME)

clean:
	#rm *.o
	rm *.mod
	ls -l $(EXEFNAME)
