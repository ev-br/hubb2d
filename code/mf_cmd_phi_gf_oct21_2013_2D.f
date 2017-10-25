!
! Determinant diagrammatic MC code for the Fermi-Hubbard model
! at finite-temperature.
!
!
!  for details see
! 
!  http://montecarlo.csi.cuny.edu/umass/fermion.html
!
! **********************************
!    This is a non-MPI version 
! **********************************
!
!
! *****************   This version [mf_cmd_phi_gf_oct21_2013_2d.f] is derived from mf_cmd_phi__dec2_2010_2D.f 
!                     by the addition of measurment of the full green's function.   
!
! *********** The difference w/respect to __nov24_2010.f version is that in ene(...) phi/L is used, not bare \phi
!
!    Need: 2. profile -- measurements etc
!          3. time --- rewrite the way m_u(:) etc are filled ?
!          4. BC twist
! 
!
!
!

!---------------------
!--- Variables
!---------------------
	module vrbls
	implicit none

	real*8,  parameter :: nul=0.0d0, un1=1.0d0
	real*8   :: pi
	integer, parameter :: i_hu = huge(1)

	integer, parameter :: nkink_m = 512  ! max # of kinks on one site
	integer, parameter :: nmnm_max = 8192 ! max total # of kinks

!------------------------------------------------------
!     Phys. parameters
!------------------------------------------------------
	real*8 :: U, U_in    ! MINUS interaction strength, initial for thermalization
	real*8 :: beta       ! inverse temperature
	real*8 :: beta_ini   ! start thermalization with this \beta
	real*8 :: beta_fin
	real*8 :: mu         ! chem. potential
	real*8 :: ef, kf     ! Fermi energy & momentum

    real*8 :: eta        ! GF- vs Z-sector weight

	real*8,parameter :: gf_eta=1.038d0   ! 1 + GF anomalous dimension, U(1) universality class

	real*8,parameter :: H_eta=0.03d0   ! eta in Heisenberg universality class

    real*8 :: alpha         ! Rubtsov's alpha


!------------------------------------------------------
!     Update related
!------------------------------------------------------

	real*8 :: bun           !  = beta*U*Nsite 

!------------------------------------------------------
!     Configuration
!------------------------------------------------------

! The data structure is purposedly redundant, for in all cases where 
! there is a memory/speed tradeoff, it is solved in favor of speed.
!
! The interaction vertices (a.k.a 'kinks')  live on the sites of a d-dimensional lattice
! and on the [0,\beta] interval of the \tau-continuum. The kinks are labelled by
! the 'names', e.g. unique numbers. The configuration consists of:
!
! 1. Overall list of kinks' names, managed by the name manager, see GetName() and DropName() subroutines
!    Given the kink name, one has access to the kink's site, temporal position (a.k.a 'time', also 'tau'),
!    and the position in TheMatrix, see below. 
!
! 2. Each site hosts a list of names of kinks which belong to this site. The order of kinks in these
!    lists is arbitrary.
!
! 3. The configuration weight matrix is not stored itself; one always deals with its inverse, 
!    a.k.a. TheMatrix :). 
!
! 4. ira/masha are not kinks. They are stored separately, see below.
!
	integer,   allocatable :: kname(:,:)   ! list of kinks's names @ a given site: kname(j,site)
	integer*2, allocatable :: nkink(:)     ! # of kinks @ site: nkink(site)

	integer, allocatable   :: ksite(:)      ! ksite(name) => site of a kink 'name'
	real*8, allocatable    :: ktau(:)       ! ktau(name) => tau of a kink 'name'

	integer, allocatable   :: row(:), clmn(:) ! row(name), clmn(name) => position of a 'name' in TheMatrix
	integer, allocatable   :: nm_row(:),nm_clmn(:)  ! nm_row(row) => name of the kink associated with the row

! NameManager data
	integer                        :: nmnm         ! total number of kinks
      INTEGER, DIMENSION(-2:nmnm_max) :: namelist, numnam

      logical :: present             ! =.true. is ira & masha are present in the configuration
      integer, parameter :: ira=-2, masha=-1

! TheMatrix
	integer             :: pm,lda      ! actual size & leading dimension
	real*8, allocatable :: matr(:,:)   ! TheMatrix, inverse
	real*8, allocatable :: m_u(:), m_v(:), m_w(:),m_z(:)       ! working arrays for rank-1 updates
	real*8              :: m_s

!	real*8, allocatable :: m_u2(:,:),m_v2(:,:),m_s2(:,:),m_c2(:,:)  ! for rank-2 updates
	real*8, allocatable :: m_u2(:,:),m_v2(:,:),m_s2(:,:),m_z2(:,:),m_w2(:,:),m_x2(:,:)

!-------------------------------------------------------
!     Measurement related
!-------------------------------------------------------

!
!   Errorbar analysis is performed via blocking method. The maximum number of blocks is
!   fixed, once this number is reached, the blocks are collated so that one has twice
!   as few blocks of twice the size. 
!
	integer, parameter :: b_n_max=100   ! max # of blocks
	integer  :: Z_b                     ! block size

	integer :: b_n, i_b    ! # of filled blocks  & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block
! Partition function
	real*8 :: Z            

! Diagonal estimators:
	real*8  :: PE, KE, ndens     ! PE, KE and density estimators 
	real*8  :: PE_stat(b_n_max), KE_stat(b_n_max), ndens_stat(b_n_max) ! files

	real*8  :: PE_2, PEKE
	real*8  :: PE_2_stat(b_n_max), PEKE_stat(b_n_max)

	!! For 3D specifically:
!	real*8  :: g_uu(0:20,0:20,0:20), g_ud(0:20,0:20,0:20)  ! Density-densiy correltaors, a.k.a 'dance-dance'
!	real*8  :: g_uu0 !! the current contribution to g_uu(0)
!
!	real*8  :: z_uu(0:20,0:20,0:20),z_ud(0:20,0:20,0:20)


	!!  For 1D ONLY
!	real*8  :: g_uu(0:20), g_ud(0:20)  ! Density-densiy correltaors, a.k.a 'dance-dance'
	real*8  :: g_uu0 !! the current contribution to g_uu(0)
!	real*8  :: z_uu(0:20),z_ud(0:20)

	real*8, allocatable :: g_uu(:,:), z_uu(:,:)
	real*8, allocatable :: g_ud(:,:), z_ud(:,:)

!! the full Green's function
	real*8, allocatable :: gf_r(:,:,:), gf_i(:,:,:), zgf(:,:,:)
	integer :: G_nfreq=64 !! number of frequencies for the Green's function

! Integrated ira-masha correlator:
	real*8  :: im               ! estimator
	real*8  :: im_stat(b_n_max) ! file

! Integrated spin-spin correlator:
	real*8 :: gss, gss0
	real*8 :: gss_stat(b_n_max)
 
! <Sz:Sz>@nn sites
        real*8 :: nn_szsz
        real*8 :: nn_szsz_stat(b_n_max)

! Service histograms:
	integer, parameter :: im_t_max=15
	real*8  :: im_t(-im_t_max:im_t_max)          ! t(ira-masha) histogram

	integer, parameter :: nmnm_sup = 7500        ! kink # distribution
	real*8, dimension(0:nmnm_sup) :: nmnm_distr

	integer, parameter :: det_distr_max=25       ! det. distribution
	real*8 :: det_distr(-det_distr_max:det_distr_max)


!------------------------------------------------------
!      Flow control and output
!------------------------------------------------------

!
!  Since configuration weight determinant ratios (a.k.a.'det') are
!  recalculated incrementally, the roundoff error buildup is possible
!  Also, if the current configuration is close to the nodal surface,
!  the det-s might get extravagantly large/small, which also might lead to
!  error acumulation. In order to deal with these, determinants are recalculated 
!  from scratch (which costs N**3 operations!) (i) every so often (once per N**4 
!  proves to be OK), and (ii) every time unusually large det pops out (the latter
!  signals that the pre-update configuration was close to a nodal surface)
!
!
	real*8 :: tolerance    ! determinant recalculation tolerance level
	real*8 :: recalc_rate  ! recalc. rate counter

! Update manager probabilities
	real*8, dimension(1:10) :: prob

      real*8 :: step, step_p, step_w, step_m, step_t, step_r
!     step    counts MC steps
!     step_p  the number of steps between printouts
!     step_w  the number of steps between writing to disk
!     step_m  the number of steps between measurements
!     step_t  the number of steps for thermolization
!     step_r  the number of steps for det. recalculation

!  Counters for printing, writing, measuring, thermolization, det recalculation
      real*8 :: i_p, i_w, i_m, i_t, i_r

      real*8  :: n_sw   ! number of sweeps for thermalization, 1 sweep = \beta*U*Nsite updates

! Address/accept counters	for various updates
	real*8 :: c_a_v, a_a_v, c_d_v, a_d_v
	real*8 :: c_a_v2, a_a_v2, a_d_v2, c_d_v2
	real*8 :: c_a_v2a, a_a_v2a, a_d_v2a, c_d_v2a

! Kink number output
      integer :: nm_min=0, nm_max=i_hu
      real*8  :: nm_av=0

! 'status' variable holds the name of the update; proves useful for debugging
	character*6 :: status, prevstatus

! Output channel #:
	integer, parameter :: OUT_UNIT=22        ! For MPI code, writing to the stdout is useless, 
											 ! thus all the output goes to file(s)

! rndm() seed, see the 'mumbers' module
	integer :: r_ij, r_kl
	integer, allocatable :: ar_ij(:), ar_kl(:)


! Walltime constraints
	character*8  :: date               ! "how much on your watch?" -- "Six o'cloch"
	character*10 :: time
	character*5  :: zone

	integer :: tvalues(8)
      
	real*8 :: time_prev, time_curr, time_elapsed, time_limit
	integer :: hr_curr, hr_prev

	character*50 :: myparfname, fnamesuffix, outputfname
	character*50 :: str,istr, anfname, rndmstr                     ! service vrbls


!  Working protocol: First, the file named 'infile' is read. The format is:
!
!    nfls
!    nc_1      par_XXX
!    nc_2      par_YYY
!    ........
!    nc_nfls   par_ZZZ
!
!  Here nfls is the number of different parameter files, and nc_i is the number
!  of clones to run per i-th parameter set. The latter number must match the 
!  corresponding entry in the parameter file.
!
!  For a parameter file par_XXX, _XXX is used as a suffix for all the output:
!    stat_XXX_#.dat, out_XXX_#.out etc, where # is the number of a clone (1...nc_i).
!
!  Parameter file also contains the walltime limit in hours. System time is
!  checked in the printouts, and as it exceeds the limit, the process 
!  saves and wraps up.

	end module vrbls




!===============================================
!=== Main block :)
!===============================================
	program MAIN
	use vrbls; use det_n2
	use mumbers; !use tree
    use lattice
    use green_functions
	implicit none 

	logical :: acpt                   ! accepted update flag
	integer :: i, j, ans1, ans2, ans3
	real*8  :: r, dp, det , dt
    character*6 :: cvoid


! ------- pi pi pi pi ------------
	pi = 4.d0*atan(1.d0)

!--------- Set the clock --------
	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	time_prev=0.d0; time_elapsed = 0.d0; hr_prev = 0


!--------- Get the par-file name from the command line --------------	
	call  GETARG(1 , myparfname)	
	
	! determine the suffix for output files
	fnamesuffix = myparfname(4:)
  
	outputfname = 'out'//trim(fnamesuffix)   ! send output here


!------- Reading parameter file ------------------

	open(OUT_UNIT,file=outputfname,position='append')      ! to be open until the MC starts
	write(OUT_UNIT,*) ' reading ', myparfname,' ...'

    open( 1, file=trim(myparfname) )

    read(1,*)d
    dd = 2*d

    allocate( N(1:d), N2(1:d) )
    Nsite=1
    do i=1,d 
        read(1,*) N(i)  ! Number of sites per given dimension
        Nsite=Nsite*N(i)
        N2(i)=N(i)/2
    end do

	read(1,*) ans1  ! 0 if new configuration, 1 if old one
	read(1,*) ans2  ! 0 if new statistics,    1 if old one
	read(1,*) ans3  ! 0 if new rndm() seed, 1 if read one
	
	read(1,*) mu            ! Chemical potential

	read(1,*) U,U_in        ! - U, interaction, & initial for thermalization
	read(1,*) beta, beta_ini     ! Inverse temperature
	    bun = beta*U*Nsite  ! Shorthand
	    step_r = bun*bun    ! step for determinant recalculation
	    beta_fin = beta; 
		
	read(1,*) eta          ! G- vs. Z-sector weight
	read(1,*) mtau         ! # of points in \tau for tabulation of GFs
						   ! Again, for inlining reasons it's necessary to have static variables here
	   if(mtau>mt_max)then; print*,'mtau > mtau_max'; call mystop; 
	   endif               
	   bmt = beta/mtau; bmt1 = 1.d0/bmt                   ! shorthands
	read(1,*) tolerance    ! det recalculation tolerance
      read(1,*) n_sw                 ! number of sweeps for thermolization
      read(1,*) step_p               ! step for printing
      read(1,*) step_w               ! step for writing to disk
      read(1,*) step_m               ! step for measuring
	read(1,*) time_limit           ! time limit [hrs]
 	    time_limit = time_limit*3600   ! seconds      
	read(1,*) cvoid                ! horizontal separator :)
! Update manager probabilities: 

      read(1,*) prob(1); prob(2)=prob(1)+prob(1)                ! add/drop
      read(1,*) dp; prob(3)=prob(2)+dp; prob(4)=prob(3)+dp      ! add/drop 2 @ same site
      read(1,*) dp; prob(5)=prob(4)+dp; prob(6)=prob(5)+dp      ! add/drop 2 @ same site
	if(abs(prob(6)-1.d0)>1.0d-10)then
	   print*,'Update manager probs. normalization problem...'
	   call mystop
	endif

!--- Random number generator:
!
!  Number of seeds provided MUST be equal to the groupsize (one seed per process)
!  Otherwise the program crashes with mysterious error message
!
!
	if(ans3==0)then     ! init rndm() state

	  read(1,*)cvoid

	  
	  read(1,*) r_ij,r_kl
	  call init_rndm(r_ij,r_kl)
	  
	else                ! read rndm() state

	  rndmstr='rndm'//trim(fnamesuffix)//'.dat'
	  call read_rndm(rndmstr)

	endif

      ! BC twist
      allocate(phi(d))
      read(1,*)phi(1)
	  read(1,*)phi(2)

	close(1)

	write(OUT_UNIT,*)'....done'


!---- Shift chem. potential w/ Rubtsov's alpha @@@alpha
	alpha = 0.5d0          ! assume 1/2 fillin'
	mu = mu + alpha*U
	print*,' shift \mu to ',mu



!---------- DONE reading, allocating memory:

!--- Lattice
	allocate (ass(dd,Nsite),back(dd),x(1:d,1:Nsite))
	CALL ASSA

!--- Configuration
	allocate( nkink(1:Nsite), kname(nkink_m,1:Nsite) )  ! site data
	allocate( ksite(-2:nmnm_max),ktau(-2:nmnm_max) )    ! global lists

	allocate( row(-2:nmnm_max),clmn(-2:nmnm_max) )      ! c/r <-> name links
	
	lda=128; allocate(nm_row(lda),nm_clmn(lda))         ! Allocate 128*128 matrix first. 
										! more memory will be allocated later if necessary	
	allocate(matr(1:lda,1:lda),m_u(1:lda),m_v(1:lda))   ! TheMatrix
	allocate(m_w(1:lda),m_z(1:lda))
	allocate(m_u2(lda,2),m_v2(lda,2),m_z2(lda,2),m_w2(lda,2),m_x2(lda,2),m_s2(2,2)) 

!	allocate( g_uu(  ), z_uu( 0:N(1)-1, 0:N(2)-1, 0:N(3)-1 ) )
!	allocate( g_ud( 0:N(1)-1, 0:N(2)-1, 0:N(3)-1 ), z_ud( 0:N(1)-1),0:N(2)-1, 0:N(3)-1 ) )
        allocate( g_uu(0:N(1)-1, 0:N(2)-1), z_uu(0:N(1)-1, 0:N(2)-1) )
        allocate( g_ud(0:N(1)-1, 0:N(2)-1), z_ud(0:N(1)-1, 0:N(2)-1) )
		
        allocate( gf_r(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1) )
		allocate( gf_i(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1) )
		allocate( zgf(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1))
	

! Name Manager is initialized in either init_cnf or rd_cnf


!------------ DONE with memory allocation, proceed to tabulations

!--- Tabulate Green functions
	write(OUT_UNIT,*)'tabulating GFs...'  
	if(mtau>mt_max)then; write(OUT_UNIT,*)'mt_max is too small'; 
	call mystop; endif

	print*,' regular tabulate'
	call TABULATE(beta, mtau, mu)
	write(OUT_UNIT,*)'...done'


!--- Structure output
      nm_av=0; nm_min=i_hu; nm_max=0
      nmnm_distr=0; det_distr=0.d0; 


! 'how much on your watch?' a.k.a initialize the clock
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! take care of the runs that span across midnight 


	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' init. done: ',dt,' sec'	
	
!--- close output
	close(OUT_UNIT)


!----------- DONE with initializations


!=========== Read / init conf & thermalize =========


	    ! 0 -- init&therm, 1-- rd, 2-- rd& retherm [w/ tuning \beta]
	select case(ans1)
        case(2)                !read & rethermalize w/ tuning \beta
	     call rd_cnf
	     call therm1
	case(1)                !read & rethermalize w/ fixed \beta
	     call rd_cnf
	     call therm2
	case default
	     call init_cnf
	     call therm1
	end select

! nullify counters
	i_p=0.d0; i_w=0.d0; i_m=0.d0; step=0.d0; 
        c_a_v=0.d0; a_a_v=0.d0; c_d_v=0.d0; a_d_v=0.d0
        c_a_v2=0.d0; a_a_v2=0.d0; c_d_v2=0.d0; a_d_v2=0.d0
        c_a_v2a=0.d0; a_a_v2a=0.d0; c_d_v2a=0.d0; a_d_v2a=0.d0
	recalc_rate=0.d0

! Construct TheMatrix
	det=recalc_matrix(pm,lda,matr)


!============== Read / init  statistics ==============
      if (ans2 == 1) then
	 call init_stat
         call RD_STAT
      else
         call INIT_STAT
      end if


!==================== main MC loop  ===============
	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'  '
	write(OUT_UNIT,*)' Start MC loop'
	close(OUT_UNIT)

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_m=i_m+un1; i_r=i_r+1.d0

!------------------------ update --------------------------
	   r = rndm()   ; 
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call add_2_same(acpt,det)
	   else if(r<prob(4))then;    call drop_2_same(acpt,det)	
	   else if(r<prob(5))then;    call add_2_any(acpt,det)
	   else if(r<prob(6))then;    call drop_2_any(acpt,det)	
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

! det. distribution
            i = floor(log10(abs(det)))
	      if(abs(i)<det_distr_max) det_distr(i) = det_distr(i)+1.d0

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-----------------------------------------------------

         nm_av=nm_av+(un1*nmnm)/step_p        
         IF(nmnm>nm_max)nm_max=nmnm
         IF(nmnm<nm_min)nm_min=nmnm

	   if(nmnm<nmnm_sup) nmnm_distr(nmnm) = nmnm_distr(nmnm) + 1.d0


         if (i_m == step_m)  then; i_m=nul; call MEASURE; end if
         if (i_p == step_p)  then; i_p=nul; call PRNT;
	       end if  
         if (i_w == step_w)  then;  i_w=nul; call WRT;    
	   end if 

	enddo


	contains  

!*************************  subroutines below ONLY ************************

!-----------------
!--- Add a kink
!-----------------
	subroutine add(acpt,det)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio
	integer :: name, site,j,nk, vova, sv
	real*8  :: tnew, tv

	c_a_v = c_a_v + un1; acpt=.false.
	prevstatus=status; status='_add__' 

	if(pm==lda)then; call resize_matrix(-1); endif

!------- check
	if(nmnm>=nmnm_max)then;
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif
!---------- end check

	site=Nsite*rndm()+1.d0
	if(site>Nsite)site=Nsite    ! select a site to play on
	tnew = beta*rndm()                    ! select tau		

!------------- determinant
	if(pm==0)then; det=g0000-alpha  !ratio = g0000**2
	else

! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
       sv = ksite(vova); tv=ktau(vova)
	   m_v(clmn(vova)) = GREENFUN(sv, tv, site, tnew)
	   m_u(row(vova))  = GREENFUN(site, tnew, sv, tv)
	enddo
	m_s = g0000-alpha
	det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*bun/(nmnm+1)


! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 

!	call CheckGetName;    ! this slows things down a lot 

	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; 
	ktau(name)=tnew                ! global list

	if(pm==0)then                                    ! update TheMatrix
	   pm=1; matr(pm,pm)=1.d0/( g0000 -alpha)
	else
	   call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	a_a_v = a_a_v + un1

	end subroutine add



!-----------------
!--- Drop a kink
!-----------------
	subroutine drop(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio 
	integer :: site,j,nk, name, num,r,c, vova

	c_d_v = c_d_v + un1;   acpt=.false.
	prevstatus=status; status='_drop_'

	if(nmnm==0)return    ! nothing to drop yet :)
 
	num=nmnm*rndm()+1.d0                           ! play a kink to be dropped
	if(num>nmnm)num=nmnm; name=namelist(num); 
	site=ksite(name); nk = nkink(site)
	
!---------- determinant
	if(pm==1)then; det=1.d0/matr(pm,pm)
	else
	    r = row(name); c=clmn(name)
	    det = det_m1(pm,lda,matr,r,c)
	endif
!----------------------

	ratio = det*det*bun/nmnm
	ratio = un1/ratio

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
	if(pm==1)then; pm=0; nkink(site)=0
	else

	   do j=1,nk; if(kname(j,site)==name)exit   ! find name
	   enddo

	   kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site

	   vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
	   vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	   call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix

	endif

!	call CheckDropName(name);
	call DropName(name)

	a_d_v = a_d_v + un1

	end subroutine drop


!---------------------
!--- insert a kink @ specified position (site, tnew)
!---------------------   & return the new name & the det ratio & the new matrix
	subroutine insert(site, tnew,det, name)
	implicit none
	integer :: site, name
	real*8  :: tnew, det

	integer :: j,vova, nk
	real*8  :: tv

!------------- determinant
	if(pm==0)then; det=g0000-alpha  !ratio = g0000**2
	else

! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
       tv=ktau(vova)
	   m_v(clmn(vova)) = GREENFUN(ksite(vova), tv, site, tnew)
	   m_u(row(vova))  = GREENFUN(site, tnew, ksite(vova), tv)
	enddo
	m_s = g0000-alpha
	det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! pm==0
!---------------------------


! update configuration 
	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; 
	ktau(name)=tnew                ! global list

	if(pm==0)then                                    ! update TheMatrix
	   pm=1; matr(pm,pm)=1.d0/(g0000-alpha)
	else
	   call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	end subroutine insert


!---------------------
!--- remove a kink 'name'
!---------------------   & return the new name & the det ratio
	subroutine remove(det,name)
	implicit none

	integer :: site, name, nk, r,c, vova
	real*8  :: det

	site=ksite(name); nk = nkink(site)
	
!---------- determinant
	if(pm==1)then; det=1.d0/matr(pm,pm)
	else
	    r = row(name); c=clmn(name)
	    det = det_m1(pm,lda,matr,r,c)
	endif
!----------------------

! update configuration 
	if(pm==1)then; pm=0; nkink(site)=0
	else

	   do j=1,nk; if(kname(j,site)==name)exit   ! find name
	   enddo

	   kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site

	   vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
	   vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	   call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix

	endif

!	call CheckDropName(name);
	call DropName(name)

	end subroutine remove



!-----------------
!--- Add TWO kinks, same site
!-----------------
	subroutine add_2_same(acpt,det)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio
	integer :: name, name2, site,j,nk, vova, sv
	real*8  :: tnew, tv, tnew2, temp2x2(2,2)

	c_a_v2 = c_a_v2 + un1; acpt=.false.
	prevstatus=status; status='_add__' 

	if(pm>lda-5)then; call resize_matrix(-1); endif

!------- check
	if(nmnm>=nmnm_max)then;
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif
	if(present)then
	   open(out_unit, file=outputfname, position='append')
	   write(out_unit,*)' add_2 can''t handle ira/masha @ step=', step
	   close(out_unit)
	   call mystop
	endif
!---------- end check

	site=Nsite*rndm()+1.d0
	if(site>Nsite)site=Nsite    ! select a site to play on
	tnew = beta*rndm()                    ! select tau
	tnew2 = beta*rndm()             


!------------------ do two rank-1 updates
!	call insert(site,tnew,det1,name)
!	call insert(site,tnew2,det2,name2)
!
!	print*,'--- ? ', pm
!	print*,det1*det2, det1,det2
!
!	call remove(det1,name)
!	call remove(det2,name)
!----------------------------------------

!------------- determinant   
	m_s2(1,1) = g0000-alpha; m_s2(2,2)=g0000-alpha;
	m_s2(1,2) = GREENFUN(site, tnew, site,tnew2)
	m_s2(2,1) = GREENFUN(site, tnew2, site, tnew)
	temp2x2=m_s2

	if(pm==0)then; det=full_inv(2,2,temp2x2)   !m_s2(1,1)*m_s2(2,2)-m_s2(1,2)*m_s2(2,1)
	else
! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
       tv=ktau(vova); sv = ksite(vova)
	   m_v2(clmn(vova),1) = GREENFUN(sv, tv, site, tnew2)
	   m_v2(clmn(vova),2) = GREENFUN(sv, tv, site, tnew)
	   m_u2(row(vova),1)  = GREENFUN(site, tnew2, sv, tv)
	   m_u2(row(vova),2)  = GREENFUN(site, tnew, sv,tv)
	enddo

	det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_z2)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*beta*beta*U*U/( (nkink(site)+2)*(nkink(site)+1) )

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
!	call CheckGetName;    ! this slows things down a lot 
	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; ktau(name)=tnew                ! global list

	call getname(name2)
	nk=nkink(site)
	nkink(site)=nk+1; kname(nk+1,site)=name2
	ksite(name2)=site ; ktau(name2)=tnew2

	if(pm==0)then                                    ! update TheMatrix
	   pm=2; 
	   matr(1:2,1:2) = temp2x2
	else
	   call inv_p2(pm,lda,matr,det,m_s2,m_z2,m_v2,m_w2,m_x2)        
	endif


	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	clmn(name2) = pm-1 ; nm_clmn(pm-1) = name2     ! careful: @ inv_p2 pm -> pm+2
	row(name2)  = pm-1 ; nm_row(pm-1)  = name2


!------------------ do two rank-1 updates
!	call remove(det1,name)

!	call remove(det2,name2)
!
!	print*,'--- ? ', pm
!	print*,det1*det2, det1,det2
!	print*,det
!
!	call insert(site,tnew,det1,name)
!	call insert(site,tnew2,det2,name2)
!
!	pause
!----------------------------------------


	a_a_v2 = a_a_v2 + un1
!	call try_drop2_same(acpt,det,name,name2,ratio)

	end subroutine add_2_same



!-----------------
!--- TRY adding TWO kinks, same site
!-----------------
	subroutine try_add2_same(site,tnew,tnew2,ratio_)

	real*8  :: ratio, ratio_
	integer :: site, j, vova, sv
	real*8  :: tnew, tv, tnew2, temp2x2(2,2)

	if(pm==lda)then; call resize_matrix(-1); endif

!------------- determinant   
	m_s2(1,1) = g0000-alpha; m_s2(2,2)=g0000-alpha;
	m_s2(1,2) = GREENFUN(site, tnew, site, tnew2)
	m_s2(2,1) = GREENFUN(site, tnew2, site, tnew)
	temp2x2=m_s2

	if(pm==0)then; det=full_inv(2,2,temp2x2)   !m_s2(1,1)*m_s2(2,2)-m_s2(1,2)*m_s2(2,1)
	else
! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
	   tv=ktau(vova); sv = ksite(vova)
	   m_v2(clmn(vova),1) = GREENFUN(sv,tv, site, tnew2)
	   m_v2(clmn(vova),2) = GREENFUN(sv,tv, site, tnew)
	   m_u2(row(vova),1)  = GREENFUN(site ,tnew2, sv, tv)
	   m_u2(row(vova),2)  = GREENFUN(site ,tnew, sv, tv)
	enddo

	det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_z2)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*beta*U*beta*U/( (nkink(site)+2)*(nkink(site)+1) )

	! Metropolis
	if(abs(ratio-1.d0/ratio_)>1d-10)then
	  print*,'nuh @ step = ', step
	  print*,'ratio = ', ratio
	  print*,'1/ratio_= ', 1.d0/ratio_
	  call mystop
	else;
!	   print*,'try_add2 ok @ step = ', step
   	endif


	end subroutine try_add2_same


!-----------------
!--- Drop TWO kinks (same site)
!-----------------
	subroutine drop_2_same(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio, temp2x2(2,2), tname, tname2
	integer :: site,j,j2,jj,jj2,nk, name,name2, r,c,r2,c2, vova
	integer :: rr1,rr2,cc1,cc2
!	real*8 :: det_1, det_2, a_1(lda,lda), a_2(lda,lda)

	c_d_v2 = c_d_v2 + un1; acpt=.false.
	prevstatus=status; status='_drop2'; 
	if(nmnm<2)return    ! nothing to drop yet :)

! play a site
	site = Nsite*rndm()+1.d0; if(site>Nsite)site=Nsite
	nk = nkink(site)
	if(nk<2)return              ! nothing to drop yet

! play the unfortunate ones
	j=nk*rndm()+1.d0; if(j>nk)j=nk
	name = kname(j,site)

1111    continue
	j2=nk*rndm()+1.d0; if(j2>nk)j2=nk
	name2 = kname(j2,site)
	if(name==name2)goto 1111   ! can't drop it twice

	tname=ktau(name); tname2=ktau(name2)
	
!---------- determinant
	if(pm==2)then; !det=1.d0/matr(pm,pm)
	         temp2x2=matr(1:2,1:2)
		 det = full_inv(2,2,temp2x2)
		 det = 1.d0/det 
	else
	    r  = row(name);  c  = clmn(name)
	    r2 = row(name2); c2 = clmn(name2)
	    det = det_m2(pm,lda,matr,r,r2,c,c2)
	endif
!----------------------
	
	ratio = det*det*beta*beta*U*U/( nk*(nk-1) )
	ratio = 1.d0/ratio

! Metropolis
        if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
	if(pm==2)then; pm=0; nkink(site)=0
	else
	   if(j==j2)then;
	      print*,'  j==j2 @ step = ', step
	      call mystop
	   endif

	   if(j>j2)then; jj=j; jj2=j2    
	           else; jj=j2; jj2=j
	   endif                         ! jj>jj2 now

	   if(jj==nk)      then; kname(jj2,site) = kname(nk-1,site)     !swap jj2 & nk-1 
           elseif(jj==nk-1)then; kname(jj2,site) = kname(nk,site)
	   else;                 
				 kname(jj,site)  = kname(nk,site)
				 kname(jj2,site) = kname(nk-1,site)
	   endif
	   nkink(site) = nk-2

	   cc1=min(c,c2); cc2=max(c,c2)
	   rr1=min(r,r2); rr2=max(r,r2)      ! need to sort rows & clmns

	   if(rr2/=pm)then
	     vova = nm_row(pm); row(vova) = rr2; nm_row(rr2) = vova      ! matrix links
	     vova = nm_clmn(pm); clmn(vova) = cc2; nm_clmn(cc2) = vova
	   endif

	   if(rr1/=pm-1)then
	     vova = nm_row(pm-1); row(vova) = rr1; nm_row(rr1) = vova      ! matrix links
	     vova = nm_clmn(pm-1); clmn(vova) = cc1; nm_clmn(cc1) = vova
	   endif

	   call inv_m2(pm,lda,matr,det,rr1,rr2,cc1,cc2,m_z2,m_w2)      ! TheMatrix
	endif

!	call CheckDropName(name);
	call DropName(name)
	call DropName(name2)

!	call try_add2_same(site, tname,tname2,ratio)

	a_d_v2=a_d_v2 + 1.d0

	end subroutine drop_2_same



!-----------------
!--- TRY dropping two kinks
!-----------------
	subroutine try_drop2_same(acpt,det,name,name2,ratio_)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio, ratio_ , temp2x2(2,2)
	integer :: site, nk, name, name2, r, c, r2, c2

	prevstatus=status; status='_drop_'

	if(nmnm==0)return    ! nothing to drop yet :)
 
!	num=nmnm*rndm()+1.d0                           ! play a kink to be dropped
!	if(num>nmnm)num=nmnm; name=namelist(num); 
	site=ksite(name); nk = nkink(site)

	if(site/=ksite(name2))then
	   print*,'how come?'
	print*,'name:  ', name, ksite(name), ktau(name)
	print*,'name2: ', name2, ksite(name2), ktau(name2)
		call mystop
	endif
	
!---------- determinant
	if(pm==2)then; !det=1.d0/matr(pm,pm)
	         temp2x2=matr(1:2,1:2)
		 det = full_inv(2,2,temp2x2)
		 det = 1.d0/det 
	else
	    r  = row(name);  c  = clmn(name)
	    r2 = row(name2); c2 = clmn(name2)
	    det = det_m2(pm,lda,matr,r,r2,c,c2)
	endif
!----------------------


	ratio = det*det*beta*U*beta*U/nk/(nk-1)  !/nmnm
!	print*,' td2 (3) ', ratio, ratio_

! Metropolis
	if(abs(ratio-ratio_)>1d-10)then
	  print*,'nuh @ step = ', step
	  print*,'ratio = ', ratio
	  print*,'ratio_= ', ratio_
	  call mystop
   	endif

	end subroutine try_drop2_same



!-----------------
!--- Add TWO kinks, any sites
!-----------------
	subroutine add_2_any(acpt,det)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio
	integer :: name, name2, site,site2,j,nk, vova, sv
	real*8  :: tnew, tv, tnew2, temp2x2(2,2)

	c_a_v2a = c_a_v2a + un1; acpt=.false.
	prevstatus=status; status='_add2a' 

	if(pm>lda-5)then; call resize_matrix(-1); endif

!------- check
	if(nmnm>=nmnm_max)then;
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif
	if(present)then
	   open(out_unit, file=outputfname, position='append')
	   write(out_unit,*)' add_2 can''t handle ira/masha @ step=', step
	   close(out_unit)
	   call mystop
	endif
!---------- end check
	
	site=Nsite*rndm()+1.d0
	if(site>Nsite)site=Nsite    ! select a site to play on
	tnew = beta*rndm()                    ! select tau	

	site2=Nsite*rndm()+1.d0
	if(site2>Nsite)site2=Nsite
	tnew2 = beta*rndm()             

!------------- determinant   
	m_s2(1,1) = g0000-alpha; m_s2(2,2)=g0000-alpha;
	m_s2(1,2) = GREENFUN(site, tnew, site2, tnew2)
	m_s2(2,1) = GREENFUN(site2, tnew2, site, tnew)
	temp2x2=m_s2

	if(pm==0)then; det=full_inv(2,2,temp2x2)   !m_s2(1,1)*m_s2(2,2)-m_s2(1,2)*m_s2(2,1)
	else
! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
	   tv=ktau(vova); sv = ksite(vova)
	   m_v2(clmn(vova),1) = GREENFUN(sv, tv, site2, tnew2)
	   m_v2(clmn(vova),2) = GREENFUN(sv, tv, site, tnew)
	   m_u2(row(vova),1)  = GREENFUN(site2, tnew2, sv, tv)
	   m_u2(row(vova),2)  = GREENFUN(site, tnew, sv, tv)
	enddo

	det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_z2)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*bun*bun/( (nmnm+2)*(nmnm+1) )

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
!	call CheckGetName;    ! this slows things down a lot 
	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; ktau(name)=tnew                ! global list

	call getname(name2)
	nk=nkink(site2)
	nkink(site2)=nk+1; kname(nk+1,site2)=name2
	ksite(name2)=site2 ; ktau(name2)=tnew2

	if(pm==0)then                                    ! update TheMatrix
	   pm=2; 
	   matr(1:2,1:2) = temp2x2
	else
	   call inv_p2(pm,lda,matr,det,m_s2,m_z2,m_v2,m_w2,m_x2)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	clmn(name2) = pm-1 ; nm_clmn(pm-1) = name2     ! careful: @ inv_p2 pm -> pm+2
	row(name2)  = pm-1 ; nm_row(pm-1)  = name2

	a_a_v2a = a_a_v2a + un1
!	call try_drop2_any(name,name2,ratio)

	end subroutine add_2_any


!-----------------
!--- TRY dropping TWO kinks (any sites)
!-----------------
	subroutine try_drop2_any(name, name2,ratio_)
!
! This is complementary to subroutine add above
!
	real*8 :: ratio, ratio_ , temp2x2(2,2), det
	integer :: site,site2, name,name2, r,c,r2,c2

! play unfortunate ones
	site=ksite(name); !nk = nkink(site)
	site2 = ksite(name2)

!---------- determinant
	if(pm==2)then; !det=1.d0/matr(pm,pm)
	         temp2x2=matr(1:2,1:2)
		 det = full_inv(2,2,temp2x2)
		 det = 1.d0/det 
	else
	    r  = row(name);  c  = clmn(name)
	    r2 = row(name2); c2 = clmn(name2)
	    det = det_m2(pm,lda,matr,r,r2,c,c2)
	endif
!----------------------
	
	ratio = det*det*bun*bun/( nmnm*(nmnm-1) )
	

! Metropolis
        if(abs(ratio-ratio_)>1.d-8)then
	  print*,'try_drop_2_any fails @ step = ', step
	  print*,'ratio = ', ratio
	  print*,'ratio_= ', ratio_
	  call mystop
	endif


	end subroutine try_drop2_any



!-----------------
!--- Drop TWO kinks (any sites)
!-----------------
	subroutine drop_2_any(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio, ratio_ , temp2x2(2,2), tname, tname2
	integer :: site,site2,j,j2,jj,jj2,nk, name,name2, num,r,c,r2,c2, vova
	integer :: rr1,rr2,cc1,cc2

	c_d_v2a = c_d_v2a + un1; acpt=.false.
	prevstatus=status; status='drop2a'; 
	if(nmnm<2)return    ! nothing to drop yet :)

! play unfortunate ones
	num=nmnm*rndm()+1.d0                           ! play a kink to be dropped
	if(num>nmnm)num=nmnm; name=namelist(num); 
	site=ksite(name); !nk = nkink(site)
	tname = ktau(name)

2222    continue
	num=nmnm*rndm()+1.d0
	if(num>nmnm)num=nmnm; name2=namelist(num);
	if(name2==name)goto 2222
	site2 = ksite(name2)
	tname2 = ktau(name2)

!---------- determinant
	if(pm==2)then; !det=1.d0/matr(pm,pm)
	         temp2x2=matr(1:2,1:2)
		 det = full_inv(2,2,temp2x2)
		 det = 1.d0/det 
	else
	    r  = row(name);  c  = clmn(name)
	    r2 = row(name2); c2 = clmn(name2)
	    det = det_m2(pm,lda,matr,r,r2,c,c2)
	endif
!----------------------
	
	ratio = det*det*bun*bun/( nmnm*(nmnm-1) )
	ratio = 1.d0/ratio

! Metropolis
        if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
	if(pm==2)then; pm=0; nkink(site)  = nkink(site) -1 
			     nkink(site2) = nkink(site2)-1
	else

	     ! find their positions
	   do j=1,nkink(site); if(kname(j,site)==name)exit   ! find name
	   enddo
	   do j2=1,nkink(site2); if(kname(j2,site2)==name2)exit   ! find name
	   enddo

	   if(site==site2)then   ! same as in drop_2_same:

	     if(j>j2)then; jj=j; jj2=j2    
	             else; jj=j2; jj2=j
	     endif                         ! jj>jj2 now

	     nk=nkink(site)
	     if(jj==nk)      then; kname(jj2,site) = kname(nk-1,site)     !swap jj2 & nk-1 
             elseif(jj==nk-1)then; kname(jj2,site) = kname(nk,site)
	     else;                 
		   		   kname(jj,site)  = kname(nk,site)
				   kname(jj2,site) = kname(nk-1,site)
	     endif
	     nkink(site) = nk-2

	   else                  ! they're on the different sites; 

	     nk=nkink(site);  kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site	
	     nk=nkink(site2); kname(j2,site2) = kname(nk,site2); nkink(site2) = nk- 1 

	   endif

	   cc1=min(c,c2); cc2=max(c,c2)
	   rr1=min(r,r2); rr2=max(r,r2)      ! need to sort rows & clmns

	   if(rr2/=pm)then
	     vova = nm_row(pm); row(vova) = rr2; nm_row(rr2) = vova      ! matrix links
	     vova = nm_clmn(pm); clmn(vova) = cc2; nm_clmn(cc2) = vova
	   endif

	   if(rr1/=pm-1)then
	     vova = nm_row(pm-1); row(vova) = rr1; nm_row(rr1) = vova      ! matrix links
	     vova = nm_clmn(pm-1); clmn(vova) = cc1; nm_clmn(cc1) = vova
	   endif

	   call inv_m2(pm,lda,matr,det,rr1,rr2,cc1,cc2,m_z2,m_w2)      ! TheMatrix
	endif

!	call CheckDropName(name);
	call DropName(name)
	call DropName(name2)

    !call check_cnf

	a_d_v2a=a_d_v2a + 1.d0

	end subroutine drop_2_any



!-----------------
!--- TRY adding TWO kinks, any sites
!-----------------
	subroutine try_add_2_any(site,tnew,site2,tnew2,ratio_)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio, ratio_
	integer :: site, site2, j, vova, sv
	real*8  :: tnew, tv, tnew2, temp2x2(2,2)

	if(pm==lda)then; call resize_matrix(-1); endif

!------------- determinant   
	m_s2(1,1) = g0000-alpha; m_s2(2,2)=g0000-alpha;
	m_s2(1,2) = GREENFUN(site, tnew, site2, tnew2)
	m_s2(2,1) = GREENFUN(site2, tnew2, site, tnew)
	temp2x2=m_s2

	if(pm==0)then; det=full_inv(2,2,temp2x2)   !m_s2(1,1)*m_s2(2,2)-m_s2(1,2)*m_s2(2,1)
	else
! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
	   tv=ktau(vova); sv = ksite(vova)
	   m_v2(clmn(vova),1) = GREENFUN(sv, tv, site2, tnew2)
	   m_v2(clmn(vova),2) = GREENFUN(sv, tv, site, tnew)
	   m_u2(row(vova),1)  = GREENFUN(site2, tnew2, sv, tv)
	   m_u2(row(vova),2)  = GREENFUN(site, tnew, sv, tv)
	enddo

	det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_z2)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*bun*bun/( (nmnm+2)*(nmnm+1) )

! Metropolis
      if(abs(ratio-1.d0/ratio_)>1d-8)then; 
	  print*,'try_add_2_any fails @ step = ', step
	  print*,'ratio = ', ratio
	  print*,'ratio_= ', ratio_
	  call mystop
      endif


	end subroutine try_add_2_any



!---------- check some consistency
	subroutine check_cnf
	implicit none

	integer :: site, name, j

! check site links
	do site=1,Nsite
	   do j=1,nkink(site)
		name = kname(j,site)
		if( site /= ksite(name) )then
		   print*; print*,' site links @ step = ', step, status
		   print*,' name:  ', name, site, ksite(name)

		   call dump_site_lists

		   call mystop
		endif
	   enddo
	enddo

	end subroutine check_cnf


	subroutine dump_site_lists
	implicit none

	integer :: site, j, name

	print*,'site,              names  '
	do site=1,Nsite
	   print*,'** ',site,' ** ', kname(1:nkink(site),site)
	enddo

	print*
	print*,' name,          site   '
	do j=1,nmnm; name=namelist(j)
	   print*,name, ksite(name)
	enddo
	print*,'-----------------------'; print*


	end subroutine dump_site_lists



!------------------
!--- Measurements
!------------------
	subroutine measure

	integer :: j
	real*8  :: this_PE, this_KE, this_nn_szsz

	if(present)then       !   GF sector

	    j = floor((ktau(masha)-ktau(ira))*im_t_max/beta)    
	    im_t(j) = im_t(j)+1.d0                       ! t(ira-masha) histogram
	    im = im + 1.d0                               ! integrated correlator

	else                  !   Z sector

	i_b = i_b + 1 ; 	Z = Z + un1
	this_PE = 1.d0*nmnm/beta + alpha*alpha*U*Nsite          ! PE

	PE = PE + this_PE
	this_KE = diag_KE()
	KE = KE + this_KE                                     ! KE  (need to multiply by TWO)
	
	PE_2 = PE_2 + this_PE*this_PE
	PEKE = PEKE + this_PE*this_KE

	g_uu0=diag_dens() + 2.d0*alpha
	ndens = ndens + g_uu0                                   ! number density 
!	call dance_dance_2(this_PE,g_uu0)                       ! dance-dance correlators

	call measure_GF
	call dance_dance_3(this_PE,g_uu0)

        this_nn_szsz = measure_nn_szsz(g_uu0) ; 
        nn_szsz = nn_szsz + this_nn_szsz

	endif  ! present


!-------------- is a current block full?
	if( i_b == Z_b )then    ! wrap it up
	    b_n = b_n + 1; i_b=0
		PE_stat(b_n) = PE; PE = 0.d0
		KE_stat(b_n) = KE; KE = 0.d0
		ndens_stat(b_n) = ndens; ndens = 0.d0

		PE_2_stat(b_n) = PE_2 ; PE_2 = 0.d0
		PEKE_stat(b_n) = PEKE ; PEKE = 0.d0

                nn_szsz_stat(b_n)=nn_szsz; nn_szsz=0.d0

! normalize the spin structure factor
		gss_stat(b_n) = gss*(N(1)**(1+ H_eta))    ! no need to scale by Nsite; 
		gss=0.d0                                  ! cf "AFM-like" test [search for "AFM-like"]
                                                          ! gss estimator already gives the sublattice magnetization;
                                                          ! this way, Staudt's S(Q) = gss*L^d
		if(b_n == b_n_max)then                  
	        call collate( im_stat, b_n)
		    call collate( gss_stat, b_n)
			call collate( PE_stat, b_n )
			call collate( KE_stat, b_n )
			call collate( ndens_stat, b_n )
			call collate( PE_2_stat, b_n )
			call collate( PEKE_stat, b_n )
                  call collate( nn_szsz_stat, b_n )
			b_n = b_n/2; Z_b = Z_b * 2.d0
		endif


	endif



	end subroutine measure


!--------------------------------
!--- measure density: insert an extra kink @ random place
!--------------------------------
	real*8 function diag_dens()

	real*8 :: det,ti, tn
	integer :: site, j, vova, sv

! select where to insert a kink
	site=Nsite*rndm()+1.d0; if(site>Nsite)site=Nsite
	tn=beta*rndm()
	
!------------- determinant
	if(nmnm==0)then; det=g0000-alpha  
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     ti=ktau(vova); sv = ksite(vova)
	     m_u(j)  = GREENFUN(site, tn, sv, ti)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     ti=ktau(vova); sv = ksite(vova)
	     m_v(j) = GREENFUN(sv, ti, site, tn)
	  enddo

	  m_s = g0000-alpha

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! nmnm==0

! estimator
	diag_dens = det*2.d0   ! a factor of 2 reflects the spin summation

	end function diag_dens



!--------------------------------
!--- measure GF: insert an extra kink @ random place
!--------------------------------
	subroutine measure_GF()

	real*8 :: det,ti, t1, t2, xi_k, phase
	integer :: site1, site2, x1(d), x2(d), j, vova, sv
	integer :: k, i1,i2

! select where to insert the kinks
	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite
	x1(:)=x(:,site1); t1=beta*rndm()
	
	site2=Nsite*rndm()+1.d0; if(site2>Nsite)site2=Nsite
	x2(:)=x(:,site2); t2=beta*rndm()
	
!------------- determinant
	if(nmnm==0)then; det=GREENFUN(site1, t1, site2, t2)  
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     ti=ktau(vova); sv = ksite(vova)
	     m_u(j)  = GREENFUN(site1, t1, sv, ti)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     ti=ktau(vova); sv = ksite(vova)
	     m_v(j) = GREENFUN(sv, ti, site2, t2)
	  enddo

	  m_s = GREENFUN(site1, t1, site2, t2)

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! nmnm==0

! estimator -------------------
	i1=abs(x1(1)-x2(1)) ; i1=min(i1,N(1)-i1)
	i2=abs(x1(2)-x2(2)) ; i2=min(i2,N(2)-i2)
	
	do k=0, G_nfreq-1
		xi_k = 2.d0*pi*(k+0.5d0)/beta
		phase=xi_k*(t1-t2)
		gf_r(i1,i2,k)=gf_r(i1,i2,k) + det*cos(phase)*beta
		gf_i(i1,i2,k)=gf_i(i1,i2,k) + det*sin(phase)*beta
		zgf(i1,i2,k)=zgf(i1,i2,k) + 1.d0
	enddo
!----------------------------------------	

	end subroutine measure_GF




!-----------------------------------------------
!--- measure kinetic energy from the Z-sector
!-----------------------------------------------
	real*8 function diag_KE()

	integer :: site1, site2, vova, j, i, sv
	real*8  :: t, det,tv
!
!  This estmator is for the tight-binding dispersion ONLY (n.n. sites)
!
	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite; 

	i=dd*rndm()+1; if(i>dd)i=dd; site2 = ass(i,site1)
	t=beta*rndm()

!------------- determinant
	m_s = GREENFUN(site2, t, site1, t)        ! a corner 

	if(pm==0) then; det = m_s        
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     tv=ktau(vova); sv = ksite(vova)
	     m_u(j)  = GREENFUN(site2, t, sv, tv)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     tv=ktau(vova); sv = ksite(vova)
	     m_v(j) = GREENFUN(sv, tv, site1, t)
	  enddo

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif

! return value
	diag_KE = dd*Nsite*det  ! d*Nsite is the # of bonds,       
				! a factor of 2 is due to the spin summation
				! BUT : THERE MUST BE ANOTHER FACTOR OF 2 BECAUSE e(k)= -2t cos(...)  !!!


	end function diag_KE



!--------------------------------
!--- measure density-density correlators
!--------------------------------
!	subroutine dance_dance()
!
! This is a very time-consuming version since it calculates 
! the contributions to all L/2 distances
!
!	real*8 :: det1,det2, t1, tv 
!	integer :: site1,site2, x1(d),x2(d),xv(d),j,vova, dir,i


! play two extra half-kinks
!	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite
!	x1(:)=x(:,site1); t1=beta*rndm()

!	          do i=0,N2(1)

!	site2=site1;  dir=1 
!	     do j=0,i-1; site2= ass(dir,site2); enddo; 
!	x2=x(:,site2)



!--- determinantz
!	m_s = g0000

!	m_s2(1,1)=g0000; m_s2(1,2)=GREENFUN(x2,t1,x1,t1)
!	m_s2(2,2)=g0000; m_s2(2,1)=GREENFUN(x1,t1,x2,t1)


!	if(pm==0)then; det1=g0000**2; det2 = det1 - m_s2(1,2)*m_s2(2,1)
!	else

! 1st for up-down
!	  do j=1,pm; vova=nm_row(j)               ! fill a column
!	     xv=x(:,ksite(vova)); tv=ktau(vova)
!	     m_u(j)  = GREENFUN(x1,t1,xv,tv)
!	  enddo
!	  m_u2(1:pm,1)=m_u(1:pm)

!	  do j=1,pm; vova=nm_clmn(j)              ! fill a row
!	     xv=x(:,ksite(vova)); tv=ktau(vova)
!	     m_v(j) = GREENFUN(xv,tv,x1,t1)
!	  enddo
!	  m_v2(1,1:pm)=m_v(1:pm)

!	  det1 = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)     

! 2nd for up-down
!	  do j=1,pm; vova=nm_row(j)               ! fill a column
!	     xv=x(:,ksite(vova)); tv=ktau(vova)
!	     m_u(j)  = GREENFUN(x2,t1,xv,tv)
!	  enddo
!	  m_u2(1:pm,2)=m_u(1:pm)

!	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
!	     xv=x(:,ksite(vova)); tv=ktau(vova)
!	     m_v(j) = GREENFUN(xv,tv,x2,t1)
!	  enddo
!	  m_v2(2,1:pm)=m_v(1:pm)
	  
!	  det2 = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)     


! up-down
!	  det1 = det1*det2


! up-up
!	  det2 = det_p2(pm,lda,matr,m_u2,m_v2,m_s2,m_c2)

!	endif


!	g_ud(i) = g_ud(i) + det1 
!	g_uu(i) = g_uu(i) + det2 

!	                    enddo   ! i



!	end subroutine dance_dance



!===========================================================================
	subroutine dance_dance_3(this_PE,g_uu0)
	real*8 :: this_PE, g_uu0   ! sent through from measure()

	real*8  :: det2, this_uu, this_ud,t0, t1, tv
	integer :: j, i,k, x1(d),x2(d), site1,site2, r,c,num,name, vova, sv

	if(d .ne. 2) then 
		print*, ' Dens-dens correlator ERROR: this one only works in 2D !!'
		call mystop
	endif


!************************** g_{up,down}
	if(nmnm>0)then

	  ! select the kink to split
	  num=floor(nmnm*rndm())+1 ; name= namelist(num)
	  site1=ksite(name)

	  r = row(name) ; c= clmn(name)  	  ! find its position in the matrix
	  x1(:)=x(:,site1)  	  ! get its coordinate for determinant calculation
	  t0 = ktau(name)      	  ! Let the measuring time be the same as that on 'name' -- for easier testing

	  site2=floor(Nsite*rndm())+1
	  x2(:)=x(:,site2)

! **** Get the second determinant ************************
	    ! row
	    do j=1,pm; vova=nm_clmn(j)
	       tv=ktau(vova); sv = ksite(vova)
	       m_v(j) = GREENFUN(sv, tv, site2, t0) - GREENFUN(sv, tv, site1, t0)
	    enddo

	    ! column
	    do j=1,pm; vova=nm_row(j)
	       tv=ktau(vova); sv = ksite(vova)
	       m_z(j) = GREENFUN(site2, t0, sv, tv) - GREENFUN(site1, t0, sv, tv)
	    enddo

	    det2 = det_rc(pm,lda,matr,r,m_v,c,m_z,m_w)
!********************************************************

        i=abs(x2(1)-x1(1)) ; !!! i=min(i,N(1)-i)
        k=abs(x2(2)-x1(2)) ; !!! k=min(k,N(2)-k)
	    
		this_ud = nmnm*det2/(U*beta*Nsite)  + alpha*alpha
 
	else
	        site1=Nsite*rndm()+1; if(site1>Nsite)site1=Nsite
       	        site2=Nsite*rndm()+1; if(site2>Nsite)site2=Nsite
		x1(:)=x(:,site1) ; x2=x(:,site2)
	
		i=abs(x2(1)-x1(1)) ; !!!; i=min(i,N(1)-i)
		k=abs(x2(2)-x1(2)) ; !!! k=min(k,N(2)-k)

		this_ud = alpha*alpha

	endif

        !AFM-type check
        !this_ud=( (-1)**i + 1 ) /2

	g_ud(i,k)=g_ud(i,k)+this_ud
	z_ud(i,k)=z_ud(i,k)+1.d0

        gss = gss + (-1)**(i+k)*this_ud/2.d0    ! 1/4 since Sz=(n_up-n_down)/2 & 2 due to up-up + down-down
        gss0 = gss0 + (-1)**(i+k)*this_ud/2.d0

!*************************************** g_{up,up} : add two extra (1/2)-kinks

	t1=beta*rndm()
	site1=Nsite*rndm()+1; if(site1>Nsite)site1=Nsite
	x1 = x(:,site1)

	site2=Nsite*rndm()+1; if(site2>Nsite)site2=Nsite
	x2 = x(:,site2)

	if(site1==site2)then
	        this_uu = g_uu0/2.d0   ! just n_{up}
		i=0 ; k=0 
	else
       i=abs(x2(1)-x1(1))  !!!; i=min(i,N(1)-i)
       k=abs(x2(2)-x1(2)) !!! ; k=min(k,N(2)-k)

	   !--- determinantz
	   m_s2(1,1)=g0000-alpha; m_s2(1,2)=GREENFUN(site2, t1, site1, t1)
	   m_s2(2,2)=g0000-alpha; m_s2(2,1)=GREENFUN(site1, t1, site2, t1)
	   if(pm==0)then;  
	        det = (g0000-alpha)**2 - m_s2(1,2)*m_s2(2,1)
	   else
	        ! 1st for up-down
	        do j=1,pm; vova=nm_row(j)               ! fill a column
	           tv=ktau(vova); sv = ksite(vova)
	           m_u(j)  = GREENFUN(site1, t1, sv, tv)
	        enddo
	        m_u2(1:pm,1)=m_u(1:pm)

	        do j=1,pm; vova=nm_clmn(j)              ! fill a row
	           tv=ktau(vova); sv = ksite(vova)
	           m_v(j) = GREENFUN(sv, tv, site1, t1)
	        enddo
	        m_v2(1:pm,1)=m_v(1:pm)

	        ! 2nd for up-down
	        do j=1,pm; vova=nm_row(j)               ! fill a column
	           tv=ktau(vova); sv = ksite(vova)
	           m_u(j)  = GREENFUN(site2, t1, sv, tv)
	        enddo
	        m_u2(1:pm,2)=m_u(1:pm)

	        do j=1,pm; vova=nm_clmn(j)             ! fill a row
	           tv=ktau(vova); sv = ksite(vova)
	           m_v(j) = GREENFUN(sv, tv, site2, t1)
	        enddo
	        m_v2(1:pm,2)=m_v(1:pm)
	  
	        det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_w2)
	   endif
	   this_uu = det + alpha*alpha

	endif    ! site1==site2


        !AFM-type check
        !this_uu=( (-1)**i + 1 ) /2

	g_uu(i,k)=g_uu(i,k)  +this_uu
	z_uu(i,k)=z_uu(i,k)+1.d0
	gss = gss + (-1)**(i+k)*this_uu/2.d0     ! 1/4 since Sz=(n_up-n_down)/2 & 2 due to up-up + down-down
        gss0 = gss0 + (-1)**(i+k)*this_uu/2.d0
 

	end subroutine dance_dance_3
!===========================================================================




!===========================================================================
	real*8 function measure_nn_szsz(g_uu0)
	real*8 :: g_uu0   ! sent through from measure()

	real*8  :: det2, this_uu, this_ud,t0, t1, tv
	integer :: j, x1(d),x2(d), site1,site2, r,c,num,name, vova,dir, sv

	if(d .ne. 2) then 
		print*, ' Dens-dens correlator ERROR: this one only works in 2D !!'
		call mystop
	endif


!************************** g_{up,down}
	if(nmnm>0)then

	  ! select the kink to split
	  num=floor(nmnm*rndm())+1 ; name= namelist(num)
	  site1=ksite(name)

	  r = row(name) ; c= clmn(name)  	  ! find its position in the matrix
	  x1(:)=x(:,site1)  	  ! get its coordinate for determinant calculation
	  t0 = ktau(name)      	  ! Let the measuring time be the same as that on 'name' -- for easier testing
	  
          dir = dd*rndm()+1.d0; if(dir>dd)dir=dd
	  site2=ass(dir,site1)
	  x2(:)=x(:,site2)

! **** Get the second determinant ************************
	    ! row
	    do j=1,pm; vova=nm_clmn(j)
	       tv=ktau(vova); sv = ksite(vova)
	       m_v(j) = GREENFUN(sv, tv, site2, t0) - GREENFUN(sv, tv, site1, t0)
	    enddo

	    ! column
	    do j=1,pm; vova=nm_row(j)
	       tv=ktau(vova); sv = ksite(vova)
	       m_z(j) = GREENFUN(site2, t0, sv, tv) - GREENFUN(site1, t0, sv, tv)
	    enddo

	    det2 = det_rc(pm,lda,matr,r,m_v,c,m_z,m_w)
!********************************************************

		this_ud = nmnm*det2/(U*beta*Nsite)  + alpha*alpha
 
	else
	        site1=Nsite*rndm()+1; if(site1>Nsite)site1=Nsite
       	        !site2=Nsite*rndm()+1; if(site2>Nsite)site2=Nsite

                dir = dd*rndm()+1.d0; if(dir>dd)dir=dd
	        site2=ass(dir,site1)

        	this_ud = alpha*alpha

	endif


!*************************************** g_{up,up} : add two extra (1/2)-kinks

	t1=beta*rndm()
	site1=Nsite*rndm()+1; if(site1>Nsite)site1=Nsite
	x1 = x(:,site1)

	!site2=Nsite*rndm()+1; if(site2>Nsite)site2=Nsite
        
        dir = dd*rndm()+1.d0; if(dir>dd)dir=dd
	site2=ass(dir,site1)
	x2 = x(:,site2)

	if(site1==site2)then
	        this_uu = g_uu0/2.d0   ! just n_{up}
	else
	   !--- determinantz
	   m_s2(1,1)=g0000-alpha; m_s2(1,2)=GREENFUN(site2, t1, site1, t1)
	   m_s2(2,2)=g0000-alpha; m_s2(2,1)=GREENFUN(site1, t1, site2, t1)
	   if(pm==0)then;  
		det = (g0000-alpha)**2 - m_s2(1,2)*m_s2(2,1)
	   else
	        ! 1st for up-down
	        do j=1,pm; vova=nm_row(j)               ! fill a column
	           tv=ktau(vova); sv = ksite(vova)
	           m_u(j)  = GREENFUN(site1, t1, sv, tv)
	        enddo
	        m_u2(1:pm,1)=m_u(1:pm)

	        do j=1,pm; vova=nm_clmn(j)              ! fill a row
	           tv=ktau(vova); sv = ksite(vova)
	           m_v(j) = GREENFUN(sv, tv, site1, t1)
	        enddo
	        m_v2(1:pm,1)=m_v(1:pm)

	        ! 2nd for up-down
	        do j=1,pm; vova=nm_row(j)               ! fill a column
	           tv=ktau(vova); sv = ksite(vova)
	           m_u(j)  = GREENFUN(site2, t1, sv, tv)
	        enddo
	        m_u2(1:pm,2)=m_u(1:pm)

	        do j=1,pm; vova=nm_clmn(j)             ! fill a row
	           tv=ktau(vova); sv = ksite(vova)
	           m_v(j) = GREENFUN(sv, tv, site2, t1)
	        enddo
	        m_v2(1:pm,2)=m_v(1:pm)
	  
	        det = det_p2_(pm,lda,matr,m_u2,m_v2,m_s2,m_w2)
	   endif
	   this_uu = det + alpha*alpha

	endif    ! site1==site2


        measure_nn_szsz = (this_uu +this_ud)/2.d0 -0.25d0


	end function measure_nn_szsz
!===========================================================================



!================  DONE with measurements; various service functions below 



!------------------------
!--- Collate the array
!------------------------
      subroutine collate(arr,n)
      real*8, dimension(*) :: arr
      integer :: n, Zb        ! array length & # of measurements per array entry

      integer :: i

      do i=1,n/2
          arr(i)=arr(2*i)+arr(2*i-1)
      enddo

      end subroutine collate


!-------------------------------
!--- Analyze block statistics: average and dispersion
!-------------------------------
	subroutine bSTAT(arr,n,Zb,av,err)
	real*8, dimension(*) :: arr
	integer              :: n, Zb
	real*8               :: av, err

	real*8 :: av2, dif

	av  = sum( arr(1:n) )/Zb/n
	av2 = sum( arr(1:n)**2 )/Zb/Zb/n

	dif = av2 - av*av; 	if(dif<0.d0)dif=0.d0

	err = sqrt( dif ) / sqrt(1.d0*n)


	end subroutine bSTAT


!------------------------------------
!--- Print out and check the runtime
!------------------------------------
	subroutine prnt
	integer :: i,j,k
	real*8 :: xxx, dt

	real*8 :: PE_av, PE_err, KE_av, KE_err
	real*8 :: PE_2_av, PE_2_err, PEKE_av, PEKE_err, dbl, dbl_err,nn_av,nn_err
	real*8 :: gss_av, gss_err
	real*8 :: dens_av, dens_err,dummy, summ

	logical :: lastone

        character*99 :: fname


! 'how much on your watch?'
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! across the midnight

	time_elapsed = time_elapsed + dt

	lastone=.false.
	if( time_elapsed > time_limit )  lastone=.true.   ! time to wrap up
!-------------------------------------------

	open(1,file=outputfname,position='append')

	write(1,*)''
      write(1,*)'-------------------------', time_elapsed/3600,' hrs'
	
	if(i_t<step_t)then
	  write(1,*)' thermalization: ',1.d2*i_t/step_t,' % done'
	  write(1,*)' current U = ', bun/beta/Nsite
	  write(1,*)'  '
	endif


      do i=1,d
      write(1,fmt=700)i,N(i)
	end do
 700  format (8x,'N(',i1,') =', i4)

      write(1,*)' '
	write(1,*)' 1/T  = ', beta, '  -U  = ', bun/beta/Nsite
	write(1,*)' \mu  = ', mu, ' nnn = ', g0000
	write(1,*)' \alpha = ', alpha
      write(1,fmt=1111)phi(1),phi(2)
 1111 format( '  phi = ', G11.5, G11.5,G11.5)
	write(1,*)'  '

      write(1,*)' MC step (mln) =', step/1.d6,'  Z(mln) = ',Z/1.d6
	write(1,*)' Z_b(mln) = ',Z_b/1d6,' b_n = ', b_n
      write(1,*)' '

	write(1,fmt=771)nmnm,nm_max,nm_min,nm_av
 771  format(' nmnm-> ',I5,' now [ max = ',I5,' min = ',I5, ' average = ',G11.5,' ]')



!---- nullify counters
      nm_max=0; nm_min=i_hu; nm_av=nul


! is there enough statistics?
	if(b_n>3)then

!--- pot. energy -----------------------------------
	call bSTAT(PE_stat(1:b_n),b_n,Z_b,PE_av,PE_err)

	write(1,*)'  '
      write(1,fmt=701) PE_av, PE_err
 701  format(8x,'PE =',g12.5,4x,' +/- ',g10.3)

      
!--- kin. energy -----------------------------------
	call bSTAT(KE_stat(1:b_n),b_n,Z_b,KE_av,KE_err)

	write(1,*)'  '
      write(1,fmt=711) KE_av, KE_err
 711  format(8x,'KE =',g12.5,4x,' +/- ',g10.3)
      

!--- number density ---------------------------------
	call bSTAT(ndens_stat(1:b_n),b_n,Z_b,dens_av,dens_err)

	write(1,*)'  '
      write(1,fmt=702) dens_av, dens_err
 702  format(8x,'dens =',g12.5,4x,' +/- ',g10.3)
      

!--- d(dbl.occ)/dT------------------------------------
	call bSTAT(PE_2_stat(1:b_n),b_n,Z_b,PE_2_av,PE_2_err)
	call bSTAT(PEKE_stat(1:b_n),b_n,Z_b,PEKE_av,PEKE_err)

	dbl = PE_2_av - PE_av**2   + 2.d0*( PEKE_av - PE_av*KE_av )   ! KE is to by multiplied by two

	dbl_err = ( (PE_err/PE_av)**2 + (KE_err/KE_av)**2  ) *(PE_av*2.d0*KE_av)**2
	dbl_err = PE_2_err**2 + 2.d0*PE_err**2 + PEKE_err**2 + dbl_err
	dbl_err = sqrt( dbl_err )

	dbl     = dbl*beta*beta/U/Nsite
	dbl_err = dbl_err*beta*beta/U/Nsite
	

	write(1,*)'  '
	write(1,fmt=703) dbl, dbl_err
 703    format(8x,'dD/dT =',g12.5,4x,' +/- ',g10.3)


! Fermi energy & momentum for the non-interacting gas of the same density:
!	kf = (3.d0*pi*pi*dens_av)**(1.d0/3.d0)
!	ef = kf**2 !/2.d0                          ! m=1
	
!	write(1,fmt=777)ef, kf
! 777  format(8x,'E_F =',g12.5,4x,' k_F = ',g12.5 )

!--- spin structure factor ---------------------------
	call bSTAT(gss_stat(1:b_n),b_n,Z_b,gss_av,gss_err)
	write(1,*)' '
      write(1,fmt=723) gss_av, gss_err
 723  format(8x,'NMSF =',g12.5,4x,' +/- ',g12.5)

!--- nn spin-spin corr  ---------------------------
	call bSTAT(nn_szsz_stat(1:b_n),b_n,Z_b,nn_av,nn_err)
	write(1,*)' '
      write(1,fmt=733) nn_av, nn_err
 733  format(8x,'<SzSz>@nn =',g12.5,4x,' +/- ',g12.5)



!-----------------------------------------------------
	write(1,*)
	write(1,*)' address/accept:'
	write(1,*)' add/drop  ',a_a_v/(c_a_v+.001),' / ',a_d_v/(c_d_v+.001)
	write(1,*)' a/d2_same ',a_a_v2/(c_a_v2+.001),' / ',a_d_v2/(c_d_v2+.001)
	write(1,*)' a/d2_any  ',a_a_v2a/(c_a_v2a+.001),' / ',a_d_v2a/(c_d_v2a+.001)
	write(1,*)' recalc / backsc ', recalc_rate/(step + 0.0001)


	endif   ! b_n>3

	close(1)


!=============================  writeouts: various service distributions

        fname='gf_ii_xi'//trim(fnamesuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
!print*, g_ud(0,0,0)/z_ud(0,0,0)
	do k=0, G_nfreq-1
	   write(1,*) k, gf_r(0,0,k)/zgf(0,0,k)  ,gf_i(0,0,k)/zgf(0,0,k)
	enddo
	close(1)
	
        fname='gf_ij_xi0'//trim(fnamesuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
	do i=0,N(1)/2
	do j=0,N(2)/2
	   write(1,*) i,j, gf_r(i,j,0)/zgf(i,j,0) ,gf_i(i,j,0)/zgf(i,j,0)
	enddo
	enddo
	close(1)


        fname='g____ud'//trim(fnamesuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
	do i=0,N(1)-1
	do j=0,N(2)-1
	   write(1,*)i,j, g_uu(i,j)/z_uu(i,j)  ,g_ud(i,j)/z_ud(i,j)
	enddo
	enddo
	close(1)

        fname='g____uu'//trim(fnamesuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
	do i=0,N(1)-1
	do j=0,N(2)-1
	   write(1,*)i,j, g_uu(i,j)/z_uu(i,j) ! ,g_ud(i,j,k)/z_ud(i,j,k)
	enddo
	enddo
	close(1)

        fname='g____szsz'//trim(fnamesuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
	summ =0.d0
	do i=0,N(1)-1
	do j=0,N(2)-1
	   dummy = 2.d0*g_uu(i,j)/z_uu(i,j) + 2.d0*g_ud(i,j)/z_ud(i,j)  -1.d0 
           dummy = dummy / 4.d0 	
           write(1,*)i,j, dummy
	enddo
	enddo
	close(1)

! write kink number distrr
      xxx=sum(nmnm_distr(1:nmnm_sup))
      if(xxx==0.d0)xxx=1.d0

      open(2,file='nmnm'//trim(fnamesuffix)//'.dat')
	do j=0,nmnm_sup; write(2,*)j,nmnm_distr(j)/xxx
      enddo
	close(2)

! write det distr
	xxx=sum(det_distr(-det_distr_max:det_distr_max))
	if(xxx==0.d0)xxx=1.d0

	open(2,file='det'//trim(fnamesuffix)//'.dat')
	do j=-det_distr_max,det_distr_max
	   write(2,*)j, det_distr(j)/xxx
	enddo
	close(2)



!---------  uncomment this if you want to see the configuration -----------

! write kink visualization
!	open(2,file=trim(workdirpath)//'viz'//trim(fnamesuffix)//'.dat')
!	do i=1,nmnm; name=namelist(i)
!		write(2,*)x(:,ksite(name)),ktau(name)
!	enddo
!	close(2)
!
!	if(present)then
!
!	open(2,file=trim(workdirpath)//'viz_i'//trim(fnamesuffix)//'.dat')
!		write(2,*)x(:,ksite(ira)), ktau(ira)
!	close(2)
!
!	open(2,file=trim(workdirpath)//'viz_m'//trim(fnamesuffix)//'.dat')
!	    write(2,*)x(:,ksite(masha)),ktau(masha)
!	close(2)
!
!	endif
!---------------------------------------------------------------------------



! time to wrap up? --- write everything, release allocated memory and wrap up.
	if(lastone)then; 
	    open(1,file=outputfname,position='append')
		write(1,*)'Time to wrap up, dude..'; 
		close(1)
	    call wrt
		call mystop
	endif


	end subroutine prnt
	


!-------------------
!--- Write data
!-------------------
	subroutine wrt
	integer :: j, name

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'write .....'

! configuration
	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    write(1,*)beta
          write(1,*)present
           if(present)then
              write(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
              write(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
           endif
	  write(1,*)nmnm,pm
	  write(1,*)nm_row(1:pm)
	  write(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	     name=namelist(j)
	     write(1,*)name, ksite(name),ktau(name), row(name), clmn(name)
	  enddo
	  write(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! rndm state ---- this nothin' but useless rigor
!	rndmstr = 'rndm'//trim(fnamesuffix)//'.dat'
!	call save_rndm(rndmstr)
!
! There is a switch in the parfile as to read 
! the generator state or to seed it anew.
!

! statistics
	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  write(1,*)d,N
	  write(1,*)beta, U, mu
        write(1,*)phi
	  write(1,*)step, Z
	  write(1,*)Z_b, b_n, i_b
	  write(1,*)PE
	  write(1,*)PE_stat(1:b_n)
	  write(1,*)KE
	  write(1,*)KE_stat(1:b_n)
	  write(1,*)ndens
	  write(1,*)ndens_stat(1:b_n)
	  write(1,*)im
	  write(1,*)im_stat(1:b_n)
	  write(1,*)PE_2
	  write(1,*)PE_2_stat(1:b_n)
	  write(1,*)PEKE
	  write(1,*)PEKE_stat(1:b_n)
	  write(1,*)'------'
	  write(1,*)g_uu(0:N(1)-1, 0:N(2)-1)
	  write(1,*)z_uu(0:N(1)-1, 0:N(2)-1)
	  write(1,*)g_ud(0:N(1)-1, 0:N(2)-1)
	  write(1,*)z_ud(0:N(1)-1, 0:N(2)-1)
	  write(1,*)gss
	  write(1,*)gss_stat(1:b_n)
          write(1,*)nn_szsz
          write(1,*)nn_szsz_stat(1:b_n)
	  write(1,*)gf_r(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	  write(1,*)gf_i(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	  write(1,*)zgf(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	
	close(1)


	write(OUT_UNIT,*)'writing done!'
	close(OUT_UNIT)

	end subroutine wrt


!---------------------
!--- Read statistics
!---------------------
	subroutine rd_stat
	integer :: ddd
	real*8  :: dummy
	character*6 :: cvoid


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading stat .....'

	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  read(1,*)ddd,N
	      if(ddd/=d) call mystop
	  read(1,*)dummy, dummy, dummy
        read(1,*)phi
	  read(1,*)step, Z
	  read(1,*)Z_b, b_n, i_b
	  read(1,*)PE
	  read(1,*)PE_stat(1:b_n)
	  read(1,*)KE
	  read(1,*)KE_stat(1:b_n)
	  read(1,*)ndens
	  read(1,*)ndens_stat(1:b_n)
	  read(1,*)im
	  read(1,*)im_stat(1:b_n)
	  read(1,*)PE_2
	  read(1,*)PE_2_stat(1:b_n)
	  read(1,*)PEKE
	  read(1,*)PEKE_stat(1:b_n)
	  read(1,*)cvoid
	  read(1,*)g_uu(0:N(1)-1, 0:N(2)-1)
	  read(1,*)z_uu(0:N(1)-1, 0:N(2)-1)
	  read(1,*)g_ud(0:N(1)-1, 0:N(2)-1)
	  read(1,*)z_ud(0:N(1)-1, 0:N(2)-1)
	  read(1,*)gss
	  read(1,*)gss_stat(1:b_n)
          read(1,*)nn_szsz
          read(1,*)nn_szsz_stat(1:b_n)
	  read(1,*)gf_r(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	  read(1,*)gf_i(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	  read(1,*)zgf(0:N(1)/2, 0:N(2)/2, 0:G_nfreq-1)
	close(1)

	det_distr=0.d0; nmnm_distr=0.d0

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_stat


!---------------------
!--- Init statistics
!---------------------
	subroutine init_stat

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init stat .....'

	Z_b = 10
	i_b = 0; b_n = 0

	Z=0.d0
	im = 0.d0; im_stat=0.d0
	gss=0.d0; gss_stat=0.d0
	PE = 0.d0; PE_stat=0.d0
	KE = 0.d0; KE_stat=0.d0
	ndens = 0.d0; ndens_stat=0.d0
        nn_szsz=0.d0; nn_szsz_stat=0.d0
	g_uu=0.d0; g_ud=0.d0
	z_uu=0.d0; z_ud=0.d0
	PE_2 = 0.d0; PE_2_stat = 0.d0
	PEKE = 0.d0; PEKE_stat = 0.d0

        gss0=0.d0
	
	gf_r=0.d0; gf_i=0.d0; zgf=0.d0



!---------------------------------

	im_t=0.d0; det_distr=0.d0; nmnm_distr=0.d0; 

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)


	end subroutine init_stat


!--------------------------
!--- Read configuration
!--------------------------
	subroutine rd_cnf
	integer :: j,name, site
	real*8 :: bbb

	character*50 :: rndmstr

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading conf.....'
	
! configuration
	namelist(ira)=ira; namelist(masha)=masha
      do j=1,nmnm_max; namelist(j)=j; numnam(j)=j
      enddo; nmnm=0; pm=0

	nkink=0; kname=0

	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    read(1,*)beta;  
	    if( abs(beta-beta_fin)>1d-8 )then
		print*, 'reading beta, beta_fin = ', beta, beta_fin
		write(OUT_UNIT,*) ' reading beta, beta_fin = ', beta, beta_fin
	    endif

	    if(ans1==2)then ; beta_ini = beta ; 
	    endif

        read(1,*)present
           if(present)then
              read(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
              read(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
           endif
	  read(1,*)nmnm,pm
	     if(nmnm>nmnm_max)then
	        print*,'rd_cnf: nmnm= ',nmnm,' while nmnm_max= ',nmnm_max
			call mystop 
	     endif
!================== allocate enough memory to fit TheMatrix
	    deallocate(matr,m_u,m_v,m_w,m_z,nm_row,nm_clmn)    ! deallocate first
	    deallocate(m_u2,m_v2,m_z2,m_w2,m_x2,m_s2) 
		lda=(int(pm/64)+1)*64+128
		allocate(nm_row(lda),nm_clmn(lda))
		allocate(matr(lda,lda))
		allocate(m_u(lda),m_v(lda),m_w(lda),m_z(lda))
	allocate(m_u2(lda,2),m_v2(lda,2),m_z2(lda,2),m_w2(lda,2),m_x2(lda,2),m_s2(2,2)) 

!=========================================
	  read(1,*)nm_row(1:pm)
	  read(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	      read(1,*)name, site,ktau(name), row(name), clmn(name)
	      namelist(j)=name; numnam(name)=j
	      ktau(name) = ktau(name)
	         ksite(name)=site
	         nkink(site)=nkink(site)+1
	         kname(nkink(site),site)=name
	  enddo
	  read(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! calculate TheMatrix
	bbb = recalc_matrix(pm,lda,matr)
	

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_cnf


!--------------------------
!--- Init configuration
!--------------------------
	subroutine init_cnf
	integer :: j

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init conf.....'

! configuration
	nkink=0; kname=0
	ksite=0; ktau=0.d0
	row=0; clmn=0; pm=0
	nm_row=0; nm_clmn=0

! name manager
	namelist(ira)=ira; namelist(masha)=masha	
      do j=1,nmnm_max
          namelist(j)=j; numnam(j)=j
      enddo; nmnm=0; 

! no ira, masha @ the beginning
      present=.false.

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine init_cnf


!-------------------------------
!--- Increase LDA (matrix storage): 
!      if lnew == -1 then just add 512
!-------------------------------
	subroutine resize_matrix(lnew)
	integer :: lnew

	integer :: lda_new 
	real*8, allocatable :: matr_new(:,:)
	real*8, allocatable :: nm_r_new(:), nm_c_new(:)

	if(lnew==-1)then; lda_new=lda+512
	else; lda_new = lnew
	endif

	allocate(matr_new(lda_new,lda_new))
	allocate(nm_r_new(lda_new),nm_c_new(lda_new))

! save TheMatrix as is
	matr_new(1:pm,1:pm)=matr(1:pm,1:pm)
	nm_r_new(1:pm)=nm_row(1:pm); nm_c_new(1:pm)=nm_clmn(1:pm)

! Resize
	deallocate(matr); deallocate(m_u,m_v,m_w,m_z)
	deallocate(m_u2,m_v2,m_s2,m_z2,m_x2,m_w2)
	deallocate(nm_row,nm_clmn)

	lda=lda_new

	allocate(matr(lda,lda),m_u(lda),m_v(lda))
	allocate(m_w(lda),m_z(lda))
	allocate(m_u2(lda,2),m_v2(lda,2),m_z2(lda,2),m_w2(lda,2),m_x2(lda,2),m_s2(2,2)) 
	allocate(nm_row(lda),nm_clmn(lda))


! Restore
	matr(1:pm,1:pm)=matr_new(1:pm,1:pm)
	nm_row(1:pm)=nm_r_new(1:pm); nm_clmn(1:pm)=nm_c_new(1:pm)

! Cleanup
	deallocate(matr_new,nm_r_new,nm_c_new)


	end subroutine resize_matrix


!-------------------------------
!--- Recalculate the inverse from scratch (N**3 operations)
!-------------------------------
	real*8 function recalc_matrix(pm,lda,a)
	integer :: pm,lda     ! actual size & leading dimension of A
	real*8  :: a(lda,lda)

	integer :: i,j, vova, lesha
	real*8  :: ti,tj, det_big

	recalc_matrix = 0.;

	if(pm<2)return   ! no use to recalculate

	recalc_rate = recalc_rate + 1.d0

! build the matrix
	do j=1,pm; do i=1,pm
	  vova=nm_row(i); lesha=nm_clmn(j)
	  ti=ktau(vova); tj=ktau(lesha)
	  a(i,j)=GREENFUN(ksite(lesha), tj, ksite(vova), ti)
	  if(i==j)a(i,j)=a(i,j)-alpha
	enddo; enddo

! invert
	det_big = full_inv(pm,lda,a)
	
	recalc_matrix = det_big       ! return the absolute value of the determinant

	end function recalc_matrix




!----------------------
!---- Getname function for the Name Manager
!----------------------
      SUBROUTINE GetName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nick

      IF(nmnm<nmnm_max)THEN; 
	  nmnm=nmnm+1
        nick=namelist(nmnm); numnam(nick)=nmnm
      ELSE
	    ! has been moved into CheckGetName()
      ENDIF   

      END SUBROUTINE GetName

!------------------------------
!--- Check for GetName: nmnm>nmnm_max
!------------------------------
	subroutine CheckGetName

	if(nmnm>=nmnm_max)then; PRINT*,'GetName-> list is over!'
	   call mystop
	endif

	end subroutine CheckGetName

       
!----------------------
!---- DropName function for Name Manager
!----------------------
      SUBROUTINE DropName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nick
      INTEGER :: nm, lastname

      nm=numnam(nick)
      IF(nm<nmnm)THEN
         lastname=namelist(nmnm); namelist(nm)=lastname
         numnam(lastname)=nm; namelist(nmnm)=nick
         numnam(nick)=nmnm; nmnm=nmnm-1
      ELSE IF(nm==nmnm)THEN
		 nmnm=nmnm-1
      ELSE
	   ! has been moved into the CheckDropName()
      ENDIF
      RETURN 

      END SUBROUTINE DropName 


!------------------------------
!--- Check whether the name to be dropped exists
!------------------------------
	subroutine CheckDropName(nick)
	integer :: nick

	if(numnam(nick)>nmnm)then; 
	PRINT*,'DropName-> No such name:',nick; call mystop
	endif

	end subroutine CheckDropName



!---------------------------------
!--- Calculate the det from scratch & compare
!---------------------------------
	subroutine check_recalc

	real*8, allocatable :: a(:,:)
	real*8 :: det_1
	integer :: i,j

	if(pm<2)return
	
      allocate(a(lda,lda))
      det_1 = recalc_matrix(pm,lda,a)

      do i=1,pm; do j=1,pm
	    if(abs(a(i,j)-matr(i,j))>1d-8)then 
	       print*; print*,'check_recalc failed'
	       print*,i,j
		   print*,matr(i,j),a(i,j),abs(a(i,j)-matr(i,j))
	       print*,status,step,pm
	       print*,det_1
	       deallocate(a)
	       call mystop
	    endif
	enddo; enddo


	deallocate(a)


	end subroutine check_recalc


!---------------------------------
!--- Thermalize from scratch:
!      - diag updates only
!      - tune U
!---------------------------------
	subroutine therm1
	logical :: acpt
	real*8 :: r, pex, det !,beta_fin
	integer :: j,ii	

	integer :: howmanytimes = 20    ! adjust \beta this many times 
	integer :: howmanystepsperbeta  ! how many MC steps with each \beta

	step_t = n_sw*bun              ! bun is set while reading the parameter file
	howmanystepsperbeta = step_t/howmanytimes

	nm_av=0.d0; nm_max=0.d0; nm_min=0.d0
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0 ; i_r = 0.d0

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [diag] will be done by ', step_t,' steps'
	write(OUT_UNIT,*),'beta_ini, _fin = ', beta_ini, beta_fin
	close(OUT_UNIT)

!---------------------------------------------------
	do j=1,howmanytimes

	  if( howmanytimes-j<3)then; beta = beta_fin         ! adjust \beta
	  else;
	       pex =  1.d0*howmanytimes/j-0.99  
   	       beta = beta_fin + (beta_ini-beta_fin)*exp(-1.d0/pex)
	  endif
	  bmt = beta/mtau ; bmt1 = 1.d0/bmt 
	  bun = beta*U*Nsite	  

	  call tabulate(beta, mtau, mu)                     ! GFs

	  det = recalc_matrix(pm,lda,matr) 

	  do ii = 1,howmanystepsperbeta
		i_t = i_t+1.d0 ; i_r = i_r + 1.d0
		i_p = i_p+1.d0 ; i_w = i_w + 1.d0 ; step = step + 1.d0

!--- Update: diagonal add/drop ONLY
	   r = rndm()   ; 
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call add_2_same(acpt,det)
	   else if(r<prob(4))then;    call drop_2_same(acpt,det)	
	   else if(r<prob(5))then;    call add_2_any(acpt,det)
	   else if(r<prob(6))then;    call drop_2_any(acpt,det)	
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         nm_av=nm_av+(un1*nmnm)/step_p        
         IF(nmnm>nm_max)nm_max=nmnm
         IF(nmnm<nm_min)nm_min=nmnm

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	  enddo	   ! ii
	enddo      ! j


!______________________________________


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6, ' mln steps'
	write(OUT_UNIT,*)'  '
	close(OUT_UNIT)

	if(n_sw>0) call wrt   

	nm_av=0.d0; nm_max=0.d0; nm_min=0.d0

	end subroutine therm1


!---------------------------------
!--- Thermalize via worm scheme
!---------------------------------
	subroutine therm2
	logical :: acpt
	real*8 :: r, det


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [worm] will be done by ',  n_sw,' sweeps'
	close(OUT_UNIT)

	step_t = n_sw*bun
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_t=i_t+1.d0; i_r = i_r + 1.d0
	   if(step>=step_t)exit


!--- Update :)
	   r = rndm()   ; 
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call add_2_same(acpt,det)
	   else if(r<prob(4))then;    call drop_2_same(acpt,det)	
	   else if(r<prob(5))then;    call add_2_any(acpt,det)
	   else if(r<prob(6))then;    call drop_2_any(acpt,det)	
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	enddo

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6, ' mln steps'
	write(OUT_UNIT,*)'  '

	if(n_sw>0) call wrt   

	end subroutine therm2


!--------------------------
!--- clean up & stop 
!--------------------------
	subroutine mystop

	if(allocated(ass))deallocate(ass)
	if(allocated(back))deallocate(back); if(allocated(x))deallocate(x)
	if(allocated(nkink))deallocate(nkink)
	if(allocated(kname))deallocate(kname)
	if(allocated(ksite))deallocate(ksite)
	if(allocated(ktau))deallocate(ktau)
	if(allocated(row))deallocate(row)
	if(allocated(clmn))deallocate(clmn)
	if(allocated(nm_row))deallocate(nm_row)
	if(allocated(nm_clmn))deallocate(nm_clmn)
	if(allocated(matr))deallocate(matr)
	if(allocated(m_u))deallocate(m_u)
	if(allocated(m_v))deallocate(m_v)
	if(allocated(m_w))deallocate(m_w)
	if(allocated(m_z))deallocate(m_z)
	if(allocated(m_u2))deallocate(m_u2)
	if(allocated(m_v2))deallocate(m_v2)
	if(allocated(m_s2))deallocate(m_s2)
	if(allocated(m_z2))deallocate(m_z2)
	if(allocated(m_w2))deallocate(m_w2)
	if(allocated(m_x2))deallocate(m_x2)

	stop

	end subroutine mystop


	end program MAIN
