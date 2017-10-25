    module green_functions
    use lattice
    implicit none

    public :: GREENFUN, TABULATE
    public :: g0000, bmt, bmt1, mt_max, mtau, phi

    private

!--- Green function array
    integer, parameter :: mt_max=2000    ! Max. number of the mesh points for tau
    integer :: mtau
    real*8, allocatable :: GR_DAT(:,:,:), GRD_DAT(:,:,:)

    real*8 :: bmt, bmt1     ! a shorthand for beta/mtau, and its inverse
    real*8 :: g0000         ! green function w/ all arg = 0, non-interacting density  
    real*8, allocatable :: phi(:)   ! BC twist

    real*8 :: beta_tabulated

    contains


!----------------------------------------------
!---  Green Function, spline interpolation 
!----------------------------------------------
    real*8 function GREENFUN(site1, tau1, site2, tau2)
    implicit none
    integer :: site1, site2
    double precision :: tau, tau1, tau2, dt, gre
    integer :: x1(1:d), x2(1:d),j, sgn

    integer :: nxy(1:d)
    integer :: nx, ny, nta  !, ntb
    double precision :: tta,ttb,ga,gb,c, gra,grb   !,p

! prepare \tau
    tau=tau1-tau2
    dt=tau; sgn=1

    if(tau < 1.d-14)then; dt=beta_tabulated + tau; sgn=-1; endif
! G(t=0) must be understood as G(t-> -0) = -G(t=\beta)
! A long way to accomplish this is below, commented out. A short way is above :).
!
!     if (abs(tau) < 1.d-14) then; dt = beta; sgn=-1
!     else if (tau > 0.d0) then; dt=tau; sgn=1
!   else; dt=beta+tau; sgn=-1
!    end if
!----------------------------------------



! prepare coords, don't forger about PBC
    x1 = x(:, site1); x2 = x(:, site2)
    nxy = distance_v(x1, x2)
    nx = nxy(1); ny = nxy(2)

!----------------------------------- spline

    nta=dt*bmt1    ! Recall, bmt=beta/mtau, bmt1=1/bmt 
    tta=dt-nta*bmt 
    ttb=tta - bmt   

!cccccccccccccccccccccccccccccccccccccc
      
    ga=GR_DAT(nta,nx,ny)
    gb=GR_DAT(nta+1,nx,ny)

    gra=GRD_DAT(nta,nx,ny)
    grb=GRD_DAT(nta+1,nx,ny)

    c=(ga-gb)*bmt1

    gre=(c+gra)*ttb + (c+grb)*tta
    gre=gre*tta*ttb*bmt1 + gb*tta-ga*ttb
    gre=gre*bmt1

    GREENFUN = gre*sgn

    end function GREENFUN


!-------------------------------------
!     Tabulates Green function and its time derivate at positive tau.
!-------------------------------------
    subroutine TABULATE(beta, mtau, mu)
    implicit none
!
! These tabulated values will be used in greenfun() for the spline interpolation. 
!
    real*8, intent(in) :: beta
    integer, intent(in) :: mtau
    real*8, intent(in) :: mu

    integer :: ntab
    integer :: i, mx, my !, mz
    integer :: nt, nx, ny !, nz
    double precision :: phase, eps, gamma
    double precision :: pp(d),term
    real*8 :: expet(0:mtau), ww
    real*8 :: kx,ky !,kz

    real*8 :: BC_twist, cos_BC_twist
    integer :: site
    real*8 :: ttt

    do i=1,d; pp(i)=4.d0*asin(1.d0)/N(i); enddo

! Green functions
	ntab=maxval(N(:)/2+1)                          ! the longest jump w/PBC is N/2
	allocate( GR_DAT(0:mt_max+1, 0:ntab-1, 0:ntab-1) )
	allocate( GRD_DAT(0:mt_max+1, 0:ntab-1, 0:ntab-1) )

 
!------------------ Cosines  ------------------
!      do i=1,d
!      do mx=0, N(i)-1;                ! The single-particle dispersion is tight-binding,
!              co(mx,i)=-2.d0*cos(pp(i)*mx);   ! thus the spectrum is cosine.
!    end do
!    enddo
!----------------------------------------------

! nullify'em
    GR_DAT=0.d0; GRD_DAT=0.d0    

! sum over momentums 1st                 ! This weird loop sequence is UNAVOIDABLE
!    do mz=0, N(3)-1            ! in order to ensure the data locality:
    do my=0, N(2)-1;           ! otherwise tabulation for L>10 takes
        do mx=0, N(1)-1;           ! hours and hours.

            kx=pp(1)*mx + phi(1)/N(1)
            ky=pp(2)*my + phi(2)/N(2)
!           kz=pp(3)*mz + phi(3)/N(3)

! spectrum
!           eps=co(mx,1)+co(my,2)+co(mz,3)-mu
!           eps = -2.d0*( cos(kx) + cos(ky) + cos(kz) ) -mu
            eps = -2.d0*( cos(kx) + cos(ky) ) -mu
!           eps = -mu             ! atomic limit, dammit!


            gamma=-eps*beta
            gamma=exp(gamma)+1.d0
            gamma=-1.d0/gamma

            ww = exp(-eps*bmt)             ! bmt=beta/mtau
            do nt=0,mtau; expet(nt)=ww**nt
            enddo

! coordinates 
           !do nz = 0, Ntab-1; 
            do ny = 0, Ntab-1; do nx = 0, Ntab-1
                !phase=pp(1)*nx*mx +pp(2)*ny*my+pp(3)*nz*mz
                phase = kx*nx + ky*ny ! + kz*nz

                do nt=0,mtau
                    term = cos(phase)*expet(nt)*gamma/Nsite
                    GR_DAT(nt,nx,ny) = GR_DAT(nt,nx,ny) + term
                    GRD_DAT(nt,nx,ny) = GRD_DAT(nt,nx,ny) -eps*term
                enddo                  ! tau
            enddo; enddo !; enddo  ! coordinates
        enddo; enddo ! ; enddo      ! momentum
!------------------------

! store beta used in tabulations for further use in GREENFUN
    beta_tabulated = beta

! fill in fictitious nt=mtau+1, see GREENFUN for explanation
    GR_DAT(mtau+1,:,:)=0.d0; GRD_DAT(mtau+1,:,:)=0.d0

! fill g0000
	site = 1; ttt=0.d0
	g0000 = GREENFUN(site,ttt,site,ttt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!	site0=1; site1=site0
!	do j=1,N(1)
!	  print*,site0,site1,GR_DAT(100,site0,site1),GRD_DAT(100,site0,site1)
!	  site1=ass(1,site1); site1=ass(2,site1)
!	enddo

    print*,'done with TABULATE'

    end subroutine TABULATE

    end module green_functions
