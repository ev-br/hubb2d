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
    integer :: site1, site2, sgn
    double precision :: tau, tau1, tau2, dt, gre


    integer :: nta  !, ntb
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

!----------------------------------- spline

    nta=dt*bmt1    ! Recall, bmt=beta/mtau, bmt1=1/bmt 
    tta=dt-nta*bmt 
    ttb=tta - bmt   

!cccccccccccccccccccccccccccccccccccccc
      
    ga=GR_DAT(nta, site1, site2)
    gb=GR_DAT(nta+1, site1, site2)

    gra=GRD_DAT(nta, site1, site2)
    grb=GRD_DAT(nta+1, site1, site2)

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

    real*8, allocatable :: ham(:,:)
    integer :: site, site1, j
    real*8 :: factor, ww,ttt,term, gamma, expet(0:mtau)
    integer :: nt

    ! lapack stuff
    character*1 :: jobz,uplo
    integer     :: ldh, lwork,info
    real*8, allocatable  :: work(:), eps(:)

    print*,'start with TABULATE'

! Green functions
    if(allocated(GR_DAT))deallocate(GR_DAT)
    if(allocated(GRD_DAT))deallocate(GRD_DAT)
	allocate( GR_DAT(0:mt_max+1, 1:Nsite, 1:Nsite) )
	allocate( GRD_DAT(0:mt_max+1, 1:Nsite, 1:Nsite) )

    GR_DAT=0.d0; GRD_DAT=0.d0


!--------- build the hamiltonian
    allocate(ham(1:Nsite, 1:Nsite)) 
    ham=0.d0
    do site=1, Nsite
       ham(site,site)=ham(site, site) - mu
       do j=1, dd
          site1 = ass(j, site)
          if(site1 > 0) ham(site, site1) = ham(site, site1) - 1.d0
       enddo
    enddo

!  compute eigenvalues; for LAPACK parameters and arguments, see
!  http://www.netlib.org/lapack/double/dsyev.f
!  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

    jobz='V'  ! compute eigenvectors
    uplo='U'  ! use upper diag of ham(:,:) --- doesn't matter really
    ldh=Nsite
    lwork=12*Nsite
    allocate( work(lwork), eps(Nsite) )

! query the optimal workspace size
    call dsyev(jobz,uplo,Nsite,ham,ldh,eps,work,-1,info)
    lwork=work(1)
    deallocate(work); allocate(work(lwork))

! diagonalize
    call dsyev(jobz, uplo, Nsite, ham, ldh, eps, work, lwork, info)

    if(info/=0)then; print*,'*** dsyev returns info = ', info
        print*,'*** check the TABULATE routine'
        stop
    endif

!	print*,'optimal lwork =', work(1),lwork
!
!	print*,'eigenvalues:'
!	do j=1,Nsite
!		print*,eps(j)
!	enddo
!
!	print*; print*,'----------- eigenvector # '
!	do j=1,Nsite;
!		!write(1,*)j,ham(j,6)
!		print*,j,ham(j,1)
!	enddo
!	print*; print*; print*

!------------- have the spectrum, proceed to the GFs
    do j=1, Nsite
        gamma=-eps(j)*beta
        gamma=exp(gamma)+1.d0
        gamma=-1.d0/gamma

        ww = exp(-eps(j)*bmt) ! bmt=beta/mtau
        do nt=0,mtau; expet(nt)=ww**nt
        enddo

        do site=1,Nsite;
            do site1=1,Nsite
                factor = ham(site, j)*ham(site1, j)
                do nt=0,mtau
                    term = factor*expet(nt)*gamma    !/Nsite
                    GR_DAT(nt, site, site1) = GR_DAT(nt, site, site1) + term
                    GRD_DAT(nt, site, site1) = GRD_DAT(nt, site, site1) -eps(j)*term
                enddo
            enddo
        enddo ! site, site1

	enddo   ! j: eigenvalues
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
