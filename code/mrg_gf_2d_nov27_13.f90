!
!   Analyze the statistics
!
	implicit none

!------------------- Phys. parameters --------------------

	real*8,parameter :: pi=3.14159265358979323846d0

      integer, parameter  :: d=2, dd=4            ! d is the global dimension, dd=d+d  
	integer, allocatable :: ddd(:)
      integer, allocatable :: N(:,:)           ! Array for Numbers of sites per dimensions
      integer, allocatable :: Nsite(:)            ! Total number of sites

	real*8,allocatable :: U(:)          ! MINUS interaction strength
	real*8,allocatable :: beta(:)       ! inverse temperature
	real*8,allocatable :: mu(:)         ! chem. potential

	real*8,parameter :: H_eta=0.03d0   ! eta in Heisenberg universality class

!--------------  Measurement related -------------------

	integer, parameter   :: b_n_max=1000
	integer,allocatable  :: Z_b(:)           ! # block size
	integer,allocatable  :: b_n(:), i_b(:)    ! filled block nbr. & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block
! partition function
	real*8,allocatable :: Z(:)

! diagonal estimators:
	real*8,allocatable  :: PE(:), KE(:), ndens(:)     ! PE, KE and density estimators 
	real*8,allocatable  :: PE_stat(:,:), KE_stat(:,:), ndens_stat(:,:) ! files

	real*8, allocatable :: PE_2(:), PEKE(:)
	real*8, allocatable :: PE_2_stat(:,:), PEKE_stat(:,:)

! integrated ira-masha correlator:
	real*8,allocatable  :: im(:)               ! estimator
	real*8,allocatable  :: im_stat(:,:) ! file

! Integrated spin-spin correlator:
	real*8, allocatable :: gss(:)
	real*8, allocatable :: gss_stat(:,:)
	
	real*8, allocatable :: gf_r(:,:,:,:), gf_i(:,:,:,:), zgf(:,:,:,:)
	real*8, allocatable :: gf_r_all(:,:,:), gf_i_all(:,:,:)
	real*8, allocatable :: gf_k_r(:,:,:), gf_k_i(:,:,:)
	real*8, allocatable :: sigma_k_r(:,:,:), sigma_k_i(:,:,:)
	real*8, allocatable :: sigma_bfe_r(:,:), sigma_bfe_i(:,:)
	
	integer, allocatable :: BF_index_vector_p(:,:)
	
	integer :: G_nfreq=64 !! number of frequencies for the Green's function
	integer :: S_nbf_p=3 !! number of basis functions per dimension for Sigma expansion

! <SzSz>@nn via a direct estimator
	real*8, allocatable :: nn_szsz(:)
	real*8, allocatable :: nn_szsz_stat(:,:)

! density-density
	real*8, allocatable :: g_uu(:,:,:), g_ud(:,:,:)
        real*8, allocatable :: z_uu(:,:,:), z_ud(:,:,:)

	real*8, allocatable :: step(:)

 
!-----------  overall stat ---------------
	integer :: b_n_all, i_all
	real*8, allocatable :: PE_all(:), KE_all(:),im_all(:),ndens_all(:)
	real*8, allocatable :: gss_all(:)
	real*8, allocatable :: chi_all(:)

	real*8, allocatable :: PE_2_all(:), PEKE_all(:)

	real*8, allocatable :: ene_sc(:)  ! (PE + KE per fermion ) / E_F
	real*8, allocatable :: tef(:)     ! T/E_F


	real*8, allocatable :: g_uu_all(:,:), g_ud_all(:,:), g_ss(:,:)
        real*8, allocatable :: z_uu_all(:,:), z_ud_all(:,:)

        real*8, allocatable :: nn_szsz_g(:)    ! <SzSz>@nn from g_uu & g_ud arrays
        integer :: nn_szsz_g_howmany

        real*8, allocatable :: nn_szsz_all(:)  ! <SzSz>@nn via direct estimator
        

	real*8 :: step_all, Z_all

	real*8 :: av, err
	real*8 :: ef, def,kf

!-------------------- discard what -----------
	integer, allocatable :: istart(:,:)
	integer :: istart_num                ! how many tries
	integer :: klm, tmp
	logical, allocatable ::  empty(:)


!-------------------- working variables -----------------

	integer :: nfls     ! # of stat. files to process
	integer :: i,j,p,k, nf, iii

	integer :: who_max,Z_b_max

	real*8 :: eef, npart

	character*50 :: fname, fsuffix

!=========================================


!-------------- read the files ----------------

	print*,' Enter # of files: '; read*,nfls

	
	allocate( ddd(nfls),N(nfls,d),Nsite(nfls) )
	allocate( U(nfls),beta(nfls),mu(nfls) ) 
	allocate( Z_b(nfls),b_n(nfls),i_b(nfls),Z(nfls),step(nfls) )
	allocate( PE(nfls),KE(nfls),ndens(nfls),im(nfls) )
	allocate( PE_2(nfls), PEKE(nfls))
	allocate( PE_stat(nfls,b_n_max),KE_stat(nfls,b_n_max) )
	allocate( PE_2_stat(nfls,b_n_max), PEKE_stat(nfls,b_n_max) )
	allocate( ndens_stat(nfls,b_n_max),im_stat(nfls,b_n_max) )
	allocate( gss(nfls),gss_stat(nfls,b_n_max) )
	allocate( nn_szsz(nfls),nn_szsz_stat(nfls,b_n_max) )
	allocate( g_uu(nfls,0:20,0:20), g_ud(nfls,0:20,0:20) )
	allocate( z_uu(nfls,0:20,0:20), z_ud(nfls,0:20,0:20) )

	do i=1,nfls; 
	   print*,' Enter filename #',i,' : '; read*,fname
	   call rd_stat(i,fname)
	   print*,i,' : ',trim(fname),' is read.'
	enddo

	read*,fsuffix


!-------------- check files ---------------------
	do i=1,nfls
	  if( any ( mu/=mu(1) ) )then; print*,'mu...',mu(:)
	  endif
	  if( any ( beta/=beta(1) ) )then; print*,'beta...',beta(:)
	  endif
	  if( any ( U/=U(1) ) )then; print*,'U...',U(:)
	  endif
	  if( any ( N(:,1)/=N(1,1) ) )then; print*,'N(,1)...',N(:,1)
	  endif
	  if( any ( N(:,2)/=N(1,2) ) )then; print*,'N(,2)...',N(:,2)
	  endif
	enddo


        if( any( N(:,:) >20 ) )then; print*,'there''s gonna be a problem with g_uu'
                                     stop
        endif

    call tabulate_basis_function_index_vectors
	
!--------  do density-density correlatorz ----------------
	allocate( g_uu_all(0:N(1,1)-1,0:N(1,2)-1), g_ud_all(0:N(1,1)-1,0:N(1,2)-1))
	allocate( z_uu_all(0:N(1,1)-1,0:N(1,2)-1), z_ud_all(0:N(1,1)-1,0:N(1,2)-1))
	allocate( g_ss(0:N(1,1)-1,0:N(1,2)-1))
	
	Z_all = sum(Z(:))
        g_uu_all=0.d0 ; g_ud_all=0.d0; g_ss=0.d0

	do i=0,N(1,1)-1;
	do j=0,N(1,2)-1;
                do iii = 1,nfls
	  g_uu_all(i,j)= g_uu_all(i,j)+  g_uu(iii,i,j)/z_uu(iii,i,j)
	  g_ud_all(i,j)= g_ud_all(i,j)+  g_ud(iii,i,j)/z_ud(iii,i,j)
                enddo
	enddo
	enddo

	open(1,file='g_uu'//trim(fsuffix)//'.dat') 
	open(2,file='g_ud'//trim(fsuffix)//'.dat')
	
	do i=0,N(1,1)-1;
	do j=0,N(1,2)-1;
	       write(1,*) i,j, g_uu_all(i,j)
	       write(2,*) i,j, g_ud_all(i,j)
	enddo
	enddo

	close(1); close(2)




	open(1,file='g_ss'//trim(fsuffix)//'.dat')
	
	do i=0,N(1,1)-1;
	do j=0,N(1,2)-1;
		g_ss(i,j)=2.0d0*g_uu_all(i,j)+2.0d0*g_ud_all(i,j)-1.d0
		write(1,*) i,j, g_ss(i,j)/4
	enddo
	enddo
		
	close(1)


! nn <Sz:Sz> correlations  from g_uu & g_ud arrays
        nn_szsz_g_howmany=nfls*d
        allocate(  nn_szsz_g(1:nn_szsz_g_howmany) )
        nn_szsz_g=0.d0 ; j=1
        do i=1,nfls
           nn_szsz_g(j) = (g_uu(i,0,1)/z_uu(i,0,1)+g_ud(i,0,1)/z_ud(i,0,1)-0.5d0)/2.d0
           nn_szsz_g(j+1) = (g_uu(i,1,0)/z_uu(i,1,0)+g_ud(i,1,0)/z_ud(i,1,0)-0.5d0)/2.d0
           j=j+2
        enddo

!        do j=1,nfls*d
!                print*,j,nn_szsz(j)
!        enddo
!        pause


!------------- Magnetic structure factor ---------------
	
!	gss(1)=0
!
!	do i=0,N(1,1)-1;
!	do j=0,N(1,2)-1;
!	do k=0,N(1,3)-1;
!
!		gss(1)=gss(1) + (-1)**(i+j+k)*g_ss(i,j,k)
!
!	enddo
!	enddo
!	enddo
!
!	print*, ' '
!	print*, 'MSF=', gss(1)
!	print*, ' '
!
!	gss=gss*(N(1,1)**(1.03d0)) !/(Nsite(1))
!
!	print*, ' '
!	print*, 'NMSF=', gss(1)
!	print*, ' '
	
!-------------------------------------------------------


	deallocate(g_uu_all, g_ud_all)

!---------------------------------------------------------

! ----- do GF --------------------------------------------

        gf_r_all=0.d0 ; gf_i_all=0.d0

	do i=0,N(1,1)/2;
	do j=0,N(1,2)/2;
		do nf=0, G_nfreq-1 
                do iii = 1,nfls
!					print*,'------------------------------'
!					print*, 'iii=', iii
!					print*, 'i,j,nf=', i,j,nf
!					print*, 'gr,gi,z=', gf_r(iii,i,j,nf), gf_i(iii,i,j,nf), zgf(iii,i,j,nf)
					
					gf_r_all(i,j,nf)= gf_r_all(i,j,nf) +  gf_r(iii,i,j,nf)/(nfls*zgf(iii,i,j,nf))
					gf_i_all(i,j,nf)= gf_i_all(i,j,nf) +  gf_i(iii,i,j,nf)/(nfls*zgf(iii,i,j,nf))
				enddo
		enddo
	enddo
	enddo


	fname='gf_ii_xi'//trim(fsuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
!print*, g_ud(0,0,0)/z_ud(0,0,0)
	do k=0, G_nfreq-1
	   write(1,*) k, gf_r_all(0,0,k),  gf_i_all(0,0,k)
	enddo
	close(1)
	
	fname='gf_ij_xi0'//trim(fsuffix)//'.dat'
	open(1,file=trim(adjustl(fname)))
!print*, g_ud(0,0,0)/z_ud(0,0,0)
	do i=0,N(1,1)/2
	do j=0,N(1,2)/2
	   write(1,*) i,j, gf_r_all(i,j,0), gf_i_all(i,j,0)
	enddo
	enddo
	close(1)
	
	call fourier_GF
	call plot_sigma_k

	print*, 'got Sigma well'

! --------------------------------------------------------


!-------------- equate block sizes ----------------

	Z_b_max = maxval( Z_b )

!c	print*,'before: ',Z_b, Z_b_max


	do i = 1,nfls
	
	  p = Z_b_max / Z_b(i); !print*,'i,p(i) =  ',i,p
	  Z_b(i) = Z_b(i)*p; b_n(i) = b_n(i)/p

	  do j = 1,b_n(i);  	 
	     PE_stat(i,j) = 1.d0*sum( PE_stat(i,p*(j-1)+1:p*j) ) / p
	     KE_stat(i,j) = 1.d0*sum( KE_stat(i,p*(j-1)+1:p*j) ) / p

	     PE_2_stat(i,j) = 1.d0*sum( PE_2_stat(i,p*(j-1)+1:p*j) ) / p
	     PEKE_stat(i,j) = 1.d0*sum( PEKE_stat(i,p*(j-1)+1:p*j) ) / p

	     im_stat(i,j) = 1.d0*sum( im_stat(i,p*(j-1)+1:p*j) ) / p
	     ndens_stat(i,j) = 1.d0*sum( ndens_stat(i,p*(j-1)+1:p*j) ) / p

             gss_stat(i,j)=1.d0*sum( gss_stat(i,p*(j-1)+1:p*j) ) / p

             nn_szsz_stat(i,j)=1.d0*sum( nn_szsz_stat(i,p*(j-1)+1:p*j) ) / p

	  enddo

	enddo

!c	print*,'after: ',Z_b, Z_b_max

! ------------  where to start from --------------------
	istart_num=3
	allocate( istart(1:nfls, 1:istart_num) )
	
	do i=1,nfls
	  istart(i,:)=(/ 1,b_n(i)/4,b_n(i)/2 /)     ! everyth, drop 1st quarter, 1st half
	enddo


	  do klm=1,istart_num


! -------------  glue all the files -----------------
	print*; print*; 
	print*,'-------------------------------------------------' 
	print*


	allocate(empty(1:nfls)); empty=.false.

  	b_n_all = 0 !sum(b_n)
	do i=1,nfls
	   tmp = b_n(i) - istart(i,klm) + 1
	   b_n_all = b_n_all + tmp
	   if(tmp==0)then; print*,'No stat. left in #', i
	                   empty(i)=.true. 
	   endif
	enddo
	if(b_n_all==0)then; print*,'No stat. left whatsoever...'
	                    stop
	endif


	allocate( PE_all(b_n_all), KE_all(b_n_all), im_all(b_n_all) )
	allocate( PE_2_all(b_n_all),  PEKE_all(b_n_all))
	allocate( ndens_all(b_n_all), ene_sc(b_n_all),tef(b_n_all) )
	allocate( gss_all(b_n_all),chi_all(b_n_all) )
        allocate( nn_szsz_all(b_n_all) )

	i_all=0; Z_all = 0

	do i=1,nfls;  
	
	   if(empty(i))cycle
	  
	  do j=istart(i,klm),b_n(i)
	     i_all = i_all + 1

	     ndens_all(i_all) = ndens_stat(i,j)

	     eef = (3.d0*pi*pi*ndens_all(i_all))**(2.d0/3.d0)
	     npart = ndens_all(i_all)*Nsite(i)

	     PE_all(i_all) = -PE_stat(i,j) / Nsite(i)

	     KE_all(i_all) = -2.d0* KE_stat(i,j) / Nsite(i)

	     ene_sc(i_all) = KE_all(i_all) + PE_all(i_all)

	     tef(i_all) = 1.d0 / beta(1) / eef

		 gss_all(i_all) = gss_stat(i,j)
		 chi_all(i_all) = gss_all(i_all)/( (N(1,1)**(1.d0+H_eta))  )

             nn_szsz_all(i_all) = nn_szsz_stat(i,j)    		 

	  enddo

	enddo





	Z_all = Z_b(1) * i_all

	step_all = sum(step); 



! ---------------  prntout ------------------
	
	print*,d, N(1,:)
	print*,'beta = ', beta(1), ' U = ', U(1), ' mu = ', mu(1)
	print*
	

	print*,'MC step (mln) = ', step_all/1.d6
	print*,' Z  = ',Z_all/1.d6


!--- pot. energy -----------------------------------
	print*; print*,'doing PE: '
	call mrg(PE_all(1:b_n_all),b_n_all,1)


      
!--- kin. energy -----------------------------------
	print*; print*,'doing KE: '
	call mrg(KE_all(1:b_n_all),b_n_all,1)



!---------  scaled energy ------------------------
	print*; print*,'doing total energy per lattice site: '
	call mrg(ene_sc(1:b_n_all),b_n_all,1)
	

!--- ndens
	print*; print*,'doing denity: '
	call mrg(ndens_all(1:b_n_all),b_n_all,1)

!------Magnetic structure factor -----------------
	print*; print*,'doing NMSF: '
	call mrg(gss_all(1:b_n_all),b_n_all,1)

!------ Magnetic Susceptibility-----------------
	print*; print*,'doing magnetic susceptibility: '
	call mrg(chi_all(1:b_n_all),b_n_all,1)

!------  nn <SzSz> -----------------
	print*; print*,'doing <SzSz>@nn : '
	call mrg(nn_szsz_all(1:b_n_all),b_n_all,1)



! ---- nn <Sz:Sz>
        call bSTAT(nn_szsz_g,nn_szsz_g_howmany,1,av,err)
        print 999, av, err
 999    format(4x,'<Sz:Sz>@nn = ',g12.5,4x,' +/- ',g12.5, '   (estimated err)')
        print*

!--- T / E_F
!	print*; print*,'doing T / E_F: '
!	call mrg(tef(1:b_n_all),b_n_all,1)



!--- condensate
!c	print*; print*,'doing g_im(w=0,k=0): '
!c	call mrg(im_all(1:b_n_all),b_n_all,1)



	deallocate(PE_all,KE_all,ndens_all,im_all, empty, ene_sc,tef)
	deallocate(PE_2_all, PEKE_all)
	deallocate(gss_all, chi_all, nn_szsz_all)

		enddo   ! klm



	contains  !========================================

!---------------------
!--- Read statistics
!---------------------
	subroutine rd_stat(i,fname)
	integer :: i
	character*50 :: fname
	character*6  :: cvoid
	double precision :: phi(1:d)


	open(1,file=trim(fname))
	  read(1,*)ddd(i),N(i,:)
	  
!---------- GF -------------------------------------------
	if(i==1) then
		allocate( gf_r(nfls,0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1) )
		allocate( gf_i(nfls, 0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
		allocate( zgf(nfls, 0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1) )

		allocate( gf_r_all(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
		allocate( gf_i_all(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
	
		allocate( gf_k_r(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
		allocate( gf_k_i(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
	
		allocate( sigma_k_r(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
		allocate( sigma_k_i(0:N(1,1)/2, 0:N(1,2)/2, 0:G_nfreq-1))
				
		
		
		allocate(BF_index_vector_p(1:d,0:S_nbf_p**d-1))
		
		allocate( sigma_bfe_r(0:S_nbf_p**d-1, 0:G_nfreq-1))
		allocate( sigma_bfe_i(0:S_nbf_p**d-1, 0:G_nfreq-1))

		
	endif

	
!----------------------------------------------------------
	  read(1,*)beta(i),U(i),mu(i)
	  read(1,*)phi
	  read(1,*)step(i), Z(i)
	  read(1,*)Z_b(i), b_n(i), i_b(i)
!	print*, 'i=', i, 'z_b=', Z_b(i)
!	print*, 'b_n=', b_n(i)
!	print*, 'step', step(i)
!	print*, 'Z=', Z(i)
	  read(1,*)PE(i)
	  read(1,*)PE_stat(i,1:b_n(i));   
	  read(1,*)KE(i)
	  read(1,*)KE_stat(i,1:b_n(i))
	  read(1,*)ndens(i)
	  read(1,*)ndens_stat(i,1:b_n(i))
	  read(1,*)im(i)
	  read(1,*)im_stat(i,1:b_n(i))
	  read(1,*)PE_2(i)
	  read(1,*)PE_2_stat(i,1:b_n(i))
	  read(1,*)PEKE(i)
	  read(1,*)PEKE_stat(i,1:b_n(i))
	  read(1,*)cvoid
	  read(1,*)g_uu(i,0:N(i,1)-1, 0:N(i,2)-1)
	  read(1,*)z_uu(i,0:N(i,1)-1, 0:N(i,2)-1)
	  read(1,*)g_ud(i,0:N(i,1)-1, 0:N(i,2)-1)
	  read(1,*)z_ud(i,0:N(i,1)-1, 0:N(i,2)-1)
	  read(1,*)gss(i)
	  read(1,*)gss_stat(i,1:b_n(i))
          read(1,*)nn_szsz(i)
          read(1,*)nn_szsz_stat(i,1:b_n(i))
	  read(1,*)gf_r(i, 0:N(i,1)/2, 0:N(i,2)/2, 0:G_nfreq-1)
	  read(1,*)gf_i(i, 0:N(i,1)/2, 0:N(i,2)/2, 0:G_nfreq-1)
	  read(1,*)zgf(i, 0:N(i,1)/2, 0:N(i,2)/2, 0:G_nfreq-1)	
	close(1)

	Nsite(i)=N(i,1)*N(i,2)

	PE_stat(i,1:b_n(i)) = PE_stat(i,1:b_n(i)) / Z_b(i)
	KE_stat(i,1:b_n(i)) = KE_stat(i,1:b_n(i)) / Z_b(i)
	
	PE_2_stat(i,1:b_n(i)) = PE_2_stat(i,1:b_n(i)) / Z_b(i)
	PEKE_stat(i,1:b_n(i)) = PEKE_stat(i,1:b_n(i)) / Z_b(i)

	ndens_stat(i,1:b_n(i)) = ndens_stat(i,1:b_n(i)) / Z_b(i)
	im_stat(i,1:b_n(i)) = im_stat(i,1:b_n(i)) / Z_b(i)
	gss_stat(i,1:b_n(i)) = gss_stat(i,1:b_n(i)) / Z_b(i)

        nn_szsz_stat(i,1:b_n(i))= nn_szsz_stat(i,1:b_n(i)) / Z_b(i)

	end subroutine rd_stat


!-------------------------------
!--- Analyze block statistics
!-------------------------------
	subroutine bSTAT(arr,n,Zb,av,err)
	integer              :: n, Zb
	real*8, dimension(1:n) :: arr
	real*8               :: av, err

	real*8 :: av2

	av  = sum( arr(1:n) )/Zb/n
	av2 = sum( arr(1:n)**2 )/Zb/Zb/n

				!av2 = av2 + (arr(j)/Zb)**2

	err = sqrt( av2 - av*av ) / sqrt(1.d0*n)


	end subroutine bSTAT


!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
	subroutine mrg(arr,n,Zb)
	integer, intent(in)              :: n, Zb
	real*8, dimension(1:n), intent(in) :: arr

	real*8  :: av, err, arr1(1:n)
	integer :: i, n1, zb1


	zb1 = zb; 	arr1(1:n) = arr(1:n); n1=n

	print*,'-----------'
	
	do;
	
! emit
	call bSTAT(arr1,n1,zb1,av,err)
      print 777, av, err,n1
 777  format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
	if(n1<3)exit

! merge
	n1=INT(n1/2); zb1=zb1*2
	do i=1,n1
	    arr1(i) =  arr1(2*i-1) + arr1(2*i)
	enddo

	enddo

	print*,'------------'; print*; 


	end subroutine mrg


	subroutine fourier_GF

      integer :: i, mx, my !, mz
      integer :: nt, nx, ny !, nz
	  integer :: nxr, nyr !, nz
      double precision :: phase
      double precision :: pp(d), factor
      real*8 :: kx,ky !,kz


	do i=1,d; pp(i)=4.d0*asin(1.d0)/N(1,i); enddo


	gf_k_r=0.d0; gf_k_i=0.d0

! sum over momenta 1st 
!	do mz=0, N(1,3)-1      
	do my=0, N(1,2)/2;     
      do mx=0, N(1,1)/2;    
	  
        kx=pp(1)*mx
        ky=pp(2)*my
!        kz=pp(3)*mz 



! coordinates 
          !do nz = 0, Ntab-1; 
		  do ny = 0, N(1,2)-1; do nx = 0, N(1,1)-1
	       !phase=pp(1)*nx*mx +pp(2)*ny*my+pp(3)*nz*mz
             phase = kx*nx + ky*ny ! + kz*nz
			 factor=cos(phase)
			 nxr=min(nx,N(1,1)-nx)
			 nyr=min(ny,N(1,2)-ny)
				
	       do nt=0, G_nfreq-1

			gf_k_r(mx,my, nt)=gf_k_r(mx,my, nt) + gf_r_all(nxr,nyr,nt)*factor
			gf_k_i(mx,my, nt)=gf_k_i(mx,my, nt) + gf_i_all(nxr,nyr,nt)*factor


		   enddo                  ! tau
	    enddo; enddo !; enddo  ! coordinates
	enddo; enddo ! ; enddo      ! momentum
!------------------------
	
	fname='gf_k'//trim(fsuffix)//'.io'
	open(1,file=trim(adjustl(fname)))
	   write(1,*) gf_k_r
	   write(1,*) gf_k_i
	close(1)
		
	
	
	end subroutine fourier_GF
	
	
	subroutine bfe_Sigma
	
      integer :: i, mx, my !, mz
      integer :: nbf, nf
      double precision :: phase
      double precision :: pp(d), factor, p0(1:d)
	  double precision :: FR(0:S_nbf_p**d-1, 0:G_nfreq-1, 0:N(1,2)/2)
	  double precision :: FI(0:S_nbf_p**d-1, 0:G_nfreq-1, 0:N(1,2)/2)
	  double precision :: F0R(0:N(1,2)/2)
	  double precision :: F0I(0:N(1,2)/2)
	  double precision :: F1R(0:N(1,2)/2)
	  double precision :: F1I(0:N(1,2)/2)
	  double precision :: sum

	print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	print*, 'Number of momentum basis functions per dimension, S_nbf_p=', S_nbf_p
	print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'



	do i=1,d; pp(i)=2.d0*pi/N(1,i); enddo
	
	FR=0.d0
	FI=0.d0
		
	do nf=0, G_nfreq-1
		do nbf=0, S_nbf_p**d-1
			
			do my=0, N(1,2)/2;     
				do mx=0, N(1,1)/2; 
					p0(1)=pp(1)*mx
					p0(2)=pp(2)*my
					
					F0R(mx)= BF_p_c(nbf, p0)*sigma_k_r(mx, my, nf)
					F0I(mx)= BF_p_c(nbf, p0)*sigma_k_i(mx, my, nf)
					
				enddo
				
				F1R(my)=Simpson(N(1,1)/2, F0R)
				F1I(my)=Simpson(N(1,1)/2, F0I)
			enddo
			
			sigma_bfe_r(nbf,nf)=Simpson(N(1,2)/2, F1R)
			sigma_bfe_i(nbf,nf)=Simpson(N(1,2)/2, F1I)
		enddo
	enddo
		

	sigma_bfe_r=sigma_bfe_r*16.d0/(N(1,1)*N(1,2))
	sigma_bfe_i=sigma_bfe_i*16.d0/(N(1,1)*N(1,2))
	
      
	end subroutine bfe_Sigma
	
	!! simpson integration (not normalized)
	double precision function Simpson(N, F)
		integer, intent(in) :: N    
		double precision, intent(in) :: F(0:N)
		!double precision, intent(in) :: x(0:N)
		double precision :: sum
		integer :: i
		
		if(mod(N,2)/=0) then
			print*, 'Simpson integration cannot work with odd N!'
			print*, 'N=', N
			stop
		endif
		
		sum=0.d0
		i=0
		do while (i/=N)
			sum=sum + (F(i) + 4.d0*F(i+1) + F(i+2))/6.d0
			i=i+2 
		enddo
	
		Simpson=sum
		return
		
	end function Simpson
	
	
	
	
	
	
!!**************************************************************************************
!!       Basis Function Expansions 
!!**************************************************************************************

!! conjugate function
double precision function BF_p_c(index, p0)
	implicit none
	integer, intent(in) :: index
	double precision, intent(in) :: p0(1:d)

	BF_p_c=1.d0
	do i=1, d
		if(BF_index_vector_p(i, index)==0) then
			BF_p_c=BF_p_c*0.5d0
		endif
		BF_p_c=2.d0*BF_p_c*cos(BF_index_vector_p(i, index)*p0(i)) 
	enddo
	
	return
end function BF_p_c 

subroutine tabulate_basis_function_index_vectors
  implicit none
  integer :: fp, fp1, fp0
  
	do fp=0, S_nbf_p**d-1
		
		fp1=fp
		do i=1, d
		fp0=fp1/(S_nbf_p**(d-i))
		fp1=fp1-fp0*(S_nbf_p**(d-i))
		BF_index_vector_p(d+1-i,fp)=fp0
		enddo
	enddo
  
end subroutine tabulate_basis_function_index_vectors

!!**************************************************************************************
	
	
	subroutine plot_sigma_k
	real*8 :: kx, ky, pp(d), sr, si, nrm, sum1, sum2, sum3
	integer :: mx, my, m1, m2, k, i 
	real*8 :: g0_1_r, g0_1_i, xi_k, eps, gr, gi
	real*8 :: KR(0:G_nfreq-1), KI(0:G_nfreq-1)
	real*8 :: K0R(0:G_nfreq-1), K0I(0:G_nfreq-1), K0
	
	
	do i=1,d; pp(i)=4.d0*asin(1.d0)/N(1,i); enddo
	
	KR=0.d0
	KE=0.d0
	K0R=0.d0
	K0I=0.d0
	K0=0.d0
	nrm=0.d0
	
	do my=0, N(1,2)/2;     
      do mx=0, N(1,1)/2;    
        kx=pp(1)*mx
        ky=pp(2)*my
		eps = -2.d0*( cos(kx) + cos(ky) ) - mu(1)
		g0_1_r=-eps
		
		do k=0, G_nfreq-1
			xi_k = 2.d0*pi*(k+0.5d0)/beta(1)
			g0_1_i=xi_k
			
			gr=gf_k_r(mx,my,k)
			gi=gf_k_i(mx,my,k)
			
			sigma_k_r(mx,my,k)= g0_1_r - gr/(gr*gr+gi*gi)
			sigma_k_i(mx,my,k)= g0_1_i + gi/(gr*gr+gi*gi)
						
		enddo

		
	  enddo
	enddo
	
	
	!nxr=min(nx,N(1,1)-nx)
	
! ---------- do kinetic energy -----------------------------------
!
	do my=0, N(1,2)-1;     
      do mx=0, N(1,1)-1;    
        kx=pp(1)*mx
        ky=pp(2)*my
		eps = -2.d0*( cos(kx) + cos(ky) ) !- mu(1)
		g0_1_r=-eps
		nrm=nrm+1.d0
	
		do k=0, G_nfreq-1
			xi_k = 2.d0*pi*(k+0.5d0)/beta(1)
			g0_1_i=xi_k
			
			m1=min(mx,N(1,1)-mx)
			m2=min(my,N(1,2)-my)
			
			gr=gf_k_r(m1,m2,k)
			gi=gf_k_i(m1,m2,k)
			
			KR(k)=KR(k)+eps*gr
			KI(k)=KI(I)+eps*gi
			K0R(k)=K0R(K)-eps*eps/(eps*eps+xi_k*xi_k)
			K0I(k)=K0I(K)-eps*xi_k/(eps*eps+xi_k*xi_k)
			
		enddo
		
		K0=K0+eps/(1.d0 + exp(beta(1)*eps))
		
	  enddo
	enddo
!-----------------------------------------------------------------
	
	
	
	K0=K0/nrm
	KR=KR/nrm
	KI=KI/nrm
	K0R=K0R/nrm
	K0I=K0I/nrm
	
	fname='Sigma_k_'//trim(fsuffix)//'.io'
	open(1,file=trim(adjustl(fname)))
	   write(1,*) sigma_k_r
	   write(1,*) sigma_k_i
	close(1)
		
	call bfe_Sigma
	
	fname='Sigma_bfe_'//trim(fsuffix)//'.io'
	open(1,file=trim(adjustl(fname)))
	   write(1,*) sigma_bfe_r
	   write(1,*) sigma_bfe_i
	close(1)
	
	fname='Sigma0_k_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   k=0
	   do my=0, N(1,2)/2; mx=my
		kx=pp(1)*mx
		write(1,*) kx, sigma_k_r(mx,my,k), sigma_k_i(mx,my,k)
	   enddo
	close(1)
	
	fname='Sigma63_k_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   k=63
	   do my=0, N(1,2)/2; mx=my
		kx=pp(1)*mx
		write(1,*) kx, sigma_k_r(mx,my,k), sigma_k_i(mx,my,k)
	   enddo
	close(1)
	
	fname='G63_k_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   k=63
	   do my=0, N(1,2)/2; mx=my
		kx=pp(1)*mx
		write(1,*) kx, gf_k_r(mx,my,k), gf_k_i(mx,my,k)
	   enddo
	close(1)
	
	fname='Sigma_loc_k_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   do k=0, G_nfreq-1
		sr=0.d0
		si=0.d0
		nrm=0.d0
		do my=0, N(1,2)-1; do mx=0, N(1,1)-1
			m1=min(mx,N(1,1)-mx)
			m2=min(my,N(1,2)-my)			
			sr=sr+sigma_k_r(m1,m2,k)
			si=si+sigma_k_i(m1,m2,k)
			nrm=nrm+1.d0
		enddo; enddo
		write(1,*) k, sr/nrm, si/nrm
	   enddo
	close(1)
	
	fname='G_loc_k_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   do k=0, G_nfreq-1
		sr=0.d0
		si=0.d0
		nrm=0.d0
		do my=0, N(1,2)-1; do mx=0, N(1,1)-1
			m1=min(mx,N(1,1)-mx)
			m2=min(my,N(1,2)-my)	
			sr=sr+gf_k_r(m1,m2,k)
			si=si+gf_k_i(m1,m2,k)
			nrm=nrm+1.d0
		enddo; enddo
		write(1,*) k, sr/nrm, si/nrm
	   enddo
	close(1)
	
	
	fname='K_xi.dat'
	open(1,file=trim(adjustl(fname)))
	   sum1=0.d0; sum2=0.d0; sum3=0.d0
	   do k=0, G_nfreq-1
		write(1,*) k, KR(k), KI(k)
		sum1=sum1+K0R(k)
		sum2=sum2+KR(k)
		sum3=sum3-2.d0*(gf_r_all(0,1,k) + gf_r_all(1,0,k))
	   enddo
	   sum1=4.d0*sum1/beta(1)
	   sum2=4.d0*sum2/beta(1)
	   sum3=4.d0*sum3/beta(1)
	   print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	   print*, 'From the Greenz function:'
	   print*, ' '
	   print*, 'KE0=', 2.d0*K0
	   print*, '2 sum_n_k e_k G0(xi_n) =', sum1
	   print*, 'KE =',(sum2-sum1 + 2.d0*K0)
	   print*, '2 sum_n_k e_k G(xi_n)  =', sum2
	   print*,  ' sum G_{i,i+1}(xi_n)  =', sum3
	   print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	close(1)


	fname='K0_xi.dat'
	open(1,file=trim(adjustl(fname)))
	   do k=0, G_nfreq-1
		write(1,*) k, K0R(k), K0I(k)
	   enddo
	close(1)

	
	fname='KmK0_xi_DDMC.dat'
	open(1,file=trim(adjustl(fname)))
	   do k=0, G_nfreq-1
		write(1,*) k, KR(k)-K0R(k), KI(k)-K0I(k)
	   enddo
	close(1)
	
	
	endsubroutine plot_sigma_k
	
	end
