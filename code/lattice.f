!------------------------------------
!--- d-dimensional lattice with PBC
!------------------------------------
    module lattice
    implicit none

    public :: d, dd, N, N2, Nsite, Nbond, ass, back, x, ASSA, distance_v, num_neighb

    private

    integer :: d, dd
    integer, allocatable :: N(:), N2(:)
    integer :: Nsite      ! Total number of sites
    integer :: Nbond      ! Total number of bonds

!------------------------------------------------------
!     Sites and Associations
!------------------------------------------------------

      integer, allocatable :: ass(:,:)
      integer, allocatable :: back(:)
      integer, allocatable :: x(:,:)  ! coordinates x(1:d,site)


!     The array ass(...) specifies the nearest-neighbor sites in both positive and
!     negative i-direction. For example, in a 2-D lattice, like the one below,
!     the associations to the io site are:  
!     (here +(-)1 = positive (negative) x-direction, 
!     +(-)2 = positive (negative) y-direction)
      
!     ass(+1,i0)=i1;    ass(-1,i0)=i3;    ass(+2,i0)=i4;    ass(-2,i0)=i2
      
!                            ...    
!                             |  
!                     ... -- i2 -- ...
!                       |     |     |
!                ...-- i3 -- i0 -- i1 --.... 
!                       |     |     |
!                     ... -- i4 -- ... 
!                             |  
!                            ...
!           
!     * The array back(..) specifies the 'opposite direction', e.g. in 3D
!           back(1) = 4, back(4)=1 (positive x is 1, negative x is 4), and so forth
!

    contains


!--------------------------------------------
!--- Arranges associations between sites
!--------------------------------------------
      subroutine ASSA  

      integer :: site, site1, i, i1, i2 !, i3 
      integer :: ic(d), NN(d+1) 
!     ic(i) is the i-coordinate of a site, 
!     the relation between site indexes and coordinates reads:
!     site = 1 + (ic(1)-1) + N(1)*(ic(2)-1) +...
!               + N(1)*...*N(d-1)*(ic(d)-1)

      allocate (ass(dd,Nsite),back(dd),x(1:d,1:Nsite))
      Nbond = d*Nsite   ! number of bonds

      NN(1)=1; do i=2,d+1; NN(i)=NN(i-1)*N(i-1); enddo
       
      do i=1,d; back(i)=i+d; back(i+d)=i; enddo 
      
      ic=1; ic(1)=0
      DO site=1, Nsite

!------------ Coordinates for site ----------
         i1=1 
         DO
         if (ic(i1) < N(i1)) then
            ic(i1)=ic(i1)+1
            DO i2=1,i1-1; ic(i2)=1; END DO
            EXIT
         else; i1=i1+1;  end if 

         END DO

!-------------------------------------------------------


         DO i=1, d
            if (ic(i) < N(i)) then
               site1=site+NN(i)
            else
               site1=site+NN(i)-NN(i+1)
            end if

            ass(i,site)=site1
            ass(back(i),site1)=site
            x(:,site)=ic(:)
                 
         END DO
      END DO
      
      end subroutine ASSA


!----------------------------------------------
!--- Distance between two coord vectors, incl BC
!----------------------------------------------
    function distance_v(x1, x2)
    implicit none
    integer, dimension(1:d) :: distance_v
    integer, dimension(1:d) :: x1, x2
    integer :: j, xx(1:d)

    j=1; xx(j) = abs(x1(j)-x2(j)); xx(j) = min(xx(j), N(j) - xx(j))
    j=2; xx(j) = abs(x1(j)-x2(j)); xx(j) = min(xx(j), N(j) - xx(j))

    distance_v = xx

    end function distance_v


!----------------------------------------------
!--- Number of neighbors
!----------------------------------------------
    function num_neighb(site)
    implicit none
    integer :: num_neighb
    integer, intent(in) :: site

    num_neighb = dd

    end function num_neighb

    end module lattice
