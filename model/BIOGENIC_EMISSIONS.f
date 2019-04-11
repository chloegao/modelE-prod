#include "rundeck_opts.h"

      module biogenic_emis
      integer, parameter :: npolynb=20, ntype=16
!@dbparam base_isopreneX factor to tune the base isoprene emissions
!@+ globally when #defined BIOGENIC_EMISSIONS
      real*8 :: base_isopreneX=1.d0
!@var baseisop (kg C / m2 /s ) base isoprene emission
      real*8, allocatable, dimension(:,:,:) :: baseisop
      real*8, parameter :: isopcoeff(npolynb) = (/
     *     -1.86E-01, 2.19E+00,  2.12E+00, -2.43E-01, -4.72E+00,
     *      1.05E+01, 3.77E-01, -3.05E+00,  3.05E-01,  5.16E-01,
     *      2.72E+00,-4.82E+00, -9.59E-04, -1.29E+00,  9.37E-01,
     *     -3.31E-01, 1.07E+00, -7.59E-02, -3.01E-01, -4.07E-01 /) 

      end module biogenic_emis


      subroutine alloc_biogenic_emis(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use biogenic_emis,     only : baseisop
      use tracers_DRYDEP,    only : nent
      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: I_0,I_1,J_0,J_1

      call getDomainBounds(grid)
      I_0=grid%I_STRT
      I_1=grid%I_STOP
      J_0=grid%J_STRT
      J_1=grid%J_STOP

      allocate( baseisop(I_0:I_1,J_0:J_1,NENT) )

      end subroutine alloc_biogenic_emis


      subroutine isoprene_emission(i,j,itype,emisop,CanTemp)

      use biogenic_emis
      use rad_com,        only : cosz1, cfrac
      use pbl_drv,        only : t_pbl_args
      use constant,       only : rgas
      use tracers_drydep, only : xylai, ijuse, nent
      use constant,       only : tf
      use ghy_com,        only : fearth
    
      implicit none
 
      real*8, intent(in) :: CanTemp !@var CanTemp canopy temperature (C)
      type (t_pbl_args) :: pbl_args
      integer :: inveg
      integer, intent(in) :: i,j,itype  
      real*8 :: tlai,embio,clight,biofit,tcorr,sfdens,tmmpk
      real*8, intent(out) :: emisop

      select case(itype)

      case(1:3)      ! ocean, ocean ice, landice:

        emisop=0.d0  ! no emissions

      case(4)        ! land:

! Temperature of canpoy 

        ! estimated surface density
!        sfdens=100.d0*pbl_args%psurf/(rgas*pbl_args%tgv) 
        tmmpk = CanTemp + tf

        emisop=0.d0
        tlai=0.d0

        do inveg=1,nent
          tlai=tlai+xylai(i,j,inveg)*baseisop(i,j,inveg)
        end do

! Light correction
        if ((cosz1(i,j) > 0.d0).and.(tlai > 0.d0)) then ! Only calculate
        ! for grid cell with sunlight and isoprene-emitting vegetation
          embio=0.d0
          do inveg=1,nent
            if (xylai(i,j,inveg)*baseisop(i,j,inveg) > 0.d0) then
              clight=biofit(isopcoeff,xylai(i,j,inveg),
     &               cosz1(i,j),cfrac(i,j))
              embio=embio+baseisop(i,j,inveg)*clight*ijuse(i,j,inveg)
            endif
          end do

! Temperature correction
          if(tmmpk > tf) THEN
            emisop=tcorr(tmmpk)*embio
          else
            emisop=0.d0
          endif

        endif  

      end select

      return                                                          
      end                              


      subroutine rdisopcf                                              
!@sum These polynomial coefficients normally should be read
!@+  in from the file 'isoprene.coef'. I am hardcoding them
!@+  here, as the quickest way to make sure this is ESMF-
!@+  compliant, as I do not suspect we will commit this code.
!@ THIS IS NOW OBSOLETE

      use biogenic_emis

      implicit none

! polynomial coefficients for biogenic isoprene emissions:
! From Guenther.

      if(npolynb /= 20)call stop_model('npolynb problem',255)

      return                                                            
      end                       


      subroutine rdisobase
!@sum These baseline emissions factors normally should be read
!@+  in from the file 'isopemis.table'. I am hardcoding them
!@+  here, as the quickest way to make sure this is ESMF-
!@+  compliant, as I do not suspect we will commit this code.
!@+  Units are atoms C cm^-2 leaf s^-1
!@+  Construct the base emission for each grid box                         
!@+  Output is baseisop in kg C m^-2 s^-1

      use biogenic_emis
      use tracers_DRYDEP,      only : nent
      use ent_drv,             only : map_ent_pfts_to_megan_pfts
      use constant,            only : byavog
      use geom,                only : imaxj
      use domain_decomp_atm,   only : getDomainBounds, grid
      use trchem_shindell_com, only : nMeganPFT      
      use ghy_com,             only : fearth

      implicit none

      integer:: i,j,k,j_0,j_1,i_0,i_1

! Baseline isoprene emissions factors for the 16 Megan types 
      real*8, dimension(0:nMeganPFT), parameter ::  convert = (/0.d0,
     & 600.d0, 1.d0, 3000.d0, 7000.d0, 10000.d0, 7000.d0, 10000.d0, 
     & 11000.d0, 2000.d0, 4000.d0, 4000.d0, 1600.d0, 800.d0, 200.d0, 
     & 50.d0, 1.d0/)

      real*8 :: factor
      real*8, dimension(nent)      :: v_dummy, h_dummy, megan_map
      real*8, dimension(nMeganPFT) :: v_m_dummy

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &               I_STRT=I_0, I_STOP=I_1)

      if(nent/=18) call stop_model('Number of Ent types incorrect',255)

! Isoprene is traced in terms of equivalent C atoms.
! Compute the baseline ISOPRENE emissions, which depend on veg type. 

! Convert from microgram/m^2/hr to C/cm^2 leaf/s
      factor = (1.d-9)/3600.d0

      do J=J_0,J_1
        do I=I_0,imaxj(J)
            call map_ent_pfts_to_megan_pfts(v_dummy, h_dummy,
     &      v_m_dummy, i, j, megan_map)
            if (fearth(i,j)>0.d0) then
               do k=1,nent
                   baseisop(i,j,k)=convert(megan_map(k))*
     &             factor*base_isopreneX
               enddo ! k
            endif ! fearth>0
        enddo     ! i
      enddo       ! j

      return
      end

! Local Air temperature correction for isoprene emissions:
! From Guenther et al. 1992

      real*8 function tcorr(temp)

      implicit none

      real*8, intent(in) :: temp
      real*8, parameter :: r=8.314d0, ct1=9.5d4, ct2=2.3d5, t1=3.03d2,
     *     t3=3.14d2 

      tcorr =    exp( ct1/(r*t1*temp)*(temp-t1) ) /
     &   (1.D0 + exp( ct2/(r*t1*temp)*(temp-t3) ))

      return
      end
