#include "rundeck_opts.h"

      module biogenic_emis
      integer, parameter :: npolynb=20, ntype=16, nvegtype=74
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
      use biogenic_emis, only:  baseisop,nvegtype

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: I_0,I_1,J_0,J_1

      call getDomainBounds(grid)
      I_0=grid%I_STRT
      I_1=grid%I_STOP
      J_0=grid%J_STRT
      J_1=grid%J_STOP

      allocate( baseisop(I_0:I_1,J_0:J_1,nvegtype) )

      end subroutine alloc_biogenic_emis


      subroutine isoprene_emission(i,j,itype,emisop,CanTemp)

      use biogenic_emis
      use rad_com, only : cosz1,cfrac
      use pbl_drv, only : t_pbl_args
      use constant, only : rgas
      use tracers_drydep, only : xylai,ijreg,ijuse
      use constant, only : tf
    
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

        do inveg=1,ijreg(i,j)
          tlai=tlai+xylai(i,j,inveg)*baseisop(i,j,inveg)
        end do

! Light correction

        if ((cosz1(i,j) > 0.).and.(tlai > 0.)) then ! Only calculate
        ! for grid cell with sunlight and isoprene-emitting vegetation

          embio=0.d0
          do inveg=1,ijreg(i,j)
            if (xylai(i,j,inveg)*baseisop(i,j,inveg) > 0.) then
              clight=biofit(isopcoeff,xylai(i,j,inveg),
     &               cosz1(i,j),cfrac(i,j))
              embio=embio+baseisop(i,j,inveg)*
     &              clight*ijuse(i,j,inveg)*1.d-3
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
      use tracers_drydep, only : ijreg,ijland
      use constant, only   : byavog
      use geom, only : imaxj
      use domain_decomp_atm, only : getDomainBounds, grid

      implicit none

      integer:: i,j,k,j_0,j_1,i_0,i_1
      real*8, dimension(0:nvegtype-1), parameter ::  converT = (/
     *     0.00E+00, 2.79E+12, 8.71E+11, 2.79E+12, 2.79E+12, 2.79E+12, 
     *     2.79E+12, 2.79E+12, 3.34E+12, 2.79E+12, 2.79E+12, 2.79E+12, 
     *     2.79E+12, 2.79E+12, 2.79E+12, 2.79E+12, 2.79E+12, 2.79E+12, 
     *     2.79E+12, 2.79E+12, 1.67E+12, 1.67E+12, 1.67E+12, 1.39E+12, 
     *     3.34E+12, 6.27E+12, 6.27E+12, 3.34E+12, 0.93E+12, 0.93E+12,
     *     8.71E+11, 8.71E+11, 2.49E+12, 1.39E+12, 2.79E+12, 2.79E+12,
     *     8.71E+11, 8.71E+11, 8.71E+11, 8.71E+11, 0.93E+12, 1.39E+12,
     *     8.71E+11, 2.79E+12, 1.67E+12, 1.67E+12, 3.34E+12, 2.79E+12,
     *     9.41E+12, 2.79E+12, 3.34E+12, 3.34E+12, 2.79E+12, 2.79E+12,
     *     3.34E+12, 2.79E+12, 4.18E+12, 4.18E+12, 1.39E+12, 3.34E+12,
     *     3.14E+12, 3.34E+12, 1.39E+12, 3.34E+12, 8.71E+11, 2.79E+12,
     *     2.79E+12, 2.79E+12, 2.79E+12, 3.34E+12, 2.79E+12, 3.34E+12,
     *     8.71E+11, 2.79E+12 /)

      real*8 :: factor                       

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &               I_STRT=I_0, I_STOP=I_1)

      if(nvegtype /= 74) call stop_model('nvegtype problem',255)

! Isoprene is traced in terms of equivalent C atoms.
! Compute the baseline ISOPRENE emissions, which depend on veg type   
! 12.d-3 is the carbon mol wt. in kg/mole.
! 1.d4 is for cm^2 -> m^2

      factor = 12.d-3*1.d4*byavog      

      do J=J_0,J_1
        do I=I_0,imaxj(J)
          do k=1,ijreg(i,j)
            baseisop(i,j,K) = 
     &      convert(ijland(i,j,k)+1)*factor*base_isopreneX
          enddo                                                       
        enddo                                                          
      enddo                                                             
                                                                        
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
