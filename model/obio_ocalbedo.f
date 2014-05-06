#include "rundeck_opts.h"      

      module obio_ocalbedo_mod

      implicit none

      integer, parameter :: nlt=33                ! number of spectral channels
      real*8, dimension(nlt), private :: wfac
      real*8, dimension(nlt), protected :: aw, bw ! absorption,scattering coefficients of water

      integer, dimension(nlt), protected :: lam   ! wavelength in nm
      logical, private :: initialized=.false.

      contains

      subroutine obio_ocalbedo(wind,solz,bocvn,xocvn,chl,
     .                         rod,ros,hycgr,i,j)

***********************************************************************
***** this routine is used by both atmosphere and ocean at each (i,j)
***** where implicitly it is assumed that
***** wind,solz,bocvn,xocvn,chl are in the atmos gird if hycgr=.false. 
***** wind,solz,                are in the ocean gird if hycgr=.true.
***** when this routine is called from within ocean, it does not compute
***** albedo coefficients. Those are computed on the atmos grid.
***** when this routine is called from within ocean, it does not pass
***** reflectances
***********************************************************************

c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
c  Derive surface reflectance as a function of solz and wind
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
      USE CONSTANT, only : radian

      implicit none

      real*8, intent(in)  :: wind, solz, chl
      logical, intent(in) :: hycgr
      integer, intent(in) :: i,j
      real*8, dimension(6), intent(out) :: bocvn, xocvn
      real*8, dimension(:), intent(out) :: rod, ros

      integer nl
      real*8 cn,rof,rosps,rospd,rtheta
      real*8 sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real*8 tanrpls,sinp,tanp,a,b

      real*8 :: sunz
 
      real*8 :: sum1, sum2, part_sum
      logical :: obio_reflectance, res
      real*8, dimension(nlt) :: refl

!!!!!!!!!!  Boris' part !!!!!!!!!!!!!!!!!
!@sum Those are the weights which were obtained by normalizing solar flux
!@sum for the Lamda's range 0 - 4000 nm. We used Landau fitting function for it.
!@sum They are used to be used for getting 6 band approximation based on 33
!@sum Watson Gregg band calculations
!@sum band_6(j) = Sum(band_33(i)*weight(i))/Sum(part_sum(i))
C**** WHY IS WEIGHT ONLY DECLARED TO BE 31 AND NOT 33?
      real*8 :: weight(31) = (/0.0158378,0.0201205,0.0241885,0.0277778,
     .                       0.0307124,0.0329082,0.0343586,0.0351143,
     .                       0.0352609,0.0349008,0.0341389,0.0330742,
     .                       0.0317941,0.0303725,0.0288696,0.0273329,
     .                       0.0500921,0.0643897,0.0686573,0.0532013,
     .                       0.0416379,0.0330341,0.0265929,0.0217156,
     .                       0.0179725,0.0150596,0.0127618,0.0158128,
     .                       0.0232875,0.0313132,0.0184843/)

C**** gband is the distribution of the 31 bands for the 6 band GISS code
      integer :: gband(7) = (/ 1, 18, 19, 23, 26, 30, 32 /)
c      real*8 :: part_sum(6) = (/0.526854,0.0643897,0.196531,0.066281,
c     .                        0.066922,0.0497974/)
        real*8 :: lam8(nlt)
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!

      real*8 :: roair, rn
      integer :: ngiss

      call init
      rn = 1.341d0  ! index of refraction of pure seawater
      roair = 1.2D3 ! density of air g/m3  SHOULD BE INTERACTIVE?

      sunz=acos(solz)/radian  !in degs

c  Foam and diffuse reflectance
      if (wind .gt. 4.0) then
        if (wind .le. 7.0) then
          cn = 6.2D-4 + 1.56D-3/wind
          rof = roair*cn*2.2D-5*wind*wind - 4.0D-4
        else
          cn = 0.49D-3 + 0.065D-3*wind
          rof = (roair*cn*4.5D-5 - 4.0D-5)*wind*wind
        endif
        rosps = 0.057d0
      else
        rof = 0.0
        rosps = 0.066d0
      endif
      
c  Direct
c   Fresnel reflectance for sunz < 40, wind < 2 m/s
      if (sunz .lt. 40.0 .or. wind .lt. 2.0) then
        if (sunz .eq. 0.0) then
          rospd = 0.0211d0
        else
          rtheta = sunz*radian
          sintr = sin(rtheta)/rn
          rthetar = asin(sintr)
          rmin = rtheta - rthetar
          rpls = rtheta + rthetar
          sinrmin = sin(rmin)
          sinrpls = sin(rpls)
          tanrmin = tan(rmin)
          tanrpls = tan(rpls)
          sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
          tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
          rospd = 0.5*(sinp + tanp)
        endif
      else
       !Empirical fit otherwise
        a = 0.0253d0
        b = -7.14D-4*wind + 0.0618d0
        rospd = a*exp(b*(sunz-40.0))
      endif

c  Reflectance totals
      do nl = 1,nlt
        ros(nl) = rosps + rof*wfac(nl)
        rod(nl) = rospd + rof*wfac(nl)
      enddo

      !lam is integer, lam8 is real8
      lam8=float(lam)

      if (hycgr) return  !we do not compute albedo coefs
                         !from within ocean, but from atmos
      
!!!!!!!!!! Boris' part !!!!!!!!!!!!!!!!!
C**** get chlorophyll term
      ! function obio_reflectance calculates reflectance 
      ! as a function of chl and wavelength (lam)

      res = obio_reflectance(refl,chl,lam8,nlt,i,j)
 
!  transition between band33 and band6 approximation

! loop over giss radiation bands
      do ngiss=1,6
        sum1 = 0.0
        sum2 = 0.0
        do nl=gband(ngiss), gband(ngiss+1)-1 
          if (refl(nl).lt.0) refl(nl) = 0.0
          ros(nl) = ros(nl) + refl(nl)
          sum1 = sum1+weight(nl)*rod(nl)
          sum2 = sum2+weight(nl)*ros(nl)
        enddo
        part_sum=sum(weight(gband(ngiss):gband(ngiss+1)-1))
        xocvn(ngiss) = sum1/part_sum
        bocvn(ngiss) = sum2/part_sum

      if (xocvn(ngiss).ge.1. .or. bocvn(ngiss).ge.1.) then
         print*, 'XOCVN/BOCVN greater than 1 at ngiss,i,j=',
     .                ngiss,i,j
         do nl=gband(ngiss), gband(ngiss+1)-1
           write(*,'(i5,4e12.4)')
     .     nl,rod(nl),ros(nl),refl(nl),weight(gband(ngiss))
         enddo
         stop
      endif
      enddo
    
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!
      return
      end subroutine obio_ocalbedo

!=======================================================================
      subroutine init

      use filemanager, only: openunit, closeunit
      implicit none

      real*8 :: saw, sbw
      real*8 :: b0, b1, b2, b3, a0, a1, a2, a3, expterm, tlog, fac, rlam
      integer :: nl,ic , iu_bio, lambda
      character title*50
      data a0,a1,a2,a3 /0.9976d0, 0.2194d0,  5.554d-2,  6.7d-3 /
      data b0,b1,b2,b3 /5.026d0, -0.01138d0, 9.552d-6, -2.698d-9/

      if (initialized) then
        return
      else
        initialized=.true.
      endif
      call openunit('cfle1',iu_bio,.false.,.true.)
      do ic = 1,6
        read(iu_bio,'(a50)')title
      enddo
      do nl = 1,nlt
        read(iu_bio,20) lambda,saw,sbw
        lam(nl) = lambda
        aw(nl) = saw
        bw(nl) = sbw
        if (lam(nl) .lt. 900) then
          expterm = exp(-(aw(nl)+0.5*bw(nl)))
          tlog = dlog(1.0D-36+expterm)
          fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
          wfac(nl) = max(0d0,min(fac,1d0))
        else
          rlam = float(lam(nl))
          fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
          wfac(nl) = max(fac,0d0)
        endif
      enddo
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)
      return

      end subroutine init
!=======================================================================

      end module obio_ocalbedo_mod
