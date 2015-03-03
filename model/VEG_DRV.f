#include "rundeck_opts.h"


      module veg_drv
!@sum veg_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang, Y. Kim

      implicit none
      private
      save

      public get_vdata, get_cropdata, get_laimaxdata, get_soil_C_total

#ifndef USE_ENT
      public init_vegetation,reset_veg_to_defaults
     &     ,veg_set_cell, veg_save_cell,upd_gh, updveg

      real*8,public :: cosday,sinday
#endif
      contains


#ifndef USE_ENT

      subroutine init_vegetation(redogh,istart)
!@sum initializes vegetation
      use Dictionary_mod
      use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      integer, intent(in) :: istart
      logical, intent(in) :: redogh

      call sync_param( "cond_scheme", cond_scheme)  !nyk 5/1/03
      call sync_param( "vegCO2X_off", vegCO2X_off)  !nyk 3/2/04
      call sync_param( "crops_yr", crops_yr)

      call read_veg_data(redogh,istart)
      call upd_gh
      end subroutine init_vegetation


      subroutine read_veg_data(redogh,istart)
!@sum reads vegetation arrays and rundeck parameters
      use Dictionary_mod
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, AM_I_ROOT
      use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      use veg_com
      use fluxes, only : focean
      use model_com, only : modelEclock
      use ghy_com, only : fearth

      implicit none

      integer, intent(in) :: istart
      logical, intent(in) :: redogh

      INTEGER :: J_1, J_0, J_1H, J_0H, I_1H, I_0H, I_1, I_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      integer iu_veg
      integer i, j, k
      logical veg_data_missing

      real*8, allocatable :: veg_c4(:,:)
      real*8 :: vc4
      integer :: read_c4_grass = 0
      integer variable_lk
      integer jyear

C****
C**** Extract useful local domain parameters from "grid"
C****
      call modelEclock%get(year=jyear)

      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO


      call get_vdata(vdata, (/""/))
c**** read rundeck parameters
      call sync_param( "read_c4_grass", read_c4_grass)
      call  get_param( "variable_lk", variable_lk )

C**** Update vegetation file if necessary (i.e. crops_yr =0 or >0)
      if(crops_yr.eq.0) call updveg(jyear,.false.)
      if(crops_yr.gt.0) call updveg(crops_yr,.false.)

      if (istart.le.2 .or. redogh) then ! initial. foliage arrays (adf)
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0
      end if
      if (istart.le.0) return   ! avoid reading unneeded files

     !!! spgsn=.1d0

c**** check whether ground hydrology data exist at this point.
      veg_data_missing = .false.
      do j=J_0,J_1
        do i=I_0,I_1
          if (variable_lk==0) then
            if ( fearth(i,j) <= 0.d0 ) cycle
          else
            if ( focean(i,j) >= 1.d0 ) cycle
          endif
          if ( sum(vdata(i,j,1:12)).eq.0 ) then
            print *,"No vegetation data: i,j=",i,j,vdata(i,j,1:12)
            veg_data_missing = .true.
          end if
        enddo
      enddo
      if ( veg_data_missing ) then
         if (AM_I_ROOT()) then
            write(6,*) 'Vegetation data is missing at some pts'
            write(6,*) 'If you have a non-standard land mask, please'
            write(6,*) 'consider using extended GH data and rfs file.'
            call stop_model(
     &           'Vegetation data is missing at some cells',255)
         end if
      endif

      end subroutine read_veg_data


      subroutine upd_gh
!@sum initializes (or re-initializes) the vegetation data
      use constant, only : twopi,one
      use TimeConstants_mod, only: DAYS_PER_YEAR
      use Dictionary_mod
      use DOMAIN_DECOMP_ATM, only : GRID
      use DOMAIN_DECOMP_ATM, only : getDomainBounds
      use fluxes, only : focean
      use geom, only : lat2d
      use veg_com !, only : vdata,Cint,Qfol
      use ghy_com, only : ngm,fearth
      use ghy_com, only : dz_ij
!      use surf_albedo, only: albvnh  !nyk !!! not used? i.a.

      implicit none

!---  local vars
      integer i, j
      real*8 dif,frdn,frup,phase,scs0,scsim,scsre,sfv,sla0
      real*8 almass0, almassre, almassim  !nyk
!      real*8 sla0f, slimf, slref !nyk
      real*8 slim,slre,svh,z
      real*8 snm,snf  ! temporary sums (adf)
      real*8 cwc_sum
      integer iv, l
      real*8 fv,dz(ngm)
      integer n
!---  parameters for different types of vegetation
!---          tundr  grass  shrub  trees  decid evrgr  rainf crops
      real*8, parameter :: alamax(11) =
     $     (/ 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0
     &     ,0.d0, 0.d0, 2.d0 /)
      real*8, parameter :: alamin(11) =
     $     (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)

      real*8, parameter :: aroot(11) =
     $     (/ 12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: broot(11) =
     $     (/  1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: rsar(11) =
     $     (/100d0, 100d0, 200d0, 200d0, 200d0, 300d0,250d0, 125d0
     &     ,0.d0, 0.d0, 100d0 /)
      real*8, parameter :: vhght(11) =
     $     (/0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0
     &     ,0.d0, 0.d0, 1.5d0 /)
! Mean canopy nitrogen (nmv; g/m2[leaf]) and Rubisco factors (nfv) for each
! vegetation type (adf)
      real*8, parameter :: nmv(11) =
     $     (/1.6d0,0.82d0,2.38d0,1.03d0,1.25d0,2.9d0,2.7d0,2.50d0
     &     ,0.d0, 0.d0, 0.82d0 /)
      real*8, parameter :: nfv(11) =
     $     (/1.4d0,1.5d0 ,1.3d0 ,1.3d0 ,1.5d0 ,0.9d0,1.1d0,1.3d0
     &     ,0.d0, 0.d0, 1.5d0 /)
      integer, parameter :: laday(11) =
     $     (/ 196,  196,  196,  196,  196,  196,  196,  196
     &     ,0, 0, 196 /)
      real*8, parameter :: can_w_coef(11) =
     &     (/ 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4
     &     ,0.d0, 0.d0, 1.d-4 /)

! Specific leaf areas (sleafa, kg[C]/m2) (adf, nyk)
! Values below 1/(m2/kg) to get kg/m2 for multiplying.
! Sources: White, M.A., et.al. (2000), Earth Interactions, 4:1-85.
!          Leonardos,E.D.,et.al.(2003), Physiologia Plantarum, 117:521+.
!               From winter wheat grown at 20 C (20 m2/kg[dry mass])
!                                      and  5 C (13 m2/kg[dry mass])
!          Francesco Tubiello, personal communication, crop 18-20 m2/kg.
      real*8, parameter :: sleafa(11) =
     $     (/ 1./30.5d0,1./49.0d0,1./30.5d0,1./40.5d0,1./32.0d0,
     $        1./8.2d0,1./32.0d0,1./18.0d0, 0.d0, 0.d0, 1./49.0d0 /)

c****             tundr grass shrub trees decid evrgr rainf crops
c****
c**** laday(veg type, lat belt) = day of peak lai
c old peak lai:  2nd line is for latitudes < 23.5 deg
c****    1  temperate latitudes
c****    2  non-temperate latitudes
c     data  laday/ 196,  196,  196,  196,  196,  196,  105,  196/
c     data  laday/ 196,  288,  288,  288,  288,  196,  105,  288/
c****
c**** contents of ala(k,i,j),  lai coefficients
c****   1  average leaf area index
c****   2  real amplitude of leaf area index
c****   3  imaginary amplitude of leaf area index
c****
c**** contents of almass(i,j),  leaf mass
c****      leaf mass = ala*sleafa = kg[C]/ground area
c****
c**** contents of acs(k,i,j),  cs coefficients
c****   1  average stomatal conductance
c****   2  real amplitude of stomatal conductance
c****   3  imaginary amplitude of stomatal conductance
c****
!@dbparam ghy_default_data if == 1 reset all GHY data to defaults
!@+ (do not read it from files)
      integer :: ghy_default_data = 0
      integer variable_lk

      INTEGER :: I_0, I_1, J_1, J_0, J_1H, J_0H
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

!      entry upd_gh ! need to redo if vdata changes
C****
C**** Extract parameters from "grid" in case we entered here
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      call  get_param( "variable_lk", variable_lk )
c****
c**** set the global arrays  ala, acs, afb, afr, anm, anf
c****
      ala(:,:,:)=0.
      alaf(:,:,:,:)=0.          !nyk lai by veg functional type
      alaif(:,:,:)=0.
      acs(:,:,:)=0.
      avh(:,:)=0.
      afb(:,:)=0.
      afr(:,:,:)=0.
      acs(1,:,:)=.01d0
      anm(:,:)=0.0D0            ! Global mean canopy N array (adf)
      anf(:,:)=0.0D0            ! Global Rubisco factors (adf)
      almass(:,:,:)=0.0D0       ! Global leaf mass at a time, nyk
!rar  aalbveg(:,:)=0.08D0 ! no need, it is set in daily_earth
      can_w_capacity(:,:) = 0.d0

      do j=J_0,J_1
        do i=I_0,I_1
          afb(i,j)=vdata(i,j,1)+vdata(i,j,10)
          if(afb(i,j).gt..999) afb(i,j)=1.

          if (variable_lk==0) then
            if ( fearth(i,j) <= 0.d0 .or. afb(i,j) >= 1.d0 ) cycle
          else
            if ( focean(i,j) >= 1.d0 .or. afb(i,j) >= 1.d0 ) cycle
          endif
          !if(focean(i,j) >= 1.d0 .or. afb(i,j) >= 1.d0) cycle
c**** calculate lai, cs coefficicents
          sfv=0.d0
          sla0=0.     !For calculating sla
          slre=0.     !For calculating sla
          slim=0.     !For calculating sla
          !sla0f=0.     !For calculating alaf
          !slref=0.     !For calculating alaf
          !slimf=0.     !For calculating alaf
          almass0=0.  !nyk For specific leaf area
          almassre=0. !nyk For specific leaf area
          almassim=0. !nyk For specific leaf area
          scs0=0.
          scsre=0.
          scsim=0.
          svh=0.
          snm=0. ! adf
          snf=0. ! adf
          cwc_sum = 0.d0
          do iv=1,11
            if ( iv==9 .or. iv==10 ) cycle
            phase=twopi*laday(iv)/DAYS_PER_YEAR
            if(lat2d(i,j).lt.0.) phase=phase+twopi/2.
            fv=vdata(i,j,iv+1)
            sfv=sfv+fv
            svh=svh+fv*vhght(iv)
            snm=snm+fv*nmv(iv) ! adf
            snf=snf+fv*nfv(iv) ! adf
            dif=(alamax(iv) - alamin(iv))
            sla0=sla0+fv*(alamax(iv) + alamin(iv))
            slre=slre+fv*dif*cos(phase)
            slim=slim+fv*dif*sin(phase)
            !nyk-------------
            !alaf and alaif
            alaf(1,iv,i,j) = 0.5*(alamax(iv) + alamin(iv))
            alaf(2,iv,i,j) = 0.5*dif*cos(phase)
            alaf(3,iv,i,j) = 0.5*dif*sin(phase)
            alaif(iv,i,j)= alaf(1,iv,i,j)+
     $           cosday*alaf(2,iv,i,j)+sinday*alaf(3,iv,i,j) !save ij
           !nyk-------------
            !almaxmin = almaxmin + sleafa(iv)*fv*(alamax(iv)-alamin(iv))
            almass0 = almass0 + sleafa(iv)*fv*(alamax(iv) - alamin(iv))
            !almass0 = almass0 + sleafa(iv)*fv*(alamax(iv) + alamin(iv))
            !almassre = almassre + sleafa(iv)*fv*dif*cos(phase)
            !almassim = almassim + sleafa(iv)*fv*dif*sin(phase)
            !----------------
            scs0=scs0+fv*(alamax(iv) + alamin(iv))/rsar(iv)
            scsre=scsre+fv*dif*cos(phase)/rsar(iv)
            scsim=scsim+fv*dif*sin(phase)/rsar(iv)
            cwc_sum = cwc_sum + fv*can_w_coef(iv)
          end do
          ala(1,i,j)=.5/sfv*sla0
          ala(2,i,j)=.5/sfv*slre
          ala(3,i,j)=.5/sfv*slim

          acs(1,i,j)=.5/sfv*scs0
          acs(2,i,j)=.5/sfv*scsre
          acs(3,i,j)=.5/sfv*scsim
          avh(i,j)=svh/sfv
          can_w_capacity(i,j) = cwc_sum/sfv
          anm(i,j)=snm/sfv ! adf
          anf(i,j)=snf/sfv ! adf
          almass(1,i,j) = almass0   !This just computes total growth for
          almass(2,i,j) = 0.   ! year via difference between max and min
          almass(3,i,j) = 0.
          !almass(1,i,j)= 0.5/sfv*almass0 !nyk
          !almass(2,i,j)= 0.5/sfv*almassre !nyk
          !almass(3,i,j)= 0.5/sfv*almassim !nyk

c**** calculate root fraction afr averaged over vegetation types
          do n=1,ngm
            dz(n)=dz_ij(i,j,n)
            if(dz(n).le.0.) go to 320
          end do
 320      n=n-1
          do iv=1,11
            if ( iv==9 .or. iv==10 ) cycle
            fv=vdata(i,j,iv+1)
            z=0.
            frup=0.
            do l=1,n
              z=z+dz(l)
              frdn=aroot(iv)*z**broot(iv)
              frdn=min(frdn,one)
              if(l.eq.n)frdn=1.
              afr(l,i,j) = afr(l,i,j) + fv*(frdn-frup)
              frup=frdn
            end do
          end do
          do l=1,n
            afr(l,i,j) = afr(l,i,j)/(1.-afb(i,j))
          end do
        end do
      end do

      ! write(98,*) afr
      return
      end subroutine upd_gh


      subroutine reset_veg_to_defaults( reset_prognostic )
      use veg_com, only: vdata
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      logical, intent(in) :: reset_prognostic
      integer i,j

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
      do i=I_0,I_1

      vdata(i,j,1:12)= (/  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.62451148d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.37548852d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00 /)
      enddo
      enddo

      end subroutine reset_veg_to_defaults


      subroutine veg_set_cell (vegcell,i0,j0,pres,ts,ghy_data_only)
      !subroutine veg_set_cell (i0,j0,ghy_data_only)
!@sum resets the vegetation module to a new cell i0,j0
      use constant, only : gasc
      use ghy_com, only : ngm
      use sle001, only : shw !,pres,ts !fr,snowm, ws,shc,
      use vegetation, only :  ! alaie,rs,nm,nf,alai,vh
     &     ! alait,vfraction 
     &     !!!Qf,!Ci, ! added by adf, nyk !fdir,parinc,vegalbedo,sbeta,!Ci,
     &     !Ca,     ! nyk
     &     cond_scheme,vegCO2X_off  !nyk
     &     !,cnc
      use rad_com, only : cosz1, CO2X, CO2ppm  !nyk
     &    ,FSRDIR,SRVISSURF  !adf, nyk
      use veg_com, only :
     &     Cint,Qfol           ! added by adf
     &     ,cnc_ij
     &     ,aalbveg    ! nyk
     &     ,afr,avh,anm,anf,ala,alaif,alaf,vdata,acs
      use vegetation, only : t_vegcell

      implicit none

      type(t_vegcell) :: vegcell
      integer, intent(in) :: i0,j0
      real*8, intent(in)  :: pres,ts
!@var ghy_data_only set only ghy data (no call to vegetatin expected)
      logical, optional :: ghy_data_only
      real*8, parameter :: spgsn=.1d0
      integer l
      real*8 aa
!      real*8 aalbveg0, sfv  !nyk
!      integer northsouth,iv  !nyk
!      real*8 alaic,vh,shtpr,alai  ! adf
      real*8 alaic, alai

      ! call stop_model("veg_set_cell called", 255)

c**** fr: root fraction in layer l  (1=fr(1)+fr(2)+...+fr(n))
      do l=1,ngm
        vegcell%fr(l)=afr(l,i0,j0)
      end do
c**** vh: vegetation height
      vegcell%vh=avh(i0,j0)
      vegcell%nm=anm(i0,j0)     ! mean canopy nitrogen (g/m2[leaf]) (adf)
      vegcell%nf=anf(i0,j0)     ! canopy nitrogen factor (adf)
      vegcell%snowm=vegcell%vh*spgsn

c**** alai: leaf area index
      alai=ala(1,i0,j0)+cosday*ala(2,i0,j0)+sinday*ala(3,i0,j0)
      alai=max(alai,1.d0)
      !lai by functional type, not normalized by cover fraction!
      alaif(:,i0,j0) = 0.d0
      vegcell%alait(:) = 0.d0
      vegcell%vfraction(:) = 0.d0
      do L=1,11
        if ( L==9 .or. L==10 ) cycle
        alaif(L,i0,j0)=alaf(1,L,i0,j0)+
     &       cosday*alaf(2,L,i0,j0)+sinday*alaf(3,L,i0,j0) !save ij
        vegcell%alait(L) = alaif(L,i0,j0)  !for VEGETATION.f
        vegcell%vfraction(L)=vdata(i0,j0,L+1) !for VEGETATION.f
      end do

      alaic=5.0
      vegcell%alaie=alaic*(1.-exp(-alai/alaic))
c**** rs: minimum stomatal resistance
      vegcell%rs=alai/
     &     (acs(1,i0,j0)+cosday*acs(2,i0,j0)+sinday*acs(3,i0,j0))
c???  cnc=alai/rs   redefined before being used (qsbal,cond)

      !---------------------------------------------------------
      !nyk vegetation albedo.  Only really updated daily, but have to
      !get it initialized somewhere after ALBVNH is calculated.
      !ALBVNH is unfortunately set *after* ground hydr is initialized.
      !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
!      aalbveg0=0.d0               !nyk
!      sfv=0.d0
!      if (j0.le.jm/2) then
!
!       northsouth=1            !southern hemisphere
!      else
!        northsouth=2            !northern hemisphere
!      end if
!      do iv=1,8
!        aalbveg0 = aalbveg0 + fv*(ALBVNH(iv+1,1,northsouth))
!        sfv = sfv + fv
!      end do
!      aalbveg(i0,j0) = aalbveg0/sfv
!      write (99,*) 'aalbveg ghinij', aalbveg(i0,j0) !nyk
      !---------------------------------------------------------

      !ws(0,2)=.0001d0*alai
      vegcell%ws_can = .0001d0*alai
!!!      ws(0,2)=can_w_capacity(i0,j0)*alai
c shc(0,2) is the heat capacity of the canopy
      aa=ala(1,i0,j0)
      !shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw
      vegcell%shc_can=(.010d0+.002d0*aa+.001d0*aa**2)*shw

! This is a hack to prevent the program from using uninitialized
! radiation arrays at the beginning of the run.
! If program returns at this point, it sets only the data needed for
! ground hydrology but not the vegetation which is ok if veg_conductance
! is not called.
      if ( present(ghy_data_only) ) then
        if ( ghy_data_only ) return
      endif

!----------------------------------------------------------------------!
      if (cond_scheme.eq.2) then  !new conductance scheme
! Sine of solar elevation (rad).
        vegcell%sbeta=cosz1(i0,j0)
! adf Fraction of solar radiation at ground that is direct beam.
        vegcell%fdir=FSRDIR(i0,j0)
! nyk Calculate incident PAR, photosynthetically active radiation.
!     *SRVISSURF is from SRDVIS:
!     = incident visible solar radiation (dir+dif) on the surface
!     = estimated as 53% of total solar flux density, so wavelength
!       range is UV through ~760 or 770 nm cutoff (not strict)
!     *SRDVIS is normalized to solar zenith = 0.
!     *From integrating solar flux density (TOA) over solar spectrum,
!       PAR(400-700 nm) is ~ 82% of the flux density of SRDVIS.
        vegcell%parinc=0.82*SRVISSURF(i0,j0)*vegcell%sbeta
! nyk Get vegetation grid albedo, temporary until canopy scheme in place
        vegcell%vegalbedo = aalbveg(i0,j0)
! Internal foliage CO2 concentration (mol/m3).
        vegcell%Ci=Cint(i0,j0)
! Foliage surface mixing ratio (kg/kg).
        vegcell%Qf=Qfol(i0,j0)
        !!!Qf=Qfol(i0,j0)
        vegcell%CNC = cnc_ij(i0,j0)
! Atmospheric CO2 concentration at land surface (mol/m3)
        if (vegCO2X_off==0) then  !Default
           vegcell%Ca = CO2X*CO2ppm*(1.0D-06)*pres*100.0/gasc/ts
        else
           vegcell%Ca = CO2ppm*(1.0D-06)*pres*100.0/gasc/ts
        endif
!        write(99,*) "vegCO2X_off,CO2X,CO2ppm,Ca,pres,gasc,ts",
!     %       vegCO2X_off,CO2X,CO2ppm,Ca,pres,gasc,ts
!        stop
      end if
!----------------------------------------------------------------------!
      vegcell%alai = alai
      return
      end subroutine veg_set_cell


      subroutine veg_save_cell(i,j)
      !Save vegetation arrays for cell.
      use vegetation, only: cond_scheme, Ci, Qf, CNC
      use veg_com, only:  Cint, Qfol, cnc_ij

      integer, intent(in):: i,j

      if (cond_scheme.eq.2) then      !new conductance scheme only
        Cint(i,j)=Ci                  ! Save last value
        Qfol(i,j)=Qf                  ! Save last value
        cnc_ij(i,j)=CNC
      end if
      end subroutine veg_save_cell



      ! end module veg_drv

!***********************************************************************

      subroutine updveg (year,reset_veg)
!@sum  reads appropriate crops data and updates the vegetation file
!@auth R. Ruedy
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, AM_I_ROOT
      use veg_com, only : vdata
      USE GEOM, only : imaxj
      !use veg_drv, only : upd_gh
      implicit none
      integer, intent(in) :: year
      logical, intent(in) :: reset_veg
      !---
      real*8 wt,crops
      integer, save :: year_old=-1
      integer :: i,j,k
      real*8, allocatable :: cropdata(:,:)
      INTEGER :: I_0, I_1, J_1, J_0, J_0H, J_1H, I_0H, I_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO =J_0H,   J_STOP_HALO =J_1H )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** check whether update is needed
      if (year.eq.year_old) return

      call get_vdata(vdata, (/""/))

c****     check whether a no-crops vege-file was used
      if( maxval( vdata(I_0:I_1, J_0:J_1,9) ) > 0.d0)
     *     call stop_model('updveg: use no_crops_VEG_file',255)

      allocate( cropdata(I_0H:I_1H, J_0H:J_1H) )
      call get_cropdata(year, cropdata)

      if (AM_I_ROOT()) write(6,*)
     *     'Using crops data from year',year
c**** Modify the vegetation fractions
      do j=J_0,J_1
      do i=I_0,I_1
        crops = cropdata(i,j)
         if ( crops > 0.d0 ) then
            do k=1,12
              vdata(i,j,k) = vdata(i,j,k)*(1.d0-crops)
            end do
              vdata(i,j,9) = crops
         end if
      end do
      end do
      deallocate( cropdata )

      if(reset_veg) call upd_gh

      year_old = year
      return
      end subroutine updveg

#endif


      subroutine get_vdata(vdata, vegnames)
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds!, AM_I_ROOT
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      !use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      !use veg_com
      !use model_com, only : jyear,focean
      !use ghy_com, only : fearth
#ifdef USE_ENT
      use ent_mod, only: N_COVERTYPES, N_OTHER,COVER_SAND !ykim - use ent_const to accomodate GISS and Ent PFTs.
#endif
      implicit none
#ifndef USE_ENT
      integer, parameter :: N_COVERTYPES = 12
      integer, parameter :: N_OTHER = 2
      integer, parameter :: COVER_SAND = 1
#endif
      real*8, intent(out) :: vdata(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      character(len=*) :: vegnames(:)
      !---
      INTEGER :: J_1, J_0, J_1H, J_0H, I_1H, I_0H, I_1, I_0
      integer :: i, j, k, fid
      real*8 :: s
      !character(len=32) :: vegnames(N_COVERTYPES-N_OTHER)
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** read land surface parameters or use defaults
      ! vegnames should be set elsewhere.  Do this once the model
      ! actually uses Ent vegetation types as input.
cddd      vegnames = (/
cddd     &     'brightsoil     ','tundra         ','grass          ',
cddd     &     'shrub_and_grass','tree_and_grass ','deciduous      ',
cddd     &     'evergreen      ','rainforest     ','cultivation    ',
cddd     &     'darksoil       '
cddd     &     /)

      if(file_exists('VEG')) then
        fid = par_open(grid,'VEG','read')
        do k=1,size(vegnames)
          vdata(:,:,k) = 0.     ! if type is absent in file, assume 0
          call read_dist_data(grid,fid,trim(vegnames(k)),vdata(:,:,k))
        enddo
        call par_close(grid,fid)
      else
        vdata(:,:,:) = 0.
        vdata(:,:,1) = COVER_SAND ! all bare soil if no input data available
      endif

c**** zero-out vdata(11) until it is properly read in
      do k=N_COVERTYPES-N_OTHER+1, N_COVERTYPES
        vdata(:,:,k) = 0.
      end do

      ! make sure that veg fractions are reasonable
      do j=J_0,J_1
        do i=I_0,I_1
          !print *, i,j
          !print *, vdata(i,j,:) 
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vdata(i,j,k) < 1.d-4 ) vdata(i,j,k) = 0.d0
          enddo
          s = sum( vdata(i,j,:) )
          if ( s > .9d0 ) then
            vdata(i,j,:) = vdata(i,j,:)/s
          else if ( s < .1d0 ) then
            print *, "missing veg data at ",i,j,"assume bare soil"
            vdata(i,j,: ) = 0.d0
            vdata(i,j,COVER_SAND) = 1.d0
          else
            print *,i,j,s
            print *, vdata(i,j,:) 
            call stop_model("Incorrect data in VEG file",255)
          endif
        enddo
      enddo

      end subroutine get_vdata

cddd      module cropdata_mod
cddd      use timestream_mod, only : timestream
cddd      implicit none
cddd!@var CROPstream interface for reading and time-interpolating the crop file
cddd!@+   See usage notes in timestream_mod
cddd      type(timestream) :: CROPstream
cddd      logical :: have_crops_file
cddd      end module cropdata_mod
      subroutine get_cropdata(year, cropdata)
!@sum get_cropdata reads timeseries file for crop fraction and
!@+   interpolates to requested year.

      use domain_decomp_atm, only : grid
      use timestream_mod, only : init_stream,read_stream,timestream
      !use cropdata_mod
      use filemanager, only : file_exists
      implicit none
      integer, intent(in) :: year
      real*8 :: cropdata(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO)
c
      logical, save :: init = .false.
      integer :: day
      type(timestream), save :: CROPstream
      logical, save :: have_crops_file

      day = 1 ! to pass a required argument

      if (.not. init) then
        init = .true.
        have_crops_file = file_exists('CROPS')
        if(have_crops_file) then
          call init_stream(grid,CROPstream,'CROPS','crops',
     &         0d0,1d30,'none',year,day)
        endif
      endif

      cropdata = 0.d0 ! make sure halo is intialized
      if(have_crops_file) then
        call read_stream(grid,CROPstream,year,day,cropdata)
        !Set min to zero, since no land mask yet -nyk 1/22/08
        cropdata = max(0d0, cropdata)
      endif

      end subroutine get_cropdata



      subroutine get_laimaxdata(laimaxdata)
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds!, AM_I_ROOT
      use DOMAIN_DECOMP_ATM, only : READT_PARALLEL
      use filemanager
      !use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      !use veg_com
      !use model_com, only : jyear,focean
      !use ghy_com, only : fearth
#ifdef USE_ENT
      use ent_mod, only: N_COVERTYPES, N_OTHER,COVER_SAND !ykim - use ent_const to accomodate GISS and Ent PFTs.
#endif
      implicit none
#ifndef USE_ENT
      integer, parameter :: N_COVERTYPES = 12
      integer, parameter :: N_OTHER = 2
      integer, parameter :: COVER_SAND = 1
#endif
      real*8, intent(out) :: laimaxdata(
     &     grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      !---
      INTEGER :: J_1, J_0, J_1H, J_0H, I_1H, I_0H, I_1, I_0
      integer :: i, j, k, iu_laimax
      real*8 :: s

      CALL getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** read land surface parameters or use defaults
      call openunit("LAImax",iu_LAIMAX,.true.,.true.)
      do k=1,N_COVERTYPES-N_OTHER
        CALL READT_PARALLEL
     *    (grid,iu_LAIMAX,NAMEUNIT(iu_LAIMAX),laimaxdata(:,:,K),1)
      end do
c**** zero-out laimaxdata(11) until it is properly read in
      do k=N_COVERTYPES-N_OTHER+1, N_COVERTYPES
        laimaxdata(:,:,k) = 0.
      end do
      call closeunit(iu_LAIMAX)

      end subroutine get_laimaxdata



      subroutine get_soil_C_total(ncasa, soil_C_total)
      use DOMAIN_DECOMP_ATM, only : GRID
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      implicit none
      integer, intent(in) :: ncasa
      real*8 :: !,intent(out) ::
     &     soil_C_total(ncasa,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      !---
      integer :: fid
      if(file_exists('SOILCARB_global')) then
        fid = par_open(grid,'SOILCARB_global','read')
        call read_dist_data(grid,fid,'soil_C_total',soil_C_total,jdim=3)
        call par_close(grid,fid)
      else
        soil_C_total = 0.
      endif
      end subroutine get_soil_C_total

      end module veg_drv
