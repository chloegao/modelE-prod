#include "rundeck_opts.h"

! A WIP version of OCN_TRACER appropriate for RUNTIME_NTM_OCEAN and
! tracers whose ICs are defined in a netcdf file containing arrays
! whose names match those specified in ocean_trname in the rundeck.

      subroutine tracer_ic_ocean(atmocn)
      use model_com, only: itime,itimei
      use ocn_tracer_com, only : tracerlist, ocn_tracer_entry
      use ocean, only : im,jm,lmo
      use ocean, only : dxypo,mo
      use ocean, only : trmo
      use ocean, only : txmo,tymo,tzmo
      use ocean, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : getDomainBounds, am_i_root
      use oceanr_dim, only : grid=>ogrid
      use exchange_types, only : atmocn_xchng_vars
      use dictionary_mod
      use pario, only : read_dist_data,par_open,par_close
      implicit none
      type(atmocn_xchng_vars) :: atmocn
c
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     tr_ic
      integer n,i,j,l,nt,fid
      integer :: j_0s, j_1s, j_0, j_1
      type(ocn_tracer_entry), pointer :: entry

      if(itime.ne.itimei) return

      call getDomainBounds(grid, j_strt_skp = j_0s, j_stop_skp = j_1s,
     &     j_strt = j_0, j_stop = j_1)

      fid = par_open(grid,'OCN_TRACER_IC','read')

! Loop over tracers, read the IC for each, convert to extensive units (kg).
! straits IC not an option yet.
      do nt=1,tracerlist%getsize()
        entry=>tracerlist%at(nt)
        call read_dist_data(grid,fid,trim(entry%trname),tr_ic)
        do l=1,lmo
        do j=j_0,j_1
        trmo(:,j,l,nt) = 0.
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          trmo(i,j,l,nt) = tr_ic(i,j,l)*mo(i,j,l)*dxypo(j)
        enddo
        enddo
        enddo
        enddo
      enddo

      call par_close(grid,fid)

      return
      end subroutine tracer_ic_ocean

      subroutine oc_tdecay(dts)
      implicit none
      real*8, intent(in) :: dts
      return
      end subroutine oc_tdecay

      SUBROUTINE OCN_TR_AGE(DTS)
!@sum OCN_TR_AGE age tracers in ocean
!@auth Gavin Schmidt/Natassa Romanou
      USE MODEL_COM, only : itime
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE OCN_TRACER_COM, only : n_age
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 age_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=0 and add 1 below
C**** this is mass*age (kg*year)
C**** age=1/(JDperY*24*3600) in years
      age_inc=dts/(INT_DAYS_PER_YEAR*SECONDS_PER_DAY)
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,1,n_age)=0 ; TXMO(I,J,1,n_age)=0 
                TYMO(I,J,1,n_age)=0 ; TZMO(I,J,1,n_age)=0
              else
                TRMO(I,J,L,n_age)= TRMO(I,J,L,n_age) +
     +                age_inc * MO(I,J,L) * oXYP(I,J)
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_age

      SUBROUTINE OCN_TR_CFC(DTS)
!@sum OCN_TR_CFC tracer in ocean
!@auth Natassa Romanou
      USE Dictionary_mod, only : get_param
      USE MODEL_COM, only : itime,modelEclock
      USE CONSTANT,   only : grav
      USE OCN_TRACER_COM, only : n_cfc,icfcyear,cfc11nh,cfc11sh
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo,dxypo,g0m,s0m,olat=>olat2d_dg ! 2D array containing lat at each i,j
      USE OFLUXES,    only : oRSI,oAPRESS,ocnatm
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid
      USE ODIAG, only : ij_cfcair,ij_kw,ij_csat,ij_cfcflux,oij=>oij_loc
      use model_com, only: modeleclock
      use runtimecontrols_mod, only: ocn_cfc

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 :: cfc_inc, Xconv,a,pres,g,s,sst,sss,temgsp,wind,pnoice,Xkw
     .              ,solub,solub_cfc,schmidtno_cfc,Sc,kw,cfcair,csat
     .              ,fluxa,flux,flux_tendency,rho_water,dp1d
     .              ,Pnorth,Psouth,trmopro,fluxb
      real*8 :: ys ! northern boundary of SH constant-value domain (deg N)
      real*8 :: yn ! southern boundary of NH constant-value domain (deg N)
      real*8 :: wt_sh ! weight for SH constant-value domain
      real*8,External   :: VOLGSP
      integer i,j,l,k
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1,year, month, dayOfYear, date
      real*8 :: cfc_conc_const

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** interpolate atmospheric values
      call modelEclock%get(year=year, month=month, date=date,
     &     dayOfYear=dayOfYear)

      !pick appropriate values for each year
      do i=1,105
        if (year == icfcyear(i)) then
         yn=10.d0
         ys=-10.d0
         Pnorth=cfc11nh(i)
         Psouth=cfc11sh(i)
        endif
      enddo


C**** at each time step set surface tracer conc=1+flux from atmos
      Xconv = 1/3.6d5
      a = 0.337
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
      IF (focean(i,j)>0) then
      wind=ocnatm%wsavg(i,j)        !owind(i,j)
      pnoice = 1.d0 - oRSI(i,j)     !1-fice
      k = 1    !surface only
      pres = oAPRESS(i,j)    !surface atm. pressure
     .     + MO(I,J,k)*GRAV*.5   !pressure at first layer
      pres = pres * 0.00000986923266716     ! atm
        g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        sst=TEMGSP(g,s,pres)     !in situ   temperature
        sss=s*1000.d0            !convert to psu (eg. ocean mean salinity=35psu)
        rho_water = 1d0/VOLGSP(g,s,pres)
        dp1d = MO(I,J,K)/rho_water   !local thickenss of each layer in meters
!     if (i.eq.50.and.j.eq.90) then
!     write(*,'(a,5i5,5e12.4)')'CFC OUTPUT:',
!    .   i,j,date,month,year,
!    .   pres,sst,sss,rho_water,dp1d
!     endif

      Xkw = Xconv * a * wind**2       ! units in m/s
      solub = solub_cfc(sst,sss,11)   !mol/m3/pptv
      Sc = schmidtno_cfc(sst,11)
      kw = pnoice*Xkw/sqrt(Sc/660)

!     if (i.eq.50.and.j.eq.90) then
!     write(*,'(a,5i5,6e12.4)')'CFC OUTPUT:',
!    .   i,j,date,month,year,
!    .   Xconv,a,wind,solub,Sc,kw   
!     endif
#ifdef OCN_CFCconst
!     cfcair = 1.          !pptv to derive green's functions -- corresponds to year=1951
      call get_param('cfc_conc_const',cfc_conc_const)
      cfcair=cfc_conc_const
#else
      if(olat(i,j) <= ys) then
          wt_sh = 1d0
      elseif(olat(i,j) >= yn) then
          wt_sh = 0d0
      else
          wt_sh = (yn-olat(i,j))/(yn-ys)
      endif
      cfcair = (1.-wt_sh) * Pnorth + wt_sh * Psouth
#endif

!mo units: kg/m2
      csat = solub * cfcair * pres/1       ! mol/m3
      fluxa = kw * csat                     ! mol/m2/s
      fluxb = kw * trmo(i,j,1,n_cfc) *1000.d0/137.37d0
     .           * rho_water/mo(i,j,1)/dxypo(j)            !mol/m2/s
      flux = fluxa - fluxb     !mol/m2/s

      flux_tendency= flux /dp1d * 1000.d0*pnoice  !mili-mol/m3/s

!     if (i.eq.50.and.j.eq.90) then
!      write(*,'(a,5i5,9e12.4)')'CFC OUTPUT 1:', 
!    . i,j,date,month,year, cfcair,pres,solub,csat,kw,fluxa,
!    . flux,flux_tendency,trmo(i,j,1,n_cfc)
!     endif

      trmopro=trmo(i,j,1,n_cfc)   !save for printout
      trmo(i,j,1,n_cfc) =  trmo(i,j,1,n_cfc)
     .                  + (flux_tendency * DTS)*1e-6*137.37
     .                                 *mo(i,j,1)*dxypo(j)/rho_water   !kg
      if (i.eq.50.and.j.eq.90) then
         write(*,'(a,5i5,12e12.4)')'CFC OUTPUT 1:',
     . date,month,year,i,j,pnoice,cfcair,pres,solub,csat,kw,fluxa,
     . fluxb,flux,dp1d,trmopro,trmo(i,j,1,n_cfc)
      endif


!     we do not set the moments for the tracer field here... 
!     maybe we need to do that for the online gasexchange though...
!     TXMO(I,J,1,n_cfc)=0
!     TYMO(I,J,1,n_cfc)=0 ; TZMO(I,J,1,n_cfc)=0

       if (ocn_cfc) then
         OIJ(I,J,IJ_cfcair) = OIJ(I,J,IJ_cfcair) + cfcair
         OIJ(I,J,IJ_kw) = OIJ(I,J,IJ_kw) +  kw
         OIJ(I,J,IJ_csat) = OIJ(I,J,IJ_csat) +  csat
         OIJ(I,J,IJ_cfcflux) = OIJ(I,J,IJ_cfcflux) +  flux
      endif
      ENDIF
      ENDDO
      ENDDO

      end SUBROUTINE OCN_TR_CFC


      SUBROUTINE read_atmcfc
  
      USE FILEMANAGER, only: openunit,closeunit
      use ocn_tracer_com, only: icfcyear,cfc11nh,cfc11sh

      implicit none
      integer iu_file,i
      character(len=80) :: first_line_dummy

      allocate(icfcyear(105),cfc11nh(105), cfc11sh(105))

      !read in the values from cfc1112.atm
      call openunit('cfcatm_data',iu_file,.false.,.true.)
      read(iu_file,*)
      do i=1,105
      read(iu_file,*)icfcyear(i), cfc11nh(i), cfc11sh(i)
      enddo
      call closeunit(iu_file)

      end SUBROUTINE read_atmcfc


      real*8 function solub_cfc(pt,ps,kn)

!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /www/html/ipsl/OCMIP/phase2/simulations/CFC/boundcond/RCS/sol_cfc.f,v $
!_ $Revision: 1.2 $
!_ $Date: 1998/07/17 07:37:02 $   ;  $State: Exp $
!_ $Author: jomce $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sol_cfc.f,v $
!_ Revision 1.2  1998/07/17 07:37:02  jomce
!_ Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
!_ following warning from Matthew Hecht at NCAR.
!_
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!     function sol_cfc=solub_cfc(pt,ps,kn)
!-------------------------------------------------------------------
!
!     CFC 11 and 12 Solubilities in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!
!     pt:       temperature (degre Celcius)
!     ps:       salinity    (o/oo)
!     kn:       11 = CFC-11, 12 = CFC-12
!     sol_cfc:  in mol/m3/pptv
!               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

!
!     J-C Dutay - LSCE
!-------------------------------------------------------------------
      implicit none
      real*8 :: a1,a2,a3,a4,b1,b2,b3,ta,d,pt,ps
      integer kn

!c coefficient for solubility in  mol/l/atm
!c ----------------------------------------
!
!     for CFC 11
!     ----------
      if (kn.eq.11) then
            a1 = -229.9261d0
            a2 =  319.6552d0
            a3 =  119.4471d0
            a4 =   -1.39165d0
            b1 =   -0.142382d0
            b2 =    0.091459d0
            b3 =   -0.0157274d0
      endif
!    
!     for CFC/12
!     ----------
      if (kn.eq.12) then
             a1 = -218.0971d0
             a2 =  298.9702d0
             a3 =  113.8049d0
             a4 =   -1.39165d0
             b1 =   -0.143566d0
             b2 =    0.091015d0
             b3 =   -0.0153924d0
       endif
 

      ta       = ( pt + 273.16d0)* 0.01d0
      d    = (b3*ta + b2)*ta + b1
 
 
      solub_cfc = exp(a1 + a2/ta + a3*log(ta) + a4*ta*ta + ps*d)

!     conversion from mol/(l * atm) to mol/(m^3 * atm) 
!     ------------------------------------------------
      solub_cfc = 1000.d0 * solub_cfc
 
!     conversion from mol/(m^3 * atm) to mol/(m3 * pptv) 
!     --------------------------------------------------
      solub_cfc = 1.0d-12 * solub_cfc

      end function solub_cfc


      real*8 function schmidtno_cfc(t,kn)

!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /home/orr/ocmip/web/OCMIP/phase2/simulations/CFC/boundcond/RCS/sc_cfc.f,v $
!_ $Revision: 1.1 $
!_ $Date: 1998/07/07 15:22:00 $   ;  $State: Exp $
!_ $Author: orr $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sc_cfc.f,v $
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!---------------------------------------------------
!     CFC 11 and 12 schmidt number 
!     as a fonction of temperature. 
!
!     ref: Zheng et al (1998), JGR, vol 103,No C1 
!
!     t: temperature (degree Celcius)
!     kn: = 11 for CFC-11,  12 for CFC-12
!
!     J-C Dutay - LSCE
!---------------------------------------------------
!
!   coefficients with t in degre Celcius
!   ------------------------------------
      implicit none
      real*8 :: a1,a2,a3,a4,t
      integer kn

      if (kn==11) then
      a1 = 3501.8d0
      a2 = -210.31d0
      a3 =    6.1851d0
      a4 =   -0.07513d0
      endif
 
      if (kn==12) then
      a1 = 3845.4d0
      a2 = -228.95d0
      a3 =    6.1908d0
      a4 =   -0.067430d0
      endif

      schmidtno_cfc = a1 + a2 * t + a3 *t*t + a4 *t*t*t;

      end function schmidtno_cfc

      SUBROUTINE OCN_TR_VENT(DTS)
!@sum OCN_VENT tracer in ocean
!@auth Natassa Romanou
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_vent
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 vent_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=1
      vent_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_vent)= vent_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_vent)= 0. 
                TYMO(I,J,L,n_vent)= 0.
                TZMO(I,J,L,n_vent)= 0.
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_vent

      SUBROUTINE OCN_TR_GASX(DTS)
!@sum OCN_GASX tracer in ocean
!@auth Natassa Romanou
! same as ventilation tracer but take ice into account
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_gasx
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid
      USE OFLUXES, only : oRSI

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 gasx_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

      gasx_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
#ifdef OTRAC_unspecified
                !under ice do not specify anything
                if (oRSI(i,j).ne.1.) then
                TRMO(I,J,L,n_gasx)= gasx_inc * MO(I,J,L) * oXYP(I,J)
     .                            * (1-oRSI(i,j))
                else
                write(*,'(a,3i5,2e12.4)')'ice frac=1',
     .          itime,i,j,orsi(i,j),trmo(i,j,l,n_gasx)
                endif
#else
        !in this formulation we specify 0 under the ice
                TRMO(I,J,L,n_gasx)= gasx_inc * MO(I,J,L) * oXYP(I,J)
     .                            * (1-oRSI(i,j))
                TXMO(I,J,L,n_gasx)= 0.
                TYMO(I,J,L,n_gasx)= 0.
                TZMO(I,J,L,n_gasx)= 0.
#endif
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_gasx

      SUBROUTINE OCN_TR_WaterMass(DTS)
!@sum OCN_WaterMass tracer in ocean
!@auth Natassa Romanou
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_wms1,n_wms2,n_wms3
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo, oLON_DG,oLAT_DG,ZOE=>ZE


      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 wms1_inc,wms2_inc,wms3_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=1
C**** this is mass*age (kg*year)
C**** age=1/(JDperY*24*3600) in years
      if (n_wms1 /= 0) then
      wms1_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          if (oLAT_DG(j,2) .le. -55.d0) then    !Antarctic
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_wms1)= wms1_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms1)= 0.
                TYMO(I,J,L,n_wms1)= 0.
                TZMO(I,J,L,n_wms1)= 0.
              end if
            end if
          ENDDO
          endif
        ENDDO
      ENDDO
      endif

      if (n_wms2 /= 0) then
      wms2_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          if (oLAT_DG(j,2) .ge. 40.d0) then    !North Atlantic
          DO I=1,IMAXJ(J)
          if (oLON_DG(i,2) .ge. -90.d0 .and. oLON_DG(i,2) .le. 60.d0)
     .      then    !North Atlantic
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_wms2)= wms2_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms2)= 0.
                TYMO(I,J,L,n_wms2)= 0.
                TZMO(I,J,L,n_wms2)= 0.
              end if
            end if
          end if
          ENDDO
          endif
        ENDDO
      ENDDO
      endif

      if (n_wms3 /= 0) then
      wms3_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (zoe(l).ge.3000.) then
                TRMO(I,J,L,n_wms3)= wms3_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms3)= 0.
                TYMO(I,J,L,n_wms3)= 0.
                TZMO(I,J,L,n_wms3)= 0.
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
      endif
C****
      return
      end subroutine ocn_tr_WaterMass

       subroutine OCN_TR_DetrSettl(DTS)
!@sum OCN_TR_DetrSettl in ocean
!@auth Natassa Romanou
      USE CONSTANT, only : grav
      USE MODEL_COM, only : itime,itimei
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      USE OCN_TRACER_COM, only : n_dets
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo,g0m,s0m,dxypo
      USE OFLUXES,    only : oAPRESS

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      Real*8,External   :: VOLGSP
      integer i,j,l
      real*8  trnd,wsdeth,wsdet,rho_water,g,s,pres
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

!---------------------------------------------------------------
! ---  detrital settling
! tracer to simulate detrital settling for ideal tracer with initial
! distribution = 1e3
!---------------------------------------------------------------
!
!!the sinking term is given in units (m/hr)*(mgr,chl/m3)
!!in order to be converted into mgr,chl/m3/hr as the tendency
!!terms are in the phytoplankton equations, 
!!we need to multiply by dz of each layer:
!!  dz(k  ) * P_tend(k  ) = dz(k  ) * P_tend(k  ) - trnd
!!  dz(k+1) * P_tend(k+1) = dz(k+1) * P_tend(k+1) + trnd
!!this way we ensure conservation of tracer after vertical adjustment
!!the /hr factor is bcz the obio timestep is in hrs.
!
      wsdeth = 20.0/24.0     !as for nitrogen (m/hr)
      wsdet = wsdeth

! adjust settling velocity based on temperature
!     do k = 1,kmax
!       wsdet(k) = wsdeth*viscfac(k)*pnoice(k)
!     enddo

! convert to m/s
      wsdet = wsdet/SECONDS_PER_HOUR

      ! initialization :-)
      if (itime.eq.itimei) then 
      DO J=J_0,J_1;  DO I=1,IMAXJ(J)
      DO L=1,LMO-1
        if (l.le.lmm(i,j)) then
          trmo(I,J,L,n_dets) = 1000.d0*MO(I,J,L)*oXYP(I,J)
        endif
      ENDDO
      ENDDO;ENDDO
      endif

      DO J=J_0,J_1;  DO I=1,IMAXJ(J)
      pres = oAPRESS(i,j)    !surface atm. pressure
      DO L=1,LMO-1
        if (l.le.lmm(i,j)) then
         pres=pres+MO(I,J,L)*GRAV*.5
         g=G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
         s=S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
         rho_water = 1d0/VOLGSP(g,s,pres)
         trnd = trmo(I,J,L,n_dets)*wsdet*DTS*rho_water/MO(I,J,L)
         trmo(I,J,L  ,n_dets)=trmo(I,J,L  ,n_dets) - trnd
         trmo(I,J,L+1,n_dets)=trmo(I,J,L+1,n_dets) + trnd

!        TXMO(I,J,L,n_dets)= 0.
!        TYMO(I,J,L,n_dets)= 0.
!        TZMO(I,J,L,n_dets)= 0.
        endif
      ENDDO
      ENDDO;ENDDO
      
!!let detritus that reaches the bottom, disappear in the sediment
!!        k = kmax
!!        trnd = det(k,nt)*wsdet(k,nt)
!!        D_tend(k,nt)   = D_tend(k,nt)   - trnd/dp1d(k)

!      !diagnostic for carbon export at compensation depth
!      cexp = 0.
!      do  k=1,kzc    
!      !term2: settling C detritus contribution
!      !dont set cexp = 0 here, because adds to before
!        nt= 1            !only the for carbon detritus
!        cexp = cexp 
!     .        + det(k,nt)*wsdet(k,nt)
!     .        * 24.d0 * 365.d0
!     .        * 1.d-15                 !ugC/l -> PgC/yr
!     .        * dxypo(j) 
!      enddo
       end subroutine OCN_TR_DetrSettl

      SUBROUTINE OCN_TR_DIC(DTS)
!@sum OCN_DIC tracer in ocean
! dic preindustrial in the ocean + gas exchange
!@auth Natassa Romanou
      USE MODEL_COM, only : itime,itime,itimei
      USE OCN_TRACER_COM, only : n_vent
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

!     ! initialization :-)
!     if (itime.eq.itimei) then 
!     filename='dic_inicond'
!     call bio_inicond_g(filename,fldo2,fldoz)
!         trmo(:,:,:,n_dets) = fldo2
!     endif

C****
      return
      end subroutine ocn_tr_dic

      subroutine diagtco (m,nt0,atmocn)
      use exchange_types, only : atmocn_xchng_vars
      implicit none
      integer, intent(in) :: m
      integer, intent(in) :: nt0
      type(atmocn_xchng_vars) :: atmocn
      return
      end subroutine diagtco
