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
C**** this is mass*age (kg*year)
C**** age=1/(JDperY*24*3600) in years
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
