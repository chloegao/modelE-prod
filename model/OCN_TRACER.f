#include "rundeck_opts.h"

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@sum  OCN_TRACER: tracer-dependent routines for GISS Ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Tracer initialisation + sources: tracer_ic_ocean
!@+
!@auth Jean Lerner/Gavin Schmidt

      subroutine tracer_ic_ocean(atmocn)
!@sum tracer_ic_ocean initialise ocean tracers
!@auth Gavin Schmidt
!@ver 1.0
      USE MODEL_COM, only: itime
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry, n_age
#ifdef TRACERS_SPECIAL_O18
      USE OCN_TRACER_COM, only : water_tracer_ic
#endif

      USE OCEAN, only : im,jm,lmo,dxypo,mo,lmm,imaxj,oXYP
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,mo,s0m,sxmo,symo,szmo,oc_tracer_mean
#endif
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst,msist,ssist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
#ifdef TRACERS_WATER
     *     ,trsist
#endif
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, haveLatitude,
     *     broadcast, GLOBALSUM
      USE OCEANR_DIM, only : grid=>ogrid
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, pack_data, unpack_data
      !USE OCEAN, only : gather_ocean
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      USE Dictionary_mod
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
c
      integer n,i,j,l,nst,i1,j1,i2,j2,ll
#ifdef TRACERS_SPECIAL_O18
      integer iu_O18ic,ip1,im1
      character*80 title
      real*4, dimension(im,jm,lmo) :: t0m4,tzm4
      real*8 fwm,afac
#endif
      real*8 t01,t02,trsum,wsum,tratio,frac_tr
      real*8, dimension(im,jm,lmo) :: mo_glob,s0m_glob,trmo_glob
      real*8 :: OTRACJ(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) 
      INTEGER :: J_0S, J_1S, J_0, J_1, J_0H, J_1H
      type(ocn_tracer_entry), pointer :: entry

      call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     *     J_STRT = J_0, J_STOP = J_1, 
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

#ifdef TRACERS_SPECIAL_O18
      if (is_set_param('water_tracer_ic')) then
        call get_param('water_tracer_ic', water_tracer_ic)
      endif
      !call sync_param ( "water_tracer_ic",water_tracer_ic  )
#endif

C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      call pack_data(grid,mo,mo_glob)
      call pack_data(grid,s0m,s0m_glob)

      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if (entry%trname.eq.'OceanAge') n_age=n
        if (itime.eq.entry%itime_tr0) then
        select case (entry%trname(1:6))

        case default
#ifdef TRACERS_OCEAN
          trmo(:,:,:,n)=0.
          txmo(:,:,:,n)=0.
          tymo(:,:,:,n)=0.
          tzmo(:,:,:,n)=0.
C**** straits
          if (am_i_root()) then
            trmst(:,:,n)=0.
            txmst(:,:,n)=0.
            tzmst(:,:,n)=0.
          end if
          CALL broadcast(grid, trmst)
          CALL broadcast(grid, txmst)
          CALL broadcast(grid, tzmst)

#endif

#ifdef TRACERS_WATER
          if (am_i_root()) trsist(:,:,n)=0.
          CALL broadcast(grid, trsist)
#endif

#if (defined TRACERS_OCEAN) && (defined TRACERS_ZEBRA)
        case ('zebraL')

           read(entry%trname(7:8),'(I2)') ll ! gets level from name
           do l=1,lmo
             do j=J_0,J_1
               do i=1,im
                 if (l.eq.ll .and. l.le.lmm(i,j)) then
                   trmo(i,j,l,n) = mo(i,j,l)*dxypo(j) ! set conc=1 for l=ll
                 else
                   trmo(i,j,l,n) = 0.   ! zero otherwise
                 end if
               end do
             end do
           end do

           if (am_i_root()) then
             do nst=1,nmst
               do l=1,lmst(nst)
                 if (l.eq.ll) then
                   trmst(nst,l,n)= mmst(nst,l)
                 else
                   trmst(nst,l,n)= 0.
                 endif
               end do
             end do
           end if

! set all gradients to zero initially
           txmo(:,:,:,n) = 0; tymo(:,:,:,n)= 0. ; tzmo(:,:,:,n)=0.
           if (am_i_root()) then
             txmst(:,:,n)=0. ; tzmst(:,:,n)=0.
          endif
          CALL broadcast(grid, trmst)
          CALL broadcast(grid, txmst)
          CALL broadcast(grid, tzmst)
#endif

        case ('Water', 'H2O18', 'HDO', 'H2O17' )
#if (defined TRACERS_OCEAN) && (defined TRACERS_SPECIAL_O18)
C**** Open ic file for isotope tracers
          if(water_tracer_ic.eq.1)
     *         call openunit("H2O18ic",iu_O18ic,.true.,.true.)
#endif

#ifdef TRACERS_OCEAN
#ifndef TRACERS_SPECIAL_O18
C**** main ocean variabiles and gradients
          do j=J_0,J_1
          do i=1,im
            do l=1,lmm(i,j)
              trmo(i,j,l,n)=entry%trw0*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
              txmo(i,j,l,n)=-entry%trw0*sxmo(i,j,l)
              tymo(i,j,l,n)=-entry%trw0*symo(i,j,l)
              tzmo(i,j,l,n)=-entry%trw0*szmo(i,j,l)
            end do
          end do
          end do
#else
C**** read in initial conditions for isotopes
C**** search through for correct tracer (since there is no guarantee
C**** that they will be in same order as tracer indices).
C**** data are now in 'per mil'
          if(water_tracer_ic.eq.1) then
            rewind (iu_O18ic)
 10         read  (iu_O18ic,err=800,end=810) title,t0m4,tzm4
            if (index(title,trim(entry%trname)).eq.0) goto 10
            write (6,*) 'Read from H2O18ic: ',title
            call closeunit(iu_O18ic)
          else
            t0m4(:,:,:)=0.
            tzm4(:,:,:)=0.
c            if(n.eq.n_Water) t0m4(:,:,:)=1.
          endif

C**** Turn per mil data into mass ratios (using current standard)
          t0m4(:,:,:)=(t0m4(:,:,:)*1d-3+1.)*entry%trw0
          tzm4(:,:,:)=0.    ! tzm4(:,:,:)*1d-3*entry%trw0 corrupted?
C****
          do l=1,lmo
            txmo(:,J_0:J_1,l,n) = 0.
            tymo(:,J_0:J_1,l,n) = 0.
C**** Define East-West horizontal gradients
            im1=im-1
            i=im
            do j=J_0S,J_1S
              do ip1=1,im
                if (lmm(i,j).ge.l) then
                  if (lmm(im1,j).ge.l) then
                    if (lmm(ip1,j).ge.l) then
                      txmo(i,j,l,n)=2.5d-1*(t0m4(ip1,j,l)-t0m4(im1,j,l))
                    else
                      txmo(i,j,l,n)=  5d-1*(t0m4(  i,j,l)-t0m4(im1,j,l))
                    end if
                  else
                    if (lmm(ip1,j).ge.l)
     *                   txmo(i,j,l,n)=5d-1*(t0m4(ip1,j,l)-t0m4(i,j,l))
                  end if
                end if
                im1=i
                i=ip1
              end do
            end do
C**** Define North-South horizontal gradients
            do j=J_0S,J_1S
              do i=1,im
                if (lmm(i,j).ge.l)  then
                  if (lmm(i,j-1).ge.l)  then
                    if (lmm(i,j+1).ge.l)  then
                      tymo(i,j,l,n)=2.5d-1*(t0m4(i,j+1,l)-t0m4(i,j-1,l))
                    else
                      tymo(i,j,l,n)=  5d-1*(t0m4(i,  j,l)-t0m4(i,j-1,l))
                    end if
                  else
                    if (lmm(i,j+1).ge.l)
     *                   tymo(i,j,l,n)=5d-1*(t0m4(i,j+1,l)-t0m4(i,j,l))
                  end if
                end if
              end do
            end do
C**** Multiply ratios by freshwater mass
            do j=J_0,J_1
            do i=1,im
              if (l.le.lmm(i,j)) then
                fwm = mo(i,j,l)*dxypo(j)-s0m(i,j,l)
                trmo(i,j,l,n)=t0m4(i,j,l)*fwm
                txmo(i,j,l,n)=txmo(i,j,l,n)*fwm-sxmo(i,j,l)*t0m4(i,j,l)
                tymo(i,j,l,n)=tymo(i,j,l,n)*fwm-symo(i,j,l)*t0m4(i,j,l)
                tzmo(i,j,l,n)=tzm4(i,j,l)  *fwm-szmo(i,j,l)*t0m4(i,j,l)
              else
                trmo(i,j,l,n)=0.
                txmo(i,j,l,n)=0.
                tymo(i,j,l,n)=0.
                tzmo(i,j,l,n)=0.
              end if
            end do
            end do
          end do
#endif

C**** Initiallise strait values based on adjacent ocean boxes
          !call gather_ocean(1)  ! mo,g0m,gx-zmo,s0m,sx-zmo,trmo,tx-zmo

          call pack_data(grid,trmo(:,:,:,n),trmo_glob)

          if(am_I_root()) then
          do nst=1,nmst
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=trmo_glob(i1,j1,l)/(mo_glob(i1,j1,l)*dxypo(j1)
     *             -s0m_glob(i1,j1,l))
              t02=trmo_glob(i2,j2,l)/(mo_glob(i2,j2,l)*dxypo(j2)
     *             -s0m_glob(i2,j2,l))
              trmst(l,nst,n) = 5d-1*(mmst(l,nst)-s0mst(l,nst))*(t01+t02)
              txmst(l,nst,n) = -5d-1*sxmst(l,nst)*(t01+t02)
              tzmst(l,nst,n) = -5d-1*szmst(l,nst)*(t01+t02)
            end do
            do l=lmst(nst)+1,lmo
              trmst(l,nst,n) = 0.
              txmst(l,nst,n) = 0.
              tzmst(l,nst,n) = 0.
            end do
#ifdef TRACERS_WATER
            trsist(n,1:2,nst) = entry%trw0*(msist(1,nst)*xsi(1:2)
     *           -ssist(1:2,nst))
            trsist(n,3:lmi,nst)=entry%trw0*(msist(2,nst)*xsi(3:lmi)
     *           -ssist(3:lmi,nst))
#endif
          end do
          end if

          call bcast_straits(.false.) ! bcst tracers

C**** Balance tracers so that average concentration is TRW0 
C**** or oc_tracer_mean
          
          CALL CONSERV_OTR(OTRACJ,N)
          CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

          if (AM_I_ROOT()) then
            if (oc_tracer_mean(n).ne.-999.) then
              tratio=entry%trw0*(oc_tracer_mean(n)*1d-3+1.)
            else
              tratio=entry%trw0
            end if
            
            if (entry%trname.eq.'Water') then
              wsum = trsum
              frac_tr = -999.
            else
              frac_tr = tratio/(trsum/wsum)
              write(6,*) "Average oceanic tracer concentration ",
     *           entry%trname,(trsum/(wsum*entry%trw0)-1d0)*1d3,frac_tr
            end if
          end if
          call broadcast(grid,  frac_tr)
         
          if (frac_tr.ne.-999.) then    
            do l=1,lmo
              do j=J_0,J_1
                do i=1,imaxj(j)
                  trmo(i,j,l,n) = trmo(i,j,l,n) * frac_tr
                  if (frac_tr.lt.1.) then
                    TXMO(I,J,L,n)=frac_tr*TXMO(I,J,L,n)
                    TYMO(I,J,L,n)=frac_tr*TYMO(I,J,L,n)
                    TZMO(I,J,L,n)=frac_tr*TZMO(I,J,L,n)
                  end if
                end do
              end do
            end do
C**** straits
            if (am_i_root()) then
              DO NST=1,NMST
                DO L=1,LMST(NST)
                  TRMST(L,NST,N)=frac_tr*TRMST(L,NST,N)
                  if (frac_tr.lt.1) then
                    TXMST(L,NST,N)=frac_tr*TXMST(L,NST,N)
                    TZMST(L,NST,N)=frac_tr*TZMST(L,NST,N)
                  end if
                END DO
              END DO
            end if
            CALL broadcast(grid, TRMST)
            CALL broadcast(grid, TXMST)
            CALL broadcast(grid, TZMST)

C**** Check
            CALL CONSERV_OTR(OTRACJ,N)
            CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

            if (AM_I_ROOT()) then
              tratio=(trsum/(wsum*entry%trw0)-1.)*1000.
              write(6,*) "New ocean tracer mean: ",tratio
     *             ,oc_tracer_mean(n)
            end if
          end if
#endif
C****
        end select
        write(6,*) entry%trname," tracer initialised in ocean"
        end if
      end do

C**** ensure that atmospheric arrays are properly updated (i.e. gtracer)
      ! call now made at higher level
      !CALL TOC2SST(atmocn)

      return
 800  write(6,*) "Error reading input file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
 810  write(6,*) "Tracer ",entry%trname," not found in file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      SUBROUTINE OC_TDECAY(DTS)
!@sum OC_TDECAY decays radioactive tracers in ocean
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only:tracerlist, ocn_tracer_entry,expDecayRate
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      IMPLICIT NONE
      real*8, intent(in) :: dts
      integer n
      type(ocn_tracer_entry), pointer :: entry

      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if (entry%trdecay.gt.0. .and. itime.ge.entry%itime_tr0) then
C**** Oceanic decay
          trmo(:,:,:,n)   = expDecayRate(n)*trmo(:,:,:,n)
          txmo(:,:,:,n)   = expDecayRate(n)*txmo(:,:,:,n)
          tymo(:,:,:,n)   = expDecayRate(n)*tymo(:,:,:,n)
          tzmo(:,:,:,n)   = expDecayRate(n)*tzmo(:,:,:,n)
        end if
      end do
C****
      return
      end subroutine oc_tdecay

      SUBROUTINE OCN_TR_AGE(DTS)
!@sum OCN_TR_AGE age tracers in ocean
!@auth Gavin Schmidt/Natassa Romanou
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_age
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj,
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
C**** age=1/(INT_DAYS_PER_YEAR*24*3600) in years
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
#endif

      SUBROUTINE DIAGTCO (M,NT0,atmocn)
!@sum  DIAGTCO Keeps track of the conservation properties of ocean tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
      USE OCEAN, only : IMO=>IM, oJ_BUDG, oJ_0B, oJ_1B
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      USE OCEANR_DIM, only : oGRID
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT0 index denoting tracer number, NT is index in tconsrv 
      INTEGER, INTENT(IN) :: nt0
!@var atmocn
      type(atmocn_xchng_vars) :: atmocn
c
      INTEGER :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(IMO,
     &                  oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(oJ_0B:oJ_1B) :: TOTALJ
      INTEGER :: nm,ni
      INTEGER :: I, J, J_0, J_1, I_0, I_1

      call getDomainBounds(ogrid, J_STRT=J_0, J_STOP=J_1)
      I_0 = oGRID%I_STRT
      I_1 = oGRID%I_STOP
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.
      nt=nt0+atmocn%natmtrcons
      if (atmocn%nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL (kg)
        call conserv_otr(total,nt0)
        nm=atmocn%nofmt(m,nt)
        ni=atmocn%nofmt(1,nt)

C**** Calculate zonal sums
        totalj(oj_0b:oj_1b)=0.
        do j=j_0,j_1
        do i=i_0,i_1
          totalj(oj_budg(i,j)) = totalj(oj_budg(i,j)) + total(i,j)
        end do
        end do

c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          atmocn%tconsrv(oJ_0b:oJ_1b,nm,nt) =
     &         atmocn%tconsrv(oJ_0b:oJ_1b,nm,nt)
     &         +(totalj(oJ_0b:oJ_1b)-atmocn%tconsrv(oJ_0b:oJ_1b,ni,nt)) 
        end if
C**** Save current value in TCONSRV(NI)
        atmocn%tconsrv(oJ_0b:oJ_1b,ni,nt)=totalj(oJ_0b:oJ_1b)
      end if
      return
      end subroutine diagtco

      SUBROUTINE conserv_OTR(OTR,ITR)
!@sum  conserv_OTR calculates zonal ocean tracer (kg) on ocean grid
!@auth Gavin Schmidt
      USE OCEAN, only : IMO=>IM,JMO=>JM,oXYP,LMM, imaxj,trmo
      USE STRAITS, only : nmst,jst,ist,lmst,trmst
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      USE OCEANR_DIM, only : ogrid

      IMPLICIT NONE
!@var OTR zonal ocean tracer (kg)
      REAL*8, DIMENSION(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: OTR 
      INTEGER, INTENT(IN) :: ITR
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds(ogrid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          OTR(I,J) = SUM(TRMO(I,J,:LMM(I,J),ITR))
        END DO
      END DO
      if (HAVE_SOUTH_POLE) OTR(2:IMO,1)   = OTR(1,1)
      if (HAVE_NORTH_POLE) OTR(2:IMO,JMO) = OTR(1,JMO)
C**** include straits variables
      DO N=1,NMST
        I=IST(N,1)
        J=JST(N,1)
        If (J >= J_0 .and. J <= J_1)
     &       OTR(I,J)=OTR(I,J)+SUM(TRMST(:LMST(N),N,ITR))
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OTR
