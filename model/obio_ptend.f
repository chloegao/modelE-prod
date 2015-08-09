#include "rundeck_opts.h"

      subroutine obio_ptend(vrbos,kmax,i,j)

c  Computes tendencies of biological particles (units/hr) 
c  tracer
c  P(1) = nitrate (uM)
c  P(2) = ammonium (uM)
c  P(3) = silica (uM)
c  P(4) = iron (nM)
c  P(5) = diatoms (mg chl m-3)
c  P(6) = chlorophytes (mg chl m-3)
c  P(7) = cyanobacteria (mg chl m-3)
c  P(8) = coccolithophores (mg chl m-3)
c  P(9) = herbivores (mg chl m-3)
 
      USE obio_dim
      USE obio_incom,only: cnratio,cfratio,remin,obio_wsh,bf,cchlratio
     .                    ,wsdeth,rkn,rks,rkf,Rm,phygross,bn,solFe
     .                    ,mgchltouMC,uMtomgm3
      USE obio_forc, only: tirrq
      USE obio_com, only : dp1d,obio_P,obio_ws,P_tend,D_tend,C_tend
     .                    ,gro,rlamz,dratez1,dratez2
     .                    ,greff,pnoice,drate,tfac,regen,Fescav
     .                    ,wshc,rikd,rmuplsr,det
     .                    ,gcmax1d,covice_ij,atmFe_ij
     .                    ,temp1d,wsdet,tzoo,p1d
     .                    ,rhs,pp2_1d,flimit
#ifdef restoreIRON
     .                    ,Iron_BC
#endif

#ifdef OBIO_RUNOFF
#ifdef POC_RUNOFF
     .                      ,rpocconc_loc
#endif
#ifdef NITR_RUNOFF
     .                      ,rnitrconc_loc
!                         ,rnitrfmlo_loc
#endif
#ifdef SILI_RUNOFF
     .                      ,rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
     .                      ,rironconc_loc
      USE obio_incom, only: estFe
#endif
#endif

#ifdef OBIO_RUNOFF
      USE OFLUXES, only:  oFLOWO
!      USE ocean, only:  oxyp
      USE MODEL_COM, only:  dtsrc 
#endif
#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo
      USE MODEL_COM, only : nstep=>itime
#else  /* HYCOM ocean */
      USE hycom_dim_glob, only : kdm
      USE hycom_scalars, only : nstep
#endif

      implicit none

      integer i,j,k,kto
      integer nt,kmax

      real rmu4(nchl)     !growth on ammonium, NH4
      real rmu3(nchl)     !growth on nitrate, NO4
      real zoo(ntyp)      !herbivores (zooplankton)
      real dphy(ntyp)     !death rate of phytoplankton
      real viscfac(kdm)



      real bs,ptot
      real Pzoo,Gzoo,Dzoo1,Dzoo2,exc
      real tirrqice,upn,upa,upf,ups

      real rnut2,tnit,tmp,rnut1,rmmn,framm,rmml,rmmlice,rmms
      real rmmf,rlim,rlimice,grate,rlimnfix,rlimrkn,rfix,gratenfix1
      real graterkn,gratenlim,gratenfix,gron,gronfix(kdm)
      real dphyt

      logical vrbos

      real term
      character*7, dimension(12) ::  bio_var = (/
     .  ' NO3   ',' NH4   ',' SiO2  ',' Iron  ','Diatom ','Chlphy ',
     .  'Cyanob ','Coccol ','Herbiv ','N/Cdet ','Silica ','Fe_det '/)


       rhs=0.0
       obio_ws = 0.0
       P_tend = 0.0
       C_tend = 0.0
       D_tend = 0.0
       wsdet = 0.0
       rmu4 = 0.0
       rmu3 = 0.0
       zoo  = 0.0
       dphy = 0.0
       gro = 0.0
       viscfac = 0.0
       pp2_1d = 0.0
       gronfix = 0.0



!define no ice points based on covice (here: covice_ij)
      pnoice(1)=1.-covice_ij
      do k=2,kdm
#ifdef change_PNOICE
         pnoice(k)=1.
#else
         pnoice(k)=pnoice(1)
#endif
      enddo

      bs = 2.0*bn

!tendency terms are computed in the mid-level (m)

c  Start Model Space Loop
!      m = indext2     !index of "past" (t-1)

!Iron + atm iron: disperse in layer and convert to nM
!we do not need to multiply this by pnoice, because iron
!is deposited over ice and presumably when ice melts will enter
!the ocean
       k = 1
       term = atmFe_ij*solFe*1.d-3/max(p1d(k+1),1.e-3)
       rhs(k,4,4) = term
       P_tend(k,4) = term

!change: March 15, 2010
!      do k=2,kmax
!         rhs(k,4,4) = 0.0
!         P_tend(k,4) = 0.0
!      enddo

!Grazing/regeneration of ammonium
       do k = 1,kmax

         ptot = 0.0
         do nt = nnut+1,ntyp-nzoo
          ptot = ptot + obio_P(k,nt)
         enddo
         ptot = max(ptot,1.0E-36)

         Pzoo = obio_P(k,ntyp)
         gzoo = tzoo*Rm*(1.0-exp(-rlamz*ptot))*Pzoo
         dzoo1 = dratez1*Pzoo
         dzoo2 = dratez2*Pzoo*Pzoo

         term = ((1.0-greff)*gzoo-dzoo1-dzoo2) * pnoice(k)
         rhs(k,ntyp,ntyp) = term    !!!!!!!this also takes into accout P(5,6,7,8)
         P_tend(k,ntyp) = term

cdiag    if (vrbos)
cdiag.   write(*,*)'ptend1: ',
cdiag.        nstep,i,j,k,ntyp,obio_P(k,ntyp),tzoo,Rm,rlamz,
cdiag.        ptot,Pzoo,dratez1,dratez2,greff,P_tend(k,ntyp)

         do nt = nnut+1,ntyp-nzoo
          !fraction of grazing for this group
          zoo(nt) = gzoo*obio_P(k,nt)/ptot
          dphy(nt) = drate*obio_P(k,nt)

          term = -zoo(nt) * pnoice(k)
          rhs(k,nt,ntyp) = term
          P_tend(k,nt) = term      !!nonlin term takes into accountP(ntyp),P(5:8)

cdiag     if (vrbos) write(*,*)
cdiag.     'obio_ptend2: ',nstep,nt,i,j,k,gzoo,obio_P(k,nt),ptot,
cdiag.      zoo(nt),P_tend(k,nt)

          term = -dphy(nt) * pnoice(k)           !death of phytoplankton
          rhs(k,nt,nt) = term
          P_tend(k,nt) = P_tend(k,nt) + term

cdiag     if (vrbos) write(*,*)
cdiag.     'obio_ptend3: ',nstep,nt,i,j,k,drate,dphy(nt),
cdiag.     obio_P(k,nt),P_tend(k,nt)

cdiag      if (vrbos) then
cdiag       if(k.eq.1 .and. nt.eq.nnut+1) write(101,'(a,a)')
cdiag.      'nstep    k   nt   dp1d   zoo(nt) dphy(nt)   obio_P(k,nt)'
cdiag.     ,'gzoo   ptot    drate  P_tend(k,nt)'
cdiag       write(101,'(3i5,8(1x,es8.2))')nstep,k,nt,dp1d(k),
cdiag.        zoo(nt),dphy(nt),obio_P(k,nt),gzoo,ptot,drate,
cdiag.        P_tend(k,nt)
cdiag      endif

         enddo

         exc = greff*gzoo

         !remineralization 
         term = tfac(k)*remin(1)*det(k,1)/cnratio * pnoice(k)
         rhs(k,1,ntyp+1) = term    !put this in diff column
         P_tend(k,1) = term

         !regeneration from zooplankton
         term = bn*(exc + regen*dzoo2) * pnoice(k)
         rhs(k,2,ntyp) = term
         P_tend(k,2) = term

         !remineralization 
         term = tfac(k)*remin(2)*det(k,2) * pnoice(k)
         rhs(k,3,ntyp+2) = term    !put this in diff column
         P_tend(k,3) = term

         !regeneration from zooplankton
         term = bf*(exc + regen*dzoo2) * pnoice(k)    
         rhs(k,4,ntyp) = term
         P_tend(k,4) = P_tend(k,4) + term

         !remineralization 
         term = tfac(k)*remin(3)*det(k,3) * pnoice(k)
         rhs(k,4,ntyp+3) = term                        !put this in diff column
         P_tend(k,4) = P_tend(k,4) + term

         term = -Fescav(k)
         rhs(k,4,13) = term             !this should actually be rhs(4,4) but we
                                        !have already defined it so let it be rhs(4,13)
         P_tend(k,4) = P_tend(k,4) + term


         term = exc*mgchltouMC*pnoice(k)
     *                  + regen*dzoo2*mgchltouMC*pnoice(k)
         rhs(k,13,ntyp) = term             
         C_tend(k,1) = term
 
cdiag    if (vrbos) then
cdiag     if(k.eq.1) write(102,'(a,a)')
cdiag.    'nstep    k  dp1d(k) tfac(k)   det(k,1) bn    exc   ',
cdiag.    ' det(3)   P_t(1)   P_t(2)   P_t(3) P_t(4) '
cdiag       write(102,'(2i5,10(1x,es8.2))')nstep,k,
cdiag.         dp1d(k),tfac(k),det(k,1),bn,exc,det(k,3),
cdiag.         P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4)
cdiag    endif
cdiag    if (vrbos) write(*,*)'ptend4: ',
cdiag.     nstep,i,j,k,dp1d(k),tfac(k),det(k,1),bn,exc,det(k,3),
cdiag.     P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4)


         dphyt = 0.0
         do nt = nnut+1,ntyp-nzoo
          dphyt = dphyt + dphy(nt)
         enddo
 
!1st detrital fraction is carbon
         term = dphyt*cchlratio * pnoice(k)
         rhs(k,ntyp+1,nnut+1) = term
         D_tend(k,1) = term

         term = dzoo1*cchlratio * pnoice(k)
     .        + (1.0-regen)*dzoo2*cchlratio * pnoice(k)
         rhs(k,ntyp+1,ntyp) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = -tfac(k)*remin(1)*det(k,1) * pnoice(k)
         rhs(k,ntyp+1,ntyp+1) = term
         D_tend(k,1) = D_tend(k,1) + term

!2nd detrital fraction is silica
         term = bs*dphy(nnut+1) * pnoice(k)
         rhs(k,ntyp+2,nnut+1) = term
         D_tend(k,2) = term

         term = bs*zoo(nnut+1) * pnoice(k)
         rhs(k,ntyp+2,ntyp) = term
         D_tend(k,2) = D_tend(k,2) + term

         term = -tfac(k)*remin(2)*det(k,2) * pnoice(k)
         rhs(k,ntyp+2,ntyp+2) = term
         D_tend(k,2) = D_tend(k,2) + term

!3rd detrital fraction is iron
         term = bf*dphyt * pnoice(k)
         rhs(k,12,nnut+1) = term
         D_tend(k,3) = term

         term = bf*dzoo1 * pnoice(k)
     .        + bf*(1.0-regen)*dzoo2 * pnoice(k)
         rhs(k,12,ntyp) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = -tfac(k)*remin(3)*det(k,3) * pnoice(k)
         rhs(k,12,ntyp+3) = term
         D_tend(k,3) = D_tend(k,3) + term

!change June 1, 2010
         term = 0.1 * Fescav(k)
!endofchange
         rhs(k,12,4) = term
         D_tend(k,3) = D_tend(k,3) + term

cdiag    if (vrbos) then
cdiag     if(k.eq.1) write(103,'(a,a)')
cdiag.    'nstep    k  dp1d(k)   dphyt   dzoo1     dzoo2   dphy(5) '
cdiag.   ,'zoo(5)   tfac(k)  D_t(1)     D_t(2) D_t(3)'
cdiag    write(103,'(2i5,10(1x,es8.2))')nstep,k,dp1d(k),
cdiag.   dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)
cdiag    endif
cdiag    if (vrbos) write(*,*)'ptend5 :',
cdiag.   nstep,i,j,k,
cdiag.   dp1d(k),dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)

       enddo !kmax

c Day: Grow
      do k = 1,kmax

cdiag if(vrbos)
cdiag.write(*,*)'obio_ptend6: ',nstep,i,j,k,tirrq(k)

CHECK if (tirrq(k) .gt. 0.0)then
        tirrqice = tirrq(k)*0.01  !reduce light in ice by half

c Light-regulated growth
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then

        ! Diatoms
        nt = 1

cdiag   if (vrbos) then
cdiag     if (k.eq.1)write(104,'(a,3i5)')
cdiag.                'ptend: diatoms growth, nstep,k,nt=',nstep,k,nt
cdiag     if (k.eq.1)write(104,'(5x,a,12x,a)')
cdiag.      'dp    tirrq    P(1)  P(2)  P(3)  P(4)','rkn    rks    rkf'
cdiag     write(104,'(9(1x,es8.2))')dp1d(k),tirrq(k),
cdiag.                       obio_P(k,1),obio_P(k,2),obio_P(k,3),
cdiag.                       obio_P(k,4),rkn(nt),rks(nt),rkf(nt)
cdiag   endif
cdiag   if (vrbos)
cdiag     write(*,'(a,4i5,9e12.4)')'ptend7 :', 
cdiag.              nstep,i,j,k,dp1d(k),tirrq(k),
cdiag.                       obio_P(k,1),obio_P(k,2),obio_P(k,3),
cdiag.                       obio_P(k,4),rkn(nt),rks(nt),rkf(nt)

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
         tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))     !nitrate
          tmp = 1.0 - rnut2
        ! Enforce preferential utilization of ammonium
          rnut1 = min(tmp,tnit)
           rmmn = rnut1 + rnut2         !nitrate limitation
          if ( rmmn .ne. 0.d0 ) then
            framm = rnut2/rmmn
          else
            framm = 0.d0
          endif
           rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))         !light limitation
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))

        rmms = obio_P(k,3)/(rks(nt)+obio_P(k,3))      !silicate limitation
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron limitation
        rlim = min(rmml,rmmn,rmms,rmmf)
        rlimice = min(rmmlice,rmmn,rmms,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml*pnoice(k)
        flimit(k,nt,2) = rmmlice*(1-pnoice(k))
        flimit(k,nt,3) = rmmn*pnoice(k)  
        flimit(k,nt,4) = rmms*pnoice(k)  
        flimit(k,nt,5) = rmmf*pnoice(k)  

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
     .          + rmuplsr(k,nt)*rlimice*(1.0-pnoice(k))
        rmu4(nt) = grate*framm
        rmu3(nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production  in mgC/m2/hr because:
        ! [gro]= mg,chl/m3/hr, [dp]= m, 
        ! [cchlratio]= mgl/mgl,[phygross]=no units 
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio

!       write(*,'(a,5i5,11e12.4)')'obio_ptend, pp:',
!    .  nstep,i,j,k,nt,gro(k,nt), phygross,p1d(k+1),cchlratio
!    . ,grate,obio_P(k,nt+nnut),rmuplsr(k,nt),rlim,rmml,rmmn,rmmf
      endif

!!#endif


!!#if NCHL_DEFINED > 1
      if (nchl > 1) then

! Chlorophytes
        nt = 2

cdiag   if (vrbos) then
cdiag     if (k.eq.1)write(105,'(a,3i5)')
cdiag.                'ptend: chlorop growth, nstep,k,nt=',nstep,k,nt
cdiag     if (k.eq.1)write(105,'(5x,a,12x,a)')
cdiag.      'dp    tirrq    P(1)  P(2)  P(3)  P(4)','rkn    rks    rkf'
cdiag     write(105,'(9(1x,es8.2))')dp1d(k),tirrq(k),
cdiag.               obio_P(k,1),obio_P(k,2),obio_P(k,4),rkn(nt),
cdiag.                       rkf(nt), gro(k,nt),obio_P(k,nt+nnut)
cdiag   endif
cdiag   if (vrbos)
cdiag.    write(*,*)'ptend8: ',nstep,i,j,k,dp1d(k),tirrq(k),
cdiag.               obio_P(k,1),obio_P(k,2),obio_P(k,4),rkn(nt),
cdiag.                       rkf(nt), gro(k,nt),obio_P(k,nt+nnut)


        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit  = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp   = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2         !nitrate limitation                
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn                           
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))     !light limitation
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron limitation
        rlim = min(rmml,rmmn,rmmf)
        rlimice = min(rmmlice,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim * pnoice(k)
     .        + rmuplsr(k,nt)*rlimice * (1.0-pnoice(k))
        rmu4(nt) = grate*framm
        rmu3(nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term  
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio

!       write(*,'(a,5i5,11e12.4)')'obio_ptend, pp:',
!    .  nstep,i,j,k,nt,gro(k,nt), phygross,p1d(k+1),cchlratio
!    . ,grate,obio_P(k,nt+nnut),rmuplsr(k,nt),rlim,rmml,rmmn,rmmf
      endif
!!#endif


!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
! Cyanobacteria
        nt = 3
        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)
        rlimnfix = min(rmml,rmmf)         !limitation for N2 fixation
        rlimrkn = min(rmml,rkn(nt),rmmf)   !limitation at kn

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
        rmu4(nt) = grate*framm
        rmu3(nt) = grate*(1.0-framm)
        rfix = 0.25*exp(-(75.0*obio_P(k,nt+nnut)))
        rfix = max(rfix,0.0)
c        rfix = min(rfix,0.2)

        gratenfix1 = rmuplsr(k,nt)*rlimnfix*rfix  !N fix
        graterkn = rmuplsr(k,nt)*rlimrkn
        gratenlim = graterkn - grate
        gratenfix = min(gratenlim,gratenfix1)  !N fix cannot exceed kn
        gratenfix = max(gratenfix,0.0) * pnoice(k)

        gron = grate*obio_P(k,nt+nnut)
        gronfix(k) = gratenfix*obio_P(k,nt+nnut)
        gro(k,nt) = gron + gronfix(k)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

cdiag   if (vrbos) write(*,*)
cdiag.   'obio_ptend9: ',nt,i,j,k,rmuplsr(k,nt),rlim,rlimnfix,rlimrkn,
cdiag.    rfix,gratenfix,gronfix,obio_P(k,nt+nnut),gron,gro(k,nt),
cdiag.    P_tend(k,nt+nnut)

cdiag   if (vrbos) write(*,*)
cdiag.  'obio_ptend10: ',nt,i,j,k,tirrq(k),rikd(k,nt),rmml,obio_P(k,4),
cdiag.   rkf(nt),rmmf,rkn(nt),rnut1,rnut2,tmp,tnit,rmmn,obio_P(k,1),
cdiag.   obio_P(k,2)

        !Net primary production
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio

!       write(*,'(a,5i5,11e12.4)')'obio_ptend, pp:',
!    .  nstep,i,j,k,nt,gro(k,nt), phygross,p1d(k+1),cchlratio
!    . ,grate,obio_P(k,nt+nnut),rmuplsr(k,nt),rlim,rmml,rmmn,rmmf
      endif
!!#endif


!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
! Coccolithophores
        nt = 4

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
        rmu4(nt) = grate*framm 
        rmu3(nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt) * pnoice(k)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term
        gcmax1d(k) = max(gcmax1d(k),grate)

        !Net primary production
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio

!       write(*,'(a,5i5,11e12.4)')'obio_ptend, pp:',
!    .  nstep,i,j,k,nt,gro(k,nt), phygross,p1d(k+1),cchlratio
!    . ,grate,obio_P(k,nt+nnut),rmuplsr(k,nt),rlim,rmml,rmmn,rmmf
      endif
!!#endif


!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
! Dinoflagellates
        nt = 5

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim
        rmu4(nt) = grate*framm * pnoice(k)
        rmu3(nt) = grate*(1.0-framm) * pnoice(k)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio 

!       write(*,'(a,5i5,11e12.4)')'obio_ptend, pp:',
!    .  nstep,i,j,k,nt,gro(k,nt), phygross,p1d(k+1),cchlratio
!    . ,grate,obio_P(k,nt+nnut),rmuplsr(k,nt),rlim,rmml,rmmn,rmmf
      endif
!!#endif

cdiag   if (vrbos) then
cdiag   do nt = 1,nchl
cdiag     write(801,'(i7,2e12.4)')
cdiag.          nstep,gro(1,nt),obio_P(1,nt+nnut)
cdiag   enddo
cdiag   endif

cdiag  if(vrbos)write(*,*)'ptend11: ',
cdiag.  nstep,i,j,k,tirrq(k),P_tend(k,5),P_tend(k,6),P_tend(k,7) 
cdiag.  ,P_tend(k,8),P_tend(k,9)


! Nutrient uptake
        upn = 0.0
        upa = 0.0
        upf = 0.0
        do nt = 1,nchl
         term = bn*(rmu3(nt)*obio_P(k,nnut+nt))
         upn = upn + term
         rhs(k,1,nnut+nt) = -upn    !use - sign here because that is how it goes in P_tend

         term = bn*(rmu4(nt)*obio_P(k,nnut+nt))
         upa = upa + term
         rhs(k,2,nnut+nt) = -upa    !use - sign here because that is how it goes in P_tend

         term = bf*gro(k,nt)
         upf = upf + term
         rhs(k,4,nnut+nt) = -upf    !use - sign here because that is how it goes in P_tend
        enddo

        term = bs*gro(k,1)
        ups = term
        rhs(k,3,nnut+1) = -ups    !use - sign here because that is how it goes in P_tend

cdiag  if(vrbos)then
cdiag     if(k.eq.1) write(106,'(a,a)')
cdiag.'nstep    k  dp1d(k)   bn   rmu3(1)  rmu3(2)   rmu3(3)',
cdiag.'rmu3(4)  P(5)   P(6)   P(7)   P(8)     P_t(1)bef P_t(1) '
cdiag   write(106,'(2i5,12(1x,es8.2))')
cdiag. nstep,k,dp1d(k),bn,rmu3(1),rmu3(2),rmu3(3),
cdiag. rmu3(4),obio_P(k,nnut+1),obio_P(k,nnut+2),obio_P(k,nnut+3),
cdiag. obio_P(k,nnut+4),P_tend(k,1),P_tend(k,1) - upn
cdiag  endif

cdiag  if(vrbos)write(*,*)'obio_ptend12: ',
cdiag. nstep,i,j,k,dp1d(k),bn,rmu3(1),rmu3(2),rmu3(3),rmu3(4),
cdiag. obio_P(k,5),obio_P(k,6),obio_P(k,7),obio_P(k,8),
cdiag. P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4),
cdiag. upn,upa,upf,ups

        P_tend(k,1) = P_tend(k,1) - upn
        P_tend(k,2) = P_tend(k,2) - upa
        P_tend(k,3) = P_tend(k,3) - ups
        P_tend(k,4) = P_tend(k,4) - upf

CHECK  endif !tirrq(k) .gt. 0.0
      enddo  !kmax

!  Fix up nitrogen uptake from N-fixation -- give back to water
!  column in reverse order of layers
      do k = 1,kdm
       kto = (kdm-k)+1
       P_tend(kto,1) = P_tend(kto,1) - bn*gronfix(k)
      enddo

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!       if (oFLOWO(i,j) .gt. 0.)then
!	 rnitr_loc = rnitrmflo_loc
!     .    / (oFLOWO(i,j) * dxypo(j)) ! kg/s => kg,N/kg,w/s
!     .    * 1.d3     ! kg,N to g,N
!     .    * (1./14.) ! g,N to mol,N
!     .    * 1.d3     ! mol,N to mmol,N
!     .    * rho_water ! kg,water to m3 water
!         if (i.eq.169.and. j.eq.59) then
!           write(*,'(/,a,2i5,6e12.4)')'i,j,rnitrmflo, rnitr, oFLOWO,
!     .       dxypo, nitr, rho_water:',i,j,rnitrmflo_loc(i,j),
!     .       rnitr_loc(i,j),
!     .       oFLOWO(i,j),dxypo(j),obio_P(1,1),rho_water
!	else
!	  rnitr_loc = 0.
!	endif
         term = rnitrconc_loc(i,j) 
     .    * oFLOWO(i,j)/dtsrc         ! kg,N/kg,w => kg,N/m2,w/s
     .    / dp1d(1)                   ! kg,N/m2,w/s => kg,N/m3,w/s
     .    * 1.d6/14.                  ! kg,N/m3,w/s => mmol,N/m3,w/s
     .    * 3600.                     ! mmol,N/m3,w/s => mmol,N/m3,w/hr
         rhs(1,1,17) = term
         P_tend(1,1) = P_tend(1,1) + term

!	if (i.eq.169 .and. j.eq.59) then
!	  write(*,'(/,a,2i5,5e12.4)')'i,j,rnitrconc,rnitr,oFLOWO,
!     .       dtsrc,dp1d(1):',i,j,rnitrconc_loc(i,j),rhs(1,1,17),
!     .       oFLOWO(i,j),dtsrc,dp1d(1) 
!                endif

#endif
#ifdef SILI_RUNOFF
	term = rsiliconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,S/kg,w => kg,S/m2,w/s
     .   / dp1d(1)                   ! kg,S/m2,w/s => kg,S/m3,w/s
     .   * 1.d6/28.055               ! kg,S/m3,w/s => mmol,S/m3,w/s
     .   * 3600.                     ! mmol,S/m3,w/s => mmol,S/m3,w/hr
        rhs(1,3,17) = term
        P_tend(1,3) = P_tend(1,3) + term
#endif
#ifdef IRON_RUNOFF
	term = rironconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,Fe/kg,w => kg,Fe/m2,w/s
     .   / dp1d(1)                   ! kg,Fe/m2,w/s => kg,Fe/m3,w/s
     .   * 1.d9/55.845              ! kg,Fe/m3,w/s => umol,Fe/m3,w/s
     .   * 3600.                     ! umol,Fe/m3,w/s => umol,Fe/m3,w/hr
     .   * estFe                     ! estuarine retention rate
        rhs(1,4,17) = term
        P_tend(1,4) = P_tend(1,4) + term
#endif
#ifdef POC_RUNOFF
	term = rpocconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,C/kg,w => kg,C/m2,w/s
     .   / dp1d(1)                   ! kg,C/m2,w/s => kg,C/m3,w/s
     .   * 1.d6                      ! kg,C/m3,w/s => mg,C/m3,w/s
     .   * 3600.                     ! mg,C/m3,w/s => mg,C/m3,w/hr
        rhs(1,10,17) = term
        D_tend(1,1) = D_tend(1,1) + term
#endif
#endif

#ifdef restoreIRON
!iron bottom sink/source
!whether sink or source is determined from Iron_BC
!Iron_BC > 0 sink of iron through sedimentation
!Iron_BC < 0 source of iron through sediment resuspension (Moore et al, 2004) 
!     k = kmax
!       if (p1d(kmax) >= 3700.) then        !for deep regions, bottom cell (lower 200m)
!           term =  - obio_P(k,4) / (200/(Iron_BC*3700))
!           term = term/(365*24)     !convert to per hr
!           rhs(k,4,14) = term
!           P_tend(k,4) = P_tend(k,4) + term
!       endif
!alternative
        if (p1d(kmax) <= 1100.) then    !for shelf regions, bottom cell
          !!do k = 1,kmax
            k=kmax
            term =  50.d0/(24.d0*dp1d(k))     ! in nano-moleFe/m3/hr
            rhs(k,4,14) = term
            P_tend(k,4) = P_tend(k,4) + term
          !!enddo
!       write(*,'(a,5i5,3e12.4)')'obio_ptend, iron source:',
!    .  nstep,i,j,k,kmax,dp1d(kmax),P_tend(k,4)-term,term
        endif
#endif


#ifdef TRACERS_Alkalinity
      call obio_alkalinity(kmax,i,j)
#endif
  
cdiag if(vrbos) then
cdiag do k=1,kdm
cdiag  write(*,*)'obio_ptend13: ',
cdiag. nstep,i,j,k,tirrq(k),gro(k,1),gro(k,2),
cdiag. gro(k,3),gro(k,4),P_tend(k,8),P_tend(k,9)
cdiag enddo
cdiag endif

      call obio_carbon(gro,vrbos,kmax,i,j)

cdiag if(vrbos) then
cdiag do k=1,kdm
cdiag  write(*,*)'obio_ptend14: ',
cdiag.  nstep,i,j,k,P_tend(k,8),P_tend(k,9)
cdiag enddo
cdiag endif

 
cdiag if (vrbos) then
cdiag  do k=1,1
cdiag   write (*,107) k,(bio_var(l),l=1,9),
cdiag.   (bio_var(nt),obio_P(k,nt),P_tend(k,nt),
cdiag.    (rhs(k,nt,l),l=1,16),nt=1,9),
cdiag.   (bio_var(nt),det(k,nt-9),D_tend(k,nt-9),
cdiag.    (rhs(k,nt,l),l=1,6),nt=10,12)
cdiag
cdiag   write (*,107) k,(bio_var(l),l=7,12),
cdiag.   (bio_var(nt),obio_P(k,nt),P_tend(k,nt),
cdiag.    (rhs(k,nt,l),l=9,16),nt=1,9),
cdiag.   (bio_var(nt),det(k,nt-9),D_tend(k,nt-9),
cdiag.    (rhs(k,nt,l),l=7,12),nt=10,12)
cdiag  end do
cdiag end if
 107  format (/'lyr',i3,4x,'amount   tndcy   ',
     .    9(2x,a7)/(a7,2es9.1,2x,6es9.1))

!compute here rates but all detritus update done inside the update routine

c Sinking rate temperature (viscosity) dependence (also convert to /hr)
      do k = 1,kmax
        viscfac(k) = 0.451 + 0.0178*temp1d(k)
      enddo
      do nt = 1,nchl
       do k = 1,kmax
         obio_ws(k,nt) = obio_wsh(nt)*viscfac(k)*pnoice(k)
       enddo
      enddo
      nt = 4
      do k = 1,kmax
        obio_ws(k,nt) = wshc(k)*viscfac(k)*pnoice(k)
      enddo
      do nt = 1,nchl
        obio_ws(kmax+1,nt)=obio_ws(kmax,nt)
      enddo
      do nt = 1,ndet
       do k = 1,kmax
         wsdet(k,nt) = wsdeth(nt)*viscfac(k)*pnoice(k)
       enddo
      enddo

cdiag  if(vrbos .and. diagno)then
cdiag   do k=1,kmax
cdiag    write(905,'(2i5,11e12.4)')
cdiag.     nstep,k,temp1d(k),viscfac(k),pnoice,
cdiag.             (obio_wsh(nt),nt=1,nchl),(obio_ws(k,nt),nt=1,nchl)
cdiag      write(906,'(2i5,4e12.4)')
cdiag.       nstep,k,wshc(k),viscfac(k),pnoice,obio_ws(k,4)
cdiag      write(*,'(a,2i5,8e12.4)')'obio_ptend, wsdet ',
cdiag.       nstep,k,(wsdeth(nt),nt=1,ndet),viscfac(k),pnoice,
cdiag.              (wsdet(k,nt),nt=1,ndet)
cdiag   enddo
cdiag  endif

c Save method for hard boundary condition (no flux)
c      srate = 0.0 - obio_wsh(n)*tracer(i,k-1,m,n)

      return
      end subroutine obio_ptend
