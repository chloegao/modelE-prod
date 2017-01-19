#include "rundeck_opts.h"
      SUBROUTINE GET_CC_CDNC_MX(L,nmodes,ncaero,MCDNL1,MCDNO1)
!@sum specific calculation to get cloud droplet number for indirect effects
!@auth Surabi Menon 
!@contains routines for calculating cloud droplet number (cm-3) for convective clouds
!@this is called in CLOUDS2_E1 if MATRIX is used to set aerosols
      USE CLOUDS_COM
      IMPLICIT NONE
      real*8 AIRM,EXPL,EXPO,WCDNO,WCDNL,rho
     *,MCDNL1,MCDNO1,amass,tams,smturb,DXYPJ,PL,TL
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSM6,SSM7,SSM8,
     *       SSMAL,SSMAO,SSML,SSMO
      integer, PARAMETER :: nt=31 ! this MUST be equal to sntm
      real*8,dimension(nt)::DSS,DSU,ncaero
      integer L,n,nmodes,nm 

      SSMAL=0.d0
      SSMAO=0.d0

*************************************************************************************************
C** Can use different approaches to convert mass to number, very sensitive to sizes assumed
C** Converting aerosol mass (kg/m3) to number (m-3) from Leon Rotstayn, based on IPCC Table 5
C** This was used in the Menon and Rotstayn (2006) Clim. Dyn. paper
C**_____________________________________________________________________________________________
c Factor of 132.14/32.06 converts from sulfur to ammonium sulfate
c 1.69e11 converts from mass (kg/m3) to number conc. (/cm3) for assumed lognormal dist.
c 1.21e11 converts from hydrophilic mass (kg/m3) to number conc.(/cm3) for
c assumed lognormal distribution for carbonaceous aerosols.
C**_____________________________________________________________________________________________
c     SSM1 = 1.69d11*(132.14d0/32.06d0)*DSU(1)         !all sulfate 
c     SSM2=(DSU(2)/2169.d0)/(0.004189d0*(.44d0**3)   ) !seasalt 0.1-1 um range
c     SSM4 = 1.21d11*((DSU(4))+ DSU(6))                !OCI,BCI: aged (1 day efolding time)
c     SSM5 = 1.21d11*((0.8*DSU(5))+(0.6*DSU(7)))       !OCB,BCB: 80% and 60% as hydrophillic) 
c     SSM8 = 1.21d11*((DSU(18))+(DSU(19))+(DSU(20))+(DSU(21))) !SOA: 100% hydrophillic) 
!     SSMD1=(DSU(10)/2.5d3)/(0.004189d0*(0.75d0**3.))  !dust(clay) coated w/ sulf
!     SSMD2=(DSU(11)/2.65d3)/(0.004189d0*(2.2d0**3.))  !dust(silt1) coated w/ sulfate
!     SSMD3=(DSU(12)/2.65d3)/(0.004189d0*(4.4d0**3.))  !dust(silt2) coated w/ sulfate
C*** Old way of getting Na in cm-3 as from Lohmann et al. 1999, JGR
C** 4.189 is the conversion for 4/3*pi 
c     SSM1=(DSU(1)/1777.d0)/(4.189*1.d12*(0.15**3.))      !sulfate 
c     SSM2=(DSU(2)/2169.d0)/(4.189*1.d12*(0.44**3.))      !seasalt 0.1-1 um range
c     SSM3=(0.1d0*DSU(3)/2169.)/(0.004189*(1.7**3.))   !seasalt 1-4 um range    
c     SSM3=(DSU(4)/1000.d0)/(0.004189*(0.20**3.))      !OCIA                   
c     SSM4=(0.80*DSU(5)/1000.d0)/(0.004189*(0.20**3.)) !OCB                    
c     SSM5=(DSU(6)/1000.d0)/(0.004189*(0.09**3.))      !BCIA                   
c     SSM6=(0.60*DSU(7)/1000.d0)/(0.004189*(0.09**3.)) !BCB                    
c     SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0     !BCIA and BCI aged use 100% & 100% 
c     SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
c     SSM3=((DSU(18)+DSU(19)+DSU(20)+DSU(21))/1000.d0)/(0.004189*(0.20**3.))!SOA                   
c     SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
c     SSMO=SSMAO/(0.004189d0*(.085d0**3.d0)) + SSM2
*************************************************************************************************
C** Get Na in cm-3 as from Lohmann et al. 1999, JGR, 104,D8, 9169-9198
C** Number = [ (Mass/den) } / [ (4/3 * pi * r^3* ) }
C** effective radius for SO4 = 0.15 um; OC = 0.20 um, BC = 0.09 um, SS1 = 0.44 um; SS2 = 1.7 um
C** 4/3*pi*r^3 for each species can be precalculated as  (in cm3)
C** SO4=1.414*1.d-14; SS1 = 3.568*1.d-13, SS2=2.058*1.d-11; OC=3.351*1.d-14; BC = 3.054*1.d-15
C**  1/density of each species /(4/3 * pi * r^3) is then precalculated as
C*** SO4 = 3.981*1.d10, SS1 = 1.401*1.d9; SS2 = 2.43*1.d7; OC = 2.98*1.d10; BC = 3.27*1.d11
C** We assume densities as: SO4=1777; SS1=SS2=2000; OC=BC=1000 (kg/m3)
C** the values are copied below
C** The two columns for Vol and 1/den/vol reflect differences from size difference
C** when using eff. rad versus a std. vol radius of .052 and .085 for cont. and mari. aerosols
C**Rad (um)    	 	Volume(cm3)     Vol(cm3)	m3/cm3		m3/cm3
C* Species		4/3*pi*r^3  rl=.052,ro=.085 um	1/den/vol	1/den/vol
C**SO4: 0.15		1.41372E-14	5.89d-16	3.98E+10	9.55d11
C**SS1:0.44		3.56818E-13	2.57d-15	1.40E+09	1.94d11
C**SS2:1.7		2.05795E-11	2.06d-11	2.43E+07	2.42d07
C**OC: 0.2		3.35103E-14	5.89d-16	2.98E+10	1.70d12
C**BC:0.09		3.05363E-15	5.89d-16	3.27E+11	1.70d12
c
c     SSM1 = 9.55d11*DSU(1)         ! all sulfate 
c     SSM2 = 1.94d11*DSU(2)         ! SS 01.-1 um 
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um 
c     SSM4 = 1.70d12*DSU(4)         ! OCIA aged industrial OC
c     SSM5 = 1.70d12*DSU(5)*fq_aer(n_OCB) ! OCB
c     SSM6 = 1.70d12*DSU(6)         ! BCIA aged industrial BC
c     SSM7 = 1.70d12*DSU(7)*fq_aer(n_BCB) ! BCB
c     SSM8 = 1.70d12*(DSU(18)*fq_aer(n_isopp1a)
c    &               +DSU(19)*fq_aer(n_isopp2a)
c    &               +DSU(20)*fq_aer(n_apinp1a)
c    &               +DSU(21)*fq_aer(n_apinp2a)
c    &               )! SOA
c
C** Land Na (cm-3) is from sulfate+OC+BC and is SSMAL
C** Ocean Na (cm-3) is from sulfate+OC+BC+Seasalt and is SSMAO
c
c     SSMAL=SSM1+SSM4+SSM5+SSM6+SSM7+SSM8     !SO4,OCIA,BCIA,BCB,OCB,SOA
c     SSMAO=SSMAL+SSM2                        ! Land aerosols (SSMAL) + Sea-salt
C*************Use matrix activated fraction for aerosol conc. 
      do nm=1,nmodes
C** here we do not want sea-salt contribution for land aerosols
!       Mode #     1     2     3     4     5     6     7     8     9    ! # to identify the mode in MODES1, etc. below.
!     DATA MNAME/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','SSS',
!    &           'OCC','BC1','BC2','BC3','OCS','DBC','BOC','BCS','MXX'/
!       Mode #     10    11    12    13    14    15    16    17    18   ! # t
       if( nm.eq.7 .or. nm.eq.8) ncaero(nm)=1.d-30
       SSMAL= SSMAL + (ncaero(nm)) 
c      if(ncaero(nm).gt.1.e-04)
c    * write(6,*)"incldmat",ncaero(nm),nm,SSMAL
      enddo
      do nm=1,nmodes
       SSMAO= SSMAO + (ncaero(nm))
      enddo
c     write(6,*)"CC Matrix",SSMAL,SSMAO,l  

      IF(SSMAL.le.100.d0) SSMAL=100.d0
      IF(SSMAO.le.50.d0) SSMAO=50.d0
C** We use data from Texas based on Segal et al. 2004 from Leon R.
      MCDNL1 = 174.8d0 + (1.51d0*SSMAL**0.886d0)
      MCDNO1 = -29.6d0 + (4.917d0*SSMAO**0.694d0)
c     if (MCDNL1.gt.100.) 
c    *write(6,*)"CDNC for MC Clds",MCDNL1,MCDNO1,SSMAL,SSMAO,L
c     endif
c
      RETURN
      END SUBROUTINE GET_CC_CDNC_MX 
C****************************************************************************
c
      SUBROUTINE GET_CC_CDNC(L,AIRM,DXYPJ,PL,TL,DSS,MCDNL1,MCDNO1)
!@sum specific calculation to get cloud droplet number for indirect effects
!@auth Surabi Menon 
!@Use for calculating cloud droplet number for convective clouds
!@when using mass based aerosols
      USE CLOUDS_COM
      USE TRACER_COM, only: n_apinp1a, n_apinp2a,n_BCB, n_isopp1a, 
     * n_isopp2a, n_OCB 
      USE OldTracer_mod, only: fq_aer
      USE CONSTANT,only:mb2kg,RGAS
      IMPLICIT NONE
      real*8 AIRM,EXPL,EXPO,WCDNO,WCDNL,rho
     *,MCDNL1,MCDNO1,amass,tams,smturb,DXYPJ,PL,TL
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSM6,SSM7,SSM8,
     *       SSMAL,SSMAO,SSML,SSMO
      real*8 SSMD1,SSMD2,SSMD3, SSM1a
      integer, PARAMETER :: nt=31 ! this MUST be equal to sntm
      real*8,dimension(nt)::DSS,DSU
      integer L,n

      do n = 1,nt
        DSU(n)=1.d-10
      end do
      SSMAL=0.d0
      SSMAO=0.d0

C** add in terms for AMASS from other program to get aerosol mass conc.
C** amass is airmass in kg
      amass=AIRM*mb2kg*DXYPJ
C** This is air density in kg/m3
      rho=1d2*PL/(RGAS*TL)
C*** DSU gives you aerosol mass in  kg/m3, DSS is in kg of species
C*** DSS/amass is mass mixing ratio of aerosol (kg/kg)
      tams=1.d0/amass*rho
      do n = 1,nt
        DSU(n) =DSS(n)*tams
C** Special case if not including dust-sulfate hetchem reactions
C** but including nitrates
        if (n.gt.9.and.n.le.13) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif 
        if (n.gt.14.and.n.le.17) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif
      enddo
      SSM1 = 9.55d11*DSU(1)         ! all sulfate
C**nitrate,  vol radius = 0.3 um and effe rad = 0.15 um
C** Choose value that works as for sulfates
      SSM1a= (DSU(14)/1700.d0)/(0.004189d0*(0.058**3))   
c     SSM2 = 1.94d11*DSU(2)         ! SS 0.01-1 um
      SSM2 = 1.89d10*DSU(2)         ! SS 0.01-1 um
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um
#ifdef TRACERS_AEROSOLS_VBS
      SSM4 = 1.70d12*sum(DSU(23:31))
      SSM5 = 0.d0
#else
      SSM4 = 1.70d12*DSU(4)         ! OCIA aged industrial OC
      SSM5 = 1.70d12*DSU(5)*fq_aer(n_OCB) ! OCB
#endif
      SSM6 = 1.70d12*DSU(6)         ! BCIA aged industrial BC
      SSM7 = 1.70d12*DSU(7)*fq_aer(n_BCB) ! BCB
#ifdef TRACERS_AEROSOLS_SOA
      SSM8 = 1.70d12*(DSU(18)*fq_aer(n_isopp1a)
     &               +DSU(19)*fq_aer(n_isopp2a)
#ifdef TRACERS_TERP
     &               +DSU(20)*fq_aer(n_apinp1a)
     &               +DSU(21)*fq_aer(n_apinp2a)
#endif
     &               )! SOA
#else
      SSM8=0.d0
#endif
c
C** Land Na (cm-3) is from sulfate+OC+BC + NO3 and is SSMAL
C** Ocean Na (cm-3) is from sulfate+OC+BC+Seasalt and is SSMAO
c
      SSMAL=SSM1+SSM4+SSM5+SSM6+SSM7+SSM8+SSM1a       !SO4,OCIA,BCIA,BCB,OCB,SOA,NO3
      SSMAO=SSMAL+SSM2                           ! Land aerosols (SSMAL) + Sea-salt
c
      IF(SSMAL.le.100.d0) SSMAL=100.d0
      IF(SSMAO.le.50.d0) SSMAO=50.d0
c     write(6,*)"AEROSOL MASS",DSU(1),DSU(2),DSU(4),DSU(5),DSU(6),DSU(7)
C** We use data from Texas based on Segal et al. 2004 from Leon R.
      MCDNL1 = 174.8d0 + (1.51d0*SSMAL**0.886d0)
      MCDNO1 = -29.6d0 + (4.917d0*SSMAO**0.694d0)
c     write(6,*)"CDNC for MC Clds",MCDNL1,MCDNO1,SSMAL,SSMAO,L
c
      RETURN
      END SUBROUTINE GET_CC_CDNC


C************************************************************************************
C** For large-scale stratus clouds
C*************************************************************************************
      SUBROUTINE GET_CDNC(L,LHX,WCONST,WMUI,AIRM,WMX,DXYPJ,
     *FCLD,CAREA,CLDSAVL,DSS,PL,TL,OLDCDL,
     *VVEL,SME,DSU,CDNL0,CDNL1)
!@sum specific calculation to get cloud droplet number for indirect effects
!@auth Surabi Menon 
!@Use for calculating cloud droplet for large-scale stratus clouds
!@when using mass based aerosols
!@input is mostly aerosol mass and a few cloud properties 
      USE CLOUDS_COM
      USE TRACER_COM, only: n_apinp1a, n_apinp2a,n_BCB, n_isopp1a, 
     * n_isopp2a, n_OCB 
      USE OldTracer_mod, only: fq_aer
      USE CONSTANT,only:mb2kg,LHE,LHS,RGAS
      IMPLICIT NONE
      real*8 CAREA,CLDSAVL,AIRM,WMX,OLDCDL,VVEL  ! VVEL is in cm/s
     *,SME,rho,PL,TL,WTURB
      integer, PARAMETER :: nt=31 ! this MUST be equal to sntm
      real*8,dimension(nt)::DSS,DSU
      real*8 EXPL,EXPO,WCDNL,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,amass,tams
     *,FCLD,WCONST,LHX,WMUI,DXYPJ
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSM6,SSM7,SSM8,SSMAL,SSML
      real*8 SSMD1,SSMD2,SSMD3,SSM1a
      real*8 term1,term2,vterm,alf
      integer L,n

      do n = 1,nt
        DSU(n)=1.d-10                       
      end do
      SSMAL=0.d0

!add in terms for AMASS from other program to get aerosol mass conc.
      amass=AIRM*mb2kg*DXYPJ  
C** This is air density in kg/m3
      rho=1d2*PL/(RGAS*TL)
C*** DSU gives you aerosol mass in  kg/m3, DSS is in kg of species
C*** DSS/amass is mass mixing ratio of aerosol (kg/kg)
      tams=1.d0/amass*rho
      do n = 1,nt
        DSU(n) =DSS(n)*tams
C** Special case if not including dust-sulfate hetchem reactions but with NO3
        if (n.gt.9.and.n.le.13) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif 
        if (n.gt.14.and.n.le.17) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif 
      enddo

*************************************************************************************************
c     SSM1 = 3.98d10*DSU(1)         ! all sulfate 
c     SSM2 = 1.40d09*DSU(2)         ! SS 01.-1 um 
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um 
c     SSM4 = 2.98d10*DSU(4)         ! OCIA aged industrial OC
c     SSM5 = 2.98d10*DSU(5)         ! OCB with 80% solubility                   
c     SSM6 = 3.27d11*DSU(6)         ! BCIA aged industrial BC
c     SSM7 = 3.27d11*DSU(7)         ! BCB with 80% solubility 
c     SSM8 = 2.98d11*(DSU(18)+DSU(19)+DSU(20)+DSU(21))! SOA with 100% solubility 

      SSM1 = 9.55d11*DSU(1)         ! all sulfate 
      SSM1a= (DSU(14)/1700.d0)/(0.004189d0*(0.058**3)) !nitrate
c     SSM2 = 1.94d11*DSU(2)         ! SS 01.-1 um 
      SSM2 = 1.89d10*DSU(2)         ! SS 01.-1 um 
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um 
#ifdef TRACERS_AEROSOLS_VBS
      SSM4 = 1.70d12*sum(DSU(23:31))
      SSM5 = 0.d0
#else
      SSM4 = 1.70d12*DSU(4)         ! OCIA aged industrial OC
      SSM5 = 1.70d12*DSU(5)*fq_aer(n_OCB) ! OCB
#endif
      SSM6 = 1.70d12*DSU(6)         ! BCIA aged industrial BC
      SSM7 = 1.70d12*DSU(7)*fq_aer(n_BCB) ! BCB
#ifdef TRACERS_AEROSOLS_SOA
      SSM8 = 1.70d12*(DSU(18)*fq_aer(n_isopp1a)
     &               +DSU(19)*fq_aer(n_isopp2a)
#ifdef TRACERS_TERP
     &               +DSU(20)*fq_aer(n_apinp1a)
     &               +DSU(21)*fq_aer(n_apinp2a)
#endif
     &               )! SOA with 100% solubility
#else
      SSM8=0.d0
#endif
c
C** Na (cm-3) is from sulfate+OC+BC+sea-salt 
      SSMAL=SSM1+SSM2+SSM4+SSM5+SSM6+SSM7+SSM8+SSM1a     !SO4,OCIA,BCIA,BCB,OCB,SS1,SOA,NO3
c     write(6,*)"SSMAL",SSMAL,SSM1,SSM4,DSU(1),DSU(4),l
      IF(SSMAL.le.100.d0) SSMAL=100.d0

C** Here we use Gultepe's paramet for CDNC = f(Na)
c     EXPL=(298.d0*log10(SSMAL))-595.d0
c     EXPO=(162.d0*log10(SSMAL))-273.d0

c     IF (EXPO.LT.10.d0) EXPO=10.d0
c     IF (EXPL.LT.10.d0) EXPL=10.d0
C**Note smalphaf is to mimic turbulence. Need to create a dependency if using Gultepe's scheme
C** use CTEI effect for CDNC as in Menon et al. 2002 JAS, but replace CTEI with WTURB
C**      smalphaf=(1.d0+ 2.d0*SMFPML)*0.5d0 !old way when we did not have ATURB
c     smalphaf = (1.d0 +2.d0*SME)*0.5d0
c     WCDNO= EXPO*smalphaf
c     WCDNL= EXPL*smalphaf
C**** use Lohmann et al. 2007, ACPD,7,371-3761 formualtion
C**Convert velocity term to cm/s
      alf=0.023   !cm^4/s
      vterm = 1.d2*( .01*VVEL + (1.33*(sqrt(SME))) )    ! SME is EGCM
C* WTURB = 0.817*sqrt(egcm). To equate to 1.33*sqrt(egcm)
c     vterm = 1.d2*( VVEL + (WTURB*1.63) )          
      if(vterm.gt.0.)  then 
       term1=SSMAL*vterm
       term2=(alf*SSMAL) + vterm
       WCDNL= ((term1/term2)**1.27)*1.d-1
c      if(WTURB.gt.0.)write(6,*)"Lohmann",vterm,WCDNL,term2,term1,SSMAL
c    *,VVEL,WTURB,"Aerosol no",SSM1,SSM2,SSM4,SSM5,SSM6,SSM7,SSM8
      else
       WCDNL=10.d0   
      endif
      CDNL0=OLDCDL   !term initialised to 10. CLOUDS_DRV                  
C** Using the new CDNC scheme where we calculate it as a function of
C** gas phase sulfate and cloud area changes
      CCLD0 = CLDSAVL   !cld frac from previous time step   was 3D
      CCLD1 = 1.d0-CAREA  !CLDSS(L)cld frac from present time step 
      DCLD = CCLD1-CCLD0     ! cloud fraction change
      dfn = 0.0d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0  
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
      endif
      IF (CDNL1.le.10.d0) CDNL1=10.d0

      RETURN
      
      END SUBROUTINE GET_CDNC 

C**************************************************************************
C** Here we calculate the autoconversion rate 
C**************************************************************************
      SUBROUTINE GET_QAUT(L,PL,TL,FCLD,WMX,SCDNCW,RCLD,RHOW,r6,r6c,
     *QCRIT,QAUT)
!@sum specific calculation to obtain autoconversion rates for the indirect effects
!@auth Surabi Menon 
!@contains various routines that may be used to get autoconversion that depends on cloud droplet number or size
!@when using mass based aerosols
      USE TRACER_COM, only: n_apinp1a, n_apinp2a,n_BCB, n_isopp1a, 
     * n_isopp2a, n_OCB 
      USE OldTracer_mod, only: fq_aer
      USE CONSTANT,only:TWOPI,GRAV,by6,by3,RGAS
      IMPLICIT NONE
      real*8 TL,WMX,SCDNCW,QAUT,RHOW,FCLD,RCLD,rho,r3c
      real*8 QCRIT,rcr,dynvis,PL,qcr,QAU,GAMA2
      real*8 epsi,epsis,betu6,betd6,bet6,cwc,r6c,r6,QAUT1,QAUT2
      real*8, PARAMETER :: bcon=1.15d23,  ! s-1 condensation rate const.
     &                     kap2=1.9d11    ! cm-3 s-1
      integer L

C**sm Tripoli and Cotton (1980), adapted from Jones Tech rpt.
C     using Qcrit for QAUT
C** This is air density in kg/m3 
       rho=1d2*PL/(RGAS*TL)
c     dynvis=2.4784d0*1.d-6*TL/(TL+393.16d0)
c     rcr=7.0d0
c     QCRIT=(2.d0/3.d0)*TWOPI*1.d9*((rcr*1.d-06)**3.d0)*SCDNCW/
c    *(.001d0*rho)
c
c     QAUT=(0.104d0*GRAV*0.55d0*((.001d0*rho)**(4.d0/3.d0))*
c    * ((WMX/(FCLD+1.d-30))**(7.d0/3.d0)))/
c    * (dynvis*((SCDNCW*1.d09)**(1.d0/3.d0)))

c     if(QAUT.le.0.) write(6,*)"QCR",QAUT,SCDNCW,WMX,l
C***Using Berrys 1967 scheme
c      GAMA2 = 0.35
c      QAUT=(GAMA2*rho*(WMX/(FCLD+1.d-20))**2)/(120.d0+((5.7d0
c    *    *SCDNCW*(FCLD+1.d-30))/(rho*WMX)))

C*** Use R6AUTO from Rotstayn and Liu (2005) GRL...
C***Using dispersion parameter for alhpa of 0.003 (liu and Daum, 2002, Nature)
       epsi=1.d0-0.7d0*exp(-0.003d0*SCDNCW)
       epsis=epsi**2
       betu6=(1.d0+3.d0*epsis)*(1.d0+4.d0*epsis)*(1.d0+5.d0*epsis)
       betd6=(1.d0+epsis)*(1.d0+2.d0*epsis)
       bet6 = (betu6/betd6)
       cwc = WMX/(FCLD+1.d-30)*1.d3*rho    ! gives in-cld water content in g/m3
       r6=RCLD*bet6   ! check to make sure it is in um
       r6c=4.09d-4*(bcon**by6)*(SCDNCW**by6)/(cwc**by3)
       r3c=r6c/(bet6**by6)*1.d-6   ! gives r in meters 
       QCRIT = (2.d0*TWOPI/3.d0)*rhow*r3c**3*SCDNCW*1.d6/rho   !Qcrit unitless: r,N in m
       QAUT1=(3.d0/(2.d0*TWOPI*rhow))**2
       QAUT2=kap2*bet6*(cwc**3)/SCDNCW
C** Units of kg m^-3 s^-1
       QAU=QAUT1*QAUT2*1.d-9
       QAUT=QAU/rho     ! QAUT is in s-1
C*** In main CLOUDS2.f if qc is gt. qcri start Qaut otherwise it is 0        
c     if (RCLD.gt.5.d0)
c     write(6,*)"QAUT",r6c,r6,RCLD,cwc,WMX,FCLD,SCDNCW,QCRIT,QAUT
      RETURN
    
      END SUBROUTINE GET_QAUT

C**************************************************************************
C** Here we calculate the updated CDNC if aerosol mass has changed
C**************************************************************************
      SUBROUTINE GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX,FCLD,
     *CLDSSL,CLDSAVL,VVEL,SME,DSU,OLDCDL,CDNL0,CDNL1)
!@sum specific calculation to get cloud droplet number for indirect effects
!@auth Surabi Menon 
!@ Use for calculating cloud droplet number for large-scale stratus clouds
!@when using mass-based aerosols that are updated after various cloud processes in CLOUDS2_E1
      USE CLOUDS_COM
      USE TRACER_COM, only: n_apinp1a, n_apinp2a,n_BCB, n_isopp1a, 
     * n_isopp2a, n_OCB 
      USE OldTracer_mod, only: fq_aer
      USE CONSTANT,only:LHE,LHS
      IMPLICIT NONE
      real*8 ::CLDSSL,CLDSAVL,WMX
     *,OLDCDL,VVEL,SME,WTURB  ! VVEL is in cm/s
      real*8 EXPL,EXPO,WCDNL,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,FCLD
     *,LHX,WMUI,WCONST
      real*8 term1,term2,vterm,alf
      integer, PARAMETER :: nt=31 ! this MUST be equal to sntm
      real*8,dimension(nt)::DSU

      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSM6,SSM7,SSM8,SSMAL,SSML
      real*8 SSMD1,SSMD2,SSMD3,SSM1a

      integer L

      SSMAL=0.d0

c     SSM1 = 3.98d10*DSU(1)         ! all sulfate 
c     SSM2 = 1.40d09*DSU(2)         ! SS 01.-1 um 
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um 
c     SSM4 = 2.98d10*DSU(4)         ! OCIA aged industrial OC
c     SSM5 = 2.98d10*DSU(5)         ! OCB with 80% solubility                   
c     SSM6 = 3.27d11*DSU(6)         ! BCIA aged industrial BC
c     SSM7 = 3.27d11*DSU(7)         ! BCB with 80% solubility 
c     SSM8 = 2.98d10*(DSU(18)+DSU(19)+DSU(20)+DSU(21))! SOA with 100% solubility                   

      SSM1 = 9.55d11*DSU(1)         ! all sulfate 
      SSM1a=(DSU(14)/1700.d0)/(0.004189d0*(0.058**3))   !nitrate
c     SSM2 = 1.94d11*DSU(2)         ! SS 01.-1 um 
      SSM2 = 1.89d10*DSU(2)         ! SS 01.-1 um 
c     SSM3 = 2.43d07*DSU(3)         ! SS in 1-4 um 
#ifdef TRACERS_AEROSOLS_VBS
      SSM4 = 1.70d12*sum(DSU(23:31))
      SSM5 = 0.d0
#else
      SSM4 = 1.70d12*DSU(4)         ! OCIA aged industrial OC
      SSM5 = 1.70d12*DSU(5)*fq_aer(n_OCB) ! OCB
#endif
      SSM6 = 1.70d12*DSU(6)         ! BCIA aged industrial BC
      SSM7 = 1.70d12*DSU(7)*fq_aer(n_BCB) ! BCB
#ifdef TRACERS_AEROSOLS_SOA
      SSM8 = 1.70d12*(DSU(18)*fq_aer(n_isopp1a)
     &               +DSU(19)*fq_aer(n_isopp2a)
#ifdef TRACERS_TERP
     &               +DSU(20)*fq_aer(n_apinp1a)
     &               +DSU(21)*fq_aer(n_apinp2a)
#endif
     &               )! SOA
#else
      SSM8=0.d0
#endif
c
C** Land Na (cm-3) is from sulfate+OC+BC+seasalt 
      SSMAL=SSM1+SSM2+SSM4+SSM5+SSM6+SSM7+SSM8+SSM1a     !SO4,OCIA,BCIA,BCB,OCB,SS1,SOA,NO3
      IF(SSMAL.le.100.d0) SSMAL=100.d0

C**** use Lohmann et al. 2007, ACPD,7,371-3761 formualtion
C**Convert velocity term to cm/s
      alf=0.023   !cm^4/s
      vterm = 1.d2*( .01*VVEL + (1.33*(sqrt(SME))) )    ! SME is EGCM
C* WTURB = 0.817*sqrt(egcm). To equate to 1.33*sqrt(egcm)
c     vterm = 1.d2*( VVEL + (WTURB*1.63) )          
      if(vterm.gt.0.)  then 
       term1=SSMAL*vterm
       term2=(alf*SSMAL) + vterm
       WCDNL= ((term1/term2)**1.27)*1.d-1
      else
       WCDNL=10.d0
      endif
      CDNL0=OLDCDL

C** Using the new CDNC scheme where we calculate it as a function of
C** aerosol and cloud area changes
      CCLD0 = CLDSAVL   ! cld frac from previous time step  was 3D
      CCLD1 = CLDSSL    !cld frac from present time step  was 3D
      DCLD = CCLD1-CCLD0   ! cloud fraction change
      dfn = 0.d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
      endif
      IF (CDNL1.le.10.d0) CDNL1=10.d0

      RETURN

      END SUBROUTINE GET_CDNC_UPD
