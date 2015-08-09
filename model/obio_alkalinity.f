#include "rundeck_opts.h"

      subroutine obio_alkalinity(kmax,i,j)

!@sum  online computation of alkalinity
!@auth Natassa Romanou

! compute source/sink term for alkalinity
! zc     compensation depth: photosynthesis = respiration
! or more precisely p1d(kzc)
! J_PO4  source/sinks of phosphate
! J_Ca   source/sinks of calcium carbonate
! Jprod  production of organic phosphorus which is linked to particulate
!        organic carbon POC production by 
!        J_POC=(1-sigma_Ca)*J_DOC=(1-sigma_Ca)*rC:P*Jprod
!        to estimate J_DOC use total productivity pp
!        Really, production of particulate organic carbon (or phosphorus)
!        above the compensation depth
! F_Ca   downward flux of CaCO3
! Fc     instantaneous downward flux of particulate organic phosphorus at 
!        compensation depth
! sigma_Ca fraction of P converted to DOP, or fraction of C converted to DOC
! 1-sigma_Ca fraction of P NOT converted to DOP, or fraction of C NOT converted to DOC
! alk tendency is also computed under ice.  this is because 
! alk changes due to sal and temp changes.
! alkalinity units should be umol/kg (uE/kg to go into obio_carbon)
! alk tendency computed here is in umol/m3/s, convert to umol/kg
! *** a note on units: I assume here that units are umol,N/m3
!-------------------------------------------------------------------------
! J_ALK = -rN:P* J_PO4 + 2*J_Ca
! J_PO4 = J_NO3/rN:P, tendency of phosphates
! J_NO3 : tendency of nitrates, P_tend(:,1)
! J_Ca = R* rC:P* (1-sigma_Ca)* Jprod, for  z<zc
! J_Ca = - dF_Ca/dz,   for z>zc
! Jprod = pp, for z<zc 
! Jprod = 0,  for z>zc
! Distinguish two cases: pp based on all species
!                        pp based on coccolithophores only
! F_Ca = R * rC:P * Fc * exp(-(z-zc)/d)
!   Fc = (1-sigma_Ca) * integral_0_zc (J_DOC*dz)
!
! note sign errors in equations 28a and 31 of notes.
!-------------------------------------------------------------------------

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,
     .      npratio,uMtomgm3,cnratio,bn,zc
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend,
     .      rhs,alk1d,caexp,kzc
#ifdef OBIO_RUNOFF
#ifdef ALK_RUNOFF
      use obio_com, only: ralkconc_loc
      USE OFLUXES, only:  oFLOWO
      USE MODEL_COM, only: dtsrc
#endif
#endif


#ifdef OBIO_ON_GARYocean
      USE OCEAN, only: dxypo
#else
      USE hycom_arrays, only: scp2
#endif

      implicit none

      integer nt,k,kmax,nchl1,nchl2,i,j
      real*8 J_PO4(kmax),pp,Jprod(kmax),Jprod_sum,Fc,zz,F_Ca(kmax+1),
     .       J_Ca(kmax),term,term1,term2,DOP
!--------------------------------------------------------------------------
!only compute tendency terms if total depth greater than conpensation depth
      if (p1d(kmax+1) .lt. p1d(kzc)) then
         A_tend(:) = 0.
         go to 100
      endif

!pp and therefore Jprod are defined at dp points (mid level)
!F_Ca is defined at pressure interfaces
!J_Ca is therefore defined mid-layer (dp points)

!compute sources/sinks of phosphate
!J_PO4 units uM/hr
      do k=1,kmax
      J_PO4(k) =  P_tend(k,1)/npratio   !approximate by nitrate conc tendency
                                        !NO3/PO4 ratio from Conkright et al, 1994
     .            /14.d0                !expressed as nitrates
      term = -1.d0*npratio * J_PO4(k)   !uM,N/hr= mili-mol,N/m3/hr
      rhs(k,15,1) = term
      A_tend(k)= term 
      enddo

!distinguish two cases: OCMIP uses total pp for Jprod, whereas here we also
!consider the case where Jprod only includes coccolithophores
#ifdef Jprod_based_on_pp
      nchl1=1
      nchl2=nchl
#endif
#ifdef Jprod_based_on_cocc
      nchl1=4
      nchl2=4   
#endif

      if (nchl1.le.0 .or. nchl2.le.0) then
          print*, nchl1, nchl2, 
     .   'MUST SET Jprod_based_on_pp or Jprod_based_on_cocc in rundeck'
          stop  
      endif

!compute total primary production
!and the Jprod term (production of organic phosphorus)
      do k=1,kmax
        if (p1d(k) .le. p1d(kzc))  then
          pp=0.
          do nt=nchl1,nchl2
             pp=pp+pp2_1d(k,nt)/dp1d(k)   !mgC/m3/hr
          enddo
          Jprod(k) = pp/cpratio  !mgP/m3/hr
        else
          Jprod(k)=0.
        endif
      enddo

!integrate net primary production down to zc
      Jprod_sum = 0.
      do k=1,kmax
      if (p1d(k) .le. p1d(kzc))  then
       Jprod_sum = Jprod_sum + Jprod(k)*dp1d(k)   ! mgP/m2/hr
      endif
      enddo

      Fc = (1.d0-sigma_Ca)* Jprod_sum     !mgP/m2/hr

!compute downward flux of CaCO3
!only below the euphotic zone (the compensation layer)
      F_Ca = 0.d0
      do k=kzc,kmax+1
           F_Ca(k) = rain_ratio*cpratio*Fc
     .             * exp(-1.d0*(p1d(k)-p1d(kzc))/d_Ca)    !mgC/m2/hr
      enddo

      !p1d(kzc) is really the compensation depth
      !F_Ca(kzc) is the CaCO3 export
      caexp = 0.d0
      do k=1,kzc
      caexp = caexp + F_Ca(k)
     .                *24.d0*365.d0
     .                *1.d-15       !PgC/m2/yr
#ifdef OBIO_ON_GARYocean
     .                * dxypo(j)    ! -> Pg,C/yr
#else
     .                * scp2(i,j)   ! -> Pg,C/yr
#endif
      enddo

!compute sources/sinks of CaCO3
      do k=1,kmax
         if (p1d(k) .le. p1d(kzc))  then
             !formation of calcium carbonate above compensation depth
             J_Ca(k) = -1.d0*rain_ratio*cpratio*(1.-sigma_Ca)*Jprod(k)   !mgC/m3/hr
         else
             !dissolution of calcium carbonate below compensation depth
             J_Ca(k) = -1.d0* (F_Ca(k+1)-F_Ca(k)) / dp1d(k)   !mgC/m3/hr
         endif
       term = 2.d0* J_Ca(k)/cnratio  ! mgC/m3/hr -> mili-mol,N/m3/hr  
                                     ! no need to multiply here by mol.weight
                                     ! because already in cnratio (see obio_init)
       rhs(k,15,5) = term
       A_tend(k) = A_tend(k) + term      
      enddo

#ifdef OBIO_RUNOFF
#ifdef ALK_RUNOFF
      term = ralkconc_loc(i,j)
     . * oFLOWO(i,j)/dtsrc           ! mol/kg,w => mol/m2,w/s
     . / dp1d(1)                     ! mol/m2,w/s => mol/m3,w/s
     . * 1.d3                        ! mol/m3,w/s => mmol/m3,w/s
     . * 3600.                       ! mmol/m3,w/s => mmol/m3/hr
      rhs(1,15,17) = term
      A_tend(1) = A_tend(1) + term
#endif
#endif 

      !for consistency, keep term that goes into rhs table in uM/hr = mili-mol,N/m3/hr
      !convert A_tend terms into uE/kg/hr, the actual units of alkalinity 
      A_tend = A_tend /1024.5d0 *1.d3     ! mili-mol,N/m3/hr -> umol/m3/hr -> umol/kg/hr

!!!!!!!!!! NEED TO ADD BOTTOM BOUNDARY CONDITIONS 

 100  continue
      end subroutine obio_alkalinity
