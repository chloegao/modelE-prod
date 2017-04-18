#include "rundeck_opts.h"


!@sum multicoag  :  performs coagulation on the aerosol size distribution
!@+  defined by Nk and Mk (number and mass).  See "An Efficient
!@+   Numerical Solution to the Stochastic Collection Equation", S.
!@+   Tzivion, G. Feingold, and Z. Levin, J Atmos Sci, 44, no 21, 3139-
!@+   3149, 1987.  Unless otherwise noted, all equation references refer
!@+   to this paper.  Some equations are taken from "Atmospheric Chemistry
!@+   and Physics: From Air Pollution to Climate Change" by Seinfeld
!@+   and Pandis (S&P).  

!@+   This routine uses a "moving sectional" approach in which the
!@+   aerosol size bins are defined in terms of dry aerosol mass.
!@+   Addition or loss of water, therefore, does not affect which bin
!@+   a particle falls into.  As a result, this routine does not
!@+   change Mk(water), although water masses are needed to compute
!@+   particle sizes and, therefore, coagulation coefficients.  Aerosol
!@+   water masses in each size bin will need to be updated later
!@+   (in another routine) to reflect changes that result from
!@+   coagulation.
!@+   The user must supply the mass and number distributions, Mk and Nk,
!@+   as well as the time step, dt.

!@auth Yunha Lee

#if (defined TOMAS_12_10NM) || (defined TOMAS_12_3NM)

      SUBROUTINE multicoag(dt)


      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE CONSTANT,   only:  pi,gasc   
      IMPLICIT NONE

      real dt !time step [s]
      integer n,c,bh,bl,dts_old,count ! for 15 size bins lumping
      integer k,j,i,jj,kk    !counters
      real*8 dNdt(ibins), dMdt(ibins,icomp-idiag)
!@var   dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
!@+   the mass boundaries of the size bins.  xbar is the average mass
!@+   of a given size bin (it varies with time in this algorithm).  phi
!@+   and eff are defined in the reference, equations 13a and b.

      real*8 xbar(ibins), phi(ibins), eff(ibins)
      real*8 Nki(ibins), Mki(ibins, icomp) 
!@var kij represents the coagulation coefficient (cm3/s) normalized by the
!@+   volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
      real kij(ibins,ibins)
      real Dpk(ibins)             !diameter (m) of particles in bin k
      real Dk(ibins)              !Diffusivity (m2/s) of bin k particles
      real ck(ibins)              !Mean velocity (m/2) of bin k particles
      real olddiff                !used to iterate to find diffusivity
      real density                !density (kg/m3) of particles
      real mu                     !viscosity of air (kg/m s)
      real mfp                    !mean free path of air molecule (m)
      real Kn                     !Knudsen number of particle
      real*8 mp         !particle mass (kg)
      real beta                   !correction for coagulation coeff.
      real*8 aerodens
      external aerodens
      real*8 Mktot      !total mass of aerosol

      !temporary summation variables
      real*8 k1m(icomp-idiag),k1mx(icomp-idiag)
      real*8 k1mx2(icomp-idiag)
      real*8 k1mtot,k1mxtot
      real*8 sk2mtot, sk2mxtot
      real*8 sk2m(icomp-idiag), sk2mx(icomp-idiag)
      real*8 sk2mx2(icomp-idiag)
      real*8 High_in,dNtot,dMtot
      real*8 mtotal

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl   

      real zeta                      !see reference, eqn 6
      real tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
      real dts  !internal time step (<dt for stability)
      real tsum !time so far
      real*8 Neps !minimum value for Nk
      real*8 mi, mf   !initial and final masses
      logical is_nan
      external is_nan
      parameter(zeta=1.28125, dtlimit=0.25, itlimit=10)
      real kB  !kB is Boltzmann constant (J/K)
      parameter (kB=1.38e-23, Neps=1.0e-3)

C-----CODE--------------------------------------------------------------

      tsum = 0.0
      dts_old=10.
C If any Nk are zero, then set them to a small value to avoid division by zero
      do k=1,ibins
         if (Nk(k) .lt. Neps) then
            Nk(k)=Neps
            Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
            do j=1,icomp
               if (j.ne.srtso4)then
                  Mk(k,j)=0.d0
               endif
            enddo
         endif
cyhl To check whether mass is conserved during coagulation. 
c         print*,'Nk',k,Nk(k),Mk(k,2)
      enddo

      Nki(:)=Nk(:)
      Mki(:,:) =Mk(:,:)

C Calculate air viscosity and mean free path

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6

      do k=1,ibins

         mso4=Mk(k,srtso4) 
         mnacl=Mk(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=Mk(k,srtnh4) !0.1875*mso4  !assume ammonium bisulfate
         mecob=Mk(k,srtecob)
         mecil=Mk(k,srtecil)
         mocil=Mk(k,srtocil)
         mocob=Mk(k,srtocob)
         mdust=Mk(k,srtdust)          
         mh2o=Mk(k,srth2o)   

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 

         Mktot=0.d0
         do j=1,icomp
            Mktot=Mktot+Mk(k,j)
         enddo
         mp=Mktot/Nk(k)
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
          ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo

C Calculate coagulation coefficients

      do i=1,ibins
         do j=1,ibins
            Kn=4.0*(Dk(i)+Dk(j))          
     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
            kij(j,i)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
            kij(j,i)=kij(j,i)*1.0e6/boxvol  !normalize by grid cell volume
         enddo
      enddo
 10   continue     !repeat process here if multiple time steps are needed

C Calculate xbar, phi and eff

      do k=1,ibins

         xbar(k)=0.0
         do j=1,icomp-idiag
            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
         enddo

       if(k.lt.ibins-1)then !from 1 to 10 bins

         eff(k)=2./9.*Nk(k)/xk(k)
     *        *(4.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
         phi(k)=2./9.*Nk(k)/xk(k)
     *        *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./3.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./3.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
      else                      ! from 11 bins to 12 bins
            eff(k)=2./31./31.*Nk(k)/xk(k)
     *           *(32.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
            phi(k)=2./31./31.*Nk(k)/xk(k)
     *           *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./31.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./31.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
      endif

      enddo

C Necessary initializations
         sk2mtot=0.0
         sk2mxtot=0.0
         do j=1,icomp-idiag
            sk2m(j)=0.0
            sk2mx(j)=0.0
            sk2mx2(j)=0.0
         enddo

C Calculate rates of change for Nk and Mk

      do k=1,ibins

         !Initialize to zero
         do j=1,icomp-idiag
            k1m(j)=0.0
            k1mx(j)=0.0
            k1mx2(j)=0.0
         enddo
         High_in=0.0
         k1mtot=0.0
         k1mxtot=0.0

         !Calculate sums
         do j=1,icomp-idiag
            if (k .gt. 1.and.k.lt.ibins) then
               do i=1,k-1
                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)*zeta
                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2.
     *                 *zeta**3.
               enddo
            elseif(k.eq.ibins)then
                  k1m(j)= sk2m(j)+kij(k,k-1)*Mk(k-1,j)
                  k1mx(j)=sk2mx(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)*4.754
                  k1mx2(j)=sk2mx2(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)**2.
     *                 *107.4365
            endif
            k1mtot=k1mtot+k1m(j)
            k1mxtot=k1mxtot+k1mx(j)
         enddo
         if(k.lt.ibins)then
            do i=k+1, ibins
               High_in=High_in+Nk(i)*kij(k,i)
            enddo
         endif

         if(k.lt.ibins-1)then

         dNdt(k)= -Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.125
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/6./xk(k)*k1mxtot)
     *    -kij(k,k)*(phi(k)/3.*xbar(k)*Nk(k)+(eff(k)-phi(k))/18.
     *        /xk(k)*zeta*xbar(k)*xbar(k)*Nk(k))

         if (k .gt. 1) then      
cyhl Nk*low_in changes to -0.5*Kij*Nk**2.
         dNdt(k)=dNdt(k)+0.625*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *        *sk2mxtot)
     * +kij(k-1,k-1)*(phi(k-1)/3.*xbar(k-1)*Nk(k-1)+(eff(k-1)
     *  -phi(k-1))/18./xk(k-1)*zeta*xbar(k-1)*xbar(k-1)*Nk(k-1))


cyhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
cyhl It results in much closer result as 30 bins. Apr.27.08
 
      call nanstop(dNdt(k),246,k,0)

         endif

         do j=1,icomp-idiag
                
            dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     *  -(phi(k)*xk(k+1)*k1m(j)+(eff(k)+2.*phi(k))/6.*k1mx(j)
     *           +(eff(k)-phi(k))/6./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)/3. ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*(phi(k)*xk(k+1)*Mk(k,j)/3.
     *           +(eff(k)+2.*phi(k))/6.*zeta*xbar(k)*Mk(k,j)/3.
     * +(eff(k)-phi(k))/6./xk(k)*zeta**3.*xbar(k)**2.*Mk(k,j)/3.)
      
cyhl  Term9(-kij(k,k)*Nk(k)*Mk(k,j)) is cancled out by term6 (k)  
          
            if (k .gt. 1) then
               dMdt(k,j)=dMdt(k,j)
     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)+2.*phi(k-1))/6.*sk2mx(j)
     *        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *  +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *    +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *    *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *    /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)

               call nanstop(dMdt(k,J),265,k,J)    
            endif
         enddo
      else if (k.eq.ibins-1)then

         dNdt(k)=0.625*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *        *sk2mxtot)
     * +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)/3.*(phi(k-1)+(eff(k-1)
     *  -phi(k-1))/6./xk(k-1)*zeta*xbar(k-1))

cyhl updated the following
         dNdt(k)=dNdt(k)-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.02
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *    -kij(k,k)*xbar(k)*Nk(k)*0.484*(phi(k)+(eff(k)-phi(k))/62.
     *        /xk(k)*4.754*xbar(k))

cyhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
cyhl It results in much closer result as 30 bins. Apr.27.08

         do j=1,icomp-idiag
            dMdt(k,j)=
     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)+2.*phi(k-1))/6.*sk2mx(j)
     *        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *  +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *    +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *    *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *    /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)


cyhl updated the following
            dMdt(k,j)= dMdt(k,j)+Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     * -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k))
     *    *k1mx(j)+(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *   +0.484*phi(k))*4.754*xbar(k)
     * +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
         enddo

      else if (k.eq.ibins)then
         dNdt(k)=-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.103226
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *    -kij(k,k)*0.484*xbar(k)*Nk(k)*(phi(k)+(eff(k)-phi(k))
     *       /62./xk(k)*4.754*xbar(k))
     *    +0.52*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mxtot)
     * +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)*0.484*(phi(k-1)+(eff(k-1)
     *  -phi(k-1))/62./xk(k-1)*4.754*xbar(k-1))

         do j=1,icomp-idiag
            dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     * -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k))*k1mx(j)
     *   +(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *   +0.484*phi(k))*4.754*xbar(k)
     * +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)

     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)/62.+0.484*phi(k-1))
     *    *sk2mx(j)+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)*0.103226
     *  +kij(k-1,k-1)*Mk(k-1,j)*0.484*(phi(k-1)*xk(k)
     * +(eff(k-1)/62.+0.484*phi(k-1))*4.754*xbar(k-1)
     * +(eff(k-1)-phi(k-1))/62./xk(k-1)*107.4365*xbar(k-1)**2.)
         enddo
      endif

         !Save the summations that are needed for the next size bin
         sk2mtot=k1mtot
         sk2mxtot=k1mxtot
         do j=1,icomp-idiag
            sk2m(j)=k1m(j)
            sk2mx(j)=k1mx(j)
            sk2mx2(j)=k1mx2(j)
         enddo
c         print*,'dNdt',Nk(k),dNdt(k)
      enddo  !end of main k loop


C Update Nk and Mk according to rates of change and time step

      !If any Mkj are zero, add a small amount to achieve finite
      !time steps
      do k=1,ibins
         do j=1,icomp-idiag
            if (Mk(k,j) .eq. 0.d0) then
               !add a small amount of mass
               mtotal=0.d0
               do jj=1,icomp-idiag
                  mtotal=mtotal+Mk(k,jj)
               enddo
               Mk(k,j)=1.d-10*mtotal
            endif
         enddo
      enddo

      !Choose time step
      dts=dt-tsum      !try to take entire remaining time step
cdbg      limit='comp'
      do k=1,ibins
         if (Nk(k) .gt. Neps) then
            !limit rates of change for this bin
            if (dNdt(k) .lt. 0.0) tlimit=dtlimit
            if (dNdt(k) .gt. 0.0) tlimit=itlimit
            if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then 
               dts=Nk(k)*tlimit/abs(dNdt(k))
!               if(dts.eq.0) print*,'dts Nk',k,dts,dNdt(k),Nk(k)
            endif
            do j=1,icomp-idiag
               if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
               if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
               if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then 
                  mtotal=0.d0
                  do jj=1,icomp-idiag
                     mtotal=mtotal+Mk(k,jj)
                  enddo         !only use this criteria if this species is significant

                  if ((Mk(k,j)/mtotal) .gt. 5.d-4) then
                     dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
 
                     if(dts.eq.0)print*,'dts MK',k,j,dts,dt,tsum,Mk(k,j)
     &                    ,nk(k),mtotal,dMdt(k,j),xk(k),xk(k+1)
                     if(dts.eq.0.and.dts_old.eq.0.) then
                       open (1044,file='debug_coag.dat',access='append',
     &                      status='unknown')
                       write(1044,*)'start',k,j,dts,dt,tsum,Mk(k,j)
     &                    ,nk(k),mtotal,dMdt(k,j),xk(k),xk(k+1)
                       write(1044,*)'dts MK, T,P',temp,pres
                       write(1044,*)'dts Mk, BOX',boxvol,boxmass,rh
                       do kk=1,ibins
                       write(1044,*)'dts MK, Nk=',Nk(kk),Nki(kk)
                       enddo
                       do jj=1,icomp
                         write(1044,*)'comp',jj
                       do kk=1,ibins
            write(1044,*)'dts MK, Mk=',Mk(kk,jj),Mki(kk,jj)
                       enddo
                       enddo
                     endif

                  else
                     if (dMdt(k,j) .lt. 0.0) then !set dmdt to 0 to avoid very small mk going negative
                        dMdt(k,j)=0.0
                     endif
                  endif
c                  print*,'mk',k,j,Mk(k,j)/mtotal,dMdt(k,j)
               endif
            enddo
c            print*,'dNdt',k,Nk(k),dNdt(k),Mk(k,1),dMdt(k,1) 
         else
            !nothing in this bin - don't let it affect time step
            Nk(k)=Neps
            Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
            if (dNdt(k) .lt. 0.0) dNdt(k)=0.0  !make sure mass/number don't go negative
            do j=1,icomp-idiag
               if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
            enddo
         endif
      enddo

c      if (dts .lt. 20.) write(*,*) 'dts<20. in multicoag',dts,tsum
!       if (dts .eq. 0.) then
!          write(*,*) 'time step is 0. dts=',dts
!          call stop_model('dts=0 in multicoag',255)
!       endif

!YUNHA (Sep, 2012) this is newly added to prevent an occasional crash. 
      if(dts.eq.0) then
! When dts is small and this is due to the mass/number is out of size range, 
! it calls mnfix. Doing this, it will help to avoid dts=0 case. 
        Nki(:)=Nk(:)
        Mki(:,:)=Mk(:,:)
        call mnfix(Nki,Mki)       !YUNHA LEE (08/28/2012) 
        count=0
        do k=1,ibins
          if(Nki(k).eq.Nk(k)) then
            count=count+1
          endif
        enddo

        if(count.lt.ibins) then
          print*,'dts=0 but mnfix worked',dts,tsum,count ! count<5 is random choice
          Nk(:)=Nki(:)
          Mk(:,:)=Mki(:,:)
        endif
      endif
     
      if(dts.eq.0.and.count.eq.ibins)then
!This case, next dts will be zero and model should be stopped. 
        write(*,*) 'time step is 0',count
        call stop_model('dts=0 in multicoag',255)
      endif

      do k=1,ibins
         Nk(k)=Nk(k)+dNdt(k)*dts
         do j=1,icomp-idiag
            Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
         enddo
      enddo

      if(dts.lt.1e-10.and.dts_old.lt.1e-10) then
! When dts is small in two sequently, then mnfix is called. It might help to get a larger dts next time.  However, I am not sure if it is helpful   !YUNHA LEE (09/05/2012) 
        Nki(:)=Nk(:)
        Mki(:,:)=Mk(:,:)
        call mnfix(Nki,Mki)     
        Nk(:)=Nki(:)
        Mk(:,:)=Mki(:,:)
      endif

      tsum=tsum+dts
      dts_old=dts
      if (tsum .lt. dt) goto 10

      RETURN
      END SUBROUTINE multicoag
#endif 


c$$$  This is multicoag for TOMAS-30 model
c$$$
c$$$!@+   **************************************************
c$$$!@+   *  multicoag                                     *
c$$$!@+   **************************************************
c$$$
c$$$!@auth   Peter Adams, June 1999
c$$$!@+   Modified to allow for multicomponent aerosols, February 2000
c$$$
c$$$!@sum   :  performs coagulation on the aerosol size distribution
c$$$!@+   defined by Nk and Mk (number and mass).  See "An Efficient
c$$$!@+   Numerical Solution to the Stochastic Collection Equation", S.
c$$$!@+   Tzivion, G. Feingold, and Z. Levin, J Atmos Sci, 44, no 21, 3139-
c$$$!@+   3149, 1987.  Unless otherwise noted, all equation references refer
c$$$!@+   to this paper.  Some equations are taken from "Atmospheric Chemistry
c$$$!@+   and Physics: From Air Pollution to Climate Change" by Seinfeld
c$$$!@+   and Pandis (S&P).  
c$$$
c$$$!@+   This routine uses a "moving sectional" approach in which the
c$$$!@+   aerosol size bins are defined in terms of dry aerosol mass.
c$$$!@+   Addition or loss of water, therefore, does not affect which bin
c$$$!@+   a particle falls into.  As a result, this routine does not
c$$$!@+   change Mk(water), although water masses are needed to compute
c$$$!@+   particle sizes and, therefore, coagulation coefficients.  Aerosol
c$$$!@+   water masses in each size bin will need to be updated later
c$$$!@+   (in another routine) to reflect changes that result from
c$$$!@+   coagulation.
c$$$
c$$$C-----INPUTS------------------------------------------------------------
c$$$
c$$$!@+   The user must supply the mass and number distributions, Mk and Nk,
c$$$!@+   as well as the time step, dt.
c$$$
c$$$C-----OUTPUTS-----------------------------------------------------------
c$$$
c$$$!@+   The program updates Nk and Mk.
c$$$
c$$$      SUBROUTINE multicoag(dt)
c$$$
c$$$      USE TOMAS_AEROSOL
c$$$      USE TRACER_COM, only : xk
c$$$      USE CONSTANT,   only:  pi,gasc   
c$$$      IMPLICIT NONE
c$$$
c$$$C-----ARGUMENT DECLARATIONS---------------------------------------------
c$$$
c$$$c      real dt         !time step (s)
c$$$      real dt
c$$$
c$$$C-----VARIABLE DECLARATIONS---------------------------------------------
c$$$
c$$$      integer k,j,i,jj,kk    !counters
c$$$      real*8 dNdt(ibins), dMdt(ibins,icomp-idiag)
c$$$      real*8 xbar(ibins), phi(ibins), eff(ibins)
c$$$
c$$$C kij represents the coagulation coefficient (cm3/s) normalized by the
c$$$C volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
c$$$      real kij(ibins,ibins)
c$$$      real Dpk(ibins)             !diameter (m) of particles in bin k
c$$$      real Dk(ibins)              !Diffusivity (m2/s) of bin k particles
c$$$      real ck(ibins)              !Mean velocity (m/2) of bin k particles
c$$$      real olddiff                !used to iterate to find diffusivity
c$$$      real density                !density (kg/m3) of particles
c$$$      real mu                     !viscosity of air (kg/m s)
c$$$      real mfp                    !mean free path of air molecule (m)
c$$$      real Kn                     !Knudsen number of particle
c$$$      real*8 mp         !particle mass (kg)
c$$$      real beta                   !correction for coagulation coeff.
c$$$      real aerodens
c$$$      external aerodens
c$$$      real*8 Mktot      !total mass of aerosol
c$$$
c$$$      !temporary summation variables
c$$$      real*8 k1m(icomp-idiag),k1mx(icomp-idiag)
c$$$      real*8 k1mx2(icomp-idiag)
c$$$      real*8 k1mtot,k1mxtot
c$$$      real*8 sk2mtot, sk2mxtot
c$$$      real*8 sk2m(icomp-idiag), sk2mx(icomp-idiag)
c$$$      real*8 sk2mx2(icomp-idiag)
c$$$      real*8 in
c$$$      real*8 mtotal
c$$$
c$$$      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
c$$$      real mecil,mecob,mocil,mocob
c$$$      real mdust,mnacl   
c$$$
c$$$      real zeta                      !see reference, eqn 6
c$$$      real tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
c$$$      real dts  !internal time step (<dt for stability)
c$$$      real tsum !time so far
c$$$      real*8  Neps              !minimum value for Nk
c$$$cdbg      character*12 limit        !description of what limits time step
c$$$
c$$$      real*8 mi, mf             !initial and final masses
c$$$      logical is_nan
c$$$      external is_nan
c$$$      
c$$$!@+   VARIABLE COMMENTS...
c$$$
c$$$!@+   dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
c$$$!@+   the mass boundaries of the size bins.  xbar is the average mass
c$$$!@+   of a given size bin (it varies with time in this algorithm).  phi
c$$$!@+   and eff are defined in the reference, equations 13a and b.
c$$$
c$$$C-----ADJUSTABLE PARAMETERS---------------------------------------------
c$$$
c$$$      parameter(zeta=1.0625, dtlimit=0.25, itlimit=10.)
c$$$      real kB  !kB is Boltzmann constant (J/K)
c$$$      parameter (kB=1.38e-23, Neps=1.0e-3)
c$$$
c$$$ 1    format(16E15.3)
c$$$
c$$$C-----CODE--------------------------------------------------------------
c$$$
c$$$      tsum = 0.0
c$$$
c$$$C If any Nk are zero, then set them to a small value to avoid division by zero
c$$$      do k=1,ibins
c$$$         if (Nk(k) .lt. Neps) then
c$$$            Nk(k)=Neps
c$$$            Mk(k,srtso4)=Neps*sqrt(xk(k+1)*xk(k)) !make the added particles SO4
c$$$            do j=1,icomp
c$$$               if (j.ne.srtso4)then
c$$$                  Mk(k,j)=0.d0
c$$$               endif
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$
c$$$C Calculate air viscosity and mean free path
c$$$
c$$$      mu=2.5277e-7*temp**0.75302
c$$$      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6
c$$$!      call nanstop(mfp,124, 0,0)
c$$$C Calculate particle sizes and diffusivities
c$$$      do k=1,ibins
c$$$
c$$$         mso4=Mk(k,srtso4) 
c$$$         mnacl=Mk(k,srtna)
c$$$         mno3=0.e0
c$$$         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
c$$$         mnh4=0.1875*mso4  !assume ammonium bisulfate
c$$$         mecob=Mk(k,srtecob)
c$$$         mecil=Mk(k,srtecil)
c$$$         mocil=Mk(k,srtocil)
c$$$         mocob=Mk(k,srtocob)
c$$$         mdust=Mk(k,srtdust)          
c$$$         mh2o=Mk(k,srth2o)   
c$$$
c$$$         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
c$$$     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
c$$$
c$$$c$$$         density=aerodens(Mk(k,srtso4),0.d0,Mk(k,srtnh4),
c$$$c$$$     &        Mk(k,srtna),Mk(k,srtecil),Mk(k,srtecob),
c$$$c$$$     &        Mk(k,srtocil),Mk(k,srtocob),Mk(k,srtdust),
c$$$c$$$     &        Mk(k,srth2o))    !assume bisulfate
c$$$!         call nanstop(density,130,k,0)
c$$$
c$$$Ckpc  Add 0.2x first for ammonium, and then add 1.0x in the loop
c$$$         Mktot=0.d0
c$$$         do j=1,icomp
c$$$            Mktot=Mktot+Mk(k,j)
c$$$         enddo
c$$$         mp=Mktot/Nk(k)
c$$$
c$$$         Dpk(k)=((mp/density)*(6./pi))**(0.333)
c$$$      call nanstop(dpk(k),139,k,0)
c$$$         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
c$$$
c$$$         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
c$$$     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
c$$$
c$$$         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
c$$$
c$$$      enddo
c$$$
c$$$C Calculate coagulation coefficients
c$$$
c$$$      do i=1,ibins
c$$$         do j=1,ibins
c$$$            Kn=4.0*(Dk(i)+Dk(j))          
c$$$     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
c$$$
c$$$            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
c$$$
c$$$            !This is S&P eqn 12.46 with non-continuum correction, beta
c$$$            kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
c$$$
c$$$            kij(i,j)=kij(i,j)*1.0e6/boxvol  !normalize by grid cell volume
c$$$
c$$$         enddo
c$$$      enddo
c$$$
c$$$ 10   continue     !repeat process here if multiple time steps are needed
c$$$
c$$$C Calculate xbar, phi and eff
c$$$
c$$$      do k=1,ibins
c$$$
c$$$         xbar(k)=0.0
c$$$         do j=1,icomp-idiag
c$$$            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
c$$$         enddo
c$$$
c$$$         eff(k)=2.*Nk(k)/xk(k)*(2.-xbar(k)/xk(k))    !eqn 13a
c$$$      call nanstop(eff(k),180,k,0)
c$$$         phi(k)=2.*Nk(k)/xk(k)*(xbar(k)/xk(k)-1.)    !eqn 13b
c$$$      call nanstop(phi(k),182,k,0)   
c$$$         !Constraints in equation 15
c$$$         if (xbar(k) .lt. xk(k)) then
c$$$            eff(k)=2.*Nk(k)/xk(k)
c$$$            phi(k)=0.0
c$$$         else if (xbar(k) .gt. xk(k+1)) then
c$$$            phi(k)=2.*Nk(k)/xk(k)
c$$$            eff(k)=0.0
c$$$         endif
c$$$      enddo
c$$$
c$$$C Necessary initializations
c$$$         sk2mtot=0.0
c$$$         sk2mxtot=0.0
c$$$         do j=1,icomp-idiag
c$$$            sk2m(j)=0.0
c$$$            sk2mx(j)=0.0
c$$$            sk2mx2(j)=0.0
c$$$         enddo
c$$$
c$$$C Calculate rates of change for Nk and Mk
c$$$
c$$$      do k=1,ibins
c$$$
c$$$         !Initialize to zero
c$$$         do j=1,icomp-idiag
c$$$            k1m(j)=0.0
c$$$            k1mx(j)=0.0
c$$$            k1mx2(j)=0.0
c$$$         enddo
c$$$         in=0.0
c$$$         k1mtot=0.0
c$$$         k1mxtot=0.0
c$$$
c$$$         !Calculate sums
c$$$         do j=1,icomp-idiag
c$$$            if (k .gt. 1) then
c$$$               do i=1,k-1
c$$$
c$$$                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
c$$$                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)
c$$$                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2
c$$$               enddo
c$$$            endif
c$$$            k1mtot=k1mtot+k1m(j)
c$$$            k1mxtot=k1mxtot+k1mx(j)
c$$$         enddo
c$$$         if (k .lt. ibins) then
c$$$            do i=k+1,ibins
c$$$               in=in+Nk(i)*kij(k,i)
c$$$            enddo
c$$$         endif
c$$$
c$$$         !Calculate rates of change
c$$$         dNdt(k)= 
c$$$     &           -kij(k,k)*Nk(k)**2
c$$$     &           -phi(k)*k1mtot
c$$$     &           -zeta*(eff(k)-phi(k))/(2*xk(k))*k1mxtot
c$$$     &           -Nk(k)*in
c$$$      call nanstop(dNdt(k),240,k,0)
c$$$         if (k .gt. 1) then
c$$$         dNdt(k)=dNdt(k)+
c$$$     &           0.5*kij(k-1,k-1)*Nk(k-1)**2
c$$$     &           +phi(k-1)*sk2mtot
c$$$     &           +zeta*(eff(k-1)-phi(k-1))/(2*xk(k-1))*sk2mxtot
c$$$      call nanstop(dNdt(k),246,k,0)
c$$$         endif
c$$$
c$$$         do j=1,icomp-idiag
c$$$            dMdt(k,j)= 
c$$$     &           +Nk(k)*k1m(j)
c$$$     &           -kij(k,k)*Nk(k)*Mk(k,j)
c$$$     &           -Mk(k,j)*in
c$$$     &           -phi(k)*xk(k+1)*k1m(j)
c$$$     &           -0.5*zeta*eff(k)*k1mx(j)
c$$$     &           +zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
c$$$      call nanstop(dMdt(k,j),257,k,j)
c$$$C      write(79,*) '243,k,j,dMdt',k,j,dMdt(k,j)
c$$$            if (k .gt. 1) then
c$$$               dMdt(k,j)=dMdt(k,j)+
c$$$     &           kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)
c$$$     &           +phi(k-1)*xk(k)*sk2m(j)
c$$$     &           +0.5*zeta*eff(k-1)*sk2mx(j)
c$$$     &           -zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
c$$$      call nanstop(dMdt(k,j),265,k,j)
c$$$C      write(79,*) '243,k,j,dMdt',k,j,dMdt(k,j)
c$$$            endif
c$$$cdbg            if (j. eq. srtso4) then
c$$$cdbg               if (k. gt. 1) then
c$$$cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
c$$$cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
c$$$cdbg     &               0.5*zeta*eff(k)*k1mx(j),
c$$$cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j),
c$$$cdbg     &               kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j),
c$$$cdbg     &               phi(k-1)*xk(k)*sk2m(j),
c$$$cdbg     &               0.5*zeta*eff(k-1)*sk2mx(j),
c$$$cdbg     &               zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
c$$$cdbg               else
c$$$cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
c$$$cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
c$$$cdbg     &               0.5*zeta*eff(k)*k1mx(j),
c$$$cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
c$$$cdbg               endif
c$$$cdbg            endif
c$$$         enddo
c$$$
c$$$cdbg         write(*,*) 'k,dNdt,dMdt: ', k, dNdt(k), dMdt(k,srtso4)
c$$$
c$$$         !Save the summations that are needed for the next size bin
c$$$         sk2mtot=k1mtot
c$$$         sk2mxtot=k1mxtot
c$$$         do j=1,icomp-idiag
c$$$            sk2m(j)=k1m(j)
c$$$            sk2mx(j)=k1mx(j)
c$$$            sk2mx2(j)=k1mx2(j)
c$$$         enddo
c$$$
c$$$      enddo  !end of main k loop
c$$$
c$$$C Update Nk and Mk according to rates of change and time step
c$$$
c$$$      !If any Mkj are zero, add a small amount to achieve finite
c$$$      !time steps
c$$$      do k=1,ibins
c$$$         do j=1,icomp-idiag
c$$$            if (Mk(k,j) .eq. 0.d0) then
c$$$               !add a small amount of mass
c$$$               mtotal=0.d0
c$$$               do jj=1,icomp-idiag
c$$$                  mtotal=mtotal+Mk(k,jj)
c$$$               enddo
c$$$               Mk(k,j)=1.d-10*mtotal
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$      !Choose time step
c$$$      dts=dt-tsum      !try to take entire remaining time step
c$$$cdbg      limit='comp'
c$$$      do k=1,ibins
c$$$         if (Nk(k) .gt. Neps) then
c$$$            !limit rates of change for this bin
c$$$            if (dNdt(k) .lt. 0.0) tlimit=dtlimit
c$$$            if (dNdt(k) .gt. 0.0) tlimit=itlimit
c$$$            if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then 
c$$$               dts=Nk(k)*tlimit/abs(dNdt(k))
c$$$
c$$$C      write(79,*) 'k,tlimit,dts',k,tlimit,dts
c$$$cdbg               limit='number'
c$$$cdbg               write(limit(8:9),'(I2)') k
c$$$cdbg               write(*,*) Nk(k), dNdt(k)
c$$$            endif
c$$$            do j=1,icomp-idiag
c$$$               if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
c$$$               if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
c$$$               if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then 
c$$$               mtotal=0.d0
c$$$               do jj=1,icomp-idiag
c$$$                  mtotal=mtotal+Mk(k,jj)
c$$$               enddo
c$$$               !only use this criteria if this species is significant
c$$$               if ((Mk(k,j)/mtotal) .gt. 1.d-5) then
c$$$                  dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
c$$$
c$$$C      write(79,*) 'k,j,tlimit,dts',k,j,tlimit,dts
c$$$               else
c$$$                  if (dMdt(k,j) .lt. 0.0) then
c$$$                     !set dmdt to 0 to avoid very small mk going negative
c$$$                     dMdt(k,j)=0.0
c$$$                  endif
c$$$               endif
c$$$cdbg                  limit='mass'
c$$$cdbg                  write(limit(6:7),'(I2)') k
c$$$cdbg                  write(limit(9:9),'(I1)') j
c$$$cdbg                  write(*,*) Mk(k,j), dMdt(k,j)
c$$$               endif
c$$$            enddo
c$$$         else
c$$$            !nothing in this bin - don't let it affect time step
c$$$            Nk(k)=Neps
c$$$            Mk(k,srtso4)=Neps*1.4*xk(k) !make the added particles SO4
c$$$            !make sure mass/number don't go negative
c$$$            if (dNdt(k) .lt. 0.0) dNdt(k)=0.0
c$$$            do j=1,icomp-idiag
c$$$               if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$c      if (dts .lt. 20.) write(*,*) 'dts<20. in multicoag'
c$$$       if (dts .eq. 0.) then
c$$$          write(*,*) 'time step is 0'
c$$$C          pause
c$$$          stop
c$$$C       go to 20
c$$$       endif
c$$$
c$$$      !Change Nk and Mk
c$$$cdbg      write(*,*) 't=',tsum+dts,' ',limit
c$$$      do k=1,ibins
c$$$         Nk(k)=Nk(k)+dNdt(k)*dts
c$$$         do j=1,icomp-idiag
c$$$            Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
c$$$         enddo
c$$$      enddo
c$$$c      print *, 'tsum=', tsum, 'dts=', dts
c$$$      !Update time and repeat process if necessary
c$$$      tsum=tsum+dts
c$$$c      print *, 'tsum=', tsum
c$$$      if (tsum .lt. dt) goto 10
c$$$
c$$$      RETURN
c$$$      END subroutine multicoag



!@sum cond_nuc   :  calculates the change in the aerosol size distribution
!@+   due to so4 condensation and binary/ternary nucleation during the
!@+   overal microphysics timestep.

!@auth   Jeff Pierce, May 2007

!@+   Initial values of
!@+   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Nnuci - number of nucleation size particles per size bin in grid cell
!@var   Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water
!@var   H2SO4rate - rate of H2SO4 chemical production [kg s^-1]
!@var   dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

!@var   Nkf, Mkf, Gcf - same as above, but final values
!@var   Nknuc, Mknuc - same as above, final values from just nucleation
!@var   Nkcond, Mkcond - same as above, but final values from just condensation
!@var   fn, fn1

      SUBROUTINE cond_nuc(Nki,Mki,Gci,Nkf,Mkf,Gcf,fnavg,fn1avg,
     &     H2SO4rate,dti,num_iter,Nknuc,Mknuc,Nkcond,Mkcond,lev)            

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      IMPLICIT NONE

      real*8 Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      real*8  Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      real*8 Nknuc(ibins), Mknuc(ibins, icomp)
      real*8 Nkcond(ibins),Mkcond(ibins,icomp)
      real*8 H2SO4rate
      real dti
      real*8 fnavg    ! nucleation rate of clusters cm-3 s-1
      real*8 fn1avg   ! formation rate of particles to first size bin cm-3 s-1
      integer, intent(in) :: lev
      real*8 dt
      integer i,j,k,c,jc           ! counters
      real*8 fn       ! nucleation rate of clusters cm-3 s-1
      real*8 fn1      ! formation rate of particles to first size bin cm-3 s-1
      real*8 CSi,CSa   ! intial and average condensation sinks
      real*8 CS1,CS2       ! guesses for condensation sink [s^-1]
      real*8 CStest   !guess for condensation sink
      real*8 Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
      real*8 Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
      real*8 Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
      logical nflg ! returned from nucleation, says whether nucleation occurred or not
      real*8 mcond,mcond1    !mass to condense [kg]
      real*8 tol      !tolerance
      real*8 eps      !small number
      real*8 sinkfrac(ibins) !fraction of condensation sink coming from bin k
      real*8 totmass  !the total mass of H2SO4 generated during the timestep
      real*8 tmass
      real*8 CSch     !fractional change in condensation sink
      real*8 CSch_tol !tolerance in change in condensation sink
      real*8 addt     !adaptive timestep time
      real*8 time_rem !time remaining
      integer num_iter !number of iteration
      real*8 sumH2SO4 !used for finding average H2SO4 conc over timestep
      integer iter ! number of iteration
      real*8 rnuc !critical radius [nm]
      real*8 gasConc  !gas concentration [kg]
      real*8 mass_change !change in mass during nucleation.f
      real*8 total_nh4_1,total_nh4_2
      real*8 min_tstep !minimum timestep [s]
      integer nuc_bin           ! the nucleation bin
      real*8 sumfn, sumfn1 ! used for getting average nucleation rates
      real*8 mcond_soa
      parameter(eps=1E-30)
      parameter(CSch_tol=0.01)
      parameter(min_tstep=1.0d0)


      dt = dble(dti)

C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gc1(j)=Gci(j)
      enddo
      do k=1,ibins
         Nk1(k)=Nki(k)
         Nknuc(k)=Nki(k)
         Nkcond(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
            Mknuc(k,j)=Mki(k,j)
            Mkcond(k,j)=Mki(k,j)
         enddo
      enddo

C     Get initial condensation sink
      CS1 = 0.d0
      call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac)

C     Get initial H2SO4 concentration guess (assuming no nucleation)
C     Make sure that H2SO4 concentration doesn't exceed the amount generated
C     during that timestep (this will happen when the condensation sink is very low)

C     get the steady state H2SO4 concentration
      call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,lev)
      Gc1(srtso4) = gasConc
      addt = min_tstep

      totmass = H2SO4rate*addt*96.d0/98.d0

C     Get change size distribution due to nucleation with initial guess  
      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin,
     &     addt,lev)          

! for so4
      mass_change = 0.d0

      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
      enddo
      mcond = totmass-mass_change ! mass of h2so4 to condense

      if (mcond.lt.0.d0)then
         tmass = 0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tmass = tmass + Mk2(k,j)
            enddo
         enddo

         if (abs(mcond).gt.totmass*1.0d-8) then
            if (-mcond.lt.Mk2(nuc_bin,srtso4)) then

               tmass = 0.d0
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(nuc_bin,j)
               enddo
               Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
               Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
               mcond = 0.d0
            else
               print*,'budget fudge 2 in cond_nuc.f'
               do k=2,ibins
                  Nk2(k) = Nk1(k)
                  Mk2(k,srtso4) = Mk1(k,srtso4)
               enddo
               Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
               Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
               mcond = 0.d0        

            endif
         else
            mcond = 0.d0
         endif
      endif
      
      mass_change = 0.d0

      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtocil)-Mk1(k,srtocil))
      enddo

      mcond_soa = SOArate*addt-mass_change ! mass of soa to condense
      tmass = 0.d0
      do k=1,ibins-1
         do j=1,icomp-idiag
            tmass = tmass + Mk2(k,j)
         enddo
      enddo
      if (mcond_soa.gt.tmass)then ! limit soa
         mcond_soa=tmass
      endif

C     Get guess for condensation
      call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3)

      if(mcond_soa.eq.0) goto 17
      do k=1,ibins
         Nk2(k)=Nk3(k)
         do j=1,icomp
            Mk2(k,j)=Mk3(k,j)
         enddo
      enddo
      call ezcond(Nk2,Mk2,mcond_soa,srtocil,Nk3,Mk3)

 17   continue

      Gc3(srtnh4) = Gc1(srtnh4)   

      call eznh3eqm(Gc3,Mk3)
      call ezwatereqm(Mk3)

! check to see how much condensation sink changed
      call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac) !problem here. 


      if(CS2.eq.CS1.or.CS1.le.0.) then
!MODELE-TOMAS : Somehow CSch becomes NaN, which mean no condensation and nucleation occurs.
! I take whole timestep in this case. Need to look this further in future. 
        time_rem=dt 
        addt=dt  !This will let it go only one. 

      else
         
         CSch = abs(CS2 - CS1)/CS1    
       
c      if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
         ! get starting adaptive timestep to not allow condensationk sink
         ! to change that much
         addt = addt*CSch_tol/CSch/2
         addt = min(addt,dt)
         addt = max(addt,min_tstep)
         time_rem = dt ! time remaining
         
      endif
         num_iter = 0
         sumH2SO4=0.d0
         sumfn = 0.d0
         sumfn1 = 0.d0
         ! do adaptive timesteps
         do while (time_rem .gt. 0.d0)
            num_iter = num_iter + 1
C     get the steady state H2SO4 concentration
            if (num_iter.gt.1)then ! no need to recalculate for first step
               call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,lev)
               Gc1(srtso4) = gasConc
            endif

            sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
            totmass = H2SO4rate*addt*96.d0/98.d0

            call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,
     &           nuc_bin,addt,lev) 
            sumfn = sumfn + fn*addt
            sumfn1 = sumfn1 + fn1*addt

            mass_change = 0.d0
            do k=1,ibins
               mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
            enddo
            mcond = totmass-mass_change ! mass of h2so4 to condense

            if (mcond.lt.0.d0)then
               tmass = 0.d0
               do k=1,ibins
                  do j=1,icomp-idiag
                     tmass = tmass + Mk2(k,j)
                  enddo
               enddo
               if (abs(mcond).gt.totmass*1.0D-8) then
                  if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
c                     if (CS1.gt.1.0D-5)then
c                        print*,'budget fudge 1 in cond_nuc.f'
c                     endif
                     tmass = 0.d0
                     do j=1,icomp-idiag
                        tmass = tmass + Mk2(nuc_bin,j)
                     enddo
                     Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                     Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                     mcond = 0.d0
                  else
                     print*,'budget fudge 2 in cond_nuc.f'
                     do k=2,ibins
                        Nk2(k) = Nk1(k)
                        Mk2(k,srtso4) = Mk1(k,srtso4)
                     enddo
                     Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
                     Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
                     mcond = 0.d0 
                  endif
               else
                  mcond = 0.d0
               endif
            endif

            mass_change = 0.d0

            do k=1,ibins
               mass_change=mass_change+(Mk2(k,srtocil)-Mk1(k,srtocil))
            enddo

            mcond_soa = SOArate*addt-mass_change ! mass of soa to condense
            tmass = 0.d0
            do k=1,ibins-1
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(k,j)
               enddo
            enddo
            if (mcond_soa.gt.tmass)then  ! limit soa
               mcond_soa=tmass
            endif
            
            do k=1,ibins
               Nknuc(k) = Nknuc(k)+Nk2(k)-Nk1(k)
               do j=1,icomp-idiag
                  Mknuc(k,j)=Mknuc(k,j)+Mk2(k,j)-Mk1(k,j)
               enddo
            enddo          
           
            call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3)
            do k=1,ibins
               Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
               do j=1,icomp-idiag
                  Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
               enddo
            enddo

            if(mcond_soa.eq.0) goto 19
            do k=1,ibins
               Nk2(k)=Nk3(k)
               do j=1,icomp
                  Mk2(k,j)=Mk3(k,j)
               enddo
            enddo

            call ezcond(Nk2,Mk2,mcond_soa,srtocil,Nk3,Mk3)
            do k=1,ibins
               Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
               do j=1,icomp-idiag
                  Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
               enddo
            enddO
 19         continue

            Gc3(srtnh4) = Gc1(srtnh4)

            call eznh3eqm(Gc3,Mk3)
            call ezwatereqm(Mk3)
            
! check to see how much condensation sink changed
            call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)  

            time_rem = time_rem - addt
            if (time_rem .gt. 0.d0) then
               
!MODELE-TOMAS : newly added to prevent NaN for CSch. 
!Somehow, mcond >0 but not nucleation and no condensation occurs. 
!same size distribution causes problems. 

               if(CS2.eq.CS1.or.CS1.le.0.) then
                  addt = time_rem
                  addt = max(addt,min_tstep)
               else
                  CSch = abs(CS2 - CS1)/CS1 !TOMAS - CS1 = ZERO?
                  
                  
                  addt = min(addt*CSch_tol/CSch,addt*1.5d0) ! allow adaptive timestep to change
                  addt = min(addt,time_rem) ! allow adaptive timestep to change
                  addt = max(addt,min_tstep)
               endif
Cjrp               endif
               CS1 = CS2
               Gc1(srtnh4)=Gc3(srtnh4)
               do k=1,ibins
                  Nk1(k)=Nk3(k)
                  do j=1,icomp
                     Mk1(k,j)=Mk3(k,j)
                  enddo
               enddo         
            endif
         enddo
         Gcf(srtso4)=sumH2SO4/dt
         fnavg = sumfn/dt
         fn1avg = sumfn1/dt

      do k=1,ibins
         Nkf(k)=Nk3(k)
         do j=1,icomp
            Mkf(k,j)=Mk3(k,j)
         enddo
      enddo      
      Gcf(srtnh4)=Gc3(srtnh4)

      return
      end


!@sum getCondSink_kerm :  calculates the condensation sink (first order loss
!@+    rate of condensing gases) from the aerosol size distribution.
!@+    This is the cond sink in kerminen et al 2004 Parameterization for 
!@+    new particle formation AS&T Eqn 6.
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nk(ibins) - number of particles per size bin in grid cell
!@var   Nnuc - number of particles per size bin in grid cell
!@var   Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mk(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   spec - number of the species we are finding the condensation sink for

!   Output 
!   =================
!@var   CS - condensation sink [s^-1]
!@var   sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getCondSink_kerm(Nko,Mko,CS,Dpmean,Dp1,dens1)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 CS
      real*8 Dpmean             ! the number mean diameter [m]
      real*8 Dp1                ! the size of the first size bin [m]
      real*8 dens1              ! the density of the first size bin [kg/m3]

      integer i,j,k,c           ! counters
      real*8 mu                  !viscosity of air (kg/m s)
      real*8 mfp                 !mean free path of air molecule (m)

      real Di                   !diffusivity of gas in air (m2/s)
      real*8 Neps               !tolerance for number
      real density              !density [kg m^-3]
      real*8 mp                 !mass per particle [kg]
      real*8 Dpk(ibins)         !diameter of particle [m]
      real*8 Kn                 !Knudson number
      real*8 beta(ibins)        !non-continuum correction factor
      real*8 Mktot              !total mass in bin [kg]
      real*8 Dtot,Ntot          ! used on getting the number mean diameter
      
      real mso4, mh2o, mno3, mnh4 !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      real gasdiff
      external aerodens         !!, gasdiff
      
      parameter(Neps=1.0d10)
      real*8 alpha(icomp)       ! accomodation coef  
      data alpha/0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65/
      real Sv(icomp)            !parameter used for estimating diffusivity
      data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88,42.88,
     &     42.88/
      
      
C     get some parameters  
!!      mu=2.5277e-7*temp**0.75302
!!      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6
!for old debugging      mfp=2.0*mu/(pres*sqrt(8.0*0.6589/(pi*gasc*temp)))  !S&P eqn 8.6


cyhl the following should be used instead of above two lines!
      Di=gasdiff(temp,pres,98.0,Sv(srtso4)) ! Di is diffusivity of condensing gas in air [m2/s]
      mfp=2.d0*Di/sqrt(8.0*gasc*temp/(pi*(molwt(srtso4)+2.)/1000.))   !0.098 for H2SO4 m.w. [kg/mol]
!     the denominator is mean speed. sqrt(8.0*R*temp/(pi*molwt(srtso4)/1000.)) !S&P2 eqn 9.2 ms=mean speed [m/s]
         
cyhl but this needs some tuning before it actually uses! 

c      Di=gasdiff(temp,pres,98.0,Sv(srtso4))
c      print*,'Di',Di

C     get size dependent values
      CS = 0.d0
      Ntot = 0.d0
      Dtot = 0.d0
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
               Mktot=Mktot+Mko(k,j)
            enddo

            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4) !0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate                  
            mp=Mktot/Nko(k)
          else
!nothing in this bin - set to "typical value"
            density=1500.
            mp=sqrt(xk(k+1)*xk(k))
          endif
          Dpk(k)=((mp/density)*(6./pi))**(0.333)
          Kn=2.0*mfp/Dpk(k)     !S&P eqn 11.35 (text)
          CS=CS+0.5d0*(Dpk(k)*Nko(k)/(boxvol*1.0D-6)*(1+Kn))/
     &         (1.d0+0.377d0*Kn+1.33d0*Kn*(1+Kn))
          Ntot = Ntot + Nko(k)
          Dtot = Dtot + Nko(k)*Dpk(k)
          if (k.eq.1)then
            Dp1=Dpk(k)
            dens1 = density
          endif
        enddo      
        
        if (Ntot.gt.1D15)then
          Dpmean = Dtot/Ntot
        else
          Dpmean = 150.d0
        endif
        
        return
        end
      
      
!@sum getCondSink :  calculates the condensation sink (first order loss
!@+   rate of condensing gases) from the aerosol size distribution.
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nk(ibins) - number of particles per size bin in grid cell
!@var   Nnuc - number of particles per size bin in grid cell
!@var   Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mk(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   spec - number of the species we are finding the condensation sink for

!   Output
!   =================
!@var   CS - condensation sink [s^-1]
!@var   sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getCondSink(Nko,Mko,spec,CS,sinkfrac)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 CS, sinkfrac(ibins)
      integer spec

      integer i,j,k,c           ! counters
      real*8 mu                 !viscosity of air (kg/m s)
      real*8 mfp                !mean free path of air molecule (m)
      real Di                   !diffusivity of gas in air (m2/s), and molecular weight (kg/mol)
      real*8 Neps               !tolerance for number
      real*8 density            !density [kg m^-3]
      real mw
      real*8 mp                 !mass per particle [kg]
      real*8 Dpk(ibins)         !diameter of particle [m]
      real*8 Kn                 !Knudson number
      real*8 beta(ibins)        !non-continuum correction factor
      real*8 Mktot              !total mass in bin [kg]

      real mso4, mh2o, mno3, mnh4 !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      real gasdiff
      external aerodens         !, gasdiff

      parameter(Neps=1.0d10)
      real*8 alpha(icomp)       ! accomodation coef  
      data alpha/0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65/
      real Sv(icomp)            !parameter used for estimating diffusivity
      data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88,42.88,
     &         42.88/


cyhl the following two lines are commented out (not accurate mfp)
C     get some parameters
c$$$      mu=2.5277e-7*temp**0.75302 
c$$$      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6 !bug?


cyhl mfp is now for the condensing gas in the air   10/17/2010
      Di=gasdiff(temp,pres,98.0,Sv(spec))  ! YHL(10/17/2010) - this is not accurate for SOA, but leave this for now. 
      mfp=2.d0*Di/sqrt(8.0*gasc*temp/(pi*(molwt(srtso4)+2.)/1000.))   
!     molwt(srtso4) is 96, so adding 2 will make 98. 
!     the denominator is mean speed. sqrt(8.0*R*temp/(pi*molwt(srtso4)/1000.)) !S&P2 eqn 9.2 ms=mean speed [m/s]

C     get size dependent values
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
                  Mktot=Mktot+Mko(k,j)
            enddo

            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate                  
            mp=Mktot/Nko(k)
          else
!nothing in this bin - set to "typical value"
            density=1500.d0
            mp=sqrt(xk(k+1)*xk(k))
          endif
          Dpk(k)=((mp/density)*(6.d0/pi))**(1.d0/3.d0)
          Kn=2.0*mfp/Dpk(k)     !S&P eqn 11.35 (text)
          beta(k)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(spec)) !S&P eqn 11.35
        enddo      
C     get condensation sink
        CS = 0.d0
        surf_area = 0.d0
        do k=1,ibins
          CS = CS + Dpk(k)*Nko(k)*beta(k)
          surf_area = surf_area+Nko(k)*pi*(Dpk(k)*1.0D6)**2
        enddo
        do k=1,ibins
          sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
        enddo
c     CS = 2.d0*pi*Di*CS/(boxvol*1D-6)      
        CS = 2.d0*pi*dble(Di)*CS/(boxvol*1D-6)
        surf_area = surf_area/boxvol
        return
        end
      

!@sum getCoagLoss  :  calculates the first order loss rate of particles
!@+   of a given size with respect to coagulation.
!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!@var   d1: diameter in [m] of particle

C-----OUTPUTS-----------------------------------------------------------

!@var   ltc: first order loss rate with respect to coagulation [s-1]
!@var   sinkfrac: fraction of sink from each bin

      SUBROUTINE getCoagLoss(d1,ltc,Nko,Mko,sinkfrac)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 d1       ! diameter of the particle [m]
      real*8 ltc    ! first order loss rate [s-1]
      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 sinkfrac(ibins)

      integer i,k,j
      real*8 density  ! density of particles [kg/m3]
      real*8 density1 ! density of particles in first bin
      real*8 MW, kB
      real*8 mu       !viscosity of air (kg/m s)
      real*8 mfp      !mean free path of air molecule (m)
      real*8 Kn       !Knudsen number of particle
      real*8 Mktot    !total mass of aerosol
      real*8 kij(ibins)
      real*8 Dpk(ibins) !diameter (m) of particles in bin k
      real*8 Dk(ibins),Dk1 !Diffusivity (m2/s) of bin k particles
      real*8 ck(ibins),ck1 !Mean velocity (m/2) of bin k particles
      real*8 neps
      real*8 meps
      real*8 mp       ! mass of the particle (kg)
      real*8 beta     !correction for coagulation coeff.
      real*8 kij_self !coagulation coefficient for self coagulation

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      external aerodens

      parameter (kB= 1.38E-23) !pi and gas constant (J/mol K)
      parameter (neps=1E8, meps=1E-8)

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6

C Calculate particle sizes and diffusivities
      do k=1,ibins
         Mktot = 0.d0
         do j=1, icomp
            Mktot=Mktot+Mko(k,j)
         enddo
         Mktot=Mktot+2.d0*Mko(k,srtso4)/96.d0-Mko(k,srtnh4)/18.d0 ! account for h+

         if (Mktot.gt.meps)then
            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
         else
            density = 1400.
         endif
         if(Nko(k).gt.neps .and. Mktot.gt.meps)then
            mp=Mktot/Nko(k)
         else
            mp=sqrt(xk(k)*xk(k+1))
         endif
         if (k.eq.1) density1 = density
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo
      
      Kn = 2.0*mfp/d1
      Dk1 = kB*temp/(3.0*pi*mu*d1)             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
      mp = 4.d0/3.d0*pi*(d1/2.d0)**3.d0*density1
      ck1=sqrt(8.0*kB*temp/(pi*mp)) !S&P Table 12.1      

      do i=1,ibins
         Kn=4.0*(Dk(i)+Dk1)          
     &        /(sqrt(ck(i)**2+ck1**2)*(Dpk(i)+d1)) !S&P eqn 12.51
         beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
                                !This is S&P eqn 12.46 with non-continuum correction, beta
         kij(i)=2.0*pi*(Dpk(i)+d1)*(Dk(i)+Dk1)*beta
         kij(i)=kij(i)*1.0e6/boxvol !normalize by grid cell volume
c         kij(i)=kij(i)*1.0e6 !cm3/s
      enddo
Cjrp
Cjrp      !self coagulation
Cjrp      Kn=4.0*(Dk1+Dk1)          
Cjrp     &     /(sqrt(ck1**2+ck1**2)*(d1+d1)) !S&P eqn 12.51
Cjrp      beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
Cjrp                                !This is S&P eqn 12.46 with non-continuum correction, beta
Cjrp      kij_self=2.0*pi*(d1+d1)*(Dk1+Dk1)*beta
Cjrp      print*,'kij_self',kij_self*1E6
Cjrp      kij_self=kij_self*1.0e6/boxvol !normalize by grid cell volume      
         
      ltc = 0.d0

      do i=1,ibins
         ltc = ltc + kij(i)*Nko(i)
      enddo
      do i=1,ibins
         sinkfrac(i) = kij(i)*Nko(i)/ltc
      enddo

      RETURN
      END


!@sum NH3_GISStoTOMAS   :  puts ammonia to the particle phase until 
!@+   there is 2 moles of ammonium per mole of sulfate and the remainder
!@+   of ammonia is left in the gas phase.
!@auth   Jeff Pierce, April 2007

      SUBROUTINE NH3_GISStoTOMAS(giss_nh3g,giss_nh4a,Gce,Mke)

      USE TOMAS_AEROSOL

      IMPLICIT NONE

      real*8 Gce(icomp-1)
      real*8 Mke(ibins,icomp)
      real*8 giss_nh3g, giss_nh4a

      integer k
      real*8 tot_nh3  !total kmoles of ammonia
      real*8 tot_so4  !total kmoles of so4
      real*8 sfrac    !fraction of sulfate that is in that bin

      ! get the total number of kmol nh3
      tot_nh3 = giss_nh3g/17.d0 + giss_nh4a/18.d0

      ! get the total number of kmol so4
      tot_so4=0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then ! no free ammonia
        Gce(srtnh4) = 0.d0      ! no gas phase ammonia
        do k=1,ibins
          sfrac = Mke(k,srtso4)/96.d0/tot_so4
          Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
        enddo
      else                      ! free ammonia
c     Mnuce(srtnh4) = Mnuce(srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
        do k=1,ibins
          Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
        enddo
        Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
      endif
      
      RETURN
      END


!@sum getNucRate  :  calls the Vehkamaki 2002 and Napari 2002 nucleation
!@+   parameterizations and gets the binary and ternary nucleation rates.
!@auth   Jeff Pierce, April 2007

!   Initial values of
!   =================

!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water

!-----OUTPUTS-----------------------------------------------------------

!@var   fn - nucleation rate [# cm-3 s-1]
!@var   rnuc - radius of nuclei [nm]
!@var   nflg - says if nucleation happend

      SUBROUTINE getNucRate(Gci,fn,mnuc,nflg,l)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi  
      USE TRACER_COM, only : xk

      IMPLICIT NONE

      integer j,i,k,l
      real*8 Gci(icomp-1)
      real*8 fn       ! nucleation rate to first bin cm-3 s-1
      real*8 mnuc     !mass of nucleating particle [kg]
      logical nflg

      real*8 nh3ppt   ! gas phase ammonia in pptv
      real*8 h2so4    ! gas phase h2so4 in molec cc-1
      real*8 gtime    ! time to grow to first size bin [s]
      real*8 ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      real*8 Mktot    ! total mass in bin
      real*8 neps
      real*8 meps
      real*8 density  ! density of particle [kg/m3]
      real*8 frac     ! fraction of particles growing into first size bin
      real*8 d1,d2    ! diameters of particles [m]
      real*8 mp       ! mass of particle [kg]
      real*8 mold     ! saved mass in first bin
      real*8 rnuc     ! critical nucleation radius [nm]
      real*8 sinkfrac(ibins) ! fraction of loss to different size bins
      real*8 nadd     ! number to add
      real*8 CS       ! kerminan condensation sink [m-2]
      real*8 Dpmean   ! the number wet mean diameter of the existing aerosol
      real*8 Dp1      ! the wet diameter of bin 1
      real*8 dens1    ! density in bin 1 [kg m-3]
      real*8 GR       ! growth rate [nm hr-1]
      real*8 gamma,eta ! used in kerminen 2004 parameterzation
      real*8 drymass,wetmass,WR
      real*8 fn_c     ! barrierless nucleation rate
      real*8 h1,h2,h3,h4,h5,h6
      parameter (neps=1E8, meps=1E-8)


      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      rnuc = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then

            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc 
            if (bin_nuc.eq.1) then

               if((actv_nuc.eq.1).and.(l.le.3))then
                  call bl_nucl(h2so4,fn,rnuc)
               else
                  call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
               endif
            endif
              
            if (fn.gt.1.0d-6)then
               nflg=.true.
             else
                fn = 0.d0
                nflg=.false.
             endif
          endif
         call cf_nucl(temp,rh,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max for ternary
         fn = min(fn,fn_c)    
      else
         nflg=.false.
      endif

      if (fn.gt.0.d0) then
         call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1)
         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mk(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mk(1,j)
         enddo
         WR = wetmass/drymass
         
         call getGrowthTime(d1,Dp1,Gci(srtso4)*(1.d0+soa_amp)*WR,temp,
     &        boxvol,dens1,gtime)
         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(temp/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR

         fn = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))

         mnuc = sqrt(xk(1)*xk(2))
      endif

      return
      end
      

!@sum getH2SO4conc :  uses newtons method to solve for the steady state 
!@+   H2SO4 concentration when nucleation is occuring.

!@+   It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)

!@+   where P is the production rate of H2SO4, CS is the condensation sink
!@+   and M(H2SO4) is the loss of mass towards making new particles.

!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
!@var   CS - condensation sink [s-1]
!@var   NH3conc - ammonium in box [kg box-1]
!@var   prev - logical flag saying if a previous guess should be used or not
!@var   gasConc_prev - the previous guess [kg/box] (not used if prev is false)

!-----OUTPUTS-----------------------------------------------------------

!@var   gasConc - gas H2SO4 [kg/box]

      SUBROUTINE getH2SO4conc(H2SO4rate,CS,NH3conc,gasConc,level)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 H2SO4rate
      real*8 CS
      real*8 NH3conc
      real*8 gasConc
      logical prev
      real*8 gasConc_prev

      integer, intent(in) :: level  ! vertical layer 
      integer i,j,k,c           ! counters
      real*8 fn, rnuc           ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      real*8 mnuc, mnuc1        ! mass of nucleated particle [kg]
      real*8 fn1, rnuc1         ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      real*8 res                ! d[H2SO4]/dt, need to find the solution where res = 0
      real*8 massnuc            ! mass being removed by nucleation [kg s-1 box-1]
      real*8 gasConc1           ! perturbed gasConc
      real*8 gasConc_hi, gasConc_lo
      real*8 res1               ! perturbed res
      real*8 res_new            ! new guess for res
      real*8 dresdgasConc       ! derivative for newtons method
      real*8 Gci(icomp-1)       !array to carry gas concentrations
      logical nflg              !says if nucleation occured

      real*8 H2SO4min           !minimum H2SO4 concentration in parameterizations (molec/cm3)
!     real*8 pi
      integer iter,iter1
      real*8 CSeps              ! low limit for CS
      real*8 max_H2SO4conc      !maximum H2SO4 concentration in parameterizations (kg/box)
      real*8 nh3ppt             !ammonia concentration in ppt
     
      parameter(H2SO4min=1.D4) !molecules cm-3
      parameter(CSeps=1.0d-20)


      do i=1,icomp-1
         Gci(i)=0.d0
      enddo
      Gci(srtnh4)=NH3conc
! make sure CS doesn't equal zero
! some specific stuff for napari vs. vehk
      if ((bin_nuc.eq.1).or.(tern_nuc.eq.1))then
         nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)
         if ((nh3ppt.gt.1.0d0).and.(tern_nuc.eq.1))then
            max_H2SO4conc=1.0D9*boxvol/1000.d0*98.d0/6.022d23
         elseif (bin_nuc.eq.1)then
            max_H2SO4conc=1.0D11*boxvol/1000.d0*98.d0/6.022d23
         else
            max_H2SO4conc = 1.0D100
         endif
      else
         max_H2SO4conc = 1.0D100
      endif
      
C     Checks for when condensation sink is very small
      if (CS.gt.CSeps) then
         gasConc = H2SO4rate/CS
      else
         if ((bin_nuc.eq.1).or.(tern_nuc.eq.1)) then
            gasConc = max_H2SO4conc
         else
            print*,'condesation sink too small in getH2SO4conc.f'
            call stop_model('CS too small getH2SO4',255)
         endif
      endif
      
      gasConc = min(gasConc,max_H2SO4conc)
      Gci(srtso4) = gasConc
      call getNucRate(Gci,fn,mnuc,nflg,level)
      
      if (fn.gt.0.d0) then      ! nucleation occured
         gasConc_lo = H2SO4min*boxvol/(1000.d0/98.d0*6.022d23) !convert to kg/box
         
C     Test to see if gasConc_lo gives a res < 0 (this means ANY nucleation is too high)
         Gci(srtso4) = gasConc_lo*1.000001d0
         call getNucRate(Gci,fn1,mnuc1,nflg,level)
         if (fn1.gt.0.d0) then
            massnuc = mnuc1*fn1*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1350.*fn1*boxvol*
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1800.*fn1*boxvol*
c     &           98.d0/96.d0

            res = H2SO4rate - CS*gasConc_lo - massnuc
            if (res.lt.0.d0) then ! any nucleation too high

               gasConc = gasConc_lo*1.000001 ! have nucleation occur and fix mass balance after
               return
            endif
         endif
         
         gasConc_hi = gasConc   ! we know this must be the upper limit (since no nucleation)
                                !take density of nucleated particle to be 1350 kg/m3
         massnuc = mnuc*fn*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
         res = H2SO4rate - CS*gasConc - massnuc
         
                                ! check to make sure that we can get solution
         if (res.gt.0.d0) then
            return
         endif
         
         iter = 0

         do while ((abs(res/H2SO4rate).gt.1.D-4).and.(iter.lt.40))
            iter = iter+1
            if (res .lt. 0.d0) then ! H2SO4 concentration too high, must reduce
               gasConc_hi = gasConc ! old guess is new upper bound
            elseif (res .gt. 0.d0) then ! H2SO4 concentration too low, must increase
               gasConc_lo = gasConc ! old guess is new lower bound
            endif

            gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
            Gci(srtso4) = gasConc
            call getNucRate(Gci,fn,mnuc,nflg,level)
            massnuc = mnuc*fn*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
            res = H2SO4rate - CS*gasConc - massnuc

            if (iter.eq.30.and.CS.gt.1.0D-5)then
               print*,'getH2SO4conc iter break'
               print*,'H2SO4rate',H2SO4rate,'CS',CS
               print*,'gasConc',gasConc,'massnuc',massnuc
               print*,'res/H2SO4rate',res/H2SO4rate
            endif
         enddo
         
      else                      ! nucleation didn't occur
      endif
      
      return
      end
      



!@sum getGrowthTime  :  calculates the time it takes for a particle to grow
!@+   from one size to the next by condensation of sulfuric acid (and
!@+   associated NH3 and water) onto particles.

!@+   This subroutine  assumes that the growth happens entirely in the kinetic
!@+   regine such that the dDp/dt is not size dependent.  The time for growth 
!@+   to the first size bin may then be approximated by the time for growth via
!@+   sulfuric acid (not including nh4 and water) to the size of the first size bin
!@+   (not including nh4 and water).

!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!@var   d1: intial diameter [m]
!@var   d2: final diameter [m]
!@var     h2so4: h2so4 ammount [kg]
!@var     temp: temperature [K]
!@var     boxvol: box volume [cm3]

C-----OUTPUTS-----------------------------------------------------------

!@var  gtime: the time it takes the particle to grow to first size bin [s]

      SUBROUTINE getGrowthTime(d1,d2,h2so4,temp,boxvol,density,gtime)

      USE CONSTANT,   only:  pi,gasc  
      IMPLICIT NONE

      real*8 d1,d2    ! initial and final diameters [m]
      real*8 h2so4    ! h2so4 ammount [kg]
      real*8 temp     ! temperature [K]
      real*8 boxvol   ! box volume [cm3]
      real*8 gtime    ! the time it will take the particle to grow 
                                ! to first size bin [s]
      real*8 density  ! density of particles in first bin [kg/m3]

      real*8 MW
      real*8 csulf    ! concentration of sulf acid [kmol/m3]
      real*8 mspeed   ! mean speed of molecules [m/s]
      real*8 alpha    ! accomidation coef

      parameter(MW=98.d0) ! density [kg/m3], mol wgt sulf [kg/kmol]
      parameter(alpha=0.65)


      csulf = h2so4/MW/(boxvol*1d-6) ! SA conc. [kmol/m3]
      mspeed = sqrt(8.d0*gasc*temp*1000.d0/(pi*MW))

C   Kinetic regime expression (S&P 11.25) solved for T
      gtime = (d2-d1)/(4.d0*MW/density*mspeed*alpha*csulf)

      RETURN
      END




!@sum nucleation :  calls the Vehkamaki 2002 and Napari 2002 nucleation
!@+   parameterizations and gets the binary and ternary nucleation rates.
!@+   The number of particles added to the first size bin is calculated
!@+   by comparing the growth rate of the new particles to the coagulation
!@+   sink.

!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!   Initial values of
!   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water
!@var   dt - total model time step to be taken (s)

!-----OUTPUTS-----------------------------------------------------------

!@var   Nkf, Mkf, Gcf - same as above, but final values
!@var   fn, fn1

      SUBROUTINE nucleation(Nki,Mki,Gci,Nkf,Mkf,Gcf,fn,fn1,totsulf,
     &     nuc_bin,dt,l)


      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi 
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      integer j,i,k,l
      real*8 Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      real*8 Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      real*8 totsulf
      integer nuc_bin
      real*8 dt
      real*8 fn       ! nucleation rate of clusters cm-3 s-1
      real*8 fn1      ! formation rate of particles to first size bin cm-3 s-1

      real*8 nh3ppt   ! gas phase ammonia in pptv
      real*8 h2so4    ! gas phase h2so4 in molec cc-1
      real*8 rnuc     ! critical nucleation radius [nm]
      real*8 gtime    ! time to grow to first size bin [s]
      real*8 ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      real*8 Mktot    ! total mass in bin
      real*8 neps
      real*8 meps
      real*8 density  ! density of particle [kg/m3]
      real*8 frac     ! fraction of particles growing into first size bin
      real*8 d1,d2    ! diameters of particles [m]
      real*8 mp       ! mass of particle [kg]
      real*8 Dpk(ibins) !diameter of particle [m]
      real*8 mold     ! saved mass in first bin
      real*8 mnuc     !mass of nucleation
      real*8 sinkfrac(ibins) ! fraction of loss to different size bins
      real*8 nadd     ! number to add
      real*8 CS       ! kerminan condensation sink [m-2]
      real*8 Dpmean   ! the number wet mean diameter of the existing aerosol
      real*8 Dp1      ! the wet diameter of bin 1
      real*8 dens1    ! density in bin 1 [kg m-3]
      real*8 GR       ! growth rate [nm hr-1]
      real*8 gamma,eta ! used in kerminen 2004 parameterzation
      real*8 drymass,wetmass,WR

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      external aerodens
      real*8 fn_c     ! barrierless nucleation rate
      real*8 h1,h2,h3,h4,h5,h6
      parameter (neps=1E8, meps=1E-8)


      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      fn1 = 0.d0
      rnuc = 0.d0
      gtime = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then

            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
         elseif (bin_nuc.eq.1) then

            if((actv_nuc.eq.1).and.(l.le.3))then
              call bl_nucl(h2so4,fn,rnuc)
            else
              call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            endif
            if (fn.lt.1.0d-6)then
               fn = 0.d0
            endif
         endif    
         call cf_nucl(temp,rh,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max
         fn = min(fn,fn_c) 
      endif

      d1 = rnuc*2.d0*1D-9
      do k=1,ibins
         if (Nki(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
               Mktot=Mktot+Mki(k,j)
            enddo

            mso4=Mki(k,srtso4) 
            mnacl=Mki(k,srtna)
            mno3=0.0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mki(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mki(k,srtecob)
            mecil=Mki(k,srtecil)
            mocil=Mki(k,srtocil)
            mocob=Mki(k,srtocob)
            mdust=Mki(k,srtdust)          
            mh2o=Mki(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
            
            mp=Mktot/Nki(k)
         else
                                !nothing in this bin - set to "typical value"
            density=1500.d0
            mp=sqrt(xk(k+1)*xk(k))
         endif
         Dpk(k)=((mp/density)*(6.d0/pi))**(1.d0/3.d0)
      enddo

C     if nucleation occured, see how many particles grow to join the first size
C     section
      if (fn.gt.0.d0.and.Dpk(1).gt.d1) then

         call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1)
         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mk(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mk(1,j)
         enddo
         WR = wetmass/drymass
         
         call getGrowthTime(d1,Dp1,Gci(srtso4)*(1.d0+soa_amp)*WR,temp,
     &        boxvol,dens1,gtime)
         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(temp/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR

         fn1 = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))

         mnuc = sqrt(xk(1)*xk(2))
 
         nadd = fn1
         
         nuc_bin = 1
         
         mold = Mki(nuc_bin,srtso4)
         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+nadd*mnuc*
     &        boxvol*dt/(1.d0+soa_amp)
         Mkf(nuc_bin,srtocil) = Mki(nuc_bin,srtocil)+nadd*mnuc*
     &        boxvol*dt*(1.d0-1.d0/(1.d0+soa_amp))
         Nkf(nuc_bin) = Nki(nuc_bin)+nadd*boxvol*dt
         Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
         Gcf(srtnh4) = Gci(srtnh4)
         
         do k=1,ibins
            if (k .ne. nuc_bin)then
               Nkf(k) = Nki(k)
               do i=1,icomp
                  Mkf(k,i) = Mki(k,i)
               enddo
            else
               do i=1,icomp
                  if (i.ne.srtso4) then
                     Mkf(k,i) = Mki(k,i)
                  endif
               enddo
            endif
         enddo
         
         do k=1,ibins
            if (Nkf(k).lt.1.d0) then
               Nkf(k) = 0.d0
               do j=1,icomp
                  Mkf(k,j) = 0.d0
               enddo
            endif
         enddo
         call mnfix(Nkf,Mkf)
         
C     there is a chance that Gcf will go less than zero because we are artificially growing
C     particles into the first size bin.  don't let it go less than zero.

         
      else if (fn.gt.0.d0.and.d1.gt.Dpk(1)) then
         fn1=fn
         k=1
         do while (d1.ge.Dpk(k+1))
            k=k+1
         enddo

         nuc_bin=k
         mnuc=sqrt(xk(nuc_bin)*xk(nuc_bin+1))

         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+fn1*mnuc*
     &        boxvol*dt/(1.d0+soa_amp)
         Mkf(nuc_bin,srtocil) = Mki(nuc_bin,srtocil)+fn1*mnuc*
     &        boxvol*dt*(1.d0-1.d0/(1.d0+soa_amp))
         Nkf(nuc_bin) = Nki(nuc_bin)+fn1*boxvol*dt
         Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
         Gcf(srtnh4) = Gci(srtnh4)

         if(k.gt.5)then
            print*,'big rnuc',Dpk(k),d1,Dpk(k+1)
         endif

         do k=1,ibins
            if (k .ne. nuc_bin)then
               Nkf(k) = Nki(k)
               do i=1,icomp
                  Mkf(k,i) = Mki(k,i)
               enddo
            else
               do i=1,icomp
                  if (i.ne.srtso4) then
                     Mkf(k,i) = Mki(k,i)
                  endif
               enddo
            endif
         enddo
         
         do k=1,ibins
            if (Nkf(k).lt.1.d0) then
               Nkf(k) = 0.d0
               do j=1,icomp
                  Mkf(k,j) = 0.d0
               enddo
            endif
         enddo
         call mnfix(Nkf,Mkf)
         
      else
         
         do k=1,ibins
            Nkf(k) = Nki(k)
            do i=1,icomp
               Mkf(k,i) = Mki(k,i)
            enddo
         enddo
         
      endif
      
      RETURN
      END


!@sum ezcond  :  takes a given amount of mass and condenses it
!@+   across the bins accordingly.  
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
!@var   mcond - mass of species to condense [kg/grid cell]
!@var   spec - the number of the species to condense

!-----OUTPUTS-----------------------------------------------------------

!     Nkf, Mkf - same as above, but final values

      SUBROUTINE ezcond(Nki,Mki,mcondi,spec,Nkf,Mkf)

      USE TOMAS_AEROSOL
      USE DOMAIN_DECOMP_ATM, only : am_i_root 
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nki(ibins), Mki(ibins, icomp)
      real*8 Nkf(ibins), Mkf(ibins, icomp)
      real*8 mcondi
      integer spec

      integer i,j,k,c           ! counters
      real*8 mcond
      real*8 CS                 ! condensation sink [s^-1]
      real*8 sinkfrac(ibins+1)  ! fraction of CS in size bin
      real*8 Nk1(ibins), Mk1(ibins, icomp)
      real*8 Nk2(ibins), Mk2(ibins, icomp)
      real*8 madd               ! mass to add to each bin [kg]
      real*8 maddp(ibins)       ! mass to add per particle [kg]
      real*8 mconds             ! mass to add per step [kg]
      integer nsteps            ! number of condensation steps necessary
      integer floor, ceil       ! the floor and ceiling (temporary)
      real*8 eps                ! small number
      real*8 tdt                !the value 2/3
      real*8 mpo,mpw            !dry and "wet" mass of particle
      real*8 WR                 !wet ratio
      real*8 tau(ibins)         !driving force for condensation
      real*8 totsinkfrac        ! total sink fraction not including nuc bin
      real*8 CSeps              ! lower limit for condensation sink
      real*8 tot_m,tot_s        !total mass, total sulfate mass
      real*8 ratio              ! used in mass correction
      real*8 fracch(ibins,icomp)
      real*8 totch
      real*8 zeros(icomp)
      real*8 tot_i,tot_f,tot_fa ! used for conservation of mass check
      
      parameter(eps=1.d-40)
      parameter(CSeps=1.d-20)


      tdt=2.d0/3.d0

      mcond=mcondi

! initialize variables
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo

      call mnfix(Nk1,Mk1)

! get the sink fractions
      call getCondSink(Nk1,Mk1,spec,CS,sinkfrac) ! set Nnuc to zero for this calc

! make sure that condensation sink isn't too small
      if (CS.lt.CSeps) then     ! just make particles in first bin
         Mkf(1,spec) = Mk1(1,spec) + mcond
         Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
         do j=1,icomp
            if (icomp.ne.spec) then
               Mkf(1,j) = Mk1(1,j)
            endif
         enddo
         do k=2,ibins
            Nkf(k) = Nk1(k)
            do j=1,icomp
               Mkf(k,j) = Mk1(k,j)
            enddo
         enddo
         return
      endif
      
! determine how much mass to add to each size bin
! also determine how many condensation steps we need
      totsinkfrac = 0.d0
      do k=1,ibins
        totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
      enddo
      nsteps = 1
      do k=1,ibins
         if (sinkfrac(k).lt.1.0D-20)then
            madd = 0.d0
         else
            madd = mcond*sinkfrac(k)/totsinkfrac
         endif
         mpo=0.0
         do j=1,icomp-idiag
            mpo=mpo + Mk1(k,j)
         enddo
         floor = int(madd*0.00001/mpo)
         ceil = floor + 1
         nsteps = max(nsteps,ceil) ! don't let the mass increase by more than 10%
      enddo

! mass to condense each step
      mconds = mcond/nsteps
      
! do steps of condensation
      do i=1,nsteps
         if (i.ne.1) then
            call getCondSink(Nk1,Mk1,spec,
     &        CS,sinkfrac)      ! set Nnuc to zero for this calculation
            totsinkfrac = 0.d0
            do k=1,ibins
              totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
            enddo
         endif      
         
         tot_m=0.d0
         tot_s=0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tot_m = tot_m + Mk1(k,j)
               if (j.eq.spec) then
                  tot_s = tot_s + Mk1(k,j)
               endif
            enddo
         enddo

         if (mcond.gt.tot_m*1.0D-3) then

            do k=1,ibins
               mpo=0.0
               mpw=0.0
!WIN'S CODE MODIFICATION 6/19/06
!THIS MUST CHANGED WITH THE NEW dmdt_int.f
               do j=1,icomp-idiag
                  mpo = mpo+Mk1(k,j) !accumulate dry mass
               enddo
               do j=1,icomp
                  mpw = mpw+Mk1(k,j) ! have wet mass include amso4
               enddo
               WR = mpw/mpo     !WR = wet ratio = total mass/dry mass
               if (Nk1(k) .gt. 0.d0) then
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
                  mpw=mpw/Nk1(k)
                  tau(k)=1.5d0*((mpw+maddp(k)*WR)**tdt-mpw**tdt) !added WR to moxid term (win, 5/15/06)
c     tau(k)=0.d0
c     maddp(k)=0.d0
               else
                                !nothing in this bin - set tau to zero
                  tau(k)=0.d0
                  maddp(k) = 0.d0
               endif
            enddo

            call mnfix(Nk1,Mk1)
                                ! do condensation

            call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,maddp)
c     call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
C     jrp         totch=0.0
C     jrp         do k=1,ibins
C     jrp            do j=1,icomp
C     jrp               fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
C     jrp               totch = totch + (Mk2(k,j)-Mk1(k,j))
C     jrp            enddo
C     jrp         enddo

         elseif (mcond.gt.tot_s*1.0D-12) then
            do k=1,ibins
               if (Nk1(k) .gt. 0.d0) then
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac
               else
                  maddp(k) = 0.d0
               endif
               Mk2(k,spec)=Mk1(k,spec)+maddp(k)
               do j=1,icomp
                  if (j.ne.spec) then
                     Mk2(k,j)=Mk1(k,j)
                  endif
               enddo
               Nk2(k)=Nk1(k)
            enddo
            call mnfix(Nk2,Mk2)
         else ! do nothing
            mcond = 0.d0
            do k=1,ibins
               Nk2(k)=Nk1(k)
               do j=1,icomp
                  Mk2(k,j)=Mk1(k,j)
               enddo
            enddo
         endif
         if (i.ne.nsteps)then
            do k=1,ibins
               Nk1(k)=Nk2(k)
               do j=1,icomp
                  Mk1(k,j)=Mk2(k,j)
               enddo
            enddo            
         endif

      enddo

      do k=1,ibins
         Nkf(k)=Nk2(k)
         do j=1,icomp
            Mkf(k,j)=Mk2(k,j)
         enddo
      enddo

! check for conservation of mass
      tot_i = 0.d0
      tot_fa = mcond
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,spec)
         tot_f=tot_f+Mkf(k,spec)
         tot_fa=tot_fa+Mki(k,spec)
      enddo

      if (mcond.gt.0.d0.and.
     &    abs((mcond-(tot_f-tot_i))/mcond).gt.0.d0) then
         if (abs((mcond-(tot_f-tot_i))/mcond).lt.1.d0) then
            ! do correction of mass
            ratio = (tot_f-tot_i)/mcond
            do k=1,ibins
               Mkf(k,spec)=Mki(k,spec)+
     &              (Mkf(k,spec)-Mki(k,spec))/ratio
            enddo
            call mnfix(Nkf,Mkf)
         else
            if(am_i_root())then
            print*,'ERROR in ezcond',spec
            print*,'Condensation error',(mcond-(tot_f-tot_i))/mcond
            print*,'mcond',mcond,'change',tot_f-tot_i
            print*,'tot_i',tot_i,'tot_fa',tot_fa,'tot_f',tot_f
            print*,'Nki',Nki
            print*,'Nkf',Nkf
            print*,'Mki',Mki
            print*,'Mkf',Mkf
            call STOP_model('ERROR in ezcond',255)
            endif
         endif
      endif

! check for conservation of mass
      tot_i = 0.d0
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,srtnh4)
         tot_f=tot_f+Mkf(k,srtnh4)
      enddo
      if(tot_i.gt.0)then !YUNHA LEE  - I added this to avoid floating invalid in MODELE-TOMAS.
      if (abs(tot_f-tot_i)/tot_i.gt.1.0D-8)then
         print*,'No N conservation in ezcond.f'
         print*,'initial',tot_i
         print*,'final',tot_f
      endif
      endif

      return
      end


!@auth  Peter Adams/Modified by Yunha Lee Mar 2008 

!@sum TMCOND :  CONDENSATION Based on Tzivion, Feingold, Levin, JAS 1989 and 
!@+             Stevens, Feingold, Cotton, JAS 1996

!@+  The supersaturation is calculated outside of the routine and assumed
!@+  to be constant at its average value over the timestep.
!@+  
!@+  The method has three basic components:
!@+  (1) first a top hat representation of the distribution is construced
!@+      in each bin and these are translated according to the analytic
!@+      solutions
!@+  (2) The translated tophats are then remapped to bins.  Here if a 
!@+      top hat entirely or in part lies below the lowest bin it is 
!@+      not counted.
!@+  

!@+   Additional notes (Peter Adams)

!@+   I have changed the routine to handle multicomponent aerosols.  The
!@+   arrays of mass moments are now two dimensional (size and species).
!@+   Only a single component (CSPECIES) is allowed to condense during
!@+   a given call to this routine.  Multicomponent condensation/evaporation
!@+   is accomplished via multiple calls.  Variables YLC and YUC are
!@+   similar to YL and YU except that they refer to the mass of the 
!@+   condensing species, rather than total aerosol mass.

!@+   I have removed ventilation variables (VSW/VNTF) from the subroutine
!@+   call.  They still exist internally within this subroutine, but
!@+   are initialized such that they do nothing.

!@+   I have created a new variable, AMKDRY, which is the total mass in
!@+   a size bin (sum of all chemical components excluding water).  I
!@+   have also created WR, which is the ratio of total wet mass to 
!@+   total dry mass in a size bin.

!@+   AMKC(k,j) is the total amount of mass after condensation of species
!@+   j in particles that BEGAN in bin k.  It is used as a diagnostic
!@+   for tracking down numerical errors.

!@+   End of my additional notes

!@var  TAU(k) ......... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
!@var  X(K) ........ Array of bin limits in mass space
!@var  AMKD(K,J) ... Input array of mass moments
!@var  ANKD(K) ..... Input array of number moments
!@var  AMK(K,J) .... Output array of mass moments
!@var  ANK(K) ...... Output array of number moments
!@var  CSPECIES .... Index of chemical species that is condensing

      SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,moxd)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE DOMAIN_DECOMP_ATM, only : am_i_root

      IMPLICIT NONE

      INTEGER L,I,J,K,IMN,CSPECIES,kk
      real*8 DN,DM,DYI,TAU(ibins),XL,XU,YL,YLC,YU,YUC
      real*8 TEPS,NEPS,EX2,ZERO
      real*8 XI,XX,XP,YM,WTH,W1,W2,WW,AVG
      real*8 VSW,VNTF(ibins)
      real*8 TAU_L, maxtau
      real*8 X(ibins+1),AMKD(ibins,icomp),ANKD(ibins)
      real*8 AMK(ibins,icomp),ANK(ibins)
      real*8 AMKDRY(ibins), AMKWET(ibins), WR(ibins)
      real*8 DMDT_INT
      real*8 AMKD_tot
      PARAMETER (TEPS=1.0d-40,NEPs=1.0d-20)
      PARAMETER (EX2=2.d0/3.d0,ZERO=0.0d0)
      real*8 moxd(ibins)        !moxid/Nact (win, 5/25/06)
      real*8 c1, c2             !correction factor (win, 5/25/06)
      real*8 xk_hi,tmpvar,xk_lo,frac_lo_n,frac_lo_m
!      external DMDT_INT

 3    format(I4,200E20.11)



C If any ANKD are zero, set them to a small value to avoid division by zero
      do k=1,ibins
         if (ANKD(k) .lt. NEPS) then
            ANKD(k)=NEPS
            AMKD(k,srtso4)=NEPS*sqrt(X(k)*X(k+1)) !make the added particles SO4
            do j=1,icomp
               if (j.ne.srtso4)then
                  AMKD(k,j)=0.d0
               endif
            enddo
         endif
      enddo

Cpja Sometimes, after repeated condensation calls, the average bin mass
Cpja can be just above the bin boundary - in that case, transfer a some
Cpja to the next highest bin

!TOMAS - If AMKD_tot/ANKD(k) >X (k+1), it will still great after this part. 
! mp is same by moving 10% mass and number.  That's the point here??
      do k=1,ibins-1
        AMKD_tot=0.0d0
        do kk=1,icomp-idiag
        AMKD_tot=AMKD_tot+AMKD(k,kk)
        enddo
         if (AMKD_tot/ANKD(k).gt.X(k+1)) then
            do j=1,icomp
               AMKD(k+1,j)=AMKD(k+1,j)+0.1d0*AMKD(k,j)
               AMKD(k,j)=AMKD(k,j)*0.9d0
            enddo
            ANKD(k+1)=ANKD(k+1)+0.1d0*ANKD(k)
            ANKD(k)=ANKD(k)*0.9d0
         endif
      enddo

Cpja Initialize ventilation variables so they don't do anything
      VSW=0.0d0
      DO L=1,ibins
         VNTF(L)=0.0d0
      ENDDO

Cpja Initialize AMKDRY and WR
      DO L=1,ibins
         AMKDRY(L)=0.0d0
         AMKWET(L)=0.0d0
         DO J=1,icomp-idiag
            AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
         ENDDO
         DO J=1,icomp
            AMKWET(L)=AMKWET(L)+AMKD(L,J)
         ENDDO
         WR(L)=AMKWET(L)/AMKDRY(L)
      ENDDO

Cpja Initialize X() array of particle masses based on xk()
      DO L=1,ibins
         X(L)=xk(L)
      ENDDO

c
c Only solve when significant forcing is available
c
      maxtau=0.0d0
      do l=1,ibins
         maxtau=max(maxtau,abs(TAU(l)))
      enddo
      IF(ABS(maxtau).LT.TEPS)THEN
         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=AMKD(L,J)
            ENDDO
            ANK(L)=ANKD(L)
         ENDDO
      ELSE
         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=0.d0
            ENDDO
            ANK(L)=0.d0
         ENDDO
         WW=0.5d0
c        IF(TAU.LT.0.)WW=.5d0
c
c identify tophats and do lagrangian growth
c
         DO L=1,ibins
            IF(ANKD(L).EQ.0.)GOTO 200

            !if tau is zero, leave everything in same bin
            IF (TAU(L) .EQ. 0.) THEN
               ANK(L)=ANK(L)+ANKD(L)
               DO J=1,icomp
                  AMK(L,J)=AMK(L,J)+AMKD(L,J)
               ENDDO
            ENDIF
            IF (TAU(L) .EQ. 0.) GOTO 200

Cpja Limiting AVG, the average particle size to lie within the size
Cpja bounds causes particles to grow or shrink arbitrarily and is
Cpja wreacking havoc with choosing condensational timesteps and
Cpja conserving mass.  I have turned them off.
c            AVG=MAX(X(L),MIN(X(L+1),AMKDRY(L)/(NEPS+ANKD(L))))
            AVG=AMKDRY(L)/ANKD(L)
            XX=X(L)/AVG
cyhl Eq below needs to be changed for 15 size bins. 
cyhl double check this equation!!!
            if(l.lt.ibins-1)then
               XI=.5d0 + XX*(2.5d0 - 2.0d0*XX)

cyhl this is the old equation for 30 bins: XI=.5d0 + XX*(1.5d0 - XX)
            if (XI .LT. 1.d0) then
               !W1 will have sqrt of negative number
               write(*,*)'ERROR: tmcond - XI<1 for bin: ',L,XK(L)
               write(*,*)'lower limit is',X(L),TAU(L)
               write(*,*)'AVG is ',AVG
               write(*,*)'Nk is ', ANKD(L)
               write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
               write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
               call stop_model('ERROR in tmcond',255)
            endif
cyhl            W1 =SQRT(12.d0*(XI-1.d0))*AVG
cyhl            W2 =MIN(X(L+1)-AVG,AVG-X(L))
cyhl            WTH=W1*WW+W2*(1.d0-WW)
            W1 =SQRT(12.d0*(XI-1.d0))*AVG/4.0d0 ! cyhl 4.0=xk(k+1)/xk(k)
            W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            
            else                   ! 32 = xk(k+1)/xk(k)
            XI=.5d0 + XX*(16.5d0 - 16.0d0*XX)
            if (XI .LT. 1.d0) then
               !W1 will have sqrt of negative number
               write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
               write(*,*)'lower limit is',X(L)
               write(*,*)'AVG is ',AVG
               write(*,*)'Nk is ', ANKD(L)
               write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
               write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
               call stop_model('ERROR in tmcond',255)
            endif
               W1 =SQRT(12.d0*(XI-1.d0))*AVG/32.0d0 ! cyhl 32.0=xk(k+1)/xk(k)
               W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            endif

            WTH=W1*WW+W2*(1.d0-WW)  ! average of w1 and w2
            IF(WTH.GT.1.) then
               write(*,*)'WTH>1 in cond.f, bin #',L,W1,W2
               call stop_model('ERROR in tmcond',255)
            ENDIF
            XU=AVG+WTH*.5d0
            XL=AVG-WTH*.5d0
c Ventilation added bin-by-bin
            TAU_L=TAU(l)*MAX(1.d0,VNTF(L)*VSW)
            IF(TAU_L/TAU(l).GT. 6.) THEN
               PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
            ENDIF
            IF(TAU_L.GT.TAU(l)) THEN 
               PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
            ENDIF
! prior to 5/25/06 (win)
!            YU=DMDT_INT(XU,TAU_L,WR(L))
!            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
!            IF (YU .GT. X(ibins+1) ) THEN
!               YUC=YUC*X(ibins+1)/YU
!               YU=X(ibins+1)
!            ENDIF
!            YL=DMDT_INT(XL,TAU_L,WR(L)) 
!            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
!add new correction factor to YU and YL (win, 5/25/06)
            YU=DMDT_INT(XU,TAU_L,WR(L))
            YL=DMDT_INT(XL,TAU_L,WR(L)) 
            if(moxd(L).eq.0d0) then
               c1=1.d0          !for so4cond call, without correction factor.
            else
               c1 = moxd(L)*2.d0/(YU+YL-XU-XL)
            endif
            c2 = c1 - (c1-1.d0)*(XU+XL)/(YU+YL)
            YU = YU*c2
            YL = YL*c2
!end part for fudging to get higher AVG 

            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
            IF (YU .GT. X(ibins+1) ) THEN
               YUC=YUC*X(ibins+1)/YU
               YU=X(ibins+1)
            ENDIF
            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
            DYI=1.d0/(YU-YL)

c            print*,'XL',XL,'YL',YL,'XU',XU,'YU',YU
c
c deal with portion of distribution that lies below lowest gridpoint
c
            IF(YL.LT.X(1))THEN
cpja Instead of the following, I will just add all new condensed
cpja mass to the same size bin
c               if ((YL/XL-1.d0) .LT. 1.d-3) then
c                  !insignificant growth - leave alone
c                  ANK(L)=ANK(L)+ANKD(L)
c                  DO J=1,icomp-1
c                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
c                  ENDDO
c                  GOTO 200
c               else
c                  !subtract out lower portion
c                  write(*,*)'ERROR in cond - low portion subtracted'
c                  write(*,*) 'Nk,Mk: ',ANKD(L),AMKD(L,1),AMKD(L,2)
c                  write(*,*) 'TAU: ', TAU_L
c                  write(*,*) 'XL, YL, YLC: ',XL,YL,YLC
c                  write(*,*) 'XU, YU, YUC: ',XU,YU,YUC
c                  ANKD(L)=ANKD(L)*MAX(ZERO,(YU-X(1)))*DYI
c                  YL=X(1)
c                  YLC=X(1)*AMKD(1,CSPECIES)/AMKDRY(1)
c                  DYI=1.d0/(YU-YL)
c               endif
               ANK(L)=ANK(L)+ANKD(L)
               do j=1,icomp
                  if (J.EQ.CSPECIES) then
                     AMK(L,J)=AMK(L,J)+(YUC+YLC)*.5d0*ANKD(L)
                  else
                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
                  endif
               enddo
               GOTO 200
            ENDIF
            IF(YU.LT.X(1))GOTO 200
c
c Begin remapping (start search at present location if condensation)
c
            IMN=1
            IF(TAU(l).GT.0.)IMN=L
            DO I=IMN,ibins
               IF(YL.LT.X(I+1))THEN
                  IF(YU.LE.X(I+1))THEN
                     DN=ANKD(L)
                     do j=1,icomp
                        DM=AMKD(L,J)
                        IF (J.EQ.CSPECIES) THEN
                           AMK(I,J)=(YUC+YLC)*.5d0*DN+AMK(I,J)
                        ELSE
                           AMK(I,J)=AMK(I,J)+DM
                        ENDIF
                     enddo
                     ANK(I)=ANK(I)+DN
                  ELSE
                     DN=ANKD(L)*(X(I+1)-YL)*DYI
                     do j=1,icomp
                        DM=AMKD(L,J)*(X(I+1)-YL)*DYI
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(I+1),-1.0d0*TAU_L,WR(L))
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                           AMK(I,J)=DN*(YM+YLC)*0.5d0+AMK(I,J)
                        ELSE
                           AMK(I,J)=AMK(I,J)+DM
                        ENDIF
                     enddo
                     ANK(I)=ANK(I)+DN
                     DO K=I+1,ibins
                        IF(YU.LE.X(K+1))GOTO 100
                        DN=ANKD(L)*(X(K+1)-X(K))*DYI
                        do j=1,icomp
                           DM=AMKD(L,J)*(X(K+1)-X(K))*DYI
                           IF (J.EQ.CSPECIES) THEN
                              XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                              YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                              AMK(K,J)=DN*1.5d0*YM+AMK(K,J)
                           ELSE
                              AMK(K,J)=AMK(K,J)+DM
                           ENDIF
                        enddo
                        ANK(K)=ANK(K)+DN
                     ENDDO
                     print*,'Trying to put stuff in bin ibins+1'
                     call stop_model('ERROR in tmcond,',255)
 100                 CONTINUE
                     DN=ANKD(L)*(YU-X(K))*DYI
                     do j=1,icomp
                        DM=AMKD(L,J)*(YU-X(K))*DYI
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                           AMK(K,J)=DN*(YUC+YM)*0.5d0+AMK(K,J)
                        ELSE
                           AMK(K,J)=AMK(K,J)+DM
                        ENDIF
                     enddo
                     ANK(K)=ANK(K)+DN
                  ENDIF  !YU.LE.X(I+1)
                  GOTO200
               ENDIF   !YL.LT.X(I+1)
            ENDDO !I loop
 200        CONTINUE
         ENDDO    !L loop
      ENDIF

      RETURN
      END


c----------------------------------------------------------------------
!@sum Function DMDT_INT:  Here we apply the analytic solution to the
!@+  droplet growth equation in mass space for a given scale length which
!@+  mimics the inclusion of gas kinetic effects
!@+ Reference: Stevens et al. 1996, Elements of the Microphysical Structure
!@+           of Numerically Simulated Nonprecipitating Stratocumulus,
!@+           J. Atmos. Sci., 53(7),980-1006. 

!@+  This calculates a solution for m(t+dt) using eqn.(A3) from the reference

!@+  Comments by Peter Adams :
!@+  I have changed the length scale.  Non-continuum effects are
!@+  assumed to be taken into account in choice of tau (in so4cond
!@+  subroutine).

!@+  I have also added another argument to the function call, WR.  This
!@+  is the ratio of wet mass to dry mass of the particle.  I use this
!@+  information to calculate the amount of growth of the wet particle,
!@+  but then return the resulting dry mass.  This is the appropriate
!@+  way to implement the condensation algorithm in a moving sectional
!@+  framework.

!@+  End of comments

!@var  M0 ......... initial mass
!@var  L0 ......... length scale
!@var  Tau ........ forcing from vapor field

      real*8 FUNCTION DMDT_INT(M0,TAU,WR)
 
      IMPLICIT NONE
      real*8 M0,TAU,X,L0,C,ZERO,WR,MH2O
      PARAMETER (C=2.d0/3.d0,L0=0.0d0,ZERO=0.0d0)
 
      MH2O=(WR-1.d0)*M0
      X=((M0+MH2O)**C+L0)
      X=MAX(ZERO,SQRT(MAX(ZERO,C*TAU+X))-L0)
!win,5/14/06      DMDT_INT=X*X*X-MH2O
      DMDT_INT = X*X*X/WR  !<step5.2> change calculation to keep WR 
                           !constant after condensation/evap (win, 5/14/06)

Cpja Perform some numerical checks on dmdt_int
      if ((tau .gt. 0.0) .and. (dmdt_int .lt. m0)) dmdt_int=m0
      if ((tau .lt. 0.0) .and. (dmdt_int .gt. m0)) dmdt_int=m0

      RETURN
      END


!@sum cf_nucl  :  calculates the barrierless nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of...
!@+   Reference : Clement and Ford (1999) Atmos. Environ. 33:489-499
!@auth   Jeff Pierce, April 2007

      SUBROUTINE cf_nucl(temp,rh,cna,nh3ppt,fn)

      IMPLICIT NONE

      real*8 cna      ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppt   ! mixing ratio of ammonia in ppt
      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 alpha1


      if (nh3ppt .lt. 0.1) then
         alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
      else
         alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
      endif
      fn = alpha1*cna**2*3600.
c sensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
      if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion

 10   return
      end




!@sum vehk_nucl    :  calculates the binary nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of 

!@+   Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck, 
!@+   M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric 
!@+   Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions." 
!@+   Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).

!@auth   Jeff Pierce, April 2007

      SUBROUTINE vehk_nucl(temp,rh,cnai,fn,rnuc)

      IMPLICIT NONE

      real*8 cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
      real*8 gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 xstar                ! mole fraction sulfuric acid in cluster
      real*8 ntot                 ! total number of molecules in cluster
      integer i                 ! counter

c     Nucleation Rate Coefficients
c
      data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598,
     $          2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
      data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796,
     $       0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
      data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, 
     $       0.00864245, -0.000407382, 0.0000377124, 0.000157098,
     $       -0.0000198417, -9.24619d-6 /
      data fb3 /0.0000722811, 0.0000404196, -4.7758d-6, 3.10801d-7,
     $       -0.0000228947, -4.01957d-7, 3.21794d-8, 4.00914d-7,
     $       2.46048d-8, 5.00427d-9 /
      data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985,
     $       0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  / 

c     Coefficients of total number of molecules in cluster     
      data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323,
     $         -0.0125211, -0.038546, -0.0183749, -0.0619974,
     $         0.0121827, 0.000320184 /
      data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104,
     $         0.00580655, -0.000672316, 0.000172072, 0.000906958,
     $         -0.00010665, -0.0000174762 /      
      data gb2 /0.00102485, 0.000192654, -0.0000154757, 2.65066d-6,
     $         -0.000101674, 2.60288d-6, -3.71766d-7, -9.11728d-7,
     $         2.5346d-7, 6.06504d-8 /
      data gb3 /-2.18646d-6, -6.7043d-7, 5.66661d-8, -3.67471d-9,
     $         2.88195d-7, 1.19416d-8, -5.14875d-10, -5.36796d-9,
     $         -3.63519d-10, -1.42177d-11 /
      data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251,
     $         0.0942243, -0.00851515, 0.00026866, -0.00774234,
     $         0.000610065, 0.000135751 /


      cna=cnai/5.

c     Respect the limits of the parameterization
      if (cna .lt. 1.d4) then ! limit sulf acid conc
         fn = 0.
         rnuc = 1.
c         print*,'cna < 1D4', cna
         goto 10
      endif
      if (cna .gt. 1.0d11) cna=1.0e11 ! limit sulfuric acid conc  
      if (temp .lt. 230.15) temp=230.15 ! limit temp
      if (temp .gt. 305.15) temp=305.15 ! limit temp
      if (rh .lt. 1d-4) rh=1d-4 ! limit rh
      if (rh .gt. 1.) rh=1. ! limit rh
c
c     Mole fraction of sulfuric acid
      xstar=0.740997-0.00266379*temp-0.00349998*log(cna)
     &   +0.0000504022*temp*log(cna)+0.00201048*log(rh)
     &   -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2.
     &   -0.0000179059*temp*(log(rh))**2.
     &   +0.000184403*(log(rh))**3.
     &   -1.50345d-6*temp*(log(rh))**3.
c 
c     Nucleation rate coefficients 
      do i=1, 10
         fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2.
     &        +fb3(i)*temp**3.+fb4(i)/xstar
      enddo
c
c     Nucleation rate (1/cm3-s)
      fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2.
     &    +fb(4)*(log(rh))**3.+fb(5)*log(cna)
     &    +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna)
     &    +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2.
     &    +fb(10)*(log(cna))**3.)

c
c   Cap at 10^6 particles/s, limit for parameterization
      if (fn.gt.1.0d6) then
         fn=1.0d6
      endif
c
c     Coefficients of total number of molecules in cluster 
      do i=1, 10
         gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2.
     &        +gb3(i)*temp**3.+gb4(i)/xstar
      enddo
c     Total number of molecules in cluster
      ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2.
     &    +gb(4)*(log(rh))**3.+gb(5)*log(cna)
     &    +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna)
     &    +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2.
     &    +gb(10)*(log(cna))**3.)

c     cluster radius
      rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

 10   return
      end




!@sum  bl_nucl :  calculates a simple binary nucleation rate of 1 nm
!@+   particles.
!@+       j_1nm = A * [H2SO4]
!@auth   Jeff Pierce, April 2007

      SUBROUTINE bl_nucl(cnai,fn,rnuc)

      IMPLICIT NONE


      real*8,intent(in) :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]

      real*8 fn                 ! nucleation rate [cm-3 s-1]
      real*8 rnuc               ! critical cluster radius [nm]

      real*8 cna                ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 A                  ! prefactor... empirical
      parameter(A=2.0D-6)


      cna=cnai

      fn=A*cna
      rnuc=0.5d0                ! particle diameter of 1 nm

      return
      end SUBROUTINE bl_nucl
            

!@sum  napa_nucl :  calculates the ternary nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of 
!@+     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
!@+    Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
!@+     Research-Atmospheres 107, no. D19 (2002).

!@auth   Jeff Pierce, April 2007

      SUBROUTINE napa_nucl(temp,rh,cnai,nh3ppti,fn,rnuc)

      IMPLICIT NONE

      real*8 cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppti              ! concentration of gas phase ammonia

      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
      real*8 fnl                  ! natural log of nucleation rate
      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppt               ! concentration of gas phase ammonia
      integer i                 ! counter

      data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916,
     $         0.31176, -0.0200738, 0.165536,
     $         6.52645, 3.68024, -0.066514, 0.65874,
     $         0.0599321, -0.732731, 0.728429, 41.3016,
     $         -0.160336, 8.57868, 0.0530167, -2.32736        /

      data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378,
     $         1.64009, -0.752115, 3.26623, -0.258002, -0.204098,
     $         -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736,
     $         -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/
     
      data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052,
     $         0.00939836, -0.00343852, 0.00525813, -0.0489703,
     $         0.00143456, 0.00106259, 0.0122938, -0.00165718,
     $         -0.0362432, 0.000147186, -0.027422, 0.000904383,
     $         -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/
     
      data a3 /-0.000824007, -9.74576d-06, -4.14699d-07, 7.61229d-06,
     $         -1.74927d-05, -1.09753d-05, -8.98038d-06, 0.000146967,
     $         -2.02036d-06, -1.2656d-06, 6.18554d-05, 3.41744d-06,
     $         4.93337d-05, -2.37711d-07, 4.93478d-05, -5.73788d-07, 
     $         8.39522d-08, -6.48365d-07, -2.93564d-05, 8.0459d-08   /


      cna=cnai
      nh3ppt=nh3ppti

c     Napari's parameterization is only valid within limited area
      if ((cna .lt. 1.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
         fn = 0.
         rnuc = 1
         goto 10
      endif  
      if (cna .gt. 1.0d9) cna=1.0d9 ! limit sulfuric acid conc
      if (nh3ppt .gt. 100.) nh3ppt=100. ! limit temp  
      if (temp .lt. 240.) temp=240. ! limit temp
      if (temp .gt. 300.) temp=300. ! limit temp
      if (rh .lt. 0.05) rh=0.05 ! limit rh 
      if (rh .gt. 0.95) rh=0.95 ! limit rh

      do i=1,20
         fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
      enddo

      fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2.
     &  +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh
     &  +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt)
     &  *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna)
     &  +fa(12)*rh
     &  *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh)
     &  *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna)
     &  *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt)
     &  +fa(18)*rh
     &  *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20)
     &  *(log(cna))**2.*(log(nh3ppt))**2.
c
c
      fn=exp(fnl)
c   Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=1.0d6
        fnl=log(fn)
      endif

      rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.

 10   return
      end



!@sum gasdiff: This function returns the diffusion constant of a species in
!@+   air (m2/s).  It uses the method of Fuller, Schettler, and
!@+   Giddings as described in Perry's Handbook for Chemical
!@+   Engineers.
!@auth   Peter Adams, May 2000

      real FUNCTION gasdiff(temp,pres,mw,Sv)

      IMPLICIT NONE

      real*8 temp, pres  !temperature (K) and pressure (Pa) of air
      real mw          !molecular weight (g/mol) of diffusing species
      real Sv          !sum of atomic diffusion volumes of diffusing species

      real mwair, Svair   !same as above, but for air
      real mwf, Svf

      parameter(mwair=28.9, Svair=20.1)

      mwf=sqrt((mw+mwair)/(mw*mwair))
      Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
      gasdiff=1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf

      RETURN
      END

!@sum mnfix  :  examines the mass and number distributions and
!@+   determines if any bins have an average mass outside their normal
!@+   range.  I have seen this happen because the GCM advection seems
!@+   to treat the mass and number tracers inconsistently.  If any bins
!@+   are out of range, I shift some mass and number to a new bin in
!@+   a way that conserves both.

!@auth   Peter Adams, September 2000/Jeff Pierce, August 2007

!@var   Nkx and Mkx are the number and mass distributions


      SUBROUTINE mnfix(Nkx,Mkx)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk

      IMPLICIT NONE

      real*8 Nkx(ibins), Mkx(ibins,icomp)

      integer k,j,L,jj          !counters
      real*8 tot_mass           ! total dry mass
      real*8 xk_hi, xk_lo       ! geometric mean mass of bins that mass is moving to
      real*8 xk_hi1, xk_hi2     ! used as tests to see if mass is way to low or high
      real*8 Neps,Meps
      real*8 tmpvar             ! temporary variable
      real*8 frac_lo_n, frac_lo_m ! fraction of the number and mass going to the lower bin
      real*8 totmass
      parameter(Neps=1.d-20,Meps=1.d0-42)


      ! first check for negative tracers
      do k=1,ibins
         if (Nkx(k).lt.Neps)then
            if (Nkx(k).gt.-1.0d10)then
               Nkx(k)=Neps
               do j=1,icomp
                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
               enddo
            else
               print*,'Negative tracer in mnfix'
               print*,'Nk',Nkx
               print*,'Mkx',Mkx
               call stop_model('- tracer in mnfix',255)
            endif
         endif
         totmass=0.d0
         do j=1,icomp-idiag
            totmass=totmass+Mkx(k,j)
         enddo
         if (totmass.lt.Meps) then
            Nkx(k)=Neps
            do jj=1,icomp
               if (jj.eq.srtso4)then
                  Mkx(k,jj)=Neps*sqrt(xk(k+1)*xk(k))
               else
                  Mkx(k,jj)=0.d0
               endif
            enddo
         endif            
         do j=1,icomp-idiag
            if (Mkx(k,j).lt.0.d0)then
               if (Mkx(k,j).gt.-1.0d0)then
                  Nkx(k)=Neps
                  do jj=1,icomp
                     if (jj.eq.srtso4)then
                        Mkx(k,jj)=Neps*sqrt(xk(k+1)*xk(k))
                     else
                        Mkx(k,jj)=0.d0
                     endif
                  enddo
               else
                  print*,'Negative tracer in mnfix'
                  print*,'Nk',Nkx
                  print*,'Mkx',Mkx
                  call stop_model('- tracer in mnfix',255)
               endif
            endif
         enddo
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
         if (tot_mass/Nkx(k).gt.xk(k+1).or.
     &        tot_mass/Nkx(k).lt.xk(k)) then
c            print*,'out of bounts in mnfix, fixing'
c            print*,k,'AVG',tot_mass/Nkx(k),'lo',xk(k),'hi',xk(k+1)
            ! figure out which bins to redistribute to
            xk_hi1 = sqrt(xk(2)*xk(1))
            xk_hi2 = sqrt(xk(ibins+1)*xk(ibins))
            if (xk_hi1.gt.tot_mass/Nkx(k)) then
               !mass per particle very low
               !conserve mass at expense of number
               tmpvar = Nkx(k)
               Nkx(k)=Neps
               Nkx(1)=Nkx(1)+tot_mass/sqrt(xk(2)*xk(1))
               do j=1,icomp
                  tmpvar=Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
                  Mkx(1,j)=Mkx(1,j)+tmpvar
               enddo
            elseif (xk_hi2.lt.tot_mass/Nkx(k)) then
               !mass per particle very high
               !conserve mass at expernse of number
               Nkx(k)=Neps
              Nkx(ibins)=Nkx(ibins)+tot_mass/sqrt(xk(ibins+1)*xk(ibins))
               do j=1,icomp
                  tmpvar=Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
                  Mkx(ibins,j)=Mkx(ibins,j)+tmpvar
               enddo               
            else ! mass of particle is somewhere within the bins
               L = 2
               xk_hi = sqrt(xk(L+1)*xk(L))
               do while (xk_hi .lt. tot_mass/Nkx(k))
                  L=L+1
                  xk_hi = sqrt(xk(L+1)*xk(L))
               enddo
               xk_lo = sqrt(xk(L)*xk(L-1))
                                ! figure out how much of the number to distribute to the lower bin
               frac_lo_n = (tot_mass - Nkx(k)*xk_hi)/
     &              (Nkx(k)*(xk_lo-xk_hi))
               frac_lo_m = frac_lo_n*Nkx(k)*xk_lo/tot_mass

               tmpvar = Nkx(k)
               Nkx(k) = Neps
               Nkx(L-1) = Nkx(L-1) + frac_lo_n*tmpvar
               Nkx(L) = Nkx(L) + (1-frac_lo_n)*tmpvar
               do j=1,icomp
                  tmpvar = Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j) = 0.d0
                  endif
                  Mkx(L-1,j) = Mkx(L-1,j) + frac_lo_m*tmpvar
                  Mkx(L,j) = Mkx(L,j) + (1-frac_lo_m)*tmpvar
               enddo
               tot_mass=0.0d0
               do j=1,icomp-idiag
                  tot_mass=tot_mass+Mkx(k,j)
               enddo
            endif 
         endif
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
         if (tot_mass/Nkx(k).gt.xk(k+1).or.
     &        tot_mass/Nkx(k).lt.xk(k)) then      
            print*,'ERROR in mnfix'
            print*,'bin',k,'lo',xk(k),'hi',xk(k+1),'avg',tot_mass/Nkx(k)
            print*,'tot_mass',tot_mass
            print*,'Nkx(k)',Nkx(k)
            print*,'Mkx'
            do j=1,icomp
               print*,Mkx(k,j)
            enddo
            call stop_model('out of range in mnfix',255)
         endif
      enddo

      return
      end





C=======================================================================
C
C *** SUBROUTINE getCCN_kappa
C *** WRITTEN BY Yunha Lee
C *** Compute CCN at 0.1, 0.2, 0.3% 
C
C=======================================================================
C
      SUBROUTINE getCCN_kappa(I,J,L)
C
      USE TOMAS_AEROSOL  
      USE TRACER_COM, only : ntm,trm,tr_mm
     &     ,nbins,xk,trpdens,n_AECIL,
     &       n_AOCIL,n_AOCOB,n_ASO4,n_ANACL,n_ADUST,
     &       n_AECOB
      USE TRDIAG_COM, only: taijls=>taijls_loc,ijlt_ccn_01
     &     ,ijlt_ccn_03,ijlt_ccn_02
      USE CONSTANT, only: pi,gasc
      implicit none 
      REAL*8 SURT,DIAM3(NBINS+1),Tvol,DENS(7),A3
c      REAL*8 Tp,BOXM,BOXV
      integer i, j, l, si, n
      integer k,kk,tracnum

      
!@var constants needed for CCN calculation 
      real*8, parameter :: Mv=18.015d-3
      real*8, parameter :: rhow= 1000.d0
!@var temporal CCN 
      real*8, dimension(nsmax) :: ccn_mod 
!@var temporal Sc 
      real*8 Scnew
!@var Sc at each size boundary
      real*8, dimension (nbins) :: Sc,kappa

!@var CCN at 0.1/0.2/0.3%    (#/m3)
c      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: CCN_TOMAS


C initialize CCN_mod
      CCN_mod(:)=0.
C get density 

      dens(1)=trpdens(n_ASO4(1))
      dens(2)=trpdens(n_ANACL(1))
      dens(3)=trpdens(n_AECOB(1))
      dens(4)=trpdens(n_AECIL(1))
      dens(5)=trpdens(n_AOCOB(1))
      dens(6)=trpdens(n_AOCIL(1))
      dens(7)=trpdens(n_ADUST(1))

C surface tension
      SURT   = 0.0761-1.55E-4*(Temp-273.)

      A3 = (4*Mv*SURT/(gasc*Temp*rhow))**3

      DO N=1,NBINS+1 
C Diameter cubed in each size boundary with assuming density =1800 kg/m3. 
        diam3(n)= xk(n)/1800.d0*6.d0/pi
      ENDDO

      DO N=1,NBINS

        Tvol =Mk(N,1)/dens(1)+Mk(N,2)/dens(2)+Mk(N,3)/dens(3)+
     &       Mk(N,4)/dens(4)+Mk(N,5)/dens(5)+Mk(N,6)/dens(6)
     &       +Mk(N,7)/dens(7)   ! total vol. of species 

        kappa(n)=(0.6*Mk(n,1)/dens(1)+1.28*Mk(n,2)/dens(2)
     &       +0.227*Mk(n,6)/dens(6))/Tvol ! average kappa in a bin 
C note that kappa is hard-coded here. 

        Sc(n) = sqrt(4.d0*A3/27.d0/Diam3(n)/kappa(n))
        Sc(n) = min(Sc(n), 10.d0) ! HACK!!! Yunha must fix this. Chances are that particles of diameter 2.55e-9 in mode 1 are just too small for this calculation?
        Sc(n) = exp(Sc(n))
        Sc(n)=(Sc(n)-1.d0)*100.d0

c        print*,'debug_kappa',n,kappa(n),Sc(n)

        if(Sc(n) .lt. 0.) Sc(n)=1.e-6
      ENDDO

      DO N=1,NBINS

C compute CCN at various Smax
 
        DO SI=1,nsmax 
          IF(SC(N) .LE. SMAX(SI)) CCN_mod(SI)=CCN_mod(SI)
     *         +Nk(n)/boxvol ! unit is cm-3 now 
          
C     INTERPOLATION :
          if(N .LT. NBINS)THEN 
            IF(SC(N+1) .lt. SMAX(si) .and. SC(N) .gt. SMAX(SI) ) THEN 
C     compute new Sc (I+1) using the upper limit Dp to determine the activation fraction  
              Scnew = sqrt(4.d0*A3/27.d0/Diam3(n+1)/kappa(n))
              Scnew = min(Scnew, 10.d0) ! HACK!!! Yunha must fix this. Chances are that particles of diameter 2.55e-9 in mode 1 are just too small for this calculation?
              Scnew = exp(Scnew)
              Scnew=(Scnew-1.d0)*100.d0

              if (Sc(n) .ne. Scnew) ! HACK! This is needed to avoid division by zero. Yunha must verify that this is indeed the correct behavior, as it appears to be the case.
     &        CCN_mod(SI)=CCN_mod(SI)+Nk(n)/boxvol*
     &         (1/(dlog(100.+SMAX(SI)))**(2)-1/(dlog(100.+Scnew))**(2))/
     &         (1/(dlog(100.+Sc(n)))**(2)-1/(dlog(100.+Scnew))**(2))
              
            ENDIF    
          ENDIF
          
        ENDDO ! SMAX

      ENDDO ! size bin
C        PRINT*,'interpolate',i,j,l,si,CCN_mod(si) 
C assing each CCN_mod to 3d array 

      CCN_TOMAS(I,J,L,:)=CCN_mod(:) ! ; unit is not cm-3 yet 

        taijls(i,j,l,ijlt_ccn_01)=taijls(i,j,l,ijlt_ccn_01)+
     &   CCN_mod(1)

       taijls(i,j,l,ijlt_ccn_02)=taijls(i,j,l,ijlt_ccn_02)+
     &   CCN_mod(2)

       taijls(i,j,l,ijlt_ccn_03)=taijls(i,j,l,ijlt_ccn_03)+
     &   CCN_mod(3)



      RETURN
      END
