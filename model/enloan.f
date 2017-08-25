#include "hycom_mpi_hacks.h"
      subroutine enloan(m,n,mm,nn,k1m,k1n)
c
c --- use as an auxiliary ice model to keep SST above freezing to counteract
c --- possible atmo => ocean interpolation errors.
c --- this is a stripped version of hycom's original enloan routine.
c
      USE SEAICE, only : fsss,tfrez,Ei
c
      USE HYCOM_DIM,only : isp,ifp,ilp,kk,idm,J_0,J_1,J_0H,J_1H
      USE HYCOM_SCALARS, only : thkmin,onem,nstep,delt1,g,spcifh
     &     ,equatn,epsil,brntop,brnbot,itest,jtest
      USE HYCOM_ARRAYS
c
      implicit none
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
c
      real tmelt,tmxl,dpth,heatfx,transf,capcty,dpgt0,old,amount,apermsq
      logical vrbos
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tmelt in winter. return this borrowed energy to the 'energy bank'
c --- in summer as quickly as conditions allow.
c
      do 10 j=J_0,J_1
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
c
      if (dp(i,j,k1n).le.0.) then
        write (*,'(i9,2i5,a)') nstep,i,j,'  zero mxlayr thickness'
        stop '(enloan error)'
      end if
c
c --- if saln < 0 due to excessive freshwater input, bring up salt from below

      if (saln(i,j,k1n).lt.0.) then
        write(*,'(i8,a,2i4,f8.4)')nstep,' warning: neg S =',i,j
     .      ,saln(i,j,k1n)
        transf=-saln(i,j,k1n)*dp(i,j,k1n)
        saln(i,j,k1n)=0.
        k=1
        do while (transf.gt.0.)
          k=k+1
          if (k.gt.kk) then
            write(*,'(i8,a)') nstep,' enloan error: column S < 0'
          end if
          kn=k+nn
          capcty=min(saln(i,j,kn)*dp(i,j,kn),transf)
          dpgt0=max(epsil,dp(i,j,kn))
          saln(i,j,kn)=saln(i,j,kn)-capcty/dpgt0
          transf=transf-capcty
        end do
      endif
c
c --- calculate hypothetical mixed-layer temp
      dpth=max(dp(i,j,k1n),thkmin*onem)
      tmxl=temp(i,j,k1n)+surflx(i,j)*delt1*g/(spcifh*dpth)
      tmelt=tfrez(saln(i,j,k1n),0.)
      old=loan_ice(i,j)

      if (tmxl.lt.tmelt-.3) then		! warm SST by forming new ice
        amount=scp2(i,j)*(tmelt-tmxl)*spcifh*dpth/g		! > 0	  [J]
        loan_ice(i,j)=loan_ice(i,j)+amount
        apermsq=amount*scp2i(i,j)					! J/m^2
        surflx(i,j)=surflx(i,j)+apermsq/delt1				! W/m^2
        if (tmxl.lt.tmelt-.5)
     .   print 100, i,j,' (enloan) -borrow- T =',temp(i,j,k1n),' =>',
     .   temp(i,j,k1n)+apermsq*g/(spcifh*dpth),
     .    'loan =',old*scp2i(i,j),' =>',loan_ice(i,j)*scp2i(i,j)
 100    format (2i5,2(a,f6.2),2x,2(a,es11.3))

      else if (loan_ice(i,j).gt.0.) then	! use warm SST to melt ice
        amount=scp2(i,j)*(tmxl-tmelt)*spcifh*dpth/g		! > 0	  [J]
        amount=max(0.,min(loan_ice(i,j),amount))
     .    * 0.2					! retard melting
        loan_ice(i,j)=loan_ice(i,j)-amount
        apermsq=amount*scp2i(i,j)					! J/m^2
        surflx(i,j)=surflx(i,j)-apermsq/delt1				! W/m^2
        if (apermsq.gt.1.e5)						! J/m^2
     .  print 100, i,j,' (enloan) -return- T =',temp(i,j,k1n),' =>',
     .   temp(i,j,k1n)-apermsq*g/(spcifh*dpth),
     .    'loan =',old,' =>',loan_ice(i,j)
      else
        if (vrbos)
     .  print 100, i,j,' (enloan)       tmxl =',tmxl,' vs',tmelt,
     .    'loan =',old*scp2i(i,j),' =>',loan_ice(i,j)*scp2i(i,j)
      endif
c
      if (vrbos) write (*,103) nstep,i,j,
     .  '  exiting enloan   ice_loan,surflx =',loan_ice(i,j)*scp2i(i,j),
     .  surflx(i,j),'       temp    saln    dens   thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
     .   dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)
 103  format (i8,2i5,a,2f9.1/a/(i3,2f8.3,f8.3,f8.2,f8.1))

 10   continue

! --- nudge loan_ice toward warmer water to avoid multiyear buildup

      call icexport(loan_ice)

      return
      end
c
c> Revision history
c>
c> July 2017 - created stripped version of hycom's enloan
