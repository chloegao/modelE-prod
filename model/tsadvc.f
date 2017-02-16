#include "hycom_mpi_hacks.h"
      subroutine tsadvc(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 1.0 -- cyclic in j
      USE HYCOM_SCALARS, only : wts1,onemm,wts2,delt1,lp,nstep,itest
     &     ,jtest,diagno,temdff,theta
      USE HYCOM_ARRAYS
      USE HYCOM_ARRAYS_GLOB, only: tempGlob => temp, salnGlob => saln

      USE HYCOM_DIM_GLOB,    only: ispGlob  => isp,  ifpGlob  => ifp
      USE HYCOM_DIM_GLOB,    only: ilpGlob  => ilp

      USE HYCOM_DIM, only : ii, jj, kk, idm, jdm, kdm,
     &     J_0, J_1, I_0H, J_0H, J_1H, ogrid,
     &     ip, iu, iv, iq, ifp, ilp, isp, ifq, ilq, isq,
     &     ifu, ilu, isu, ifv, ilv, isv, jchunk

      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT,HALO_UPDATE,NORTH,SOUTH,
     &                          GLOBALMIN,GLOBALMAX

      implicit none
c
      integer i,j,k,l,m,n,mm,nn,km,kn,k1m,k1n,ja,jb
c
      real smin,smax,tmin,tmax,sminn,smaxx,tminn,tmaxx,posdef,flxdiv
     .    ,offset,factor,q,pold,pmid,pnew,snew,tnew,val
      real uflxn(idm,J_0H:J_1H,kdm),vflxn(idm,J_0H:J_1H,kdm),
     .     sign(idm,J_0H:J_1H,kdm),pn(idm,J_0H:J_1H,kdm+1)

      integer kp
      real sigocn,hfharm
      external sigocn,hfharm
c
      do 81 k=1,kk
      km=k+mm
      kn=k+nn
      kp=min(k+1,kk)
c
      smin=999.
      tmin=999.
      smax=-999.
      tmax=-999.
      sminn=999.
      tminn=999.
      smaxx=-999.
      tmaxx=-999.
c
c --- --------------------------------------
c --- advection of thermodynamic variable(s)
c --- --------------------------------------
c
      posdef=100.
c
      call cpy_p_par(dp(I_0H,J_0H,km))
      call cpy_p_par(dp(I_0H,J_0H,kn))
      CALL HALO_UPDATE(ogrid,vflx(:,:,k), FROM=NORTH)

c
      do 49 j=J_0,J_1
      jb= PERIODIC_INDEX(j+1, jj)
      do 49 l=1,isp(j)
      do 49 i=ifp(j,l),ilp(j,l)
c
c --- time smoothing of thermodynamic variable(s) (part 1)
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      temp(i,j,km)=temp(i,j,km)*(wts1*pmid+onemm)+
     .             temp(i,j,kn)* wts2*pold
      saln(i,j,km)=saln(i,j,km)*(wts1*pmid+onemm)+
     .             saln(i,j,kn)* wts2*pold
c
c --- before calling 'advem', make sure (a) mass fluxes are consistent
c --- with layer thickness change, and (b) all fields are positive-definite
      flxdiv=(uflx(i+1,j,k)-uflx(i,j,k)
     .       +vflx(i,jb ,k)-vflx(i,j,k))*delt1*scp2i(i,j)
      util2(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)-flxdiv)
      util1(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)+flxdiv)
      offset=min(0.,util1(i,j),util2(i,j))
      util2(i,j)=util2(i,j)-offset
      util1(i,j)=util1(i,j)-offset
c
      temp(i,j,kn)=temp(i,j,kn)+posdef
      smin=min(smin,saln(i,j,kn))
      smax=max(smax,saln(i,j,kn))
      tmin=min(tmin,temp(i,j,kn))
      tmax=max(tmax,temp(i,j,kn))
 49   continue
c
      call GLOBALMIN(ogrid,smin,val); smin=val;
      call GLOBALMAX(ogrid,smax,val); smax=val;
      call GLOBALMIN(ogrid,tmin,val); tmin=val;
      call GLOBALMAX(ogrid,tmax,val); tmax=val;

      if (tmin.lt.0. .or. smin.lt.0.) then

        ! gather only temp,saln for kn
        call gathPrvTsadvc(temp(1,J_0H,kn),saln(1,J_0H,kn),
     &                     tempGlob(1,1,kn),salnGlob(1,1,kn))

        if (AM_I_ROOT()) then
        do 490 j=1,jj
        do 490 l=1,ispGlob(j)
        do 490 i=ifpGlob(j,l),ilpGlob(j,l)
        if ((tmin.lt.0. and. tempGlob(i,j,kn).eq.tmin) .or.
     .      (smin.lt.0. and. salnGlob(i,j,kn).eq.smin)) then
        write (lp,101) nstep,i,j,k,' neg. temp/saln bfore advem call ',
     .  tempGlob(i,j,kn)-posdef,salnGlob(i,j,kn)
 101    format (i9,' i,j,k =',2i5,i3,a,2f8.2)
        itest=i
        jtest=j
        end if
 490    continue
        end if !AM_I_ROOT

!       call stencl(kn)

      end if
c
      call advem(2,temp(1,J_0H,kn),uflx(1,J_0H,k),
     .           vflx(1,J_0H,k),
     .           scp2(1,J_0H),scp2i(1,J_0H),
     .           delt1,util1(1,J_0H),util2(1,J_0H))
c
      call advem(2,saln(1,J_0H,kn),uflx(1,J_0H,k),
     .           vflx(1,J_0H,k),
     .           scp2(1,J_0H),scp2i(1,J_0H),
     .           delt1,util1(1,J_0H),util2(1,J_0H))
c
      do 46 j=J_0,J_1
      do 46 l=1,isp(j)
      do 46 i=ifp(j,l),ilp(j,l)
      temp(i,j,kn)=temp(i,j,kn)-posdef
      sminn=min(sminn,saln(i,j,kn))
      smaxx=max(smaxx,saln(i,j,kn))
      tminn=min(tminn,temp(i,j,kn))
      tmaxx=max(tmaxx,temp(i,j,kn))
c
c --- time smoothing of thickness field
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      pnew=max(0.,dp(i,j,kn))
      dp(i,j,km)=pmid*wts1+(pold+pnew)*wts2
c --- time smoothing of thermodynamic variable(s) (part 2)
      pmid=max(0.,dp(i,j,km))
      temp(i,j,km)=(temp(i,j,km)+temp(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      saln(i,j,km)=(saln(i,j,km)+saln(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      th3d(i,j,km)=sigocn(temp(i,j,km),saln(i,j,km))
c
c --- build up time integral of mass field variables
      pmid=max(0.,dp(i,j,km))
      dpav (i,j,k)=dpav (i,j,k)+pmid
      temav(i,j,k)=temav(i,j,k)+temp(i,j,km)*pmid
      salav(i,j,k)=salav(i,j,k)+saln(i,j,km)*pmid
      th3av(i,j,k)=th3av(i,j,k)+th3d(i,j,km)*pmid
 46   continue
c
      call GLOBALMIN(ogrid,sminn,val); sminn=val;
      call GLOBALMAX(ogrid,smaxx,val); smaxx=val;
      call GLOBALMIN(ogrid,tminn,val); tminn=val;
      call GLOBALMAX(ogrid,tmaxx,val); tmaxx=val;


      if (tminn+posdef.lt.0. .or. sminn.lt.0.) then

        if (AM_I_ROOT()) then
        call gathPrvTsadvc(temp(1,J_0H,kn),saln(1,J_0H,kn),
     &                     tempGlob(1,1,kn),salnGlob(1,1,kn))

        do 492 j=1,jj
        do 492 l=1,isp(j)
        do 492 i=ifp(j,l),ilp(j,l)
        if (tempGlob(i,j,kn).eq.tminn .or. salnGlob(i,j,kn).eq.sminn)
     .  write (lp,101) nstep,i,j,k,' neg. temp/saln after advem call ',
     .  tempGlob(i,j,kn),salnGlob(i,j,kn)
 492    continue

        endif ! AM_I_ROOT

      end if
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' th,s,dp after advection  '',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      if (diagno) then
        if (AM_I_ROOT()) then
        write (lp,'(i9,i3,'' min/max of s after advection:'',4f7.2)')
     .  nstep,k,sminn,smaxx
        endif ! AM_I_ROOT
      end if
c
c
c --- --------------------------------------
c --- diffusion of thermodynamic variable(s)
c --- --------------------------------------
c
      !call scatter_hycom_arrays

      call cpy_p_par(temp(I_0H,J_0H,kn))
      call cpy_p_par(saln(I_0H,J_0H,kn))
      call cpy_p_par(th3d(I_0H,J_0H,kn))

      CALL HALO_UPDATE(ogrid,saln(:,:,kn), FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,temp(:,:,kn), FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,  dp(:,:,kn), FROM=SOUTH+NORTH)

c
      do 145 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 144 l=1,isu(j)
      do 144 i=ifu(j,l),ilu(j,l)
      factor=scuy(i,j)*2.*hfharm(max(dp(i-1,j,kn),onemm)
     .                          ,max(dp(i  ,j,kn),onemm))
      uflux (i,j)=factor*(temp(i-1,j,kn)-temp(i,j,kn))
 144  uflux2(i,j)=factor*(saln(i-1,j,kn)-saln(i,j,kn))
c
      do 145 l=1,isv(j)
      do 145 i=ifv(j,l),ilv(j,l)
      factor=scvx(i,j)*2.*hfharm(max(dp(i,ja ,kn),onemm)
     .                          ,max(dp(i,j  ,kn),onemm))
      vflux (i,j)=factor*(temp(i,ja ,kn)-temp(i,j,kn))
 145  vflux2(i,j)=factor*(saln(i,ja ,kn)-saln(i,j,kn))
c
      CALL HALO_UPDATE(ogrid,vflux, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,vflux2, FROM=NORTH)

      do 146 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 146 l=1,isp(j)
      do 146 i=ifp(j,l),ilp(j,l)
      factor=-temdff*delt1/(scp2(i,j)*max(dp(i,j,kn),onemm))
      util1(i,j)=(uflux (i+1,j)-uflux (i,j)
     .           +vflux (i,jb )-vflux (i,j))*factor
      util2(i,j)=(uflux2(i+1,j)-uflux2(i,j)
     .           +vflux2(i,jb )-vflux2(i,j))*factor
      temp(i,j,kn)=temp(i,j,kn)+util1(i,j)
      saln(i,j,kn)=saln(i,j,kn)+util2(i,j)
      th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag. write (lp,100) nstep,i,j,k,'t,s,dt,ds,dsigdt,dsigds,cabbl =',
cdiag. temp(i,j,kn),saln(i,j,kn),util1(i,j),util2(i,j),
cdiag. dsigdt(temp(i,j,kn),saln(i,j,kn))*util1(i,j),
cdiag. dsigds(temp(i,j,kn),saln(i,j,kn))*util2(i,j),cabbl(i,j,k)
c
 146  continue
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' t,s,dp after isopyc.mix.'',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      call cpy_p_par(temp(I_0H,J_0H,km))
      call cpy_p_par(temp(I_0H,J_0H,kn))
      call cpy_p_par(saln(I_0H,J_0H,km))
      call cpy_p_par(saln(I_0H,J_0H,kn))
      call cpy_p_par(th3d(I_0H,J_0H,km))
      call cpy_p_par(th3d(I_0H,J_0H,kn))

 81   continue

c
c --- convert mass fluxes to density coord. prior to time integration
c
      call reflux_th(uflx(1,J_0H,1) ,vflx(1,J_0H,1) ,
     .               th3d(1,J_0H,k1m),
     .               p(1,J_0H,1),
     .               uflxn(1,J_0H,1),vflxn(1,J_0H,1),
     .               sign(1,J_0H,1),
     .               pn(1,J_0H,1),theta,kdm,kdm)

c --- activate this loop if -reflux- is   n o t   called
ccc      do k=1,kk
ccc      do j=1,jj
ccc      do i=1,ii
ccc      uflxn(i,j,k)=uflx(i,j,k)
ccc      vflxn(i,j,k)=vflx(i,j,k)
ccc      end do
ccc      end do
ccc      end do
c
      do 155 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 155 k=1,kk
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
      uflxav(i,j,k)=uflxav(i,j,k)+uflxn(i,j,k)		!  uflx time integral
     .   *.5*min(nstep,2)
 154  continue
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
      vflxav(i,j,k)=vflxav(i,j,k)+vflxn(i,j,k)		!  vflx time integral
     .   *.5*min(nstep,2)
 155  continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- convert mass fluxes to pressure coord. prior to time integration
c
      call reflux_pr(uflx(1,J_0H,1) ,vflx(1,J_0H,1) ,
     .               p(1,J_0H,1),
     .               uflxn(1,J_0H,1),vflxn(1,J_0H,1),
     .               pn(1,J_0H,1),kdm,kdm)
c
      do 153 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 153 k=1,kk
c
      do 152 l=1,isu(j)
      do 152 i=ifu(j,l),ilu(j,l)
      ufxavp(i,j,k)=ufxavp(i,j,k)+uflxn(i,j,k)		!  uflx time integral
     .   *.5*min(nstep,2)
 152  continue
c
      do 153 l=1,isv(j)
      do 153 i=ifv(j,l),ilv(j,l)
      vfxavp(i,j,k)=vfxavp(i,j,k)+vflxn(i,j,k)		!  vflx time integral
     .   *.5*min(nstep,2)
 153  continue
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      return
      end
c
c
c  Revision history:
c
c> June 1995 - eliminated setting of salinity in massless layers (loop 46)
c>             (this is now done in mxlayr.f)
c> Aug. 1995 - added array -cabbl- to transmit cabbeling info to -diapfl-
c> Aug. 1995 - omitted t/s/dp time smoothing in case of abrupt mxlayr.thk.change
c> Sep. 1995 - increased temdff if mixed layer occupies >90% of column
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine gathPrvTsadvc(temp,saln,tempGlob,salnGlob)
! gather variables for printing

      USE HYCOM_DIM, only: ogrid, idm, jdm, J_0H, J_1H
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      implicit none
      real temp(idm,J_0H:J_1H), saln(idm,J_0H:J_1H)
      real tempGlob(idm,jdm), salnGlob(idm,jdm)

      call pack_data( ogrid,  temp,     tempGlob )
      call pack_data( ogrid,  saln,     salnGlob )

      end subroutine gathPrvTsadvc
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

