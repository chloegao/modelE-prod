#include "rundeck_opts.h"

!global ?
      subroutine geopar(iniOCEAN)
c
c --- set up model parameters related to geography
c
c --- hycom version 0.9 -- cyclic in j
css   USE GEOM, only : dxyp
c
      USE CONSTANT, only: radian  ! radian=pi/180.
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast
cddd      USE HYCOM_DIM_GLOB, only : ii,jj,kk,ii1,isp,ifp,ilp,ip,isq,ifq,ilq
cddd     &     ,isu,ifu,ilu,jsv,jfv,jlv,ntrcr,jsp,jfp,jlp,msk,iio,jjo
cddd     &     ,iia,jja,idm,jdm, iu,iv,iq
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : lp,pi,area,avgbot,huge,flnmlat,flnmdep
     &   ,flnmbas,ipacn,ipacs,jpac,iatln,iatls,jatl,beropn
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS
      USE HYCOM_CPLER
      USE HYCOM_DYNSI_CPLER
      use filemanager, only : findunit
      use hycom_dim, only : ogrid

      implicit none
      integer i,j,k,l,n,nn,ia,ib,ja,jb,jp,iu1,iu2,iu3
c
      logical, intent(in) :: iniOCEAN
      real realat,sphdis,q,loncor(4),latcor(4)
      integer idim,jdim,length,iz,jz,nt
      character util(idm*jdm+14)*2,preambl(5)*79
      real*4 real4(idm,jdm),lat4(idm,jdm,4),lon4(idm,jdm,4)
c --- 'glufac' = regional viscosity enhancement factor
      real, parameter :: glufac=3., zero=0.

c
c --- read basin depth array
      write (lp,'(2a)') ' reading bathymetry file from ',flnmdep
      call findunit(iu1)
      open (unit=iu1,file=flnmdep,form='unformatted',status='old'
     .     ,action='read')
      read (iu1) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (lp,'(2(a,2i5))') 'depth file dimensions',iz,jz,
     .   '  should be',idm,jdm
        stop '(geopar)'
      end if
      rewind (iu1)
      read (iu1) iz,jz,((real4(i,j),i=1,iz),j=1,jz)
      close (unit=iu1)
      do 9 j=1,jj
      do 9 i=1,ii
 9    depths(i,j)=real4(i,j)
c
c     do 7 j=1,jj
c     do 7 i=1,ii
c7    if (depths(i,j).gt.0.) depths(i,j)=max(botmin.,depths(i,j))
c
c --- reset the Denmark Strait - done in advance
c
c     depths(69,168)=798.
c     depths(69,169)=798.
c     depths(70,167)=798.
c     depths(70,168)=798.
c     depths(71,167)=798.
c
c --- reset the Iceland-Faeroes ridge - done in advance
c     do i=72,74
c     depths(i,175)=798.
c     end do
c
      !write(0,*) "ok ",__FILE__,__LINE__
      if (AM_I_ROOT()) then ! print only on root
        write (lp,*) 'shown below: bottom depth'
        call zebra(depths,idm,ii1,jj)
      endif

      !write(0,*) "ok ",__FILE__,__LINE__
c
c --- determine do-loop limits for u,v,p,q points
      call bigrid(depths)

      !write(0,*) "ok ",__FILE__,__LINE__

!! copy hycom_dim arrays to global grid
      call gather_hycom_dim

      if (AM_I_ROOT()) then
ccc      do 3 i=1,ii1
ccc 3    write (lp,'('' i='',i3,'' jfp,jlp='',7(1x,2i5))') i,
ccc     . (jfp(i,l),jlp(i,l),l=1,jsp(i))
ccc      do 5 j=1,jj
ccc 5    write (lp,'('' j='',i3,'' ifp,ilp='',7(1x,2i5))') j,
ccc     . (ifp(j,l),ilp(j,l),l=1,isp(j))
c
c --- smooth bottom topography (optional)
ccc      call psmoo(depths)
c
c     call prtmsk(ip,depths,util1,idm,ii1,jj,0.,1.,
c    .     'bottom depth (m)')
c
      write (lp,'(2a)') 'read lat/lon from ',flnmlat
      call findunit(iu2)
      open (iu2,file=flnmlat,form='unformatted',status='old')
      read (iu2) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (lp,'(2(a,2i5))') 'error - idm,jdm =',iz,jz,
     .   ' in lat/lon file should be',idm,jdm
        stop '(geopar)'
      end if
      rewind iu2
      read (iu2) iz,jz,lat4,lon4
      close(iu2)
c
      do 8 j=1,jj
      do 8 n=1,4
      do 8 i=1,ii
      latij(i,j,n)=lat4(i,j,n)
 8    lonij(i,j,n)=lon4(i,j,n)
c
c     write (lp,*) 'shown below: latitude of vorticity points'
c     call zebra(latij(1,1,4),idm,ii,jj)
c     write (lp,*) 'shown below: longitude of vorticity points'
c     call zebra(lonij(1,1,4),idm,ii,jj)
c
c --- define coriolis parameter and grid size
      do 56 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 56 i=1,ii
c
      corio(i,j)=sin(latij(i,j,4)*radian)*4.*pi/86164.        !  86400 * 365 / 366
c
      scpy(i,j)=sphdis(latij(i,j ,2),lonij(i,j ,2),
     .                 latij(i,jb,2),lonij(i,jb,2))
c
      if (i.lt.ii) then
      scpx(i,j)=sphdis(latij(i  ,j,1),lonij(i  ,j,1),
     .                 latij(i+1,j,1),lonij(i+1,j,1))
#ifdef CUBED_SPHERE
c Define cell areas using great-circle polygons.
c In future, will be done for all cases.
      loncor(1:2) = (/ lonij(i:i+1:+1,j ,4) /) ! List the 4 lons/lats
      latcor(1:2) = (/ latij(i:i+1:+1,j ,4) /) ! of cell corners in
      loncor(3:4) = (/ lonij(i+1:i:-1,jb,4) /) ! counterclockwise order
      latcor(3:4) = (/ latij(i+1:i:-1,jb,4) /) ! for the area calculation
      n = 4
      do nn=1,4  ! check whether this cell is a triangle
        if(abs(loncor(nn)-loncor(1+mod(nn,4))).lt.1d-3 .and.
     &     abs(latcor(nn)-latcor(1+mod(nn,4))).lt.1d-3) then
          loncor = cshift(loncor,nn)
          latcor = cshift(latcor,nn)
          n = 3
          exit
        endif
      enddo
      call gc_polyarea(loncor,latcor,n,scp2(i,j))
#else
      scp2(i,j)=scpx(i,j)*scpy(i,j)
#endif
      scp2i(i,j)=1./scp2(i,j)
      end if
c
      scuy(i,j)=sphdis(latij(i,j ,4),lonij(i,j ,4),
     .                 latij(i,jb,4),lonij(i,jb,4))
c
      if (i.gt.1) then
      scux(i,j)=sphdis(latij(i  ,j,3),lonij(i  ,j,3),
     .                 latij(i-1,j,3),lonij(i-1,j,3))
      scu2(i,j)=scux(i,j)*scuy(i,j)
      scuxi(i,j)=1./scux(i,j)
      end if
c
      scvy(i,j)=sphdis(latij(i,j ,3),lonij(i,j ,3),
     .                 latij(i,ja,3),lonij(i,ja,3))
c
      if (i.lt.ii) then
      scvx(i,j)=sphdis(latij(i  ,j,4),lonij(i  ,j,4),
     .                 latij(i+1,j,4),lonij(i+1,j,4))
      scv2(i,j)=scvx(i,j)*scvy(i,j)
      scvyi(i,j)=1./scvy(i,j)
      end if
c
      scqy(i,j)=sphdis(latij(i,j ,1),lonij(i,j ,1),
     .                 latij(i,ja,1),lonij(i,ja,1))
c
      if (i.gt.1) then
      scqx(i,j)=sphdis(latij(i  ,j,2),lonij(i  ,j,2),
     .                 latij(i-1,j,2),lonij(i-1,j,2))
      scq2(i,j)=scqx(i,j)*scqy(i,j)
      scq2i(i,j)=1./scq2(i,j)
      end if
c
 56   continue
c
      if (beropn .and. scu2(ipacs,jpac).ne.scu2(iatln,jatl))
     .  write(*,'(a,6f13.5)') ' chk WRONG scu2'
     . ,scu2(ipacs,jpac),scu2(iatln,jatl)
c
      if ( latij(ipacs,jpac,3).ne.latij(iatls,jatl,3)
     . .or.latij(ipacn,jpac,3).ne.latij(iatln,jatl,3)
     . .or.lonij(ipacs,jpac,3).ne.lonij(iatls,jatl,3)
     . .or.lonij(ipacn,jpac,3).ne.lonij(iatln,jatl,3))
     .  write(*,'(a,8f9.2)') ' chk WRONG lat/lon '
     .,latij(ipacs,jpac,3),latij(iatls,jatl,3)
     .,latij(ipacn,jpac,3),latij(iatln,jatl,3)
     .,lonij(ipacs,jpac,3),lonij(iatls,jatl,3)
     .,lonij(ipacn,jpac,3),lonij(iatln,jatl,3)
c
      write(*,'(a,3f8.2)') 'lat/lon/depth of Bering Strait:'
     . ,latij(ipacs,jpac,3),lonij(ipacs,jpac,3),depths(ipacs,jpac)
c     write (lp,'('' shown below: coriolis parameter'')')
c     call zebra(corio,idm,ii,jj)
c     write (lp,'('' shown below: grid cell size'')')
c     call zebra(scp2,idm,ii,jj)
c
      area=0.
      avgbot=0.
c
      do 57 j=1,jj
      do 57 l=1,isp(j)
      do 57 i=ifp(j,l),ilp(j,l)
      avgbot=avgbot+depths(i,j)*scp2(i,j)
 57   area=area+scp2(i,j)
      avgbot=avgbot/area
      write (lp,100) avgbot,area
 100  format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .       -12p,f9.1)
c
c --- initialize some arrays
c
      ! uncommented by IA
      !if (nstep0.eq.0) then
      if (iniOCEAN) then
      write (lp,*) 'laying out arrays in memory ...'
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,1)=huge
      if (ip(i,j).eq.1) p(i,j,1)=zero
      pu(i,j,1)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,1)=huge
      ubavg(i,j,2)=huge
      ubavg(i,j,3)=huge
      vbavg(i,j,1)=huge
      vbavg(i,j,2)=huge
      vbavg(i,j,3)=huge
      utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=zero
      dpmixl(i,j,:)= 1.0  ! TNL: avoid NaN on the first step
      oice(i,j)=zero
      taux(i,j)=zero
      tauy(i,j)=zero
      oflxa2o(i,j)=zero
      osalt(i,j)=zero
      oemnp(i,j)=zero
      ustar(i,j)=zero
      sswflx(i,j)=zero
c
      ubavav(i,j)=zero
      vbavav(i,j)=zero
      pbavav(i,j)=zero
      sfhtav(i,j)=zero
      dpmxav(i,j)=zero
      oiceav(i,j)=zero
      eminpav(i,j)=zero
      surflav(i,j)=zero
      salflav(i,j)=zero
      brineav(i,j)=zero
      tauxav(i,j)=zero
      tauyav(i,j)=zero
c
      do 209 k=1,kk
      u  (i,j,k   )=huge
      u  (i,j,k+kk)=huge
      v  (i,j,k   )=huge
      v  (i,j,k+kk)=huge
      uflx(i,j,k)=huge
      vflx(i,j,k)=huge
      ufxcum(i,j,k)=huge
      vfxcum(i,j,k)=huge
      dpinit(i,j,k)=huge
      dpold (i,j,k)=huge
      dp (i,j,k   )=huge
      dp (i,j,k+kk)=huge
      dpu(i,j,k   )=huge
      dpu(i,j,k+kk)=huge
      dpv(i,j,k   )=huge
      dpv(i,j,k+kk)=huge
      p (i,j,k+1)=huge
      pu(i,j,k+1)=huge
      pv(i,j,k+1)=huge
c
      th3d(i,j,k+kk)=huge
      th3d(i,j,k   )=huge
      thstar(i,j,k)=huge
      do nt=1,ntrcr
        tracer(i,j,k,nt)=zero
      end do
      uav(i,j,k)=zero
      vav(i,j,k)=zero
      dpuav(i,j,k)=zero
      dpvav(i,j,k)=zero
      dpav (i,j,k)=zero
      temav(i,j,k)=zero
      salav(i,j,k)=zero
      th3av(i,j,k)=zero
      uflxav(i,j,k)=zero
      vflxav(i,j,k)=zero
      diaflx(i,j,k)=zero
 209  continue
c
      do 210 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 210 l=1,isq(j)
      do 210 i=ifq(j,l),ilq(j,l)
      pbot(i  ,j  )=0.
      pbot(i-1,j  )=0.
      pbot(i  ,ja )=0.
      pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      p(i-1,j  ,1)=0.
      p(i  ,ja ,1)=0.
      p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,k   )=0.
      dp(i  ,j  ,k+kk)=0.
      dp(i-1,j  ,k   )=0.
      dp(i-1,j  ,k+kk)=0.
      dp(i  ,ja ,k   )=0.
      dp(i  ,ja ,k+kk)=0.
      dp(i-1,ja ,k   )=0.
 210  dp(i-1,ja ,k+kk)=0.
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
      do 156 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 156 l=1,isu(j)
      do 156 i=ifu(j,l),ilu(j,l)
      pu(i,j,1)=0.
c
      depthu(i,ja)=0.
      utotn (i,ja)=0.
      pgfx  (i,ja)=0.
c
      depthu(i,jb)=0.
      utotn (i,jb)=0.
      pgfx  (i,jb)=0.
c
      do 156 k=1,kk
      dpu(i,ja,k   )=0.
      dpu(i,ja,k+kk)=0.
c
      dpu(i,jb,k   )=0.
      dpu(i,jb,k+kk)=0.
 156  continue
c
      do 158 j=1,jj
      do 158 l=1,isp(j)
      do 158 i=ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      do 158 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
      uflx(i,j,k)=0.
      ufxcum(i,j,k)=0.
      u(i,j,k   )=0.
 158  u(i,j,k+kk)=0.
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      do 166 i=1,ii1
      ia=mod(i-2+ii,ii)+1
      ib=i+1
      do 166 l=1,jsv(i)
      do 166 j=jfv(i,l),jlv(i,l)
      pv(i,j,1)=0.
c
      depthv(ia,j)=0.
      vtotn (ia,j)=0.
      pgfy  (ia,j)=0.
c
      depthv(ib,j)=0.
      vtotn (ib,j)=0.
      pgfy  (ib,j)=0.
c
      do 166 k=1,kk
      dpv(ia,j,k   )=0.
      dpv(ia,j,k+kk)=0.
c
      dpv(ib,j,k   )=0.
      dpv(ib,j,k+kk)=0.
 166  continue
c
      do 168 i=1,ii1
      do 168 l=1,jsp(i)
      do 168 jp=jfp(i,l),jlp(i,l)+1
      j=mod(jp-1+jj,jj)+1
      depthv(i,j)=0.
      vtotn (i,j)=0.
      pgfy  (i,j)=0.
      vbavg(i,j,1)=0.
      vbavg(i,j,2)=0.
      vbavg(i,j,3)=0.
      vtotm (i,j)=0.
      vflux (i,j)=0.
      vflux2(i,j)=0.
      vflux3(i,j)=0.
      via(i,j)=0.
      vib(i,j)=0.
c
      do 168 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
      vflx(i,j,k)=0.
      vfxcum(i,j,k)=0.
      v(i,j,k   )=0.
 168  v(i,j,k+kk)=0.
      write (lp,*) '... array layout completed'
      ! uncommented by IA
      endif                    ! end of nstep=0
c
c --- set 'glue' to values > 1 in regions where extra viscosity is needed
c
      do 154 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 154 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ib=mod(i     ,ii)+1
      glue(i,j)=1.
c --- add glue in coastal areas
      if (depths(i,j).lt.200.
     .    .or. depths(ia,j).lt.200. .or. depths(i,ja).lt.200.
     .    .or. depths(ib,j).lt.200. .or. depths(i,jb).lt.200.)
     .    glue(i,j)=max(glue(i,j),glufac)
      if (depths(i,j).lt.100.
     .    .or. depths(ia,j).lt.100. .or. depths(i,ja).lt.100.
     .    .or. depths(ib,j).lt.100. .or. depths(i,jb).lt.100.)
     .    glue(i,j)=max(glue(i,j),2.*glufac)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- add glue to mediterranean:
#ifdef HYCOM2deg
        if (i.ge.  91 .and. i.le.  98 .and. j.le.  25)
     .        glue(i,j)=glufac
#endif
#ifdef HYCOM1deg
        if ((i.ge. 180 .and. i.le. 198 .and. j.le.  37)
     . .or. (i .ge.188 .and. i.le. 191 .and. j.ge. 356))
     .        glue(i,j)=glufac
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 154  continue
c
c --- 1:9 represent NAT, SAT, NIN, SIN, NPA, SPA, ARC, SO, MED
      call findunit(iu3)
      open (iu3,file=flnmbas,form='formatted',status='old')
#ifdef HYCOM2deg
        do n=1,2
        read(iu3,*)
        read(iu3,'(90i1)') ((msk(i,j),j=(n-1)*jj/2+1,n*jj/2),i=1,ii)
        enddo
#endif
#ifdef HYCOM1deg
        do n=1,3
        read(iu3,*)
        read(iu3,'(4x,120i1)') ((msk(i,j),j=(n-1)*jj/3+1,n*jj/3),i=1,ii)
        enddo
#endif
      close(iu3)
c
      do i=1,ii
      do j=1,jj
      ijlist(i,j)=1000*i+j
      enddo
      enddo
c
      write(*,'(a,20i12)') 'ijlist ',((ijlist(i,j),i=30,32),j=4,5)
      endif ! AM_I_ROOT
c
      call cpl_wgt                      ! read in weights for coupler
      call init_hycom_dynsi_cpler
c
      call broadcast(ogrid,area)

      return
      end
c
      function sphdis(x1,y1,x2,y2)
c --- dist.(m) between 2 points on sphere, lat/lon (x1,y1) and lat/lon (x2,y2)
      USE CONSTANT, only: radius,radian    ! radian = pi/180.
      implicit none
      real x1,y1,x2,y2,sphdis,ang
      real x1_rad,x2_rad,ang_rad
c
      ang=mod(y2-y1+540.,360.)-180.
      x1_rad=(90-x1)*radian    ! convert degrees to radian
      x2_rad=(90-x2)*radian    ! convert degrees to radian
      ang_rad=ang*radian       ! convert degrees to radian

      sphdis=radius*acos(min(1.,cos(x1_rad)*cos(x2_rad)
     .                         +sin(x1_rad)*sin(x2_rad)*cos(ang_rad)))
      if (sphdis.eq.0.)  sphdis=radius*
     .  sqrt((x2-x1)**2+(ang*cos(.5*(x1_rad+x2_rad)))**2)*radian
cdiag if (sphdis.eq.0.) write (*,'(a,2f8.3,2x,2f8.3)')
cdiag.  'warning - zero distance between lat/lon points',x1,y1,x2,y2
      sphdis=max(sphdis,1.)
      return
      end

      subroutine gc_polyarea(lon,lat,n,area)
!@sum gc_polyarea calculates the area of a polygon on a sphere
!@+   whose edges are great circles
!@+   M. Kelley
      USE CONSTANT, only: radius,twopi
      implicit none
      integer :: n
      real*8, dimension(n) :: lon,lat ! input: lon,lat in degrees
      real*8 :: area
      real*8, dimension(3) :: vi,vip1,vn,cri,crip1,crn
      integer :: i
      real*8 :: twopibyn,cc

      twopibyn = twopi/n
      vn = v3d(lon(n),lat(n))
      vi = v3d(lon(1),lat(1))
      crn = cross3d(vn,vi)
      cri = crn
      area = 0.
      do i=1,n-1
        vip1 = v3d(lon(i+1),lat(i+1))
        crip1 = cross3d(vi,vip1)
        cc = sum(vi*cross3d(cri,crip1))
        area = area + (twopibyn-atan2(cc,sum(cri*crip1)))
        vi = vip1
        cri = crip1
      enddo
      cc = sum(vi*cross3d(cri,crn))
      area = area + (twopibyn-atan2(cc,sum(cri*crn)))
      area = area*radius*radius
      return
      contains
      function v3d(lon,lat)
      USE CONSTANT, only: radian   ! radian = pi/180.
      real*8 :: lon,lat,v3d(3)
      real*8 :: rlon,rlat  ! convert to radians
      rlon=lon*radian
      rlat=lat*radian
      v3d(1:2) = cos(rlat)*(/cos(rlon),sin(rlon)/); v3d(3) = sin(rlat)
      end function v3d
      function cross3d(v1,v2)
      real*8, dimension(3) :: v1,v2,cross3d
      cross3d = cshift(v1,1)*cshift(v2,-1)-cshift(v1,-1)*cshift(v2,1)
      end function cross3d
      end subroutine gc_polyarea
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - changed loop 56 from i/j to j/i to improve memory layout
c> Oct. 2000 - added code to compute 'glue'
c> Apr. 2001 - eliminated stmt_funcs.h
c> Dec. 2001 - added clause to assure p=0 in single-width channels (loop 209)
