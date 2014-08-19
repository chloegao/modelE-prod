!#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      subroutine archiv(n,nn)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE JulianCalendar_mod, only : jdendofm
      USE MODEL_COM, only : modelEclock,
     *  itime,iyear1,nday,aMON,xlabel,lrunid,monthi,datei
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE HYCOM_SCALARS, only : nstep,time,lp,theta,huge,baclin,onem
     &     ,thref,nhr,g
      USE HYCOM_DIM_GLOB, only : ii1,jj,JDM,kk,isp,ifp,ilp,ntrcr,isu
     &     ,ifu,ilu,isv,ifv,ilv,ii,idm,kdm
      USE HYCOM_ARRAYS_GLOB
#ifdef TRACERS_OceanBiology
      USE obio_com, only : pCO2av,ao_co2fluxav,diag_counter
     .                    ,cexpav,pp2tot_dayav
#ifdef TRACERS_Alkalinity
     .                    ,caexpav
#endif
#endif

      use filemanager, only : findunit
c
      implicit none
      integer i,j,k,l,n,nn,kn
c
      integer no,nop,length,nt
      real factor,vol,tts2,temavg,sst,dpsmo(idm,jdm,kdm)
     .    ,icearea,icevol,icearean,icevoln,iceareas,icevols
      character flnm*40,intvl*4
      character what*16
      real*4 real4(idm,jdm)
     .   ,time4,watcum4,empcum4,thref4,theta4(kdm),unused
      integer*4 length4,idm4,jdm4,kdm4,nstep4
      integer*4 irecl ! specific record lenth, machine dependent
      logical, parameter :: smooth = .false.     ! smooth fields before saving
      data unused/0./
      integer, parameter :: 
     . mon_date(13)=(/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: year, month, dayOfYear, date, hour
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)
      call getdte(Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,amon)
      print *,'stamp',Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,
     .  amon
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.out',xlabel(1:lrunid)
      elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
        write(*,*) 'mismatching archive date in agcm/ogcm=',
     .     (itime+1.)/nday,time
        stop 'mismatching archive date'   
      else
        write(flnm,'(a3,i4.4,2a)') amon,year,'.out',xlabel(1:lrunid)
      endif
c
c     write (lp,*) 'shown below: sea surface height'
c     call zebra(srfhgt,idm,ii1,jj)
      write (lp,*) 'shown below: sea surface temperature'
      call zebra(temp(1,1,1+nn),idm,ii1,jj)
c
      if (date.le.9999) then
        write (intvl,'(i4.4)') date
      else
        stop ' wrong date > 9999'
      endif
c
      no=4096 
      inquire (iolength=irecl)  real4(1,1) ! length of an unformatted real*4  
                                           ! irecl=1 on COMPAQ, irecl=4 on SGI
      length=((irecl*idm*JDM+no+15)/no)*no
c
c
      call findunit(nop)
      open (unit=nop,file=flnm,status='unknown',form='unformatted',
     .      access='direct',recl=length)
      no=1
      length4=length
      idm4=idm
      jdm4=jdm
      kdm4=kdm
      nstep4=nstep
c --- time4: model integration time in day; time: starts as Judian Date 
      time4 = time - mon_date(monthi) - datei + 1
      do k=1,kk
        theta4(k)=theta(k)
      end do
      write(flnm(1:17),'(1x,2(i2.2,a),i4.4,a,i2)')
     .             month,'/',date,'/',year,' hr ',hour+nhr
      write (lp,'(a/9x,2a,f7.1)') 'storing history data in',flnm(1:17)
     .               ,' date=',time4
      write (nop,rec=no) length4,idm4,jdm4,kdm4,nstep4,time4
     .      ,unused,theta4,flnm(1:17)
c
      if (smooth) then
c
        do j=1,jj
         do k=1,kk
          do l=1,isp(j)
           do i=ifp(j,l),ilp(j,l)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
           end do
          end do
        end do
       end do
c
       do k=2,kk
        call psmo1(p(1,1,k),pbot)
       end do
c
      dpsmo=huge
       do j=1,jj
        do k=1,kk
         do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            dpsmo(i,j,k)=p(i,j,k+1)-p(i,j,k)
          end do
         end do
        end do
       end do
c
      end if                            !  smooth
c
      no=no+1
      call r8tor4(ubavg(1,1,n),real4)
      if (smooth) call usmoo4(real4)
      write (nop,rec=no) 'ubavg           ',0,real4
      write (lp,100)     'ubavg           ',0,no
      no=no+1
      call r8tor4(vbavg(1,1,n),real4)
      if (smooth) call vsmoo4(real4)
      write (nop,rec=no) 'vbavg           ',0,real4
      write (lp,100)     'vbavg           ',0,no
      no=no+1
      call r8tor4(srfhgt,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'srfhgt (m)      ',0,real4
      write (lp,100)     'srfhgt (m)      ',0,no
      no=no+1
      call r8tor4(dpmixl,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'mix_dpth        ',0,real4
      write (lp,100)     'mix_dpth        ',0,no
      no=no+1
      call r8tor4(oice,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'icecover(%)     ',0,real4
      write (lp,100)     'icecover(%)     ',0,no
c
      do 75 k=1,kk
      kn=k+nn
      no=no+1
      call r8tor4(u(1,1,kn),real4)
      if (smooth) call usmoo4(real4)
      write (nop,rec=no) 'u               ',k,real4
      write (lp,100)     'u               ',k,no
      no=no+1
      call r8tor4(v(1,1,kn),real4)
      if (smooth) call vsmoo4(real4)
      write (nop,rec=no) 'v               ',k,real4
      write (lp,100)     'v               ',k,no
      no=no+1
      if (smooth) then
        call r8tor4(dpsmo(1,1,k),real4)
      else
        call r8tor4(dp(1,1,kn),real4)
      endif
      write (nop,rec=no) 'dp              ',k,real4
      write (lp,100)     'dp              ',k,no
      no=no+1
      call r8tor4(temp(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'temp            ',k,real4
      write (lp,100)     'temp            ',k,no
      no=no+1
      call r8tor4(saln(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'saln            ',k,real4
      write (lp,100)     'saln            ',k,no
      no=no+1
      call r8tor4(th3d(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'th3d            ',k,real4
      write (lp,100)     'th3d            ',k,no
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- ifort miscomputes 'no' in the following loop:
ccc      do nt=1,ntrcr
ccc        no=no+1
ccc        call r8tor4(tracer(1,1,k,nt),real4)
ccc        if (smooth) call psmoo4(real4)
ccc        write (what,'(a6,i2,4x)') 'tracer',nt
ccc        write (nop,rec=no) what,k,real4
ccc      write (lp,100)       what,k,no
ccc      end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- temporary: store diffusivity as tracer 1
c     tracer(:,:,:,1)=diadff*1.e4
c     write (what,'(a12)') 'diapyc.diffu'
c
c --- code around compiler glitch:
      do nt=1,ntrcr
        call r8tor4(tracer(1,1,k,nt),real4)
        if (smooth) call psmoo4(real4)
        write (what,'(a6,i2,5x)') 'tracer',nt
        write (nop,rec=no+nt) what,k,real4
      write (lp,100)       what,k,no+nt
      end do
      no=no+ntrcr
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 75   continue
c
c --- output time-averaged fields
c
      factor=baclin/(date*SECONDS_PER_DAY)
c
      do 55 j=1,jj
      do 55 l=1,isp(j)
      do 55 i=ifp(j,l),ilp(j,l)
      eminpav(i,j)=eminpav(i,j)*factor
      surflav(i,j)=surflav(i,j)*factor
      salflav(i,j)=salflav(i,j)*factor
      brineav(i,j)=brineav(i,j)*factor
      tauxav(i,j)=tauxav(i,j)*factor
      tauyav(i,j)=tauyav(i,j)*factor
 55   continue
c
#ifdef TRACERS_OceanBiology
      if (diag_counter .ne. 0.d0) then
      do j=1,jj; do l=1,isp(j); do i=ifp(j,l),ilp(j,l)
        ao_co2fluxav(i,j)=ao_co2fluxav(i,j)/diag_counter
        pco2av(i,j)=pco2av(i,j)/diag_counter
        pp2tot_dayav(i,j)=pp2tot_dayav(i,j)/diag_counter
        cexpav(i,j)=cexpav(i,j)/diag_counter
#ifdef TRACERS_Alkalinity
        caexpav(i,j)=caexpav(i,j)/diag_counter
#endif
       enddo; enddo; enddo
        i=243;j=1;
        write(*,'(a,3i5,f10.3,4e12.4)')'2222222222',
     .   nstep,i,j,diag_counter,ao_co2fluxav(i,j),
     .   pco2av(i,j),pp2tot_dayav(i,j),cexpav(i,j)
      endif
#endif

      do 58 k=1,kk
c
      do 59 j=1,jj
      do 591 l=1,isu(j)
      do 591 i=ifu(j,l),ilu(j,l)
      if (dpuav(i,j,k).gt.0) uav(i,j,k)=uav(i,j,k)/dpuav(i,j,k)
 591  uflxav(i,j,k)=uflxav(i,j,k)*baclin
      do 592 l=1,isv(j)
      do 592 i=ifv(j,l),ilv(j,l)
      if (dpvav(i,j,k).gt.0) vav(i,j,k)=vav(i,j,k)/dpvav(i,j,k)
 592  vflxav(i,j,k)=vflxav(i,j,k)*baclin
      do 59 l=1,isp(j)
      do 59 i=ifp(j,l),ilp(j,l)
      if (dpav(i,j,k).gt.0.) then
        temav(i,j,k)=temav(i,j,k)/dpav(i,j,k)
        salav(i,j,k)=salav(i,j,k)/dpav(i,j,k)
        th3av(i,j,k)=th3av(i,j,k)/dpav(i,j,k)
      end if
      dpav(i,j,k)=dpav(i,j,k)*factor
c
      diaflx(i,j,k)=-diaflx(i,j,k)/(2.*onem)
c --- convert diapycnal thickness changes into actual interface fluxes
      if (k.gt.1) diaflx(i,j,k)=diaflx(i,j,k)+diaflx(i,j,k-1)
 59   continue
ccc      write (lp,'(a,i3)') 'shown below: N.Atl. diaflx, bottm of layer',k
ccc      call zebra(diaflx(1,int(.8*jdm),k),idm,idm/3,idm/3)
c
      no=no+1
      call r8tor4(uflxav(1,1,k),real4)
      write (nop,rec=no) '     uflxav_'//intvl,k,real4
      write (lp,100)     '     uflxav_'//intvl,k,no
      no=no+1
      call r8tor4(vflxav(1,1,k),real4)
      write (nop,rec=no) '     vflxav_'//intvl,k,real4
      write (lp,100)     '     vflxav_'//intvl,k,no
      no=no+1
      call r8tor4(diaflx(1,1,k),real4)
      write (nop,rec=no) '     diaflx_'//intvl,k,real4
      write (lp,100)     '     diaflx_'//intvl,k,no
 58   continue
c
      do 56 j=1,jj
      do 561 l=1,isu(j)
      do 561 i=ifu(j,l),ilu(j,l)
 561  ubavav(i,j)=ubavav(i,j)*factor
      do 562 l=1,isv(j)
      do 562 i=ifv(j,l),ilv(j,l)
 562  vbavav(i,j)=vbavav(i,j)*factor
      do 56 l=1,isp(j)
      do 56 i=ifp(j,l),ilp(j,l)
      pbavav(i,j)=pbavav(i,j)*factor
c     sfhtav(i,j)=sfhtav(i,j)*factor+thref*pbavav(i,j)/g  ! meter
      sfhtav(i,j)=sfhtav(i,j)*factor
      dpmxav(i,j)=dpmxav(i,j)*factor
 56   oiceav(i,j)=oiceav(i,j)*factor
c
c     write (lp,'(3a)') 'shown below: ',intvl,'- day SSH average'
c     call zebra(sfhtav,idm,ii1,jj)
      write (lp,'(3a)') 'shown below: ',intvl,'- day SST average'
      call zebra(temav,idm,ii1,jj)
c
      no=no+1
      call r8tor4(ubavav,real4)
      write (nop,rec=no) '     ubavav_'//intvl,0,real4
      write (lp,100)     '     ubavav_'//intvl,0,no
      no=no+1
      call r8tor4(vbavav,real4)
      write (nop,rec=no) '     vbavav_'//intvl,0,real4
      write (lp,100)     '     vbavav_'//intvl,0,no
      no=no+1
      call r8tor4(sfhtav,real4)
      write (nop,rec=no) '     sfhtav_'//intvl,0,real4
      write (lp,100)     '     sfhtav_'//intvl,0,no
      no=no+1
      call r8tor4(dpmxav,real4)
      write (nop,rec=no) '     dpmxav_'//intvl,0,real4
      write (lp,100)     '     dpmxav_'//intvl,0,no
      no=no+1
      call r8tor4(oiceav,real4)
      write (nop,rec=no) '     oiceav_'//intvl,1,real4
      write (lp,100)     '     oiceav_'//intvl,0,no
c
#ifdef TRACERS_OceanBiology
      no=no+1
      call r8tor4(ao_co2fluxav,real4)
      write (nop,rec=no) '  ao_co2flux'//intvl,1,real4
      write (lp,100)     '  ao_co2flux'//intvl,0,no

      no=no+1
      call r8tor4(pco2av,real4)
      write (nop,rec=no) '      pco2av'//intvl,1,real4
      write (lp,100)     '      pco2av'//intvl,0,no

      no=no+1
      call r8tor4(pp2tot_dayav,real4)
      write (nop,rec=no) 'pp2tot_dayav'//intvl,1,real4
      write (lp,100)     'pp2tot_dayav'//intvl,0,no

      no=no+1
      write(*,*)'archyb1: ',cexpav(243,1)
      call r8tor4(cexpav,real4)
      write(*,*)'archyb2: ',real4(243,1)
      write (nop,rec=no) '      cexpav'//intvl,1,real4
      write (lp,100)     '      cexpav'//intvl,0,no

#ifdef TRACERS_Alkalinity
      no=no+1
      call r8tor4(caexpav,real4)
      write (nop,rec=no) '     caexpav'//intvl,1,real4
      write (lp,100)     '     caexpav'//intvl,0,no
#endif
#endif
c
      do 57 k=1,kk
      no=no+1
      call r8tor4(uav(1,1,k),real4)
      write (nop,rec=no) '        uav_'//intvl,k,real4
      write (lp,100)     '        uav_'//intvl,k,no
      no=no+1
      call r8tor4(vav(1,1,k),real4)
      write (nop,rec=no) '        vav_'//intvl,k,real4
      write (lp,100)     '        vav_'//intvl,k,no
      no=no+1
      call r8tor4(dpav(1,1,k),real4)
      write (nop,rec=no) '       dpav_'//intvl,k,real4
      write (lp,100)     '       dpav_'//intvl,k,no
      no=no+1
      call r8tor4(temav(1,1,k),real4)
      write (nop,rec=no) '      temav_'//intvl,k,real4
      write (lp,100)     '      temav_'//intvl,k,no
      no=no+1
      call r8tor4(salav(1,1,k),real4)
      write (nop,rec=no) '      salav_'//intvl,k,real4
      write (lp,100)     '      salav_'//intvl,k,no
      no=no+1
      call r8tor4(th3av(1,1,k),real4)
      write (nop,rec=no) '     th3dav_'//intvl,k,real4
      write (lp,100)     '     th3dav_'//intvl,k,no
 57   continue
c
c --- time-averaged surface fluxes:
      no=no+1
      call r8tor4(eminpav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) '  eminp m/s_'//intvl,0,real4
      write (lp,100)     '  eminp m/s_'//intvl,0,no
      no=no+1
      call r8tor4(surflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' htflx W/m2_'//intvl,0,real4
      write (lp,100)     ' htflx W/m2_'//intvl,0,no
      no=no+1
      call r8tor4(salflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' sflx g/m2s_'//intvl,0,real4
      write (lp,100)     ' sflx g/m2s_'//intvl,0,no
      no=no+1
      call r8tor4(brineav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'brine g/m2s_'//intvl,0,real4
      write (lp,100)     'brine g/m2s_'//intvl,0,no
      no=no+1
      call r8tor4(tauxav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' Tau_x N/m2_'//intvl,0,real4
      write (lp,100)     ' Tau_x N/m2_'//intvl,0,no
      no=no+1
      call r8tor4(tauyav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' Tau_y N/m2_'//intvl,0,real4
      write (lp,100)     ' Tau_y N/m2_'//intvl,0,no
c
      close (unit=nop)
 100  format (9x,a,' (layer',i3,') archived as record',i5)
      write (lp,*) no,' records archived'
c


      do 60 j=1,jj
c
      do 601 i=1,ii
      eminpav(i,j)=0.
      surflav(i,j)=0.
       salflav(i,j)=0.
      brineav(i,j)=0.
      tauxav(i,j)=0.
      tauyav(i,j)=0.
c
#ifdef TRACERS_OceanBiology
        diag_counter =0 
        ao_co2fluxav(i,j)=0.
        pco2av(i,j)=0.
        pp2tot_dayav(i,j)=0.
        cexpav(i,j)=0.
#ifdef TRACERS_Alkalinity
        caexpav(i,j)=0.
#endif
#endif

      ubavav(i,j)=0.
      vbavav(i,j)=0.
      pbavav(i,j)=0.
      dpmxav(i,j)=0.
      sfhtav(i,j)=0.
 601  oiceav(i,j)=0.
c
      do 60 k=1,kk
      do 602 i=1,ii
      uav(i,j,k)=0.
      vav(i,j,k)=0.
      dpuav(i,j,k)=0.
      dpvav(i,j,k)=0.
      dpav (i,j,k)=0.
      temav(i,j,k)=0.
      salav(i,j,k)=0.
      th3av(i,j,k)=0.
      uflxav(i,j,k)=0.
      vflxav(i,j,k)=0.
 602  diaflx(i,j,k)=0.
 60   continue
c
      return
      end
c
      subroutine r8tor4(real8,real4)
c
      USE HYCOM_DIM_GLOB
      implicit none
      integer i,j
c
      real real8(idm,jdm)
      real*4 real4(idm,jdm)
c
      do 1 j=1,jdm
      do 1 i=1,idm
 1    real4(i,j)=real8(i,j)
c
      return
      end
c
c
c> Revision history:
c>
c> June 2001 - corrected sign error in diaflx
c> Feb. 2005 - added multiple tracer capability
c> June 2005 - reordered output fields, added time-av'ged mxlyr & ice temp/th
