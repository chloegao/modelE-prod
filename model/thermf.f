      subroutine thermf(m,n,mm,nn,k1m,k1n
     .   ,sss_restore_dt,sss_restore_dtice)
c
      USE SEAICE, only : fsss
c
c --- hycom version 0.9
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : baclin,thref,watcum,empcum,nstep,nstep0
     &     ,diagno,lp,area,spcifh,avgbot,g,onem,slfcum,delt1,itest,jtest
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, GLOBALSUM
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
c
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
c
      real thknss,radfl,radflw,radfli,vpmx,prcp,prcpw,prcpi,
     .     evap,evapw,evapi,exchng,target,old,
     .     rmean,tmean,smean,vmean,boxvol,fxbias,
     &     slfcol(J_0H:J_1H),watcol(J_0H:J_1H),empcol(J_0H:J_1H),
     &     rhocol(J_0H:J_1H),temcol(J_0H:J_1H),salcol(J_0H:J_1H),
     &     sf1col(J_0H:J_1H),sf2col(J_0H:J_1H),clpcol(J_0H:J_1H),
     &     numcol(J_0H:J_1H),fxbiasj(J_0H:J_1H)
      integer iprime,ktop
      real qsatur,totl,eptt,salrlx,sf1cum,sf2cum,bias,clpcum,numcum
      external qsatur
      data ktop/3/
ccc      data ktop/2/                        !  normally set to 3
css   data salrlx/0.3215e-7/          !  1/(1 yr)
      data salrlx/0.6430e-8/          !  1/(5 yr)
      real*8 :: sss_restore_dt,sss_restore_dtice
      real :: piston
      logical sss_relax
#ifdef STANDALONE_OCEAN
      data sss_relax/.true./
      real*8, dimension(:,:), pointer :: sssobs,rsiobs
#else
      data sss_relax/.false./
#endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: weak salinity restoring via corrective surface flux
c
ccc      if (mod(time+.0001,30.).lt..0002) then            !  once a month
ccc        totl=0.
ccc        eptt=0.
ccc        do 8 j=1,jj
ccc        do 8 k=1,kk
ccc        do 8 l=1,isp(j)
ccc        do 8 i=ifp(j,l),ilp(j,l)
ccc        eptt=eptt+oemnp(i,j)*scp2(i,j)
ccc 8      totl=totl+saln(i,j,k)*dp(i,j,k)*scp2(i,j)
ccc        totl=totl/g                                     !  10^-3 kg
cccc
cccc --- initialize total salt content
ccc        if (saltot.le.0.) saltot=totl
cccc --- determine corrective flux using multi-year relaxation time scale
ccc        pcpcor=(totl-saltot)*salrlx                     !  10^-3 kg/sec
ccc     .   *thref/35.                                     !  => m^3/sec
ccc     .   /eptt                                          !  => rel. units
ccc        write (lp,'(i9,a,f9.5,a,f7.2)') nstep,
ccc     .  '  overall salt gain (psu):',(totl-saltot)*thref/(area*avgbot),
ccc     .  '  => precip bias(%):',100.*pcpcor
ccc      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do 66 j=J_0,J_1
      do 66 k=1,kk
      kn=k+nn
      do 66 l=1,isp(j)
      do 66 i=ifp(j,l),ilp(j,l)
      if (glue(i,j).gt.1. .and. saln(i,j,kn).gt.40.)                !  Med fudge
     .  saln(i,j,kn)=saln(i,j,kn)+(40.-saln(i,j,kn))*baclin*3.e-8
 66   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
c --- for conservation reasons, (E-P)-induced global virtual salt flux must
c --- be proportional to global E-P. this requires a global corrrection.
c
      do 82 j=J_0,J_1
      sf2col(j)=0.
      do 82 l=1,isp(j)
      do 82 i=ifp(j,l),ilp(j,l)
 82   sf2col(j)=sf2col(j)+(saln(i,j,k1n)-34.7)*oemnp(i,j)/thref
     .                                                   *scp2(i,j)
c
      call GLOBALSUM(ogrid,sf2col,sf2cum, all=.true.)
      sf2cum=sf2cum/area
c
      do 85 j=J_0,J_1
c
      watcol(j)=0.
      empcol(j)=0.
      slfcol(j)=0.
      sf1col(j)=0.
      rhocol(j)=0.
      temcol(j)=0.
      salcol(j)=0.
      clpcol(j)=0.
      numcol(j)=0.
c
      do 85 l=1,isp(j)
c
      do 85 i=ifp(j,l),ilp(j,l)
c     ustar(i,j)=sqrt(sqrt(taux(i,j)**2+tauy(i,j)**2)*thref)
c    .                             *(1.-oice(i,j))
c
      surflx(i,j)=oflxa2o(i,j)
c
c --- oemnp = evaporation minus precipitation over open water (m/sec)
c --- salflx = salt flux (10^-3 kg/m^2/sec) in +p direction, salt from SI added
c     salflx(i,j)=saln(i,j,k1n)*(osalt(i,j)*fsss-oemnp(i,j)/thref)
      salflx(i,j)=-saln(i,j,k1n)*oemnp(i,j)/thref+osalt(i,j)*1.e3
     .            +sf2cum			! global E-P constraint
      sf1col(j)=sf1col(j)+salflx(i,j)*scp2(i,j)
c
c --- clip neg.(outgoing) salt flux to prevent S < 0 in top layer
      old=salflx(i,j)
      salflx(i,j)=max(salflx(i,j),-saln(i,j,k1n)*dp(i,j,k1n)/(g*delt1)
     .   * .7)
      if (old.lt.salflx(i,j)) then			! diagnostic use
        clpcol(j)=clpcol(j)+(salflx(i,j)-old)*scp2(i,j)
        numcol(j)=numcol(j)+1.
        print '(2i5,2(a,es14.6))',i,j,' salflx',old,' =>',salflx(i,j)
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if (abs(i-itest).le.2.and.abs(j-jtest).le.2) 
     .  write (*,100) nstep,i,j,
     .  '  surflx    salflx     salin     oemnp     osalt      oice',
     .  surflx(i,j),salflx(i,j),saln(i,j,k1n),oemnp(i,j),osalt(i,j),
     .  oice(i,j)
100   format(i9,2i5,a/18x,6es10.3/(19x,a/18x,6es10.3))
c
      watcol(j)=watcol(j)+surflx(i,j)*scp2(i,j)
      empcol(j)=empcol(j)+oemnp(i,j)*scp2(i,j)
      slfcol(j)=slfcol(j)+salflx(i,j)*scp2(i,j)
      rhocol(j)=rhocol(j)+th3d(i,j,k1n)*scp2(i,j)
      temcol(j)=temcol(j)+temp(i,j,k1n)*scp2(i,j)
      salcol(j)=salcol(j)+saln(i,j,k1n)*scp2(i,j)
 85   continue
c
      call GLOBALSUM(ogrid,watcol,watcum, all=.true.)
      call GLOBALSUM(ogrid,empcol,empcum, all=.true.)
      call GLOBALSUM(ogrid,slfcol,slfcum, all=.true.)
      call GLOBALSUM(ogrid,sf1col,sf1cum, all=.true.)
      call GLOBALSUM(ogrid,clpcol,clpcum, all=.true.)
      call GLOBALSUM(ogrid,numcol,numcum, all=.true.)
      if( AM_I_ROOT() )
     .  print '(i9,a,2es11.3)',nstep,' avg.salt flux, E-P adjustment:',
     .   sf1cum/area,sf2cum
c
c --- correct salt flux globally for local clipping done to prevent S < 0
      bias=(slfcum-sf1cum)/area
      if (numcum.gt.0.) then

        do 83 j=J_0,J_1
        sf2col(j)=0.
        do 83 l=1,isp(j)
        do 83 i=ifp(j,l),ilp(j,l)
        salflx(i,j)=salflx(i,j)-bias
        sf2col(j)=sf2col(j)+salflx(i,j)*scp2(i,j)
 83     continue
c --- optional, diagnostic use only:
        call GLOBALSUM(ogrid,sf2col,sf2cum, all=.true.)
        if( AM_I_ROOT() ) then
          print '(a,3es15.7)','orig/clipped/restored salt flux:',
     .      sf1cum,slfcum,sf2cum
          print '(a,i4,a,es11.3)','salt flux clipped at',nint(numcum),
     .      ' points => required global flux corr:',bias
        end if
      end if

#ifdef STANDALONE_OCEAN
      if (sss_relax) then
        sssobs => atmocn%sssobs  ! units:  psu/1000
        rsiobs => atmocn%rsiobs  ! sea ice fraction [0-1]
c --- surface salinity restoration
        do 84 j=J_0,J_1
        fxbiasj(j)=0.
        do 84 l=1,isp(j)
        do 84 i=ifp(j,l),ilp(j,l)
          piston=((1.-rsiobs(i,j))/sss_restore_dt     ! open water
     .          +     rsiobs(i,j) /sss_restore_dtice) ! under ice
     .          *12./SECONDS_PER_DAY                  ! 12m depth scale
          old=salflx(i,j)
          fxbias=(sssobs(i,j)*1000.-saln(i,j,k1n))*piston/thref
          salflx(i,j)=salflx(i,j)+fxbias
          fxbiasj(j)=fxbiasj(j)+fxbias*scp2(i,j)
          if( AM_I_ROOT() ) then
            if (abs(i-itest).le.2 .and. abs(j-jtest).le.2) 
     .   write (*,'(i9,2i5,a,2f7.2/19x,a,2es10.3,a,f6.2)') nstep,i,j,
     . '  actual & reference salinity:',saln(i,j,k1n),sssobs(i,j)*1.e3,
     . '  salflx before & after relax:',old,salflx(i,j)
          end if
 84     continue
      end if
#endif

      if (nstep.eq.nstep0+1 .or. diagno) then

      call GLOBALSUM(ogrid,rhocol,rmean, all=.true.)
      call GLOBALSUM(ogrid,temcol,tmean, all=.true.)
      call GLOBALSUM(ogrid,salcol,smean, all=.true.)
c
      if( AM_I_ROOT() ) then
      write (lp,'(i9,''energy residual (w/m^2)'',f9.2)') nstep,
     .  watcum/(area*(nstep-nstep0))
      write (lp,'(9x,''resulting temp drift (deg/century):'',f9.3)')
     .  watcum*thref/(spcifh*avgbot*area*(nstep-nstep0)) *
     &  36500.*SECONDS_PER_DAY
      write (lp,'(i9,''e - p residual (mm/year)'',f9.2)') nstep,
     .  empcum/(area*(nstep-nstep0))*SECONDS_PER_DAY*365000.
      write (lp,'(9x,''saln drift resulting from e-p (psu/century):''
     .                                                     ,f9.3)')
     .  -empcum*35./(avgbot*area*(nstep-nstep0)) *36500.*SECONDS_PER_DAY
css   write (lp,'(i9,''salt residual (T/year)'',3f9.2)') nstep,
css  .  slfcum*365.*SECONDS_PER_DAY*g/onem
      write (lp,'(7x,''saln drift resulting from salfl (psu/century):''
     .                                                     ,f9.3)')
     .  slfcum*36500.*SECONDS_PER_DAY*g/(avgbot*area*onem)
      write (lp,'(i9,a,3f9.3)') nstep,' mean surf. sig,temp,saln:',
     .    rmean/area,tmean/area,smean/area
      end if ! AM_I_ROOT
c
      rmean=0.
      smean=0.
      tmean=0.
      vmean=0.
      ! reusing these arrays for next sums 
      rhocol=0.; temcol=0.; salcol=0.; watcol=0.;
      do 84 j=J_0,J_1
      do 84 k=kk,1,-1
      kn=k+nn
      if (nstep.eq.nstep0+1) kn=k+mm
      do 84 l=1,isp(j)
      do 84 i=ifp(j,l),ilp(j,l)
      boxvol=dp(i,j,kn)*scp2(i,j)
      !--changing following to facilitate GLOBALSUM with bitwise reprocibility
      ! rmean=rmean+th3d(i,j,kn)*boxvol
      ! smean=smean+saln(i,j,kn)*boxvol
      ! tmean=tmean+temp(i,j,kn)*boxvol
      ! vmean=vmean+boxvol
      !-----------------------------------------------
      rhocol(j)=rhocol(j)+th3d(i,j,kn)*boxvol
      salcol(j)=salcol(j)+saln(i,j,kn)*boxvol
      temcol(j)=temcol(j)+temp(i,j,kn)*boxvol
      watcol(j)=watcol(j)+boxvol
 84   continue
      call GLOBALSUM(ogrid,rhocol,rmean, all=.true.) 
      call GLOBALSUM(ogrid,salcol,smean, all=.true.) 
      call GLOBALSUM(ogrid,temcol,tmean, all=.true.) 
      call GLOBALSUM(ogrid,watcol,vmean, all=.true.) 
c
      if( AM_I_ROOT() )
     .    write (lp,'(i9,a,3f9.3)') nstep,' mean basin sig,temp,saln:',
     .    rmean/vmean,tmean/vmean,smean/vmean
      end if                                !  diagno = .true.
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> Apr. 2000 - added diagnostics of t/s trends implied by srf.flux residuals
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> Nov. 2000 - added global density diagnostis to global t/s diagnostics
c> Apr. 2001 - eliminated stmt_funcs.h
c> Apr. 2001 - added diapycnal flux forcing at lateral boundaries
c> June 2001 - eliminated -buoyfl- calculation (now done in mxlayr)
c> Aug. 2001 - improved handling of flux calculations over water and ice
c> Dec. 2001 - fixed bug in -evapw- calculation (replaced max by min)
