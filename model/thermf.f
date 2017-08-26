      subroutine thermf(m,n,mm,nn,k1m,k1n
     .   ,sss_restore_dt,sss_restore_dtice)
c
      USE SEAICE, only : fsss
c
c --- hycom version 0.9
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : baclin,thref,slfcum,watcum,empcum,nstep
     & ,nstep0,diagno,area,spcifh,avgbot,g,onem,delt1,itest,jtest
     & ,stdsal
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, GLOBALSUM
      use TimeConstants_mod, only: SECONDS_PER_DAY
c
      implicit none
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
      logical vrbos
      real thknss,radfl,radflw,radfli,vpmx,prcp,prcpw,prcpi,
     .     evap,evapw,evapi,exchng,target,old,
     .     rmean,tmean,smean,vmean,boxvol,fxbias,
     &     slfcol(J_0H:J_1H),watcol(J_0H:J_1H),empcol(J_0H:J_1H),
     &     rhocol(J_0H:J_1H),temcol(J_0H:J_1H),salcol(J_0H:J_1H),
     &     sf1col(J_0H:J_1H),sf2col(J_0H:J_1H),clpcol(J_0H:J_1H),
     &     numcol(J_0H:J_1H),fxbiasj(J_0H:J_1H)
      integer iprime,ktop
      real qsatur,totl,eptt,salrlx,sf1cum,sf2cum,bias,numcum
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
ccc        write (*,'(i9,a,f9.5,a,f7.2)') nstep,
ccc     .  '  overall salt gain (psu):',(totl-saltot)*thref/(area*avgbot),
ccc     .  '  => precip bias(%):',100.*pcpcor
ccc      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
c --- for conservation reasons, (E-P)-induced global virtual salt flux must
c --- be proportional to global E-P. this requires a global corrrection.
c --- sf1cum,sf2cum are the uncorrected and corrected global salt fluxes.
c
      do 82 j=J_0,J_1
      sf2col(j)=0.
      fxbiasj(j)=0.
      do 82 l=1,isp(j)
      do 82 i=ifp(j,l),ilp(j,l)
 82   sf2col(j)=sf2col(j)+(-stdsal*oemnp(i,j)/thref+osalt(i,j)*1.e3)
     .  *scp2(i,j)
c
      call GLOBALSUM(ogrid,sf2col,sf2cum, all=.true.)	! desired glob.saltflx
c
      do 85 j=J_0,J_1
c
      watcol(j)=0.
      empcol(j)=0.
      sf1col(j)=0.
      rhocol(j)=0.
      temcol(j)=0.
      salcol(j)=0.
      numcol(j)=0.
c
      do 85 l=1,isp(j)
      do 85 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
c     vrbos=nstep.gt.23000 .and. i.eq.357 .and. j.eq.167
c     vrbos=i.eq.154 .and. j.eq.198
c     vrbos=i.eq.354 .and. j.eq.315
      surflx(i,j)=oflxa2o(i,j)				! heat flux
c
c --- oemnp = evaporation minus precipitation over open water (m/sec)
c --- salflx = salt flux (10^-3 kg/m^2/sec) in +p direction
      salflx(i,j)=-saln(i,j,k1n)*oemnp(i,j)/thref+osalt(i,j)*1.e3
c
c -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
c --- local modifications to salt flux go here:

#ifdef STANDALONE_OCEAN
      if (sss_relax) then
        sssobs => atmocn%sssobs  ! units:  psu/1000
        rsiobs => atmocn%rsiobs  ! sea ice fraction [0-1]
c --- surface salinity restoration
        piston=((1.-rsiobs(i,j))/sss_restore_dt		! open water
     .        +     rsiobs(i,j) /sss_restore_dtice)	! under ice
     .        *12./SECONDS_PER_DAY			! 12m depth scale
        old=salflx(i,j)
        fxbias=(sssobs(i,j)*1000.-saln(i,j,k1n))*piston/thref
        salflx(i,j)=salflx(i,j)+fxbias
        fxbiasj(j)=fxbiasj(j)+fxbias*scp2(i,j)
        if (vrbos)
     .    write (*,'(i9,2i5,a,2f7.2/19x,a,2es10.3,a,f6.2)') nstep,
     .    i,j,'  actual & reference salinity:',saln(i,j,k1n),
     .    sssobs(i,j)*1.e3,'  salflx before & after relax:',old,
     .    salflx(i,j)
      end if
#endif

c --- clip neg.(outgoing) salt flux to prevent S < 0 in top layer
      old=salflx(i,j)
      salflx(i,j)=max(salflx(i,j),
     .  -saln(i,j,k1n)*dp(i,j,k1n)/(g*delt1) * .7)
      if (old.lt.salflx(i,j)) then			! diagnostic use
        numcol(j)=numcol(j)+1.
        print '(2i5,a,f5.1,3x,a,2es14.6)',i,j,' S =',saln(i,j,k1n),
     .    'orig.,modified salflx =',old,salflx(i,j)
      end if
c
c --- clip pos.(incoming) salt flux to prevent S > 40 in top layer
      if (glue(i,j).gt.1.) then
       old=salflx(i,j)
       salflx(i,j)=min(salflx(i,j),
     .   (40.-saln(i,j,k1n))*dp(i,j,k1n)/(g*delt1) * .1)
       if (old.gt.salflx(i,j)) then			! diagnostic use
         numcol(j)=numcol(j)+1.
         print '(2i5,a,f5.1,3x,a,2es14.6)',i,j,' S =',saln(i,j,k1n),
     .    'orig.,modified salflx =',old,salflx(i,j)
       end if
      end if
c
      if (vrbos) write (*,103) nstep,i,j,
     .  'hflx',surflx(i,j),'sflx',salflx(i,j),'sst',temp(i,j,k1n),
     .  'sss',saln(i,j,k1n),'oemp',oemnp(i,j),'oslt',osalt(i,j),
     .  'ice%',oice(i,j)*100.,'loan',loan_ice(i,j)
 103  format (i8,' (thermf) i,j ='2i5/(8(a5,'=',es9.2)))
c -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
c
      sf1col(j)=sf1col(j)+salflx(i,j)*scp2(i,j)
      watcol(j)=watcol(j)+surflx(i,j)*scp2(i,j)
      empcol(j)=empcol(j)+oemnp(i,j)*scp2(i,j)
      rhocol(j)=rhocol(j)+th3d(i,j,k1n)*scp2(i,j)
      temcol(j)=temcol(j)+temp(i,j,k1n)*scp2(i,j)
      salcol(j)=salcol(j)+saln(i,j,k1n)*scp2(i,j)
 85   continue
c
      call GLOBALSUM(ogrid,watcol,watcum, all=.true.)
      call GLOBALSUM(ogrid,empcol,empcum, all=.true.)
      call GLOBALSUM(ogrid,sf1col,sf1cum, all=.true.)
      call GLOBALSUM(ogrid,numcol,numcum, all=.true.)
c
c --- now apply global salt flux constraint to compensate for (a) use of
c --- local salinity, (b) local clipping to prevent S < 0 or S > Smax
      bias=(sf2cum-sf1cum)/area

c --- optional, diagnostic use only:
      if( AM_I_ROOT() ) then
        print '(a,2es15.7,a,es11.3)','actual/desired avg.salt flux:',
     .    sf1cum/area,sf2cum/area,' => bias =',bias
        print '(a,i4,a)','salt flux clipped at',nint(numcum),' points'
      end if
c
      do 83 j=J_0,J_1
      do 83 l=1,isp(j)
      do 83 i=ifp(j,l),ilp(j,l)
 83   salflx(i,j)=salflx(i,j)+bias
c
      if (nstep.eq.nstep0+1 .or. diagno) then

      call GLOBALSUM(ogrid,rhocol,rmean, all=.true.)
      call GLOBALSUM(ogrid,temcol,tmean, all=.true.)
      call GLOBALSUM(ogrid,salcol,smean, all=.true.)
c
      if( AM_I_ROOT() ) then
      write (*,'(i9,''energy residual (w/m^2)'',f9.2)') nstep,
     .  watcum/(area*(nstep-nstep0))
      write (*,'(9x,''resulting temp drift (deg/century):'',f9.3)')
     .  watcum*thref/(spcifh*avgbot*area*(nstep-nstep0)) *
     &  36500.*SECONDS_PER_DAY
      write (*,'(i9,''e - p residual (mm/year)'',f9.2)') nstep,
     .  empcum/(area*(nstep-nstep0))*SECONDS_PER_DAY*365000.
      write (*,'(9x,''saln drift resulting from e-p (psu/century):''
     .                                                     ,f9.3)')
     .  -empcum*35./(avgbot*area*(nstep-nstep0)) *36500.*SECONDS_PER_DAY
css   write (*,'(i9,''salt residual (T/year)'',3f9.2)') nstep,
css  .  sf1cum*365.*SECONDS_PER_DAY*g/onem
      write (*,'(7x,''saln drift resulting from salfl (psu/century):''
     .                                                     ,f9.3)')
     .  sf1cum*36500.*SECONDS_PER_DAY*g/(avgbot*area*onem)
      write (*,'(i9,a,3f9.3)') nstep,' mean surf. sig,temp,saln:',
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
     .    write (*,'(i9,a,3f9.3)') nstep,' mean basin sig,temp,saln:',
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
