#include "hycom_mpi_hacks.h"
      subroutine advem(iord,fld,u,v,scal,scali,dt,fco,fc)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT,HALO_UPDATE,NORTH,SOUTH,
     &                          haveLatitude,GLOBALSUM
      USE HYCOM_SCALARS, only: itest,jtest

      implicit none
c
c combined monotone scheme, for details see section 3.3 (eqs. 34 to 37)
c in smolarkiewicz and clark, 1986, j.comput.phys.,67,no 2, p. 396-438
c and smolarkiewicz and grabowski, 1989, j.comput.phys.
c  fld    - transported mixing ratio, e.g., salinity or temperature
c  u,v    - mass fluxes satisfying continuity equation
c  scal   - spatial increments (squared)
c  scali  - inverse of scal
c  dt     - temporal increment
c  fco,fc - depth of the layer at previous and new time step
c
      integer i,j,l,n,ia,ib,ja,jb
c
      real fld(idm,J_0H:J_1H),u(idm,J_0H:J_1H),
     .     v(idm,J_0H:J_1H),scal(idm,J_0H:J_1H),
     .     scali(idm,J_0H:J_1H),fco(idm,J_0H:J_1H),
     .     fc(idm,J_0H:J_1H)
      real fmx(idm,J_0H:J_1H),fmn(idm,J_0H:J_1H),
     .     flp(idm,J_0H:J_1H),fln(idm,J_0H:J_1H),
     .     flx(idm,J_0H:J_1H),fly(idm,J_0H:J_1H)
      real u1(idm,J_0H:J_1H),v1(idm,J_0H:J_1H),
     .     flxdiv(idm,J_0H:J_1H),clipj(J_0H:J_1H),
     .     vlumj(J_0H:J_1H)

      real dt,onemu,q,clip,vlume,amount,bfore,after
      integer iord,ip1,im1,jp1,jm1
      logical wrap,recovr
      data recovr/.false./
c
c --- if iord=1, scheme reduces to simple donor cell scheme.
      parameter (onemu=.0098)          !  SI units
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional code for checking conservation properties
ccc      bfore=0.
ccc      do 14 j=1,jj
ccc      do 14 l=1,isp(j)
ccc      do 14 i=ifp(j,l),ilp(j,l)
ccc 14   bfore=bfore+fld(i,j)*fco(i,j)*scal(i,j)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- compute low-order and part of antidiffusive fluxes
c
      fmx=0.d0; fmn=0.d0; fly=0.d0;
      flp=0.d0; fln=0.d0;

      call cpy_p_par(fld)

      CALL HALO_UPDATE(ogrid,fld, FROM=SOUTH+NORTH)
c
      do 11 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      u1(i,j)=.5*abs(u(i,j))*(fld(i,j)-fld(i-1,j))
      if (u(i,j).ge.0.) then
        q=fld(i-1,j)
      else
        q=fld(i  ,j)
      end if
    2 flx(i,j)=u(i,j)*q
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      v1(i,j)=.5*abs(v(i,j))*(fld(i,j)-fld(i,ja ))
      if (v(i,j).ge.0.) then
        q=fld(i,ja )
      else
        q=fld(i,j  )
      end if
    3 fly(i,j)=v(i,j)*q
c
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      ja = PERIODIC_INDEX(j-1, jj)
      if (ip(i,ja).eq.0) ja=j
      jb = PERIODIC_INDEX(j+1, jj)
      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
   11 fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
c
      do 22 j=J_0,J_1
      do 22 l=1,isp(j)
      flx(ifp(j,l)  ,j)=0.
      flx(ilp(j,l)+1,j)=0.
  22  continue
c
      do 33 i=1,ii1
      wrap=jfv(i,1).eq.1	! true if j=1 and j=jj are both water points
      do 33 l=1,jsp(i)
      j=jfp(i,l)
      if (haveLatitude(ogrid, J=j)) then
        if (j.gt.1 .or. .not.wrap) fly(i,j)=0.
      endif
      j=mod(jlp(i,l),jj)+1
      if (haveLatitude(ogrid, J=j)) then
        if (j.gt.1 .or. .not.wrap) fly(i,j)=0.
      endif
   33 continue
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,'(''advem (1)''2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3,
cdiag.1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fld(i-1,j),u(i,j),fld(i,j-1),
cdiag.v(i,j),fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 61 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      vlumj(j)=0.
      clipj(j)=0.
      do 61 l=1,isp(j)
      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      q=fld(i,j)*fco(i,j)-flxdiv(i,j)
      amount=max(fmn(i,j)*fc(i,j),min(q,fmx(i,j)*fc(i,j)))
      if (recovr) then
        vlumj(j)=vlumj(j)+scal(i,j)*fc(i,j)
        clipj(j)=clipj(j)+(q-amount)*scal(i,j)
      end if
   61 fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
      if (iord.le.1) go to 100

      CALL HALO_UPDATE(ogrid,flxdiv, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,fco,    FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,fc,     FROM=SOUTH)
c
c --- finish computation of antidiffusive fluxes
c
      do 8 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 7 l=1,isu(j)
      do 7 i=ifu(j,l),ilu(j,l)
    7 flx(i,j)=u1(i,j)-u(i,j)*(flxdiv(i,j)+flxdiv(i-1,j))
     .   /(fco(i,j)+fco(i-1,j)+fc(i,j)+fc(i-1,j)+onemu)
c
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
    8 fly(i,j)=v1(i,j)-v(i,j)*(flxdiv(i,j)+flxdiv(i,ja ))
     .   /(fco(i,j)+fco(i,ja )+fc(i,j)+fc(i,ja )+onemu)
c
c---- limit antidiffusive fluxes
c

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 16 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 16 l=1,isp(j)
      do 16 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j))*fc(i,j)*scal(i,j)/( (onemu
     .  -min(0.,flx(i+1,j))+max(0.,flx(i,j))
     .  -min(0.,fly(i,jb ))+max(0.,fly(i,j)) )*dt)
      fln(i,j)=(fld(i,j)-fmn(i,j))*fc(i,j)*scal(i,j)/( (onemu
     .  +max(0.,flx(i+1,j))-min(0.,flx(i,j))
     .  +max(0.,fly(i,jb ))-min(0.,fly(i,j)) )*dt)
   16 continue
c
      call cpy_p_par(flp)
      call cpy_p_par(fln)

      CALL HALO_UPDATE(ogrid,fln, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,flp, FROM=SOUTH)
c
      do 18 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 17 l=1,isu(j)
      do 17 i=ifu(j,l),ilu(j,l)
      flx(i,j)=max(0.,flx(i,j))*min(1.,flp(i,j),fln(i-1,j))
     .        +min(0.,flx(i,j))*min(1.,flp(i-1,j),fln(i,j))
   17 continue
c
      do 18 l=1,isv(j)
      do 18 i=ifv(j,l),ilv(j,l)
      fly(i,j)=max(0.,fly(i,j))*min(1.,flp(i,j),fln(i,ja ))
     .        +min(0.,fly(i,j))*min(1.,flp(i,ja ),fln(i,j))
   18 continue

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,'(''advem (2)''2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3,
cdiag.1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fld(i-1,j),u(i,j),fld(i,ja ),
cdiag.v(i,j),fld(i,j),v(i,jb ),fld(i,jb ),u(i+1,j),fld(i+1,j)
c
      do 62 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 62 l=1,isp(j)
      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      q=fld(i,j)*fc(i,j)-flxdiv(i,j)
      amount=max(fmn(i,j)*fc(i,j),min(q,fmx(i,j)*fc(i,j)))
      if (recovr) clipj(j)=clipj(j)+(q-amount)*scal(i,j)
   62 fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
  100 continue
c
c --- revover 'clipped' amount and return to field
c
      if (recovr) then
        vlume=0.
        clip=0.
c
        call GLOBALSUM(ogrid,vlumj,vlume, all=.true.)
        call GLOBALSUM(ogrid,clipj,clip,  all=.true.)

        if (vlume.ne.0.) then
          clip=clip/vlume
cdiag     write (lp,'(a,1pe11.3)') 'tracer drift in advem:',-clip
          do 13 j=J_0,J_1
          do 13 l=1,isp(j)
          do 13 i=ifp(j,l),ilp(j,l)
   13     fld(i,j)=fld(i,j)+clip
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional code for checking conservation properties
ccc      after=0.
ccc      do 15 j=1,jj
ccc      do 15 l=1,isp(j)
ccc      do 15 i=ifp(j,l),ilp(j,l)
ccc 15   after=after+fld(i,j)*fc(i,j)*scal(i,j)
ccc      write (lp,'(a,1p,3e14.6,e11.1)') 'advem conservation:',
ccc     .  bfore,after,after-bfore,(after-bfore)/bfore
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
c
c
c> Revision history:
c>
c> Mar. 2000 - removed 'cushn' and added logic to assure global conservation
c> Apr. 2000 - conversion to SI units
c> Apr. 2000 - changed i/j loop nesting to j/i
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Sep. 2000 - fixed cyclicity problem in loop 33
