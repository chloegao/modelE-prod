#include "rundeck_opts.h"
      SUBROUTINE cheminit
!@sum cheminit initialize model chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls jplrts,phtlst,inphot,wave,reactn

C**** GLOBAL parameters and variables:
      USE FILEMANAGER, only: openunit,closeunit
      USE MODEL_COM, only: Itime, ItimeI
      USE DOMAIN_DECOMP_ATM, only : getDomainBounds, grid
      USE TRACER_COM, only: oh_live,no3_live
      USE TRCHEM_Shindell_COM, only:nc,ny,numfam,JPPJ,nn,ks,nps,nds,
     &    ndnr,kps,kds,kpnr,kdnr,nnr,kss,nr,npnr,nr2,nr3,nmm,nhet,
     &    prnls,prnrts,prnchg,lprn,jprn,iprn,ay,nss,pHOx,pOx,pNOx,
     &    yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yNO3,yRXPAR,yXO2N,acetone,
     &    allowSomeChemReinit
     &    ,pCLOx,pCLx,pOClOx,pBrOx,yCl2,yCl2O2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var iu_data temporary unit number
!@var i,l loop dummy
      integer :: iu_data,i,l,j
      integer :: J_0,J_1,J_0S,J_1S,J_1H,J_0H,I_0,I_1
         
      call getDomainBounds(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C Read chem diagnostics parameters and molecule names
C from MOLEC file:
      call openunit('MOLEC',iu_data,.false.,.true.)
      read(iu_data,100)prnls,prnrts,prnchg,lprn,jprn,iprn
      read(iu_data,110)ay
      call closeunit(iu_data)

C Read JPL chemical reactions/rates from unit JPLRX:
      call jplrts

C Read photolysis parameters and reactions from unit JPLPH:
      call phtlst

c fastj initialization routine:
      call inphot

c Set up wavelengths 200 -730 nm & O2 and O3 cross sections:
      call wave

c Set up arrays of reaction numbers involving each molecule:
      call reactn

C Initialize a few (IM,JM,LM) arrays, first hour only:
      IF(Itime == ItimeI .and. allowSomeChemReinit == 1) THEN
        ! allowSomeChemReinit condition b/c these are in RSF files:
        pHOx(I_0:I_1,J_0:J_1,:)     =1.d0
        pOx(I_0:I_1,J_0:J_1,:)      =1.d0
        pNOx(I_0:I_1,J_0:J_1,:)     =1.d0
        yCH3O2(I_0:I_1,J_0:J_1,:)   =1.d0 ! 1.d5 ??
        yC2O3(I_0:I_1,J_0:J_1,:)    =0.d0
        yROR(I_0:I_1,J_0:J_1,:)     =0.d0
        yXO2(I_0:I_1,J_0:J_1,:)     =0.d0
        yAldehyde(I_0:I_1,J_0:J_1,:)=0.d0
        yNO3(I_0:I_1,J_0:J_1,:)     =0.d0
        yXO2N(I_0:I_1,J_0:J_1,:)    =0.d0
        yRXPAR(I_0:I_1,J_0:J_1,:)   =0.d0
        oh_live(I_0:I_1,J_0:J_1,:)  =0.d0
        no3_live(I_0:I_1,J_0:J_1,:) =0.d0
        acetone(I_0:I_1,J_0:J_1,:)  =0.d0
        pClOx(I_0:I_1,J_0:J_1,:)    =1.d0
        pClx(I_0:I_1,J_0:J_1,:)     =0.d0
        pOClOx(I_0:I_1,J_0:J_1,:)   =0.d0
        pBrOx(I_0:I_1,J_0:J_1,:)    =1.d0
        yCl2(I_0:I_1,J_0:J_1,:)     =0.d0
        yCl2O2(I_0:I_1,J_0:J_1,:)   =0.d0
      END IF

 100  format(/3(50x,l1/),3(50x,i8/))
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_TERP
 110  format(6(///10(a8)),(///2(a8)))
#else
 110  format(6(///10(a8)),(///1(a8)))
#endif
#else
#ifdef TRACERS_TERP
 110  format(5(///10(a8)),(///4(a8)))
#else
 110  format(5(///10(a8)),(///3(a8)))
#endif
#endif  /* TRACERS_AEROSOLS_SOA */


      return
      END SUBROUTINE cheminit



      SUBROUTINE jplrts
!@sum jplrts read/set up chemical reaction rates from JPL
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls lstnumc

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: nr,nr2,nr3,nmm,nhet,pe,ea,nst,ro,
     &                               r1,sn,sb,nn,nnr,nc

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
C
!@var ate temporary reactants names array
!@var i,ii,j dummy loop variable
!@var iu_data temporary unit number
      CHARACTER*8, DIMENSION(4) :: ate
      character(len=300) :: out_line
      INTEGER :: i,ii,j,iu_data

C Read in the number of each type of reaction:
      call openunit('JPLRX',iu_data,.false.,.true.)
      read(iu_data,124)nr,nr2,nr3,nmm,nhet
      write(out_line,*)' '
      call write_parallel(trim(out_line))
      write(out_line,*) 'Chemical reactions used in the model: '
      call write_parallel(trim(out_line))

      do i=1,nr               ! >>> begin loop over total reactions <<<
        if(i <= nr-nhet) then !non-hetero
          if(i <= nr2) then   !mono or bi
            if(i > nr2-nmm) then ! read monomolecular reactions
              if(i == nr2-nmm+1) read(iu_data,22)ate
              read(iu_data,16)ate,pe(i),ea(i),nst(i-nr2+nmm)
              write(out_line,30) i,ate(1),' + ',ate(2),
     &        ' --> ',ate(3),' + ',ate(4)
              call write_parallel(trim(out_line))
            else                  ! read bimolecular reactions
   5          read(iu_data,16)ate,pe(i),ea(i)
              write(out_line,30) i,ate(1),' + ',ate(2),
     &        ' --> ',ate(3),' + ',ate(4)
              call write_parallel(trim(out_line))
            end if
          else                    ! read trimolecular reactions
 20         if(i == nr2+1) read(iu_data,22)ate
            ii=i-nr2
            read(iu_data,21)ate,ro(ii),sn(ii),r1(ii),sb(ii)
            write(out_line,30) i,ate(1),' + ',ate(2),
     *      ' --> ',ate(3),' + ',ate(4)
            call write_parallel(trim(out_line))
          end if
        else                     ! read heterogeneous reactions
          if(i == nr-(nhet-1)) read(iu_data,22)ate
          read(iu_data,31)ate
          write(out_line,30) i,ate(1),' + ',ate(2),
     *    ' --> ',ate(3),' + ',ate(4)
          call write_parallel(trim(out_line))
        end if ! (i <= nr-nhet)
c
        do j=1,2
          call lstnum(ate(j),nn(j,i))
          call lstnum(ate(j+2),nnr(j,i))
        end do
      end do                ! >>> end loop over total reactions <<<

 124  format(///5(/43x,i3)///)
  27  format(/(30x,i2))
  21  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f5.2,e9.2,f4.1)
  22  format(/10x,4a8/)
  25  format(//32x,2f7.1,i6)
  16  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f8.0,i4)
  31  format(4x,a8,1x,a8,3x,a8,1x,a8)
  30  format(1x,i3,2x,a8,a3,a8,a5,a8,a3,a8)
      call closeunit(iu_data)
      return
      end SUBROUTINE jplrts



      SUBROUTINE lstnum(at,ks)
!@sum lstnum find molecule number in param list of molecules
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nc,ay

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var at local copy of species name
!@var ks local variable to be passed back to jplrts nnr or nn array.
!@var j dummy loop variable
      INTEGER                  :: j
      INTEGER,     INTENT(OUT) :: ks
      CHARACTER*8, INTENT(IN)  :: at
      
      j=1
      do while(j <= nc)
        if(at == ay(j))then
          ks = j
          return
        else
          j = j + 1
          cycle
        endif
      enddo 
      ks = nc + 1
      if (at /= 'N2' .and. at /= 'H')
     &  call stop_model('ERROR: Tracer '//trim(at)//
     &    ' does not exist in the MOLEC file',255)

      return
      end SUBROUTINE lstnum



      SUBROUTINE phtlst
!@sum phtlst read Photolysis Reactions and parameters
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls lstnum

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel 
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: nss,ks,kss,nc,JPPJ

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var ate species name
!@var nabs,al,nll,nhu,o2up,o3up currently read from JPLPH, but not used
!@var i,j dummy loop variables
!@var iu_data temporary unit number
      INTEGER                   :: nabs,nll,nhu,i,j,iu_data
      CHARACTER*8, DIMENSION(3) :: ate
      character(len=300)        :: out_line
      REAL*8                    :: al,o2up,o3up

C Read in photolysis parameters:
      call openunit('JPLPH',iu_data,.false.,.true.)
      read(iu_data,121)nss,nabs,al,nll,nhu,o2up,o3up
C Check on the number of photolysis reactions:
      IF(nss /= JPPJ)
     &call stop_model('WARNING: nss /= JPPJ, check # photo rxns',255)

c Assign ks and kss gas numbers of photolysis reactants from list:
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*) 'Photolysis reactions used in the model: '
      call write_parallel(trim(out_line))
      do i=1,JPPJ
        read(iu_data,112)ate
        write(out_line,172) i,ate(1),' + hv   --> ',ate(2),' + ',ate(3)
        call write_parallel(trim(out_line))
        call lstnum(ate(1),ks(i))
        do j=2,3
           call lstnum(ate(j),kss(j-1,i))
        end do
      end do
 121  format(//2(45x,i2/),43x,f4.2/44x,i3/45x,i2/2(40x,e7.1/))
 112  format(4x,a8,3x,a8,1x,a8)
 172  format(1x,i2,2x,a8,a12,a8,a3,a8)
      call closeunit(iu_data)

      return
      end SUBROUTINE phtlst



      SUBROUTINE reactn
!@sum reactn read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls guide,printls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nps,nds,kps,kds,nn,nnr,nr,ks,kss,
     &                               JPPJ,npnr,ndnr,kpnr,kdnr,prnls

      IMPLICIT NONE

c Chemical reaction lists:
      call guide(npnr,ndnr,kpnr,kdnr,nn,nnr,2,nr)
c Photolysis reaction lists:
      call guide( nps, nds, kps, kds,ks,kss,1,JPPJ)
C Print out some diagnostics:
      if(prnls) call printls
      
      return
      end SUBROUTINE reactn



      SUBROUTINE guide(npr,ndr,kpr,kdr,xx,nnn,ns,nre)
!@sum guide read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls calcls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: p_1,p_2,p_3,p_4

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var xx   = either nn  or ks   from reactn sub
!@var nnn  = either nnr or kss  from reactn sub
!@var kpr  = either kps or kpnr from reactn sub
!@var kdr  = either kds or kdnr from reactn sub
!@var npr  = either nps or npnr from reactn sub
!@var ndr  = either nds or ndnr from reactn sub
!@var ns   = either 1   or    2 from reactn sub
!@var nre number of reactions
      INTEGER,  DIMENSION(p_4)     :: kpr, kdr
      INTEGER,  DIMENSION(p_3)     :: npr, ndr
      INTEGER, DIMENSION(p_1,p_2)  :: xx, nnn
      INTEGER                      :: ns, nre

c Chemical and photolytic destruction:
      call calcls(xx,ns,nnn,2,ndr,kdr,nre)
c Chemical and photolytic production:
      call calcls(nnn,2,xx,ns,npr,kpr,nre)
      
      return
      end SUBROUTINE guide



      SUBROUTINE calcls(nn,ns,nnn,nns,ndr,kdr,nre)
!@sum calcls Set up reaction lists for calculated gases (1 to ny)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: ny, numfam, p_2, p_3, p_4, nfam,
     &                               prnls

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var kdr  = either kdr or kpr from guide sub
!@var ndr  = either ndr or npr from guide sub
!@var nre  number of reactions
!@var nns number of partic_ on opposite side of reaction
!@var ns   = either ns or   2 from guide sub
!@var nn   = either xx or nnn from guide sub
!@var nnn  = either xx or nnn from guide sub
!@var ii,k,j,i,ij,i2,newfam,ifam dummy variables
      INTEGER, DIMENSION(p_4)    :: kdr
      INTEGER, DIMENSION(p_3)    :: ndr
      INTEGER :: nre, nns, ns, k, j, i, ij, i2, newfam, ifam, ii
      INTEGER, DIMENSION(ns,p_2) :: nn 
      INTEGER, DIMENSION(nns,p_2):: nnn
      character(len=300) :: out_line

      k=1
      nfam(numfam+1)=ny+1
      do j=1,numfam      !families, list only interfamily reactions
        newfam=0
        kdr(j)=k
        i_loop: do i=1,nre    ! 1 to # chem or phot reactions
          ij_loop: do ij=1,ns !ns # partic (prod & chem dest=2,phot dest=1)
            ! check if molecule # nn() is element of family j:
            if(nn(ij,i) >= nfam(j).and.nn(ij,i) < nfam(j+1))then
              ! check if reaction is intrafamily:
              do i2=1,nns  ! nns # partic on opposite side of reac.
                if(nnn(i2,i) >= nfam(j).and.nnn(i2,i) < nfam(j+1))
     &          cycle i_loop
              enddo
              ! don't write same reaction twice:
              if(k /= 1)then
                if(ndr(k-1) == i.and.newfam /= 0) cycle ij_loop
              endif
              ndr(k)=i
              k=k+1
              newfam=1
            endif
          enddo ij_loop
        enddo i_loop
      enddo

      do j=numfam+1,nfam(1)-1     ! individual non-family molecules
        kdr(j)=k
        do i=1,nre                ! 1 to # chem or phot reactions
          do ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j) cycle ! nn is mol # of participant
            ndr(k)=i
            k=k+1
          enddo
        enddo
      enddo

      do 100 j=nfam(1),ny !indiv family mols.,list only intrafamily
        do ii=1,numfam-1
          if(j < nfam(ii+1))then
            ifam=ii
            goto 110
          endif
        enddo
        ifam=numfam
 110    kdr(j)=k
        do 100 i=1,nre          ! 1 to # chem or phot reactions
          do 100 ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j)goto100       ! nn is mol # of participant
c           check that reaction is intrafamily
            do i2=1,nns  ! nns # participants on opposite side of reac.
              if(nnn(i2,i) >= nfam(ifam).and.nnn(i2,i) < nfam(ifam+1))
     &        then
                 ndr(k)=i
                 k=k+1
                 goto100
              endif
            enddo
 100  continue
      kdr(ny+1)=k

      if(prnls)then
        write(*,*) 'nn array size :',k
        write(out_line,*) 'nn array size :',k
        call write_parallel(trim(out_line))
      endif
      
      return
      end SUBROUTINE calcls



      SUBROUTINE printls
!@sum printls print out some chemistry diagnostics (reaction lists)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: kpnr,npnr,kdnr,ndnr,kps,nps,
     &                         ny,nn,nnr,ks,kss,ay,kds,nds,nc

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var ireac,igas,ichange,ii dummy variables
      INTEGER :: ireac,igas,ichange,ii
      character(len=300) :: out_line

c Print reaction lists:
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kpnr(igas+1)-kpnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            if (nnr(2,npnr(ireac)) > nc) then
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ','X'
              call write_parallel(trim(out_line))
            else
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ',ay(nnr(2,npnr(ireac)))
              call write_parallel(trim(out_line))
            end if
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kdnr(igas+1)-kdnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20)
     &      ' Reaction # ',ndnr(ireac),' destroys ',
     *      ay(nn(1,ndnr(ireac))),' and  ',ay(nn(2,ndnr(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*)
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ PHOTOLYTIC PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kps(igas+1)-kps(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20) ' Reaction # ',nps(ireac),' produces ',
     *      ay(kss(1,nps(ireac))),' and  ', ay(kss(2,nps(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     & '______________ PHOTOLYTIC DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kds(igas+1)-kds(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,30) ' Reaction # ',nds(ireac),' destroys ',
     *      ay(ks(nds(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
  10  format(1x,a8)
  20  format(a12,i4,a10,a8,a6,a8)
  30  format(a12,i4,a10,a8) 
       
      return
      end SUBROUTINE printls



      SUBROUTINE wave
!@sum wave Set up Wavelengths 200-730 nm, and O2 & O3 Cross Sections
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: n_bnd3, wlt, sech, sO3, sO2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var lw dummy variable
      INTEGER  :: lw

      do lw=1,n_bnd3
        wlt(lw)=195.d0+5.d0*real(lw)
        ! sech(2,x)=O3 cross sections, sech(1,x)=O2 cross sections
        sech(2,lw)=sO3(lw)
        if(lw <= 9) sech(1,lw)=sO2(lw)
      end do
      
      return
      end SUBROUTINE wave



      SUBROUTINE inphot
!@sum inphot initialise photolysis rate data, called directly from the
!@+   cinit routine in ASAD. Currently use to read the JPL spectral data
!@+   and standard O3 and T profiles and to set the appropriate reaction
!@+   index.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls RD_TJPL,RD_PROF

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: jfacta, jlabel,jppj
     &                  ,MXFASTJ,MIEDX2,title_aer_pf,NAA

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var ipr Photolysis reaction counter
!@var cline dummmy text
!@var i dummy loop variable
!@var iu_data temporary unit number
!@var temp1 temp variable to read in jfacta
!@var temp2 temp variable to read in jlabel
      integer            :: iu_data, ipr, i, L
      character*120      :: cline
      character(len=300) :: out_line
      character*7        :: temp2
      real*8             :: temp1

c Reread the ratj_GISS.d file to map photolysis rate to reaction
c Read in quantum yield jfacta and fastj label jlabel
      ipr=0
      call openunit('RATJ',iu_data,.false.,.true.)
 10   read(iu_data,'(a)',err=20) cline
      if(cline(2:5) == '9999') then
        go to 20
      elseif(cline(1:1) == '#' .or. cline(5:5) == '$') then
        go to 10
      else
        ipr=ipr+1
        backspace iu_data
        read(iu_data,'(78x,f5.1,2x,a7)',err=20) temp1,temp2
        jfacta(ipr) = temp1
        jlabel(ipr) = temp2
        jfacta(ipr)=jfacta(ipr)*1.d-2 
        go to 10
      endif
 20   call closeunit(iu_data)
      if(ipr /= jppj) then
        write(out_line,1000) ipr,jppj
        call write_parallel(trim(out_line),crit=.true.)
        call stop_model('problem with # photolysis reactions',255)
      endif

c Print details:
      write(out_line,1100) ipr
      call write_parallel(trim(out_line))
      do i=1,ipr
        write(out_line,1200) i, jlabel(i), jfacta(i)
        call write_parallel(trim(out_line))
      enddo

c Read in JPL spectral data set:
      call openunit('SPECFJ',iu_data,.false.,.true.)
      call RD_TJPL(iu_data)
      call closeunit(iu_data)

c Read in T & O3 climatology:
      call openunit('ATMFJ',iu_data,.false.,.true.)
      call RD_PROF(iu_data)
      call closeunit(iu_data)

 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,i2,': ',a7,' (Q.Y. ',f6.3,') ')
      return
      end SUBROUTINE inphot



      SUBROUTINE RD_TJPL(NJ1)
!@sum RD_TJPL Read wavelength bins, solar fluxes, Rayleigh parameters,
!@+   T-dependent cross sections and Rayleigh/aerosol scattering phase
!@+   functions with temperature dependences. Current data originates
!@+   from JPL'97.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: NJVAL,WBIN,WL,NWWW,FL,QRAYL,QBC,
     &    Q1D,TQQ,QQQ,NAA,QAAFASTJ,NK,WAAFASTJ,PAA,zpdep,npdep,jpdep,
     &    lpdep,NS,TITLE0,NW1,NW2,TITLEJ,JPPJ,jind,jlabel,jfacta,
     &    title_aer_pf,QO2,QO3,DUMMY,rad_FL
     &    ,SSA,RAA,NP,SF3_fact,SF2_fact,bin4_1991,bin4_1988,bin5_1988

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var NJ1 local copy of unit number to read
!@var i,j,k,iw dummy loop variables
!@var jj dummy variable
!@var nQQQ minus no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
!@var NJVAL2 temporary test for NJVAL= its constant value...
C     bin4_1988 fastj2 bin#4 photon flux for year 1988
C     bin4_1991 fastj2 bin#4 photon flux for year 1991
C     bin5_1988 fastj2 bin#5 photon flux for year 1988
      INTEGER, INTENT(IN) :: NJ1
      INTEGER             :: i,j,k,iw,jj,nqqq,NJVAL2
      character(len=300)  :: out_line

      TQQ = 0.d0

      if(rad_FL == 0)then
        bin4_1991 = 9.431d+11
        bin4_1988 = 9.115E+11
        bin5_1988 = 5.305E+12
      endif

C Read in spectral data:
      READ(NJ1,'(A)') TITLE0
      WRITE(out_line,'(1X,A)') TITLE0
      call write_parallel(trim(out_line))
      READ(NJ1,'(10X,4I5)') NJVAL2,NWWW,NW1,NW2
      IF(NJVAL /= NJVAL2) THEN
        WRITE(out_line,*)'NJVAL (constant)= ',NJVAL,' but it is ',
     &  NJVAL2,'when read in from SPECFJ file.  Please reconcile.'
        call write_parallel(trim(out_line),crit=.true.)
        call stop_model('NJVAL problem in RD_TJPL',255)
      END IF
      NQQQ = NJVAL-3
      READ(NJ1,102) (WBIN(IW),IW=1,NWWW)
      READ(NJ1,102) (WBIN(IW+1),IW=1,NWWW)
      READ(NJ1,102) (WL(IW),IW=1,NWWW)
      if(rad_FL == 0)then
        READ(NJ1,102) (FL(IW),IW=1,NWWW)
      else
        READ(NJ1,102) (DUMMY(IW),IW=1,NWWW)
      endif
      READ(NJ1,102) (QRAYL(IW),IW=1,NWWW)
      READ(NJ1,102) (QBC(IW),IW=1,NWWW)   !From Loiusse et al[JGR,96]

C Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields(each at 3 temps):
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,1),TQQ(K,1), (QO2(IW,K),IW=1,NWWW)
      ENDDO
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,2),TQQ(K,2), (QO3(IW,K),IW=1,NWWW)
      ENDDO
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,3),TQQ(K,3), (Q1D(IW,K),IW=1,NWWW)
      ENDDO
      do k=1,3
        write(out_line,200) titlej(1,k),(tqq(i,k),i=1,3)
        call write_parallel(trim(out_line))
      enddo

C Read remaining species:  X-sections at 2 T's :
      DO J=1,NQQQ
        READ(NJ1,103) TITLEJ(1,J+3),TQQ(1,J+3),(QQQ(IW,1,J),IW=1,NWWW)
        READ(NJ1,103) TITLEJ(2,J+3),TQQ(2,J+3),(QQQ(IW,2,J),IW=1,NWWW)
        write(out_line,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
        call write_parallel(trim(out_line))
      ENDDO
      READ(NJ1,'(A)') TITLE0

C (Don't) read pressure dependencies:
      npdep=0

c Zero index arrays:
      jind=0
      jpdep=0

C Set mapping index:
      do j=1,NJVAL
        do k=1,jppj
          if(jlabel(k) == titlej(1,j)) jind(k)=j
        enddo
        do k=1,npdep
          if(lpdep(k) == titlej(1,j)) jpdep(j)=k
        enddo
      enddo
      do k=1,jppj
        if(jfacta(k) == 0.d0) then
          write(out_line,*) 'Not using photolysis reaction ',k
          call write_parallel(trim(out_line))
        endif
        if(jind(k) == 0) then
          if(jfacta(k) == 0.d0) then
            jind(k)=1
          else
            write(out_line,*)
     &      'Which J-rate for photolysis reaction ',k,' ?'
            call write_parallel(trim(out_line),crit=.true.)
            call stop_model('J-rate problem in RD_TJPL',255)
          endif
        endif
      enddo

C Read aerosol phase functions:
      read(NJ1,'(A10,I5,/)') TITLE0,NAA
      if(NAA > NP)then 
        write(out_line,350) NAA
        call write_parallel(trim(out_line),crit=.true.)
        call stop_model('NAA too large in RD_TJPL',255)
      endif
      NK=4        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) title_aer_pf(j)
        do k=1,NK
          read(NJ1,'(A5,F8.4,F7.3,F8.4,1x,8F6.3)') WAAFASTJ(k,j),
     &    QAAFASTJ(k,j),RAA(k,j),SSA(k,j),(PAA(i,k,j),i=1,8)
        enddo
      enddo

      write(out_line,*) 'Aerosol phase functions & wavelengths'
      call write_parallel(trim(out_line))
      DO J=1,NAA
        write(out_line,'(1x,A8,I2,A,9F8.1)')
     $  title_aer_pf(J),J,'  wavel=',(WAAFASTJ(K,J),K=1,NK)
        call write_parallel(trim(out_line))
        write(out_line,'(9x,I2,A,9F8.4)') J,'  Qext =',
     &  (QAAFASTJ(K,J),K=1,NK)
        call write_parallel(trim(out_line))
      ENDDO   

      if(rad_FL == 0)then
        SF2_fact=FL(5)/bin5_1988
        SF3_fact=0.1d-6*(FL(4)-bin4_1988)/(bin4_1991-bin4_1988)
      endif

  101 FORMAT(8E10.3)
  102 FORMAT((10X,6E10.3)/(10X,6E10.3)/(10X,6E10.3))
  103 FORMAT(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))
  104 FORMAT(13x,i2)
  105 FORMAT(A7,3x,7E10.3)
  110 format(3x,a5)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(1x,' pr.dep:',a10,7(1pE10.3))
  350 format(' Too many phase functions supplied; increase NP to ',i2)
      RETURN
      END SUBROUTINE RD_TJPL



      SUBROUTINE READ_FL(end_of_day)
!@sum READ_FL Instead of reading the photon fluxes (FL) once from the 
!@+   SPECFJ file, this now varyies year-to-year as read from FLTRAN file.
!@+   Format is like the SPECFJ file, data should be consistent with the   
!@+   RADN9 file.
!@auth Greg Faluvegi (based on RD_TJPL above)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: NWWW,FL,FLX,DUMMY,rad_FL
     & ,SF2_fact,SF3_fact,bin4_1991,bin4_1988,bin5_1988
      use model_com, only: modelEclock
      USE RAD_COM, only: s0_yr
      USE RADPAR, only: icycs0,icycs0f

      IMPLICIT NONE
      
C**** Local parameters and variables and arguments:
C bin4_1988 fastj2 bin#4 photon flux for year 1988
C bin4_1991 fastj2 bin#4 photon flux for year 1991
C bin5_1988 fastj2 bin#5 photon flux for year 1988
      integer :: yearx,iunit,i,iw,wantYear,firstYear,lastYear,icyc
      logical, intent(in) :: end_of_day
      character(len=300) :: out_line
      logical :: found1988, found1991
 
      ! only for start of years and restarts:
      if(.not. end_of_day .or. modelEclock%getDayOfYear() == 1) then

        ! set year we are looking for based on rad code s0_yr:
        if(s0_yr==0)then 
          wantYear=modelEclock%getYear()
        else
          wantYear=s0_yr
        end if

        if(wantYear > 2000)then
          icyc=icycs0f
        else
          icyc=icycs0
        end if

        ! scan the file to make sure needed years are there
        ! and to see whether we need to cycle based on the 
        ! initial few or last few years:

        found1988=.false. ; found1991=.false.
        CALL openunit('FLTRAN',iunit,.false.,.true.)
        READ(iunit,*) ! 1 line of comments
        i=0
        scanLoop: do
          i=i+1
          READ(iunit,102,end=100) yearx,(FLX(IW),IW=1,NWWW)
          if(i==1)firstYear=yearx
          if(yearx==1988)found1988=.true.
          if(yearx==1991)found1991=.true.
        end do scanLoop
 100    lastYear=yearx    
        rewind(iunit)
        if(.not.found1988)call stop_model('1988 problem READ_FL',13)
        if(.not.found1991)call stop_model('1991 problem READ_FL',13)
       
        if(lastYear-firstYear+1 < icyc)
     &  call stop_model('years in FLTRAN file < icyc',13)
        if(wantYear < firstYear)then
          write(out_line,*)'READ_FL year ',wantYear,' outside of file.'
          call write_parallel(trim(out_line))
          ! next line depends on integer arithmatic:
          wantYear=wantYear+icyc*((firstYear-wantYear+icyc-1)/icyc)
          write(out_line,*)'Using: ',wantYear,' instead.'
          call write_parallel(trim(out_line))
        else if(wantYear > lastYear)then
          write(out_line,*)'READ_FL year ',wantYear,' outside of file.'
          call write_parallel(trim(out_line))
          ! next line depends on integer arithmatic:
          wantYear=wantYear-icyc*((wantYear-lastYear+icyc-1)/icyc)
          write(out_line,*)'Using: ',wantYear,' instead.'
          call write_parallel(trim(out_line))
        end if

        ! now read file with appropriate (safe) target year:
        READ(iunit,*) ! 1 line of comments
        readLoop: do
          READ(iunit,102,end=101) yearx,(FLX(IW),IW=1,NWWW)
          if(yearx == wantYear) then
            FL(1:NWWW)=FLX(1:NWWW)
          else
            DUMMY(1:NWWW)=FLX(1:NWWW)
          endif
          if(yearx == 1988)then
            if(yearx == wantYear)then
              bin4_1988=FL(4); bin5_1988=FL(5)
            else      
              bin4_1988=DUMMY(4); bin5_1988=DUMMY(5)
            endif
          else if(yearx == 1991)then
            if(yearx == wantYear)then
              bin4_1991=FL(4)
            else
              bin4_1991=DUMMY(4)
            endif
          endif
          if(yearx >= wantYear.and.yearx >= 1991) exit readLoop
        end do readLoop

        write(out_line,*)'READ_FL Using year ',wantYear,
     &  ' bin4_now/1988/1991= ',FL(4),bin4_1988,bin4_1991,
     &  ' bin5_now/1988= ',FL(5),bin5_1988
        call write_parallel(trim(out_line))

        call closeunit(iunit)
      endif

      if(rad_FL > 0)then
        SF2_fact=FL(5)/bin5_1988
        SF3_fact=0.1d-6*(FL(4)-bin4_1988)/(bin4_1991-bin4_1988)
      endif
  102 FORMAT((I4,6X,6E10.3)/(10X,6E10.3)/(10X,6E10.3))
      RETURN 

 101  CONTINUE ! This should no longer be reached.        
      call stop_model("READ_FL end of file problem.",13)
      RETURN 
      END SUBROUTINE READ_FL  


      SUBROUTINE rd_prof(nj2)
!@sum rd_prof input T & O3 reference profiles, define Black Carbon prof.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: TITLE0,
     & TREF2, OREF2, BREF2, ZZHT

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nj2 local unit number
!@var ia,i,m,l,lat,mon,ntlats,ntmons,n216 local dummy variables
      INTEGER, INTENT(IN) :: nj2
      integer :: ia, i, m, l, lat, mon, ntlats, ntmons, n216
      character(len=300) :: out_line
      REAL*8 :: ofac, ofak

      READ(NJ2,'(A)') TITLE0
      WRITE(out_line,'(1X,A)') TITLE0
      call write_parallel(trim(out_line))
      READ(NJ2,'(2I5)') NTLATS,NTMONS
      WRITE(out_line,1000) NTLATS,NTMONS
      call write_parallel(trim(out_line))
      N216 = MIN0(216, NTLATS*NTMONS)
      DO IA=1,N216
        READ(NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = MIN(12, MAX(1, MON))
        L = MIN(18, MAX(1, (LAT+95)/10))
        READ(NJ2,'(3X,11F7.1)') (TREF2(I,L,M), I=1,41)
        READ(NJ2,'(3X,11F7.4)') (OREF2(I,L,M), I=1,31)
      ENDDO
  
c Extend climatology to 100 km:
      ofac=exp(-2.d5/ZZHT)
      do i=32,51
        ofak=ofac**(i-31)
        do m=1,ntmons
          do l=1,ntlats
            oref2(i,l,m)=oref2(31,l,m)*ofak
          enddo
        enddo
      enddo
      do l=1,ntlats
        do m=1,ntmons
          do i=42,51
            tref2(i,l,m)=tref2(41,l,m)
          enddo
        enddo
      enddo

c Approximate Black Carbon up to 10 km; surface 200 ng/m3 (Liousse et
c al) Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-11 g/m2/cm as BREF is in
c cm))
      do i=1,6;  BREF2(i) =10.d0*1.0d-11; end do
      do i=7,51; BREF2(i) =0.d0         ; end do

      return
 1000 format(1x,'Data: ',i3,' Lats x ',i2,' Months')

      end SUBROUTINE rd_prof
