#include "rundeck_opts.h"

      subroutine obio_kpar(kmax,vrbos,i,j,im,jm,kdm,nstep,dtsrc,dxypo,
     &                     ogrid,mo,g0m,s0m,grav)
 
      USE obio_dim
      USE obio_com,   only : npst,npnd,p1d,Kd,Kpar
     .                      ,delta_temp1d,temp1d

      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      USE SW2OCEAN, only : lsrpd,fsr
      use ocalbedo_mod, only: nlt


      implicit none

      type(DIST_GRID),intent(in) :: ogrid
      integer, intent(in) :: im,jm,kdm,nstep
      real, intent(in) :: dtsrc,dxypo(jm),grav,     
     &                mo(im,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm),
     &                g0m(im,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm),
     &                s0m(im,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) 
      !real, intent(inout) :: fsr(lsrpd)


      integer i,j,k
      integer nl,ih,icd,ich,ntr,kmax

      real Ebotq,actot,bctot,bbctot,a,bt,bb
      real acdom450,bbc(10),Etopq,zd,zirrq,chl,chlm,fac
      real*8 temgsp


      real Edz(nlt,kdm),Esz(nlt,kdm)
      real Euz(nlt,kdm)
      real Edtop(nlt),Estop(nlt)
      real fchl(nchl)
      real g,s,pres,delta_g

      logical vrbos

      data bbc / 0.002, 0.00071, 0.0032, 0.00071, 0.0029,
     .           0.0,   0.0,     0.0,    0.0,     0.0/


          !integrate kd to get kpar
          Kpar(k) = 0.0
          delta_temp1d(k) = 0.0
          do nl = npst,npnd
             Kpar(k) = Kpar(k) + Kd(nl,k)   !in W/m2
          enddo !nl
          pres=pres+MO(I,J,k)*GRAV*.5
          g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
          s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
          !temperature change due to Kpar
          delta_g =      Kpar(k)              ! W/m2
     .                  * 1.                  !  -> Joules/s/m2
     .                  / mo(i,j,k)           !  -> Joules/kg/s
     .                  * dtsrc               !  -> Joules/kg
          delta_temp1d(k) = temp1d(k) - TEMGSP(g+delta_g,s,pres)
          !add missing pressure to get to the bottom of layer k
          pres=pres+MO(I,J,k)*GRAV*.5
!!!!!!!!!!!! need to write hycom implementation

          !compute fractions
          if (i.eq.40.and.j.eq.40) then
          write(*,'(a,i9,3i5,3e12.4,i5,2e12.4)')'kpar=',
     .        nstep,i,j,k,p1d(k),kpar(k),delta_temp1d(k),
     .        lsrpd,fsr(k),kpar(k)/kpar(1)
          endif
          if (vrbos) then
          write(*,'(a,i9,3i5,3e12.4,i5,2e12.4)')'kpar=',
     .        nstep,i,j,k,p1d(k),kpar(k),delta_temp1d(k),
     .        lsrpd,fsr(k),kpar(k)/kpar(1)
          endif


      !compute par ratios
      do k = 1,kmax
          fsr(k) = kpar(k) / kpar(1)
      enddo

      return
      end
