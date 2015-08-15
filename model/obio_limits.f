#include "rundeck_opts.h"

      subroutine obio_limits(name)

! compute min max for each tracer and check if in bounds

      USE obio_dim


#ifdef OBIO_ON_GARYocean
      USE OCEANRES,       only : idm=>imo,jdm=>jmo,kdm=>lmo
      USE OCN_TRACER_COM, only : ntrcr=>ntm
      USE OCEAN,          only : ip=>focean
      USE obio_com,       only : tracer
#else
      USE hycom_dim_glob, only : kdm,jj,isp,ifp,ilp,ip,ntrcr,idm,jdm
      USE hycom_arrays_glob, only : tracer
#endif

      implicit none
 
      real trminima(11)
      real trmaxima(11)
      data trminima/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.2,1913./
      data trmaxima/40.,0.6,108.,3.5,3.,0.6,0.4,0.3,0.8,16.3,2400./

      integer i,j,k,l

      integer nt,ineg1,jneg1,ipos1,jpos1
      integer kneg1,kpos1
      real tracer_min,tracer_max
      character name*(*)

      do nt=1,ntrcr
       tracer_min= 1.e30
       tracer_max=-1.e30
        ineg1=-1
        jneg1=-1
        kneg1=-1
        ipos1=-1
        ipos1=-1
        kpos1=-1
         do k=1,kdm
         do 1000 j=1,jdm
         do 1000 i=1,idm
         if (ip(i,j).gt.0 .and. tracer(i,j,k,nt).le.tracer_min) then
             tracer_min = tracer(i,j,k,nt)
             ineg1=i
             jneg1=j
             kneg1=k
         endif
         if (ip(i,j).gt.0 .and. tracer(i,j,k,nt).ge.tracer_max) then
             tracer_max = tracer(i,j,k,nt)
             ipos1=i
             jpos1=j
             kpos1=k
         endif
 1000    continue
         enddo
       write(*,'(a,a,i3,a,2(es9.2,1x,3i4))')
     . name,', tracer=',nt,', min,max=',
     . tracer_min,ineg1,jneg1,kneg1,
     . tracer_max,ipos1,jpos1,kpos1
      enddo

      return
      end
