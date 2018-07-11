      subroutine obio_rivers

      use MODEL_COM, only: dtsrc
      USE obio_incom, only: estFe
      use obio_com, only: ralkconc_loc
      use obio_com, only:  rpocconc_loc,rnitrconc_loc,rsiliconc_loc
     .                    ,rironconc_loc,rdocconc_loc,rdicconc_loc
!                         ,rnitrfmlo_loc

!       if (oFLOWO(i,j) .gt. 0.)then
!        rnitr_loc = rnitrmflo_loc
!     .    / (oFLOWO(i,j) * dxypo(j)) ! kg/s => kg,N/kg,w/s
!     .    * 1.d3     ! kg,N to g,N
!     .    * (1./14.) ! g,N to mol,N
!     .    * 1.d3     ! mol,N to mmol,N
!     .    * rho_water ! kg,water to m3 water
!         if (i.eq.169.and. j.eq.59) then
!           write(*,'(/,a,2i5,6e12.4)')'i,j,rnitrmflo, rnitr, oFLOWO,
!     .       dxypo, nitr, rho_water:',i,j,rnitrmflo_loc(i,j),
!     .       rnitr_loc(i,j),
!     .       oFLOWO(i,j),dxypo(j),obio_P(1,1),rho_water
!       else
!         rnitr_loc = 0.
!       endif

         term = rnitrconc_loc(i,j)
     .    * river_runoff/dtsrc         ! kg,N/kg,w => kg,N/m2,w/s
     .    / dp1d(1)                   ! kg,N/m2,w/s => kg,N/m3,w/s
     .    * 1.d6/14.                  ! kg,N/m3,w/s => mmol,N/m3,w/s
     .    * 3600.                     ! mmol,N/m3,w/s => mmol,N/m3,w/hr
         rhs(1,1,17) = term
         P_tend(1,1) = P_tend(1,1) + term
!       if (i.eq.169 .and. j.eq.59) then
!         write(*,'(/,a,2i5,5e12.4)')'i,j,rnitrconc,rnitr,oFLOWO,
!     .       dtsrc,dp1d(1):',i,j,rnitrconc_loc(i,j),rhs(1,1,17),
!     .       oFLOWO(i,j),dtsrc,dp1d(1) 
!                endif

        term = rsiliconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,S/kg,w => kg,S/m2,w/s
     .   / dp1d(1)                   ! kg,S/m2,w/s => kg,S/m3,w/s
     .   * 1.d6/28.055               ! kg,S/m3,w/s => mmol,S/m3,w/s
     .   * 3600.                     ! mmol,S/m3,w/s => mmol,S/m3,w/hr
        rhs(1,3,17) = term
        P_tend(1,3) = P_tend(1,3) + term

        term = rironconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,Fe/kg,w => kg,Fe/m2,w/s
     .   / dp1d(1)                   ! kg,Fe/m2,w/s => kg,Fe/m3,w/s
     .   * 1.d9/55.845              ! kg,Fe/m3,w/s => umol,Fe/m3,w/s
     .   * 3600.                     ! umol,Fe/m3,w/s => umol,Fe/m3,w/hr
     .   * estFe                     ! estuarine retention rate
        rhs(1,4,17) = term
        P_tend(1,4) = P_tend(1,4) + term

        term = rpocconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,C/kg,w => kg,C/m2,w/s
     .   / dp1d(1)                   ! kg,C/m2,w/s => kg,C/m3,w/s
     .   * 1.d6                      ! kg,C/m3,w/s => mg,C/m3,w/s
     .   * 3600.                     ! mg,C/m3,w/s => mg,C/m3,w/hr
        rhs(1,10,17) = term
        D_tend(1,1) = D_tend(1,1) + term

#ifdef OBIO_RUNOFF
        term = rdocconc_loc(i,j)
     .    * river_runoff/dtsrc          ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
     .    * 3600.                      ! mmol,C/m3,w/s => mmol,C/m3,w/hr
        rhs(1,13,17) = term
        C_tend(1,1) = C_tend(1,1) + term

        term = rdicconc_loc(i,j)
     .    * river_runoff/dtsrc          ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
     .    *3600.                       ! mmol,C/m3,w/s => mmol,C/m3,w/hr
        rhs(1,14,17) = term
        C_tend(1,2) = C_tend(1,2) + term

#ifdef TRACERS_Alkalinity
      term = ralkconc_loc(i,j)
     . * river_runoff/dtsrc           ! mol/kg,w => mol/m2,w/s
     . / dp1d(1)                     ! mol/m2,w/s => mol/m3,w/s
     . * 1.d3                        ! mol/m3,w/s => mmol/m3,w/s
!    . * 3600.                       ! mmol/m3,w/s => mmol/m3/hr    commented out July 2016
      rhs(1,15,17) = term
      A_tend(1) = A_tend(1) + term
#endif 

      end subroutine obio_rivers
