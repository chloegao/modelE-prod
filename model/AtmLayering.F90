#include "rundeck_opts.h"

module VerticalRes
!@sum Vertical Resolution and Layering file
  use constant, only : grav,mb2kg
  implicit none

!@var PSF global mean surface pressure  (mb)
  real*8, parameter :: PSF = 984d0

  integer, private :: iii ! iterator
  integer, parameter, private :: lmbig=200 ! needed until PDTs are supported

  type layering_t
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of constant-pressure domain
    integer :: lm,ls1
!@var PLbot nominal pressure levels at bottom of layers (mb)
    real*8, dimension(lmbig) :: plbot
  end type layering_t

!
! Define a set of layering options
!

! 12 layers, top at 10 mb
  type(layering_t), parameter, private :: L12 = layering_t( &
    lm = 12, &
    ls1 = 9, &
    plbot = (/  &
      PSF, 934d0, 854d0, 720d0, 550.d0,   &  ! Pbot L=1,5
      390d0, 285d0, 210d0,                &  !      L=...
      150d0,                              &  !      L=LS1
      100d0,  60d0,  30d0, 10d0           &  !      L=..,LM+1
      , (0d0, iii=12+2,lmbig) /) &
      )

! 20 layers, top at .1 mb
  type(layering_t), parameter, private :: L20 = layering_t( &
    lm = 20, &
    ls1 = 11, &
    plbot = (/  &
      PSF, 964d0, 934d0, 884d0, 810d0,   &  ! Pbot L=1,..
      710d0, 550d0, 390d0, 285d0, 210d0, &  !      L=...
      150d0,                             &  !      L=LS1
      110d0,  80d0,  55d0,  35d0,  20d0, &  !      L=...
      10d0,   3d0,   1d0,  .3d0, .1d0    &  !      L=..,LM+1
      , (0d0, iii=20+2,lmbig) /) &
      )

! 40 layers, top at .1 mb
  type(layering_t), parameter, private :: L40 = layering_t( &
    lm = 40, &
    ls1 = 24, &
    plbot = (/  &
      PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0, &  !  L=1,..
      785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0, &  !  L=...
      441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0, &
      197d0, 173d0,                                    &
      150d0,                                           &  !  L=LS1
      128d0, 108d0,  90d0,  73d0,  57d0,  43d0,  31d0, &  !  L=...
      20d0,  10d0,5.62d0,3.16d0,1.78d0,  1.d0,         &
      .562d0,.316d0,.178d0, .1d0                       &  !  L=..,LM+1
      , (0d0, iii=40+2,lmbig) /) &
      )

! Select layering using rundeck CPP symbol ATM_LAYERING
  type(layering_t), parameter, private :: layering = ATM_LAYERING

! Expose the components of the layering selection
  integer, parameter :: &
    lm = layering%lm, &
    ls1 = layering%ls1
  real*8, parameter :: &
    plbot(1:lm+1) = layering%plbot(1:lm+1)

!@var delp nominal pressure thicknesses of layers (mb)
  real*8, dimension(lm), parameter, private :: delp = plbot(1:lm)-plbot(2:lm+1)

!@var PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere

  real*8, parameter :: &
    PTOP = plbot(ls1), &
    PMTOP = plbot(lm+1), &
    PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP

!@var tropomask unity from layers 1-ls1-1, zero above
!@var stratmask zero from layers 1-ls1-1, unity above
  real*8, dimension(lm), parameter, private :: &
    tropomask = (/ (1d0, iii=1,ls1-1), (0d0, iii=ls1,lm) /), &
    stratmask = 1d0-tropomask

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  real*8, parameter :: MDRYA = psf*mb2kg, MTOP = pmtop*mb2kg, MFIXs = pstrat*mb2kg, &
    MFIX(LM)  = delp*stratmask*mb2kg, &
    MFRAC(LM) = delp*tropomask/psfmpt

End Module VerticalRes
