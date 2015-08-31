#include "rundeck_opts.h"

!
! Defines a set of layering options and selects one
! via rundeck CPP symbol OCN_LAYERING
!

module olayers

  implicit none
  private

  public :: lmo,lmo_min,dZO

  integer :: iii ! iterator
  integer, parameter :: lmobig=100
  type layering_t
!@param lmo maximum number of ocean layers in a column
!@param lmo_min minimum number of ocean layers in a column
     integer :: lmo,lmo_min
!@param dZO  nominal thicknesses (m) of ocean layers
     real*8, dimension(lmobig) :: dZO
  end type layering_t

  ! L32
  ! L  dZO   ZOE   ZOC      L  dZO   ZOE   ZOC
  ! =  ===   ===   ===      =  ===   ===   ===
  ! 1   12    12     6     17  204  1900  1798
  ! 2   18    30    21     18  206  2100  2003
  ! 3   26    56    43     19  206  2312  2209
  ! 4   36    92    74     20  206  2518  2415
  ! 5   48   140   116     21  206  2724  2621
  ! 6   62   202   171     22  206  2930  2827
  ! 7   78   280   241     23  206  3136  3033
  ! 8   96   376   328     24  206  3342  3239
  ! 9  116   492   434     25  206  3548  3445
  !10  134   626   559     26  206  3754  3651
  !11  150   776   701     27  206  3960  3857
  !12  164   940   858     28  206  4166  4063
  !13  176  1116  1028     29  206  4372  4269
  !14  186  1302  1209     30  206  4578  4475
  !15  194  1496  1399     31  206  4784  4681
  !16  200  1696  1596     32  206  4990  4887
  type(layering_t), parameter, private :: L32 = layering_t( &
       lmo=32, &
       lmo_min=2, &
       dZO = (/ &
         12, 18, 26, 36,  48, 62, 78, 96, &
         116,134,150,164, 176,186,194,200, &
         204,206,206,206, 206,206,206,206, &
         206,206,206,206, 206,206,206,206 &
         , (0, iii=32+1,lmobig) /) &
         )

  ! L13
  ! L  dZO    ZE   ZOC      L  dZO    ZE   ZOC
  ! =  ===    ==   ===      =  ===    ==   ===
  ! 1   12    12     6      7  137   386   318
  ! 2   18    30    21      8  205   591   488
  ! 3   27    57    44      9  308   899   745
  ! 4   41    98    77     10  461  1360  1129
  ! 5   61   158   128     11  692  2052  1706
  ! 6   91   249   204     12 1038  3090  2571
  !                        13 1557  4647  3868
  type(layering_t), parameter, private :: L13 = layering_t( &
       lmo=13, &
       lmo_min=2, &
       dZO = (/ &
         12d0, 18d0, 27d0, 40.5d0, 60.75d0, 91.125d0,  &
         136.6875d0, 205.03125d0, 307.546875d0,  &
         461.3203125d0, 691.98046875d0,  &
         1037.970703125d0, 1556.9560546875d0  &
         , (0d0, iii=13+1,lmobig) /) &
       )

! Select layering using rundeck CPP symbol OCN_LAYERING
  type(layering_t), parameter :: layering = OCN_LAYERING

! Expose the components of the layering selection
  integer, parameter :: &
       lmo     = layering%lmo, &
       lmo_min = layering%lmo_min
  real*8, parameter :: dZO(1:lmo) = layering%dZO(1:lmo)

end module olayers
