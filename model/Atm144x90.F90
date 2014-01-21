#include "rundeck_opts.h"
module HorizontalRes
!@sum Horizontal Resolution file, 2x2.5 Lat-Lon Grid
!@auth NCCS ASTG
  Implicit None
!@var IM,JM = longitudinal and latitudinal number of grid cells
#ifdef SCM
  Integer*4,Parameter :: IM=1,JM=1
#else
  Integer*4,Parameter :: IM=144,JM=90
#endif
end module HorizontalRes 
