#include "../model/include/rundeck_opts.h"
!#include "model/testSuites.inc"
#include "shared/testSuites.inc"
#include "profiler/testSuites.inc"
#ifdef TRACERS_ON
#include "tracers/testSuites.inc"
#endif
#ifdef USE_MPI
#include "MPI_Support/testSuites.inc"
#endif
