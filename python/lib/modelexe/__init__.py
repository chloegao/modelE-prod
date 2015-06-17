# Import from the Rundeck-specific wrapped library.
# This gives us access to the ModelE Fortran code.
from f90wrap_modele import *

# Set up fexception-based stop_model(), etc.
modele_python_mod.init_for_python()
