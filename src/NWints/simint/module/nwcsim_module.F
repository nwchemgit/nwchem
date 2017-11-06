#ifdef USE_SIMINT
      module nwcsim_shell
      use SimintFortran
#include "nwc_const.fh"
C     shells 
      type(c_simint_shell), target :: smnt_sh(nw_max_shells,3)
c dummy shell for 3c- 
      type(c_simint_shell), target :: zero_sh(nw_max_shells,3)
      integer :: nwcsim_noshell(3), nwcsim_bas(3)
      integer :: nwcsim_nbas
      logical :: nwcsim_initialized = .false.
      end module  nwcsim_shell
#else 
c
c     stubs when simint is not used
c     empty module
c
#endif
      
