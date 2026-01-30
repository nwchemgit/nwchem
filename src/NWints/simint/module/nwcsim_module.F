#ifdef USE_SIMINT
      module nwcsim_shell
c      use ISO_FORTRAN_ENV
      use SimintFortran
#include "nwc_const.fh"
#define SIMINT_SCREEN_NONE         0
#define SIMINT_SCREEN_SCHWARZ      1
#define SIMINT_SCREEN_FASTSCHWARZ  2
C     shells 
      integer, parameter :: max_bases=4
      type(c_simint_shell), target :: smnt_sh(nw_max_shells,max_bases)
c dummy shell for 3c- 
      type(c_simint_shell), target :: zero_sh(nw_max_shells,max_bases)
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
      integer :: nwcsim_noshell(max_bases), nwcsim_bas(max_bases)
      integer :: nwcsim_nbas
      logical :: nwcsim_initialized = .false.
      double precision :: smnt_screen_tol 
c      integer (int64) :: smnt_screen_method=SIMINT_SCREEN_FASTSCHWARZ
      integer*8 smnt_screen_method
      integer*4 smnt_deriv
      parameter(smnt_deriv=1)
      end module  nwcsim_shell
#else 
c
c     stubs when simint is not used
c     empty module
c
#endif
      
