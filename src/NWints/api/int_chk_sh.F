      logical function int_chk_sh(basisin,shell)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "geobasmapP.fh"
c::passed
      integer basisin, shell
c::local
      integer basis
c
      int_chk_sh = bas_check_handle(basisin,'int_chk_sh')
      if (.not.int_chk_sh) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
      int_chk_sh = shell.gt.0 .and.
     &             shell.lt.ncont_tot_gb(basis)
      if (.not.int_chk_sh) then
        write(6,*)'invalid contraction/shell index'
        write(6,*)' shell range is: <1:',ncont_tot_gb(basis),'>'
        write(6,*)' input shell is: ',shell
      endif
      end

