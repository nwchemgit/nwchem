      subroutine scf_movecs_read
      implicit none
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
      logical movecs_read
      external movecs_read
c
      if (.not. movecs_read(movecs_in, 1,
     $     dbl_mb(k_occ), dbl_mb(k_eval),
     $     g_movecs))
     $     call errquit('scf_movecs_read failed',0)
c
      if (scftype .eq. 'UHF') then
         if (.not. movecs_read(movecs_in, 2,
     $        dbl_mb(k_occ+nbf), dbl_mb(k_eval+nbf),
     $        g_movecs(2))) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) '  Duplicating RHF/ROHF vectors for UHF'
            endif
            call ga_copy(g_movecs(1), g_movecs(2))
            call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
            call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),1)
         endif
      endif
c
      end
      subroutine scf_movecs_write
      implicit none
#include "crohf.fh"
#include "mafdecls.fh"
#include "cscf.fh"
      logical movecs_write
      external movecs_write
c
      integer nset
c
      if (scftype .eq. 'UHF') then
         nset = 2
      else
         nset = 1
      endif
c
      if (.not. movecs_write(movecs_out, title, trans,
     $     nbf, nset, nbf,
     $     dbl_mb(k_occ), nbf, dbl_mb(k_eval),
     $     nbf, g_movecs)) call errquit
     $     ('scf_movecs_write: movecs_write failed', 0)
c
      end
