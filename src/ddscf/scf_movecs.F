      subroutine scf_movecs_read
      implicit none
#include "cscf.fh"
#include "mafdecls.fh"
      logical movecs_read
      external movecs_read
c
      if (.not. movecs_read(movecs_in, 1,
     $     dbl_mb(k_occ), dbl_mb(k_eval),
     $     g_movecs))
     $     call errquit('scf_movecs_read failed',0)
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
      if (.not. movecs_write(movecs_out, title, trans,
     $     nbf, 1, nbf,
     $     dbl_mb(k_occ), nbf, dbl_mb(k_eval),
     $     nbf, g_movecs)) call errquit
     $     ('scf_movecs_write: movecs_write failed', 0)
c
      end
