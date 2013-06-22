        subroutine tce_mrcc_macpy(k_a,k_b,size)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"

      double precision k_a(*),k_b(*)
      integer size
      integer i,j

      do i=1,size 
        k_b(i) = k_a(i)
      enddo

      return
      end

c $Id$
