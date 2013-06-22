      subroutine tce_mrcc_refs_sym()
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
#include "geom.fh"

      character*4 irrepname
      integer i,j,irrep
      integer iref
      integer k
      logical nodezero

       nodezero = (ga_nodeid().eq.0)

        if(nodezero)write(LuOut,"(/,'Symmetry of references',/)")

       do iref=1,nref
        irrep = 0
        do j = 1,nblcks(1,iref)+nblcks(2,iref)
         do k=1,int_mb(k_rangem(iref)+j-1)
            irrep = ieor(irrep,int_mb(k_symm(iref)+j-1))
         enddo
        enddo
        call sym_irrepname(geom,irrep+1,irrepname)
          if(nodezero) then
            write(LuOut,"('Ref.',I4,' sym:',A4)")iref,irrepname
          endif
        enddo

        return
        end

c $Id$
