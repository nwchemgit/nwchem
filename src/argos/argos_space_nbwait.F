      subroutine argos_space_nbwait(ibbl,lnew,lhandl,lself,lpbcs,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,isl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,lhandl
      integer isl(msa,mis2),iwl(mwm,miw2)
      logical lself,lpbcs,lnew
c
      lpbcs=lpbc9
      call argos_space_nbwait2(ibbl,lnew,int_mb(i_bb),lhandl,lself,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,
     + iwl,int_mb(i_packw),isl,int_mb(i_pack))
c
      return
      end
c $Id$
