c
c $Id: qmmm_cons.F,v 1.2 2005-06-30 00:28:39 marat Exp $
c
c
      subroutine qmmm_cons_init(irtdb)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
c
      if (.not. rtdb_get(irtdb,'qmmm:fixsolvent',mt_log,1,fixw))
     $     fixw=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixsolute',mt_log,1,fixs))
     $     fixs=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqmlink',mt_log,1,freeqmlink))
     $     freeqmlink=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqm',mt_log,1,freeqm))
     $     freeqm=.false.
     
      if(fixw) call mm_fix_solvent()
      if(fixs) call mm_fix_solute(mm_classical)
      if(freeqmlink) call mm_free_solute(mm_quant_link)
      if(freeqm) call mm_free_solute(mm_quant)

      return
      end
c
c



