      logical function analyz(irtdb0)
c
c $Id$
c
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
c
      integer irtdb0
c
c     NWChem MD analysis module
c
c     Author: T.P.Straatsma
c       High Performance Computational Chemistry
c       Environmental Molecular Sciences Laboratory
c       Pacific Northwest National Laboratory
c     Copyright 1998 Pacific Northwest National Laboratory
c
c     me : number current node
c     np : number of nodes
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      if(me.eq.0) then
      write(6,*)
      write(6,*)
      call util_print_centered(6,'Analysis Module',40,.true.)
      write(6,*)
      write(6,*)
      call ecce_print_module_entry('analyze')
      endif
c
      irtdb=irtdb0
c
      call ana_task()
c
      analyz=.true.
c
      return
      end
