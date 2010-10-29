      logical function diana(irtdb0)
c
c $Id$
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
c
      integer irtdb0
c
c     NWChem Data Intensive Analysis module
c
c     Author: T.P.Straatsma
c             Computational Biology and Bioinformatics Group
c             Computational Sciences and Mathematics Division
c             Pacific Northwest National Laboratory
c
c     Copyright 2006 Pacific Northwest National Laboratory
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
      call util_print_centered
     + (6,'Data Intensive Analysis Module',40,.true.)
      write(6,*)
      write(6,*)
      call ecce_print_module_entry('diana')
      endif
c
      irtdb=irtdb0
c
      call dia_task()
c
      diana=.true.
c
      if(me.eq.0) then
      write(6,*)
      write(6,*)
      call util_print_centered
     + (6,'Data Intensive Analysis Module End',40,.true.)
      write(6,*)
      write(6,*)
      endif
c
      return
      end
