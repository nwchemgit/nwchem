      subroutine argos_cafe_topol_init(lfnout)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer numbi,lfnout
c
      numbi=ma_sizeof(mt_int,1,mt_byte)
c
      if(me.eq.0) then
      iffld=1
      if(ffield(1:6).eq.'charmm') iffld=2
      if(ffield(1:6).eq.'gromos') iffld=3
      endif
      call ga_brdcst(mcf_69,iffld,numbi,0)
c
c      if(iffld.eq.1) then
c      if(me.eq.0) write(lfnout,1001)
c 1001 format(' Force field is amber')
c      endif
c      if(iffld.eq.2) then
c      if(me.eq.0) write(lfnout,1002)
c 1002 format(' Force field is charmm')
c      endif
c      if(iffld.eq.3) then
c      if(me.eq.0) write(lfnout,1003)
c 1003 format(' Force field is gromos')
c      endif
c
      return
      end
c $Id$
