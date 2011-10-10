      subroutine argos_cafe_print_top(lfnout,npatom,nptopw,nptops)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnout,npatom,nptopw,nptops
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' TOPOLOGY')
c
      if(npatom.ne.0) then
      call argos_cafe_print_atypes(lfnout,byte_mb(i_nam),int_mb(i_num),
     + dbl_mb(i_mas),dbl_mb(i_vdw))
      call argos_cafe_print_ctypes(lfnout,dbl_mb(i_chg))
      endif
c      call argos_cafe_print_solvent(lfnout)
      if(nptops.ne.0) then
      call argos_cafe_print_solute(lfnout,int_mb(i_imul))
      endif
      endif
c
      return
      end
c $Id$
