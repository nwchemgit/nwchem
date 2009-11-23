      subroutine argos_prop_step(mdstep,stime,eww,esw,ess,fss,esk,
     + epme)
c
c $Id: argos_prop_main.F,v 1.106 2008-07-03 18:05:48 d3j191 Exp $
c
      implicit none
c
#include "argos_prop.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical frequency
      external frequency
c
      integer mdstep
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2),epme(3)
      real*8 fss(msf,msf,3,2)
      real*8 stime,esk(msf)
c
      if(me.eq.0) call argos_prop_stat(mdstep,stime,eww,esw,
     + ess,esk,epme(iset))
c
      if(frequency(mdstep,npener)) then
      call argos_cafe_print_energy(lfnout)
      endif
c
      if(frequency(mdstep,nfprop)) call argos_prop_record()
c
      return
      end
