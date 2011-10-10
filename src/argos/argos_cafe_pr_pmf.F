      subroutine argos_cafe_pr_pmf(lfnout,ipmf,upmf)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfnout,ipmf(numpmf,8)
      real*8 upmf(numpmf)
c
      integer i
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,'Solute potentials of mean force',//,
     + '     i     U',/)
      do 1 i=1,npmfi
      write(lfnout,1001) i,upmf(i)
 1001 format(i7,e12.5)
    1 continue
      endif
      return
      end
c $Id$
