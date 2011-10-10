      subroutine argos_cafe_print_ctypes(lfnout,chg)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfnout
      real*8 chg(mqt,mqp,mset)
c
      integer i,j
c
      write(lfnout,1000)
 1000 format(/,' CHARGE TYPES',//,t10,'Charge',/)
c
      do 1 i=1,mqt
      write(lfnout,1001) i,(chg(i,1,j)/qfac,j=1,nparms)
 1001 format(i5,10f12.6)
    1 continue
c
      write(lfnout,1002)
 1002 format(//,t10,'Polarization',/)
c
      do 2 i=1,mqt
      write(lfnout,1003) i,(chg(i,2,j),j=1,nparms)
 1003 format(i5,10f12.6)
    2 continue
 
c
      return
      end
c $Id$
