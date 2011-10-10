      subroutine argos_cafe_print_solute(lfnout,imul)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfnout
      integer imul(mmult,4)
      integer i
c
      write(lfnout,1000)
 1000 format(/,' SOLUTE',//)
c
      write(lfnout,1002) totchg
 1002 format(' Total solute charge is ',f20.10,/)
c
      if(nmult(1)+nmult(2)+nmult(3)+nmult(4).le.0) return
c
      write(lfnout,1001) (imul(i,1),i=1,nmul)
 1001 format(' Atom involved in multi-node bonded interaction',i5)
c
      return
      end
c $Id$
