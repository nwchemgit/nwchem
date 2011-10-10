      subroutine argos_pme_coeff(bmod,barr)
c
      implicit none
c
#include "argos_pme_common.fh"
c
      real*8 bmod(ngmax,3),barr(ngmax)
c
      real*8 w,arr(25),darr(25)
      integer i
c
      w=zero
c
      call argos_pme_splfil(w,arr,darr)
c
      do 1 i=1,ngmax
      barr(i)=zero
    1 continue
c
      do 2 i=2,morder+1
      barr(i)=arr(i-1)
    2 continue
c
      call argos_pme_dftmod(bmod(1,1),barr,ngx)
      call argos_pme_dftmod(bmod(1,2),barr,ngy)
      call argos_pme_dftmod(bmod(1,3),barr,ngz)
c
      return
      end
c $Id$
