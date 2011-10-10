      subroutine argos_space_print()
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer i_lcnt,l_lcnt
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' DOMAIN DECOMPOSITION',/)
c
      write(lfnout,1001) np,npx,npy,npz
 1001 format(' Processor count ',i5,' =',i5,' x',i5,' x',i5)
      write(lfnout,1002) nbx*nby*nbz,nbx,nby,nbz
 1002 format(' Cell count      ',i5,' =',i5,' x',i5,' x',i5)
c
      if(mod(nbx,npx)+mod(nby,npy)+mod(nbz,npz).ne.0) then
      write(lfnout,1003)
 1003 format(/,' WARNING: Inefficient distribution of cells over ',
     + 'processors')
      endif
c
      write(lfnout,1004) bxmin,bymin,bzmin
 1004 format(/,' Minimum cell size ',f12.6,2(' x',f12.6))
c
      if(nred(1)+nred(2)+nred(3).gt.0) then
      write(lfnout,1005)
 1005 format(/,' Warning: ',/)
      if(nred(1).gt.0) write(lfnout,1006) 'x',nred(1)
      if(nred(2).gt.0) write(lfnout,1006) 'y',nred(2)
      if(nred(3).gt.0) write(lfnout,1006) 'z',nred(3)
 1006 format(' Reduced number of cells in ',a,'-dimension: ',i5)
      endif
c
      if(nable.eq.1) write(lfnout,1007)
 1007 format(/,' Read previous box pair list',/)
c
      if(nable.eq.2) write(lfnout,1008)
 1008 format(/,' Unable to read previous box pair list',/)
c
      endif
c
      if(.not.ma_push_get(mt_int,3*np,'lcnt',l_lcnt,i_lcnt))
     + call md_abort('Failed to allocate memory for lcnt',0)
c
      call argos_space_prtcnt(int_mb(i_lcnt))
c
      if(.not.ma_pop_stack(l_lcnt))
     + call md_abort('Failed to deallocate lcnt',0)
c
      return
      end
c $Id$
