      subroutine argos_cafe_fpme(fw,nwmloc,fs,isq,isgm,nsaloc,grid,mgz,
     + knode,kfrom,lsldo,node)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc,node,mgz
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
      integer knode(*),kfrom(*)
      logical lsldo(*)
      real*8 grid(ngrx,ngry,mgz)
c
      call argos_cafe_forpme(dbl_mb(i_pmew),fw,int_mb(i_iwq),nwmloc,
     + dbl_mb(i_pmes),fs,isq,isgm,nsaloc,
     + dbl_mb(i_chg),dbl_mb(i_theta),dbl_mb(i_dtheta),
     + grid,mgz,knode,kfrom,lsldo,node,int_mb(i_lseq))
c
      return
      end
c $Id$
