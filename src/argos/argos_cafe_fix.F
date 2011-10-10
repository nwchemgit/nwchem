      subroutine argos_cafe_fix(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,xwm,xs)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnout,lfn,numfix,numw,nums
      integer iwdt(mwm),iwgmn(mwm),isgan(msa),isatt(msa),isdt(msa)
      integer issgm(msa)
      real*8 xwm(mwm,3),xs(msa,3)
c
      integer i_temp,l_temp
c
      if(.not.ma_push_get(mt_int,nsatot,'itemp',l_temp,i_temp))
     + call md_abort
     + ('Unable to allocate temp array in argos_cafe_fix',0)
c
      call argos_cafe_fix_sw(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,byte_mb(i_snam),xwm,xs,
     + int_mb(i_ibnd(2)),mbt(2),numb(2),int_mb(i_temp))
c
      if(.not.ma_pop_stack(l_temp))
     + call md_abort
     + ('Unable to deallocate temp array in argos_cafe_fix',0)
      return
      end
c $Id$
