      subroutine argos_cafe_add_esa(esa)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 esa(nsa,2)
c
      if(numb(2).gt.0) call argos_cafe_add_ebonds(numb(2),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),esa)
c
      if(numh(2).gt.0) call argos_cafe_add_eangles(numh(2),mht(2),
     + mhp(2),int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),
     + dbl_mb(i_rub(2)),esa)
c
      if(numd(2).gt.0) call argos_cafe_add_edihes(numd(2),mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),esa)
c
      if(numi(2).gt.0) call argos_cafe_add_eimprs(numi(2),mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),esa)
c
      return
      end
c $Id$
