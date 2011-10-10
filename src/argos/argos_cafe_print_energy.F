      subroutine argos_cafe_print_energy(lfnout)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      if(numb(2).gt.0) call argos_cafe_pr_bonds(lfnout,numb(2),
     + mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + byte_mb(i_snam))
c
      if(numh(2).gt.0) call argos_cafe_pr_angles(lfnout,numh(2),
     + mht(2),mhp(2),
     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),
     + dbl_mb(i_rub(2)),
     + byte_mb(i_snam))
c
      if(numd(2).gt.0) call argos_cafe_pr_dihes(lfnout,numd(2),
     + mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),
     + byte_mb(i_snam))
c
      if(numi(2).gt.0) call argos_cafe_pr_imprs(lfnout,numi(2),
     + mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),
     + byte_mb(i_snam))
c
      if(numpmf.gt.0) call argos_cafe_pr_pmf(lfnout,int_mb(i_ipmf),
     + dbl_mb(i_upmf))
c
      return
      end
c $Id$
