      subroutine argos_em()
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
c
      integer i_pcgw,l_pcgw,i_pcgs,l_pcgs
c
      llong=ltwin
c
c     center of mass coordinates
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     periodic boundary conditions
c
      call argos_fold(int_mb(i_iw),int_mb(i_is),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     atom redistribution
c
      call argos_space_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),
     + dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
c     center of mass coordinates
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
      call argos_cafe_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call argos_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
      if(msdit.gt.0) call argos_stdesc(int_mb(i_iw+(lwdyn-1)*mwm),
     + dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs),dbl_mb(i_fs),
     + dbl_mb(i_wws),dbl_mb(i_wws+mwa),int_mb(i_mm),dbl_mb(i_fm),
     + dbl_mb(i_xsm))
c
      if(mcgit.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'pcgw',l_pcgw,i_pcgw))
     + call md_abort('Failed to allocate memory for pcgw',me)
      if(.not.ma_push_get(mt_dbl,3*msa,'pcgs',l_pcgs,i_pcgs))
     + call md_abort('Failed to allocate memory for pcgs',me)
      call argos_congra(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),dbl_mb(i_pcgw),
     + int_mb(i_is+(lsdyn-1)*msa),dbl_mb(i_xs),dbl_mb(i_ys),
     + dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_pcgs),dbl_mb(i_wws),
     + dbl_mb(i_wws+mwa))
      if(.not.ma_pop_stack(l_pcgs))
     + call md_abort('Failed to deallocate memory for pcgs',me)
      if(.not.ma_pop_stack(l_pcgw))
     + call md_abort('Failed to deallocate memory for pcgw',me)
      endif
c
      call argos_cafe_print_energy(lfnout)
c
      call argos_space_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
c
      return
      end
c $Id$
