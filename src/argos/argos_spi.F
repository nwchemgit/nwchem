      subroutine argos_spi()
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
c
      external argos_space_rdtrj,argos_space_skip,frequency
      logical argos_space_rdtrj,argos_space_skip,frequency
c
      integer import
c
      if(impfr.gt.1) then
      if(.not.argos_space_skip(lfntri,impfr-1)) return
      endif
c
      do 1 import=impfr,impto
c
      if(.not.frequency(import+1-impfr,nftri)) then
      if(.not.argos_space_skip(lfntri,1)) return
      goto 1
      endif
c
      lpair=.true.
      lload=.true.
      lhop=.false.
      llong=ltwin
c
c     read frame from trajectory file
c
      if(.not.argos_space_rdtrj(lfntri,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs),nwmloc,nsaloc)) return
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
      call argos_cafe_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call argos_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
c     atomic forces and potential energies
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call argos_prop_proper(import,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      call argos_prop_step(import,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme)
c
cx      call argos_prop_print()
c
      call rtdb_put(irtdb,'md:energy',mt_dbl,1,epot)
c
c     print energies
c
      if(me.eq.0.and.frequency(import,npener))
     + call argos_cafe_print_energy(lfnout)
c      if(me.eq.0.and.frequency(import,nfprop)) call argos_prop_record()
c
      if(me.eq.0.and.frequency(import,npforc))
     + call argos_space_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
c
      if(ifidi.ne.0) then
      call argos_fd(int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),dbl_mb(i_fs),
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw))
      endif
c
      if(itest.eq.1) call argos_test()
c
    1 continue
c
      return
      end
c $Id$
