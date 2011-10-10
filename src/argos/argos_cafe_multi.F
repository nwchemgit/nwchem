      subroutine argos_cafe_multi(nums,xs,fs,isgan,isatt,
     + isfrc,isdyn,isct1,ess,
     + fss,lfnpmf,lrec,iprpmf,lpmfci)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums,lfnpmf,iprpmf
      real*8 xs(msa,3),fs(msa,3),ess(msf,msf,mpe,2)
      real*8 fss(msf,msf,3,2)
      integer isgan(msa),isatt(msa),isfrc(msa),isdyn(msa),isct1(msa)
      logical lrec,lpmfci
c
      lpmfc=lpmfci
      edrs=zero
      epmf=zero
      ubias=zero
c
      if(mmult.gt.0) then
      call argos_cafe_mlist(nums,xs,fs,isgan,isfrc,isdyn,isct1,
     + mbt(2),int_mb(i_ibnd(2)),mht(2),int_mb(i_iang(2)),
     + mdt(2),int_mb(i_idih(2)),mit(2),int_mb(i_iimp(2)),
     + int_mb(i_ixmul),int_mb(i_imul),dbl_mb(i_xmul),dbl_mb(i_fmul),ess)
      endif
c
      if(ndrs.gt.0) then
      call argos_cafe_dlist(nums,xs,fs,isgan,int_mb(i_idrs),
     + dbl_mb(i_rdrs),dbl_mb(i_xdrs))
      endif
c
      if(lpmf) then
      call argos_cafe_plist(nums,xs,fs,isgan,isatt,int_mb(i_ipmf),
     + int_mb(i_jpmf),dbl_mb(i_rpmf),dbl_mb(i_xpmf),dbl_mb(i_wpmf),
     + dbl_mb(i_upmf),lfnpmf,lrec,iprpmf)
      endif
c
      return
      end
c $Id$
