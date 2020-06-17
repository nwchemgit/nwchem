      subroutine sp_setup(node,iwl,xw,xwcr,vw,gw,numwm,
     + isl,xs,xscr,vs,gs,numsa,lpack)
c
c $Id$
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),xscr(msm,3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      logical lpack
c
      integer i_qlst,l_qlst
c
      if(iguide.gt.0) then
      call sp_gagetixvf(node,iwl,int_mb(i_packw),xw,xwcr,vw,gw,numwm,
     + isl,int_mb(i_pack),xs,vs,gs,numsa,int_mb(i_ipl))
      else
      call sp_gagetixv(node,iwl,int_mb(i_packw),xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
      endif
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      call sp_qatoms(isl,numsa)
      if(.not.ma_push_get(mt_int,nsa,'qlst',l_qlst,i_qlst))
     + call md_abort('Failed to allocate memory for qlst',0)
      call sp_latoms(isl,int_mb(i_qlst),numsa)
      if(.not.ma_pop_stack(l_qlst))
     + call md_abort('Failed to deallocate memory for qlst',0)
      call sp_putix(me,iwl,xw,numwm,isl,xs,numsa)
      endif
c
      if(nsm.gt.0) call sp_getxscr(xscr,dbl_mb(i_xscr))
c
      return
      end
      subroutine sp_getxscr(xscr,xscrs)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 xscr(msm,3),xscrs(msm,3)
c
      integer i,j
c
      do 1 j=1,3
      do 2 i=1,nsm
      xscr(i,j)=xscrs(i,j)
    2 continue
    1 continue
c
      return
      end
      subroutine sp_pack_finish()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      if(npackw.gt.0) then
      if(.not.ma_pop_stack(l_packw))
     + call md_abort('Failed to deallocate memory for packw',0)
      endif
      if(npack.gt.0) then
      if(.not.ma_pop_stack(l_pack))
     + call md_abort('Failed to deallocate memory for pack',0)
      endif
c
      return
      end
      subroutine sp_pack_init(isl,numsa,iwl,numwm)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer isl(msa,mis2),numsa,iwl(mwm,miw2),numwm
c
      integer i,j,k,intsiz
c
      intsiz=8*ma_sizeof(mt_int,1,mt_byte)-1
c
      do 1 i=1,mis2
      misl(i)=0
    1 continue
c
      do 2 j=1,mis2
      do 3 i=1,numsa
      misl(j)=max(misl(j),isl(i,j))
    3 continue
    2 continue
      misl(lsdyn)=max(15,misl(lsdyn))
      misl(lsbox)=nbtot
      misl(lsnod)=np
c
      call ga_igop(msp_24,misl,mis2,'max')
c
      do 4 i=1,mis2
      k=misl(i)
      do 5 j=1,intsiz
      k=ishft(k,-1)
      if(k.eq.0) then
      nbits(i)=j+1
      goto 4
      endif
    5 continue
      nbits(i)=intsiz
    4 continue
c
      npack=1
      ipack(1)=1
      k=nbits(1)
      do 6 i=2,mis2
      if(k+nbits(i).lt.intsiz) then
      ipack(npack)=ipack(npack)+1
      k=k+nbits(i)
      else
      npack=npack+1
      ipack(npack)=1
      k=nbits(i)
      endif
    6 continue
c
      do 7 i=1,mis2
      misl(i)=2**nbits(i)-1
    7 continue
c
      do 8 i=1,miw2
      miwl(i)=0
    8 continue
c
      do 9 j=1,miw2
      do 10 i=1,numwm
      miwl(j)=max(miwl(j),iwl(i,j))
   10 continue
    9 continue
      miwl(lwdyn)=15
      miwl(lwbox)=nbtot
      miwl(lwnod)=np
c
      call ga_igop(msp_25,miwl,miw2,'max')
c
      do 11 i=1,miw2
      k=miwl(i)
      do 12 j=1,intsiz
      k=ishft(k,-1)
      if(k.eq.0) then
      nbitw(i)=j+1
      goto 11
      endif
   12 continue
      nbitw(i)=intsiz
   11 continue
c
      npackw=1
      ipackw(1)=1
      k=nbitw(1)
      do 13 i=2,miw2
      if(k+nbitw(i).lt.intsiz) then
      ipackw(npackw)=ipackw(npackw)+1
      k=k+nbitw(i)
      else
      npackw=npackw+1
      ipackw(npackw)=1
      k=nbitw(i)
      endif
   13 continue
c
      do 14 i=1,miw2
      miwl(i)=2**nbitw(i)-1
   14 continue
c
      if(me.eq.0) then
      write(lfnout,1001)
 1001 format(/,' DATA PACKING',/)
      write(lfnout,1002) ' Packing solute in  ',npack,
     + ' integers : ',(ipack(i),i=1,npack)
 1002 format(a,i2,a,15i4)
      write(lfnout,1003) ' significant bits : ',
     + (nbits(i),i=1,mis2)
 1003 format(14x,a,15i4)
      write(lfnout,1002) ' Packing solvent in ',npackw,
     + ' integers : ',(ipackw(i),i=1,npackw)
      write(lfnout,1003) ' significant bits : ',
     + (nbitw(i),i=1,miw2)
      endif
c
      if(npack.eq.mis2) npack=0
      if(npackw.eq.miw2) npackw=0
c      npackw=0
c
      if(npack.gt.0) then
      if(.not.ma_push_get(mt_int,msa*npack,'lpack',l_pack,i_pack))
     + call md_abort('Failed to allocate memory for pack',0)
      endif
      if(npackw.gt.0) then
      if(.not.ma_push_get(mt_int,mwm*npackw,'lpackw',l_packw,i_packw))
     + call md_abort('Failed to allocate memory for packw',0)
      endif
c
      call sp_pack_ga(numsa,isl,int_mb(i_pack),
     + numwm,iwl,int_mb(i_packw))
c
      return
      end
      subroutine sp_pack_ga(numsa,isl,islp,numwm,iwl,iwlp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numsa,isl(msa,mis2),islp(msa,npack)
      integer numwm,iwl(mwm,miw2),iwlp(mwm,npackw)
c
      integer il,ih,jl,jh
c
      if(numsa.gt.0.and.npack.gt.0) then
      call ga_distribution(ga_is,me,il,ih,jl,jh)
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)  
      call sp_pack(numsa,isl,islp)      
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
c
      if(numwm.gt.0.and.npackw.gt.0) then
      call ga_distribution(ga_iw,me,il,ih,jl,jh)
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      call sp_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
c
      call ga_sync()
c
      return
      end
      subroutine sp_pack(numsa,isl,islp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numsa,isl(msa,mis2),islp(msa,npack)
c
      integer i,j,k,ipck
c
      ipck=0
      do 1 j=1,npack
      do 2 i=1,numsa
      islp(i,j)=0
    2 continue
      do 3 k=1,ipack(j)
      ipck=ipck+1
      do 4 i=1,numsa
      islp(i,j)=ishft(islp(i,j),nbits(ipck))+isl(i,ipck)
    4 continue
    3 continue
    1 continue
c
      return
      end
      subroutine sp_unpack(numsa,isl,islp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numsa,isl(msa,mis2),islp(msa,npack)
c
      integer i,j,k,ipck
c
      ipck=mis2+1
      do 1 j=npack,1,-1
      do 2 k=ipack(j),1,-1
      ipck=ipck-1
      do 3 i=1,numsa
      isl(i,ipck)=iand(islp(i,j),misl(ipck))
      islp(i,j)=ishft(islp(i,j)-isl(i,ipck),-nbits(ipck))
    3 continue
    2 continue
    1 continue
c
      return
      end
      subroutine sp_packw(numwm,iwl,iwlp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numwm,iwl(mwm,miw2),iwlp(mwm,npackw)
c
      integer i,j,k,ipck
c
      ipck=0
      do 1 j=1,npackw
      do 2 i=1,numwm
      iwlp(i,j)=0
    2 continue
      do 3 k=1,ipackw(j)
      ipck=ipck+1
      do 4 i=1,numwm
      iwlp(i,j)=ishft(iwlp(i,j),nbitw(ipck))+iwl(i,ipck)
    4 continue
    3 continue
    1 continue
c
      return
      end
      subroutine sp_unpackw(numwm,iwl,iwlp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numwm,iwl(mwm,miw2),iwlp(mwm,npackw)
c
      integer i,j,k,ipck
c
      ipck=miw2+1
      do 1 j=npackw,1,-1
      do 2 k=ipackw(j),1,-1
      ipck=ipck-1
      do 3 i=1,numwm
      iwl(i,ipck)=iand(iwlp(i,j),miwl(ipck))
      iwlp(i,j)=ishft(iwlp(i,j)-iwl(i,ipck),-nbitw(ipck))
    3 continue
    2 continue
    1 continue
c
      return
      end
      subroutine sp_latoms(isl,lst,nums)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isl(msa,mis2),lst(nsa)
c
      integer i
c
      do 1 i=1,nsa
      lst(i)=0
    1 continue
c
      call cf_links(lquant,nums,isl(1,lsgan),isl(1,lsdyn),lst)
c
      call ga_igop(msp_21,lst,nsa,'+')
c
      do 2 i=1,nums
      if(iand(isl(i,lsdyn),lquant).eq.0.and.lst(isl(i,lsgan)).gt.0)
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lqlink)
    2 continue
c
      do 3 i=1,nsa
      lst(i)=0
    3 continue
c
      call cf_links(lqlink,nums,isl(1,lsgan),isl(1,lsdyn),lst)
c
      call ga_igop(msp_26,lst,nsa,'+')
c
      do 4 i=1,nums
      if(iand(isl(i,lsdyn),lqlink).eq.0.and.lst(isl(i,lsgan)).gt.0)
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lqsoft)
    4 continue
c
      return
      end
      subroutine sp_qatoms(isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
c
      integer cf_quantuma
      external cf_quantuma
c
      integer numsa
      integer isl(msa,mis2)
c
      integer i,itype
c
      do 1 i=1,numsa
      itype=cf_quantuma(isl(i,lsatt))
      if(itype.eq.1.or.itype.eq.2) isl(i,lsdyn)=ior(isl(i,lsdyn),lquant)
      if(itype.eq.2) isl(i,lsdyn)=ior(isl(i,lsdyn),lqhigh)
    1 continue
c
      return
      end
      subroutine sp_getixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
      call sp_gagetixv(node,iwl,int_mb(i_packw),xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
      subroutine sp_putp(node,pw,pwp,numwm,ps,psp,numsa,lpp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl+12*mwa+3,jl+18*mwa+2,pw,mwm)
      if(lpp) call ga_put(ga_w,il,ih,jl+18*mwa+3,jl+30*mwa+2,pwp,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl+12,jl+17,ps,msa)
      if(lpp) call ga_put(ga_s,il,ih,jl+18,jl+30,psp,msa)
      endif
c
      return
      end
      subroutine sp_getp(node,pw,pwp,numwm,ps,psp,numsa,lpp,ndx)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa,ndx
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      if(ndx.eq.1) then
      call ga_get(ga_w,il,ih,jl+12*mwa+3,jl+15*mwa+2,pw,mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      elseif(ndx.eq.2) then
      call ga_get(ga_w,il,ih,jl+15*mwa+3,jl+18*mwa+2,pw(1,1,1,2),mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      else
      call ga_get(ga_w,il,ih,jl+12*mwa+3,jl+18*mwa+2,pw,mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      endif
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      if(ndx.eq.1) then
      call ga_get(ga_s,il,ih,jl+12,jl+14,ps,msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      elseif(ndx.eq.2) then
      call ga_get(ga_s,il,ih,jl+15,jl+17,ps(1,1,2),msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      else
      call ga_get(ga_s,il,ih,jl+12,jl+17,ps,msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      endif
      endif
c
      return
      end
      subroutine sp_putix(node,iwl,xw,numwm,isl,xs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      call sp_gaputix(node,iwl,int_mb(i_packw),xw,numwm,
     + isl,int_mb(i_pack),xs,numsa,int_mb(i_ipl))
c
      return
      end
      subroutine sp_update(node,vlatt,
     + iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),vlatt(3,3)
c
      integer i,j
c
      call sp_gaputixv(node,iwl,int_mb(i_packw),xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
c
      do 1 j=1,3
      box(j)=vlat(j,j)
      boxh(j)=half*box(j)
      do 2 i=1,3
      vlat(i,j)=vlatt(i,j)
      vlati(i,j)=vlatt(i,j)
    2 continue
    1 continue
c
      call matinv(vlati,3,3)
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
      subroutine sp_gagetixv(node,iwl,iwlp,xw,xwcr,vw,numwm,isl,islp,
     + xs,vs,numsa,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,ipl,mbox)
      numwm=ipl(1,2)
      numsa=ipl(2,2)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call ga_get(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      call sp_unpackw(numwm,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_get(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_get(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_get(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call ga_get(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      call sp_unpack(numsa,isl,islp)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_get(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_get(ga_s,il,ih,jl+3,jl+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_gagetixvf(node,iwl,iwlp,xw,xwcr,vw,fw,
     + numwm,isl,islp,xs,vs,fs,numsa,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,ipl,mbox)
      numwm=ipl(1,2)
      numsa=ipl(2,2)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call ga_get(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      call sp_unpackw(numwm,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_get(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_get(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_get(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      call ga_get(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call ga_get(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      call sp_unpack(numsa,isl,islp)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_get(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_get(ga_s,il,ih,jl+3,jl+5,vs,msa)
      call ga_get(ga_s,il,ih,jl+9,jl+11,fs,msa)
      endif
c
      return
      end
      subroutine sp_gaputixv(node,iwl,iwlp,xw,xwcr,vw,numwm,isl,islp,
     + xs,vs,numsa,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      ipl(1,2)=numwm
      ipl(2,2)=numsa
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call sp_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_put(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_put(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call sp_pack(numsa,isl,islp)
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_put(ga_s,il,ih,jl+3,jl+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_gaputf(node,fw,numwm,fs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      endif
c
      return
      end
      subroutine sp_gagetf(node,fw,numwm,fs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_get(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_get(ga_s,il,ih,jl+6,jl+8,fs,msa)
      endif
c
      return
      end
      subroutine sp_gaputix(node,iwl,iwlp,xw,numwm,
     + isl,islp,xs,numsa,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      ipl(1,2)=numwm
      ipl(2,2)=numsa
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call sp_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call sp_pack(numsa,isl,islp)
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl,jl+2,xs,msa)
      endif
c
      return
      end
      subroutine sp_owner(xw,iwl,numwm,xs,isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
c     determine solvent ownership
c
      call sp_ownerw(xw,iwl,numwm,dbl_mb(i_boxs),int_mb(i_iown))
c
c     determine solute ownership
c
      call sp_owners(xs,isl,numsa,dbl_mb(i_boxs),int_mb(i_iown))
c
      return
      end
      subroutine sp_ownerw(xw,iwl,numwm,boxsiz,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 xw(mwm,3,mwa)
      real*8 boxsiz(maxbox,3)
      integer iwl(mwm,miw2),ibownr(maxbox,3)
      integer numwm
c
      integer iwm,ibx,iby,ibz,i,k,ipx,ipy,ipz
      real*8 cgx,cgy,cgz,xt(3)
c
c     determine the box and owning node for each solvent molecule
c
      do 1 iwm=1,numwm
      cgx=zero
      cgy=zero
      cgz=zero
      do 2 k=1,nwa
      cgx=cgx+xw(iwm,1,k)
      cgy=cgy+xw(iwm,2,k)
      cgz=cgz+xw(iwm,3,k)
    2 continue
      ibx=0
      iby=0
      ibz=0
      if(nbxtyp.ne.1) then
      xt(1)=cgx
      xt(2)=cgy
      xt(3)=cgz
      else
      xt(1)=box(1)*(vlati(1,1)*cgx+vlati(1,2)*cgy+vlati(1,3)*cgz)
      xt(2)=box(2)*(vlati(2,1)*cgx+vlati(2,2)*cgy+vlati(2,3)*cgz)
      xt(3)=box(3)*(vlati(3,1)*cgx+vlati(3,2)*cgy+vlati(3,3)*cgz)
      endif
      do 3 i=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(i,1)) ibx=i
    3 continue
      do 4 i=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(i,2)) iby=i
    4 continue
      do 5 i=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(i,3)) ibz=i
    5 continue
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      iwl(iwm,lwbox)=(ibz*nby+iby)*nbx+ibx
      iwl(iwm,lwnod)=(ipz*npy+ipy)*npx+ipx
    1 continue
c
      return
      end
      subroutine sp_owners(xs,isl,numsa,boxsiz,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
c
      real*8 xs(msa,3)
      real*8 boxsiz(maxbox,3)
      integer ibownr(maxbox,3)
      integer isl(msa,mis2)
      integer numsa,nfold
c
      integer isa,isaf,isal,jsa,ibox,inod,k
      integer ibx,iby,ibz,i,ipx,ipy,ipz
      real*8 xcgx,xcgy,xcgz,boxi(3),xcg(3)
      real*8 xcgmax,xcgmin,factor
c
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
      nfold=0
c
      isaf=1
      isal=0
      do 1 isa=1,numsa
c
      if(isa.lt.numsa) then
      if(isl(isa+1,lssgm).ne.isl(isaf,lssgm)) isal=isa
      else
      isal=isa
      endif
c
      if(isal.gt.0) then
c
      do 322 k=1,3
      xcgmin=xs(isaf,k)
      xcgmax=xs(isaf,k)
      do 323 jsa=isaf,isal
      xcgmin=min(xcgmin,xs(jsa,k))
      xcgmax=max(xcgmax,xs(jsa,k))
  323 continue
      xcg(k)=0.5d0*(xcgmax+xcgmin)
  322 continue
c
c      xcgx=xcg(1)
c      xcgy=xcg(2)
c      xcgz=xcg(3)
      xcgx=zero
      xcgy=zero
      xcgz=zero
      do 2 jsa=isaf,isal
      xcgx=xcgx+xs(jsa,1)
      xcgy=xcgy+xs(jsa,2)
      xcgz=xcgz+xs(jsa,3)
    2 continue
      factor=one/dble(isal-isaf+1)
      xcgx=factor*xcgx
      xcgy=factor*xcgy
      xcgz=factor*xcgz
cx
      xcgx=xcg(1)
      xcgy=xcg(2)
      xcgz=xcg(3)
c
      if(npbtyp.ne.0) then
      if(abs(xcgx).gt.boxh(1)) then
      xcgx=xcgx-nint(xcgx*boxi(1))*box(1)
      nfold=1
      endif
      if(abs(xcgy).gt.boxh(2)) then
      xcgy=xcgy-nint(xcgy*boxi(2))*box(2)
      nfold=1
      endif
      if(abs(xcgz).gt.boxh(3)) then
      xcgz=xcgz-nint(xcgz*boxi(3))*box(3)
      nfold=1
      endif
      endif
c
      ibx=0
      iby=0
      ibz=0
      do 3 i=1,nbx-1
      if(xcgx+boxh(1).gt.boxsiz(i,1)) ibx=i
    3 continue
      do 4 i=1,nby-1
      if(xcgy+boxh(2).gt.boxsiz(i,2)) iby=i
    4 continue
      do 5 i=1,nbz-1
      if(xcgz+boxh(3).gt.boxsiz(i,3)) ibz=i
    5 continue
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inod=(ipz*npy+ipy)*npx+ipx
      do 6 jsa=isaf,isal
      isl(jsa,lsbox)=ibox
      isl(jsa,lsnod)=inod
    6 continue
c
      isaf=isa+1
      isal=0
c
      endif
c
    1 continue
c
      return
      end
      subroutine sp_travel(bx,xw,vw,xwcr,gw,iwl,numwm,
     + xs,vs,gs,isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),bx(3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
      integer lenx,i_ndx,l_ndx,i_itmp,l_itmp,i_rtmp,l_rtmp
c
      call ga_sync()
c
      box(1)=bx(1)
      box(2)=bx(2)
      box(3)=bx(3)
c
      call sp_bscale(dbl_mb(i_boxs))
c
      call sp_owner(xw,iwl,nwmloc,xs,isl,nsaloc)
c
      lenx=max(nwm,nsa)
      if(.not.ma_push_get(mt_int,lenx,'ndx',l_ndx,i_ndx))
     + call md_abort('Failed to allocate ndx',0)
      if(.not.ma_push_get(mt_int,lenx,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',0)
      if(.not.ma_push_get(mt_dbl,lenx,'rtmp',l_rtmp,i_rtmp))
     + call md_abort('Failed to allocate rtmp',0)
c
      call sp_trvl(xw,vw,xwcr,gw,iwl,int_mb(i_packw),
     + xs,vs,gs,isl,int_mb(i_pack),
     + dbl_mb(i_boxs),int_mb(i_iown),int_mb(i_ipl),
     + int_mb(i_ndx),int_mb(i_itmp),dbl_mb(i_rtmp),lenx)
c
      call ga_sync()
c
      if(.not.ma_pop_stack(l_rtmp))
     + call md_abort('Failed to de-allocate rtmp',0)
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to de-allocate itmp',0)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Failed to de-allocate ndx',0)
c
      call sp_gaputixv(me,iwl,int_mb(i_packw),xw,xwcr,vw,nwmloc,
     + isl,int_mb(i_pack),xs,vs,nsaloc,int_mb(i_ipl))
c
      call sp_lbbl_indices()
c
      numwm=nwmloc
      numsa=nsaloc
c
      return
      end
      subroutine sp_trvl(xw,vw,xwcr,gw,iwl,iwlp,xs,vs,gs,isl,islp,
     + boxsiz,ibownr,ipl,ndx,itmp,rtmp,lenx)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      integer lenx
      real*8 boxsiz(maxbox,3)
      integer ipl(mbox,mip2),ibownr(maxbox,3)
      integer ndx(lenx),itmp(lenx)
      real*8 rtmp(lenx)
      logical lrec(27)
c
      integer i,indexw,indexs,j,k,ibx,iby,ibz,ipx,ipy,ipz
      integer isbox,isnod,nrbox
      integer ilp,ihp,jlp,jhp
      integer il,ih,jl,jh,ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer iliw,ihiw,jliw,jhiw
      integer ilis,ihis,jlis,jhis
      integer iwm,iwstay,jwstay,lwstay,nwgo,nwgosm
      integer nwgtsm
      integer nwnew,nwstay,iwmloc,jwmloc,lwmloc,irw
      integer isa,jsa,isstay,jsstay,icsgm,ifsgm,ilsgm
      integer nsnew,nsstay,lsstay,isaloc,jsaloc,lsaloc,irs
      integer nsgo,jnode,iwfr,iwto,isfr,isto
      real*8 factor,xscx,xscy,xscz,boxi(3)
      integer itemps,nfold
      logical lend
      character*255 string
c
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
      nfold=0
      lpbc9=.false.
c
      nwstay=0
c
c     order the solvent molecules
c
      if(nwmloc.gt.0) then
      do 1 i=1,nwmloc
      ndx(i)=i
    1 continue
      endif
      if(nwmloc.gt.1) then
      lwmloc=nwmloc/2+1
      irw=nwmloc
    2 continue
      if(lwmloc.gt.1) then
      lwmloc=lwmloc-1
      itemps=ndx(lwmloc)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 3
      endif
      endif
      iwmloc=lwmloc
      jwmloc=lwmloc+lwmloc
    4 continue
      if(jwmloc.le.irw) then
      if(jwmloc.lt.irw) then
      if((iwl(ndx(jwmloc),lwnod).eq.iwl(ndx(jwmloc+1),lwnod).and.
     + iwl(ndx(jwmloc),lwbox).le.iwl(ndx(jwmloc+1),lwbox)).or.
     + ((iwl(ndx(jwmloc),lwnod).eq.me.or.
     + (iwl(ndx(jwmloc),lwnod).ne.me.and.
     + iwl(ndx(jwmloc),lwnod).le.iwl(ndx(jwmloc+1),lwnod))).and.
     + iwl(ndx(jwmloc+1),lwnod).ne.me)) jwmloc=jwmloc+1
      endif
      if((iwl(itemps,lwnod).eq.iwl(ndx(jwmloc),lwnod).and.
     + iwl(itemps,lwbox).le.iwl(ndx(jwmloc),lwbox)).or.
     + ((iwl(itemps,lwnod).eq.me.or. (iwl(itemps,lwnod).ne.me.and.
     + iwl(itemps,lwnod).le.iwl(ndx(jwmloc),lwnod))).and.
     + iwl(ndx(jwmloc),lwnod).ne.me)) then
      ndx(iwmloc)=ndx(jwmloc)
      iwmloc=jwmloc
      jwmloc=jwmloc+jwmloc
      else
      jwmloc=irw+1
      endif
      goto 4
      endif
      ndx(iwmloc)=itemps
      goto 2
    3 continue
c
      do 5 k=1,3
      do 8 i=1,nwmloc
      rtmp(i)=xwcr(i,k)
    8 continue
      do 9 i=1,nwmloc
      xwcr(i,k)=rtmp(ndx(i))
    9 continue
      do 10 j=1,nwa
      do 11 i=1,nwmloc
      rtmp(i)=xw(i,k,j)
   11 continue
      do 12 i=1,nwmloc
      xw(i,k,j)=rtmp(ndx(i))
   12 continue
      do 13 i=1,nwmloc
      rtmp(i)=vw(i,k,j)
   13 continue
      do 14 i=1,nwmloc
      vw(i,k,j)=rtmp(ndx(i))
   14 continue
      if(iguide.gt.0) then
      do 113 i=1,nwmloc
      rtmp(i)=gw(i,k,j)
  113 continue
      do 114 i=1,nwmloc
      gw(i,k,j)=rtmp(ndx(i))
  114 continue
      endif
   10 continue
    5 continue
      do 18 k=1,miw2
      do 19 i=1,nwmloc
      itmp(i)=iwl(i,k)
   19 continue
      do 20 i=1,nwmloc
      iwl(i,k)=itmp(ndx(i))
   20 continue
   18 continue
      endif
c
      if(nwmloc.gt.0) then
      do 21 iwm=1,nwmloc
      if(iwl(iwm,lwnod).eq.me) then
      nwstay=iwm
      else
c
c     check if moving atoms go to neighboring processor
c
      do 222 k=1,27
      if(iwl(iwm,lwnod).eq.neighb(k,1)) goto 223
  222 continue
      write(string,'(a,i4,a,i4,a,i4,9f6.2)')
     +  'sp_travel: solvent molecule ',
     + iwl(iwm,lwgmn),' moving to non-neighbor ',iwl(iwm,lwnod),
     + ' from ',me,((xw(iwm,i,j),i=1,3),j=1,3)
      call md_abort(string,me)
  223 continue
c
      endif
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      if(iwl(iwm,lwnod).ne.me) write(lfndbg,'(a,3i5)')
     +  'Travel w fnd ',me,iwl(iwm,lwnod),iwl(iwm,lwgmn)
      endif
c
c     end test code
c
   21 continue
      endif
c
c     order the solute atoms
c
c     isl(isa,lsbox) : box
c     isl(isa,lsnod) : node
c     isl(isa,lssgm) : segment
c
      nsstay=0
      if(nsaloc.gt.0) then
      do 22 i=1,nsaloc
      ndx(i)=i
   22 continue
      endif
c
      if(nsaloc.gt.1) then
      lsaloc=nsaloc/2+1
      irs=nsaloc
   23 continue
      if(lsaloc.gt.1) then
      lsaloc=lsaloc-1
      itemps=ndx(lsaloc)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 24
      endif
      endif
      isaloc=lsaloc
      jsaloc=lsaloc+lsaloc
   25 continue
      if(jsaloc.le.irs) then
      if(jsaloc.lt.irs) then
      if((isl(ndx(jsaloc),lsnod).eq.isl(ndx(jsaloc+1),lsnod).and.
     + (isl(ndx(jsaloc),lsbox).lt.isl(ndx(jsaloc+1),lsbox).or.
     + (isl(ndx(jsaloc),lsbox).eq.isl(ndx(jsaloc+1),lsbox).and.
     + isl(ndx(jsaloc),lssgm).le.isl(ndx(jsaloc+1),lssgm)))).or.
     + ((isl(ndx(jsaloc),lsnod).eq.me.or.
     + (isl(ndx(jsaloc),lsnod).ne.me.and.
     + isl(ndx(jsaloc),lsnod).le.isl(ndx(jsaloc+1),lsnod))).and.
     + isl(ndx(jsaloc+1),lsnod).ne.me)) jsaloc=jsaloc+1
      endif
      if((isl(itemps,lsnod).eq.isl(ndx(jsaloc),lsnod).and.
     + (isl(itemps,lsbox).lt.isl(ndx(jsaloc),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsaloc),lsbox).and.
     + isl(itemps,lssgm).le.isl(ndx(jsaloc),lssgm)))).or.
     + ((isl(itemps,lsnod).eq.me.or. (isl(itemps,lsnod).ne.me.and.
     + isl(itemps,lsnod).le.isl(ndx(jsaloc),lsnod))).and.
     + isl(ndx(jsaloc),lsnod).ne.me)) then
      ndx(isaloc)=ndx(jsaloc)
      isaloc=jsaloc
      jsaloc=jsaloc+jsaloc
      else
      jsaloc=irs+1
      endif
      goto 25
      endif
      ndx(isaloc)=itemps
      goto 23
   24 continue
c
      do 26 k=1,3
      do 27 i=1,nsaloc
      rtmp(i)=xs(i,k)
   27 continue
      do 28 i=1,nsaloc
      xs(i,k)=rtmp(ndx(i))
   28 continue
      do 29 i=1,nsaloc
      rtmp(i)=vs(i,k)
   29 continue
      do 30 i=1,nsaloc
      vs(i,k)=rtmp(ndx(i))
   30 continue
      if(iguide.gt.0) then
      do 2129 i=1,nsaloc
      rtmp(i)=gs(i,k)
 2129 continue
      do 2130 i=1,nsaloc
      gs(i,k)=rtmp(ndx(i))
 2130 continue
      endif
   26 continue
      do 40 k=1,mis2
      do 41 i=1,nsaloc
      itmp(i)=isl(i,k)
   41 continue
      do 42 i=1,nsaloc
      isl(i,k)=itmp(ndx(i))
   42 continue
   40 continue
      endif
c
      if(nsa.gt.0) then
      do 43 isa=1,nsaloc
      if(isl(isa,lsnod).eq.me) then
      nsstay=isa
      else
c
c     check if moving atoms go to neighboring processor
c
      do 444 k=1,27
      if(isl(isa,lsnod).eq.neighb(k,1)) goto 445
  444 continue
      write(string,'(a,i4,a,i4,a,i4,3f6.2)')
     +  'sp_travel: solute segment ',
     + isl(isa,lssgm),' moving to non-neighbor ',isl(isa,lsnod),
     + ' from ',me,(xs(isa,i),i=1,3)
      call md_abort(string,me)
  445 continue
c
      endif
   43 continue
      endif
c
c     make packages ready for shipment
c
c     loop over all neighboring nodes
c
      call ga_distribution(ga_iw,me,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_distribution(ga_is,me,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
c
      indexw=0
      indexs=0
      nwgosm=0
c
      do 70 i=1,27
      jnode=neighb(i,1)
      if(jnode.ge.0.and.jnode.ne.me) then
c
c     for the solvent
c
      iwfr=0
      iwto=0
      do 71 iwm=nwstay+1,nwmloc
      if(iwl(iwm,lwnod).eq.jnode) then
      if(iwfr.eq.0) iwfr=iwm
      iwto=iwm
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      if(iwl(iwm,lwnod).ne.me) write(lfndbg,'(a,3i5)')
     +  'Travel w snd ',me,iwl(iwm,lwnod),iwl(iwm,lwgmn)
      endif
c
c     end test code
c
      endif
   71 continue
c
c     if molecules need to travel copy coordinates etc into global array
c
      nwgo=iwto-iwfr+1
      if(iwfr.eq.0) nwgo=0
      ipl(1,1)=0
      ipl(1,2)=0
c
      if(nwgo.gt.0) then
      nwgosm=nwgosm+nwgo
      il=iliw+indexw
      ih=il+nwgo-1
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,ih,jliw,jhiw,iwl(iwfr,1),mwm)
      else
      call sp_packw(ih-il+1,iwl(iwfr,1),iwlp(iwfr,1))
      call ga_put(ga_iw,il,ih,jliw,jliw+npackw-1,iwlp(iwfr,1),mwm)
      endif
      il=ilw+indexw
      ih=il+nwgo-1
      call ga_put(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(iwfr,1),mwm)
      if(iguide.gt.0) then
      call ga_put(ga_w,il,ih,jlw+6*mwa+3,jlw+9*mwa+2,gw(iwfr,1,1),mwm)
      endif
      ipl(1,1)=indexw+1
      ipl(1,2)=indexw+nwgo
      indexw=indexw+nwgo
      endif
c
c     for the solute
c
      isfr=0
      isto=0
      do 72 isa=nsstay+1,nsaloc
      if(isl(isa,lsnod).eq.jnode) then
      if(isfr.eq.0) isfr=isa
      isto=isa
      endif
   72 continue
      nsgo=isto-isfr+1
      if(isfr.eq.0) nsgo=0
      ipl(1,3)=0
      ipl(1,4)=0
      if(nsgo.gt.0) then
      il=ilis+indexs
      ih=il+nsgo-1
      if(npack.eq.0) then
      call ga_put(ga_is,il,ih,jlis,jhis,isl(isfr,1),msa)
      else
      call sp_pack(ih-il+1,isl(isfr,1),islp(isfr,1))
      call ga_put(ga_is,il,ih,jlis,jlis+npack-1,islp(isfr,1),msa)
      endif
      call ga_put(ga_s,il,ih,jls,jls+2,xs(isfr,1),msa)
      call ga_put(ga_s,il,ih,jls+3,jls+5,vs(isfr,1),msa)
      if(iguide.gt.0) then
      call ga_put(ga_s,il,ih,jls+6,jls+8,gs(isfr,1),msa)
      endif
      ipl(1,3)=indexs+1
      ipl(1,4)=indexs+nsgo
      indexs=indexs+nsgo
      endif
c
c     inform other node of number of molecules to get
c
      if(ipl(1,1).gt.0.or.ipl(1,3).gt.0) then
      call ga_distribution(ga_ip,jnode,ilp,ihp,jlp,jhp)
      ilp=ilp+2+i
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      endif
      endif
      lrec(i)=.false.
   70 continue
c
      call ga_sync()
c
c     receive molecules from other nodes
c
      nwgtsm=0
c
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
      do 74 i=1,27
      jnode=neighb(i,2)
      if(jnode.ge.0.and.jnode.ne.me.and..not.lrec(i)) then
c
      iwfr=ipl(3+i,1)
      iwto=ipl(3+i,2)
      isfr=ipl(3+i,3)
      isto=ipl(3+i,4)
c
      nwnew=iwto-iwfr+1
      nsnew=isto-isfr+1
c
      if(iwfr.eq.0) nwnew=0
      if(isfr.eq.0) nsnew=0
c
      if(nwstay+nwnew.gt.mwm) then
      write(string,'(a,i7,a,i7)')
     + 'Travel: mwm needs increase with ',nwnew,' to ',nwstay+nwnew
      call md_abort(string,me)
      endif
      if(nsstay+nsnew.gt.msa) then
      write(string,'(a,i7,a,i7)')
     + 'Travel: msa needs increase with ',nsnew,' to ',nsstay+nsnew
      call md_abort(string,me)
      endif
c
      lrec(i)=.true.
c
      if(iwfr.gt.0) then
      nwgtsm=nwgtsm+nwnew
      iwto=ipl(3+i,2)
      call ga_distribution(ga_iw,jnode,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,jnode,ilw,ihw,jlw,jhw)
c
c     get data for additional molecules
c
      il=iliw+iwfr-1
      ih=iliw+iwto-1
      if(npackw.eq.0) then
      call ga_get(ga_iw,il,ih,jliw,jhiw,iwl(nwstay+1,1),mwm)
      else
      call ga_get(ga_iw,il,ih,jliw,jliw+npackw-1,iwlp(nwstay+1,1),mwm)
      call sp_unpackw(ih-il+1,iwl(nwstay+1,1),iwlp(nwstay+1,1))
      endif
      call ga_get(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(nwstay+1,1),mwm)
      if(iguide.gt.0) then
      call ga_get(ga_w,il,ih,jlw+6*mwa+3,jlw+9*mwa+2,
     + gw(nwstay+1,1,1),mwm)
      endif
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      write(lfndbg,'(a,3i5)')
     +  ('Travel w rcv ',me,jnode,iwl(nwstay+k,lwgmn),k=1,nwnew)
      endif
c
c     end test code
c
c
c     update number of local solvent molecules
c
      nwstay=nwstay+nwnew
c
      endif
c
c     for the solute
c
      if(isfr.gt.0) then
      call ga_distribution(ga_is,jnode,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,jnode,ils,ihs,jls,jhs)
      il=ilis+isfr-1
      ih=ilis+isto-1
      jl=jlis
      jh=jhis
      if(npack.eq.0) then
      call ga_get(ga_is,il,ih,jlis,jhis,isl(nsstay+1,1),msa)
      else
      call ga_get(ga_is,il,ih,jlis,jlis+npack-1,islp(nsstay+1,1),msa)
      call sp_unpack(ih-il+1,isl(nsstay+1,1),islp(nsstay+1,1))
      endif
      call ga_get(ga_s,il,ih,jls,jls+2,xs(nsstay+1,1),msa)
      call ga_get(ga_s,il,ih,jls+3,jls+5,vs(nsstay+1,1),msa)
      if(iguide.gt.0) then
      call ga_get(ga_s,il,ih,jls+6,jls+8,gs(nsstay+1,1),msa)
      endif
c
      nsstay=nsstay+nsnew
      endif
c
      endif
c
c     reset the pointers to zero
c
      ipl(3+i,1)=0
      ipl(3+i,2)=0
      ipl(3+i,3)=0
      ipl(3+i,4)=0
c
   74 continue
c
c     reset ipl in global array
c
      call ga_put(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
c     order the solvent molecules according to subbox and
c     store indices into ip
c
c     ip(1,1)    : number of boxes on this node
c     ip(1,2)    : number of solvent molecules on this node
c     ip(2,2)    : number of solute atoms on this node
c
c     ip(3+i,1)  : index for solvents to be moved to the i-th neighbor
c
c     ip(30+i,1) : number of i-th box on this node
c     ip(30+i,2) : index to first solvent in i-th box
c     ip(30+i,3) : index to lasst solvent in i-th box
c
      if(nwstay.gt.0.and.(nwgosm.gt.0.or.nwgtsm.gt.0)) then
      do 81 i=1,nwstay
      ndx(i)=i
   81 continue
      if(nwstay.gt.1) then
      lwstay=nwstay/2+1
      irw=nwstay
   82 continue
      if(lwstay.gt.1) then
      lwstay=lwstay-1
      itemps=ndx(lwstay)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 83
      endif
      endif
      iwstay=lwstay
      jwstay=lwstay+lwstay
   84 continue
      if(jwstay.le.irw) then
      if(jwstay.lt.irw) then
      if(iwl(ndx(jwstay),lwbox).le.iwl(ndx(jwstay+1),lwbox))
     + jwstay=jwstay+1
      endif
      if(iwl(itemps,lwbox).le.iwl(ndx(jwstay),lwbox)) then
      ndx(iwstay)=ndx(jwstay)
      iwstay=jwstay
      jwstay=jwstay+jwstay
      else
      jwstay=irw+1
      endif
      goto 84
      endif
      ndx(iwstay)=itemps
      goto 82
   83 continue
c
      do 85 k=1,3
      do 88 i=1,nwstay
      rtmp(i)=xwcr(i,k)
   88 continue
      do 89 i=1,nwstay
      xwcr(i,k)=rtmp(ndx(i))
   89 continue
      do 90 j=1,nwa
      do 91 i=1,nwstay
      rtmp(i)=xw(i,k,j)
   91 continue
      do 92 i=1,nwstay
      xw(i,k,j)=rtmp(ndx(i))
   92 continue
      do 93 i=1,nwstay
      rtmp(i)=vw(i,k,j)
   93 continue
      do 94 i=1,nwstay
      vw(i,k,j)=rtmp(ndx(i))
   94 continue
      if(iguide.gt.0) then
      do 193 i=1,nwstay
      rtmp(i)=gw(i,k,j)
  193 continue
      do 194 i=1,nwstay
      gw(i,k,j)=rtmp(ndx(i))
  194 continue
      endif
   90 continue
   85 continue
      do 98 k=1,miw2
      do 99 i=1,nwstay
      itmp(i)=iwl(i,k)
   99 continue
      do 100 i=1,nwstay
      iwl(i,k)=itmp(ndx(i))
  100 continue
   98 continue
c
      endif
      endif
c
c     order the solute according to segment
c
      if(nsstay.gt.0) then
      do 122 i=1,nsstay
      ndx(i)=i
  122 continue
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  123 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 124
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  125 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lssgm).le.isl(ndx(jsstay+1),lssgm))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lssgm).le.isl(ndx(jsstay),lssgm)) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 125
      endif
      ndx(isstay)=itemps
      goto 123
  124 continue
      endif
c
c     for each segment : 1. determine box number
c                        2. assign box number to each atom
c                        3. when box not owned by node:
c                           a. assign box number
c                           b. assign correct node number
c
      goto 666
      icsgm=isl(ndx(1),lssgm)
      ifsgm=1
      ilsgm=1
      do 126 isa=2,nsstay+1
c
c     if isa is first atom of a new segment or very last atom
c
      
      if(isa.le.nsstay) then
      lend=isl(ndx(isa),lssgm).ne.icsgm
      else
      lend=.true.
      endif
      if(lend) then
      if(isa.gt.nsstay) ilsgm=nsstay
      if(ifsgm.gt.0.and.ilsgm.ge.ifsgm) then
      xscx=zero
      xscy=zero
      xscz=zero
      do 127 jsa=ifsgm,ilsgm
      xscx=xscx+xs(ndx(jsa),1)
      xscy=xscy+xs(ndx(jsa),2)
      xscz=xscz+xs(ndx(jsa),3)
  127 continue
      factor=one/dble(ilsgm-ifsgm+1)
      xscx=factor*xscx
      xscy=factor*xscy
      xscz=factor*xscz
      if(npbtyp.ne.0) then
      if(abs(xscx).gt.boxh(1)) then
      xscx=xscx-nint(xscx*boxi(1))*box(1)
      nfold=1
      endif
      if(abs(xscy).gt.boxh(2)) then
      xscy=xscy-nint(xscy*boxi(2))*box(2)
      nfold=1
      endif
      if(abs(xscz).gt.boxh(3)) then
      xscz=xscz-nint(xscz*boxi(3))*box(3)
      nfold=1
      endif
      endif
c
c     determine the box number
c
      ibx=0
      iby=0
      ibz=0
      do 128 i=1,nbx-1
      if(xscx+boxh(1).gt.boxsiz(i,1)) ibx=i
  128 continue
      do 129 i=1,nby-1
      if(xscy+boxh(2).gt.boxsiz(i,2)) iby=i
  129 continue
      do 1130 i=1,nbz-1
      if(xscz+boxh(3).gt.boxsiz(i,3)) ibz=i
 1130 continue
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      isbox=(ibz*nby+iby)*nbx+ibx
      isnod=(ipz*npy+ipy)*npx+ipx
c
c     assign box and node numbers
c
      do 1131 jsa=ifsgm,ilsgm
      isl(ndx(jsa),lsbox)=isbox
      isl(ndx(jsa),lsnod)=isnod
 1131 continue
c
      endif
      if(isa.le.nsstay) icsgm=isl(ndx(isa),lssgm)
      ifsgm=isa
      else
      ilsgm=isa
      endif
  126 continue
  666 continue
c
c     order solute according to box, segment, charge group, atom number
c
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  132 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 133
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  134 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lsbox).lt.isl(ndx(jsstay+1),lsbox).or.
     + (isl(ndx(jsstay),lsbox).eq.isl(ndx(jsstay+1),lsbox).and.
     + (isl(ndx(jsstay),lssgm).lt.isl(ndx(jsstay+1),lssgm).or.
     + (isl(ndx(jsstay),lssgm).eq.isl(ndx(jsstay+1),lssgm).and.
     + (isl(ndx(jsstay),lsgrp).lt.isl(ndx(jsstay+1),lsgrp).or.
     + (isl(ndx(jsstay),lsgrp).eq.isl(ndx(jsstay+1),lsgrp).and.
     + isl(ndx(jsstay),lsgan).le.isl(ndx(jsstay+1),lsgan)))))))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lsbox).lt.isl(ndx(jsstay),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsstay),lsbox).and.
     + (isl(itemps,lssgm).lt.isl(ndx(jsstay),lssgm).or.
     + (isl(itemps,lssgm).eq.isl(ndx(jsstay),lssgm).and.
     + (isl(itemps,lsgrp).lt.isl(ndx(jsstay),lsgrp).or.
     + (isl(itemps,lsgrp).eq.isl(ndx(jsstay),lsgrp).and.
     + isl(itemps,lsgan).le.isl(ndx(jsstay),lsgan))))))) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 134
      endif
      ndx(isstay)=itemps
      goto 132
  133 continue
      endif
c
      do 135 k=1,3
      do 136 i=1,nsstay
      rtmp(i)=xs(i,k)
  136 continue
      do 137 i=1,nsstay
      xs(i,k)=rtmp(ndx(i))
  137 continue
      do 138 i=1,nsstay
      rtmp(i)=vs(i,k)
  138 continue
      do 139 i=1,nsstay
      vs(i,k)=rtmp(ndx(i))
  139 continue
      if(iguide.gt.0) then
      do 1138 i=1,nsstay
      rtmp(i)=gs(i,k)
 1138 continue
      do 1139 i=1,nsstay
      gs(i,k)=rtmp(ndx(i))
 1139 continue
      endif
  135 continue
      do 149 k=1,mis2
      do 150 i=1,nsstay
      itmp(i)=isl(i,k)
  150 continue
      do 151 i=1,nsstay
      isl(i,k)=itmp(ndx(i))
  151 continue
  149 continue
c
      endif
c
      do 200 i=1,ipl(1,1)
      ipl(30+i,2)=0
      ipl(30+i,3)=0
      ipl(30+i,4)=0
      ipl(30+i,5)=0
  200 continue
c
      do 201 i=1,ipl(1,1)
      nrbox=ipl(30+i,1)
      if(nwstay.gt.0) then
      do 202 iwm=1,nwstay
      if(iwl(iwm,lwbox).eq.nrbox) then
      if(ipl(30+i,2).eq.0) ipl(30+i,2)=iwm
      ipl(30+i,3)=iwm
      endif
  202 continue
      endif
      if(nsstay.gt.0) then
      do 203 isa=1,nsstay
      if(isl(isa,lsbox).eq.nrbox) then
      if(ipl(30+i,4).eq.0) ipl(30+i,4)=isa
      ipl(30+i,5)=isa
      endif
  203 continue
      endif
  201 continue
c
      nwmloc=nwstay
      ipl(1,2)=nwmloc
      nsaloc=nsstay
      ipl(2,2)=nsaloc
c
      call ga_igop(msp_23,nfold,1,'+')
      lpbc9=nfold.gt.0
c
      return
      end
      subroutine sp_nbbl(n)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer n
c
      n=nbbl
      if(lqmd) n=0
c
      return
      end
      subroutine sp_prefetch_all(n,iwl,xw,isl,xs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(2),ps(2)
      real*8 pwp(2),psp(2)
c
      n=nbbl
      if(lqmd) n=0
c
      call sp_lbbl_prefetch(int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),
     + iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),xs,ps,psp,
     + .false.)
c
      icbbl=nbbl
c
      return
      end
      subroutine sp_prefetch(n,iwl,xw,isl,xs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(2),ps(2)
      real*8 pwp(2),psp(2)
      logical lfetch,lnext
      integer ibbl,nfetched
c
      n=nbbl
      if(lqmd) n=0
c
      nfetched=0
      if(ibget.lt.0) nbget=nbbl+1
c
      do 1 ibbl=1,nbbl
      call sp_ibbl_prefetch(ibbl,int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl),iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),
     + xs,ps,psp,.false.,lfetch,lnext)
      icbbl=ibbl
      if(lfetch) nfetched=nfetched+1
      if(nfetched.ge.nbget.and..not.lnext) return
 1    continue
c
      return
      end
      subroutine sp_prefetch_next(iwl,xw,isl,xs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(2),ps(2)
      real*8 pwp(2),psp(2)
      logical lfetch,lnext
      integer ibbl,nfetched
c
      if(icbbl.ge.nbbl) return
c
      nfetched=0
      do 1 ibbl=icbbl+1,nbbl
      call sp_ibbl_prefetch(ibbl,int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl),iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),
     + xs,ps,psp,.false.,lfetch,lnext)
      icbbl=ibbl
      if(lfetch) nfetched=nfetched+1
      if(nfetched.gt.0.and..not.lnext) return
 1    continue
c
      return
      end

      subroutine sp_prefetch_p(n,iwl,xw,pw,pwp,isl,xs,ps,psp)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
c
      n=nbbl
      if(lqmd) n=0
c
      call sp_lbbl_prefetch(int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),
     + iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),xs,ps,psp,
     + .true.)
c
      return
      end
      subroutine sp_lbbl_prefetch(lbbl,ipl,jpl,
     + iwl,iwlp,xw,pw,pwp,isl,islp,xs,ps,psp,lpp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer lbbl(mbbl,mbb2),ipl(mbox,mip2),jpl(mbox,mip2)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp
c
      integer ibbl,i,j,indexw,indexs,indx,jndx,nndexw,nndexs
      integer jproc,jcell
      integer jccell,iccell
      integer ih,il,jh,jl
      integer nsnon,nwnon
      integer jwfr,jwto,jsfr,jsto
      integer handle_is,handle_iw,handle_s,handle_w
      integer handle_ps,handle_pw,handle_psp,handle_pwp
      character*80 string
c
      jccell=-1
      indexw=nwmloc+1
      indexs=nsaloc+1
      nndexw=indexw
      nndexs=indexs
      do 1 ibbl=1,nbbl
c
      handle_iw=0
      handle_is=0
      handle_w=0
      handle_s=0
      handle_pw=0
      handle_ps=0
      handle_pwp=0
      handle_psp=0
c
      jproc=lbbl(ibbl,1)
      jcell=lbbl(ibbl,2)
c
      if(jproc.eq.me) then
c
      jccell=jcell
c
      else
c
      if(jcell.ne.jccell) then
c
      indexw=nndexw
      indexs=nndexs
c
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
c
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(indexw+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',indexw+nwnon)
      call ga_distribution(ga_iw,jproc,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,
     + iwl(indexw,1),mwm,handle_iw)
      else
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jl+npackw-1,
     + iwlp(indexw,1),mwm,handle_iw)
      endif
      call ga_distribution(ga_w,jproc,il,ih,jl,jh)
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(indexw,1,1),mwm,handle_w)
      if(lpp) then
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+18*mwa+2,
     + pw(indexw,1,1,1),mwm,handle_pw)
      if(lfree) call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,
     + jl+18*mwa+3,jl+30*mwa+2,pwp(indexw,1,1,1,1),mwm,handle_pwp)
      endif
      nndexw=indexw+nwnon
      endif
c
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(indexs+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (2)',indexs+nsnon)
      call ga_distribution(ga_is,jproc,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jh,
     + isl(indexs,1),msa,handle_is)
      else
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(indexs,1),msa,handle_is)
      endif
      call ga_distribution(ga_s,jproc,il,ih,jl,jh)
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(indexs,1),msa,handle_s)
      if(lpp) then
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(indexs,1,1),msa,handle_ps)
      if(lfree) call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(indexs,1,1,1),msa,handle_psp)
      endif
      nndexs=indexs+nsnon
      endif
c
      endif
c
      jccell=jcell
c
      endif
c
      lbbl(ibbl, 5)=indexw
      lbbl(ibbl, 6)=indexs
      lbbl(ibbl,15)=handle_iw
      lbbl(ibbl,16)=handle_w
      lbbl(ibbl,17)=handle_is
      lbbl(ibbl,18)=handle_s
      lbbl(ibbl,19)=handle_pw
      lbbl(ibbl,20)=handle_ps
      lbbl(ibbl,21)=handle_pwp
      lbbl(ibbl,22)=handle_psp
c
    1 continue
c
      if(iand(idebug,16).eq.16) then
      write(lfndbg,2000)
 2000 format(' Cell-cell list lbbl',/)
      write(lfndbg,2001) ((lbbl(i,j),j=1,18),i=1,nbbl)
 2001 format(18i6)
      endif
c
      return
      end
      subroutine sp_ibbl_prefetch(ibbl,lbbl,ipl,jpl,
     + iwl,iwlp,xw,pw,pwp,isl,islp,xs,ps,psp,lpp,lfetch,lnext)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl
      integer lbbl(mbbl,mbb2),ipl(mbox,mip2),jpl(mbox,mip2)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp,lfetch,lnext
c
      integer i,j,indexw,indexs,indx,jndx,nndexw,nndexs
      integer jproc,jcell
      integer jccell,iccell
      integer ih,il,jh,jl
      integer nsnon,nwnon
      integer jwfr,jwto,jsfr,jsto
      integer handle_is,handle_iw,handle_s,handle_w
      integer handle_ps,handle_pw,handle_psp,handle_pwp
      character*80 string
c
      lfetch=.false.
c
      if(ibbl.eq.1) then
      jccell=-1
      else
      jccell=lbbl(ibbl-1,2)
      endif
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
c
      handle_iw=0
      handle_is=0
      handle_w=0
      handle_s=0
      handle_pw=0
      handle_ps=0
      handle_pwp=0
      handle_psp=0
c
      jproc=lbbl(ibbl,1)
      jcell=lbbl(ibbl,2)
c
      if(jproc.eq.me) then
c
      jccell=jcell
c
      else
c
      if(jcell.ne.jccell) then
c
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
c
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(indexw+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',indexw+nwnon)
      call ga_distribution(ga_iw,jproc,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,
     + iwl(indexw,1),mwm,handle_iw)
      else
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jl+npackw-1,
     + iwlp(indexw,1),mwm,handle_iw)
      endif
      call ga_distribution(ga_w,jproc,il,ih,jl,jh)
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(indexw,1,1),mwm,handle_w)
      if(lpp) then
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+18*mwa+2,
     + pw(indexw,1,1,1),mwm,handle_pw)
      if(lfree) call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,
     + jl+18*mwa+3,jl+30*mwa+2,pwp(indexw,1,1,1,1),mwm,handle_pwp)
      endif
      lfetch=.true.
      endif
c
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(indexs+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (2)',indexs+nsnon)
      call ga_distribution(ga_is,jproc,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jh,
     + isl(indexs,1),msa,handle_is)
      else
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(indexs,1),msa,handle_is)
      endif
      call ga_distribution(ga_s,jproc,il,ih,jl,jh)
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(indexs,1),msa,handle_s)
      if(lpp) then
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(indexs,1,1),msa,handle_ps)
      if(lfree) call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(indexs,1,1,1),msa,handle_psp)
      endif
      lfetch=.true.
      endif
c
      endif
c
      jccell=jcell
c
      endif
c
      lbbl(ibbl,15)=handle_iw
      lbbl(ibbl,16)=handle_w
      lbbl(ibbl,17)=handle_is
      lbbl(ibbl,18)=handle_s
      lbbl(ibbl,19)=handle_pw
      lbbl(ibbl,20)=handle_ps
      lbbl(ibbl,21)=handle_pwp
      lbbl(ibbl,22)=handle_psp
c
      lnext=.false.
      if(ibbl.lt.nbbl) then
      lnext=lbbl(ibbl+1,2).eq.lbbl(ibbl,2)
      endif
c
      return
      end
      subroutine sp_nbwait(ibbl,lnew,lhandl,lself,lpbcs,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,isl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,lhandl
      integer isl(msa,mis2),iwl(mwm,miw2)
      logical lself,lpbcs,lnew
c
      lpbcs=lpbc9
      call sp_nbwait2(ibbl,lnew,int_mb(i_bb),lhandl,lself,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,
     + iwl,int_mb(i_packw),isl,int_mb(i_pack))
c
      return
      end
      subroutine sp_nbwait2(ibbl,lnew,lbbl,lhandl,lself,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,iwlp,isl,islp)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,lhandl
      integer lbbl(mbbl,mbb2)
      integer isl(msa,mis2),islp(msa,npack)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      logical lself,lnew
      integer indexs,indexw
      integer handle_iw,handle_w,handle_is,handle_s
      integer handle_pw,handle_ps,handle_pwp,handle_psp
      integer i
c
      lnew=.false.
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
      iwfr=lbbl(ibbl, 7)
      iwto=lbbl(ibbl, 8)
      isfr=lbbl(ibbl, 9)
      isto=lbbl(ibbl,10)
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
      if(lbbl(ibbl,1).ne.me) then
      if(lbbl(ibbl,11).gt.0) then
      jwfr=lbbl(ibbl,5)
      jwto=jwfr+lbbl(ibbl,12)-lbbl(ibbl,11)
      endif
      if(lbbl(ibbl,13).gt.0) then
      jsfr=lbbl(ibbl,6)
      jsto=jsfr+lbbl(ibbl,14)-lbbl(ibbl,13)
      endif
      endif
      lself=lbbl(ibbl,2).eq.lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
c
      handle_iw=lbbl(ibbl,15)
      handle_w=lbbl(ibbl,16)
      handle_is=lbbl(ibbl,17)
      handle_s=lbbl(ibbl,18)
      handle_pw=lbbl(ibbl,19)
      handle_ps=lbbl(ibbl,20)
      handle_pwp=lbbl(ibbl,21)
      handle_psp=lbbl(ibbl,22)
c
      if(handle_iw.ne.0.or.handle_w.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_iw)
      call ga_nbwait(handle_w)
      if(npackw.ne.0) then
      call sp_unpackw(jwto-jwfr+1,iwl(indexw,1),iwlp(indexw,1))
      endif
      endif
      if(handle_is.ne.0.or.handle_s.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_is)
      call ga_nbwait(handle_s)
      if(npackw.ne.0) then
      call sp_unpack(jsto-jsfr+1,isl(indexs,1),islp(indexs,1))
      endif
      endif
      if(handle_pw.ne.0.or.handle_ps.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_pw)
      call ga_nbwait(handle_ps)
      endif
      if(handle_pwp.ne.0.or.handle_psp.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_pwp)
      call ga_nbwait(handle_psp)
      endif
c
      return
      end
      subroutine sp_nbwaitf()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      call sp_nbwait3(int_mb(i_bb))
c
      return
      end
      subroutine sp_nbwait3(lbbl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl
      integer lbbl(mbbl,mbb2)
      integer handle_iw,handle_fw,handle_is,handle_fs
      integer handle_fwl,handle_fsl
c
      do 1 ibbl=1,nbbl
c
      handle_iw=lbbl(ibbl,15)
      handle_fw=lbbl(ibbl,16)
      handle_is=lbbl(ibbl,17)
      handle_fs=lbbl(ibbl,18)
      handle_fwl=lbbl(ibbl,19)
      handle_fsl=lbbl(ibbl,20)
c
      if(handle_iw.ne.0) then
      call ga_nbwait(handle_iw)
      endif
      if(handle_fw.ne.0) then
      call ga_nbwait(handle_fw)
      endif
      if(handle_is.ne.0) then
      call ga_nbwait(handle_is)
      endif
      if(handle_fs.ne.0) then
      call ga_nbwait(handle_fs)
      endif
      if(handle_fwl.ne.0) then
      call ga_nbwait(handle_fwl)
      endif
      if(handle_fsl.ne.0) then
      call ga_nbwait(handle_fsl)
      endif
c
    1 continue
c
      return
      end
      subroutine sp_lbbl_indices()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      call sp_lbbl_ind(int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl))
c
      return
      end
      subroutine sp_lbbl_ind(lbbl,ipl,jpl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer lbbl(mbbl,mbb2),ipl(mbox,mip2),jpl(mbox,mip2)
c
      integer ibbl,jproc,jcell,icell,jcproc,jccell,iccell
      integer i,j,indx,jndx,il,ih,jl,jh,indexw,indexs,nndexw,nndexs
      character*80 string
c
      jcproc=-1
      jccell=-1
      iccell=-1
      indexw=nwmloc+1
      indexs=nsaloc+1
      nndexw=indexw
      nndexs=indexs
c
c      if(iand(idebug,16).eq.16) then
c      write(lfndbg,3000) ((lbbl(ibbl,i),i=1,3),ibbl=1,nbbl)
c 3000 format(3i5)
c      endif
      do 1 ibbl=1,nbbl
c
      jproc=lbbl(ibbl,1)
      jcell=lbbl(ibbl,2)
      icell=lbbl(ibbl,3)
      lbbl(ibbl, 5)=0
      lbbl(ibbl, 6)=0
      if(iand(idebug,16).eq.16) then
      write(lfndbg,222) ibbl,nbbl,jproc,jcell,icell,jccell,iccell
  222 format(7i5)
      endif
      do 234 j=1,3
  234 continue
c
c     get indices for icell on the local node
c
      if(iccell.ne.icell) then
      indx=0
      do 2 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.icell) indx=30+i
    2 continue
      if(indx.eq.0) then
      do 232 j=1,3
  232 continue
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Cells on proc',i5,' :',t20,20i5,/,(t20,20i5))
      write(string,1001) icell,me
 1001 format('SP0001: Could not find icell',i5,' on proc',i5)
      call md_abort(string,0)
      endif
      iccell=icell
      endif
      lbbl(ibbl, 7)=ipl(indx,2)
      lbbl(ibbl, 8)=ipl(indx,3)
      lbbl(ibbl, 9)=ipl(indx,4)
      lbbl(ibbl,10)=ipl(indx,5)
c
c     get indices for jcell if on the local node
c
      if(jproc.eq.me) then
c
      if(jccell.ne.jcell) then
      jndx=0
      do 3 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.jcell) jndx=30+i
    3 continue
      if(jndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
      write(string,1002) jcell,me
 1002 format('SP0001: Could not find jcell',i5,' on proc',i5)
      call md_abort(string,0)
      endif
      jccell=jcell
      endif
      lbbl(ibbl,11)=ipl(jndx,2)
      lbbl(ibbl,12)=ipl(jndx,3)
      lbbl(ibbl,13)=ipl(jndx,4)
      lbbl(ibbl,14)=ipl(jndx,5)
c
c     get indices for jcell if on the remote node
c
      else
c
      if(jproc.ne.jcproc) then
      call ga_distribution(ga_ip,jproc,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,jpl,mbox)
      jcproc=jproc
      endif
c
      if(jcell.ne.jccell) then
      jndx=0
      do 4 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jcell) jndx=30+i
    4 continue
      if(jndx.eq.0) then
      if(iand(idebug,16).eq.16) then
      write(lfndbg,1003) jcell,jproc
 1003 format('Could not find cell',i5,' on proc',i5,' in sp_lbbl')
      write(lfndbg,1004) jproc,(jpl(30+i,1),i=1,jpl(1,1))
 1004 format('Cell list obtained from proc',i5,' is',/,(20i5))
      call util_flush(lfndbg)
      endif
c
      write(*,1005) jproc,(jpl(30+i,1),i=1,jpl(1,1))
 1005 format(/,'Cells on proc',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1006) jcell,jproc
 1006 format('SP0002: Could not find remote jcell',i5,' on proc',i5)
      call md_abort(string,me)
      endif
      jccell=jcell
c
      if(nbget.ne.0) then
      indexw=nndexw
      indexs=nndexs
      if(jpl(jndx,2).gt.0) nndexw=indexw+jpl(jndx,3)-jpl(jndx,2)+1
      if(jpl(jndx,4).gt.0) nndexs=indexs+jpl(jndx,5)-jpl(jndx,4)+1
      endif
c
      endif
c
      lbbl(ibbl, 5)=indexw
      lbbl(ibbl, 6)=indexs
      lbbl(ibbl,11)=jpl(jndx,2)
      lbbl(ibbl,12)=jpl(jndx,3)
      lbbl(ibbl,13)=jpl(jndx,4)
      lbbl(ibbl,14)=jpl(jndx,5)
c
      endif
c
    1 continue
c
      if(nndexw.gt.mwm) then
      write(string,1007) mwm,nndexw
 1007 format('Increase dimension mwm from ',i7,' to ',i7)
      call md_abort(string,me)
      endif
      if(nndexs.gt.msa) then
      write(string,1008) msa,nndexs
 1008 format('Increase dimension msa from ',i7,' to ',i7)
      call md_abort(string,me)
      endif
c
      return
      end
      subroutine sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      logical lself
c
      call sp_gthdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto,
     + int_mb(i_bb),int_mb(i_ipl))
c
      return
      end
      subroutine sp_gthdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto,
     + lbbl,ipl)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2)
      logical lself
c
      integer ndx,i,ibox
      character*255 string
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
      lself=ibox.eq.lbbl(ibbl,2)
c
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1001) ibox,me
 1001 format('SP0001: Could not find box',i5,' on node',i5)
      call md_abort(string,0)
      endif
c
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
c
      return
      end
      subroutine sp_getxpbl(ibbl,lhandl,
     + iwl,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,xs,ps,psp,isfr,isto,jsfr,jsto,lself,lpbcs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lself,lpbcs
c
      call sp_gtxpbl(ibbl,lhandl,
     + iwl,int_mb(i_packw),xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,int_mb(i_pack),xs,ps,psp,isfr,isto,jsfr,jsto,
     + int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),lself,.true.)
c
      lpbcs=lpbc9
c
      return
      end
      logical function sp_local(ibbl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl
      logical local
c
      call sp_gtlocal(ibbl,local,int_mb(i_bb))
c
      sp_local=local
c
      return
      end
      subroutine sp_gtlocal(ibbl,local,lbbl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl
      logical local
      integer lbbl(mbbl,mbb2)
c
      local=lbbl(ibbl,1).eq.me
c
      return
      end
      subroutine sp_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lself,lpbcs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      logical lself,lpbcs
      real*8 pw(2),pwp(2),ps(2),psp(2)
c
      call sp_gtxpbl(ibbl,lhandl,
     + iwl,int_mb(i_packw),xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,int_mb(i_pack),xs,ps,psp,isfr,isto,jsfr,jsto,
     + int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),lself,.false.)
c
      lpbcs=lpbc9
c
      return
      end
      subroutine sp_getda(jnode,jsfr,jsto,isl,isga,isq3,xs,ndxda)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer jnode,jsfr,jsto,ndxda
      integer isga(msa),isq3(msa),isl(msa,mis2)
      real*8 xs(msa,3)
c
      integer il,ih,jl,jh
      integer l_islt,i_islt
c
      if(.not.ma_push_get(mt_int,msa*mis2,'islt',l_islt,i_islt))
     + call md_abort('Failed to allocate memory for islt',0)
      call sp_gtda(jnode,jsfr,jsto,int_mb(i_islt),int_mb(i_pack),
     + isga,isq3,xs,ndxda)
      if(.not.ma_pop_stack(l_islt))
     + call md_abort('Failed to deallocate memory for islt',0)
c
      return
      end
      subroutine sp_gtda(jnode,jsfr,jsto,isl,islp,isga,isq3,xs,ndxda)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer jnode,jsfr,jsto,ndxda
      integer isl(msa,mis2),islp(msa,npack)
      integer isga(msa),isq3(msa)
      real*8 xs(msa,3)
c
      integer i,il,ih,jl,jh
c
      if(npack.eq.0) then
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl,
     + isga(ndxda+1),msa)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl+lsct3-1,jl+lsct3-1,
     + isq3(ndxda+1),msa)
      else
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(nsaloc+1,1),msa)
      call sp_unpack(jsto-jsfr+1,isl(nsaloc+1,1),islp(nsaloc+1,1))
      do 1 i=1,jsto-jsfr
      isga(ndxda+i)=isl(nsaloc+i,lsgan)
      isq3(ndxda+i)=isl(nsaloc+i,lsct3)
    1 continue
      endif
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(ndxda+1,1),msa)
c
      return
      end
      subroutine sp_gtxpbl(ibbl,lhandl,
     + iwl,iwlp,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,islp,xs,ps,psp,isfr,isto,jsfr,jsto,
     + lbbl,ipl,jpl,lself,lpp)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2),jpl(mbox,mip2)
      logical lself,lpp
c
      integer i,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      character*80 string
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
c
      lself=ibox.eq.jbox
c
      if(ibbl.eq.1) then
      jcnode=-1
      jcbox=-1
      icbox=-1
      endif
c
      if(ibox.ne.icbox) then
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1001) ibox,me
 1001 format('SP0001: Could not find local ibox',i5,' on node',i5)
      call md_abort(string,me)
      endif
c
      icbox=ibox
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
      endif
c
      if(jnode.eq.me) then
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 2 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.jbox) ndx=30+i
    2 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
      write(string,1002) jbox,me
 1002 format('SP0001: Could not find local jbox',i5,' on node',i5)
      call md_abort(string,me)
      endif
c
      jcbox=jbox
      jwfr=ipl(ndx,2)
      jwto=ipl(ndx,3)
      jsfr=ipl(ndx,4)
      jsto=ipl(ndx,5)
      endif
c
      else
c
      if(jnode.ne.jcnode) then
      call ga_distribution(ga_ip,jnode,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,jpl,mbox)
      jcnode=jnode
      endif
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0) then
      if(iand(idebug,4).eq.4) then
      write(lfndbg,8000) jbox,jnode
 8000 format('Could not find box',i5,' on node',i5,' in sp_gtxpbl')
      write(lfndbg,8001) jnode,(jpl(30+i,1),i=1,jpl(1,1))
 8001 format('Box list obtained from node',i5,' is',/,(20i5))
      call util_flush(lfndbg)
      endif
c
      write(*,1003) jnode,(jpl(30+i,1),i=1,jpl(1,1))
 1003 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1004) jbox,jnode
 1004 format('SP0002: Could not find remote jbox',i5,' on node',i5)
      call md_abort(string,me)
c
      endif
      jcbox=jbox
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_iw,jnode,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_get(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,iwl(nwmloc+1,1),mwm)
      else
      call ga_get(ga_iw,il+jwfr-1,il+jwto-1,jl,jl+npackw-1,
     + iwlp(nwmloc+1,1),mwm)
      call sp_unpackw(jwto-jwfr+1,iwl(nwmloc+1,1),iwlp(nwmloc+1,1))
      endif
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_get(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(nwmloc+1,1,1),mwm)
      if(lpp) then
      call ga_get(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+18*mwa+2,
     + pw(nwmloc+1,1,1,1),mwm)
      if(lfree) call ga_get(ga_w,il+jwfr-1,il+jwto-1,
     + jl+18*mwa+3,jl+30*mwa+2,pwp(nwmloc+1,1,1,1,1),mwm)
      endif
      jwfr=nwmloc+1
      jwto=nwmloc+1+jpl(ndx,3)-jpl(ndx,2)
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (3)',nsaloc+nsnon)
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jh,isl(nsaloc+1,1),msa)
      else
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(nsaloc+1,1),msa)
      call sp_unpack(jsto-jsfr+1,isl(nsaloc+1,1),islp(nsaloc+1,1))
      endif
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(nsaloc+1,1),msa)
      if(lpp) then
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(nsaloc+1,1,1),msa)
      if(lfree) call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(nsaloc+1,1,1,1),msa)
      endif
      jsfr=nsaloc+1
      jsto=nsaloc+1+jpl(ndx,5)-jpl(ndx,4)
      endif
      endif
c
      endif
c
      return
      end
      subroutine sp_accpbl(ibbl,lhandl,pw,pwp,ps,psp,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call sp_acpbl(ibbl,lhandl,pw,pwp,ps,psp,int_mb(i_bb),
     + int_mb(i_jpl))
c
      return
      end
      subroutine sp_acpbl(ibbl,lhandl,pw,pwp,ps,psp,lbbl,jpl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      integer lbbl(mbbl,mbb2)
      integer jpl(mbox,mip2)
c
      integer i,j,k,l,m,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      integer jwfr,jwto,jsfr,jsto
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lbbl(ibbl,4)=lhandl
c
      if(jnode.eq.me) return
      if(ibbl.lt.nbbl) then
      if(lbbl(ibbl+1,1).eq.jnode.and.lbbl(ibbl+1,2).eq.jbox) return
      endif
c
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0) then
      if(iand(idebug,4).eq.4) then
      write(lfndbg,8000) jbox,jnode
 8000 format('Could not find box',i5,' on node',i5,' in sp_acpbl')
      write(lfndbg,8001) jnode,(jpl(30+i,1),i=1,jpl(1,1))
 8001 format('Box list obtained from node',i5,' is',/,(20i5))
      call util_flush(lfndbg)
      endif
      call md_abort('Remote processor cell not found (0) ',0)
      endif
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+15*mwa+2,
     + pw(nwmloc+1,1,1,1),mwm,one)
      if(lfree) call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+15*mwa+3,
     + jl+21*mwa+2,pwp(nwmloc+1,1,1,1,1),mwm,one)
c      if(llong) call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
c     + jl+12*mwa+2,fw(nwmloc+1,1,1,2),mwm,one)
c      if(lpair) then
c      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
c      call ga_acc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,iwz(nwmloc+1),mwm,1)
c      endif
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (4)',nsaloc+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(nsaloc+1,1,1),msa,one)
      if(lfree) call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(nsaloc+1,1,1,1),msa,one)
c      if(lpair) then
c      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
c      call ga_acc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,isz(nsaloc+1),msa,1)
c      endif
      endif
c
      m=1
      if(llong) m=2
c
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=nwmloc+1,mwm
      pw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
c
      do 7 j=1,3
      do 8 i=nsaloc+1,msa
      ps(i,j,l)=zero
    8 continue
    7 continue
    9 continue
c
      return
      end
      subroutine sp_accfbl(ibbl,lhandl,fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call sp_acfbl(ibbl,lhandl,fw,fs,int_mb(i_bb),int_mb(i_jpl),
     + lpair,iwz,isz)
c
      return
      end
      subroutine sp_acfbl(ibbl,lhandl,fw,fs,lbbl,jpl,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer lbbl(mbbl,mbb2)
      integer jpl(mbox,mip2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,j,k,l,m,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      integer jwfr,jwto,jsfr,jsto
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lbbl(ibbl,4)=lhandl
c
      if(jnode.eq.me) return
      if(ibbl.lt.nbbl) then
      if(lbbl(ibbl+1,1).eq.jnode.and.lbbl(ibbl+1,2).eq.jbox) return
      endif
c
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0)
     + call md_abort('Remote processor cell not found (2) ',0)
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+6*mwa+3,jl+9*mwa+2,
     + fw(nwmloc+1,1,1,1),mwm,one)
      if(llong) call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
     + jl+12*mwa+2,fw(nwmloc+1,1,1,2),mwm,one)
      if(lpair) then
      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
      call ga_acc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,iwz(nwmloc+1),mwm,1)
      endif
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (5)',nsaloc+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+6,jl+8,
     + fs(nsaloc+1,1,1),msa,one)
      if(llong) call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+9,jl+11,
     + fs(nsaloc+1,1,2),msa,one)
      if(lpair) then
      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
      call ga_acc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,isz(nsaloc+1),msa,1)
      endif
      endif
c
      m=1
      if(llong) m=2
c
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=nwmloc+1,mwm
      fw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
c
      do 7 j=1,3
      do 8 i=nsaloc+1,msa
      fs(i,j,l)=zero
    8 continue
    7 continue
    9 continue
c
      do 10 i=nwmloc+1,mwm
      iwz(i)=0
   10 continue
      do 11 i=nsaloc+1,msa
      isz(i)=0
   11 continue
c
      return
      end
      subroutine sp_nbaccfbl(ibbl,lhandl,fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call sp_nbacfbl(ibbl,lhandl,fw,fs,int_mb(i_bb),int_mb(i_jpl),
     + lpair,iwz,isz)
c
      return
      end
      subroutine sp_nbacfbl(ibbl,lhandl,fw,fs,lbbl,jpl,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer lbbl(mbbl,mbb2)
      integer jpl(mbox,mip2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,j,k,l,m,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      integer jwfr,jwto,jsfr,jsto,indexw,indexs
      integer handle_iw,handle_fw,handle_is,handle_fs
      integer handle_fwl,handle_fsl
c
      handle_iw=0
      handle_is=0
      handle_fw=0
      handle_fs=0
      handle_fwl=0
      handle_fsl=0
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lbbl(ibbl,4)=lhandl
c
      if(jnode.eq.me) return
      if(ibbl.lt.nbbl) then
      if(lbbl(ibbl+1,1).eq.jnode.and.lbbl(ibbl+1,2).eq.jbox) return
      endif
c
c      ndx=0
c      do 3 i=1,jpl(1,1)
c      if(jpl(30+i,1).eq.jbox) ndx=30+i
c    3 continue
c      if(ndx.eq.0)
c     + call md_abort('Remote processor cell not found (3) ',0)
c
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
c
c      if(lbbl(ibbl,1).ne.me) then
c      if(lbbl(ibbl,12).gt.0) then
c      jwfr=lbbl(ibbl,5)
c      jwto=jwfr+lbbl(ibbl,12)-lbbl(ibbl,11)
c      endif
c      if(lbbl(ibbl,14).gt.0) then
c      jsfr=lbbl(ibbl,6)
c      jsto=jsfr+lbbl(ibbl,14)-lbbl(ibbl,13)
c      endif
c      endif
c
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
c
      m=1
      if(llong) m=2
c
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(indexw+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',indexw+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_w,il+jwfr-1,il+jwto-1,jl+6*mwa+3,jl+9*mwa+2,
     + fw(indexw,1,1,1),mwm,one,handle_fw)
      if(llong) call ga_nbacc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
     + jl+12*mwa+2,fw(indexw,1,1,2),mwm,one,handle_fwl)
      if(lpair) then
      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,
     + iwz(indexw),mwm,1,handle_iw)
      endif
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=indexw,indexw+nwnon-1
      fw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
    9 continue
      do 10 i=indexw,indexw+nwnon-1
      iwz(i)=0
   10 continue
      endif
c
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(indexs+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (5)',indexs+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_s,il+jsfr-1,il+jsto-1,jl+6,jl+8,
     + fs(indexs,1,1),msa,one,handle_fs)
      if(llong) call ga_nbacc(ga_s,il+jsfr-1,il+jsto-1,jl+9,jl+11,
     + fs(indexs,1,2),msa,one,handle_fsl)
      if(lpair) then
      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,
     + isz(indexs),msa,1,handle_is)
      endif
      do 12 l=1,m
      do 7 j=1,3
      do 8 i=indexs,indexs+nsnon-1
      fs(i,j,l)=zero
    8 continue
    7 continue
   12 continue
      do 11 i=indexs,indexs+nsnon-1
      isz(i)=0
   11 continue
      endif
c
      lbbl(ibbl,15)=handle_iw
      lbbl(ibbl,16)=handle_fw
      lbbl(ibbl,17)=handle_is
      lbbl(ibbl,18)=handle_fs
      lbbl(ibbl,19)=handle_fwl
      lbbl(ibbl,20)=handle_fsl
c
      return
      end
      subroutine sp_bscale(boxsiz)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 boxsiz(maxbox,3)
c
      integer ibx,iby,ibz
      real*8 boxscl
c
c     scale the subbox boundaries
c
      boxscl=box(1)/boxsiz(nbx,1)
      do 1 ibx=1,nbx
      boxsiz(ibx,1)=boxscl*boxsiz(ibx,1)
    1 continue
      boxscl=box(2)/boxsiz(nby,2)
      do 2 iby=1,nby
      boxsiz(iby,2)=boxscl*boxsiz(iby,2)
    2 continue
      boxscl=box(3)/boxsiz(nbz,3)
      do 3 ibz=1,nbz
      boxsiz(ibz,3)=boxscl*boxsiz(ibz,3)
    3 continue
c
      return
      end
      subroutine sp_balanc(stime,syntim,waltim,lpsyn)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      logical sp_diffbb
      external sp_diffbb
c
      real*8 stime,syntim,waltim
      logical lpsyn
c
      integer l_dlb,i_dlb,i_nod,l_nod,i_nen,l_nen,i_lst,l_lst
c
c      if(loadb.eq.0) return
c
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in sp_balanc at 1',me)
      if(.not.ma_push_get(mt_int,np,'nod',l_nod,i_nod))
     + call md_abort('Failed to allocate nod',me)
      if(.not.ma_push_get(mt_int,2*(np+1),'lst',l_lst,i_lst))
     + call md_abort('Failed to allocate lst',me)
      if(.not.ma_push_get(mt_int,mbbl*np,'nen',l_nen,i_nen))
     + call md_abort('Failed to allocate nen',me)
      if(.not.ma_push_get(mt_dbl,np+1,'dlb',l_dlb,i_dlb))
     + call md_abort('Failed to allocate dlb',me)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in sp_balanc at 2',me)
c
      call sp_dldbal(stime,syntim,waltim,int_mb(i_nod),int_mb(i_lst),
     + int_mb(i_nen),int_mb(i_bb),dbl_mb(i_dlb),lpsyn,int_mb(i_bindex),
     + int_mb(i_buren),np+1)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in sp_balanc at 3',me)
c
      if(.not.ma_pop_stack(l_dlb))
     + call md_abort('Failed to deallocate dlb',me)
      if(.not.ma_pop_stack(l_nen))
     + call md_abort('Failed to deallocate nen',me)
      if(.not.ma_pop_stack(l_lst))
     + call md_abort('Failed to deallocate lst',me)
      if(.not.ma_pop_stack(l_nod))
     + call md_abort('Failed to deallocate nod',me)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in sp_balanc at 4',me)
c
      if(sp_diffbb(dbl_mb(i_boxs),int_mb(i_rng))) then
      call sp_numbb(int_mb(i_iown),dbl_mb(i_boxs))
      call sp_listbb(int_mb(i_iown),dbl_mb(i_boxs),int_mb(i_bb))
      endif
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in sp_balanc at 5',me)
c
      return
      end
      subroutine sp_dldbal(stime,syntim,waltim,nod,nlst,nen,lbbl,dlb,
     + lpsyn,ibindx,iburen,npp)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer sp_btop
      external sp_btop
c
      integer npp
      integer nod(np),nlst(npp,2),nen(np,mbbl)
      integer lbbl(mbbl,mbb2)
      real*8 stime,syntim,synsum,waltim,dlb(*),factor,facx,facy,facz
      logical lpsyn
      logical lbsize
c
      integer i,j,k,n,ibusy,least,iproc,icell,jcell,itemp,jtemp
      integer intsiz,ibindx(np),iburen(np,27,2)
      logical lnb(27)
      integer nfr,nto,node
c
c     no load balancing needed if run on single processor
c
      if(np.eq.1) return
c
c     collect the synchronization times converted into integers
c
      do 1 i=1,np
      nlst(i,1)=0
      nlst(i,2)=0
      nod(i)=i-1
    1 continue
      nlst(me+1,1)=int(1.0d6*syntim)
      nlst(npp,1)=int(1.0d6*waltim)
      nlst(npp,2)=0
c
      do 2 i=1,nbbl
      iproc=lbbl(i,1)
      if(iproc.ne.me) then
      n=ibindx(iproc+1)
      if(n.gt.0) nlst(me+1,2)=ior(nlst(me+1,2),2**(n-1))
      endif
    2 continue
c
c     get the size of an integer on this machine in bits
c
      intsiz=ma_sizeof(mt_int,1,mt_byte)
c
c      if(intsiz.le.32) then
      call ga_igop(msp_14,nlst,2*npp,'+')
c      else
c      do 3 i=1,np
c      nlst(i,1)=ishft(nlst(i,1),32)+nlst(i,2)
c    3 continue
c      call ga_igop(msp_14,nlst,np,'+')
c      do 4 i=1,np
c      nlst(i,2)=iand(nlst(i,1),8589934591)
c      nlst(i,1)=ishft(nlst(i,1),-32)
c    4 continue
c      endif
c
      synsum=0.0d0
      do 5 i=1,np
      dlb(i)=dble(nlst(i,1))*1.0d-6
      synsum=synsum+dlb(i)
    5 continue
      waltim=(dble(nlst(npp,1))*1.0d-6/dble(np))
c
      if(nldup.lt.0) then
      tslow=synsum
      nldup=nldup+1
      else
      if(synsum.gt.tslow) then
      nldup=nldup+1
      else
      tslow=synsum
      nldup=0
      endif
      endif
c
c     if the accumulated synchronization time is less than 0.0001 sec
c     no load balancing is needed
c
      if(synsum.lt.1.0d-4) return
c
c     determine the busiest processor
c     if the busiest processor has no off-processor box-box pairs
c     load balancing will be done by resizing
c
      ibusy=1
      do 6 i=2,np
      if(abs(dlb(i)).lt.abs(dlb(ibusy))) ibusy=i
    6 continue
c
c     if the busiest processor has no off-processor pairs in the
c     cell-cell list loadbalancing will be done by resizing
c
      lbsize=nlst(ibusy,2).eq.0
      if(nldup.ge.lbpair) lbsize=.true.
c
c     order node list
c
      do 7 i=1,np-1
      do 8 j=i+1,np
      if(dlb(nod(i)+1).gt.dlb(nod(j)+1)) then
      n=nod(i)
      nod(i)=nod(j)
      nod(j)=n
      endif
    8 continue
    7 continue
c
c     load balancing based on cell resizing
c
c     1  if specified as only method in input
c     2. if busiest processor has no off-processor cell-cell pairs (lbsize true)
c     3. if busiest alternates between two processors (lpipo true)
c
      if(loadb.eq.1.or.(loadb.eq.3.and.(lbsize.or.lpipo))) then
      factor=(1.0d0-factld*(synsum/dble(np)-dlb(ibusy))/waltim)**third
      facx=1.0d0
      facy=1.0d0
      facz=1.0d0
      call sp_resize(ibusy-1,int_mb(i_iown),dbl_mb(i_boxs),factor,
     + facx,facy,facz)
      if(me.eq.0.and.lpsyn) then
      write(lfnsyn,3000) 1,stime,synsum,waltim,ibusy-1,
     + factor,facx,facy,facz,lpipo
 3000 format('synchronization',/,i5,3f12.6,/,i5,4f12.6,4x,l1)
      do 9 i=1,np
      write(lfnsyn,3001) nod(i),dlb(nod(i)+1)
 3001 format(i5,f12.6)
    9 continue
      call sp_wrtbxsz(lfnsyn,dbl_mb(i_boxs))
      endif
      ipairf=-3
      ipairt=-4
      lpipo=.false.
      tsyncp=synsum
      nldup=0
      tslow=synsum
      return
      endif
c
c     cascade implementation
c
      if(me.eq.0) then
c
      do 12 i=1,np
      nen(i,1)=-1
   12 continue
c
c     for each processor find least busy neighbor in cell-cell list
c
      do 10 i=1,np
      node=nod(i)
      n=nlst(node+1,2)
      do 11 j=1,27
      lnb(j)=iand(n,1).eq.1
      n=ishft(n,-1)
   11 continue
      least=-1
      do 13 j=1,27
      iproc=iburen(node+1,j,1)
      if(nen(iproc+1,1).eq.-1) then
      if(iproc.ne.node.and.lnb(j)) then
      if(least.ge.0) then
      if(dlb(iproc+1).gt.dlb(least+1)) least=iproc
      else
      if(dlb(iproc+1).gt.dlb(node+1)) least=iproc
      endif
      endif
      endif
   13 continue
c
c     nen(i,1) contains processor id that processor i-1 will receive from
c
      if(least.ge.0) then
      if(nen(least+1,1).lt.0) then
      nen(least+1,1)=node
      endif
      endif
c
   10 continue
c
      endif
c
c
c     broadcast the list
c
      call ga_brdcst(msp_16,nen(1,1),np*intsiz,0)
c
c     determine cell pair transfer list
c
c     the list contains the box-box pairs to be moved as follows:
c
c        nen(i,1) = ibox+1
c        nen(1,2) = jbox+1
c
c     the pair will move from iproc to jproc
c
      do 14 i=1,np
      if(me.eq.nen(i,1)) then
      do 15 j=1,nbbl
      if(lbbl(j,1).eq.i-1) then
      nen(i,1)=lbbl(j,3)+1
      nen(i,2)=lbbl(j,2)+1
      goto 14
      endif
   15 continue
      else
      nen(i,1)=0
      nen(i,2)=0
      endif
   14 continue
c
      call ga_igop(msp_17,nen,2*np,'+')
c
c     cascading
c
      call sp_cascad(nen,int_mb(i_bb))
c
      if(me.eq.0.and.lpsyn) then
      write(lfnsyn,3002) 2,stime,synsum,waltim,nldup,tslow
 3002 format('synchronization',/,i5,3f12.6,i5,f12.6)
      do 16 i=1,np
      write(lfnsyn,3003) nod(i),dlb(nod(i)+1)
 3003 format(i5,f12.6,i7,i5)
   16 continue
      endif
c
      return
      end
      subroutine sp_resize(ipmin,ibownr,boxsiz,factor,facx,facy,facz)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ipmin,ibownr(maxbox,3)
      real*8 factor,facx,facy,facz,boxsiz(maxbox,3)
c
      integer ibx,iby,ibz,ipminx,ipminy,ipminz
      real*8 bxi,bxr,byi,byr,bzi,bzr,rfact
c
c     load balancing based on resizing of node domains
c     ------------------------------------------------
c
      ipminx=mod(ipmin,npx)
      ipminy=mod((ipmin-ipminx)/npx,npy)
      ipminz=((ipmin-ipminx)/npx-ipminy)/npy
c
c     resizing boxes in x-dimension
c
      if(npx.gt.1.and.isload.eq.0) then
      do 7 ibx=nbx,2,-1
      boxsiz(ibx,1)=boxsiz(ibx,1)-boxsiz(ibx-1,1)
    7 continue
      bxi=zero
      bxr=zero
      do 8 ibx=1,nbx
      if(ibownr(ibx,1).eq.ipminx) then
      facx=1.0d0/boxsiz(ibx,1)
      boxsiz(ibx,1)=factor*boxsiz(ibx,1)
      if(boxsiz(ibx,1).lt.bxmin) boxsiz(ibx,1)=bxmin
      facx=facx*boxsiz(ibx,1)
      bxi=bxi+boxsiz(ibx,1)
      else
      bxr=bxr+boxsiz(ibx,1)
      endif
    8 continue
      rfact=(box(1)-bxi)/bxr
      do 9 ibx=1,nbx
      if(ibownr(ibx,1).ne.ipminx) boxsiz(ibx,1)=rfact*boxsiz(ibx,1)
      if(ibx.gt.1) boxsiz(ibx,1)=boxsiz(ibx,1)+boxsiz(ibx-1,1)
    9 continue
      endif
c
c     resizing boxes in y-dimension
c
      if(npy.gt.1.and.isload.eq.0) then
      do 10 iby=nby,2,-1
      boxsiz(iby,2)=boxsiz(iby,2)-boxsiz(iby-1,2)
   10 continue
      byi=zero
      byr=zero
      do 11 iby=1,nby
      if(ibownr(iby,2).eq.ipminy) then
      facy=1.0d0/boxsiz(iby,2)
      boxsiz(iby,2)=factor*boxsiz(iby,2)
      if(boxsiz(iby,2).lt.bymin) boxsiz(iby,2)=bymin
      facy=facy*boxsiz(iby,2)
      byi=byi+boxsiz(iby,2)
      else
      byr=byr+boxsiz(iby,2)
      endif
   11 continue
      rfact=(box(2)-byi)/byr
      do 12 iby=1,nby
      if(ibownr(iby,2).ne.ipminy) boxsiz(iby,2)=rfact*boxsiz(iby,2)
      if(iby.gt.1) boxsiz(iby,2)=boxsiz(iby,2)+boxsiz(iby-1,2)
   12 continue
      endif
c
c     resizing boxes in z-dimension
c
      if(npz.gt.1) then
      do 13 ibz=nbz,2,-1
      boxsiz(ibz,3)=boxsiz(ibz,3)-boxsiz(ibz-1,3)
   13 continue
      bzi=zero
      bzr=zero
      do 14 ibz=1,nbz
      if(ibownr(ibz,3).eq.ipminz) then
      facz=1.0d0/boxsiz(ibz,3)
      boxsiz(ibz,3)=factor*boxsiz(ibz,3)
      if(boxsiz(ibz,3).lt.bzmin) boxsiz(ibz,3)=bzmin
      facz=facz*boxsiz(ibz,3)
      bzi=bzi+boxsiz(ibz,3)
      else
      bzr=bzr+boxsiz(ibz,3)
      endif
   14 continue
      rfact=(box(3)-bzi)/bzr
      do 15 ibz=1,nbz
      if(ibownr(ibz,3).ne.ipminz) boxsiz(ibz,3)=rfact*boxsiz(ibz,3)
      if(ibz.gt.1) boxsiz(ibz,3)=boxsiz(ibz,3)+boxsiz(ibz-1,3)
   15 continue
      endif
c
      call ga_brdcst
     + (msp_19,boxsiz,ma_sizeof(mt_dbl,3*maxbox,mt_byte),0)
c
      return
      end
      subroutine sp_wrtbxsz(lfn,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 boxsiz(maxbox,3)
      integer lfn
c
      integer i
c
      write(lfn,1000) 1,nbx,(boxsiz(i,1),i=1,nbx)
      write(lfn,1000) 2,nby,(boxsiz(i,2),i=1,nby)
      write(lfn,1000) 3,nbz,(boxsiz(i,3),i=1,nbz)
 1000 format(i3,i5,5f10.5,/,(8x,5f10.5))
c
      return
      end
      logical function sp_diffbb(boxsiz,ibxrng)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibxrng(maxbox,3,2)
      real*8 boxsiz(maxbox,3)
      integer i,ib,j,k,num
      real*8 rd
c
c     temporarily skip for non-rectangular boxes
c
      if(nbxtyp.eq.1) then
      sp_diffbb=.true.
      return
      endif
c
c     initialize
c
      if(nbbdif.lt.0) then
      do 1 k=1,2
      do 2 j=1,3
      do 3 i=1,maxbox
      ibxrng(i,j,k)=-1
    3 continue
    2 continue
    1 continue
      endif
c
      nbbdif=0
c
      do 4 i=1,nbx
      num=0
      rd=zero
      do 5 j=1,nbx
      ib=i+1
      if(ib.gt.nbx) ib=ib-nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 6
      num=num+1
    5 continue
    6 continue
      if(num.ne.ibxrng(i,1,1)) then
      ibxrng(i,1,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 7 j=1,nbx
      ib=i-1
      if(ib.lt.1) ib=ib+nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 8
      num=num+1
    7 continue
    8 continue
      if(num.ne.ibxrng(i,1,2)) then
      ibxrng(i,1,2)=num
      nbbdif=nbbdif+1
      endif
    4 continue
c
      do 9 i=1,nby
      num=0
      rd=zero
      do 10 j=1,nby
      ib=i+1
      if(ib.gt.nby) ib=ib-nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 11
      num=num+1
   10 continue
   11 continue
      if(num.ne.ibxrng(i,2,1)) then
      ibxrng(i,2,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 12 j=1,nby
      ib=i-1
      if(ib.lt.1) ib=ib+nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 13
      num=num+1
   12 continue
   13 continue
      if(num.ne.ibxrng(i,2,2)) then
      ibxrng(i,2,2)=num
      nbbdif=nbbdif+1
      endif
    9 continue
c
      do 14 i=1,nbz
      num=0
      rd=zero
      do 15 j=1,nbz
      ib=i+1
      if(ib.gt.nbz) ib=ib-nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 16
      num=num+1
   15 continue
   16 continue
      if(num.ne.ibxrng(i,3,1)) then
      ibxrng(i,3,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 17 j=1,nbz
      ib=i-1
      if(ib.lt.1) ib=ib+nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 18
      num=num+1
   17 continue
   18 continue
      if(num.ne.ibxrng(i,3,2)) then
      ibxrng(i,3,2)=num
      nbbdif=nbbdif+1
      endif
   14 continue
c
      sp_diffbb=nbbdif.ne.0
c
      return
      end
      subroutine sp_cascad(ltran,lbbl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer sp_btop
      external sp_btop
c
      integer ltran(np,3),lbbl(mbbl,mbb2)
c
      integer i,j,k,ltemp
      integer icell,jcell,iproc,jproc
      logical lfirst
c
c     load balancing based on cascade box-list redistribution
c     -------------------------------------------------------
c
      do 1 i=1,np
      ltran(i,1)=ltran(i,1)-1
      ltran(i,2)=ltran(i,2)-1
      ltran(i,3)=-1
      if(ltran(i,2).ge.0) ltran(i,3)=sp_btop(ltran(i,2),int_mb(i_iown))
    1 continue
c
c     remove transfers involving double processors
c
      do 2 i=2,np
      if(ltran(i,1).ge.0) then
      iproc=sp_btop(ltran(i,1),int_mb(i_iown))
      jproc=sp_btop(ltran(i,2),int_mb(i_iown))
      do 3 j=1,i-1
      if(sp_btop(ltran(j,1),int_mb(i_iown)).eq.jproc.or.
     + sp_btop(ltran(j,2),int_mb(i_iown)).eq.jproc) then
      ltran(i,1)=-1
      ltran(i,2)=-1
      goto 2
      endif
    3 continue
      endif
    2 continue
c
c     process the transfer list
c
      lfirst=.true.
      do 4 i=1,np
c
      iproc=-1
      jproc=-1
      icell=ltran(i,1)
      jcell=ltran(i,2)
      if(icell.ge.0) then
      iproc=sp_btop(icell,int_mb(i_iown))
      jproc=sp_btop(jcell,int_mb(i_iown))
      if(lfirst) then
      lpipo=(iproc.eq.ipairt).and.(jproc.eq.ipairf)
      ipairf=iproc
      ipairt=jproc
      lfirst=.false.
      endif
c
c     remove pair
c
      if(iproc.eq.me) then
      k=0
      do 5 j=1,nbbl
      if(lbbl(j,2).ne.jcell.or.lbbl(j,3).ne.icell) then
      k=k+1
      lbbl(k,1)=lbbl(j,1)
      lbbl(k,2)=lbbl(j,2)
      lbbl(k,3)=lbbl(j,3)
      endif
    5 continue
      nbbl=k
      endif
c
c     add pair
c
      if(jproc.eq.me) then
      nbbl=nbbl+1
      lbbl(nbbl,1)=iproc
      lbbl(nbbl,2)=icell
      lbbl(nbbl,3)=jcell
      endif
c
      endif
    4 continue
c
c     order the new box-box list
c
      nbbloc=0
      do 6 i=1,nbbl-1
      do 7 j=i+1,nbbl
      if((lbbl(i,1).ne.me.and.lbbl(j,1).eq.me).or.
     + (lbbl(i,1).gt.lbbl(j,1).and.lbbl(i,1).ne.me).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).gt.lbbl(j,2)).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).eq.lbbl(j,2).and.
     + lbbl(i,3).gt.lbbl(j,3))) then
      ltemp=lbbl(i,1)
      lbbl(i,1)=lbbl(j,1)
      lbbl(j,1)=ltemp
      ltemp=lbbl(i,2)
      lbbl(i,2)=lbbl(j,2)
      lbbl(j,2)=ltemp
      ltemp=lbbl(i,3)
      lbbl(i,3)=lbbl(j,3)
      lbbl(j,3)=ltemp
      endif
    7 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    6 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      return
      end
