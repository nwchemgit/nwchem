      subroutine sp_start(lfntop,filtop,lfnrst,filrst,
     + lsyn,filsyn,nsyn,
     + rsht,rlng,rinp,
     + npxi,npyi,npzi,nbxi,nbyi,nbzi,
     + npbt,nbxt,boxt,vlatt,lpbc,
     + nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai,
     + ldb,lbp,fld,lpol,lfre,temp,tempw,temps,lqd,iguidi,
     + ldbg,idbg,prjct,mbbi,nseri,isld,irset)
c
c
c $Id: sp_init.F,v 1.98 2001-11-21 17:27:44 d3j191 Exp $
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical sp_rdrest,sp_diffbb
      external sp_rdrest,sp_diffbb
c
      integer lfntop,lfnrst,lsyn,nsyn,ldbg,idbg
      character*255 filtop,filrst,filsyn
      real*8 rsht,rlng,rinp,boxt(3),vlatt(3,3),rtemp(3)
      real*8 temp,tempw,temps,fld
      integer itemp(2)
      integer npxi,npyi,npzi,nbxi,nbyi,nbzi,isld,irset
      integer nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai
      integer lenscr,ldb,lbp,npbt,nbxt,iguidi,mbbi,nseri
      logical lpol,lfre,ignore,lpbc,lqd
      character*80 prjct
c
      integer i,j,l_f,i_f
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      project=prjct
c
      idebug=idbg
      lfndbg=ldbg
c
      lfnsyn=lsyn
      nfsync=nsyn
      isload=isld
      ireset=irset
c
      rshort=rsht
      rlong=max(rsht,rlng)
      rbox=rinp
      ntwin=1
      if(rlng.gt.rsht) ntwin=2
      ltwin=rlng.gt.rsht
      lnode0=lqmd
      lqmd=lqd
c
      loadb=ldb
      lbpair=lbp
      factld=fld
c
      lpola=lpol
      lfree=lfre
c
      nbx=nbxi
      nby=nbyi
      nbz=nbzi
      nbxin=nbxi
      nbyin=nbyi
      nbzin=nbzi
c
      npx=npxi
      npy=npyi
      npz=npzi
c
      mbbl=0
      nbbdif=-1
c
      mwm=mwmi
      msa=max(msai,msmi)
      mbbreq=mbbi
      mbblp=0
      nserie=nseri
c
      iguide=iguidi
c
      npack=0
      npackw=0
c
      call sp_nrnode()
c
      call sp_dimens(lfnrst,filrst)
c
      call sp_alloc()
c
      call sp_decomp(int_mb(i_iown),dbl_mb(i_boxs))
c
      call sp_initip(int_mb(i_iown),int_mb(i_ipl))
c
      call sp_numbb(int_mb(i_iown),dbl_mb(i_boxs))
c
      if(sp_diffbb(dbl_mb(i_boxs),int_mb(i_rng)))
     + call sp_listbb(int_mb(i_iown),dbl_mb(i_boxs),int_mb(i_bb))
c
      if(ireset.eq.0) then
      ignore=sp_rdrest(lfnrst,filrst,dbl_mb(i_boxs))
      endif

      if(np.gt.1) then
      call ga_brdcst(msp_02,nsm,ma_sizeof(mt_int,1,mt_byte),0)
      endif
      msm=max(1,nsm)
      if(.not.ma_push_get(mt_dbl,msm*3,'xscr',l_xscr,i_xscr))
     + call errquit('Failed to allocate xscr',0)
c
      lenscr=3*max(mwm*mwa,msa)
      if(.not.ma_push_get(mt_dbl,lenscr,'x',l_x,i_x))
     + call errquit('Failed to allocate x',0)
      if(.not.ma_push_get(mt_dbl,lenscr,'v',l_v,i_v))
     + call errquit('Failed to allocate v',0)
      if(.not.ma_push_get(mt_dbl,lenscr,'f',l_f,i_f))
     + call errquit('Failed to allocate f',0)
      if(.not.ma_push_get(mt_dbl,3*mwm*mwa,'r',l_r,i_r))
     + call errquit('Failed to allocate r',0)
      lenscr=max(miw2*mwm,mis2*msa)
      if(.not.ma_push_get(mt_int,lenscr,'i',l_i,i_i))
     + call errquit('Failed to allocate i',0)
      lenscr=ma_inquire_avail(mt_byte)/
     + ((9*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bf',l_bf,i_bf))
     + call errquit('Failed to allocate bf',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call ga_sync()
      call sp_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,int_mb(i_ipl),
     + dbl_mb(i_x),dbl_mb(i_v),dbl_mb(i_f),dbl_mb(i_r),int_mb(i_i),
     + dbl_mb(i_x),dbl_mb(i_v),dbl_mb(i_f),dbl_mb(i_xscr),int_mb(i_i),
     + dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_bf),dbl_mb(i_br),
     + int_mb(i_bi),lenscr,
     + dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_bf),int_mb(i_bi),
     + lenscr,int_mb(i_n),
     + int_mb(i_iown),dbl_mb(i_boxs))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bf))
     + call errquit('Failed to deallocate bf',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
      if(.not.ma_pop_stack(l_i))
     + call errquit('Failed to deallocate i',0)
      if(.not.ma_pop_stack(l_r))
     + call errquit('Failed to deallocate r',0)
      if(.not.ma_pop_stack(l_f))
     + call errquit('Failed to deallocate f',0)
      if(.not.ma_pop_stack(l_v))
     + call errquit('Failed to deallocate v',0)
      if(.not.ma_pop_stack(l_x))
     + call errquit('Failed to deallocate x',0)
c
      if(np.gt.1) then
      rtemp(1)=temp
      rtemp(2)=tempw
      rtemp(3)=temps
      call ga_brdcst(msp_01,rtemp,ma_sizeof(mt_dbl,3,mt_byte),0)
      temp=rtemp(1)
      tempw=rtemp(2)
      temps=rtemp(3)
      itemp(1)=nsm
      itemp(2)=nsf
      call ga_brdcst(msp_02,itemp,ma_sizeof(mt_int,2,mt_byte),0)
      nsm=itemp(1)
      nsf=itemp(2)
      endif
c
      nwmi=nwm
      nwai=nwa
      nsmi=nsm
      nsai=nsa
c
      mwmi=mwm
      mwai=mwa
      msai=msa
c
      npbt=npbtyp
      nbxt=nbxtyp
c
      do 1 i=1,3
      boxt(i)=box(i)
      do 2 j=1,3
      vlatt(i,j)=vlat(i,j)
    2 continue
    1 continue
c
      nsfi=nsf
      msfi=max(1,nsf)
      msmi=max(1,nsm)
c
      third=1.0d0/3.0d0
      nldup=-1
      ipairf=-1
      ipairt=-1
      lpipo=.false.
c
      lpbc=npbtyp.gt.0
c
      if(me.eq.0) then
      if(nfsync.gt.0) then 
      open(unit=lfnsyn,file=filsyn(1:index(filsyn,' ')-1),
     + status='unknown')
      write(lfnsyn,3000) np
 3000 format(i5)
      endif
      endif
c
      call ga_sync()
c
      return
      end
      subroutine sp_dimens(lfnrst,filrst)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
c
      integer lfnrst
      character*255 filrst
c
      character*1 cdum
      integer i,j,jdum,ibx,iby,ibz,itemp(10)
      real*8 rtemp(4)
c
      mbbl=0
c
      if(me.eq.0) then
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfnrst)
c
      do 2 i=1,3
      read(lfnrst,1001) cdum
 1001 format(a1)
    2 continue
      read(lfnrst,1006) nhist
 1006 format(32x,i5)
      do 6 i=1,nhist
      read(lfnrst,1007) hist(i)
 1007 format(a80)
    6 continue
      read(lfnrst,1002) npbtyp,nbxtyp
 1002 format(2i5)
      read(lfnrst,1004) ((vlat(i,j),j=1,3),i=1,3)
 1004 format(3f12.6)
      box(1)=vlat(1,1)
      box(2)=vlat(2,2)
      box(3)=vlat(3,3)
      read(lfnrst,1003) jdum
 1003 format(40x,i5)
      read(lfnrst,1001) cdum
      if(jdum.ne.0) then
      read(lfnrst,1001) cdum
      endif
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,nsf
 1005 format(6i10)
      close(unit=lfnrst,status='keep')
c
      bsize=1.05d0*max(rshort,half*rlong,rbox)
c
      if(nbx*nby*nbz.lt.np) then
      nbx=int(box(1)/bsize)
      nby=int(box(2)/bsize)
      nbz=int(box(3)/bsize)
      endif
c
      nbx=max(1,nbx,npx)
      nby=max(1,nby,npy)
      nbz=max(1,nbz,npz)
c
c      if(nbx.eq.0) nbx=1
c      if(nby.eq.0) nby=1
c      if(nbz.eq.0) nbz=1
c
c      bsize=
c     + min(bsize,box(1)/dble(nbx),box(2)/dble(nby),box(3)/dble(nbz))
c
cx      bxmin=bsize
cx      bymin=bsize
cx      bzmin=bsize
c
c      
c      if(nbx.lt.npx) then
c      i=npx/nbx
c      if(mod(npx,nbx).gt.0) i=i+1
c      bxmin=min(bsize,box(1)/dble(i))
c      nbx=npx
c      endif
c      if(nby.lt.npy) then
c      i=npy/nby
c      if(mod(npy,nby).gt.0) i=i+1
c      bymin=min(bsize,box(2)/dble(i))
c      nby=npy
c      endif
c      if(nbz.lt.npz) then
c      i=npz/nbz
c      if(mod(npz,nbz).gt.0) i=i+1
c      bzmin=min(bsize,box(3)/dble(i))
c      nbz=npz
c      endif
c
      if(nbxin.eq.0) then
      nred(1)=0
      nred(2)=0
      nred(3)=0
      if(nbx.gt.npx.and.mod(nbx,npx).gt.0) then
      nbx=(nbx/npx)*npx
      nred(1)=nbx
      endif
      endif
      if(nbyin.eq.0) then
      if(nby.gt.npy.and.mod(nby,npy).gt.0) then
      nby=(nby/npy)*npy
      nred(2)=nby
      endif
      endif
      if(nbzin.eq.0) then
      if(nbz.gt.npz.and.mod(nbz,npz).gt.0) then
      nbz=(nbz/npz)*npz
      nred(3)=nbz
      endif
      endif
c
      bxmin=bsize/dble(int((dble(nbx)*bsize)/box(1))+1)
      bymin=bsize/dble(int((dble(nby)*bsize)/box(2))+1)
      bzmin=bsize/dble(int((dble(nbz)*bsize)/box(3))+1)
c
      endif
c
      if(np.gt.1) then
      itemp(1)=nwm
      itemp(2)=nwa
      itemp(3)=nsm
      itemp(4)=nsa
      itemp(5)=nbx
      itemp(6)=nby
      itemp(7)=nbz
      itemp(8)=npbtyp
      itemp(9)=nbxtyp
      itemp(10)=nsf
      rtemp(1)=bsize
      rtemp(2)=1.001d0*bxmin
      rtemp(3)=1.001d0*bymin
      rtemp(4)=1.001d0*bzmin
      call ga_brdcst(msp_03,itemp,ma_sizeof(mt_int,10,mt_byte),0)
      call ga_brdcst(msp_04,rtemp,ma_sizeof(mt_dbl,4,mt_byte),0)
      call ga_brdcst(msp_05,box,ma_sizeof(mt_dbl,3,mt_byte),0)
      call ga_brdcst(msp_06,vlat,ma_sizeof(mt_dbl,9,mt_byte),0)
      nwm=itemp(1)
      nwa=itemp(2)
      nsm=itemp(3)
      nsa=itemp(4)
      nbx=itemp(5)
      nby=itemp(6)
      nbz=itemp(7)
      npbtyp=itemp(8)
      nbxtyp=itemp(9)
      nsf=itemp(10)
      bsize=rtemp(1)
      bxmin=rtemp(2)
      bymin=rtemp(3)
      bzmin=rtemp(4)
      endif
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      maxbox=max(nbx,nby,nbz)
      nbtot=nbx*nby*nbz
      lpbc0=nbx.eq.1.or.nby.eq.1.or.nbz.eq.1.or.
     + npx.eq.1.or.npy.eq.1.or.npz.eq.1.or.lpbc9
c
      mbox=30
      do 3 ibx=1,nbx
      do 4 iby=1,nby
      do 5 ibz=1,nbz
      if(me.eq.((((ibz-1)*npz)/nbz)*npy+(((iby-1)*npy)/nby))*npx
     + +((ibx-1)*npx)/nbx) mbox=mbox+1
    5 continue
    4 continue
    3 continue
      mbxloc=mbox-30
c
      if(np.gt.1) call ga_igop(msp_07,mbox,1,'max')
c
      mwa=max(1,nwa)
      msa=max(1,msa,(mbox-30+madbox)*((nwm*nwa+nsa)/nbtot+1)+1)
      mwm=max(1,mwm,msa/mwa+1)
c
      msa=min(msa,2*nsa+1)
      mwm=min(mwm,2*nwm+1)
c
      if(lnode0) then
      msa=nsa+1
      mwm=nwm+1
      endif
c
      return
 9999 call errquit('Failed to open restart file',0)
      return
      end
      subroutine sp_alloc
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer isize
c
      if(.not.ma_push_get(mt_int,3*maxbox,'owner',l_iown,i_iown))
     + call errquit('Failed to allocate owner',0)
      if(.not.ma_push_get(mt_dbl,3*maxbox,'bxsiz',l_boxs,i_boxs))
     + call errquit('Failed to allocate bxsiz',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'ibxrg',l_boxr,i_boxr))
     + call errquit('Failed to allocate ibxrg',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'rng',l_rng,i_rng))
     + call errquit('Failed to allocate rng',me)
c
      if(.not.ma_push_get(mt_int,mip2*mbox,'ipl',l_ipl,i_ipl))
     + call errquit('Failed to allocate ipl',0)
      if(.not.ma_push_get(mt_int,mip2*mbox,'jpl',l_jpl,i_jpl))
     + call errquit('Failed to allocate jpl',0)
c
      if(.not.ga_create(mt_int,np*mbox,mip2,'ip',mbox,mip2,ga_ip))
     + call errquit('Failed to create global array ip',0)
      if(.not.ga_create(mt_int,np*mwm,miw2,'iw',mwm,miw2,ga_iw))
     + call errquit('Failed to create global array iw',0)
      isize=6+12*mwa
      if(lpola) isize=6+18*mwa
      if(lpola.and.lfree) isize=6+30*mwa
      if(.not.ga_create(mt_dbl,np*mwm,isize,'w',mwm,isize,ga_w))
     + call errquit('Failed to create global array w',0)
      if(.not.ga_create(mt_int,np*msa,mis2,'is',msa,mis2,ga_is))
     + call errquit('Failed to create global array is',0)
      isize=39
      if(lpola) isize=45
      if(lpola.and.lfree) isize=57
      if(.not.ga_create(mt_dbl,np*msa,isize,'s',msa,isize,ga_s))
     + call errquit('Failed to create global array s',0)
c
      if(.not.ga_create(mt_int,np*mwm,1,'iwz',mwm,1,ga_iwz))
     + call errquit('Failed to create global array iwz',0)
      if(.not.ga_create(mt_int,np*msa,1,'isz',msa,1,ga_isz))
     + call errquit('Failed to create global array isz',0)
c
      return
      end
      subroutine sp_free()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      if(.not.ga_destroy(ga_isz))
     + call errquit('Error ga_destroy isz',ga_isz)
      if(.not.ga_destroy(ga_iwz))
     + call errquit('Error ga_destroy iwz',ga_iwz)
c
      if(.not.ga_destroy(ga_s))
     + call errquit('Error ga_destroy s',ga_s)
      if(.not.ga_destroy(ga_is))
     + call errquit('Error ga_destroy is',ga_is)
      if(.not.ga_destroy(ga_w))
     + call errquit('Error ga_destroy w',ga_w)
      if(.not.ga_destroy(ga_iw))
     + call errquit('Error ga_destroy iw',ga_iw)
      if(.not.ga_destroy(ga_ip))
     + call errquit('Error ga_destroy iw',ga_ip)
c
      if(.not.ma_pop_stack(l_jpl))
     + call errquit('Failed to deallocate jpl',0)
      if(.not.ma_pop_stack(l_ipl))
     + call errquit('Failed to deallocate ipl',0)
      if(.not.ma_pop_stack(l_rng))
     + call errquit('Failed to deallocate rng',0)
      if(.not.ma_pop_stack(l_boxr))
     + call errquit('Failed to deallocate boxr',0)
      if(.not.ma_pop_stack(l_boxs))
     + call errquit('Failed to deallocate boxs',0)
      if(.not.ma_pop_stack(l_iown))
     + call errquit('Failed to deallocate iown',0)
c
      return
      end
      subroutine sp_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,ipl,xw,vw,fw,xwcr,iwl,xs,vs,fs,xscr,isl,
     + bxw,bvw,bfw,brw,ibw,nw,bxs,bvs,bfs,ibs,ns,ndx,
     + ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
c
      integer sp_btop
      external sp_btop
c
      integer lfnrst,lfntop
      character*255 filrst,filtop
      integer nw,ns
      integer ibownr(maxbox,3)
      integer ipl(mbox,mip2),iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),xscr(msm,3)
      real*8 fw(mwm,3,mwa),fs(msa,3)
      real*8 bxw(nw,3,mwa),bxs(ns,3),brw(nw,3)
      real*8 bvw(nw,3,mwa),bvs(ns,3)
      real*8 bfw(nw,3,mwa),bfs(ns,3)
      integer ibw(nw,2),ibs(ns,mis2)
      real*8 boxsiz(maxbox,3)
      integer ndx(nw)
      real*8 temp,tempw,temps
c
      character*1 cdum
      real*8 rdum,cgx,cgy,cgz
      integer i,j,k,idum,jdum,kdum,number,ncyc,numw,nums
      integer icyc,ibx,iby,ibz,ipx,ipy,ipz,node,new,nold
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ili,ihi,jli,jhi,ilp,ihp,jlp,jhp
c
      integer nat,nqt,naw,nbw,nhw,ndw,now,ntw,nnw,ipp,ips,ix,nsmr
      real*8 boxi(3)
      logical lforces
c
      integer l,m,ib(3),nbox,nb(3),joff
      real*8 xtmin,xtmax,xtx,xt(3)
c
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      real*8 util_nint
      external util_nint
#endif
c
      lpbc9=.true.
      nsmr=0
c
      if(me.eq.0) then
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfnrst)
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + status='old',form='formatted',err=9899)
      rewind(lfntop)
c
      do 101 i=1,5
      read(lfntop,2001,end=9897,err=9898) cdum
 2001 format(a1)
  101 continue
      read(lfntop,2002,end=9897,err=9898) nat,nqt
 2002 format(i5)
      do 102 i=1,nat
      read(lfntop,2001,end=9897,err=9898) cdum
      do 103 j=i,nat
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  103 continue
  102 continue
      do 104 i=1,nqt
      read(lfntop,2001,end=9897,err=9898) cdum
  104 continue
      read(lfntop,2003,end=9897,err=9898) naw,nbw,nhw,ndw,now,ntw,nnw
 2003 format(5i7,2i10)
      read(lfntop,2001,end=9897,err=9898) cdum
      do 105 i=1,naw
      read(lfntop,2001,end=9897,err=9898) cdum
  105 continue
      do 106 i=1,nbw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  106 continue
      do 107 i=1,nhw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  107 continue
      do 108 i=1,ndw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  108 continue
      do 109 i=1,now
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  109 continue
      if(ntw.gt.0) then
      read(lfntop,2004,end=9897,err=9898) (idum,i=1,ntw)
      read(lfntop,2004,end=9897,err=9898) (idum,i=1,ntw)
 2004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,2005,end=9997,err=9998) (idum,i=1,nnw)
      read(lfntop,2005,end=9997,err=9998) (idum,i=1,nnw)
 2005 format(11i7)
      endif
      read(lfntop,2001,end=9897,err=9898) cdum
c
      do 1 i=1,3
      read(lfnrst,1001,end=9997,err=9998) cdum
 1001 format(a1)
    1 continue
      read(lfnrst,1016) nhist,lforces
 1016 format(32x,i5,4x,l1)
      if(nhist.gt.0) then
      do 21 i=1,nhist
      read(lfnrst,1017) hist(i)
 1017 format(a80)
   21 continue
      endif
      if(nhist.lt.mxhist) then
      nhist=nhist+1
      else
      do 22 i=1,nhist-1
      hist(i)=hist(i+1)
   22 continue
      endif
      do 23 i=1,80
      hist(nhist)(i:i)=' '
   23 continue
      read(lfnrst,1002,end=9997,err=9998) npbtyp,nbxtyp
 1002 format(i5,i5)
      read(lfnrst,1003,end=9997,err=9998) ((vlat(i,j),j=1,3),i=1,3)
 1003 format(3f12.6)
      read(lfnrst,1004,end=9997,err=9998) jdum
 1004 format(40x,i5)
      read(lfnrst,1005,end=9997,err=9998) temp,tempw,temps
 1005 format(3f12.6)
      do 2 i=1,3
      box(i)=vlat(i,i)
      boxh(i)=half*box(i)
      boxi(i)=one/box(i)
      do 3 j=1,3
      vlati(i,j)=vlat(i,j)
    3 continue
    2 continue
      call matinv(vlati,3,3)
      if(jdum.ne.0) then
      read(lfnrst,1001,end=9997,err=9998) cdum
      endif
      read(lfnrst,1006,end=9997,err=9998) idum
 1006 format(70x,i5)
      if(idum.gt.0) then
      read(lfnrst,1007,end=9997,err=9998) idum,jdum,kdum
 1007 format(3i5)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,idum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,jdum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,kdum)
 1008 format(4e20.12)
      endif
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nw+1
      numw=nw
      do 4 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
      do 44 i=1,numw
      read(lfnrst,1009,end=9997,err=9998)
     + ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1009 format(2x,6f13.8)
      if(lforces) read(lfnrst,1109,end=9997,err=9998)
     + ((bfw(i,j,k),j=1,3),k=1,nwa)
 1109 format(2x,6e13.6)
      read(lfnrst,1010,end=9997,err=9998) ibw(i,2),(brw(i,k),k=1,3)
 1010 format(i1,1x,3f13.8)
   44 continue
      do 5 i=1,numw
      cgx=zero
      cgy=zero
      cgz=zero
      do 6 k=1,nwa
      cgx=cgx+bxw(i,1,k)
      cgy=cgy+bxw(i,2,k)
      cgz=cgz+bxw(i,3,k)
      if(.not.lforces) then
      bfw(i,1,k)=zero
      bfw(i,2,k)=zero
      bfw(i,3,k)=zero
      endif
    6 continue
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
      do 7 j=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(j,1)) ibx=j
    7 continue
      do 8 j=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(j,2)) iby=j
    8 continue
      do 9 j=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(j,3)) ibz=j
    9 continue
c
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
c
      ndx(i)=(ipz*npy+ipy)*npx+ipx
      ibw(i,1)=(ibz*nby+iby)*nbx+ibx
      if(lnode0) then
      ndx(i)=0
      ibw(i,1)=0
      endif
    5 continue
      do 10 node=0,np-1
      new=0
      do 11 i=1,numw
      if(ndx(i).eq.node) then
      new=new+1
      iwl(new,lwgmn)=number+i
      iwl(new,lwnod)=node
      do 124 j=1,3
      do 12 k=1,nwa
      xw(new,j,k)=bxw(i,j,k)
      vw(new,j,k)=bvw(i,j,k)
      fw(new,j,k)=bfw(i,j,k)
   12 continue
      xwcr(new,j)=brw(i,j)
  124 continue
      iwl(new,lwbox)=ibw(i,1)
      iwl(new,lwdyn)=5*ibw(i,2)
      endif
   11 continue
c
      if(new.gt.0) then
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      nold=ipl(1,2)
      if(nold+new.gt.mwm) call errquit('Dimension mwm too small',0)
      ipl(1,2)=ipl(1,2)+new
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_iw,ili,ihi,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      ilw=ilw+nold
      ihw=ilw+new-1
      call ga_put(ga_w,ilw,ihw,jlw,jlw+3*mwa-1,xw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+3*mwa,jlw+6*mwa-1,vw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+6*mwa,jlw+6*mwa+2,xwcr(1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+6*mwa+3,jlw+9*mwa+2,fw(1,1,1),mwm)
      endif
c
   10 continue
      number=number+numw
    4 continue
      endif
c
      if(nsa.gt.0) then
      nb(1)=nbx
      nb(2)=nby
      nb(3)=nbz
      nbox=0
      k=1
      joff=0
      do 13 i=1,nsa
      if(k.gt.ns) call errquit('Increase memory for sp_rdrst buffer',0)
      read(lfnrst,1011,end=9997,err=9998)
     + ibs(k,11),(bxs(k,j),j=1,3),(bvs(k,j),j=1,3)
 1011 format(i1,1x,6f13.8)
      if(lforces) then
      read(lfnrst,1111,end=9997,err=9998) (bfs(k,j),j=1,3)
 1111 format(2x,6e13.6)
      else
      bfs(k,1)=zero
      bfs(k,2)=zero
      bfs(k,3)=zero
      endif
      read(lfntop,2009,end=9897,err=9898) (ibs(k,j),j=1,10)
 2009 format(16x,10i5)
      if(nsmr.lt.ibs(k,2)) nsmr=ibs(k,2)
      if(nserie.eq.0) then
      if(nsf.lt.ibs(k,1)) nsf=ibs(k,1)
      else
      if(nsf.lt.ibs(k,1))
     + call errquit('Error in number of solute fractions',nsf)
      endif
c
c     if segment of this atom is different from the segment of previous atom 
c     then distribute all previous atoms in the list
c
      if(k.gt.1) then
      if(ibs(k,3).ne.ibs(k-1,3)) then
      new=k-1
      goto 14
      endif
      endif
c
c     if this is the last atom distribute
c
      if(i.eq.nsa) then
      new=k
      goto 14
      endif
c
c     read next atom
c
      k=k+1
      goto 13
c
c     distribute atoms 1 through new
c
   14 continue
c
c     determine the center of geometry
c
      do 15 l=1,3
      xtmax=bxs(1,l)
      xtmin=bxs(1,l)
      do 16 j=1,new
      xtmax=max(xtmax,bxs(j,l))
      xtmin=min(xtmin,bxs(j,l))
   16 continue
      xtx=0.5d0*(xtmax+xtmin)
      if(npbtyp.ne.0) then
      if(abs(xtx).gt.boxh(l)) then
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      xtx=xtx-util_nint(xtx*boxi(l))*box(l)
#else
      xtx=xtx-nint(xtx*boxi(l))*box(l)
#endif
      endif
      endif
      ib(l)=0
      do 17 m=1,nb(l)-1
      if(xtx+boxh(l).gt.boxsiz(m,l)) ib(l)=m
   17 continue
   15 continue
c
c     periodic boundaries
c
      if(npbtyp.gt.0) then
      m=2
      if(npbtyp.eq.1) m=3
      do 18 l=1,m
      if(ib(l).ge.nb(l)) ib(l)=ib(l)-nb(l)
      if(ib(l).lt.0) ib(l)=ib(l)+nb(l)
   18 continue
      if(npbtyp.gt.1) then
      if(ib(3).ge.nb(3)) ib(3)=nb(3)-1
      if(ib(3).lt.0) ib(3)=0
      endif
      else
      do 19 l=1,3
      if(ib(l).ge.nb(l)) ib(l)=nb(l)-1
      if(ib(l).lt.0) ib(l)=0
   19 continue
      endif
c
c     determine owning node
c
      if(.not.lnode0) nbox=(ib(3)*nb(2)+ib(2))*nb(1)+ib(1)
      node=sp_btop(nbox,ibownr)
c
      do 120 j=1,new
      isl(j,lsgan)=joff+j
      isl(j,lsfrc)=ibs(j,1)
      isl(j,lsmol)=ibs(j,2)
      isl(j,lssgm)=ibs(j,3)
      isl(j,lsgrp)=ibs(j,4)
      isl(j,lspgr)=ibs(j,5)
      isl(j,lsatt)=ibs(j,6)
      isl(j,lsct1)=ibs(j,7)
      isl(j,lsct2)=ibs(j,8)
      isl(j,lsct3)=ibs(j,9)
      isl(j,lssss)=ibs(j,10)
      isl(j,lsdyn)=ibs(j,11)
      isl(j,lsbox)=nbox
      isl(j,lsnod)=node
      do 121 l=1,3
      xs(j,l)=bxs(j,l)
      vs(j,l)=bvs(j,l)
      fs(j,l)=bfs(j,l)
  121 continue
  120 continue
      joff=joff+new
c
c     communicate data to node
c
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      nold=ipl(2,2)
      if(nold+new.gt.msa)
     + call errquit('Dimension msa too small',nold+new)
      ipl(2,2)=ipl(2,2)+new
      call ga_put(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_is,ili,ihi,jli,jhi,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      ils=ils+nold
      ihs=ils+new-1
      call ga_put(ga_s,ils,ihs,jls,jls+2,xs(1,1),msa)
      call ga_put(ga_s,ils,ihs,jls+3,jls+5,vs(1,1),msa)
      call ga_put(ga_s,ils,ihs,jls+9,jls+11,fs(1,1),msa)
c
c     make first atom of next segment first in the list
c
      if(k.gt.new) then
      do 122 j=1,11
      ibs(1,j)=ibs(k,j)
  122 continue
      do 123 j=1,3
      bxs(1,j)=bxs(k,j)
      bvs(1,j)=bvs(k,j)
      bfs(1,j)=bfs(k,j)
  123 continue
      k=2
      if(i.eq.nsa) then
      new=1
      k=1
      goto 14
      endif
      endif
c
   13 continue
      endif
c
      if(nsm.gt.0) then
      do 31 i=1,nsm
      xscr(i,1)=zero
      xscr(i,2)=zero
      xscr(i,3)=zero
   31 continue
      do 32 i=1,nsm
      read(lfnrst,1012,end=99,err=99) (xscr(i,j),j=1,3)
 1012 format(2x,3f13.8)
   32 continue
      endif
   99 continue
c
      close(unit=lfnrst,status='keep')
      close(unit=lfntop,status='keep')
c
      endif
c
      return
 9897 call errquit('EOF encountered on topology file',0)
 9898 call errquit('Error reading topology file',0)
 9899 call errquit('Error opening topology file',0)
 9997 call errquit('EOF encountered on restart file',0)
 9998 call errquit('Error reading restart file',0)
 9999 call errquit('Error opening restart file',0)
      return
      end
      integer function sp_btop(ibox,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibox,ibownr(maxbox,3)
      integer iboxx,iboxy,iboxz
c
      iboxx=mod(ibox,nbx)
      iboxy=mod((ibox-iboxx)/nbx,nby)
      iboxz=((ibox-iboxx)/nbx-iboxy)/nby
      sp_btop=(ibownr(iboxz+1,3)*npy+ibownr(iboxy+1,2))*npx+
     + ibownr(iboxx+1,1)
c
      return
      end
      subroutine sp_decomp(ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
c
      integer ibx,iby,ibz,i,ix,jx,iy,jy,iz,jz,jnode,nrnod,nbiown
c
c     check dimensions of ibownr
c
      if(nbx.gt.maxbox.or.nby.gt.maxbox.or.nbz.gt.maxbox)
     + call errquit('Dimension maxbox too small',0)
c
c     determine the node dimension for each sub box
c
      do 1 ibx=1,nbx
      ibownr(ibx,1)=((ibx-1)*npx)/nbx
      boxsiz(ibx,1)=(box(1)*dble(ibx))/dble(nbx)
    1 continue
      do 2 iby=1,nby
      ibownr(iby,2)=((iby-1)*npy)/nby
      boxsiz(iby,2)=(box(2)*dble(iby))/dble(nby)
    2 continue
      do 3 ibz=1,nbz
      ibownr(ibz,3)=((ibz-1)*npz)/nbz
      boxsiz(ibz,3)=(box(3)*dble(ibz))/dble(nbz)
    3 continue
      if(idebug.eq.1) then
      write(lfndbg,8003) nbx,nby,nbz,maxbox,npx,npy,npz
 8003 format('ibownr in sp_decomp',7i5)
      write(lfndbg,8004) (ibownr(ibx,1),ibx=1,nbx)
 8004 format('ibownr x',/,(20i5))
      write(lfndbg,8005) (ibownr(iby,2),iby=1,nby)
 8005 format('ibownr y',/,(20i5))
      write(lfndbg,8006) (ibownr(ibz,3),ibz=1,nbz)
 8006 format('ibownr z',/,(20i5))
      write(lfndbg,8002)
 8002 format('boxlist')
      endif
c
c     determine neighboring nodes and store in neighb(27,2)
c     such that:
c
c     neighb(n,1) is the n-th neighbor of node me
c     neighb(n,2) is the node of which node me is the n-th neighbor
c
c     a value of -1 indicates that such node does not exist
c
      do 10 i=1,27
      neighb(i,1)=-1
      neighb(i,2)=-1
 10   continue
c
      do 11 ix=1,3
      jx=mex+ix-2
      if(npbtyp.gt.0) then
      if(npx.gt.2.and.jx.lt.0) jx=jx+npx
      if(npx.gt.2.and.jx.ge.npx) jx=jx-npx
      endif
      if(jx.ge.0.and.jx.lt.npx) then
      do 12 iy=1,3
      jy=mey+iy-2
      if(npbtyp.gt.0) then
      if(npy.gt.2.and.jy.lt.0) jy=jy+npy
      if(npy.gt.2.and.jy.ge.npy) jy=jy-npy
      endif
      if(jy.ge.0.and.jy.lt.npy) then
      do 13 iz=1,3
      jz=mez+iz-2
      if(npbtyp.eq.1) then
      if(npz.gt.2.and.jz.lt.0) jz=jz+npz
      if(npz.gt.2.and.jz.ge.npz) jz=jz-npz
      endif
      if(jz.ge.0.and.jz.lt.npz) then
      jnode=((jz*npy)+jy)*npx+jx
      neighb(3*(3*(ix-1)+(iy-1))+iz,1)=jnode
      neighb(3*(3*(3-ix)+(3-iy))+4-iz,2)=jnode
      endif
   13 continue
      endif
   12 continue
      endif
   11 continue
c
      nbiown=30
      do 20 ibx=1,nbx
      do 21 iby=1,nby
      do 22 ibz=1,nbz
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
      if(me.eq.nrnod) nbiown=nbiown+1
   22 continue
   21 continue
   20 continue
c
c     broadcast the maximum number of sub-boxes per node to all nodes
c
      if(nbiown.gt.mbox) call errquit('Error in mbox',0)
c
      return
      end
      subroutine sp_nrnode
c
      implicit none
c
#include "sp_common.fh"
c
      integer ix,iy,iz,npt,i,j,k,l
c
c     this routine distributes the available processes in the
c     Cartesian directions
c
c     npx : number of processes in x direction
c     npy : number of processes in y direction
c     npz : number of processes in z direction
c
c     determine node dimensions
c
      if(npx*npy*npz.ne.np) then
      if(npx+npy+npz.gt.0) 
     + call errquit('Specified npx*npy*npz ne np',0)
c
      npt=0
      do 1 i=1,np
      do 2 j=i,np
      do 3 k=j,np
      l=i*j*k
      if(l.eq.np) then
      if(l.eq.npt) then
      if(k.gt.npz) goto 3
      if(i+j+k.lt.npx+npy+npz) then
      npx=i
      npy=j
      npz=k
      endif
      goto 3
      else
      npt=np
      npx=i
      npy=j
      npz=k
      endif
      endif
    3 continue
    2 continue
    1 continue
c
      if(npx*npy*npz.ne.np) call errquit('nrnode: code error',0)
      endif
c
c     determine processor location of me
c
      do 4 ix=1,npx
      do 5 iy=1,npy
      do 6 iz=1,npz
      if(me.eq.((iz-1)*npy+(iy-1))*npx+(ix-1)) then
      mex=ix-1
      mey=iy-1
      mez=iz-1
      endif
    6 continue
    5 continue
    4 continue
c
      return
      end
      subroutine sp_initip(ibownr,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3),ipl(mbox,mip2)
      integer i,j,ibx,iby,ibz,nrbox,nrnod,nbiown
      integer il,ih,jl,jh
c
      do 1 j=1,mip2
      do 2 i=1,30
      ipl(i,j)=0
    2 continue
    1 continue
c
      if(idebug.gt.0) then
      write(lfndbg,8003) nbx,nby,nbz,maxbox
 8003 format('ibownr in sp_initip',4i5)
      write(lfndbg,8004) (ibownr(ibx,1),ibx=1,nbx)
 8004 format('ibownr x',/,(20i5))
      write(lfndbg,8005) (ibownr(iby,2),iby=1,nby)
 8005 format('ibownr y',/,(20i5))
      write(lfndbg,8006) (ibownr(ibz,3),ibz=1,nbz)
 8006 format('ibownr z',/,(20i5))
      call util_flush(lfndbg)
      if(idebug.gt.3) write(lfndbg,8002)
 8002 format('boxlist')
      endif
      nbiown=0
      do 3 ibx=1,nbx
      do 4 iby=1,nby
      do 5 ibz=1,nbz
      nrbox=((ibz-1)*nby+iby-1)*nbx+ibx-1
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
      if(idebug.eq.1) then
      write(lfndbg,8001) ibx,iby,ibz,nrbox,nrnod
 8001 format(5i5)
      call util_flush(lfndbg)
      endif
      if(me.eq.nrnod) then
      nbiown=nbiown+1
      ipl(30+nbiown,1)=nrbox
      ipl(30+nbiown,2)=0
      ipl(30+nbiown,3)=0
      ipl(30+nbiown,4)=0
      ipl(30+nbiown,5)=0
      endif
    5 continue
    4 continue
    3 continue
      ipl(1,1)=nbiown
c
      call ga_distribution(ga_ip,me,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(idebug.eq.3) then
      write(lfndbg,8000) (i,ipl(30+i,1),i=1,ipl(1,1))
 8000 format('ipl',/,(2i5))
      call util_flush(lfndbg)
      endif
c
      return
      end
      subroutine sp_update_i(numsa,isl,numwm,iwl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numsa,isl(msa,mis2)
      integer numwm,iwl(mwm,miw2)
c
      integer il,ih,jl,jh
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,me,il,ih,jl,jh)
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)  
      endif
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,me,il,ih,jl,jh)
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      endif
c
      call ga_sync()
c
      return
      end
      subroutine sp_numbb(ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer jbox,jnode,mbblb,i
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
      character*255 string
c
      if(idebug.gt.0) then
      write(lfndbg,'(a,i6)') 'boxsiz in sp_numbb'
      write(lfndbg,'(6f12.6)') (boxsiz(i,1),i=1,nbx)
      write(lfndbg,'(6f12.6)') (boxsiz(i,2),i=1,nby)
      write(lfndbg,'(6f12.6)') (boxsiz(i,3),i=1,nbz)
      call util_flush(lfndbg)
      endif
c
c     determine size of box-box pairlist
c     ----------------------------------
c
      mbbl=0
      mbblb=0
      do 1 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 2 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 3 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 4 jbx=0,nbx-1
      kbx=jbx-ibx
      jpx=ibownr(jbx+1,1)
c
      if(ibx.le.jbx) then
      ilx=ibx
      ihx=jbx
      else
      ilx=jbx
      ihx=ibx
      endif
c
      dx=zero
      if(ibx.ne.jbx) then
      dx=boxsiz(ihx,1)-boxsiz(ilx+1,1)
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=boxsiz(ilx,1)
      if(ihx.lt.nbx-1) dxtmp=dxtmp-boxsiz(ihx+1,1)+box(1)
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 5 jby=0,nby-1
      kby=jby-iby
      jpy=ibownr(jby+1,2)
c
      if(iby.le.jby) then
      ily=iby
      ihy=jby
      else
      ily=jby
      ihy=iby
      endif
c
      dy=zero
      if(iby.ne.jby) then
      dy=boxsiz(ihy,2)-boxsiz(ily+1,2)
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=boxsiz(ily,2)
      if(ihy.lt.nby-1) dytmp=dytmp-boxsiz(ihy+1,2)+box(2)
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 6 jbz=0,nbz-1
      kbz=jbz-ibz
      jpz=ibownr(jbz+1,3)
c
      if(ibz.le.jbz) then
      ilz=ibz
      ihz=jbz
      else
      ilz=jbz
      ihz=ibz
      endif
c
      dz=zero
      if(ibz.ne.jbz) then
      dz=boxsiz(ihz,3)-boxsiz(ilz+1,3)
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=boxsiz(ilz,3)
      if(ihz.lt.nbz-1) dztmp=dztmp-boxsiz(ihz+1,3)+box(3)
      if(dztmp.lt.dz) dz=dztmp
      if(kbz.gt.0.and.kbz.gt.iabs(kbz-nbz)) kbz=kbz-nbz
      if(kbz.lt.0.and.-kbz.gt.iabs(kbz+nbz)) kbz=kbz+nbz
      endif
      endif
c
      jbox=(jbz*nby+jby)*nbx+jbx
      jnode=(jpz*npy+jpy)*npx+jpx
c
c     determine orientation jbox in relation to ibox
c
c     lside is true if
c
c     i: 0  j: 0  k:  +
c     i: 0  j:  + k:-0+
c     i:  + j:-0+ k:-0+
c
      lside=(kbx.eq.0.and.kby.eq.0.and.kbz.ge.0)
     + .or.(kbx.eq.0.and.kby.gt.0) .or. kbx.gt.0
c
c     determine if ibox is identical to jbox
c
c     lsame=kbx.eq.0.and.kby.eq.0.and.kbz.eq.0
c
c     determine if difference in box numbers is even or odd
c
      leven=2*(iabs(ibox-jbox)/2).eq.iabs(ibox-jbox)
c
c     calculate the distance between the two boxes
c
      if(nbxtyp.eq.1) then
      dist2=
     + (vlat(1,1)*dx/box(1)+vlat(1,2)*dy/box(2)+vlat(1,3)*dz/box(3))**2+
     + (vlat(2,1)*dx/box(1)+vlat(2,2)*dy/box(2)+vlat(2,3)*dz/box(3))**2+
     + (vlat(3,1)*dx/box(1)+vlat(3,2)*dy/box(2)+vlat(3,3)*dz/box(3))**2
      else
      dist2=dx*dx+dy*dy+dz*dz
      endif
c
c     keep half of the box pairs
c
c     this test also appears in sp_numbb
c     any changes need to be made in both routines
c
      if((inode.eq.jnode.and.ibox.ge.jbox).or.(inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rlong*rlong.gt.dist2) mbbl=mbbl+1
      endif
      if(rlong*rlong.gt.dist2) mbblb=mbblb+1
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
c
      mbbl=mbbl+1
      if(mbblb+1.gt.mbbl) mbbl=mbblb+1
c
c     allocate memory for the box-box list
c
      if(np.gt.1) call ga_igop(msp_20,mbbl,1,'max')
      if(mbblp.eq.0) then
      mbbl=max(mbbl,mbbreq)
      if(.not.ma_push_get(mt_int,15*mbbl,'bb',l_bb,i_bb))
     + call errquit('Failed to allocate memory for bb',0)
      mbblp=mbbl
      else
      if(mbbl.gt.mbblp) then
      write(string,1111) mbblp,mbbl,mbbreq
 1111 format('error: lbbl increase from ',i6,' to ',i6,'(',i6,')')
      call errquit(string,me)
c      call errquit('lbbl increased beyond allocated memory',me)
      endif
      endif
c
      return
      end
      subroutine sp_listbb(ibownr,boxsiz,lbbl)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibownr(maxbox,3),lbbl(mbbl,mbb2)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer i,j,jbox,jnode,ltemp
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
c
c     Construction of the box-box pairlist
c
      nbbl=0
      do 1 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 2 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 3 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 4 jbx=0,nbx-1
      kbx=jbx-ibx
      jpx=ibownr(jbx+1,1)
c
      if(ibx.le.jbx) then
      ilx=ibx
      ihx=jbx
      else
      ilx=jbx
      ihx=ibx
      endif
c
      dx=zero
      if(ibx.ne.jbx) then
      dx=boxsiz(ihx,1)-boxsiz(ilx+1,1)
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=boxsiz(ilx,1)
      if(ihx.lt.nbx-1) dxtmp=dxtmp-boxsiz(ihx+1,1)+box(1)
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 5 jby=0,nby-1
      kby=jby-iby
      jpy=ibownr(jby+1,2)
c
      if(iby.le.jby) then
      ily=iby
      ihy=jby
      else
      ily=jby
      ihy=iby
      endif
c
      dy=zero
      if(iby.ne.jby) then
      dy=boxsiz(ihy,2)-boxsiz(ily+1,2)
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=boxsiz(ily,2)
      if(ihy.lt.nby-1) dytmp=dytmp-boxsiz(ihy+1,2)+box(2)
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 6 jbz=0,nbz-1
      kbz=jbz-ibz
      jpz=ibownr(jbz+1,3)
c
      if(ibz.le.jbz) then
      ilz=ibz
      ihz=jbz
      else
      ilz=jbz
      ihz=ibz
      endif
c
      dz=zero
      if(ibz.ne.jbz) then
      dz=boxsiz(ihz,3)-boxsiz(ilz+1,3)
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=boxsiz(ilz,3)
      if(ihz.lt.nbz-1) dztmp=dztmp-boxsiz(ihz+1,3)+box(3)
      if(dztmp.lt.dz) dz=dztmp
      if(kbz.gt.0.and.kbz.gt.iabs(kbz-nbz)) kbz=kbz-nbz
      if(kbz.lt.0.and.-kbz.gt.iabs(kbz+nbz)) kbz=kbz+nbz
      endif
      endif
c
      jbox=(jbz*nby+jby)*nbx+jbx
      jnode=(jpz*npy+jpy)*npx+jpx
c
c     determine orientation jbox in relation to ibox
c
c     lside is true if
c
c     i: 0  j: 0  k:  +
c     i: 0  j:  + k:-0+
c     i:  + j:-0+ k:-0+
c
      lside=(kbx.eq.0.and.kby.eq.0.and.kbz.ge.0)
     + .or.(kbx.eq.0.and.kby.gt.0) .or. kbx.gt.0
c
c     determine if ibox is identical to jbox
c
c     lsame=kbx.eq.0.and.kby.eq.0.and.kbz.eq.0
c
c     determine if difference in box numbers is even or odd
c
      leven=2*(iabs(ibox-jbox)/2).eq.iabs(ibox-jbox)
c
c     calculate the distance between the two boxes
c
      if(nbxtyp.eq.1) then
      dist2=
     + (vlat(1,1)*dx/box(1)+vlat(1,2)*dy/box(2)+vlat(1,3)*dz/box(3))**2+
     + (vlat(2,1)*dx/box(1)+vlat(2,2)*dy/box(2)+vlat(2,3)*dz/box(3))**2+
     + (vlat(3,1)*dx/box(1)+vlat(3,2)*dy/box(2)+vlat(3,3)*dz/box(3))**2
      else
      dist2=dx*dx+dy*dy+dz*dz
      endif
c
c     keep half of the box pairs
c
c     this test also appears in sp_numbb
c     any changes need to be made in both routines
c
      if((inode.eq.jnode.and.ibox.ge.jbox).or. (inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rlong*rlong.gt.dist2) then
      nbbl=nbbl+1
      if(nbbl.gt.mbbl) call errquit('Box-box list too small',mbbl)
      lbbl(nbbl,1)=jnode
      lbbl(nbbl,2)=jbox
      lbbl(nbbl,3)=ibox
      lbbl(nbbl,4)=0
      endif
      endif
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
      npprev=0
c
      nbbloc=0
      do 7 i=1,nbbl-1
      do 8 j=i+1,nbbl
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
    8 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    7 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      if(idebug.gt.2) then
      write(lfndbg,8000) (i,(lbbl(i,j),j=1,3),i=1,nbbl)
 8000 format('lbbl',/,(4i5))
      call util_flush(lfndbg)
      endif
c
      return
      end

      subroutine sp_finish()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      if(.not.ma_pop_stack(l_xscr))
     + call errquit('Failed to deallocate xscr',0)
      if(.not.ma_pop_stack(l_bb))
     + call errquit('Failed to deallocate memory for bb',0)
c
      call sp_free()
c
      return
      end
      subroutine sp_initf(fw,fs,llng,iwz,isz,lpair)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer iwz(mwm),isz(msa)
      logical llng,lpair
c
      integer i,j,k,l,m,il,ih,jl,jh
c
      llong=llng
c
      m=1
      if(llong) m=2
c
      do 1 l=1,m
      if(nwm.gt.0) then
      do 2 k=1,mwa
      do 3 j=1,3
      do 4 i=1,mwm
      fw(i,j,k,l)=zero
    4 continue
    3 continue
    2 continue
      endif
      if(nsa.gt.0) then
      do 5 j=1,3
      do 6 i=1,msa
      fs(i,j,l)=zero
    6 continue
    5 continue
      endif
    1 continue
c
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(llong) call ga_put(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(llong) call ga_put(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa)
      endif
c
      if(lpair) then
      do 7 i=1,mwm
      iwz(i)=0
    7 continue
      do 8 i=1,msa
      isz(i)=0
    8 continue
      call ga_zero(ga_iwz)
      call ga_zero(ga_isz)
      endif
c
      return
      end
      subroutine sp_copyg(fw,fs)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa),fs(msa,3)
c
      integer il,ih,jl,jh
c
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      endif
c
      return
      end
      subroutine sp_final(fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,il,ih,jl,jh
c
      if(np.gt.0) then
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_acc(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm,one)
      if(llong) call ga_acc(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm,one)
      call ga_get(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(ltwin) call ga_get(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_acc(ga_s,il,ih,jl+6,jl+8,fs,msa,one)
      if(llong) call ga_acc(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa,one)
      call ga_get(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(ltwin) call ga_get(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa)
      endif
      endif
c
      if(lpair) then
      if(nwm.gt.0) then
      call ga_distribution(ga_iwz,me,il,ih,jl,jh)
      call ga_acc(ga_iwz,il,ih,1,1,iwz,mwm,1)
      call ga_get(ga_iwz,il,ih,1,1,iwz,mwm)
      do 1 i=1,nwmloc
      iwz(i)=min(1,iwz(i))
    1 continue
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_isz,me,il,ih,jl,jh)
      call ga_acc(ga_isz,il,ih,1,1,isz,msa,1)
      call ga_get(ga_isz,il,ih,1,1,isz,msa)
      do 2 i=1,nsaloc
      isz(i)=min(1,isz(i))
    2 continue
      endif
      endif
c
      return
      end
      subroutine sp_wrtrst(lfnrst,filrst,lveloc,
     + pres,temp,tempw,temps,iwl,xw,vw,fw,xwcr,isl,xs,vs,fs,xscr,prjct)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst
      character*255 filrst
      logical lveloc
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3),xscr(msm,3)
      character*80 prjct
c
      integer lenscr
c
      project=prjct
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((9*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bf',l_bf,i_bf))
     + call errquit('Failed to allocate bf',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call sp_wtrst(lfnrst,filrst,lveloc,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,fw,xwcr,isl,int_mb(i_pack),xs,vs,fs,
     + xscr,int_mb(i_ipl),lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + dbl_mb(i_bf),dbl_mb(i_br),int_mb(i_bi),dbl_mb(i_bx),
     + dbl_mb(i_bv),dbl_mb(i_bf))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bf))
     + call errquit('Failed to deallocate bf',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
c
      return
      end
      subroutine sp_wtrst(lfnrst,filrst,lveloc,pres,temp,tempw,temps,
     + iwl,iwlp,xw,vw,fw,xwcr,isl,islp,xs,vs,fs,xscr,
     + ipl,nb,ibw,bxw,bvw,bfw,brw,ibs,bxs,bvs,bfs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst,nb
      character*255 filrst
      logical lveloc
      real*8 pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3),xscr(msm,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),bfw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3),bfs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
      integer idyn,idynp
      logical lforces
c
      lforces=iguide.ne.0
c
      if(ga_nodeid().eq.0) then
c
      call swatch(rdate,rtime)
c
      rewind(lfnrst)
      write(lfnrst,1000)
 1000 format('Restart file',/,' ',/,' ')
      write(lfnrst,1001) 3.3,rdate,rtime,nhist,lforces
 1001 format(f12.6,2a10,i5,4x,l1)
      hist(nhist)(1:10)=rdate
      hist(nhist)(11:20)=rtime
      hist(nhist)(21:80)=project(1:60)
      do 10 i=1,nhist
      write(lfnrst,1009) hist(i)
 1009 format(a80)
   10 continue
      write(lfnrst,1002) npbtyp,nbxtyp,((vlat(i,j),j=1,3),i=1,3)
 1002 format(2i5,/,(3f12.6))
      write(lfnrst,1003) pres
 1003 format(1pe12.5)
      write(lfnrst,1004) temp,tempw,temps
 1004 format(3f12.6)
      write(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,nsf,0,0,0
 1005 format(7i10,2i5)
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1112 i=1,nb
      ibw(i)=-1
 1112 continue
c
c     end test code
c
      do 2 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmn=ipl(1,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+lwdyn-1,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call sp_unpackw(nwmn,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lveloc)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      if(lforces)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa+3,jlw+9*mwa+2,fw,mwm)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa,jlw+6*mwa+2,xwcr,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      bxw(j,1,k)=xw(i,1,k)
      bxw(j,2,k)=xw(i,2,k)
      bxw(j,3,k)=xw(i,3,k)
      bvw(j,1,k)=vw(i,1,k)
      bvw(j,2,k)=vw(i,2,k)
      bvw(j,3,k)=vw(i,3,k)
      if(lforces) then
      bfw(j,1,k)=fw(i,1,k)
      bfw(j,2,k)=fw(i,2,k)
      bfw(j,3,k)=fw(i,3,k)
      endif
    4 continue
      brw(j,1)=xwcr(i,1)
      brw(j,2)=xwcr(i,2)
      brw(j,3)=xwcr(i,3)
      ibw(j)=iwl(i,lwdyn)
      endif
    3 continue
      endif
    2 continue
      do 5 i=1,numw
      if(lveloc) then
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
      else
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(zero,j=1,3),k=1,nwa)
      endif
 1006 format(2x,6f13.8)
      if(lforces) write(lfnrst,1106) ((bfw(i,j,k),j=1,3),k=1,nwa)
 1106 format(2x,6e13.6)
      idyn=iand(ibw(i),12)/4
      idynp=iand(ibw(i),3)
      write(lfnrst,1007) idynp,idyn,(brw(i,k),k=1,3)
 1007 format(2i1,3f13.8)
c
c     begin test code 10/31/2001
c     check if al atoms have been received
c
      if(ibw(i).lt.0)
     + call errquit('Missing solvent in wtrst',i)
c
c     end test code
c
    5 continue
      number=number+numw
    1 continue
      endif
c
      if(nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 6 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1117 i=1,nb
      ibs(i)=-1
 1117 continue
c
c     end test code
c
      do 7 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lveloc) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      if(lforces) call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 8 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      bvs(j,1)=vs(i,1)
      bvs(j,2)=vs(i,2)
      bvs(j,3)=vs(i,3)
      if(lforces) then
      bfs(j,1)=fs(i,1)
      bfs(j,2)=fs(i,2)
      bfs(j,3)=fs(i,3)
      endif
      ibs(j)=isl(i,lsdyn)
      endif
    8 continue
      endif
    7 continue
      do 9 i=1,nums
      idyn=iand(ibs(i),12)/4
      idynp=iand(ibs(i),3)
      if(lveloc) then
      write(lfnrst,1008) idynp,idyn,(bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
      else
      write(lfnrst,1008) idynp,idyn,(bxs(i,j),j=1,3),(zero,j=1,3)
      endif
 1008 format(2i1,6f13.8)
      if(lforces) write(lfnrst,1108) (bfs(i,j),j=1,3)
 1108 format(2x,3e13.6)
c
c     begin test code 10/31/2001
c     check if al atoms have been received
c
      if(ibs(i).lt.0)
     + call errquit('Missing solute atom in wtrst',i)
c
c     end test code
c
    9 continue
      number=number+nums
    6 continue
      endif
c
      if(nsm.gt.0) then
      do 21 i=1,nsm
      write(lfnrst,1109) (xscr(i,j),j=1,3)
 1109 format(2x,3f13.8)
   21 continue
      endif
c
      endif
c
      return
 9999 continue
      call errquit('Failed to open restart file',me)
      return
      end
      subroutine sp_print(lfnout)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
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
     + call errquit('Failed to allocate memory for lcnt',0)
c
      call sp_prtcnt(lfnout,int_mb(i_lcnt))
c
      if(.not.ma_pop_stack(l_lcnt))
     + call errquit('Failed to deallocate lcnt',0)
c
      return
      end
      subroutine sp_prtcnt(lfnout,lcnt)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnout,lcnt(np,3)
c
      integer i,j
c
      do 1 i=1,np
      lcnt(i,1)=0
      lcnt(i,2)=0
      lcnt(i,3)=0
    1 continue
c
      lcnt(me+1,1)=mbxloc
      lcnt(me+1,2)=nwmloc*nwa
      lcnt(me+1,3)=nsaloc
c
      if(np.gt.1) call ga_igop(msp_08,lcnt,3*np,'+')
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' Initial distribution p:b(w+s)',/)
      write(lfnout,1001) (i-1,(lcnt(i,j),j=1,3),i=1,np)
 1001 format(4(3x,i4,':',i5,'(',i7,'+',i7,')'))
      endif
c
      return
      end
      subroutine sp_printf(lfnout,filtop,lfntop,isl,xs,fs,
     + npener,esa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnout,lfntop
      character*70 filtop
      integer isl(msa,mis2),npener
      real*8 xs(msa,3),fs(msa,3,2),esa(nsa,2)
c
      integer lenscr
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((6*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call sp_prt_s(lfnout,filtop,lfntop,
     + int_mb(i_ipl),isl,int_mb(i_pack),xs,fs,
     + lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + npener,esa)
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
c
      return
      end
      subroutine sp_prt_s(lfnout,filtop,lfntop,
     + ipl,isl,islp,xs,fs,nb,ibs,bxs,bfs,
     + npener,esa)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnout,lfntop,nb,npener
      character*70 filtop
      integer ipl(mbox,mip2),isl(msa,mis2),islp(msa,npack),ibs(nb)
      real*8 xs(msa,3),fs(msa,3)
      real*8 bxs(nb,3),bfs(nb,3)
      real*8 esa(nsa,2)
c
      integer i,j,number,icyc,ncyc,node,nums,nsan,idyn,ism,iss
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ils,ihs,jls,jhs
      integer naw,nbw,nhw,ndw,now,ntw,nnw,nat,nqt,idum
      character*1 cdum
      character*16 cat
c
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
c
      if(nsan.gt.0) then
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_put(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      endif
c
      call ga_sync()
c
      if(me.ne.0) return
c
      if(npener.le.0) then
      write(lfnout,1007)
 1007 format(//,' Solute coordinates and forces',//,
     + '   mol  segment     atom  dt     x       y       z ',
     + '             fx          fy          fz',/)
      else
      write(lfnout,1008)
 1008 format(//,' Solute coordinates, forces and energies',//,
     + '   mol  segment     atom  dt     x       y       z ',
     + '             fx          fy          fz        enb',
     + '          eb          ep',/)
      endif
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + status='old',form='formatted')
      rewind(lfntop)
      do 101 i=1,5
      read(lfntop,2001) cdum
 2001 format(a1)
  101 continue
      read(lfntop,2002) nat,nqt
 2002 format(i5)
      do 102 i=1,nat
      read(lfntop,2001) cdum
      do 103 j=i,nat
      read(lfntop,2001) cdum
      read(lfntop,2001) cdum
  103 continue
  102 continue
      do 104 i=1,nqt
      read(lfntop,2001) cdum
  104 continue
      read(lfntop,2003) naw,nbw,nhw,ndw,now,ntw,nnw
 2003 format(5i7,2i10)
      read(lfntop,2001) cdum
      do 105 i=1,naw
      read(lfntop,2001) cdum
  105 continue
      do 106 i=1,nbw
      read(lfntop,2001) cdum
      read(lfntop,2001) cdum
  106 continue
      do 107 i=1,nhw
      read(lfntop,2001) cdum
      read(lfntop,2001) cdum
  107 continue
      do 108 i=1,ndw
      read(lfntop,2001) cdum
      read(lfntop,2001) cdum
  108 continue
      do 109 i=1,now
      read(lfntop,2001) cdum
      read(lfntop,2001) cdum
  109 continue
      if(ntw.gt.0) then
      read(lfntop,2004) (idum,i=1,ntw)
      read(lfntop,2004) (idum,i=1,ntw)
 2004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,2005) (idum,i=1,nnw)
      read(lfntop,2005) (idum,i=1,nnw)
 2005 format(11i7)
      endif
      read(lfntop,2001) cdum
c
      if(nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 6 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
      do 7 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 8 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      bfs(j,1)=fs(i,1)
      bfs(j,2)=fs(i,2)
      bfs(j,3)=fs(i,3)
      ibs(j)=isl(i,lsdyn)
      endif
    8 continue
      endif
    7 continue
      do 9 i=1,nums
      read(lfntop,2009) cat,ism,iss
 2009 format(a16,5x,2i5)
      idyn=iand(ibs(i),3)
      if(npener.le.0) then
      write(lfnout,1009) ism,iss,cat,idyn,
     + (bxs(i,j),j=1,3),(bfs(i,j),j=1,3)
 1009 format(2i5,':',a16,i1,1x,3f8.4,3x,3f12.3)
      else
      write(lfnout,1010) ism,iss,cat,idyn,
     + (bxs(i,j),j=1,3),(bfs(i,j),j=1,3),esa(i,1),esa(i,2),
     + esa(i,1)+esa(i,2)
 1010 format(2i5,':',a16,i1,1x,3f8.4,3x,3f12.3,3(2x,1pe10.3))
      endif
    9 continue
    6 continue
      endif
c
      close(unit=lfntop,status='keep')
c
      return
      end
      subroutine sp_wtrest(lfn)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfn
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call errquit('Failed to allocate memory for ltmp',0)
      call sp_wrest(lfn,int_mb(i_bb),int_mb(i_ltmp),mbbl+1,
     + dbl_mb(i_boxs))
      if(.not.ma_pop_stack(l_ltmp))
     + call errquit('Failed to deallocate ltmp',0)
c
      return
      end
      subroutine sp_wrest(lfn,lbbl,ltemp,mdim,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfn,mdim
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
      real*8 boxsiz(maxbox,3)
c
      integer i,j,k
c
      if(me.eq.0) then
      write(lfn,1000)
 1000 format('restart space')
      write(lfn,1001) np,mbbl,npx,npy,npz,nbx,nby,nbz,np
 1001 format(9i7)
      write(lfn,1002) (boxsiz(i,1),i=1,nbx)
      write(lfn,1002) (boxsiz(i,2),i=1,nby)
      write(lfn,1002) (boxsiz(i,3),i=1,nbz)
 1002 format(4e20.12)
      endif
c
      do 1 i=1,np
      if(i.eq.me+1) then
      do 2 k=1,mbb2
      ltemp(1,k)=0
      do 3 j=1,nbbl
      ltemp(j+1,k)=lbbl(j,k)
    3 continue
    2 continue
      ltemp(1,1)=nbbl
      else
      do 4 k=1,mbb2
      do 5 j=1,mdim
      ltemp(j,k)=0
    5 continue
    4 continue
      endif
c
      call ga_igop(msp_09,ltemp,mdim*mbb2,'+')
c
      if(me.eq.0) then
      write(lfn,1003) i-1,ltemp(1,1)
 1003 format(2i7)
      do 6 j=1,ltemp(1,1)
      write(lfn,1004) (ltemp(j+1,k),k=1,mbb2)
 1004 format(8i10)
    6 continue
      endif
c
    1 continue
c
      return
      end
      logical function sp_rdrest(lfn,fil,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfn
      character*255 fil
      real*8 boxsiz(maxbox,3)
c
      logical sp_rrest
      external sp_rrest
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call errquit('Failed to allocate memory for ltmp',0)
      sp_rdrest=sp_rrest(lfn,fil,int_mb(i_bb),int_mb(i_ltmp),
     + mbbl+1,boxsiz)
      if(.not.ma_pop_stack(l_ltmp))
     + call errquit('Failed to deallocate ltmp',0)
c
      return
      end
      logical function sp_rrest(lfn,fil,lbbl,ltemp,mdim,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfn,mdim
      character*255 fil
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
c
      integer i,j,k,node,npp,nbl,nbytes,nbxs
      integer npxp,npyp,npzp,nbxp,nbyp,nbzp,mbblx
      character*13 string
      real*8 boxsiz(maxbox,3)
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      npp=0
      read(lfn,1000,end=9997) string
 1000 format(a13)
      if(string.ne.'restart space') goto 1
      read(lfn,1001) npp,mbblx,npxp,npyp,npzp,nbxp,nbyp,nbzp,nbxs
 1001 format(9i7)
      if(mbblp.lt.mbblx) npp=0
      if(npxp.ne.npx) npp=0
      if(npyp.ne.npy) npp=0
      if(npzp.ne.npz) npp=0
      if(nbxp.ne.nbx) npp=0
      if(nbyp.ne.nby) npp=0
      if(nbzp.ne.nbz) npp=0
      if(nbxs.gt.0.and.npp.gt.0) then
      read(lfn,1002) (boxsiz(i,1),i=1,nbxp)
      read(lfn,1002) (boxsiz(i,2),i=1,nbyp)
      read(lfn,1002) (boxsiz(i,3),i=1,nbzp)
 1002 format(4e20.12)
      endif
 9997 continue
      do 2 i=1,mbb2
      ltemp(1,i)=0
    2 continue
      endif
c
      nbytes=ma_sizeof(mt_int,1,mt_byte)
      call ga_brdcst(msp_10,npp,nbytes,0)
c
      if(np.ne.npp) goto 9998
c
      nbytes=ma_sizeof(mt_dbl,3*maxbox,mt_byte)
      call ga_brdcst(msp_10,boxsiz,nbytes,0)
c
      nbytes=mdim*mbb2*ma_sizeof(mt_int,1,mt_byte)
c
      do 3 i=1,np
c
      if(me.eq.0) then
      read(lfn,1003) node,nbl
 1003 format(2i7)
      ltemp(1,1)=node
      ltemp(1,2)=nbl
      do 4 j=1,nbl
      read(lfn,1004) (ltemp(j+1,k),k=1,mbb2)
 1004 format(8i10)
    4 continue
      endif
c
      call ga_brdcst(msp_11,ltemp,nbytes,0)
c
      if(ltemp(1,1).eq.me) then
      nbbl=ltemp(1,2)
      do 5 k=1,mbb2
      do 6 j=1,nbbl
      lbbl(j,k)=ltemp(j+1,k)
    6 continue
    5 continue
      endif
c
    3 continue
c
      nable=1
      sp_rrest=.true.
      return
 9998 continue
      nable=2
      sp_rrest=.false.
      return
 9999 continue
      call errquit('Failed to open restart file',0)
      sp_rrest=.false.
      return
      end
      subroutine sp_wrttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,iwl,xw,vw,fw,xwcr,isl,xs,vs,fs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfntrj
      logical lxw,lvw,lfw,lxs,lvs,lfs
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3)
c
      integer lenscr
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((9*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bf',l_bf,i_bf))
     + call errquit('Failed to allocate bf',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call sp_wttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,fw,xwcr,isl,int_mb(i_pack),xs,vs,fs,
     + int_mb(i_ipl),lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + dbl_mb(i_bf),dbl_mb(i_br),int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + dbl_mb(i_bf))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bf))
     + call errquit('Failed to deallocate bf',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
c
      return
      end
      subroutine sp_wttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,iwl,iwlp,xw,vw,fw,xwcr,
     + isl,islp,xs,vs,fs,ipl,nb,ibw,bxw,bvw,bfw,brw,ibs,bxs,bvs,bfs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfntrj,nb
      logical lxw,lvw,lfw,lxs,lvs,lfs
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),bfw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3),bfs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
c
      if(me.eq.0) then
c
      call swatch(rdate,rtime)
c
      write(lfntrj,1000)
 1000 format('frame')
      write(lfntrj,1001) stime,temp,pres,rdate,rtime
 1001 format(2f12.6,1pe12.5,1x,2a10)
      write(lfntrj,1002) ((vlat(i,j),j=1,3),i=1,3)
 1002 format(3f12.6)
      write(lfntrj,1003) lxw,lvw,lfw,lxs,lvs,lfs,nwm,nwa,nsa
 1003 format(6l1,3i10)
c
      if((lxw.or.lvw.or.lfw).and.nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1112 i=1,nb
      ibw(i)=-1
 1112 continue
c
c     end test code
c
      do 2 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmn=ipl(1,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+lwdyn-1,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call sp_unpackw(nwmn,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lvw)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      if(lfw)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa+3,jlw+9*mwa+2,fw,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      bxw(j,1,k)=xw(i,1,k)
      bxw(j,2,k)=xw(i,2,k)
      bxw(j,3,k)=xw(i,3,k)
    4 continue
      if(lvw) then
      do 5 k=1,nwa
      bvw(j,1,k)=vw(i,1,k)
      bvw(j,2,k)=vw(i,2,k)
      bvw(j,3,k)=vw(i,3,k)
    5 continue
      endif
      if(lfw) then
      do 51 k=1,nwa
      bfw(j,1,k)=fw(i,1,k)
      bfw(j,2,k)=fw(i,2,k)
      bfw(j,3,k)=fw(i,3,k)
   51 continue
      endif
      ibw(j)=iwl(i,lwdyn)
      endif
    3 continue
      endif
    2 continue
      if(lfw) then
      if(lvw) then
      do 61 i=1,numw
      write(lfntrj,1014) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),
     + (bfw(i,j,k),j=1,3),k=1,nwa)
 1014 format(6f8.3,3f8.1)
      if(ibw(i).lt.0) call errquit('Missing solvent in wttrj',i)
   61 continue
      else
      do 71 i=1,numw
      write(lfntrj,1015) ((bxw(i,j,k),j=1,3),(bfw(i,j,k),j=1,3),k=1,nwa)
 1015 format(3f8.3,3f8.1)
      if(ibw(i).lt.0) call errquit('Missing solvent in wtrst',i)
   71 continue
      endif
      else
      if(lvw) then
      do 6 i=1,numw
      write(lfntrj,1004) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1004 format(6f8.3)
      if(ibw(i).lt.0) call errquit('Missing solvent in wttrj',i)
    6 continue
      else
      do 7 i=1,numw
      write(lfntrj,1005) ((bxw(i,j,k),j=1,3),k=1,nwa)
 1005 format(3f8.3)
      if(ibw(i).lt.0) call errquit('Missing solvent in wtrst',i)
    7 continue
      endif
      endif
      number=number+numw
    1 continue
      endif
c
      if((lxs.or.lvs.or.lfs).and.nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 8 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1117 i=1,nb
      ibs(i)=-1
 1117 continue
c
c     end test code
c
      do 9 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lvs) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      if(lfs) call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 10 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      if(lvs) then
      bvs(j,1)=vs(i,1)
      bvs(j,2)=vs(i,2)
      bvs(j,3)=vs(i,3)
      endif
      if(lfs) then
      bfs(j,1)=fs(i,1)
      bfs(j,2)=fs(i,2)
      bfs(j,3)=fs(i,3)
      endif
      ibs(j)=isl(i,lsdyn)
      endif
   10 continue
      endif
    9 continue
      if(lfs) then
      if(lvs) then
      do 111 i=1,nums
      write(lfntrj,1016) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3),
     + (bfs(i,j),j=1,3)
 1016 format(6f8.3,3f8.1)
      if(ibs(i).lt.0) call errquit('Missing solute atom in wtrst',i)
  111 continue
      else
      do 121 i=1,nums
      write(lfntrj,1017) (bxs(i,j),j=1,3),(bfs(i,j),j=1,3)
 1017 format(3f8.3,3f8.1)
      if(ibs(i).lt.0) call errquit('Missing solute atom in wtrst',i)
  121 continue
      endif
      else
      if(lvs) then
      do 11 i=1,nums
      write(lfntrj,1006) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
 1006 format(6f8.3)
      if(ibs(i).lt.0) call errquit('Missing solute atom in wtrst',i)
   11 continue
      else
      do 12 i=1,nums
      write(lfntrj,1007) (bxs(i,j),j=1,3)
 1007 format(3f8.3)
      if(ibs(i).lt.0) call errquit('Missing solute atom in wtrst',i)
   12 continue
      endif
      endif
      number=number+nums
    8 continue
      endif
c
      endif
c
      return
      end
      subroutine sp_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      call sp_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,xwcr,isl,int_mb(i_pack),xs,vs,
     + int_mb(i_ipl))
c
      return
      end
      subroutine sp_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,iwlp,xw,vw,xwcr,isl,islp,xs,vs,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
c
      write(lfnmro) nwm,nwa,nsa,stime,temp,pres,vlat
c
      do 1 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      write(lfnmro) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call sp_unpackw(nwmn,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      write(lfnmro) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      write(lfnmro) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      write(lfnmro) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      endif
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      write(lfnmro) ((isl(j,k),j=1,nsan),k=1,mis2)
      write(lfnmro) ((xs(j,k),j=1,nsan),k=1,3)
      write(lfnmro) ((vs(j,k),j=1,nsan),k=1,3)
      endif
    1 continue
c
      return
      end
      subroutine sp_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmri
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      call sp_rmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,xwcr,isl,int_mb(i_pack),xs,vs,
     + int_mb(i_ipl))
c
      return
      end
      subroutine sp_rmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,iwlp,xw,vw,xwcr,isl,islp,xs,vs,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnmri
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer i,j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ltemp(3)
      real*8 rtemp(12)
c
      if(me.eq.0) then
      read(lfnmri) ltemp,rtemp
      do 1 node=np-1,0,-1
      read(lfnmri) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_put(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      if(nwmn.gt.0) then
      read(lfnmri) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      read(lfnmri) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      read(lfnmri) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_put(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      else
      call sp_packw(nwmn,iwl,iwlp)
      call ga_put(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsan.gt.0) then
      read(lfnmri) ((isl(j,k),j=1,nsan),k=1,mis2)
      read(lfnmri) ((xs(j,k),j=1,nsan),k=1,3)
      read(lfnmri) ((vs(j,k),j=1,nsan),k=1,3)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_put(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      else
      call sp_pack(nsan,isl,islp)
      call ga_put(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_put(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_put(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      endif
    1 continue
      endif
c
      call ga_brdcst(msp_12,ltemp,3*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(msp_13,rtemp,12*ma_sizeof(mt_dbl,1,mt_byte),0)
      nwm=ltemp(1)
      nwa=ltemp(2)
      nsa=ltemp(3)
      stime=rtemp(1)
      temp=rtemp(2)
      pres=rtemp(3)
      k=3
      do 2 i=1,3
      do 3 j=1,3
      k=k+1
      vlat(i,j)=rtemp(k)
    3 continue
    2 continue
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmloc=ipl(1,2)
      nsaloc=ipl(2,2)
      if(nwmloc.gt.0) then
      call ga_distribution(ga_iw,me,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmloc-1,jli,jhi,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmloc-1,jli,jli+npackw-1,iwlp,mwm)
      call sp_unpackw(nwmloc,iwl,iwlp)
      endif
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsaloc.gt.0) then
      call ga_distribution(ga_is,me,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsaloc-1,jli,jhi,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsaloc-1,jli,jli+npack-1,islp,msa)
      call sp_unpack(nsaloc,isl,islp)
      endif
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls+3,jls+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_fix()
c
      implicit none
c
      return
      end
