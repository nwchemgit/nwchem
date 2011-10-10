      subroutine argos_space_dldbal(stime,syntim,waltim,nod,nlst,
     + nen,lbbl,dlb,lpsyn,ibindx,iburen,npp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer argos_space_btop
      external argos_space_btop
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
c      nlst(i,1)=lshift(nlst(i,1),32)+nlst(i,2)
c    3 continue
c      call ga_igop(msp_14,nlst,np,'+')
c      do 4 i=1,np
c      nlst(i,2)=iand(nlst(i,1),8589934591)
c      nlst(i,1)=rshift(nlst(i,1),32)
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
      call argos_space_resize(ibusy-1,int_mb(i_iown),
     + dbl_mb(i_boxs),factor,facx,facy,facz)
      if(me.eq.0.and.lpsyn) then
      write(lfnsyn,3000) 1,stime,synsum,waltim,ibusy-1,
     + factor,facx,facy,facz,lpipo
 3000 format('synchronization',/,i5,3f12.6,/,i5,4f12.6,4x,l1)
      do 9 i=1,np
      write(lfnsyn,3001) nod(i),dlb(nod(i)+1)
 3001 format(i5,f12.6)
    9 continue
      call argos_space_wrtbxsz(lfnsyn,dbl_mb(i_boxs))
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
      n=rshift(n,1)
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
      call argos_space_cascad(nen,int_mb(i_bb))
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
c $Id$
