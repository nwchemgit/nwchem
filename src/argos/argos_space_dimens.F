      subroutine argos_space_dimens(lfnrst,filrst)
c
      implicit none
c
#include "argos_space_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "util.fh"
c
      integer lfnrst
      character*255 filrst
c
      character*1 cdum
      integer i,j,jdum,ibx,iby,ibz,itemp(11)
      real*8 rtemp(5),rsgmr
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
 1007 format(a)
    6 continue
      read(lfnrst,1002) npbtyp,nbxtyp,rsgmr
 1002 format(2i5,f12.6)
      if(rsgm.lt.0.0d0) rsgm=rsgmr
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
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,nsf,nseq
 1005 format(7i10)
      close(unit=lfnrst,status='keep')
c
      bsize=max(rshort+half*rsgm,half*(rlong+half*rsgm),rbox)
c
      if(util_print('distribution',print_default)) then
      if(me.eq.0) write(lfnout,2001)
 2001 format(/,' Distribution information',/)
      if(me.eq.0) write(lfnout,2002) rshort,rsgm,rlong,rbox,box,bsize
 2002 format(' Short range cutoff',t35,f12.6,/,
     + ' Segment size',t35,f12.6,/,' Long range cutoff',t35,f12.6,/,
     + ' Box size rbox',t35,f12.6,//,' Box dimension',t35,3f12.6,//,
     + ' Initial cell size',t35,f12.6,/)
      endif
c
cx      if(nbx*nby*nbz.lt.np) then
      nbx=int(box(1)/bsize)
      nby=int(box(2)/bsize)
      nbz=int(box(3)/bsize)
cx      endif
c
cx      nbx=max(1,nbx,npx)
cx      nby=max(1,nby,npy)
cx      nbz=max(1,nbz,npz)
c
      if(util_print('distribution',print_default)) then
      if(me.eq.0) then
      write(lfnout,2003) nbx,nby,nbz
 2003 format(' Cell distribution',t35,3i5)
      endif
      endif
c
cx      if(nbxin.eq.0) then
cx      nred(1)=0
cx      nred(2)=0
cx      nred(3)=0
cx      if(nbx.gt.npx.and.mod(nbx,npx).gt.0) then
cx      nbx=(nbx/npx)*npx
cx      nred(1)=nbx
cx      endif
cx      endif
cx      if(nbyin.eq.0) then
cx      if(nby.gt.npy.and.mod(nby,npy).gt.0) then
cx      nby=(nby/npy)*npy
cx      nred(2)=nby
cx      endif
cx      endif
cx      if(nbzin.eq.0) then
cx      if(nbz.gt.npz.and.mod(nbz,npz).gt.0) then
cx      nbz=(nbz/npz)*npz
cx      nred(3)=nbz
cx      endif
cx      endif
c
cx      if(util_print('distribution',print_default)) then
cx      if(me.eq.0) then
cx      write(lfnout,2004) nbx,nby,nbz
cx 2004 format(' Final cell distribution',t35,3i5,/)
cx      endif
cx      endif
c
      bxmin=bsize/dble(int((dble(nbx)*bsize)/box(1))+1)
      bymin=bsize/dble(int((dble(nby)*bsize)/box(2))+1)
      bzmin=bsize/dble(int((dble(nbz)*bsize)/box(3))+1)
c
      if(util_print('distribution',print_default)) then
      if(me.eq.0) then
      write(lfnout,2005) bxmin,bymin,bzmin
 2005 format(' Minimum cell size',t35,3f12.6)
      endif
      endif
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
      itemp(11)=nseq
      rtemp(1)=bsize
      rtemp(2)=1.001d0*bxmin
      rtemp(3)=1.001d0*bymin
      rtemp(4)=1.001d0*bzmin
      rtemp(5)=rsgm
      call ga_brdcst(msp_03,itemp,ma_sizeof(mt_int,11,mt_byte),0)
      call ga_brdcst(msp_04,rtemp,ma_sizeof(mt_dbl,5,mt_byte),0)
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
      nseq=itemp(11)
      bsize=rtemp(1)
      bxmin=rtemp(2)
      bymin=rtemp(3)
      bzmin=rtemp(4)
      rsgm=rtemp(5)
      endif
c
      rbbl=rlong+half*rsgm
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      maxbox=max(nbx,nby,nbz)
      nbtot=nbx*nby*nbz
c
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
      msag=max(1,msa,(mbox-30+madbox)*((nwm*nwa+nsa)/nbtot+1)+1)
      msa=max(1,msa,(mbox-30+madbox)*((nwm*nwa+nsa)/nbtot+1)+1)
      mwmg=max(1,msag/mwa+1)
      mwm=max(1,mwm,msa/mwa+1)
c
      msa=min(msa,2*nsa+1)
      mwm=min(mwm,2*nwm+1)
      msag=min(msag,2*nsa+1)
      mwmg=min(mwmg,2*nwm+1)
c
      if(lnode0) then
      msa=nsa+1
      mwm=nwm+1
      msag=nsa+1
      mwmg=nwm+1
      endif
c
      if(util_print('distribution',print_default)) then
      if(me.eq.0) then
      write(lfnout,2006) mbox-30,madbox,nbtot
 2006 format(/,' ARRAY DIMENSION INFORMATION',//,
     + ' Number cells per processor:  ',i7,/,
     + ' Number of buffer cells:      ',i7,/,
     + ' Total number of cells:       ',i7)
      endif
      endif
c
      return
 9999 call md_abort('Failed to open restart file',0)
      return
      end
c $Id$
