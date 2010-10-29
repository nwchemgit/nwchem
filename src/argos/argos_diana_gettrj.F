      subroutine argos_diana_gettrj
c
c $Id$
c
c     read the trajectory
c     -------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      external argos_diana_rdfram
      logical argos_diana_rdfram
c
      integer nread,nframe,itmp(4),ifil,nfram
      integer meio,meframf,meframl,mefilf,mefill,meskip,melast
      integer ifram,iffram,ilfram,lq,ifr
      real*8 t(6)
      character*255 fname
c
c     deallocate any existing global arrays
c
      if(nfrtot.gt.0) then
      if(.not.ga_destroy(ga_trj))
     + call md_abort('Failed to destroy ga_trj',0)
      if(lslvnt) then
      if(.not.ga_destroy(ga_trw))
     + call md_abort('Failed to destroy ga_trw',0)
      endif
      if(.not.ga_destroy(ga_trt))
     + call md_abort('Failed to destroy ga_trt',0)
      endif
c
      read(card(8:38),'(3i10,l1)') ifrfr,ifrto,ifrsk,lslvnt
c
      if(me.eq.0) then
      write(*,'(/,3(a,i10))') ' Selected frames ',ifrfr,' to ',ifrto,
     + ' by',ifrsk
      endif
      timoff=0.0d0
      time=0.0d0
      timr=0.0d0
c
c     allocate global array
c
c     nfrtot : total number of frames
c     nfrdim : maximum number of frames per processor
c     nfrme  : actual number of frames on processor
c     nfrndx : global index of first element local list
c     ntrj   : number of solute atoms per frame in memory resident frames
c
      nfrtot=int((ifrto-ifrfr+1)/ifrsk)
      nfrdim=nfrtot/np
      if(mod(nfrtot,np).ne.0) nfrdim=nfrdim+1
      nfrme=min(nfrdim,nfrtot-me*nfrdim)
      if(nfrme.lt.0) nfrme=0
      nfrndx=me*nfrdim+1
      ntrj=nsa
c
      if(me.eq.0) then
      write(*,6001) np,nfrtot,nfrdim,nsel
 6001 format(/,' Number of processors',t45,i5,/,
     + ' Total number of selected frames',t45,i5,/,
     + ' Number of frames per processor',t45,i5,//,
     + ' Number of selected atoms',t45,i5,//)
      endif
c
      if(.not.ga_create(mt_dbl,ntrj*3,nfrtot,'trj',ntrj*3,nfrdim,
     + ga_trj)) call md_abort('Failed to create ga_trj',0)
      if(lslvnt) then
      if(.not.ga_create(mt_dbl,nwm*nwa*3,nfrtot,'trw',nwm*nwa*3,nfrdim,
     + ga_trw)) call md_abort('Failed to create ga_trw',0)
      endif
      if(.not.ga_create(mt_dbl,6,nfrtot,'trt',6,nfrdim,
     + ga_trt)) call md_abort('Failed to create ga_trt',0)
c
c      write(*,'(i5,a,2i5)') me,' allocated ga for frames ',
c     + me*nfrdim+1,min((me+1)*nfrdim,nfrtot)
c
      if(iomode.eq.0) then
c
c     read trajectory files sequentially
c     ----------------------------------
c
c     read file header
c
      call argos_diana_rdhdr(byte_mb(i_snam))
      nread=0
      nwrit=0
      nframe=0
 1010 continue
      if(argos_diana_rdfram(dbl_mb(i_xdat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 1010
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      ndata=ndata+1
      if(me.eq.0) then
      call ga_put(ga_trj,1,3*nsa,ndata,ndata,dbl_mb(i_xdat),nsa)
      if(lslvnt) call ga_put(ga_trw,1,3*nwm*nwa,ndata,ndata,
     + dbl_mb(i_wdat),nwm*nwa)
      endif
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 1010
      endif
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      else
c
c     read trajectory files in parallel
c     ---------------------------------
c
c     read file header
c
      call argos_diana_rdhdr(byte_mb(i_snam))
      call argos_diana_getofffr()
      iofffr=iofffr-ioffhd
c
      itmp(1)=nfilf
      itmp(2)=ioffhd
      itmp(3)=iofffr
      itmp(4)=0
      call ga_brdcst(mag_d01,itmp,4*ma_sizeof(mt_int,1,mt_byte),0)
      nfilf=itmp(1)
      ioffhd=itmp(2)
      iofffr=itmp(3)
c
c      write(*,'(i5,a,i10)') me,' number of frames per file is ',nfilf
c      write(*,'(i5,a,i10)') me,' trajectory header offset  is ',ioffhd
c      write(*,'(i5,a,i10)') me,' trajectory frame size     is ',iofffr
c
c     meio   : processor carrying out i/o for current processor
c     meframf: first frame to be read from trajectory for this processor
c     meframl: last frame to be read from trajectory for this processor
c     mefilf : file that has first frame for this processor
c     mefill : file that has last frame for this processor
c     meskip : number of frames to be skipped on first file for this processor
c     melast : last frame to be read from the last file for this processor
c
      meio=(me/iomode)*iomode
      meframf=ifrfr+me*nfrme*ifrsk
      meframl=ifrfr+(me+1)*nfrme*ifrsk-1
      mefilf=((meframf-1)/nfilf)+1
      mefill=((meframl-1)/nfilf)+1
      meskip=meframf-(mefilf-1)*nfilf-1
      melast=meframl-(mefill-1)*nfilf
c
cc      write(*,'(i5,a,i5,a,2i5,a,2i5,a,2i5)')
cc     + me,' p :',meio,':',meframf,meframl,':',mefilf,mefill,':',
cc     + meskip,melast
c
c     since a single processor reads the trajectory for iomode processors
c
      if(me.eq.meio) then
c
c     meframl: last frame to be read from trajectory by the i/o processor
c     mefill : file that has last frame for the i/o processor
c     melast : last frame to be read from the last file for the i/o processor
c
      meframl=min(ifrfr+(me+iomode)*nfrme*ifrsk-1,ifrto)
      mefill=((meframl-1)/nfilf)+1
      melast=meframl-(mefill-1)*nfilf
c
c      write(*,'(i5,a,i5,a,2i5,a,2i5,a,2i5)')
c     + me,' n :',meio,':',meframf,meframl,':',mefilf,mefill,':',
c     + meskip,melast
c
c     loop ov er the number of frames to be stored by this processor
c
c
      nfram=me*nfrme
      lq=index(filtrj,'?')
c
c     loop over all files for this processor
c
      do 101 ifil=mefilf,mefill
c
c     construct file name for file ifil
c
      write(cnum,'(i3.3)') ifil
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
c
c     determine first frame to read from this file
c
      iffram=1
      ilfram=nfilf
      if(ifil.eq.mefilf) iffram=meskip+1
      if(ifil.eq.mefill) ilfram=melast
c
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',form='formatted')
      do 102 ifr=iffram,ilfram
      if(mod(ifr-iffram+1,ifrsk).eq.0) then
      nfram=nfram+1
      call fseek(lfntrj,ioffhd+(ifr-1)*iofffr,0)
c      call fseek(lfntrj,ioffhd+(ifr-1)*iofffr,0,9999)
      call argos_diana_preadframe(dbl_mb(i_xdat),dbl_mb(i_wdat),t)
      call ga_put(ga_trj,1,3*nsa,nfram,nfram,dbl_mb(i_xdat),nsa)
      if(lslvnt) call ga_put(ga_trw,1,3*nwm*nwa,ndata,ndata,
     + dbl_mb(i_wdat),nwm*nwa)
      call ga_put(ga_trt,1,6,nfram,nfram,t,6)
      endif
  102 continue
  103 continue
c
      close(unit=lfntrj)
c
  101 continue
      endif
c
      endif
c
      return
 9999 continue
      call md_abort('Error in fseek',0)
      return
      end
