      subroutine dia_gettrj
c
c $Id: dia_gettrj.F,v 1.1 2007-01-27 02:26:03 d3j191 Exp $
c
c     read the trajectory
c     -------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      external dia_rdfram
      logical dia_rdfram
c
      integer nread,nframe
c
      read(card(8:37),'(3i10)') ifrfr,ifrto,ifrsk
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
      if(nfrtot.gt.0) then
      if(.not.ga_destroy(ga_trj))
     + call md_abort('Failed to destroy ga_trj',0)
      endif
c
      nfrtot=int((ifrto-ifrfr+1)/ifrsk)
      nfrdim=int((ifrto-ifrfr+1)/ifrsk)/np
      nfrme=int((ifrto-ifrfr+1)/ifrsk)/np
c
      if(me.eq.0) then
      write(*,6001) np,nfrtot,nfrme,nsel
 6001 format(' Number of processors',t40,i5,/,
     + ' Total number of selected frames',t45,i5,/,
     + ' Number of frames per processor',t45,i5,//,
     + ' Number of selected atoms',t45,i5,//)
      endif
c
      if(.not.ga_create(mt_dbl,nsa*3,nfrtot,'trj',nsa*3,nfrdim,
     + ga_trj)) call md_abort('Failed to create ga_trj',0)
c
      write(*,'(i5,a,2i5)') me,' allocated ga for frames ',
     + me*nfrdim+1,min((me+1)*nfrme,nfrtot)
c
c     read trajectory files
c
c     read file header
c
      call dia_rdhdr(byte_mb(i_snam))
      nread=0
      nwrit=0
      nframe=0
 1010 continue
      if(dia_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 1010
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      ndata=ndata+1
      if(me.eq.0) then
      call ga_put(ga_trj,1,3*nsa,ndata,ndata,dbl_mb(i_xdat),nsa)
c      write(*,'(3f12.6)') (dbl_mb(i_xdat+i-1),i=1,3)
      endif
      endif
      nframe=nframe+1
c      if(me.eq.0) then
c      write(*,'(a,i5,a,i5)') ' Frame ',nframe,', Data ',ndata
c      endif
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 1010
      endif
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      return
      end
      subroutine dia_puttrj()
c
c     read the trajectory
c     -------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      return
      end
      subroutine dia_wrttrj(ga_t,fname,fmt,nbatch,logw)
c
c     read the trajectory
c     -------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ga_t,nbatch
      character*255 fname
      real*8 x(msa,3)
      logical logw
c
      integer ifil,ibatch
      character*3 fmt
      integer i,ipd,ifr
c
      if(nbatch.le.0) return
      if(me.ne.0) return
c
      ibatch=0
      ifil=0
      ipd=index(fname,'.')
      if(ipd.eq.0) ipd=index(fname,' ')
      lsonly=.true.
c
      do 1 ifr=1,nfrtot
c
      if(ibatch.eq.0) then
      ifil=ifil+1
      write(filcop,'(a,i5.5,a,a)') fname(1:ipd-1),ifil,'.',fmt
      call dia_wthdr(fmt,byte_mb(i_snam),byte_mb(i_tag),int_mb(i_isel),
     + logw)
      endif
      ibatch=ibatch+1
c
c     write frame
c
      call ga_get(ga_t,1,3*nsa,ifr,ifr,dbl_mb(i_xdat),nsa)
      call dia_wtfram(fmt,byte_mb(i_snam),dbl_mb(i_xdat),dbl_mb(i_wdat),
     + int_mb(i_isel),byte_mb(i_tag),dbl_mb(i_val),int_mb(i_wsel),
     + int_mb(i_ndxw))
c
      if(ibatch.eq.nbatch.or.ifr.eq.nfrtot) then
      ibatch=0
      close(unit=lfncop)
      endif
c
c
    1 continue
c
      return
      end

