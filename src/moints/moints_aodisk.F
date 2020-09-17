*
* $Id$
*
c
c  Set of routines to cache AO integrals for the 4-index
c  module using exclusive-access local disk storage.
c  Integrals are retrieved in subsequent passes through
c  AO list.
c
c  Maintain state via two common block variables.
c
c       moao_ipass  < 0      (default) disk caching disabled
c                   = 0      saving to disk, initial pass
c                   > 0      retrieving from disk
c
c       moao_fd     <  0     io not initialized
c                   >= 0     io is happening
c
c  Unless explicitly initialized, there will be no disk
c  caching.
c
c  Since exclusive-access files are used, the same subset of 
c  AO integrals that was initially generated on a processor 
c  will be retrieved in subsequent passes. Thus, the task scheduling 
c  from the initial run must be emulated.
c  Use a wrapper around nxtask() and record the task number on
c  each disk record.
c
c
c
c
c
c  Following routines use the common block
c
      block data moints_moao_block
      implicit none
#include "cmointsmoao.fh"
      data moao_ipass/-1/
      data moao_fd/-1/
      data moao_lbuf/-1/
      data moao_fname/' '/
      end





c
c  Wrapper around nxtask so we can save task numbers
c  with the IO records and emulate task scheduling
c  when retrieved
c
      integer function moints_nxttask( numnodes )
      implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "cmointsmoao.fh"
      integer numnodes
      integer nxtask
      external nxtask
      
      if (moao_ipass.gt.0) then
        moints_nxttask = moao_tasknum
        if (numnodes.lt.0) call ga_sync()
      else
        moints_nxttask = nxtask( numnodes, 1 )
        if (moao_ipass.eq.0) then
          moao_tasknum = moints_nxttask
        endif
      endif
      return
      end

      

c
c  Return complete SSBB block of integrals if cached on disk
c
      logical function moints_gblk_fromdisk( blkid, ish, jsh, 
     $                                       kshlo, lshlo,
     $                                       ilo, ihi, jlo, jhi,
     $                                       kblo, kbhi, lblo, lbhi,
     $                                       ssbb )
      implicit none
#include "errquit.fh"
      integer blkid
      integer ish, jsh
      integer kshlo
      integer lshlo
      integer ilo, ihi
      integer jlo, jhi
      integer kblo, kbhi
      integer lblo, lbhi
      double precision ssbb( lblo:lbhi, kblo:kbhi, jlo:jhi, ilo:ihi )
c
#include "mafdecls.fh"
#include "cmointsmoao.fh"
c
#ifdef OLD_AODISK
      character*8 buidstr
      integer recnum
      logical moints_iorec_next
      external moints_iorec_next
#endif
      integer ssbblen
c$$$      DOUBLE PRECISION DABSSUM
c$$$      EXTERNAL DABSSUM

      moints_gblk_fromdisk = .false.
      if (moao_ipass.lt.0) return
      if ((moao_fd.ge.0).and.(moao_ipass.gt.0)) then

#ifdef OLD_AODISK
        recnum = 0
        ssbblen = (ihi-ilo+1)*(jhi-jlo+1)*(kbhi-kblo+1)*(lbhi-lblo+1)
        call dfill( ssbblen, 0.d0, ssbb, 1 )
        do while (moints_iorec_next( blkid, buidstr  ))
          if (moao_issparse.eq.1) then
            call moints_aodisk_iorec2sprs_old( ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, moao_reclen, moao_lwidth,
     $                             dbl_mb(moao_klabrec),
     $                             dbl_mb(moao_kvalrec) )
            moints_gblk_fromdisk = .true.
          else
            call errquit('moints_gblk_disk: dense write not ready',0,
     &       INT_ERR)
          endif
          recnum = recnum + 1
        enddo
#else
        ssbblen = (ihi-ilo+1)*(jhi-jlo+1)*(kbhi-kblo+1)*(lbhi-lblo+1)
        call dfill( ssbblen, 0.d0, ssbb, 1 )
        call moints_aodisk_iorec2sprs( moao_fd, moao_fptr, 
     $                            moao_tasknum,
     $                            ilo, ihi, jlo, jhi,
     $                            kblo, kbhi, lblo, lbhi,
     $                            ssbb, moao_hdrp, moao_buflen,
     $                            dbl_mb(moao_kbuf), moao_nrec )
        moints_gblk_fromdisk = .true.
#endif
C        MOAO_CUMUL = MOAO_CUMUL + DABSSUM( SSBBLEN, SSBB )
c$$$       WRITE(6,771) ISH, JSH, KSHLO, LSHLO, MOAO_TASKNUM, GA_NODEID(),
c$$$     $              BLKID, DABSSUM(SSBBLEN,SSBB)
c$$$ 771   FORMAT('qqqq-',4I4,3X,I5,I3,5X,I6,5X,F20.6)
      endif
      return
      end




c
c  Save SSBB block to disk if caching enabled.
c
      logical function moints_gblk_todisk( blkid, ish, jsh, 
     $                                     kshlo, lshlo,
     $                                     ilo, ihi, jlo, jhi,
     $                                     kblo, kbhi, lblo, lbhi,
     $                                     ssbb )
      implicit none
#include "errquit.fh"
      integer blkid
      integer ish, jsh
      integer kshlo, lshlo
      integer ilo, ihi
      integer jlo, jhi
      integer kblo, kbhi
      integer lblo, lbhi
      double precision ssbb( lblo:lbhi, kblo:kbhi, jlo:jhi, ilo:ihi )
c
#include "mafdecls.fh"
#include "cmointsmoao.fh"
c
c
      logical sparse                    ! for the moment only sparse case
      data sparse/.true./
c
c$$$      INTEGER SSBBLEN
c$$$      DOUBLE PRECISION DABSSUM
c$$$      EXTERNAL DABSSUM

      moints_gblk_todisk = .false.
      if (moao_ipass.lt.0) return
      if ((moao_fd.ge.0).and.(moao_ipass.eq.0)) then
        if (sparse) then
#ifdef OLD_AODISK
          call moints_aodisk_sprs2iorec_old( moao_fd, moao_fptr,
     $                           moao_tasknum, blkid,
     $                           ilo, ihi, jlo, jhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, moao_spreclen, moao_lwidth,
     $                           dbl_mb(moao_klabrec),
     $                           dbl_mb(moao_kvalrec),
     $                           moao_buflen, dbl_mb(moao_kbuf) )
#else
          call moints_aodisk_sprs2iorec( moao_fd, moao_fptr,
     $                           moao_tasknum, blkid, 
     $                           ilo, ihi, jlo, jhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, moao_recptr, moao_buflen, 
     $                           dbl_mb(moao_kbuf), moao_nrec )
#endif
        else
          call errquit('moints_gblk_todisk: dense read not ready',0,
     &       INT_ERR)
        endif
       moints_gblk_todisk = .true.
c
C        MOAO_CUMUL = MOAO_CUMUL + DABSSUM( SSBBLEN, SSBB )
c$$$       SSBBLEN = (IHI-ILO+1)*(JHI-JLO+1)*(KBHI-KBLO+1)*(LBHI-LBLO+1)
c$$$       WRITE(6,771) ISH, JSH, KSHLO, LSHLO, MOAO_TASKNUM, GA_NODEID(),
c$$$     $              BLKID, DABSSUM(SSBBLEN,SSBB)
c$$$ 771   FORMAT('qqqq+',4I4,3X,I5,I3,5X,I6,5X,F20.6)
c
      endif
      return
      end






c     
c     Called by application
c     Enable disk caching & check for existing file and open file
c     
      logical function moints_aodisk_init( odisk, oreuse )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      logical odisk             ! [input] toggle caching
      logical oreuse            ! [input] toggle reuse of existing file
c     
      logical fexist
      integer stat
c     
c     first time through --- enable caching TO disk
c     otherwise must be caching FROM disk
c     
      moints_aodisk_init = .false.
      if (.not.(odisk)) then
         moao_ipass = -1
         moints_aodisk_init = .true.
         return
      endif
      if (moao_ipass.lt.0) moao_ipass = 0
c     
c     if AO file exists AND reuse enabled then set ipass > 0
c     
      if (moao_fname.eq.' ') then
         call util_file_name( 'moao', .true., .true., moao_fname )
      endif
c
      inquire(file=moao_fname, exist=fexist) ! Equivalent is eaf_stat ?
      if (oreuse) then
         stat=0
         if (fexist) stat=1
! igop('*') detects the presence of one or more stat=0 inputs
#if NWCHEM_USE_IGOP_PROD
         call ga_igop( 481, stat, 1, '*')
#else
         call ga_igop( 481, stat, 1, 'min')
#endif
         if (stat.eq.1) then
            moao_ipass = 1
            if (ga_nodeid().eq.0) write(6,331)
 331        format(10x,'Existing AO integral file will be reused')
         else if (fexist) then
            call util_file_unlink(moao_fname) ! Only some nodes had it
         endif
      else if (fexist) then
         call util_file_unlink(moao_fname) ! No reuse but file there.
      endif
c     
c     open AO integral disk cache
c     
      if (moao_fd.lt.0) then
         stat = eaf_open( moao_fname, EAF_RW, moao_fd )
         if (stat.ne.0) call errquit(
     $        'moints_aodisk_init: cannot open ao file',0, DISK_ERR)
      endif
      moints_aodisk_init = .true.
      return
      end







c
c  This is called internally by moints
c  Initialize disk caching for both saving and retrieving
c  Allocate buffer space and read in first record
c
      logical function moints_aodisk_prep( )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      integer stat
      integer blkinfo(MOINTS_NBLKINFO)
      integer vlsize, dhdrsize, vallen, bufbytes
      integer nblk,nexpand
      double precision sparsesize, densesize
      common/aodisk_stat/nblk,nexpand,sparsesize,densesize

      moints_aodisk_prep = .false.
      if (moao_ipass.lt.0) then
        moints_aodisk_prep = .true.
        return
      endif
      moao_nrec = 1
c
c  allocate io buffer
c  header comes first
c  values and labels start after offset of header length
c
      if (moao_lbuf.lt.0) then
        moao_lwidth = 8/ma_sizeof(MT_INT, 1, MT_BYTE)
        vlsize = ma_sizeof(MT_INT, moao_lwidth, MT_BYTE) + 
     $           ma_sizeof(MT_DBL, 1, MT_BYTE)                              ! width of value + labels
        dhdrsize = ma_sizeof(MT_INT, MOINTS_NBLKINFO, MT_DBL )
        if (.not.ma_alloc_get(MT_DBL, MOINTS_IOBUFLEN, 'moints iobuf',
     $                        moao_lbuf, moao_kbuf)) call errquit(
     $       'moints_aodisk_prep: cannot allocate io buffer',0,
     &       DISK_ERR)
        moao_buflen   = MOINTS_IOBUFLEN                                     ! bufflen (double words)
        vallen        = moao_buflen - dhdrsize                              ! value len -= hdr length
        moao_spreclen = vallen/(vlsize/ma_sizeof(MT_DBL,1,MT_BYTE))         ! leng
        moao_kvalrec  = moao_kbuf + dhdrsize                                ! values after header
        moao_klabrec  = moao_kbuf + dhdrsize + moao_spreclen                ! labels after values
      endif
#ifndef OLD_AODISK
c
c  Initialize first IO record 
c
      if (moao_ipass.eq.0) then
        moao_recptr   = 2
        call icopy(1,moao_recptr,1,dbl_mb(moao_kbuf),1)
      endif
#endif
c
c  rewind AO integral cache unit
c
      moao_fptr = 1.d0
      moao_eof = .false.
c
c  if READ mode -- retrieve first record
c
      if (moao_ipass.gt.0) then
        bufbytes = ma_sizeof(MT_DBL,moao_buflen,MT_BYTE)
        stat =  eaf_read( moao_fd, moao_fptr, dbl_mb(moao_kbuf), 
     $                    bufbytes )
        moao_fptr = moao_fptr + bufbytes
#ifdef OLD_AODISK
        if (stat.eq.0) then
          call icopy(MOINTS_NBLKINFO, dbl_mb(moao_kbuf), 1, blkinfo, 1)
          moao_touch   = .false.
          moao_tasknum  = blkinfo(1)
          moao_reclen   = blkinfo(2)
          moao_issparse = blkinfo(3)
          moao_blkid    = blkinfo(4)
c$$$          WRITE(6,324) moao_tasknum, moao_reclen, moao_blkid
c$$$ 324      FORMAT(' First record: tasknum=',I5,'   RecLen=',I5,
c$$$     $           '    Blkid=',I5)
        else
          moao_ipass = 0                ! this is peculiar but not an error, just reset
        endif
#else
        if (stat.eq.0) then
          call icopy(MOINTS_NBLKINFO,dbl_mb(moao_kbuf+1),1,blkinfo,1)
          moao_touch    = .false.
          moao_tasknum  = blkinfo(4)
          moao_blkid    = blkinfo(5)
          moao_issparse = blkinfo(6)
          moao_recptr   = 2
          moao_hdrp     = 2
        else
          moao_ipass    = 0
        endif
#endif
      endif
      moints_aodisk_prep = .true.
      MOAO_CUMUL = 0.d0
      NEXPAND = 0
      NBLK = 0
      SPARSESIZE = 0.d0
      DENSESIZE = 0.d0
      return
      end




c
c  Called internally by moints between passes
c  Free IO buffer and increment pass count
c  Keep NO state information except moao_ipass and moao_fd
c  between passes.
c
      subroutine moints_aodisk_tidy()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      double precision fcompress
      integer nblk, nexpand
      double precision sparsesize, densesize
      common/aodisk_stat/nblk,nexpand,sparsesize,densesize

      if (moao_ipass.lt.0) return
#ifndef OLD_AODISK
      if (moao_ipass.eq.0) then
        call moints_aodisk_flushiorec( moao_fd, moao_fptr,                    ! flush last record
     $                      moao_buflen, dbl_mb(moao_kbuf) )     

        fcompress = (sparsesize/densesize)*100.d0
        if (util_print('compress stats',print_high).and.
     $    ga_nodeid().eq.0) then
          write(6,911) nblk, nexpand, fcompress
 911      format(' Total blocks        = ',I8,/,
     $           ' Expanded blocks     = ',I8,/,
     $           ' Compression factor  = ',F8.1,'%')
        endif

c$$$        CALL MOINTS_AODISK_CHECKFILE( MOAO_FD, MOAO_BUFLEN,       ! verify file for debugging
c$$$     $                             DBL_MB(MOAO_KBUF) )
      endif
#endif
      if (moao_lbuf.ge.0) then
        if (.not.ma_free_heap(moao_lbuf)) call errquit(
     $    'moints_closeaodisk: cannot free io buffer',0, MA_ERR)
        moao_lbuf = -1
      endif
      moao_ipass = moao_ipass + 1

c$$$      WRITE(6,991) MOAO_IPASS, MOAO_CUMUL, MOAO_NREC
c$$$ 991  format(' File pass=',i3,/,
c$$$     $       ' Cumulative value:',f20.6,
c$$$     $       ' Records read:',i5)
      
      return
      end




c
c  Called by application
c  Close I/O unit with option to save
c  Reset ipass count
c
      subroutine moints_aodisk_close( osave ) 
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "util.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      logical osave
c
      integer stat

      if (moao_ipass.lt.0) return
      moao_ipass = 0
      if (moao_fd.ge.0) then
        if (util_print('ao disk stats',print_high))
     $     call eaf_print_stats(moao_fd)
        stat = eaf_close(moao_fd)
        if (stat.ne.0) call errquit(
     $    'moints_closeaodisk: cannot close file',0, DISK_ERR)
        if (.not.(osave)) then
          stat = eaf_delete(moao_fname)
          if (stat.ne.0) call errquit(
     $    'moints_closeaodisk: cannot delete file',0, DISK_ERR)
        endif
        moao_fd = -1
      endif
      
      return
      end








c
c ======================================================================
c
c
c                    Utility routines
c
c
c ======================================================================
c
c
#ifdef OLD_AODISK
c
c
c  Pack the dense SSBB block into sparse form
c  with 16 bits per label
c
c
      subroutine moints_aodisk_sprs2iorec_old( fd, fptr, tasknum, blkid,
     $                             ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, reclen, lwidth,
     $                             iolab, ioval,
     $                             iobuflen, iobuf )
      implicit none
      integer fd
      double precision fptr
      integer tasknum
      integer blkid
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer reclen
      integer lwidth
      integer iolab(lwidth,reclen)
      double precision ioval(reclen)
      integer iobuflen
      double precision iobuf(*)

      integer i, j, k, l
      integer iir, issparse, lab1, lab2
      integer recnum
      double precision xx
#include "bitops_decls.fh"
#include "bitops_funcs.fh"      

      issparse = 1
      recnum = 0
      iir = 0
      call ifill( lwidth*reclen, 0, iolab, 1 )
      call dfill( reclen, 0.d0, ioval, 1 )
      do i=ilo,ihi
        do j=jlo,jhi
          lab1 = ior(ishft(i,16),j)
          do k=kblo,kbhi
            do l=lblo,lbhi
              xx = ssbb(l,k,j,i)
              if (abs(xx).gt.1.d-12) then
                iir = iir + 1
                ioval(iir)   = xx
                lab2 = ior(ishft(k,16),l)
                if (lwidth.eq.1) then
                  iolab(1,iir) = ior(ishft(lab1,32),lab2)
                elseif (lwidth.eq.2) then
                  iolab(1,iir) = lab1
                  iolab(2,iir) = lab2
                endif
                if (iir.eq.reclen) then
                  call moints_iorec_flush(fd, fptr, tasknum, iir, 
     $                                    issparse, blkid, iobuflen, 
     $                                    iobuf )
                  recnum = recnum + 1
                  call ifill( lwidth*reclen, 0, iolab, 1 )
                  call dfill( reclen, 0.d0, ioval, 1 )
                  iir = 0
                endif
              endif
            enddo  
          enddo  
        enddo  
      enddo
      if (iir.gt.0) then
        call moints_iorec_flush( fd, fptr, tasknum, iir, issparse,
     $                          blkid, iobuflen, iobuf )
        recnum = recnum + 1
        iir = 0
      endif
                
      return
      end

      

      subroutine moints_aodisk_iorec2sprs_old( ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, reclen, lwidth, iolab,
     $                             ioval )
      implicit none
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer reclen
      integer lwidth
      integer iolab(lwidth,reclen)
      double precision ioval(reclen)
      integer q, i, j, k, l
      integer i16mask
      integer onbitmask
      external onbitmask
#include "bitops_decls.fh"
#include "bitops_funcs.fh"      
      
      i16mask = onbitmask(16)
      if (lwidth.eq.1) then
        do q=1,reclen
          i = iand(ishft(iolab(1,q),-48),i16mask)
          j = iand(ishft(iolab(1,q),-32),i16mask)
          k = iand(ishft(iolab(1,q),-16),i16mask)
          l = iand(iolab(1,q),i16mask)
          ssbb(l,k,j,i) = ioval(q)
        enddo
      else if (lwidth.eq.2) then
        do q=1,reclen
          i = iand(ishft(iolab(1,q),-16),i16mask)
          j = iand(iolab(1,q),i16mask)
          k = iand(ishft(iolab(2,q),-16),i16mask)
          l = iand(iolab(2,q),i16mask)
          ssbb(l,k,j,i) = ioval(q)
        enddo
      endif
      return
      end






      subroutine moints_iorec_flush( fd, fptr, tasknum, reccnt, 
     $                               issparse, blkid, buflen, iobuf )
      implicit none
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
      integer fd
      double precision fptr
      integer tasknum
      integer reccnt
      integer issparse
      integer blkid
      integer buflen
      double precision iobuf(buflen)
c
      integer blkinfo(MOINTS_NBLKINFO)
      integer stat, bufbytes
c
      call ifill( MOINTS_NBLKINFO, 0, blkinfo, 1 )
      blkinfo(1) = tasknum
      blkinfo(2) = reccnt
      blkinfo(3) = issparse
      blkinfo(4) = blkid
      call icopy( MOINTS_NBLKINFO, blkinfo, 1, iobuf, 1 )
      bufbytes = ma_sizeof(MT_DBL,buflen,MT_BYTE)
      stat = eaf_write(fd, fptr, iobuf, bufbytes )
      fptr = fptr + bufbytes
c
c$$$      WRITE(6,772) blkid
c$$$ 772  FORMAT('---->',I8)
c

      return
      end

      
      



      logical function moints_iorec_next( blkid, buidstr )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      integer blkid
      character*8 buidstr
      integer blkinfo(MOINTS_NBLKINFO)
      integer stat
      integer bufbytes

      moints_iorec_next = .false.
      if (moao_eof) return
c
c  check if io buffer is already pre-read
c
      if (.not.(moao_touch)) then
        if (blkid.eq.moao_blkid) then
          moao_touch = .true.
          moints_iorec_next = .true.
        endif
        return
      endif
c
c  otherwise read in a buffer
c
      bufbytes = ma_sizeof(MT_DBL,moao_buflen,MT_BYTE)
      stat =  eaf_read( moao_fd, moao_fptr, dbl_mb(moao_kbuf), 
     $                  bufbytes )
      moao_fptr = moao_fptr + bufbytes
      MOAO_NREC = MOAO_NREC + 1
c
c  check for EOF and/or update info
c
      if (stat.eq.0) then
        call icopy( MOINTS_NBLKINFO, dbl_mb(moao_kbuf), 1, blkinfo, 1 )
        moao_reclen   = blkinfo(2)
        moao_issparse = blkinfo(3)
        moao_touch    = (moao_tasknum.eq.blkinfo(1)).and.
     $                  (blkid.eq.blkinfo(4))
        if (moao_touch) then
          moints_iorec_next = .true.
          return
        endif
        moao_tasknum  = blkinfo(1)
        moao_blkid    = blkinfo(4)
c$$$        WRITE(6,772) moao_blkid
c$$$ 772    FORMAT('<----',I8)
      elseif (stat.gt.0) then
        call errquit('moints_io disk io error',0, DISK_ERR)
      endif
c
c  reach here if tasknumber has changed or EOF
c  return FALSE
c
      moao_eof = (stat.lt.0)
      return
      end
        
      
c
c  Create a unique identifier string for
c  4 shell labels - 2 bytes per label
c
      subroutine moints_uid( str, i, j, k, l)
      implicit none
      character*8 str
      integer i, j, k, l
      character*2 ci, cj, ck, cl
      integer*2 ii, jj, kk, ll
      equivalence (ci,ii),(cj,jj),(ck,kk),(cl,ll)
      
      ii = i
      jj = j
      kk = k
      ll = l
      str(1:2) = ci
      str(3:4) = cj
      str(5:6) = ck
      str(7:8) = cl

      return
      end




c
c end OLD_AODISK section
c
#endif      
                         



c
c       ==============================================
c       ==============================================
c
c
c
c                      New Version 
c
c
c
c       ==============================================
c       ==============================================
c
c





#ifndef OLD_AODISK
c
c  Alternate version of sparse packing
c  Labels/indices require integer of storage, 2*NCOL + NNZ/2,
c  cf. 2*NNZ for standard version.
c
c
c  IO Record structure:
c    The first integer of each IO record points to the 
c    (double word) location of the first header. The space 
c    between is assumed to the continuation of the 
c    previous entry. ie.
c
c       __________________________
c       |                        |
c       |                        v
c       --------------------------------------------
c       |                        | header....
c       --------------------------------------------
c         |<-- continuation -->|
c              prev record
c
c
c  The header structure is:
c
c         hdr(1)    magic I    valid header
c                   magic II   invalid header, read next record (rare)
c
c         hdr(2)    >0         location of next header
c                   -1         last header for this IO record
c
c
c

      subroutine moints_aodisk_sprs2iorec( fd, fptr, tasknum, blkid,
     $                                     ilo, ihi, jlo, jhi,
     $                                     kblo, kbhi, lblo, lbhi,
     $                                     ssbb, recptr, reclen, iorec,
     $                                     nrec )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
      integer fd
      double precision fptr
      integer tasknum
      integer blkid
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer recptr
      integer reclen
      double precision iorec(reclen)
      integer nrec
c
      integer i, j, ilen, jlen, llen, klen
      integer swordlen, nir, nnz, dbilen, dijlen, nir1, nnz1
      integer rp, hdrp, ijp, rleft, bblen, nssbb
      integer irp, irrp, drp, hoff
      integer ijidx, issparse
      integer l_ij, k_ij
      integer blkinfo(MOINTS_NBLKINFO)
      integer iminus, magic_cookie, magic_cookie2
      logical iscont, ojustfit, ohdrwrite, ojust1blk
      integer labint
      double precision sprstol
      integer nblk,nexpand
      double precision sparsesize, densesize
      common/aodisk_stat/nblk,nexpand,sparsesize,densesize
c
      integer onbitmask
      external onbitmask
      integer dnonzero_cnt
      external dnonzero_cnt
c      
#include "bitops_decls.fh"

*:: data statements must come after all declarations

      data issparse/1/
      data iminus/-1/
      data sprstol/1.d-12/

#include "bitops_funcs.fh"      
c
      blkinfo(1)  = -1
      blkinfo(3)  = 0
      blkinfo(4)  = tasknum
      blkinfo(5)  = blkid
      blkinfo(6)  = issparse
      blkinfo(7)  = ilo
      blkinfo(8)  = ihi
      blkinfo(9)  = jlo
      blkinfo(10) = jhi
      blkinfo(11) = kblo 
      blkinfo(12) = kbhi
      blkinfo(13) = lblo
      blkinfo(14) = lbhi
      ohdrwrite = .false.
c
c
C      WRITE(6,910) BLKINFO(4),  BLKINFO(5),  BLKINFO(7),  
C     $             BLKINFO(8),  BLKINFO(9),  BLKINFO(10), 
C     $             BLKINFO(11), BLKINFO(12), BLKINFO(13), 
C     $             BLKINFO(14)
C 910  FORMAT('Task#:',I5,5x,' ID =',I5,3x,4(i2,'-',i2,3x))
c
      llen = lbhi - lblo + 1
      klen = kbhi - kblo + 1
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      swordlen = ma_sizeof(MT_DBL,1,MT_BYTE) + 2        ! double word + 2 byte label = 10
      magic_cookie = onbitmask(17)
      magic_cookie2 = onbitmask(21)
      labint = ma_sizeof(MT_INT,1,MT_BYTE)/2            ! labels per integer word, 2 bytes per label
c
c  Minimium header length must be remaining on IO record
c  Otherwise flush IO record
c
      dbilen = ma_sizeof(MT_INT,MOINTS_NBLKINFO,MT_DBL)
      dijlen = ma_sizeof(MT_INT,(3*ilen*jlen),MT_DBL)
      rleft = reclen - recptr + 1
      if ((dbilen+dijlen).gt.rleft) then
        call moints_aodisk_flushiorec( fd, fptr, reclen, iorec )
        recptr = 2
        call icopy( 1, recptr, 1, iorec(1), 1)
      endif
c
c  Header and record pointers 
c
      ojust1blk = recptr.eq.2
      hdrp  = recptr
      ijp   = hdrp + dbilen
      rp    = hdrp + dbilen + dijlen
      nssbb = 0
c
c  Allocate temp space
c
      if (.not. ma_push_get(MT_INT, (ilen*jlen), 'ij ',l_ij, k_ij))
     $  call errquit('moints_aodisk_sprs2iorec: no memory',0, MA_ERR)
      call ifill((ilen*jlen), 0, int_mb(k_ij), 1 )
c
c
c
      iscont = .false.
      do i=ilo,ihi
        do j=jlo,jhi

c
c  Count non-zeros -- compute space requirements 
c
          rleft = reclen - rp + 1
          nnz = dnonzero_cnt((llen*klen),sprstol,ssbb(lblo,kblo,j,i))
          nir = nnz/labint
          if (mod(nnz,labint).ne.0) nir = nir + 1
          bblen = ma_sizeof(MT_INT,(2*klen),MT_DBL) +
     $            ma_sizeof(MT_INT,nir,MT_DBL) + nnz
c
c  Compression statistics
c
          nblk = nblk + 1
          if (bblen.gt.(klen*llen)) nexpand = nexpand + 1
          sparsesize = sparsesize + bblen
          densesize = densesize + llen*klen
c
c  Flush IO record if insufficient space
c 
          
          if (bblen.gt.rleft) then
            if (.not.(ohdrwrite)) then
              blkinfo(1) = magic_cookie
              blkinfo(2) = -1                                             ! flag -- last entry in this record
              call icopy( MOINTS_NBLKINFO, blkinfo, 1, iorec(hdrp), 1 )
              ohdrwrite = .true.
            endif
            call icopy( (ilen*jlen), int_mb(k_ij), 1, iorec(ijp), 1 )
            call moints_aodisk_flushiorec( fd, fptr, reclen, iorec )
c
c  Reset pointers and arrays
c
            call icopy( 1, iminus, 1, iorec(1), 1)
            call ifill((ilen*jlen), 0, int_mb(k_ij), 1 )
            ijp = 2
            rp  = ijp + dijlen
            nssbb = 0
            iscont = .true.
            ojust1blk = .true.
          endif
c
c  Pack sparse structure into IO record
c    
          ijidx = (i-ilo)*jlen + j - jlo
          int_mb(k_ij+ijidx) = ior(ishft(rp,16),nnz)
          nssbb = nssbb + 1
          irp  = rp
          irrp = irp  + ma_sizeof(MT_INT,(2*klen),MT_DBL)
          drp  = irrp + ma_sizeof(MT_INT,nir,MT_DBL)
          call moints_sparse2d_pack( llen, klen, ssbb(lblo,kblo,j,i),
     $                               sprstol, labint, iorec(irp), nir1, 
     $                               iorec(irrp), iorec(drp), nnz1 )
          rp = drp + nnz1
          if ((nnz1.ne.nnz).or.(nir1.ne.nir)) 
     $      call errquit('moints_aodisk_sprs2iorec: internal error 2',0,
     &       DISK_ERR)
        enddo

      enddo
      if (rp.gt.(reclen+1))
     $  call errquit('moints_aodisk_sprs2iorec: internal error 1',0,
     &       DISK_ERR)
c
c  Bookkeeping for some odd cases
c
      ojustfit = (rp.eq.(reclen+1))                     ! special case
      if (ojustfit) then
        rp = 0
        if (ojust1blk) then
          hoff = 0
          if (.not.(ohdrwrite)) hoff = hdrp
          call icopy( 1, hoff, 1, iorec(1), 1 )
        endif
      endif
c
c  Copy header info & map to IO record before returning
c
      if (.not.(ohdrwrite)) then
        blkinfo(1) = magic_cookie
        blkinfo(2) = rp
        call icopy( MOINTS_NBLKINFO, blkinfo, 1, iorec(hdrp), 1 )
      endif
      call icopy( (ilen*jlen), int_mb(k_ij), 1, iorec(ijp), 1 )
c
c  Special case -- flush IO record
c
      if (ojustfit) then
        call moints_aodisk_flushiorec( fd, fptr, reclen, iorec)
        rp = 2
        call icopy( 1, rp, 1, iorec(1), 1)
      endif
c
c  Put magic cookie II at start of next header.
c  Usually overwritten by magic cookie next time around 
c  (unless we skip to next record)
c
      if (rp.le.reclen) then
        call icopy(1, magic_cookie2, 1, iorec(rp), 1 )
      endif
c
c  If SSBB block spans multiple records must encode offset
c  of end of block
c
      if (iscont) then
        call icopy( 1, rp, 1, iorec(1), 1 )
      endif
      recptr = rp
c
c  Clean up
c
      if (.not. ma_pop_stack(l_ij))
     $  call errquit('moints_aodisk_sprs2iorec: failed to pop', l_ij,
     &       MA_ERR)
      
      return
      end






c
c  Read and interpret AO file contents - one block
c  at a time -- structure as above
c
c
      subroutine moints_aodisk_iorec2sprs( fd, fptr, tasknum,
     $                                   ilo, ihi, jlo, jhi,
     $                                   kblo, kbhi, lblo, lbhi,
     $                                   ssbb, hdrp, reclen, iorec,
     $                                   nrec )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
      integer fd
      double precision fptr
      integer tasknum
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer hdrp
      integer reclen
      double precision iorec(reclen)
      integer nrec
c
      integer cookie, cookie2
      integer dbilen, ilen, jlen, ijp
      integer nexthdrp, stat, bufbytes, magic
      INTEGER IJ
      integer l_ij, k_ij
      integer blkinfo(MOINTS_NBLKINFO)
      logical ocont
      integer onbitmask
      external onbitmask
c
c
c
      cookie  = onbitmask(17)
      cookie2 = onbitmask(21)
      bufbytes = ma_sizeof(MT_DBL,reclen,MT_BYTE)
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      dbilen = ma_sizeof(MT_INT,MOINTS_NBLKINFO,MT_DBL)
      ijp = hdrp + dbilen
      if (.not. ma_push_get(MT_INT, (ilen*jlen), 'ij ',l_ij, k_ij))
     $  call errquit('moints_aodisk_iorec2sprs: no memory',0,
     &       MA_ERR)
c
c  Sanity check on current pointer and IO record
c
      call icopy( MOINTS_NBLKINFO, iorec(hdrp), 1, blkinfo, 1)
      if (blkinfo(1).ne.cookie) 
     $  call errquit('moints_aodisk_iorec2sprs: AOfile corrupt I?',0,
     &       DISK_ERR)

c$$$      WRITE(6,910) BLKINFO(4),  BLKINFO(5),  BLKINFO(7),  
c$$$     $             BLKINFO(8),  BLKINFO(9),  BLKINFO(10), 
c$$$     $             BLKINFO(11), BLKINFO(12), BLKINFO(13), 
c$$$     $             BLKINFO(14)
c$$$ 910  FORMAT('Task#:',I5,5x,' ID =',I5,3x,4(i2,'-',i2,3x))
c
c  Check consistency between args and blockinfo
c
      if ((ilo.ne.blkinfo(7)).or.(ihi.ne.blkinfo(8)).or.
     $    (jlo.ne.blkinfo(9)).or.(jhi.ne.blkinfo(10)).or.
     $    (kblo.ne.blkinfo(11)).or.(kbhi.ne.blkinfo(12)).or.
     $    (lblo.ne.blkinfo(13)).or.(lbhi.ne.blkinfo(14))) then
        write(6,623) nrec,(blkinfo(ij),ij=7,14)
 623    format(' Rec#',i6,3x,' B ',8i6)
        write(6,622) nrec,ilo,ihi,jlo,jhi,kblo,kbhi,lblo,lbhi
 622    format(' Rec#',i6,3x,' E ',8i6)
        call errquit('moints_aodisk_sprs2iorec: blk label mismatch',0,
     &       INPUT_ERR)
      endif
c
c  Unpack data structure from IO record(s)
c  -- possibly spanning multiple records
c
      ocont = .true.
      nexthdrp = blkinfo(2)

      do while (ocont)

        call icopy( (ilen*jlen), iorec(ijp), 1, int_mb(k_ij), 1)

c$$$        CALL MOINTS_AODISK_IJPRINT( ILEN, JLEN, INT_MB(K_IJ) )

        call moints_aodisk_sprs2dense_a( ilo, ihi, jlo, jhi,
     $                                   kblo, kbhi, lblo, lbhi,
     $                                   int_mb(k_ij), reclen, iorec,
     $                                   ssbb )
        ocont = nexthdrp.eq.-1
        if (nexthdrp.le.0) then
          stat = eaf_read(fd, fptr, iorec, bufbytes )
          nrec = nrec + 1
          if (stat.ne.0) goto 150
          fptr = fptr + bufbytes
          call icopy( 1, iorec(1), 1, nexthdrp, 1 )
          ijp  = 2
        endif
      enddo
c
c  Check next block -- skip ahead if invalid
c
 101  call icopy( 1, iorec(nexthdrp), 1, magic, 1)
      if (magic.eq.cookie2) then
        stat = eaf_read(fd, fptr, iorec, bufbytes )
        nrec = nrec + 1
        if (stat.ne.0) goto 150
        fptr = fptr + bufbytes
        nexthdrp = 2
        goto 101
      endif
      if (magic.ne.cookie) 
     $  call errquit('moints_aodisk_iorec2sprs: AOfile corruptII?',0,
     &       DISK_ERR)
c
c  Prefetch tasknumber from header
c
      call icopy( 4, iorec(nexthdrp), 1, blkinfo, 1 )
      tasknum = blkinfo(4)
      hdrp = nexthdrp
c
c  Cleanup
c
 150  continue
      if (.not. ma_pop_stack(l_ij))
     $  call errquit('moints_aodisk_iorec2sprs: failed to pop', l_ij,
     &       MA_ERR)
      return
      end










      subroutine moints_aodisk_flushiorec( fd, fptr, n, r )
      implicit none
#include "mafdecls.fh"
#include "eaf.fh"
#include "cmointsmoao.fh"
      integer fd
      double precision fptr
      integer n
      double precision r(n)
      integer bufbytes
      integer stat

      moao_nrec = moao_nrec + 1
      bufbytes = ma_sizeof(MT_DBL,n,MT_BYTE)
      stat = eaf_write(fd, fptr, r, bufbytes )
      fptr = fptr + bufbytes
      return
      end





      subroutine moints_aodisk_sprs2dense_a( ilo, ihi, jlo, jhi,
     $                                      kblo, kbhi, lblo, lbhi,
     $                                      ijmap, reclen, iorec, ssbb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      integer ijmap(jlo:jhi,ilo:ihi)
      integer reclen
      double precision iorec(*)
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
c
      integer labint
      integer i, j, ijp, nnz, bit16, ilen, jlen, klen, llen
      integer maxiirlen, iirlen, nir, dp
      integer l_iv, k_iv, k_irv
      integer onbitmask
      external onbitmask
#include "bitops_decls.fh"
#include "bitops_funcs.fh"      

      ilen = ihi  - ilo  + 1
      jlen = jhi  - jlo  + 1
      klen = kbhi - kblo + 1
      llen = lbhi - lblo + 1
      bit16 = onbitmask(16)
      maxiirlen = 2*klen + (klen*llen)/2 + 1
      if (.not. ma_push_get(MT_INT, maxiirlen, 'i io pack ',l_iv, k_iv))
     $  call errquit('moints_aodisk_sprs2dense_a: no memory',0, MA_ERR)
      k_irv = k_iv + 2*klen
      labint = ma_sizeof(MT_INT,1,MT_BYTE)/2

      do i=ilo,ihi
        do j=jlo,jhi
          ijp = iand(ishft(ijmap(j,i),-16),bit16)
          nnz = iand(ijmap(j,i),bit16)
          if (nnz.gt.(klen*llen)) then
            call moints_aodisk_dumpiorec( reclen, iorec )
            call errquit('moints_sprs2dense_a: internal error 3',0,
     &       INT_ERR)
          endif
          if (ijp.gt.0) then
            nir = nnz/labint
            if (mod(nnz,labint).ne.0) nir = nir + 1
            iirlen = 2*klen + nir
            if (iirlen.gt.maxiirlen) then
              stop 1190
            endif
            dp  = ijp + ma_sizeof(MT_INT,iirlen,MT_DBL)
            call icopy( iirlen, iorec(ijp), 1, int_mb(k_iv), 1)
            call moints_sparse2d_unpack( llen, klen, int_mb(k_iv),
     $                                   nnz, int_mb(k_irv), iorec(dp),
     $                                   ssbb(lblo,kblo,j,i) )
          endif
        enddo
      enddo
      if (.not. ma_pop_stack(l_iv))
     $  call errquit('moints_aodisk_sprs2iorec: failed to pop', l_iv,
     &       MA_ERR)
      
      end

c
c  Define this macro assume everyone supports integer*2
c  Undefine if compiler complains
c
#define HAVE_FORTRAN_INTEGER2

c
c   This routine packs a dense 2d matrix (m x n)
c   into a standard sparse structure,
c   column pointer with non-zero rows, row-indices
c
c
c    COLPLO[]: 1 2 3 4 5 .... n
c              | |
c              | -----------
c              |            |
c              V            V
c    V[]:      1 ....       9 ...  
c    IR[]:     2 4 5 8...   2 4 9
c
c
c
#ifdef HAVE_FORTRAN_INTEGER2
      subroutine moints_sparse2d_pack( m, n, x, tol, labint, 
     $                                 colp, nir, ir, v, nnz )
      implicit none
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      double precision x(m,n)                ! [input]  dense 2D
      double precision tol                   ! [input]  tolerance
      integer labint                         ! [input]  2byte labels per integer word
      integer colp(2,n)                      ! [output] column ptr hi,lo
      integer nir                            ! [output] row index length
      integer*2 ir(*)                        ! [output] row index 16-bit packed
      double precision v(*)                  ! [output] packed values
      integer nnz                            ! [output] number non-zeroes
c      
      integer s, t
      double precision xx


      nnz = 0
      do s=1,n
        colp(1,s) = nnz + 1
        do t=1,m
          xx = x(t,s)
          if (abs(xx).ge.tol) then
            nnz  = nnz + 1
            v(nnz)  = xx
            ir(nnz) = t
          endif
        enddo
        if (colp(1,s).gt.nnz) colp(1,s) = 0
        colp(2,s) = nnz 
      enddo
      nir = nnz/labint
      if (mod(nnz,labint).ne.0) nir = nir + 1

      return
      end





      subroutine moints_sparse2d_unpack( m, n, colp, nnz, ir, v, x )
      implicit none
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      integer colp(2,n)                      ! [input]  column ptr hi,lo
      integer nnz                            ! [input]  number non-zeroes
      double precision v(nnz)                ! [input]  packed values
      double precision x(m,n)                ! [output] dense 2D
      integer*2 ir(*)                        ! [input]  row index 16-bit packed
      integer s, rp

      do s=1,n
        if (colp(1,s).gt.0) then
          do rp=colp(1,s),colp(2,s)
            x(ir(rp),s) = v(rp)
          enddo
        endif
      enddo
      return
      end



c
c   ----------------------------------------------
c    Alternative code without integer*2 data type 
c   ----------------------------------------------
c
#else
      subroutine moints_sparse2d_pack( m, n, x, tol, labint, colp, 
     $                                 nir, ir, v, nnz )
      implicit none
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      double precision x(m,n)                ! [input]  dense 2D
      double precision tol                   ! [input]  tolerance
      integer labint                         ! [input]  2byte labels per integer word
      integer colp(2,n)                      ! [output] column ptr hi,lo
      integer nir                            ! [output] row index length
      integer ir(*)                          ! [output] row index 16-bit packed
      double precision v(*)                  ! [output] packed values
      integer nnz                            ! [output] number non-zeroes
c      
      integer s, t, bp, pack
      integer ibit16
      double precision xx
#include "bitops_decls.fh"
#include "bitops_funcs.fh"      
      integer onbitmask
      external onbitmask

      ibit16 = onbitmask(16)
      nnz = 0
      nir = 0
      pack = 0
      bp = 0
      do s=1,n
        colp(1,s) = nnz + 1
        do t=1,m
          xx = x(t,s)
          if (abs(xx).gt.tol) then
            nnz  = nnz + 1
            v(nnz) = xx
            bp = mod((labint-mod(nnz,labint)),labint)*16
C            bp = mod((2-mod(nnz,2)),2)*16
            pack = ior(pack,ishft(t,bp))
            if (bp.eq.0) then
              nir = nir + 1
              ir(nir) = pack
              pack = 0
              if (nnz.ne.(labint*nir)) stop 3331
            endif
          endif
        enddo
        if (colp(1,s).gt.nnz) colp(1,s) = 0
        colp(2,s) = nnz 
      enddo
      if (bp.ne.0) then
        nir = nir + 1
        ir(nir) = pack
      endif
      return
      end





      subroutine moints_sparse2d_unpack( m, n, colp, nnz, ir, v, x )
      implicit none
#include "mafdecls.fh"
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      integer colp(2,n)                      ! [input]  column ptr hi,lo
      integer nnz                            ! [input]  number non-zeroes
      integer ir(*)                          ! [input]  row index 16-bit packed
      double precision v(nnz)                ! [input]  packed values
      double precision x(m,n)                ! [output] dense 2D
      integer s, t, rp, rrp, parity, nshft
      integer bit16, labint
#include "bitops_decls.fh"
#include "bitops_funcs.fh"      
      integer onbitmask
      external onbitmask

      labint = ma_sizeof(MT_INT,1,MT_BYTE)/2
      bit16 = onbitmask(16)
      do s=1,n
        if (colp(1,s).gt.0) then
          do rp=colp(1,s),colp(2,s)
            parity = mod(rp,labint)
            rrp = rp/labint
            if (parity.ne.0) rrp = rrp + 1
            nshft = 16*(labint - parity)
            if (parity.eq.0) nshft = 0
            t = iand(ishft(ir(rrp),-nshft),bit16)
            x(t,s) = v(rp)
          enddo
        endif
      enddo
      return
      end

#endif











      integer function dnonzero_cnt( n, tol, v )
      implicit none
      integer n
      double precision tol
      double precision v(n)
      integer i, nnz
      
      nnz = 0
      do i=1,n
        if (abs(v(i)).ge.tol) nnz = nnz + 1
      enddo
      dnonzero_cnt = nnz
      return
      end





c
c  =========================================================
c
c    Debugging and diagnostic codes
c
c  =========================================================
c

      
      
      
c
c  Debugging code to check file consistency
c
      subroutine moints_aodisk_checkfile( fd, reclen, iorec ) 
      implicit none
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
      integer fd
      integer reclen
      double precision iorec(reclen)
      integer stat, recbytes, rp, hdrp, ijp
      integer ilen, jlen
      integer dbilen
      double precision fptr
      integer blkinfo(MOINTS_NBLKINFO)
      integer magic_cookie, magic_cookie2
      integer ijmap(1000)
      integer onbitmask
      external onbitmask

      magic_cookie  = onbitmask(17)
      magic_cookie2 = onbitmask(21)
      dbilen = ma_sizeof(MT_INT,MOINTS_NBLKINFO,MT_DBL)
      stat     = 0
      recbytes = ma_sizeof(MT_DBL,reclen,MT_BYTE)
      fptr     = 1.d0
      rp       = 1
      hdrp     = 2
      blkinfo(2) = 0
      ilen = 0
      jlen = 0

 100  stat = eaf_read(fd, fptr, iorec, recbytes )
      if (stat.eq.0) then
        fptr = fptr + recbytes
        if ((blkinfo(2).eq.-1).and.(hdrp.ne.0)) then
          call icopy((ilen*jlen), iorec(2), 1, ijmap, 1 )
C          call moints_aodisk_ijprint( ilen, jlen, ijmap )
        endif
          
        call icopy(1, iorec(1), 1, hdrp, 1 )
        if (hdrp.le.0) goto 100
 101    call icopy(MOINTS_NBLKINFO, iorec(hdrp), 1, blkinfo, 1 )
        ijp = hdrp + dbilen
        if (blkinfo(1).eq.magic_cookie2) then
c$$$          write(6,912)
c$$$ 912      format('header invalid *')
          goto 100
        endif
        if (blkinfo(1).ne.magic_cookie) stop 6669


        ilen = blkinfo(8) - blkinfo(7) + 1
        jlen = blkinfo(10) - blkinfo(9) + 1
        call icopy((ilen*jlen), iorec(ijp), 1, ijmap, 1 )
c$$$        call moints_aodisk_ijprint( ilen, jlen, ijmap )

        if (blkinfo(2).gt.0) then
          hdrp = blkinfo(2)
          goto 101
        endif
        goto 100
      endif

      return
      end




      subroutine moints_aodisk_dumpiorec( n, r )
      implicit none
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
      integer n
      double precision r(n)

      integer hdrp, ijp, i, ilen, jlen
      integer bi(MOINTS_NBLKINFO)

      call icopy(1, r(1), 1, hdrp, 1)
      do while (hdrp.gt.0) 
        call icopy(MOINTS_NBLKINFO, r(hdrp), 1, bi, 1)
        write(6,901) (bi(i),i=1,MOINTS_NBLKINFO)
 901    format(16i8)
        ilen = bi(8)  - bi(7) + 1
        jlen = bi(10) - bi(9) + 1
        ijp = hdrp + ma_sizeof(MT_INT,MOINTS_NBLKINFO,MT_DBL)
        call moints_aodisk_ijprint( ilen, jlen, r(ijp) )

        hdrp = bi(2)
      enddo
      return
      end










      logical function moints_aodisk_verifyiorec(ilen, jlen, klen, 
     $                                           llen, n, r )
      implicit none
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
      integer ilen, jlen, klen, llen
      integer n
      double precision r(n)

      integer hdrp, cookie
      integer bi(MOINTS_NBLKINFO)
      integer itmp(1000)
      integer irec, ijp
      logical ao_ijverify
      external ao_ijverify
      integer onbitmask
      external onbitmask

      moints_aodisk_verifyiorec = .true.
      cookie = onbitmask(17)
      irec = 1
      ijp = 2
      call icopy(1, r(1), 1, hdrp, 1)
      do while ((irec.eq.1).or.(hdrp.gt.0))
        if ((hdrp.eq.2).or.(irec.gt.1)) then
          call icopy(MOINTS_NBLKINFO, r(hdrp), 1, bi, 1)
          if (bi(1).ne.cookie) return
          ijp = hdrp + ma_sizeof(MT_INT,MOINTS_NBLKINFO,MT_DBL)
        endif
        call icopy( (ilen*jlen), r(ijp), 1, itmp, 1)
        if (.not.(ao_ijverify( ilen, jlen, klen, llen, 
     $                         itmp, n, r))) then
          write(6,881) 
 881      format('Failed to verify')
          moints_aodisk_verifyiorec = .false.
          return
        endif
        if ((irec.gt.1).or.(hdrp.eq.2)) hdrp = bi(2)
        irec = irec + 1
      enddo
      return
      end






      subroutine moints_aodisk_ijprint( jlen, ilen, ijmap )
      implicit none
      integer ilen, jlen
      integer*2 ijmap(2,jlen,ilen)
      integer i, j

      do j=1,jlen
        write(6,771) (ijmap(1,j,i),i=1,ilen)
 771    format(16i5)
      enddo
      return
      end

          




      logical function ao_ijverify( ilen, jlen, klen, llen, 
     $                               ijmap, n, r )
      implicit none
      integer ilen, jlen, klen, llen
      integer*2 ijmap(2,jlen,ilen)
      integer n
      double precision r(n)
c
      integer i, j
      integer ijp, nnz, nir, iirlen
      integer itmp(1000)
      logical sprs_verify
      external sprs_verify

      ao_ijverify = .true.
      do i=1,ilen
	do j=1,jlen
          if (ijmap(1,j,i).ne.0) then
            ijp = ijmap(1,j,i)
            nnz = ijmap(2,j,i)
            nir = nnz/2 + mod(nnz,2)
            iirlen = 2*klen + nir
            call icopy( iirlen, r(ijp), 1, itmp, 1)
            if (.not.(sprs_verify(llen,klen,itmp(1),
     $                            itmp(2*klen+1)))) then
              print*,' --> CHECKING AT ',ijp
              write(6,771) i,j
 771          format('Failed to verify: I =',i5,'   J =',I5)
              write(6,881) 
 881          format(//,'IJ MAP dump')
              call moints_aodisk_ijprint(jlen,ilen,ijmap)
              ao_ijverify = .false.
              return
            endif
          endif
        enddo
      enddo
      return

      end







      subroutine sprs_print( m, n, colp, ir )
      implicit none
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      integer colp(2,n)                      ! [input]  column ptr hi,lo
      integer*2 ir(*)
      integer s, rp

      write(6,711) (colp(1,s),colp(2,s),s=1,n)
 711  format(//,16(i4,2x,i4,5x))
      do s=1,n
        if (colp(1,s).gt.0) then
          write(6,772) s,(ir(rp),rp=colp(1,s),colp(2,s))
 772      format(i5,5x,16i4)
        endif
      enddo
      return
      end


      logical function sprs_verify( m, n, colp, ir )
      implicit none
      integer m                              ! [input]  rows
      integer n                              ! [input]  columns
      integer colp(2,n)                      ! [input]  column ptr hi,lo
      integer*2 ir(*)
      integer s, t

      sprs_verify = .true.
      do s=1,n
        if (colp(1,s).gt.0) then
          do t=colp(1,s),colp(2,s)
            if ((ir(t).le.0).or.(ir(t).gt.m)) then
              write(6,811) s,ir(t)
 811          format('Failed to verify:  S=',i4,'  T=',i4)
              call sprs_print( m, n, colp, ir )
              sprs_verify = .false.
              return
            endif
          enddo
        endif
      enddo
      return
      end

c  
c End !OLD_AODISK section
c
#endif
