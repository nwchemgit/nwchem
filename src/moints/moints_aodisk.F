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
c  Use a wrapper around nxtval() and record the task number on
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
      data moao_fname/''/
      end





c
c  Wrapper around nxtval so we can save task numbers
c  with the IO records and emulate task scheduling
c  when retrieved
c
      integer function moints_nxttask( numnodes )
      implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "cmointsmoao.fh"
      integer numnodes
      
      if (moao_ipass.gt.0) then
        moints_nxttask = moao_tasknum
        if (numnodes.lt.0) call ga_sync()
      else
        moints_nxttask = nxtval( numnodes )
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
      character*8 buidstr
      integer recnum, ssbblen
      logical moints_iorec_next
      external moints_iorec_next
c$$$      DOUBLE PRECISION DABSSUM
c$$$      EXTERNAL DABSSUM

      moints_gblk_fromdisk = .false.
      if (moao_ipass.lt.0) return
      if ((moao_fd.ge.0).and.(moao_ipass.gt.0)) then
        call moints_uid( buidstr, ish, jsh, kshlo, lshlo )
        recnum = 0
        ssbblen = (ihi-ilo+1)*(jhi-jlo+1)*(kbhi-kblo+1)*(lbhi-lblo+1)
        call dfill( ssbblen, 0.d0, ssbb, 1 )
        do while (moints_iorec_next( blkid, buidstr  ))
          if (moao_issparse.eq.1) then
            call moints_aodisk_iorec2sprs( ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, moao_reclen, moao_lwidth,
     $                             dbl_mb(moao_klabrec),
     $                             dbl_mb(moao_kvalrec) )
            moints_gblk_fromdisk = .true.
          else
            call errquit('moints_gblk_disk: dense write not ready',0)
          endif
          recnum = recnum + 1
        enddo
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

c$$$          call moints_aodisk_sprs2iorec_N( moao_lun, moao_tasknum, 
c$$$     $                           blkid, moao_buflen,
c$$$     $                           ilo, ihi, jlo, jhi,
c$$$     $                           kblo, kbhi, lblo, lbhi,
c$$$     $                           ssbb, moao_buflen, dbl_mb(moao_kbuf) )

          call moints_aodisk_sprs2iorec( moao_fd, moao_fptr,
     $                           moao_tasknum, blkid,
     $                           ilo, ihi, jlo, jhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, moao_spreclen, moao_lwidth,
     $                           dbl_mb(moao_klabrec),
     $                           dbl_mb(moao_kvalrec),
     $                           moao_buflen, dbl_mb(moao_kbuf) )
        else
          call errquit('moints_gblk_todisk: dense read not ready',0)
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
c  Called by application
c  Enable disk caching & check for existing file and open file
c
      logical function moints_aodisk_init( odisk, oreuse )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      logical odisk                          ! [input] toggle caching
      logical oreuse                         ! [input] toggle reuse of existing file
c
      logical fexist
      integer stat
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      external moints_moao_block
#endif
c
c  first time through --- enable caching TO disk
c  otherwise must be caching FROM disk
c
      moints_aodisk_init = .false.
      if (.not.(odisk)) then
        moao_ipass = -1
        moints_aodisk_init = .true.
        return
      endif
      if (moao_ipass.lt.0) moao_ipass = 0
c
c  if AO file exists AND reuse enabled then set ipass > 0
c      
      if (moao_fname.eq.'') then
        call util_file_name( 'moao', .true., .true., moao_fname )
      endif
      inquire(file=moao_fname, exist=fexist)
      if (oreuse) then
        stat = 0
        if (fexist) stat=1
        call ga_igop( 481, stat, 1, '*')
        if (stat.eq.1) then
          moao_ipass = 1
          if (ga_nodeid().eq.0) write(6,331)
 331      format(10x,'Existing AO integral file will be reused')
        endif
      endif
c
c  open AO integral disk cache
c
      if (moao_fd.lt.0) then
c$$$        open(unit=MOINTS_AOUNIT, file=moao_fname, status='unknown',
c$$$     $       form='unformatted',iostat=stat)
        stat = eaf_open( moao_fname, EAF_RW, moao_fd )
        if (stat.ne.0) call errquit(
     $                'moints_aodisk_init: cannot open ao file',0)
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
#include "mafdecls.fh"
#include "eaf.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      integer stat
      integer blkinfo(MOINTS_NBLKINFO)
      integer vlsize, dhdrsize, vallen, bufbytes

      moints_aodisk_prep = .false.
      if (moao_ipass.lt.0) then
        moints_aodisk_prep = .true.
        return
      endif
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
     $       'moints_aodisk_prep: cannot allocate io buffer',0)
        moao_buflen   = MOINTS_IOBUFLEN                                     ! bufflen (double words)
        vallen        = moao_buflen - dhdrsize                              ! value len -= hdr length
        moao_spreclen = vallen/(vlsize/ma_sizeof(MT_DBL,1,MT_BYTE))         ! leng
        moao_kvalrec  = moao_kbuf + dhdrsize                                ! values after header
        moao_klabrec  = moao_kbuf + dhdrsize + moao_spreclen                ! labels after values
      endif
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
      endif
      moints_aodisk_prep = .true.
      MOAO_CUMUL = 0.d0
      MOAO_NREC = 0
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
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"

      if (moao_ipass.lt.0) return
      if (moao_lbuf.ge.0) then
        if (.not.ma_free_heap(moao_lbuf)) call errquit(
     $    'moints_closeaodisk: cannot free io buffer',0)
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
#include "mafdecls.fh"
#include "eaf.fh"
#include "util.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      logical osave
c
      integer stat

      moao_ipass = 0
      if (moao_fd.ge.0) then
        if (util_print('ao disk stats',print_high))
     $     call eaf_print_stats(moao_fd)
        stat = eaf_close(moao_fd)
        if (stat.ne.0) call errquit(
     $    'moints_closeaodisk: cannot close file',0)
        if (.not.(osave)) then
          stat = eaf_delete(moao_fname)
          if (stat.ne.0) call errquit(
     $    'moints_closeaodisk: cannot delete file',0)
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
c
c
c  Pack the dense SSBB block into sparse form
c  with 16 bits per label
c
c
      subroutine moints_aodisk_sprs2iorec( fd, fptr, tasknum, blkid,
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
          lab1 = ior(lshift(i,16),j)
          do k=kblo,kbhi
            do l=lblo,lbhi
              xx = ssbb(l,k,j,i)
              if (abs(xx).gt.1.d-12) then
                iir = iir + 1
                ioval(iir)   = xx
                lab2 = ior(lshift(k,16),l)
                if (lwidth.eq.1) then
                  iolab(1,iir) = ior(lshift(lab1,32),lab2)
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

      

      subroutine moints_aodisk_iorec2sprs( ilo, ihi, jlo, jhi,
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
          i = iand(rshift(iolab(1,q),48),i16mask)
          j = iand(rshift(iolab(1,q),32),i16mask)
          k = iand(rshift(iolab(1,q),16),i16mask)
          l = iand(iolab(1,q),i16mask)
          ssbb(l,k,j,i) = ioval(q)
        enddo
      else if (lwidth.eq.2) then
        do q=1,reclen
          i = iand(rshift(iolab(1,q),16),i16mask)
          j = iand(iolab(1,q),i16mask)
          k = iand(rshift(iolab(2,q),16),i16mask)
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
        call errquit('moints_io disk io error',0)
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












#ifdef NOCOMPILE
c
c  Alternate version of sparse packing
c  Labels/indices require integer of storage, 2*NCOL + NNZ/2,
c  cf. 2*NNZ for standard version.
c
c
      subroutine moints_aodisk_sprs2iorec_N( tasknum, 
     $                             blkid, iobuflen,
     $                             ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, reclen, iorec )
      implicit none
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
      integer tasknum
      integer blkid
      integer iobuflen
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer reclen
      double precision iorec(reclen)
c
      integer i, j, ilen, jlen, llen, klen
      integer swordlen, nir, hdrlen
      integer nnz
      integer rp, rleft, rlen
      integer issparse, recnum
      logical st
      integer iirlen, drlen
      integer l_iv, k_iv, k_irv
      integer l_dv, k_dv
c      
      logical moints_sparse2d_pack
      external moints_sparse2d_pack
c
      issparse = 1
      recnum = 0
      llen = lbhi - lblo + 1
      klen = kbhi - kblo + 1
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      swordlen = ma_sizeof(MT_DBL,1,MT_BYTE) + 2       ! num of bytes for double precision value + 2 byte label = 10

      hdrlen = 3*ilen*jlen + MOINTS_NBLKINFO
      hdrlen = hdrlen + mod(hdrlen,ma_sizeof(MT_DBL,1,MT_INT))
      rp = hdrlen + mod(hdrlen,2)
      rlen = reclen - hdrlen
c
c  Allocate temp space
c
      iirlen = 2*klen + (klen*llen)/2 + 1
      drlen = klen*llen
      if (.not. ma_push_get(MT_INT, iirlen, 'i io pack ',l_iv, k_iv))
     $  call errquit('moints_aodisk_sprs2iorec: no memory',0)
      if (.not. ma_push_get(MT_DBL, drlen, 'd io pack ',l_dv, k_dv))
     $  call errquit('moints_aodisk_sprs2iorec: no memory',0)
      k_irv = k_iv + 2*klen
c
c
c
      do i=ilo,ihi
        do j=jlo,jhi

C          CALL MOINTS_MATPRINTX( LLEN, KLEN, SSBB(LBLO, KBLO, J, I))

          rleft = rlen - rp + 1

          st = moints_sparse2d_pack( klen, llen, ssbb(lblo,kblo,j,i),
     $                               int_mb(k_iv), nir, int_mb(k_irv), 
     $                               dbl_mb(k_dv), nnz )

          call icopy( (2*klen+nir), int_mb(k_iv), 1, iorec(rp), 1 )
          rp = rp + ma_sizeof(MT_INT,(2*klen+nir),MT_DBL)
          call dcopy( nnz, dbl_mb(k_dv), 1, iorec(rp), 1 )
          rp = rp + nnz

          call print_sparsemat( klen, llen, int_mb(k_iv),
     $                          int_mb(k_irv), nnz, dbl_mb(k_dv) )

        enddo
      enddo
c      
      if (.not. ma_pop_stack(l_dv))
     $  call errquit('moints_aodisk_sprs2iorec: failed to pop', l_dv)
      if (.not. ma_pop_stack(l_iv))
     $  call errquit('moints_aodisk_sprs2iorec: failed to pop', l_iv)
      return
      end



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
#ifdef INTEGER64

#else
      logical function moints_sparse2d_pack( n, m, x, colp, nir, 
     $                                       ir, v, nnz )
      implicit none
      integer n                              ! [input]  columns
      integer m                              ! [input]  rows
      double precision x(m,n)                ! [input]  dense 2D
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

      ibit16 = 2**16 - 1
      moints_sparse2d_pack = .false.
      nnz = 0
      nir = 0
      pack = 0
      do s=1,n
        colp(1,s) = nnz + 1
        do t=1,m
          xx = x(t,s)
          if (abs(xx).gt.1.d-12) then
            nnz  = nnz + 1
            v(nnz) = xx
            bp = mod((2-mod(nnz,2)),2)*16
            pack = ior(pack,lshift(t,bp))
            if (bp.eq.0) then
              nir = nir + 1
              ir(nir) = pack
              pack = 0
              if (nnz.ne.(2*nir)) stop 3331
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
      moints_sparse2d_pack = .true.
      return
      end

#endif

#endif 
c  NOCOMPILE





