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
c       moao_lun    < 0      io not initialized
c                   > 0      io is happening
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
      data moao_lun/-1/
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
c  Return complete SSBB block of integrals
c  if cached on disk
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
      integer recnum, ssbblen
      logical moints_iorec_next
      external moints_iorec_next

      moints_gblk_fromdisk = .false.
      if (moao_ipass.lt.0) return
      if ((moao_lun.gt.0).and.(moao_ipass.gt.0)) then
        recnum = 0
        ssbblen = (ihi-ilo+1)*(jhi-jlo+1)*(kbhi-kblo+1)*(lbhi-lblo+1)
        call dfill( ssbblen, 0.d0, ssbb, 1 )
        do while (moints_iorec_next( blkid  ))
          if (moao_issparse.eq.1) then
            call moints_aodisk_iorec2sprs( ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, moao_reclen, 
     $                             dbl_mb(moao_klabrec),
     $                             dbl_mb(moao_kvalrec) )
            moints_gblk_fromdisk = .true.
          else
            call errquit('moints_gblk_disk: dense read not ready',0)
          endif
          recnum = recnum + 1
        enddo
C        MOAO_CUMUL = MOAO_CUMUL + DABSSUM( SSBBLEN, SSBB )
      endif
      return
      end




c
c  Save block of SSBB block to disk if caching enabled.
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

      moints_gblk_todisk = .false.
      if (moao_ipass.lt.0) return
      if ((moao_lun.gt.0).and.(moao_ipass.eq.0)) then
        if (sparse) then
          call moints_aodisk_sprs2iorec( moao_lun, moao_tasknum, 
     $                           blkid, moao_buflen,
     $                           ilo, ihi, jlo, jhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, moao_spreclen, 
     $                           dbl_mb(moao_klabrec), 
     $                           dbl_mb(moao_kvalrec) )
        else
          call errquit('moints_gblk_todisk: dense read not ready',0)
        endif
       moints_gblk_todisk = .true.
c
C        MOAO_CUMUL = MOAO_CUMUL + DABSSUM( SSBBLEN, SSBB )
c
      endif
      return
      end






c
c  This must be called by application
c  Enable disk caching
c  Check for existing file and open file
c
      logical function moints_aodisk_init( odisk, oreuse )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
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
      if (moao_lun.lt.0) then
        open(unit=MOINTS_AOUNIT, file=moao_fname, status='unknown',
     $       form='unformatted',iostat=stat)
        if (stat.ne.0) call errquit(
     $                'moints_aodisk_init: cannot open ao file',0)
        moao_lun = MOINTS_AOUNIT
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
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      integer stat
      integer blkinfo(MOINTS_NBLKINFO)

      moints_aodisk_prep = .false.
      if (moao_ipass.lt.0) then
        moints_aodisk_prep = .true.
        return
      endif
c
c  allocate io buffers
c
      if (moao_lbuf.lt.0) then
        if (.not.ma_alloc_get(MT_DBL, MOINTS_IOBUFLEN, 'moints iobuf',
     $                        moao_lbuf, moao_kbuf)) call errquit(
     $       'moints_aodisk_prep: cannot allocate io buffer',0)
        moao_buflen   = MOINTS_IOBUFLEN                                     ! bufflen (double words)
        moao_spreclen = moao_buflen/(ma_sizeof(MT_INT,4,MT_DBL)+1)
        moao_kvalrec  = moao_kbuf                                           ! values
        moao_klabrec  = moao_kbuf + moao_spreclen                           ! labels after values
      endif
c
c  rewind AO integral cache unit
c
      rewind(moao_lun, iostat=stat)
      moao_eof = .false.
      if (stat.ne.0) call errquit(
     $  'moints_aodisk_prep: failed to rewind unit',0)
c
c  if READ mode -- retrieve first record
c
      if (moao_ipass.gt.0) then
        call moints_iorec_op( 'r', moao_lun, blkinfo, moao_buflen,
     $                     dbl_mb(moao_kbuf), stat )
        if (stat.eq.0) then
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
c  Keep NO state information except moao_ipass and moao_lun
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
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      logical osave
c
      integer stat

      moao_ipass = 0
      if (moao_lun.gt.0) then
        if (osave) then
          close(unit=moao_lun, iostat=stat, status='keep')
        else
          close(unit=moao_lun, iostat=stat, status='delete')
        endif
        if (stat.ne.0) call errquit(
     $    'moints_closeaodisk: cannot close unit',0)
        moao_lun = -1
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
      subroutine moints_aodisk_sprs2iorec( lunao, tasknum, 
     $                             blkid, iobuflen,
     $                             ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, reclen, 
     $                             iolab, ioval )
      implicit none
      integer lunao
      integer tasknum
      integer blkid
      integer iobuflen
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer reclen
      integer iolab(4,reclen)
      double precision ioval(reclen)

      integer i, j, k, l
      integer iir, issparse
      integer recnum
      double precision xx

      issparse = 1
      recnum = 0
      iir = 0
      call ifill( 4*reclen, 0, iolab, 1 )
      call dfill( reclen, 0.d0, ioval, 1 )
      do i=ilo,ihi
        do j=jlo,jhi
          do k=kblo,kbhi
            do l=lblo,lbhi
              xx = ssbb(l,k,j,i)
              if (abs(xx).gt.1.d-12) then
                iir = iir + 1
                ioval(iir)   = xx
                iolab(1,iir) = i
                iolab(2,iir) = j
                iolab(3,iir) = k
                iolab(4,iir) = l
                if (iir.eq.reclen) then
                  call moints_iorec_flush(lunao, tasknum, iir, issparse,
     $                                   blkid, iobuflen, ioval )
                  recnum = recnum + 1
                  call ifill( 4*reclen, 0, iolab, 1 )
                  call dfill( reclen, 0.d0, ioval, 1 )
                  iir = 0
                endif
              endif
            enddo  
          enddo  
        enddo  
      enddo
      if (iir.gt.0) then
        call moints_iorec_flush( lunao, tasknum, iir, issparse,
     $                          blkid, iobuflen, ioval )
        recnum = recnum + 1
        iir = 0
      endif
                
      return
      end

      


    
      subroutine moints_aodisk_iorec2sprs( ilo, ihi, jlo, jhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb, reclen, iolab,
     $                             ioval )
      implicit none
      integer ilo, ihi, jlo, jhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlo:jhi,ilo:ihi)
      integer reclen
      integer iolab(4,reclen)
      double precision ioval(reclen)
      integer q, i, j, k, l
      
      do q=1,reclen
        i = iolab(1,q)
        j = iolab(2,q)
        k = iolab(3,q)
        l = iolab(4,q)
        ssbb(l,k,j,i) = ioval(q)
      enddo
      return
      end




      subroutine moints_iorec_flush( lunao, tasknum, reccnt, issparse, 
     $                               blkid, buflen, iobuf )
      implicit none
#include "mointsmoaodef.fh"
      integer lunao
      integer tasknum
      integer reccnt
      integer issparse
      integer blkid
      integer buflen
      double precision iobuf(buflen)
c
      integer blkinfo(MOINTS_NBLKINFO)
      integer stat
c
      call ifill( MOINTS_NBLKINFO, 0, blkinfo, 1 )
      blkinfo(1) = tasknum
      blkinfo(2) = reccnt
      blkinfo(3) = issparse
      blkinfo(4) = blkid
c$$$      WRITE(6,772) blkid
c$$$ 772  FORMAT('---->',I8)
c
      call moints_iorec_op( 'w', lunao, blkinfo, buflen, iobuf, stat )

      return
      end

      
      



      logical function moints_iorec_next( blkid )
      implicit none
#include "mafdecls.fh"
#include "mointsmoaodef.fh"
#include "cmointsmoao.fh"
      integer blkid
      integer blkinfo(MOINTS_NBLKINFO)
      integer stat

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
      call moints_iorec_op( 'r', moao_lun, blkinfo, moao_buflen, 
     $                  dbl_mb(moao_kbuf), stat )
      MOAO_NREC = MOAO_NREC + 1
c
c  check for EOF and/or update info
c
      if (stat.eq.0) then
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
        
      
      
                         




      subroutine moints_iorec_op( op, lun, blkinfo, buflen, buf, stat )
      implicit none
#include "mointsmoaodef.fh"
      character*1 op
      integer lun
      integer blkinfo(MOINTS_NBLKINFO)
      integer buflen
      double precision buf(buflen)
      integer stat
      
      if ((op.eq.'w').or.(op.eq.'W')) then
        write(lun,iostat=stat) blkinfo, buf
      else
        read(lun,iostat=stat) blkinfo, buf
      endif
      return
      end

      

