      block data cgridfdata
*
* $Id: grid_file.F,v 1.7 1999-06-15 20:16:25 edo Exp $
*
      implicit none
#include "cgridfile.fh"
      data file_opened /.false./
      data l_cache /-1/
      end
      logical function grid_file_open(
     $     filename,
     $     cachesize, filesize)
      implicit none
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "eaf.fh"
      character*(*) filename
      integer cachesize, filesize
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      external cgridfdata
#endif
c
      integer availkb, sdlen
      character*20 fstype
      character*(nw_max_path_len) scratchdir
      integer util_scratch_dir_avail_for_me
      external util_scratch_dir_avail_for_me
c
c     Open file returning CHEMIO filehandle (not yet)
c     Return true on success, false on failure
c
c     Currently only allows for restart if entire file is disk resident
c     ... needs mods to allow writing of memory cached grid_pts
c     (actually gets simpler once the file routines do the caching)
c
      if (filename.eq.' ') call errquit('grid_file_open: empty name',0)
      call util_pname(filename, actualname)
c
      otrace = util_print('griddebug',print_debug)
      oinfo  = util_print('gridinfo', print_default) .and.
     $     ga_nodeid().eq.0
c
      if (otrace) then
         write(6,*) ' grid_file: opening ',
     $        actualname(1:inp_strlen(actualname))
         call util_flush(6)
      end if
c
c     Initialize common
c
      file_opened = .false.
      cur_rec = 1
      n_rec_in_file = 0
      rec_len_in_dbl = 4*buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
      call grid_buf_clear()
      call dfill(4*buffer_size, 0.0d0, values, 1)
      cur_rec_rad=0d0
#ifdef NOAIO
      oasyncio = .false.
#else
      oasyncio = .true.
#endif
c
c     Figure out sizes.  The -1 is to leave room for the end record
c     which for semidirect may be written out after EOF has been reported.
c
      if (cachesize .gt. 0) then
         max_cache_rec = cachesize/rec_len_in_dbl 
      else
         max_cache_rec = 0
      end if
      if (filesize .ge. 0) then
c
c     Filesize=0 means use default. 
c     Filesize>0 means user has specified value but need to use
c     min of this value and actual available space.
c
c     If the file is in the scratch directory then use the piece
c     of it available for this process, otherwise use all that
c     is available from eaf_stat.
c
         call util_directory_name(scratchdir, .true., ga_nodeid())
         sdlen = inp_strlen(scratchdir)
         if (scratchdir(1:sdlen) .eq. actualname(1:sdlen)) then
            availkb = util_scratch_dir_avail_for_me()
         else
            if (eaf_stat(actualname, availkb, fstype) .ne. 0)
     $           call errquit
     $           ('grid_file_open:unable to stat available disk',0)
            write(6,*) ' non-SCRATCH ', availkb
         endif
c
         max_file_rec = nint((0.95*1024.0d0*availkb) / 
     $        (8.0d0*rec_len_in_dbl))
         max_file_rec = max(0,max_file_rec)
c
         if (filesize.gt.0)
     $        max_file_rec = min(max_file_rec,filesize/rec_len_in_dbl)
      else 
c
c     Filesize < 0 means user forced no filespace at all.
c
         max_file_rec = 0
      end if
      if (max_file_rec .gt. 0) then ! For the semidirect logic
         max_file_rec = max_file_rec - 1
      else if (max_cache_rec .gt. 0) then
         max_cache_rec = max_cache_rec - 1
      endif
c
      if (max_file_rec.eq.0 .and. max_cache_rec.eq.0) then
         grid_file_open = .false.
         return
      endif
c
c     defer allocation of memory until the first write otherwise
c     this memory is sitting around during the initial guess routines
c
      l_cache = -1
      k_cache = -1
c
      grid_file_open = .true.
c
      end
      logical function grid_file_close(okeep)
      implicit none
#include "cgridfile.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "util.fh"
#include "global.fh"
      logical okeep  ! [input] If true keep, otherwise delete
c
      if (file_opened) then
         if (oinfo .or. (util_print('iostats',print_high)
     $           .and.ga_nodeid().eq.0)) 
     $        call eaf_print_stats(fd)
         if (eaf_close(fd) .ne. 0) call errquit
     $        ('grid_file_close: eaf_close failed',0)
         if (.not. okeep) then
            if (eaf_delete(actualname) .ne. 0) call errquit
     $           ('grid_file_close: delete failed', 0)
         endif
      endif
c
      if (l_cache .ne. -1) then
         if (.not. ma_free_heap(l_cache))
     $        call errquit('grid_file_close: free heap?',0)
      end if
c
      if (file_opened .or. l_cache.ne.-1) then
         call ga_igop(msg_grid_stat1, n_rec_in_file, 1, '+')
         if (oinfo) then
            write(6,1) n_rec_in_file
 1          format(/
     $           ' Parallel grid_pts file used',i8,' records'/)
            call util_flush(6)
         end if
         if (otrace) then
            write(6,*) ' grid_file: closing ',
     $           actualname(1:inp_strlen(actualname))
            call util_flush(6)
         end if
      endif
c
      l_cache = -1
      file_opened = .false.
      cur_rec = 0
      n_rec_in_file = 0
      call grid_buf_clear
      max_cache_rec = 0
c     
      grid_file_close = .true.
c
      end
      logical function grid_file_rewind()
      implicit none
#include "cgridfile.fh"
#include "inp.fh"
c
c     rewind the file ready to read grid_pts ... hide any header.
c
      cur_rec = 1
      call grid_buf_clear
c
      if (otrace) then
         write(6,*) ' grid_file: rewinding ', 
     $        actualname(1:inp_strlen(actualname))
         call util_flush(6)
      end if
c
      grid_file_rewind = .true.
c
      end
      logical function grid_buf_write(oend)
      implicit none
#include "cgridfile.fh"
      logical oend
c
      logical grid_packed_buf_write
      external grid_packed_buf_write
c
      if (n_in_rec .gt. 1) then      ! First entry is just the counter
        call grid_buf_pack()
        grid_buf_write = grid_packed_buf_write(oend)
      else
        grid_buf_write = .true.
        call grid_buf_clear()
      endif
c
      end
      logical function grid_packed_buf_write(oend)
      implicit none
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "eaf.fh"
      logical oend           ! [input]
c
c     Write the grid_pts buffer to disk.  For all but the
c     end record it is called with oend=.false. meaning
c     that size constraints should be enforced leaving
c     room for the end record.  The last record is output with
c     oend=.true. meaning just write it.
c
      integer nbitv             ! For pretty output
      integer maxrec
      double precision offset
      integer ierr
      character*80 errmsg
c
c     test for overflowing the file and memory cache
c
      maxrec = max_cache_rec + max_file_rec
      if ((.not. oend) .and. (cur_rec .gt. maxrec)) then
         grid_packed_buf_write = .false.
         if (otrace) then
            write(6,12) cur_rec
 12         format(' grid_file: file full at record ',i8)
            call util_flush(6)
         end if
         return
      endif
c
c     Opening of file, allocation of buffer and writing of
c     info message are all now deferred until the first
c     write operation
c
      if (cur_rec .eq. 1) then
         if (oinfo) then
            nbitv = 64
            write(6,20) actualname(1:inp_strlen(actualname)),
     $           rec_len_in_dbl,n_per_rec,max_cache_rec,max_file_rec
 20         format(/
     $           ' Grid_pts file          = ', a/
     $           ' Record size in doubles = ', i6,7x,
     $           ' No. of grid_pts per rec  = ', i6/
     $           ' Max. records in memory = ', i6,7x,
     $           ' Max. records in file   = ', i6/)
            call util_flush(6)
         end if
      end if
c
 10   continue
      if (cur_rec .le. max_cache_rec) then
         if (l_cache .eq. -1) then
            if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $           'grid', l_cache, k_cache)) then
               write(6,*) ' grid_file: Failed to allocate cache '//
     $              '... halving'
               max_cache_rec = max_cache_rec/2
               goto 10          ! Dirty ... go back and try again
            end if
         end if
         call dcopy(rec_len_in_dbl, buf, 1, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1)
      else
         if (.not. file_opened) then
            if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     $           call errquit('grid_file_open: eaf_open failed', 0)
            file_opened = .true.
         end if
c
c     Should really detect failures in writing and return false
c
         offset = 8.0d0*rec_len_in_dbl*(cur_rec-max_cache_rec-1)
         ierr = eaf_write(fd, offset, buf, 8*rec_len_in_dbl)
         if (ierr .ne. 0) then
            call eaf_errmsg(ierr, errmsg)
            write(6,*) ' IO offset ', offset
            write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
            call errquit('grid_packed_buf_write: write failed',0)
         endif
      end if
c
      if (otrace) then
         write(6,2) n_in_rec, cur_rec, n_rec_in_file+1
 2       format(' grid_file: wrote ',i6,' grid_pts in record',i5
     $        ,' of ',i5,'.')
c         write(6,*)' buf: '
c         call output(buf, 1, 4, 1, n_in_rec, 4, rec_len_in_dbl, 1)
         call util_flush(6)
      end if
c
      call grid_buf_clear
      n_rec_in_file   = cur_rec ! was n_rec_in_file + 1
      cur_rec         = cur_rec + 1
      grid_packed_buf_write = .true.
c
      end
      logical function grid_buf_read()
      implicit none
#include "cgridfile.fh"
      logical grid_packed_buf_read
      external grid_packed_buf_read
c
      grid_buf_read = grid_packed_buf_read()
      if (.not. grid_buf_read) return
c
      call grid_buf_unpack()
c
*      if (n_in_rec .gt. 1) then
*         grid_buf_read = .true.
*      else
*         grid_buf_read = .false.
*      end if
c     
      if (otrace) then
            write(6,2) n_in_rec, cur_rec-1, n_rec_in_file
 2          format(' grid_file: read ',i6,' grid_pts in record',i6
     $           ,' of ',i6,'.')
c         write(6,*)' buf: '
c         call output(buf, 1, 4, 1, n_in_rec, 4, rec_len_in_dbl, 1)
         call util_flush(6)
      end if
c
      end
      logical function grid_packed_buf_read()
      implicit none
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "global.fh"
c
      double precision offset
      integer ierr
      character*80 errmsg
c
      call grid_buf_clear
      if (cur_rec .gt. n_rec_in_file) then
         grid_packed_buf_read = .false.
         if (otrace) then
            write(6,1) cur_rec, n_rec_in_file
 1          format(' grid_file: EOF at record',i5,' of ',i5,'?')
            call util_flush(6)
         end if
         return
      end if
c
      grid_packed_buf_read = .true.
      if (cur_rec .le. max_cache_rec) then
c
c     Record is cached in memory
c
         call dcopy(rec_len_in_dbl, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1,
     $        buf, 1)
      else
c
c     Record is on disk
c
         offset = 8.0d0*rec_len_in_dbl*(cur_rec-max_cache_rec-1)
c
         if ((.not. oasyncio) .or. (.not. oiopending)) then
c
c     Either synchronous IO, or the record has not been read already
c
            ierr = eaf_read(fd, offset, buf, 8*rec_len_in_dbl)
            if (ierr.ne.0) then
               call eaf_errmsg(ierr,errmsg)
               write(6,*) ga_nodeid(),errmsg
               call errquit('grid_buf_read:premature end of file',
     $              cur_rec)
            endif
         else
c
c     The record has been read asynchronously ... copy it over
c
            ierr = eaf_wait(fd, aioreq)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr,errmsg)
               write(6,*) ga_nodeid(),errmsg
               call errquit('grid_buf_read:eaf_wait failed',cur_rec)
            endif
            oiopending = .false.
            call dcopy(rec_len_in_dbl, buf2, 1, buf, 1)
         endif
c         
         if (oasyncio .and. ((cur_rec+1).le.n_rec_in_file)) then
c
c     Post a read for the next record
c
            offset = 8.0d0*rec_len_in_dbl*(cur_rec-max_cache_rec)
            ierr = eaf_aread(fd, offset, buf2, 8*rec_len_in_dbl, aioreq)
            oiopending = .true.
            if (ierr.ne.0) then
               call eaf_errmsg(ierr,errmsg)
               write(6,*) ga_nodeid(),errmsg
               call errquit('grid_buf_read:asynch IO failed',cur_rec+1)
            endif
         endif
      end if
c
      cur_rec = cur_rec + 1
c
      end
      subroutine grid_buf_clear()
      implicit none
#include "cgridfile.fh"
c
c     empty the grid_pts buffer ... only the atom center and count of grid_pts
c     in the buffer is present.
c
      n_in_rec = 1
      next_value = 2
c
      end
      subroutine grid_buf_unpack()
      implicit none
#include "cgridfile.fh"
#include "mafdecls.fh"
      integer ncopy
c
c
c     Copy count
c
      n_in_rec = n_in_buf
c     
c     Copy values
c
c     Total dbls to copy is 4*(n_per_rec + 1)
c
      ncopy = 4*(n_per_rec + 1)
      call dcopy(ncopy, buf, 1, values, 1)
      cur_rec_ctr = nint(values(1,1))
      n_in_rec = nint(values(2,1))
      cur_rec_rad = values(3,1)
c
      end
      subroutine grid_buf_pack()
      implicit none
#include "cgridfile.fh"
#include "mafdecls.fh"
      integer ncopy
c
c     Copy count
c
      n_in_buf = n_in_rec
c
      values(1,1) = dble(cur_rec_ctr)
      values(2,1) = dble(n_in_rec)
      values(3,1) = cur_rec_rad
c
c     Total dbls to copy is 4*(n_per_rec + 1)
c
      ncopy = 4*(n_per_rec + 1)
      call dcopy(ncopy, values, 1, buf, 1)
c
      end
      subroutine grid_packing_info
      implicit none
#include "cgridfile.fh"
c
c     determine n_per_rec; each grid pt requires x,y,z and weight
c     NOTE: use the 1st two places to write atom center and nqpts for this buffer
c
      n_per_rec = buffer_size - 1
      n_per_rec = n_per_rec - mod(n_per_rec,2) ! Make it even
c
      end
      logical function grid_test_mem(rtdb, geom, nq_tot)
      implicit none
#include "util.fh"
#include "bas.fh"
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "sym.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input]
      integer geom              ! [input]
      integer nq_tot  ! [input]
c      integer basis             ! [input]
c     
c     By default memsize=filesize=0 implying if the calculation fits into
c     memory, then use memory, otherwise use all available disk space.
c
c     If either size is negative then that medium is not to be used.
c     If either is size is positive then it has been user set.
c     
      logical grid_file_open
      external grid_file_open
      integer avail, localmemsize, nops
      integer filesize, memsize
      character*255 gridfilename
      double precision nq
      logical ofirst
      data ofirst /.true./      ! If first call to this routine
c     
      grid_test_mem = .false.
c
      if (file_opened) return
      if (.not. rtdb_get(rtdb, 'grid:filesize',
     $     MT_INT, 1, filesize)) filesize = 0
      if (.not. rtdb_get(rtdb, 'grid:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
      if (filesize.lt.0 .and. memsize.lt.0) return
c
c     Determine if we should automatically run an in-core calculation
c
      oinfo  = util_print('gridinfo', print_default) .and.
     $     ga_nodeid().eq.0
      nops = sym_number_ops(geom) + 1
c
      rec_len_in_dbl = 4*buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
c
c     How much memory?  Leave space for GAs, for the grid_pts code,
c     and for 50,000 of misc. junk. 
c
      if (memsize .eq. 0) then
         avail = MA_inquire_avail(mt_dbl)
*      write(6,*) ' avail 1 ', avail
         if (ga_uses_ma() .and. ga_memory_limited()) 
     $        avail = avail - ga_memory_avail()
*      write(6,*) ' avail 2 ', avail
c         call intb_mem_2e4c(max2e, mem2)
*      write(6,*) ' avail 3 ', avail
c         avail = avail - 5*max(max2e,maxeri) - mem2
*      write(6,*) ' avail 4 ', avail
         avail = avail - 100000
*      write(6,*) ' avail 5 ', avail
         call ga_igop(msg_grid_file, avail, 1, 'min')
*      write(6,*) ' avail 6 ', avail
c     
c     Assume that for the small systems for which can hold the
c     grid_pts, that there is no sparsity except for symmetry.
c     Also, assume that the grid_pts will be evenly distributed
c     apart from a 10% fluctuation
c     
c         nq = 0.125d0 * dble(nbf+1)**4 
         nq=nq_tot
*      write(6,*) ' nq 1 ', nq
         if (nops .ne. 1) nq = (1.3 * nq) / nops ! 1.3 empircal from C6H6
*      write(6,*) ' nq 1 ', nq
         nq = nq / ga_nnodes()
*      write(6,*) ' nq 1 ', nq
         nq = nq * 1.1d0 + 2*rec_len_in_dbl ! Fudge factors for safety
*      write(6,*) ' nq 1 ', nq
c
         if (nq .lt. avail) then
            localmemsize = max(8*buffer_size,int(nq))
         else
            localmemsize = 0
         endif
      else
         localmemsize = 0
      endif
c
      if (.not. rtdb_cget(rtdb, 'grid:filename', 1, gridfilename)) 
     $        call errquit('grid_file: opening gridpts?', 0)
      if (localmemsize.gt.0 .or. filesize.ge.0) then
         if (.not. grid_file_open(gridfilename, 
     $        localmemsize, filesize))
     $        call errquit('grid_test_mem: grid_file_open failed',0)
         grid_test_mem = .true.
      endif
*      if (localmemsize.eq.0 .and. oinfo .and. ofirst .and. memsize.eq.0) then
*         write(6,1) nq-avail
* 1       format(/' Provide',1p,d9.2,' more words/process to ',
*     $        'cache all grid_pts in memory'/)
*      end if
c
      ofirst = .false.
c
      end
      logical function grid_file_write(nqpts, nctr, rad, aqbuf)
c      logical function grid_file_write(nints, i, j, k, l, g)
      implicit none
#include "cgridfile.fh"
c
c     insert grid_pts into internal buffer writing to disk if necessary.
c     return true on success. false if insufficient memory/file space.
c     hard fail for any other errors.
c
      integer nqpts
      integer nctr
      double precision aqbuf(4,nqpts)
      double precision rad
      logical grid_buf_write
      external grid_buf_write
c
      integer ind
c
      grid_file_write = .true.
c
      if (otrace) then
         write(6,*) ' grid_file write: nqpts,nctr,n_in_rec,next_value', 
     &                                 nqpts,nctr,n_in_rec,next_value
         call util_flush(6)
      end if
c
      n_in_rec      = next_value - 1
      n_rec_in_file = cur_rec - 1
      if(n_in_rec .eq. 1) cur_rec_ctr = nctr
      cur_rec_rad = max(rad,cur_rec_rad)
c
c     Check to see if buffer needs dumping before adding more pts.
c
      if ((nctr .eq. 0) .or. 
     &    (nctr. ne. cur_rec_ctr)) then
         grid_file_write = grid_buf_write(.false.)
         cur_rec_ctr = nctr
         if (nctr.eq.0)nctr =  n_rec_in_file
         if (.not. grid_file_write) return
         cur_rec_rad=0d0
      end if
c
c     This copy can be made a lot more efficient
c     
      do ind = 1, nqpts
         if (n_in_rec .eq. n_per_rec) then
            grid_file_write = grid_buf_write(.false.)
            if (.not. grid_file_write) return
         end if
c     
         n_in_rec   = next_value
         next_value = next_value + 1
c     
         values(1,n_in_rec)   = aqbuf(1,ind)
         values(2,n_in_rec)   = aqbuf(2,ind)
         values(3,n_in_rec)   = aqbuf(3,ind)
         values(4,n_in_rec)   = aqbuf(4,ind)
c
      end do
c
      end
      logical function grid_file_read(maxqpts, nqpts, nctr,rad,aqbuf)
      implicit none
c
#include "cgridfile.fh"
c
      integer maxqpts
      integer nqpts
      integer nctr
      double precision aqbuf(4,maxqpts),rad
c
c     Return grid pts corresponding to the same atom center.
c     Return in nqpts the number of grid pts read.  
c
      integer ncopy
      logical grid_buf_read
      external grid_buf_read
c
      grid_file_read = .true.
      nqpts = 0
c
      grid_file_read = grid_buf_read()
      if (.not. grid_file_read) return
c
c     Total dbls to copy is 4*(n_per_rec)
c
      ncopy = 4*(n_per_rec)
      call dcopy(ncopy, values(1,2), 1, aqbuf(1,1), 1)
      nctr = cur_rec_ctr
      nqpts = n_in_rec - 1
      rad =cur_rec_rad
c
      if (otrace) then
         write(6,*) ' grid_file read: nqpts,nctr,n_in_rec,next_value', 
     &                                nqpts,nctr,n_in_rec,next_value
         call util_flush(6)
      end if
c
      return
c
      end
      integer function grid_num_rec()
      implicit none
#include "cgridfile.fh"
      grid_num_rec=n_rec_in_file
      return
      end
