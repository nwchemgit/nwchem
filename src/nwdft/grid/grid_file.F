      Block data cgridfdata
*
* $Id$
*
      implicit none
#include "cgridfile.fh"
      data file_opened /.false./
      data filedump_opened /.false./
      data oreopen /.false./
      data l_cache /-1/
      data fdg /-1/
      data fd /-1/
      end
      logical function grid_file_open(enuc_in,
     $     filename,
     $     cachesize, filesize, eaf_size_in_dbl)
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "zora.fh"
      character*(*) filename
      integer cachesize, filesize
      double precision enuc_in ! [in]
      double precision enuc
      integer eaf_size_in_dbl
      integer  sdlen
      double precision availkb
      character*20 fstype
      character*(nw_max_path_len) scratchdir
      integer ierr
      integer me
      double precision util_scratch_dir_avail_for_me
      character*80 errmsg
      character*256 grinfo_in
      external util_scratch_dir_avail_for_me
      logical grid_file_rewind,grid_file_close
      logical oreafail
      external grid_file_rewind,grid_file_close
      logical grid_grinfoinit,grid_grinfo2,grid_ptsalloc
      external grid_grinfoinit,grid_grinfo2,grid_ptsalloc
#ifdef EAFHACK
      logical grid_ptsdel
      external grid_ptsdel
#endif
c     Open file returning CHEMIO filehandle (not yet)
c     Return true on success, false on failure
c
c     Currently only allows for restart if entire file is disk resident
c     ... needs mods to allow writing of memory cached grid_pts
c     (actually gets simpler once the file routines do the caching)
c
      if (filename.eq.' ') call errquit('grid_file_open: empty name',0,
     &       DISK_ERR)
      call util_pname(filename, actualname)
c
      otrace = util_print('griddebug',print_debug)
      oinfo  = util_print('gridinfo', print_default) .and.
     $     ga_nodeid().eq.0
      oreafail = util_print('grinforead', print_high)
c
      if (otrace) then
         write(LuOut,*) ' grid_file: opening ',
     $        actualname(1:inp_strlen(actualname))
         call util_flush(LuOut)
      end if
c
c     Initialize common
c
      me=ga_nodeid()
      file_opened = .false.
      cur_rec = 1
      cur_cub = -99
      nsubb=0
      ntot_sub=0
      availkb=0
      rec_len_in_dbl = 4*buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
      call grid_buf_clear()
      call dfill(4*buffer_size, 0.0d0, values, 1)
      cur_rec_rad=0d0
      eaftype=1000000+eaf_size_in_dbl
c
c     defer allocation of memory until the first write otherwise
c     this memory is sitting around during the initial guess routines
c
      call util_file_name('grinfo',.true.,.false.,grinfo_in)
      call util_pname(grinfo_in, grinfo)
      if(.not.oreopen) then
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
            if(.not.largenode) then
            call util_directory_name(scratchdir, .true., ga_nodeid())
            sdlen = inp_strlen(scratchdir)
            if (scratchdir(1:sdlen) .eq. actualname(1:sdlen)) then
               availkb = util_scratch_dir_avail_for_me()
            else
               ierr=eaf_stat(actualname, availkb, fstype)
c
c     now eaf_stat returns Mb instead of kb
c
               availkb=availkb*1024d0
               if (ierr .ne. 0) then
                  call eaf_errmsg(ierr,errmsg)
                  write(LuOut,*) ga_nodeid(),errmsg
                   call errquit
     $              ('grid_file_open:unable to stat available disk',0,
     &       DISK_ERR)
                endif
             endif
             endif
c     
c     dont use disk if less than  1MB available
c
             if(availkb.le.1024d0.or.largenode) then
c
c     dont use disk for grid pts
c     
                if (do_zora) then
                  if (ga_nodeid().eq.0) then
                    write(*,*)'WARNING: ZORA *MUST* have gridpts file '
     &                      //'Good luck!'
                  endif
                  largenode=.false.
                else
                  largenode=.true.
                  grid_file_open=.true.
                  return
                endif
             endif
            max_file_rec = nint((0.95*1024.0d0*availkb) / 
     $           (8.0d0*rec_len_in_dbl))
            max_file_rec = max(0,max_file_rec)
c     
            if (filesize.gt.0) max_file_rec = 
     =           min(max_file_rec,filesize/rec_len_in_dbl)
         else 
             call errquit(
     C           'gridfileopen: filesize lt 0',911, DISK_ERR)
            
         endif
c     
         if (max_file_rec.eq.0 .and. max_cache_rec.eq.0) then
            grid_file_open = .false.
            return
         endif
c        l_cache = -1
         k_cache = -1
         n_rec_in_file = 0
         grid_file_open = .true.
         if(.not.grid_grinfoinit(fdg,grinfo,enuc_in,'write')) call
     C        errquit(' opening1 grinfo file failed',0,0)
         filedump_opened=.true.
      else
         grid_file_open=.false.
         if (otrace) then
            write(LuOut,*) ' grid_file: grinfo opening ',grinfo
            call util_flush(LuOut)
         end if
cdbg         enuc=0d0
         if(.not.grid_grinfoinit(fdg,grinfo,enuc,'read')) goto 199

         filedump_opened=.true.
         if (otrace) then
            write(LuOut,*) ' grid_file: enucs ',enuc_in,enuc
            call util_flush(LuOut)
         end if
         if(abs(enuc-enuc_in).gt.1d-7) then
            if (otrace)
     T           write(Luout,*) ' GRIDFILE: wrong enucs ',enuc,enuc_in
c     de-allocate heap
#ifdef EAFHACK
            if(.not.grid_ptsdel(fd)) call errquit(
     A           ' grfileopen: grid_ptsdel failed fd ',fd,HEAP_ERR)
#endif
            return
	 endif
         if (otrace) then
            write(LuOut,*) ' readin mxfilerec from ',grinfo
            call util_flush(LuOut)
         end if
         if(.not.grid_grinfo2(fdg,grinfo,'read',
     I        max_file_rec, max_cache_rec,n_rec_in_file,
     L        file_opened)) goto 199
         grid_file_open=.true.
         if (.not. grid_file_rewind())
     $        call errquit('grfile: rewinding gridpts?',1, DISK_ERR)
cbugfix for idle procs
         if(n_rec_in_file.le.0) return

         file_opened=.false.
         filedump_opened=.false.
#ifndef EAFHACK
         ierr=eaf_open(actualname, eaf_rw, fd)
#else
            ierr=-1
            if(grid_ptsalloc(fd,eaf_size_in_dbl,
     C           'gridpts'))ierr=0
#endif
         if (ierr .ne. 0) then
            call eaf_errmsg(ierr,errmsg)
            write(LuOut,*) ga_nodeid(),errmsg
            call errquit('grid_file_open: eaf_open failed', ierr,
     &       DISK_ERR)
         endif
         file_opened = .true.
      endif
      return
 199  continue
      if(oreafail)
     W  write(luout,*) ga_nodeid(),' grinfo reading failed',ierr
      return
c
c
      end
      logical function grid_file_close(okeep)
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      logical okeep  ! [input] If true keep, otherwise delete
      logical oprint,lstatus
      integer ierr
      double precision ll
      logical grid_grinfo2,grid_ptsdel,file_closedone
      external grid_grinfo2,grid_ptsdel
      logical grid_grinfoinit
      external grid_grinfoinit
      oprint=util_print('iostats',print_high).and.ga_nodeid().eq.0
c
      if (file_opened) then
         if (oprint) 
     $        call eaf_print_stats(fd)
#ifndef EAFHACK
         if (eaf_close(fd) .ne. 0) call errquit
     $        ('grid_file_close: eaf_close failed',0, DISK_ERR)
#endif
         file_closedone=.true.
      endif
c
      if(okeep) then
         if(file_opened) then
            if(.not.grid_grinfo2(fdg,grinfo,'write',
     I           max_file_rec, max_cache_rec,n_rec_in_file,
     L           file_opened)) goto 199
cbad        for acr cam3lyp cdfit
         fdg=-1
         endif
      else
c invalidate enuc when deleting grid file
         if(.not.grid_grinfoinit(fdg,grinfo,-123d0,'write')) call
     C        errquit(' opening2 grinfo file failed',0,0)
#ifndef EAFHACK
         if(fd.ne.-1) then
         if(.not.file_closedone) then
         if (eaf_close(fd) .ne. 0) call errquit
     $        ('grid_file_close: eaf_close2 for fd ',fd, DISK_ERR)
      endif
            if (eaf_delete(actualname) .ne. 0) call errquit
     $           ('grid_file_close: delete failed', 0, DISK_ERR)
         endif
#else
            if (.not.grid_ptsdel(fd)) call errquit
     $           ('grid_file_close: delete failed', 1, MA_ERR)
#endif
         fdg=-1
      endif
c
#ifdef OLD
         call ga_igop(msg_grid_stat1, n_rec_in_file, 1, '+')
#endif
         if (oprint) then
            write(LuOut,1) n_rec_in_file
 1          format(/
     $           ' Parallel grid_pts file used',i8,' records'/)
            call util_flush(LuOut)
         end if
 2001    continue

         if (otrace) then
            write(LuOut,*) ' grid_file: closing ',
     $           actualname(1:inp_strlen(actualname))
            call util_flush(LuOut)
         end if
c
c     l_cache = -1
      file_opened = .false.
      cur_rec = 0
      n_rec_in_file = 0
      call grid_buf_clear
      max_cache_rec = 0
cdbg      max_file_rec = 0
c     
      grid_file_close = .true.
c
      return
 199  write(0,*) ' write failure on grinfo: ierr=',ierr
      grid_file_close=.false.
      return
      end
      logical function grid_file_rewind()
      implicit none
#include "cgridfile.fh"
#include "inp.fh"
#include "stdio.fh"
c
c     rewind the file ready to read grid_pts ... hide any header.
c
      cur_rec = 1
      call grid_buf_clear
c
      if (otrace) then
         write(LuOut,*) ' grid_file: rewinding ', 
     $        actualname(1:inp_strlen(actualname))
         call util_flush(LuOut)
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
#include "errquit.fh"
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "global.fh"
      logical oend           ! [input]
      logical grid_opencache,status,grid_ptsalloc,grid_file_close
      external grid_opencache,grid_ptsalloc,grid_file_close
c
c     Write the grid_pts buffer to disk.  For all but the
c     end record it is called with oend=.false. meaning
c     that size constraints should be enforced leaving
c     room for the end record.  The last record is output with
c     oend=.true. meaning just write it.
c
      integer maxrec
      double precision offset
      integer ierr
      character*80 errmsg
      integer grid_ptswrite
      external grid_ptswrite
c
c     test for overflowing the file and memory cache
c
      maxrec = max_cache_rec + max_file_rec
      if ((.not. oend) .and. (cur_rec .gt. maxrec)) then
         grid_packed_buf_write = .false.
         if (otrace) then
            write(LuOut,12) cur_rec
 12         format(' grid_file: file full at record ',i8)
            call util_flush(LuOut)
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
            write(LuOut,20) actualname(1:inp_strlen(actualname)),
     $           rec_len_in_dbl,n_per_rec,max_cache_rec,max_file_rec
 20         format(/
     $           ' Grid_pts file          = ', a/
     $           ' Record size in doubles = ', i6,7x,
     $           ' No. of grid_pts per rec  = ', i6/
     $           ' Max. records in memory = ', i6,7x,
     $           ' Max. recs in file   = ', i9/)
            call util_flush(LuOut)
         end if
      end if
c
 10   continue
         if (.not. file_opened) then
            file_opened = .true.
#ifndef EAFHACK
            ierr=eaf_open(actualname, eaf_rw, fd)
#else
            ierr=-1
            if(grid_ptsalloc(fd,eaftype-1000000,
     C           'gridpts')) ierr=0
#endif
         if(ierr.ne.0)  call 
     E        errquit('grid_pckbufw: eaf_open failed', 0,
     &        DISK_ERR)
         endif
       if (cur_rec .le. max_cache_rec) then
          if(.not.grid_opencache()) goto 10 ! Dirty ... go back and try again
         call util_memcpy(
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), buf,
     .        rec_len_in_dbl*8)
c     
c     dump cache to disk
c
               offset = 8.0d0*rec_len_in_dbl*(cur_rec-1)
#ifndef EAFHACK
               ierr = eaf_write(fd, offset, dbl_mb(k_cache+
     +              (cur_rec-1)*rec_len_in_dbl),
     .              8*rec_len_in_dbl)
#else
               ierr=grid_ptswrite(fd, offset, dbl_mb(k_cache+
     +              (cur_rec-1)*rec_len_in_dbl),
     .              8*rec_len_in_dbl)
#endif
               if (ierr .ne. 0) then
                  call eaf_errmsg(ierr, errmsg)
                  write(luout,*) ' availdbl ',MA_inquire_avail(mt_dbl)
                  write(LuOut,*) ' IO offset ', offset
                  write(LuOut,*) ' IO error message ',
     $                 errmsg(1:inp_strlen(errmsg)),
     $                 ' nrecdbl ',8*rec_len_in_dbl
                  call util_flush(luout)
                  call MA_summarize_allocated_blocks()
                  call util_flush(luout)
                  call util_flush(0)
                  call errquit('grid_packed_buf_write: write failed',0,
     &       DISK_ERR)
               endif
      else
c
c     Should really detect failures in writing and return false
c
         offset = 8.0d0*rec_len_in_dbl*(cur_rec-1)
#ifndef EAFHACK
         ierr = eaf_write(fd, offset, buf, 8*rec_len_in_dbl)
#else
         ierr = grid_ptswrite(fd, offset, buf, 8*rec_len_in_dbl)
#endif
         if (ierr .ne. 0) then
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO offset ', offset
            write(LuOut,*) ' IO error message ',
     $                     errmsg(1:inp_strlen(errmsg)),
     $                 '2nrecdbl ',8*rec_len_in_dbl
            call errquit('grid_packed_buf_write: write failed',0,
     &       DISK_ERR)
         endif
      end if
c
      if (otrace) then
         write(LuOut,2) n_in_rec, cur_rec, n_rec_in_file+1
 2       format(' grid_file: wrote ',i6,' grid_pts in record',i5
     $        ,' of ',i5,'.')
c         write(LuOut,*)' buf: '
c         call output(buf, 1, 4, 1, n_in_rec, 4, rec_len_in_dbl, 1)
         call util_flush(LuOut)
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
#include "stdio.fh"
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
            write(LuOut,2) n_in_rec, cur_rec-1, n_rec_in_file
 2          format(' grid_file: read ',i6,' grid_pts in record',i6
     $           ,' of ',i6,'.')
c         write(LuOut,*)' buf: '
c         call output(buf, 1, 4, 1, n_in_rec, 4, rec_len_in_dbl, 1)
         call util_flush(LuOut)
      end if
c
      end
      logical function grid_packed_buf_read()
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "global.fh"
#include "stdio.fh"
c
      double precision offset
      integer ierr
      integer grid_ptsread
      external grid_ptsread
      character*80 errmsg
c
      call grid_buf_clear
      grid_packed_buf_read = .false.
      if (cur_rec .gt. n_rec_in_file) then
cedo         if (otrace) then
            write(LuOut,1) cur_rec, n_rec_in_file
 1          format(' grid_file: EOF at record',i5,' of ',i5,'?')
cedo            call util_flush(LuOut)
cedo         end if
            call errquit(' grid_packed_buf_read: EOF reached ',0,
     &       DISK_ERR)
         return
      end if
c
      grid_packed_buf_read = .true.
      if (cur_rec .le. max_cache_rec.and.(.not.oreopen)) then
c
c     Record is cached in memory
c
         call util_memcpy(buf,
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 
     $        rec_len_in_dbl*8)
      else
c
c     Record is on disk
c
         offset = 8.0d0*rec_len_in_dbl*(cur_rec-1)
c
#ifndef EAFHACK
         ierr = eaf_read(fd, offset, buf, 8*rec_len_in_dbl)
#else
         ierr = grid_ptsread(fd, offset, buf, 8*rec_len_in_dbl)
#endif
         if (ierr.ne.0) then
            call eaf_errmsg(ierr,errmsg)
            write(LuOut,*) ga_nodeid(),errmsg
            grid_packed_buf_read = .false.
            return
         endif
c     
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
      nsubb=0
      cur_cub=-99
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
      call util_memcpy(values,buf, ncopy*8)
      cur_rec_ctr = nint(values(1,1))
      n_in_rec = nint(values(2,1))
      cur_rec_rad = values(3,1)
      nsubb = nint(values(4,1))
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
      values(4,1) = dble(nsubb)
      ntot_sub=ntot_sub+nsubb
c
c     Total dbls to copy is 4*(n_per_rec + 1)
c
      ncopy = 4*(n_per_rec + 1)
      call util_memcpy(buf, values,   ncopy*8)
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
      logical function grid_test_mem(rtdb,  nq_tot)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "cdft.fh"
#include "bas.fh"
#include "cgridfile.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "sym.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer rtdb              ! [input]
      integer nq_tot  ! [input]
      double precision enuc
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
      integer avail, localmemsize
      integer filesize, memsize,eaf_size_in_dbl
      character*255 gridfilename
      double precision nq
      logical ofirst
      data ofirst /.true./      ! If first call to this routine
c     
      grid_test_mem = .false.
c
cbad      if (file_opened) return
      oreopen=.false.
      if (.not. rtdb_get(rtdb, 'grid:filesize',
     $     MT_INT, 1, filesize)) filesize = 0
      if (.not. rtdb_get(rtdb, 'grid:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
      if (.not. rtdb_get(rtdb, 'grid:eaf_size_in_dbl',
     $     MT_INT, 1, eaf_size_in_dbl)) eaf_size_in_dbl = 0
      if (filesize.lt.0 .and. memsize.lt.0) return
c
c     Determine if we should automatically run an in-core calculation
c
      oinfo  = util_print('gridinfo', print_default) .and.
     $     ga_nodeid().eq.0
c
      rec_len_in_dbl = 4*buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
c
c     How much memory?  If not input by user, then try to 
c     allocate memory for 250K gridpts; 250K*4 dble words or 250K*4*8 bytes = 8MB
c
      if (memsize.eq.0)then
         memsize = 10000000
      else
         memsize = memsize/8
      endif
      avail = MA_inquire_avail(mt_dbl)
      if (ga_uses_ma() .and. ga_memory_limited()) 
     $   avail = avail - ga_memory_avail()
      call ga_igop(msg_grid_file, avail, 1, 'min')
      avail = min(avail, memsize)
c     
c     Assume that for the small systems for which can hold the
c     grid_pts, that there is no sparsity except for symmetry.
c     Also, assume that the grid_pts will be evenly distributed
c     apart from a 10% fluctuation
c
c      write(LuOut,*)' nq_tot = ', nq_tot     
c
c     need 4*nq_tot (x,y,z,w)
c
      nq=dble(4*nq_tot)
      nq = nq / dble(ga_nnodes())
      nq = nq * 1.1d0 + dble(2*rec_len_in_dbl) ! Fudge factors for safety
c
      if (int(nq) .lt. avail) then
         localmemsize = max(8*buffer_size,int(nq))
      else
         localmemsize = 0
      endif
c
      if (localmemsize.gt.0 .or. filesize.ge.0) then
         if(.not.geom_nuc_rep_energy(geom, enuc)) 
     .        call errquit(' gridtestmem: enucrep failed ',0, GEOM_ERR)
        call util_file_name('gridpts',.true.,.false.,gridfilename)
         if (.not. grid_file_open(enuc,gridfilename, 
     $        localmemsize, filesize,eaf_size_in_dbl))
     $        call errquit('grid_test_mem: grid_file_open failed',0,
     &       DISK_ERR)
         grid_test_mem = .true.
      endif
c
      ofirst = .false.
c
      end
      logical function grid_file_write(nqpts, nctr, rad, aqbuf,ncube)
c      logical function grid_file_write(nints, i, j, k, l, g)
      implicit none
#include "cgridfile.fh"
#include "stdio.fh"
c
c     insert grid_pts into internal buffer writing to disk if necessary.
c     return true on success. false if insufficient memory/file space.
c     hard fail for any other errors.
c
      integer nqpts
      integer nctr
      double precision aqbuf(4,nqpts)
      double precision rad
      integer ncube ! [in] cube
      logical grid_buf_write
      external grid_buf_write
c
      integer ind,ptr_cube
      save ptr_cube
c
      grid_file_write = .true.
      grid_written=.true.
c
      if (otrace) then
         write(LuOut,'(A,A,4I6)') 
     ,        ' grid_file write: nqpts,nctr,n_in_rec',
     &                  ',next_value', 
     &                                 nqpts,nctr,n_in_rec,next_value
         call util_flush(LuOut)
      end if
c
      n_in_rec      = next_value - 1
      n_rec_in_file = cur_rec - 1
      if(n_in_rec .eq. 1) cur_rec_ctr = nctr
c
c     Check to see if buffer needs dumping before adding more pts.
c
      if ((nctr .eq. 0) .or. 
     &    (nctr. ne. cur_rec_ctr).or.
     .     ((n_in_rec+nqpts.gt.n_per_rec)
c     to avoid excessive fragmentation of a single cube
     .     .and.(n_in_rec.gt.n_per_rec*3/4)
     .     )) then
         grid_file_write = grid_buf_write(.false.)
         cur_rec_ctr = nctr
         if (nctr.eq.0)nctr =  n_rec_in_file
         if (.not. grid_file_write) return
         cur_rec_rad=0d0
      end if
      cur_rec_rad = max(rad,cur_rec_rad)
c     This copy can be made a lot more efficient
      if (n_in_rec .ge. n_per_rec-1) then
         grid_file_write = grid_buf_write(.false.)
         if (.not. grid_file_write) return
      endif
      if(cur_cub.ne.ncube) then
         nsubb=nsubb+1
         cur_cub=ncube
         n_in_rec   = next_value
         next_value = next_value + 1
         values(1,n_in_rec)   = 0d0
         ptr_cube=n_in_rec
         values(2,n_in_rec)   = rad
      endif
c after 1st shell contribution to cube
         do ind = 1, nqpts
            if (n_in_rec .eq. n_per_rec) then
               grid_file_write = grid_buf_write(.false.)
               if (.not. grid_file_write) return
               nsubb=1
               cur_cub=ncube
               n_in_rec   = next_value
               next_value = next_value + 1
               values(1,n_in_rec)   = 0d0
               ptr_cube=n_in_rec
               values(2,n_in_rec)   = rad
            end if
c     
            
         n_in_rec   = next_value
         next_value = next_value + 1
c     
         values(1,ptr_cube)=values(1,ptr_cube)+1d0
         values(1,n_in_rec)   = aqbuf(1,ind)
         values(2,n_in_rec)   = aqbuf(2,ind)
         values(3,n_in_rec)   = aqbuf(3,ind)
         values(4,n_in_rec)   = aqbuf(4,ind)
c
      end do
c
      end
      logical function grid_file_read(maxqpts, nqpts, nctr,rad,aqbuf,
     .     nsub_out)
      implicit none
c
#include "cgridfile.fh"
#include "stdio.fh"
c
      integer maxqpts
      integer nqpts
      integer nctr,nsub_out
      double precision aqbuf(4,maxqpts),rad
      integer ga_nodeid
      external ga_nodeid
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
      call util_memcpy(aqbuf(1,1), values(1,2),   ncopy*8)
      nctr = cur_rec_ctr
      nqpts = n_in_rec - 1
      rad =cur_rec_rad
      nsub_out=nsubb
c
      if (otrace) then
         write(LuOut,*) ' grid_file read: nqpts,nctr,n_in_rec,',
     &                  'next_value', 
     &                                nqpts,nctr,n_in_rec,next_value
         call util_flush(LuOut)
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
      logical function grid_getnsubb(buf_in)
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
c
      double precision buf_in(4,*)
      logical grid_file_rewind,grid_file_read
      external grid_file_rewind,grid_file_read
      integer nqpts_per_buf,nqpts,ictr_buf,iqsh
      double precision rad
c
      grid_getnsubb=.false.
      ntot_sub=0
      if (.not. grid_file_rewind())
     $   call errquit('grid_getnsub: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
      do  iqsh = 1, n_rec_in_file
         if(.not.grid_file_read(nqpts_per_buf, nqpts, ictr_buf, 
     &        rad,buf_in,nsubb))
     .        call errquit('gridgnsub:fileread failed',1, GEOM_ERR)
         ntot_sub=ntot_sub+nsubb
      enddo
      grid_getnsubb=.true.

      return
      end
      logical function grid_reopen(geom)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "cgridfile.fh"
c
      integer geom ![in]
c
      double precision enuc
      integer filesize, memsize,eaf_size_in_dbl
      logical grid_file_open
      external grid_file_open
      character*255 gridfilename
      if(.not.geom_nuc_rep_energy(geom, Enuc)) 
     .     call errquit(' gridquad: enucrep failed ',0, GEOM_ERR)
      grid_reopen=.false.
        call util_file_name('gridpts',.true.,.false.,gridfilename)
      oreopen=.true.
c     extract eaf_size_in_dbl from eaftype
      eaf_size_in_dbl=(eaftype-1000000)
#if 0
cfix for atom case (enuc=0!!) and cg
      if(enuc.le.1d-6) then
         grid_reopen=.false.
         return
      endif
#endif
cdbg      filesize=0
cdbg      memsize=0
      grid_reopen=
     .     grid_file_open(enuc,
     .     gridfilename, memsize, filesize, eaf_size_in_dbl)
      return
      end
