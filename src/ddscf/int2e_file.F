      block data cint2efdata
*
* $Id$
*
      implicit none
#include "cint2efile.fh"
      data file_opened /.false./
      data l_cache /-1/
      end
      logical function int2e_file_open(
     $     filename,
     $     cachesize, filesize, precision_arg, orestart)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "eaf.fh"
      character*(*) filename
      integer cachesize, filesize
      double precision precision_arg
      logical orestart
#if defined(CRAY)
      external cint2efdata
#endif
c
      double precision length
      double precision availkb,  util_scratch_dir_avail_for_me
      integer  sdlen
      character*20 fstype
      character*(nw_max_path_len) scratchdir
      external util_scratch_dir_avail_for_me
c
c
c     Open file returning CHEMIO filehandle (not yet)
c     Return true on success, false on failure
c
c     Currently only allows for restart if entire file is disk resident
c     ... needs mods to allow writing of memory cached integrals
c     (actually gets simpler once the file routines do the caching)
c
      if (filename.eq.' ') call errquit('int2e_file_open: empty name',0,
     &       INPUT_ERR)
      call util_pname(filename, actualname)
c
      otrace = util_print('int2edebug',print_debug)
      oinfo  = util_print('int2einfo', print_default) .and.
     $     ga_nodeid().eq.0
c
      if (otrace) then
         write(6,*) ' int2e_file: opening ',
     $        actualname(1:inp_strlen(actualname))
         call util_flush(6)
      end if
c
c     Initialize common
c
      lab_pack_len = 8
      file_opened = .false.
      nlarge = 0
      cur_rec = 1
      n_rec_in_file = 0
      rec_len_in_dbl = buffer_size + ma_sizeof(mt_int, 1, mt_dbl) + 1 !!
      call int2e_buf_clear()
      call dfill(buffer_size, 0.0d0, values, 1)
#ifdef NOAIO
      oasyncio = .false.
#else
      oasyncio = .true.
#endif
c
c     work out dynamic ranges and scaling for integral compression
c
      precision = precision_arg
      ocompress = precision .gt. 1d-9 ! Was 1d-10 but not accurate enuf
**      ocompress = .false.
*      write(6,*) ' precision, ocompress ', precision, ocompress
c
      if (ocompress) then
         scale = 1.0d0/precision
         maxvalue = 2147483647*precision ! Max value compreses to 2^31-1
      else
         scale = 1.0d0
         maxvalue = 1.0d300
      end if
c
      nint_per_dbl = ma_sizeof(mt_dbl, 1, mt_int)
      call int2e_packing_info()
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
     $           ('int2e_file_open:unable to stat available disk',0,
     &       DISK_ERR)
            write(6,*) ' non-SCRATCH ', availkb
c
c     now eaf_stat returns Mb instead of kb
c
            availkb=availkb*1024
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
         int2e_file_open = .false.
         return
      endif
c
c     Figure out the size of an existing file
c
      if (orestart) then
         if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     $        call errquit('int2e_file_open: eaf_open failed', 0,
     &       INT_ERR)
         file_opened = .true.
         if (eaf_length(fd, length) .ne. 0) call errquit
     $        ('int2e_file_open: unable to determine file length',0,
     &       INT_ERR)
c
         n_rec_in_file = length / (8*rec_len_in_dbl)
         if ((length - 8.0d0*dble(n_rec_in_file*rec_len_in_dbl))
     $        .ne. 0)
     $        call errquit('int2e_file: bad file length ', 0, INT_ERR)
c
         call ga_igop(555, n_rec_in_file, 1, '+')
         if (oinfo) then
            write(6,1020) n_rec_in_file
 1020       format(/' Opened integral file with',i8,' records'/)
            call util_flush(6)
         end if
         n_rec_in_file = length / (8*rec_len_in_dbl) ! recompute local len
      end if
c
c     defer allocation of memory until the first write otherwise
c     this memory is sitting around during the initial guess routines
c
      l_cache = -1
      k_cache = -1
c
      int2e_file_open = .true.
c
      end
      subroutine int2e_packing_info
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
c     Given lab_pack_len, nint_per_dbl, and ocompress
c     determine n_per_rec and laboff
c
      if (lab_pack_len .eq. 8) then
         if (ocompress) then
c     value=32bits labels=32bits
            n_per_rec = buffer_size
         else
c     value=64bits labels=32bits
            n_per_rec = 2*(buffer_size/3)
         end if
      else if (lab_pack_len .eq. 16) then
         if (ocompress) then
c     value=32bits labels=64bits
            n_per_rec = 2*(buffer_size/3)
         else
c     value=64bits labels=64bits
            n_per_rec = (buffer_size/2)
         end if
      else
         call errquit('int2e_file_open: bad lab_pack_len', lab_pack_len,
     &       INT_ERR)
      endif
      n_per_rec = n_per_rec - mod(n_per_rec,2) ! Make it a multiple of 2
c
      if (nint_per_dbl .eq. 2) then
         if (ocompress) then
            laboff = n_per_rec + 1
         else
            laboff = 2*n_per_rec + 1
         endif
      else if (nint_per_dbl .eq. 1) then
         if (ocompress) then
            laboff = n_per_rec/2 + 1
         else
            laboff = n_per_rec + 1
         endif
      else
         call errquit('int2e_mem_size: nint_per_dbl?', nint_per_dbl,
     &       INT_ERR)
      endif
c
      end
      logical function int2e_file_close(okeep)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
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
     $        ('int2e_file_close: eaf_close failed',0, UNKNOWN_ERR)
         if (.not. okeep) then
            if (eaf_delete(actualname) .ne. 0) call errquit
     $           ('int2e_file_close: delete failed', 0, INPUT_ERR)
         endif
      endif
c
      if (l_cache .ne. -1) then
         if (.not. ma_free_heap(l_cache))
     $        call errquit('int2e_file_close: free heap?',0, MA_ERR)
      end if
c
c     if (file_opened .or. l_cache.ne.-1) then
c
      call ga_igop(msg_int2e_stat1, n_rec_in_file, 1, '+')
      call ga_igop(msg_int2e_stat2, nlarge, 1, '+')
c
c     Now checking for # of records instead of l_cache
c
      if (n_rec_in_file.gt.0) then
         if (oinfo) then
            write(6,1) n_rec_in_file, nlarge
 1          format(/
     $           ' Parallel integral file used',i8,' records with',
     $           i8,' large values'/)
            call util_flush(6)
         end if
         if (otrace) then
            write(6,*) ' int2e_file: closing ',
     $           actualname(1:inp_strlen(actualname))
            call util_flush(6)
         end if
      endif
c
      l_cache = -1
      file_opened = .false.
      cur_rec = 0
      n_rec_in_file = 0
      call int2e_buf_clear
      max_cache_rec = 0
c     
      int2e_file_close = .true.
c
      end
      logical function int2e_file_rewind()
      implicit none
#include "cint2efile.fh"
#include "inp.fh"
c
c     rewind the file ready to read integrals ... hide any header.
c
      cur_rec = 1
      call int2e_buf_clear
c
      if (otrace) then
         write(6,*) ' int2e_file: rewinding ', 
     $        actualname(1:inp_strlen(actualname))
         call util_flush(6)
      end if
c
      int2e_file_rewind = .true.
c
      end
      logical function int2e_set_bf_range(
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c
c     Subsequent integrals will have labels in these ranges only
c     until another call to this routine
c     Return true if success, false if insufficient memory/file space.
c     Hard fail for any other errors.
c
c     Reads and writes interleave correctly on the buffer.
c     Thus, both use next_value as the place to start writing and
c     reading.  Writes always set n_in_rec to next_value-1.
c     Reads do not change n_in_rec.
c
      logical int2e_buf_write
      external int2e_buf_write
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
      int2e_set_bf_range = .true.
c
c     There will be 5 entries (plus space for one integral)
c
c     1-2 The hi bits of the range
c     3-4 The lo bits of the range
c     5-5 The count of values in this range
c
      if (next_value .gt. (n_per_rec-5)) then
         int2e_set_bf_range = int2e_buf_write(.false.)
         if (.not. int2e_set_bf_range) return
      end if
      if (lab_pack_len .eq. 8) then
         if ((ihi-ilo).ge.255 .or. (jhi-jlo).ge.255 .or.
     $        (khi-klo).ge.255 .or. (lhi-llo).ge.255) then
c
c     The default of 8 bits is not enough.  If first write reset to 16
c
            if (cur_rec.eq.1 .and. n_rec_in_file.eq.0 .and.
     $           next_value.eq.2) then
               lab_pack_len = 16
               call int2e_packing_info()
            else
               write(6,*) cur_rec, n_rec_in_file, next_value
               call errquit('int2e_set_bf_range: more than 255',0,
     &       UNKNOWN_ERR)
            endif
         endif
      endif
c
c     Leave packing here at 8 bits even if are doing 16 bit
c     packing ... just for simplicity
c
c     Pack the previous counter into three 8 bit integers
c
      call int2e_buf_cntr_pack()
c
c     Store ranges so that can survive packing to 1 byte
c
      labels(1, next_value  ) = iand(ishft(ilo,-8),255)
      labels(2, next_value  ) = iand(ishft(jlo,-8),255)
      labels(3, next_value  ) = iand(ishft(klo,-8),255)
      labels(4, next_value  ) = iand(ishft(llo,-8),255)
      labels(1, next_value+1) = iand(ishft(ihi,-8),255)
      labels(2, next_value+1) = iand(ishft(jhi,-8),255)
      labels(3, next_value+1) = iand(ishft(khi,-8),255)
      labels(4, next_value+1) = iand(ishft(lhi,-8),255)
      labels(1, next_value+2) = iand(ilo,255)
      labels(2, next_value+2) = iand(jlo,255)
      labels(3, next_value+2) = iand(klo,255)
      labels(4, next_value+2) = iand(llo,255)
      labels(1, next_value+3) = iand(ihi,255)
      labels(2, next_value+3) = iand(jhi,255)
      labels(3, next_value+3) = iand(khi,255)
      labels(4, next_value+3) = iand(lhi,255)
c
c     Zero counter of # integrals in this block and save pointer to counter
c
      cntr_ptr = next_value+4
      labels(1, cntr_ptr) = 0
c
      next_value = next_value + 5
      n_in_rec   = next_value - 1
c
      if (otrace) then
         write(6,*) ' int2e_file: set range ',
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi, cntr_ptr
         call util_flush(6)
      end if
c
c     Save the range for packing
c      
      range(1,1) = ilo
      range(2,1) = jlo
      range(3,1) = klo
      range(4,1) = llo
      range(1,2) = ihi
      range(2,2) = jhi
      range(3,2) = khi
      range(4,2) = lhi
c
      end
      logical function int2e_get_bf_range(
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c
c     Get range of functions from next logical record which must have been
c     inserted with int2e_file_set_bf_range.  
c     Return true on success, false if end of file.  
c     Hard fail for any other errors.
c
      logical int2e_buf_read
      external int2e_buf_read
c
      int2e_get_bf_range = .true.
c
      if (next_value .gt. (n_in_rec-5)) then
         int2e_get_bf_range = int2e_buf_read()
         if (.not. int2e_get_bf_range) return
      end if
      ilo = labels(1, next_value  )*256 + labels(1, next_value+2) 
      jlo = labels(2, next_value  )*256 + labels(2, next_value+2) 
      klo = labels(3, next_value  )*256 + labels(3, next_value+2) 
      llo = labels(4, next_value  )*256 + labels(4, next_value+2) 
      ihi = labels(1, next_value+1)*256 + labels(1, next_value+3) 
      jhi = labels(2, next_value+1)*256 + labels(2, next_value+3) 
      khi = labels(3, next_value+1)*256 + labels(3, next_value+3) 
      lhi = labels(4, next_value+1)*256 + labels(4, next_value+3)
      cntr_ptr = next_value + 4
c
c     Form the count from 8 bit intgers
c
      call int2e_buf_cntr_unpack()
c
      next_value = next_value + 5
c
      if (otrace) then
         write(6,*) ' int2e_file: get range ',
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi, 
     $        cntr_ptr, nleft_in_range
         call util_flush(6)
      end if
c
      if (ilo.le.0 .or. jlo.le.0 .or. klo.le.0 .or. llo.le.0 .or.
     $     ihi.lt.ilo .or. jhi.lt.jlo .or. khi.lt.klo .or. lhi.lt.llo)
     $     then
         write(6,*) ' int2e_file: get range ',
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi
         write(6,*) ' int2e_file: get range next_value', next_value
         write(6,*) ' int2e_file: get range cntr_ptr', cntr_ptr
         call util_flush(6)
         call errquit('int2e_file_get_bf_range: invalid', 0,
     &       INT_ERR)
      end if
c
      range(1,1) = ilo
      range(2,1) = jlo
      range(3,1) = klo
      range(4,1) = llo
      range(1,2) = ihi
      range(2,2) = jhi
      range(3,2) = khi
      range(4,2) = lhi
c
      end
      logical function int2e_file_write(nints, i, j, k, l, g)
      implicit none
#include "cint2efile.fh"
c
c     insert integrals into internal buffer writing to disk if necessary.
c     return true on success. false if insufficient memory/file space.
c     hard fail for any other errors.
c
c     Reads and writes interleave correctly on the buffer.
c     Thus, both use next_value as the place to start writing and
c     reading.  Writes always set n_in_rec consistent with next_value-1.
c     Reads do not change n_in_rec.
c
      integer nints
      integer i(nints), j(nints), k(nints), l(nints)
      double precision g(nints)
      logical int2e_buf_write, int2e_file_write_big
      external int2e_buf_write, int2e_file_write_big
c
      integer ind
#include "itri.fh"
c
      int2e_file_write = .true.
c
      if (otrace) then
         write(6,*) ' int2e_file: write file ', nints, cntr_ptr
**         call print_integ_list(nints, i, j, k, l, g)
         call util_flush(6)
      end if
c
c     This invalidates all data beyond this point to reflect the
c     truncating effect of writing in the middle of a sequential file
c
c     THIS WILL ONLY WORK FOR FILE TRUNCATION BETWEEN RANGES
c
      n_in_rec      = next_value - 1
      n_rec_in_file = cur_rec - 1
c
c     This copy can be made a lot more efficient
c     
      do ind = 1, nints
         if (ocompress .and. (abs(g(ind)) .ge. maxvalue)) then
            int2e_file_write = int2e_file_write_big
     $           (g(ind), i(ind), j(ind), k(ind), l(ind))
            if (.not. int2e_file_write) return
         else
            if (n_in_rec .eq. n_per_rec) then
               int2e_file_write = int2e_buf_write(.false.)
               if (.not. int2e_file_write) return
            end if
c     
            n_in_rec   = next_value
            next_value = next_value + 1
            labels(1,cntr_ptr) = labels(1,cntr_ptr) + 1
c     
            labels(1,n_in_rec) = i(ind) - range(1,1)
            labels(2,n_in_rec) = j(ind) - range(2,1)
            labels(3,n_in_rec) = k(ind) - range(3,1)
            labels(4,n_in_rec) = l(ind) - range(4,1)
            values(n_in_rec)   = g(ind)
         end if
      end do
c
      end
      logical function int2e_file_write_big(gg, i, j, k, l)
      implicit none
#include "cint2efile.fh"
      double precision gg
      integer i, j, k, l
c
c     g is too large to compress to a 32 bit integer while
c     preserving the requested precision.  Store as multiple
c     pieces each less than maxvalue
c
      double precision part, absg, g
      logical int2e_buf_write
      external int2e_buf_write
c
      g = gg                    ! Do not modify arguments
      int2e_file_write_big = .true.
c
 10   absg = abs(g)
      if (absg .gt. precision) then
         part = min(absg,maxvalue)
         if (g .lt. 0.0d0) part = -part
         g = g - part
         nlarge = nlarge + 1
c         
         if (n_in_rec .eq. n_per_rec) then
            int2e_file_write_big = int2e_buf_write(.false.)
            if (.not. int2e_file_write_big) return
         end if
c     
         n_in_rec   = next_value
         next_value = next_value + 1
         labels(1,cntr_ptr) = labels(1,cntr_ptr) + 1
c     
         labels(1,n_in_rec) = i - range(1,1)
         labels(2,n_in_rec) = j - range(2,1)
         labels(3,n_in_rec) = k - range(3,1)
         labels(4,n_in_rec) = l - range(4,1)
         values(n_in_rec)   = part
c
         goto 10
      end if
c
      end
      logical function int2e_file_read(maxints, nints, i, j, k, l, g)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
      integer maxints
      integer nints
      integer i(maxints), j(maxints), k(maxints), l(maxints)
      double precision g(maxints)
c
c     Return up to maxints integrals corresponding to the same block
c     of atoms.  Return in nints the number of integrals read.  
c     The routine is called repeatedly until all integrals from the 
c     block of functions have been consumed when nints is returned as zero.
c
c     Reads and writes interleave correctly on the buffer.
c     Thus, both use next_value as the place to start writing and
c     reading.  Writes always set n_in_rec consistent with next_value-1.
c     Reads do not change n_in_rec.
c
c     Return true if success, false if no more integrals in the
c     block (nints also returned as zero), hard fail for any other errors.
c
      integer loop, ncopy, nfree, nv4
      integer ibase, jbase, kbase, lbase
      logical int2e_buf_read
      external int2e_buf_read
c
      int2e_file_read = .true.
      nints = 0
      nfree = maxints
c
 10   if (nfree .gt. 0) then    ! while (nfree > 0)
         if (next_value .gt. n_in_rec) then
            if (.not. int2e_buf_read()) goto 1000 ! EOF
         end if
         if (next_value .gt. n_in_rec) call errquit('i2efr: uh?', 0,
     &       UNKNOWN_ERR)
         if (nleft_in_range .le. 0) goto 1000 ! end of range
c
	 ibase = range(1,1)
	 jbase = range(2,1)
	 kbase = range(3,1)
	 lbase = range(4,1)
         ncopy = min(nfree,nleft_in_range)
         call dcopy(ncopy, values(next_value), 1, g(1+nints), 1)
         nv4 = (next_value-1)*4 + 1
         do loop = 1, ncopy
            i(loop+nints) = labels4(  nv4) + ibase
            j(loop+nints) = labels4(1+nv4) + jbase
            k(loop+nints) = labels4(2+nv4) + kbase
            l(loop+nints) = labels4(3+nv4) + lbase
            nv4= nv4 + 4
         end do
         nints = nints + ncopy
         nfree = nfree - ncopy
         nleft_in_range = nleft_in_range - ncopy
         next_value = next_value + ncopy
c
         goto 10
      end if
c
 1000 if (nints .eq. 0) int2e_file_read = .false.
c
      if (otrace) then
         write(6,*) ' int2e_file: read file ', nints
         call util_flush(6)
**         if (nints .gt. 0) call print_integ_list(nints, i, j, k, l, g)
      end if
c
      return
c
      end
      subroutine int2e_file_fock_block(nfock, tol2e,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
c     Add the integrals directly from the IO buffer into the
c     fock matrices.  This routine eliminates repeated calls
c     to int2e_file_read which causes unecessary cache flushes
c     from redundent data movement.  Works for the current
c     range of indices.  See fock_2e_file.F for usage.
c
      integer nfock
      double precision tol2e
      double precision dij(*),dik(*),dli(*),djk(*),dlj(*),dlk(*)
      double precision fij(*),fik(*),fli(*),fjk(*),flj(*),flk(*)
c
      integer neri
      logical int2e_buf_read
      external int2e_buf_read
c
c     While (integrals available) read buffers
c
 10   if (next_value .gt. n_in_rec) then
         if (.not. int2e_buf_read()) goto 1000 ! EOF
      end if
      if (next_value .gt. n_in_rec) call errquit('i2effb: uh?', 0,
     &       UNKNOWN_ERR)
      if (nleft_in_range .le. 0) goto 1000 ! end of range
c
      neri = nleft_in_range
c
      call fock_2e_mod_label(nfock, tol2e, neri, 
     $     labels(1,next_value), values(next_value),
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c
      nleft_in_range = nleft_in_range - neri
      next_value = next_value + neri
c     
      goto 10
c
 1000 continue
c
      end
      subroutine int2e_file_rep_fock(nfock, nbf, jfac, kfac,
     $     tol2e, dens, fock)
     $     
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c
c     Replicated data version of above routine
c
c     Add the integrals directly from the IO buffer into the
c     fock matrices.  This routine eliminates repeated calls
c     to int2e_file_read which causes unecessary cache flushes
c     from redundent data movement.  Works for the current
c     range of indices.  See fock_2e_file.F for usage.
c
      integer nbf, nfock
      double precision tol2e, jfac(*), kfac(*)
      double precision dens(nfock,nbf*nbf), fock(nfock,nbf*nbf)
c
      integer neri
      logical int2e_buf_read
      external int2e_buf_read
c
c     While (integrals available) read buffers
c
 10   if (next_value .gt. n_in_rec) then
         if (.not. int2e_buf_read()) goto 1000 ! EOF
      end if
      if (next_value .gt. n_in_rec) call errquit('i2effb: uh?', 0,
     &       UNKNOWN_ERR)
      if (nleft_in_range .le. 0) goto 1000 ! end of range
c
      neri = nleft_in_range
c
      call fock_2e_rep_mod_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels(1,next_value), values(next_value), dens, fock)
c
      nleft_in_range = nleft_in_range - neri
      next_value = next_value + neri
c     
      goto 10
c
 1000 continue
c
      end
      logical function int2e_buf_write(oend)
      implicit none
#include "cint2efile.fh"
      logical oend
c
      logical int2e_packed_buf_write
      external int2e_packed_buf_write
c
      if (n_in_rec .gt. 1) then      ! First entry is just the counter
        call int2e_buf_cntr_pack()
        call int2e_buf_pack()
        int2e_buf_write = int2e_packed_buf_write(oend)
      else
        int2e_buf_write = .true.
        call int2e_buf_clear()
      endif
c
      end
      logical function int2e_packed_buf_write(oend)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "eaf.fh"
      logical oend           ! [input]
c
c     Write the integral buffer to disk.  For all but the
c     end record it is called with oend=.false. meaning
c     that size constraints should be enforced leaving
c     room for the end record.  The last record is output with
c     oend=.true. meaning just write it.
c
c     This is done to simplify the logic for semidirect ...
c     we can test for a full file without having to worry
c     about the partial buffer that might be left.
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
         int2e_packed_buf_write = .false.
         if (otrace) then
            write(6,12) cur_rec
 12         format(' int2e_file: file full at record ',i8)
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
            if (ocompress) nbitv = 32
            if (rec_len_in_dbl .gt. 999999 .or.
     $          n_per_rec .gt. 999999 .or.
     $          max_cache_rec .gt. 999999 .or.
     $          max_file_rec .gt. 999999) then
              write(6,21) actualname(1:inp_strlen(actualname)),
     $           rec_len_in_dbl,n_per_rec,max_cache_rec,max_file_rec,
     $           lab_pack_len, nbitv
            else
              write(6,20) actualname(1:inp_strlen(actualname)),
     $           rec_len_in_dbl,n_per_rec,max_cache_rec,max_file_rec,
     $           lab_pack_len, nbitv
            endif
 20         format(/
     $           ' Integral file          = ', a/
     $           ' Record size in doubles = ', i6,7x,
     $           ' No. of integs per rec  = ', i6/
     $           ' Max. records in memory = ', i6,7x,
     $           ' Max. records in file   = ', i6/
     $           ' No. of bits per label  = ', i6,7x,
     $           ' No. of bits per value  = ', i6/)
 21         format(/
     $           ' Integral file          = ', a/
     $           ' Record size in doubles = ', i8,3x,
     $           ' No. of integs per rec  = ', i8/
     $           ' Max. records in memory = ', i8,3x,
     $           ' Max. records in file   = ', i8/
     $           ' No. of bits per label  = ', i8,3x,
     $           ' No. of bits per value  = ', i8/)
            call util_flush(6)
         end if
      end if
c
 10   continue
      if (cur_rec .le. max_cache_rec) then
         if (l_cache .eq. -1) then
            if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $           'int2e', l_cache, k_cache)) then
               write(6,*) ' int2e_file: Failed to allocate cache '//
     $              '... halving'
               max_cache_rec = max_cache_rec/2
               goto 10          ! Dirty ... go back and try again
            end if
         end if
*         call dcopy(rec_len_in_dbl, buf, 1, 
*     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1)
         call util_memcpy(dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 
     $        buf, rec_len_in_dbl*8)
      else
         if (.not. file_opened) then
            ierr=eaf_open(actualname, eaf_rw, fd)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr, errmsg)
               write(6,*) ' IO error message ',
     ,           errmsg(1:inp_strlen(errmsg))
               call errquit('int2e_file_open: eaf_open failed', 0,
     &         UNKNOWN_ERR)
             else
               file_opened = .true.
             endif
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
            call errquit('int2e_packed_buf_write: write failed',0,
     &       UNKNOWN_ERR)
         endif
      end if
c
      if (otrace) then
         write(6,2) n_in_rec, cur_rec, n_rec_in_file+1
 2       format(' int2e_file: wrote ',i6,' integrals in record',i5
     $        ,' of ',i5,'.')
         call util_flush(6)
      end if
c
      call int2e_buf_clear
      n_rec_in_file   = cur_rec ! was n_rec_in_file + 1
      cur_rec         = cur_rec + 1
      int2e_packed_buf_write = .true.
c
      end
      logical function int2e_buf_read()
      implicit none
#include "cint2efile.fh"
      logical int2e_packed_buf_read
      external int2e_packed_buf_read
c
      int2e_buf_read = int2e_packed_buf_read()
      if (.not. int2e_buf_read) return
c
      call int2e_buf_unpack()
      call int2e_buf_cntr_unpack()
c
*      if (n_in_rec .gt. 1) then
*         int2e_buf_read = .true.
*      else
*         int2e_buf_read = .false.
*      end if
c     
      if (otrace) then
            write(6,2) n_in_rec, cur_rec-1, n_rec_in_file
 2          format(' int2e_file: read ',i6,' integrals in record',i6
     $           ,' of ',i6,'.')
         call util_flush(6)
      end if
c
      end
      logical function int2e_packed_buf_read()
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "global.fh"
c
      double precision offset
      integer ierr
      character*80 errmsg
c
      call int2e_buf_clear
      if (cur_rec .gt. n_rec_in_file) then
         int2e_packed_buf_read = .false.
         if (otrace) then
            write(6,1) cur_rec, n_rec_in_file
 1          format(' int2e_file: EOF at record',i5,' of ',i5,'?')
            call util_flush(6)
         end if
         return
      end if
c
      int2e_packed_buf_read = .true.
      if (cur_rec .le. max_cache_rec) then
c
c     Record is cached in memory
c
*         call dcopy(rec_len_in_dbl, 
*     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1,
*     $        buf, 1)
         call util_memcpy(buf, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl),
     $        rec_len_in_dbl*8)
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
               call errquit('int2e_buf_read:premature end of file',
     $              cur_rec, INPUT_ERR)
            endif
         else
c
c     The record has been read asynchronously ... copy it over
c
            ierr = eaf_wait(fd, aioreq)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr,errmsg)
               write(6,*) ga_nodeid(),errmsg
               call errquit('int2e_buf_read:eaf_wait failed',cur_rec,
     &       UNKNOWN_ERR)
            endif
            oiopending = .false.
**            call dcopy(rec_len_in_dbl, buf2, 1, buf, 1)
            call util_memcpy(buf, buf2, rec_len_in_dbl*8)
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
               call errquit('int2e_buf_read:asynch IO failed',cur_rec+1,
     &       UNKNOWN_ERR)
            endif
         endif
      end if
c
      cur_rec = cur_rec + 1
c
      end
      subroutine int2e_buf_clear()
      implicit none
#include "cint2efile.fh"
c
c     empty the integral buffer ... only the count of integrals
c     in the current range in the buffer is present and it is zero
c
      n_in_rec = 1
      cntr_ptr = 1
      nleft_in_range = 0
      labels(1,1) = 0
      labels(2,1) = 0
      labels(3,1) = 0
      labels(4,1) = 0
      next_value = 2
c
      end
      subroutine int2e_file_record_position()
      implicit none
#include "cint2efile.fh"
c
c     Save enuf info so that the file can be repositioned
c     and TRUNCATED at the current location if a subsequent
c     write operations fails.
c
      posinfo(1)  = n_in_rec
      posinfo(2)  = next_value
      posinfo(3)  = cur_rec
      posinfo(4)  = cntr_ptr
      posinfo(5)  = labels(1,cntr_ptr)
      posinfo(6)  = nleft_in_range
      posinfo(7)  = n_rec_in_file
      posinfo(8)  = range(1,1)
      posinfo(9)  = range(2,1)
      posinfo(10) = range(3,1)
      posinfo(11) = range(4,1)
      posinfo(12) = range(1,2)
      posinfo(13) = range(2,2)
      posinfo(14) = range(3,2)
      posinfo(15) = range(4,2)
      posinfo(16) = nlarge
c
      if (otrace) then
         write(6,*) ' int2e_file: recorded position at ', cur_rec
         call util_flush(6)
      endif
c
      end
      subroutine int2e_file_reposition_truncate
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "eaf.fh"
c
c     Truncate file at the saved position.
c
      integer target_rec
      logical int2e_buf_read
      external int2e_buf_read
      double precision newfilelen
c
c     First load the desired record
c     
      target_rec = posinfo(3)
      if (target_rec.lt.1 .or. target_rec.gt. cur_rec) 
     $     call errquit('int2e_reposition : bad target', target_rec,
     &       UNKNOWN_ERR)
c
      if (target_rec .ne. cur_rec) then
c
c     Truncate the physical file as necessary
c
         if (cur_rec .gt. max_cache_rec) then
            if (file_opened) then
               newfilelen = max(0.0d0,
     $              8.0d0*rec_len_in_dbl*(target_rec-max_cache_rec))
               if (eaf_truncate(fd, newfilelen) .ne. 0) call errquit
     $              ('int2e_reposition: file truncate failed', 0,
     &       UNKNOWN_ERR)
            endif
            n_rec_in_file = target_rec
         endif
c
c     Read the record and leave the file positioned to write over it
c
         cur_rec = target_rec
         if (.not. int2e_buf_read()) call errquit
     $        ('int2e_reposition: failed reading record ', cur_rec,
     &       UNKNOWN_ERR)
         cur_rec = target_rec
      endif
c
      if (otrace) then
         write(6,1) n_rec_in_file, target_rec
 1       format('int2e_file: repositioned from ',i5,' to ',i5)
      endif
c
c     target record is now in the buffer. reset everything.
c
      n_in_rec            =  posinfo(1) 
      next_value          =  posinfo(2) 
      cur_rec             =  posinfo(3) 
      cntr_ptr            =  posinfo(4) 
      labels(1,cntr_ptr)  =  posinfo(5) 
      nleft_in_range      =  posinfo(6) 
      n_rec_in_file       =  posinfo(7) 
      range(1,1)          =  posinfo(8) 
      range(2,1)          =  posinfo(9) 
      range(3,1)          =  posinfo(10)
      range(4,1)          =  posinfo(11)
      range(1,2)          =  posinfo(12)
      range(2,2)          =  posinfo(13)
      range(3,2)          =  posinfo(14)
      range(4,2)          =  posinfo(15)
      nlarge              =  posinfo(16)
c
      end
      subroutine int2e_buf_unpack()
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
c
      integer i, ind, v1, v2, mask, two31
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
c
c     Copy count
c
      n_in_rec = n_in_buf
c     
c     Copy compressed values and labels 
c
      if (.not. ocompress) then
         call dcopy(n_in_rec, buf, 1, values, 1)   ! was n_per_rec
      else if (nint_per_dbl .eq. 2) then
         do i = 1, n_in_rec ! was n_per_rec
            values(i) = precision * dble(ibuf(i))
         end do
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         mask  = ishft(1,32)-1  ! 2^32 - 1
         two31 = ishft(1,31)
         two31 = two31-1
         do i = 1, n_per_rec/2
            v1 = iand(ishft(ibuf(i),-32),mask) - two31
            v2 = iand(ibuf(i),mask) - two31
            values(ind  ) = precision*dble(v1)
            values(ind+1) = precision*dble(v2)
            ind = ind + 2
         end do
      else
         call errquit('int2e_buf_unpack: wierd nint_per_dbl',
     $        nint_per_dbl, UNKNOWN_ERR)
      end if
c
      if (lab_pack_len .eq. 8) then
         call util_unpack_8(4*n_per_rec, ibuf(laboff), labels)
      else if (lab_pack_len .eq. 16) then
         call util_unpack_16(4*n_per_rec, ibuf(laboff), labels)
      else
         call errquit('int2e_buf_unpack: lab_pack_len?', lab_pack_len,
     &       UNKNOWN_ERR)
      endif
c
      end
      subroutine int2e_buf_pack()
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
c
      integer i, ind, v1, v2, two31
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
c     Copy count
c
      n_in_buf = n_in_rec
c
c     Optionally compress values
c
      if (.not. ocompress) then
         call dcopy(n_in_rec, values, 1, buf, 1)   ! was n_per_rec
      else if (nint_per_dbl .eq. 2) then
         do i = 1, n_in_rec   ! was n_per_rec
            ibuf(i) = nint(values(i)*scale)
         end do
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         two31 = ishft(1,31)
         two31 = two31-1
         do i = 1, n_per_rec/2
            v1 = nint(values(ind  )*scale) + two31
            v2 = nint(values(ind+1)*scale) + two31
            ibuf(i) = ior(ishft(v1,32),v2)
            ind = ind + 2
         end do
      else
         call errquit('int2e_buf_pack: wierd nint_per_dbl',nint_per_dbl,
     &       UNKNOWN_ERR)
      end if
c
c     Always compress labels
c
      if (lab_pack_len .eq. 8) then
         call util_pack_8(4*n_per_rec, ibuf(laboff), labels)
      else if (lab_pack_len .eq. 16) then
         call util_pack_16(4*n_per_rec, ibuf(laboff), labels)
      else
         call errquit('int2e_buf_pack: bad lab_pack_len', lab_pack_len,
     &       UNKNOWN_ERR)
      endif
c
      end
      subroutine int2e_buf_cntr_pack()
      implicit none
#include "cint2efile.fh"
c
      integer n
c
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
      n = labels(1,cntr_ptr)
      labels(1,cntr_ptr) = iand(ishft(n,-16),255)
      labels(2,cntr_ptr) = iand(ishft(n,-8),255)
      labels(3,cntr_ptr) = iand(n,255)
      labels(4,cntr_ptr) = 0
      values(cntr_ptr)   = 0.0d0
c
      end
      subroutine int2e_buf_cntr_unpack()
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
c     
      integer i, j, k, l
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
      i = labels(1,cntr_ptr)
      j = labels(2,cntr_ptr)
      k = labels(3,cntr_ptr)
      l = labels(4,cntr_ptr)
c
      nleft_in_range = ior(ishft(ior(ishft(i,8),j),8),k)
      if (i.gt.255 .or. j.gt.255 .or. k.gt.255 .or. l.ne.0 .or.
     $     nleft_in_range.gt.n_in_rec .or.
     $     values(cntr_ptr).ne.0.0d0) then
         write(6,*) i, j, k, l, cntr_ptr, 
     $        nleft_in_range, values(cntr_ptr)
         call util_flush(6)
         call errquit
     $     ('int2e_buf_cntr_unpack: invalid count', 
     $     labels(1,cntr_ptr),
     &       UNKNOWN_ERR)
      end if
c
      end
      logical function int2e_test_mem(geom, basis, tol2e)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "bas.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "sym.fh"
#include "cfock.fh"
#include "msgids.fh"
c     
      integer geom              ! [input]
      integer basis             ! [input]
      double precision tol2e    ! [input]
c     
c     By default memsize=filesize=0 implying if the calculation fits into
c     memory, then use memory, otherwise use all available disk space.
c
c     If either size is negative then that medium is not to be used.
c     If either is size is positive then it has been user set.
c     
      logical int2e_file_open
      external int2e_file_open
      integer nbf, avail, mem2, max2e, localmemsize, nops
      double precision nq
      logical ofirst
      data ofirst /.true./      ! If first call to this routine
c     
      int2e_test_mem = .false.
c
      if (file_opened) return
      if (filesize.lt.0 .and. memsize.lt.0) return
c
c     Determine if we should automatically run an in-core calculation
c
      oinfo  = util_print('int2einfo', print_default) .and.
     $     ga_nodeid().eq.0
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('int2e: bad bas',0, BASIS_ERR)
      nops = sym_number_ops(geom) + 1
c
      rec_len_in_dbl = buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
c
c     How much memory?  Leave space for GAs, for the integral code,
c     and for 50,000 of misc. junk. 
c
      if (memsize .eq. 0) then
         avail = MA_inquire_avail(mt_dbl)
*      write(6,*) ' avail 1 ', avail
         if (ga_uses_ma() .and. ga_memory_limited()) 
     $        avail = avail - ga_memory_avail()
*      write(6,*) ' avail 2 ', avail
         call intb_mem_2e4c(max2e, mem2)
*      write(6,*) ' avail 3 ', avail
         avail = avail - 5*max(max2e,maxeri) - mem2
*      write(6,*) ' avail 4 ', avail
         avail = avail - 100000
*      write(6,*) ' avail 5 ', avail
         call ga_igop(msg_int2e_file, avail, 1, 'min')
*      write(6,*) ' avail 6 ', avail
c     
c     Assume that for the small systems for which can hold the
c     integrals, that there is no sparsity except for symmetry.
c     Also, assume that the integrals will be evenly distributed
c     apart from a 10% fluctuation
c     
         nq = 0.125d0 * dble(nbf+1)**4 
*      write(6,*) ' nq 1 ', nq
         if (nops .ne. 1) nq = (1.3 * nq) / nops ! 1.3 empircal from C6H6
*      write(6,*) ' nq 1 ', nq
         nq = nq / ga_nnodes()
*      write(6,*) ' nq 1 ', nq
         nq = nq * 2d0 + 2*rec_len_in_dbl ! Allow for labels and EOF blocks
*      write(6,*) ' nq 1 ', nq
c
         if (nq .lt. avail) then
            localmemsize = max(2*buffer_size,int(nq))
         else
            localmemsize = 0
         endif
      else
         localmemsize = 0
      endif
c
      if (localmemsize.gt.0 .or. filesize.ge.0) then
         if (.not. int2e_file_open(int2efilename, 
     $        localmemsize, filesize, tol2e, .false.))
     $        call errquit('int2e_test_mem: int2e_file_open failed',0,
     &       UNKNOWN_ERR)
         int2e_test_mem = .true.
      endif
*      if (localmemsize.eq.0 .and. oinfo .and. ofirst .and. memsize.eq.0) then
*         write(6,1) nq-avail
* 1       format(/' Provide',1p,d9.2,' more words/process to ',
*     $        'cache all integrals in memory'/)
*      end if
c
      ofirst = .false.
c
      end
