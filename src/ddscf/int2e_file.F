      logical function int2e_file_open(
     $     filename,
     $     cachesize, filesize, precision_arg, orestart)
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      character*(*) filename
      integer cachesize, filesize
      double precision precision_arg
      logical orestart
#include "cint2esave.fh"
c
c     Open file returning CHEMIO filehandle (not yet)
c     Return true on success, false on failure
c
c     Currently only allows for restart if entire file is disk resident
c     ... needs mods to allow writing of memory cached integrals
c     (actually gets simpler once the file routines do the caching)
c
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
      file_opened = .false.
      nlarge = 0
      cur_rec = 1
      n_rec_in_file = 1
      rec_len_in_dbl = buffer_size + ma_sizeof(mt_int, 1, mt_dbl)
      call int2e_buf_clear()
      call dfill(rec_len_in_dbl, 0.0d0, values, 1)
c
c     work out dynamic ranges and scaling for integral compression
c
      precision = precision_arg
      ocompress = precision .gt. 1d-10
      if (ocompress) then
         scale = 1.0d0/precision
         maxvalue = 2147483647*precision ! Max value compreses to 2^31-1
         n_per_rec = buffer_size
      else
         scale = 1.0d0
         maxvalue = 1.0d99
         n_per_rec = 2*(buffer_size/3)
         n_per_rec = n_per_rec - mod(n_per_rec,2) ! Make it even
      endif
c
c     Figure out sizes
c
      if (cachesize .gt. 0) then
         max_cache_rec = cachesize/rec_len_in_dbl
      else
         max_cache_rec = 0
      endif
      if (filesize .gt. 0) then
         max_file_rec = filesize/rec_len_in_dbl
      else if (filesize .eq. 0) then
         max_file_rec = 2048    ! Default max size 
      else 
         max_file_rec = 0
      endif
c
c     Figure out the size of an existing file
c
      if (orestart) then
         open(unitno,
     $        file=actualname,
     $        status='old',
     $        form='unformatted',
     $        err=1000)
         goto 1001
 1000    call errquit('int2e_file_open: restart failed',0)
 1001    file_opened = .true.
c
         n_rec_in_file = 0
 1010    read(unitno, err=1011, end=1012)
         n_rec_in_file = n_rec_in_file + 1
         goto 1010
 1011    call errquit('int2e_file_open: error reading file',0)
 1012    rewind unitno
         if (oinfo) then
            write(6,1020) n_rec_in_file
 1020       format(/' Opened integral file with',i6,' records'/)
            call util_flush(6)
         endif
      endif
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
      logical function int2e_file_close()
      implicit none
#include "cint2efile.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cint2esave.fh"
c
      if (file_opened) close(unitno)
c
      call ga_igop(msg_int2e_stat1, n_rec_in_file, 1, '+')
      call ga_igop(msg_int2e_stat2, nlarge, 1, '+')
      if (oinfo) then
         write(6,1) n_rec_in_file, nlarge
 1       format(/
     $        ' Closing integral file ... used ',i6,' records with ',
     $        i6,' large values'/)
         call util_flush(6)
      endif
      cur_rec = 0
      n_rec_in_file = 0
      call int2e_buf_clear
c
      if (max_cache_rec .gt. 0) then
         if (.not. ma_free_heap(l_cache))
     $        call errquit('int2e_file_close: free heap?',0)
      endif
      max_cache_rec = 0
c
      if (otrace) then
         write(6,*) ' int2e_file: closing ',
     $        actualname(1:inp_strlen(actualname))
         call util_flush(6)
      end if
c
      int2e_file_close = .true.
c
      end
      logical function int2e_file_rewind()
      implicit none
#include "cint2efile.fh"
#include "inp.fh"
#include "cint2esave.fh"
c
c     rewind the file ready to read integrals ... hide any header.
c
      rewind unitno
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
#include "cint2esave.fh"
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
         int2e_set_bf_range = int2e_buf_write()
         if (.not. int2e_set_bf_range) return
      end if
c
c     Pack the previous counter into three 8 bit integers
c
      call int2e_buf_cntr_pack()
c
c     Store ranges so that can survive packing to 1 byte
c
      labels(1, next_value  ) = iand(rshift(ilo,8),255)
      labels(2, next_value  ) = iand(rshift(jlo,8),255)
      labels(3, next_value  ) = iand(rshift(klo,8),255)
      labels(4, next_value  ) = iand(rshift(llo,8),255)
      labels(1, next_value+1) = iand(rshift(ihi,8),255)
      labels(2, next_value+1) = iand(rshift(jhi,8),255)
      labels(3, next_value+1) = iand(rshift(khi,8),255)
      labels(4, next_value+1) = iand(rshift(lhi,8),255)
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
#include "cint2esave.fh"
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
         call errquit('int2e_file_get_bf_range: invalid', 0)
      endif
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
#include "cint2esave.fh"
#include "itri.fh"
c
      int2e_file_write = .true.
c
      if (otrace) then
         write(6,*) ' int2e_file: write file ', nints, cntr_ptr
         call print_integ_list(nints, i, j, k, l, g)
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
               int2e_file_write = int2e_buf_write()
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
         endif
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
#include "cint2esave.fh"
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
            int2e_file_write_big = int2e_buf_write()
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
      endif
c
      end
      logical function int2e_file_read(maxints, nints, i, j, k, l, g)
      implicit none
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
#include "cint2esave.fh"
c
      int2e_file_read = .true.
      nints = 0
      nfree = maxints
c
 10   if (nfree .gt. 0) then    ! while (nfree > 0)
         if (next_value .gt. n_in_rec) then
            if (.not. int2e_buf_read()) goto 1000 ! EOF
         end if
         if (next_value .gt. n_in_rec) call errquit('i2efr: uh?', 0)
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
         enddo
         nints = nints + ncopy
         nfree = nfree - ncopy
         nleft_in_range = nleft_in_range - ncopy
         next_value = next_value + ncopy
c
         goto 10
      endif
c
 1000 if (nints .eq. 0) int2e_file_read = .false.
c
      if (otrace) then
         write(6,*) ' int2e_file: read file ', nints
         call util_flush(6)
*         if (nints .gt. 0) call print_integ_list(nints, i, j, k, l, g)
      end if
c
      return
c
      end
      subroutine int2e_file_fock_block(nfock, tol2e,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk)
      implicit none
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
#include "cint2esave.fh"
c
c     While (integrals available) read buffers
c
 10   if (next_value .gt. n_in_rec) then
         if (.not. int2e_buf_read()) goto 1000 ! EOF
      end if
      if (next_value .gt. n_in_rec) call errquit('i2effb: uh?', 0)
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
#include "cint2efile.fh"
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
#include "cint2esave.fh"
c
c     While (integrals available) read buffers
c
 10   if (next_value .gt. n_in_rec) then
         if (.not. int2e_buf_read()) goto 1000 ! EOF
      end if
      if (next_value .gt. n_in_rec) call errquit('i2effb: uh?', 0)
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
      logical function int2e_buf_write()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "cint2esave.fh"
c
c     Opening of file, allocation of buffer and writing of
c     info message are all now deferred until the first
c     write operation
c
      if (cur_rec .eq. 1) then
         if (oinfo) then
            write(6,20) actualname(1:inp_strlen(actualname)),
     $           rec_len_in_dbl,max_cache_rec,max_file_rec, ocompress
 20         format(/
     $           ' Integral file          = ', a/
     $           ' Record size in doubles = ', i6/
     $           ' Max. records in memory = ', i6/
     $           ' Max. records in file   = ', i6/
     $           ' Compression            = ', l1/)
            call util_flush(6)
         endif
      endif
c
c     should test here for overflowing the file/memory cache
c
      call int2e_buf_cntr_pack()
      call int2e_buf_pack()
c
 10   continue
      if (cur_rec .le. max_cache_rec) then
c
         if (l_cache .eq. -1) then
            if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $           'int2e', l_cache, k_cache)) then
               write(6,*) ' Failed to allocate cache ... disabling'
               max_cache_rec = 0
               goto 10          ! Dirty ... go back and try again
            endif
         endif
         call dcopy(rec_len_in_dbl, buf, 1, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1)
      else
         if (.not. file_opened) then
            open(unitno,
     $           file=actualname,
     $           status='unknown',
     $           form='unformatted')
            file_opened = .true.
         endif
         call swrite(unitno, buf, rec_len_in_dbl)
      endif
c
      if (otrace) then
         write(6,2) n_in_rec, cur_rec, n_rec_in_file+1
 2       format(' int2e_file: wrote ',i4,' integrals in record',i4
     $        ,' of ',i4,'.')
         call util_flush(6)
      end if
c
      call int2e_buf_clear
      n_rec_in_file   = n_rec_in_file + 1
      cur_rec         = cur_rec + 1
      int2e_buf_write = .true.
c
      end
      logical function int2e_buf_read()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "cint2esave.fh"
c
      call int2e_buf_clear
      if (cur_rec .gt. n_rec_in_file) then
         int2e_buf_read = .false.
         if (otrace) then
            write(6,1) cur_rec, n_rec_in_file
 1          format(' int2e_file: EOF at record',i4,' of ',i4,'?')
            call util_flush(6)
         end if
         return
      end if
c
      if (cur_rec .le. max_cache_rec) then
         call dcopy(rec_len_in_dbl, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1,
     $        buf, 1)
      else
         call sread(unitno, buf, rec_len_in_dbl)
      endif
      call int2e_buf_unpack()
c
c     Extract the count of values in the current range in this buffer
c
      call int2e_buf_cntr_unpack()
c
      cur_rec = cur_rec + 1
      if (n_in_rec .gt. 1) then
         int2e_buf_read = .true.
      else
         int2e_buf_read = .false.
      end if
c     
      if (otrace) then
            write(6,2) n_in_rec, cur_rec, n_rec_in_file
 2          format(' int2e_file: read ',i4,' integrals in record',i4
     $           ,' of ',i4,'.')
         call util_flush(6)
      end if
c
      end
      subroutine int2e_buf_clear()
      implicit none
#include "cint2efile.fh"
#include "cint2esave.fh"
c
c     empty the integral buffer ... only the null count of integrals
c     in the currnt range in the buffer is present and it is zero
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
      subroutine int2e_buf_unpack()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
c
      integer i, nint_per_dbl, ind, v1, v2, mask, two31
#include "bitops_decls.fh"
#include "cint2esave.fh"
#include "bitops_funcs.fh"
c
      nint_per_dbl = ma_sizeof(mt_dbl, 1, mt_int)
c
c     Copy count
c
      n_in_rec = n_in_buf
c     
c     Copy compressed values and labels 
c
      if (.not. ocompress) then
         call dcopy(n_per_rec, buf, 1, values, 1)
         call util_unpack_8(4*n_per_rec, buf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 2) then
         do i = 1, n_per_rec
            values(i) = precision * ibuf(i)
         enddo
         call util_unpack_8(4*n_per_rec, ibuf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         mask  = lshift(1,32)-1  ! 2^32 - 1
         two31 = lshift(1,31)-1
         do i = 1, n_per_rec/2
            v1 = iand(rshift(ibuf(i),32),mask) - two31
            v2 = iand(ibuf(i),mask) - two31
            values(ind  ) = precision*v1
            values(ind+1) = precision*v2
            ind = ind + 2
         enddo
         call util_unpack_8(4*n_per_rec, ibuf(n_per_rec/2+1), labels)
      else
         call errquit('int2e_buf_unpack: wierd nint_per_dbl',
     $        nint_per_dbl)
      endif
c
      end
      subroutine int2e_buf_pack()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
c
      integer i, nint_per_dbl, ind, v1, v2, two31
#include "bitops_decls.fh"
#include "cint2esave.fh"
#include "bitops_funcs.fh"
c
      nint_per_dbl = ma_sizeof(mt_dbl, 1, mt_int)
c
c     Copy count
c
      n_in_buf = n_in_rec
c
c     Copy compressed values and labels 
c
      if (.not. ocompress) then
         call dcopy(n_per_rec, values, 1, buf, 1)
         call util_pack_8(4*n_per_rec, buf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 2) then
         do i = 1, n_per_rec
            ibuf(i) = nint(values(i)*scale)
         enddo
         call util_pack_8(4*n_per_rec, ibuf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         two31 = lshift(1,31)-1
         do i = 1, n_per_rec/2
            v1 = nint(values(ind  )*scale) + two31
            v2 = nint(values(ind+1)*scale) + two31
            ibuf(i) = ior(lshift(v1,32),v2)
            ind = ind + 2
         enddo
         call util_pack_8(4*n_per_rec, ibuf(n_per_rec/2+1), labels)
      else
         call errquit('int2e_buf_pack: wierd nint_per_dbl',nint_per_dbl)
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
#include "cint2esave.fh"
#include "bitops_funcs.fh"
c
      n = labels(1,cntr_ptr)
      labels(1,cntr_ptr) = iand(rshift(n,16),255)
      labels(2,cntr_ptr) = iand(rshift(n,8),255)
      labels(3,cntr_ptr) = iand(n,255)
      labels(4,cntr_ptr) = 0
      values(cntr_ptr)   = 0.0d0
c
      end
      subroutine int2e_buf_cntr_unpack()
      implicit none
#include "cint2efile.fh"
c     
      integer i, j, k, l
#include "bitops_decls.fh"
#include "cint2esave.fh"
#include "bitops_funcs.fh"
c
      i = labels(1,cntr_ptr)
      j = labels(2,cntr_ptr)
      k = labels(3,cntr_ptr)
      l = labels(4,cntr_ptr)
c
      nleft_in_range = ior(lshift(ior(lshift(i,8),j),8),k)
      if (i.gt.255 .or. j.gt.255 .or. k.gt.255 .or. l.ne.0 .or.
     $     nleft_in_range.gt.n_in_rec .or.
     $     values(cntr_ptr).ne.0.0d0) then
         write(6,*) i, j, k, l, cntr_ptr, 
     $        nleft_in_range, values(cntr_ptr)
         call util_flush(6)
         call errquit
     $     ('int2e_buf_cntr_unpack: invalid count', 
     $     labels(1,cntr_ptr))
      endif
c
      end
      logical function int2e_test_mem(geom, basis, tol2e)
      implicit none
#include "util.fh"
#include "bas.fh"
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "cfock.fh"
#include "msgids.fh"
c     
      integer geom              ! [input]
      integer basis             ! [input]
      double precision tol2e    ! [input]
c     
c     By default memsize=filesize=0, i.e., a fully direct calculation.
c     However, if we can store the integrals in memory it is much
c     better to do this by default.
c     
c     memsize<0 or filesize<0 means don't use this medium.
c     memsize>0 or filesize>0 means already enabled.
c     
      logical int2e_file_open
      external int2e_file_open
      integer nbf, avail, mem2, max2e, natom, natom_unique
      double precision nq
      logical ofirst
      data ofirst /.true./      ! If first call to this routine
#include "cint2esave.fh"
c     
      int2e_test_mem = .false.
      if (memsize.ne.0 .or. filesize.ne.0) return
      oinfo  = util_print('int2einfo', print_default) .and.
     $     ga_nodeid().eq.0
c
c     How much memory?  Leave space for GAs, for the integral code,
c     and for 50,000 of misc. junk. 
c
      avail = MA_inquire_avail(mt_dbl)
      if (ga_uses_ma()) avail = avail - ga_memory_avail()
      call int_mem_2e4c(max2e, mem2)
      mem2 = max(mem2,max2e)
      avail = avail - 5*max2e - mem2
      avail = avail - 100000
      call ga_igop(msg_int2e_file, avail, 1, 'min')
c
c     Assume that for the small systems for which can hold the
c     integrals, that there is no sparsity except for symmetry.
c     Also, assume that the integrals will be evenly distributed
c     apart from a 10% fluctuation
c
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('int2e: bad bas',0)
      if (.not. geom_ncent_unique(geom, natom_unique)) call errquit
     $     ('int2e_test_mem: bad geom handle?',0)
      if (.not. geom_ncent(geom, natom)) call errquit
     $     ('int2e_test_mem: bad geom handle?',0)
c     
      nq = 0.125d0 * dble(nbf+1)**4 
      nq = nq * (dble(natom_unique) / dble(natom))
      nq = nq / ga_nnodes()
      nq = nq * 1.1d0 + 2*(rec_len_in_dbl+1) ! Fudge factors for safety
c
      if (nq .lt. avail) then
         memsize = int(nq)      ! Set these so that scf_tidy cleans up
         filesize = nq*0.3 + 2*(buffer_size+1) ! Allow for small overflow
         if (.not. int2e_file_open(int2efilename, memsize, filesize, 
     $        tol2e, .false.))
     $        call errquit('int2e_test_mem: int2e_file_open failed',0)
         int2e_test_mem = .true.
      else if (oinfo .and. ofirst) then
         write(6,1) nq-avail
 1       format(/' Provide',1p,d9.2,' more words/process to ',
     $        'cache all integrals in memory'/)
      endif
c
      ofirst = .false.
c
      end
