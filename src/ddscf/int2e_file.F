      logical function int2e_file_open(
     $     filename,
     $     cachesize, filesize, precision_arg)
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      character*(*) filename
      integer cachesize, filesize
      double precision precision_arg
c
c     Open file returning CHEMIO filehandle (not yet)
c     Return true on success, false on failure
c
c     Currently does not allow for restart ... but will soon.
c
      call util_pname(filename, actualname)
      open(unitno,
     $     file=actualname,
     $     status='unknown',
     $     form='unformatted')
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
      cur_rec = 1
      n_rec_in_file = 1
      rec_len_in_dbl = n_per_rec + ma_sizeof(mt_int, 1, mt_dbl)
      call int2e_buf_clear()
      call dfill(rec_len_in_dbl, 0.0d0, values, 1)
c
c     work out dynamic ranges and scaling for integral compression
c
c     Need eventually to be able to diable compression for max precision
c
      if (precision .lt. 1d-11) precision = 1d-11
      precision = precision_arg
      scale = 1.0d0/precision
      maxvalue = 2147483647*precision ! Max value compreses to 2^31-1
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
      else
         max_file_rec = 2048    ! Default max size of about 2 Gb
      endif
c
      if (max_cache_rec .gt. 0) then
         if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $        'int2e', l_cache, k_cache)) then
            write(6,*) ' Failed to allocate cache ... setting to zero'
            max_cache_rec = 0
         endif
      endif
c
      if (oinfo) then
         write(6,20) actualname(1:inp_strlen(actualname)),
     $        rec_len_in_dbl,max_cache_rec,max_file_rec
 20      format(/
     $        ' Integral file          = ', a/
     $        ' Record size in doubles = ', i6/
     $        ' Max. records in memory = ', i6/
     $        ' Max. records in file   = ', i6/)
         call util_flush(6)
      endif
c
      int2e_file_open = .true.
c
      end
      logical function int2e_file_close()
      implicit none
#include "cint2efile.fh"
#include "inp.fh"
#include "mafdecls.fh"
c
      close(unitno)
c
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
c     Reads don't change n_in_rec.
c
      logical int2e_buf_write
      external int2e_buf_write
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
      labels(1, next_value  ) = ishft(ilo,-8)
      labels(2, next_value  ) = ishft(jlo,-8)
      labels(3, next_value  ) = ishft(klo,-8)
      labels(4, next_value  ) = ishft(llo,-8)
      labels(1, next_value+1) = ishft(ihi,-8)
      labels(2, next_value+1) = ishft(jhi,-8)
      labels(3, next_value+1) = ishft(khi,-8)
      labels(4, next_value+1) = ishft(lhi,-8)
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
c     Reads don't change n_in_rec.
c
      integer nints
      integer i(nints), j(nints), k(nints), l(nints)
      double precision g(nints)
      logical int2e_buf_write, int2e_file_write_big
      external int2e_buf_write, int2e_file_write_big
c
      integer ind
c
      int2e_file_write = .true.
c
      if (otrace) then
         write(6,*) ' int2e_file: write file ', nints, cntr_ptr
         call util_flush(6)
      end if
c
c     This invalidates all data beyond this point to reflect the
c     truncating effect of writing in the middle of a sequential file
c
c     THIS WILL ONLY WORK FOR TRUNCATION BETWEEN RANGES
c
      n_in_rec      = next_value - 1
      n_rec_in_file = cur_rec - 1
c
c     This copy can be made a lot more efficient
c     
      do ind = 1, nints
c
         if (abs(g(ind)) .ge. maxvalue) then
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
c
      g = gg                    ! Don't modify arguments
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
c     Reads don't change n_in_rec.
c
c     Return true if success, false if no more integrals in the
c     block (nints also returned as zero), hard fail for any other errors.
c
      integer loop, ncopy, nfree
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
         if (next_value .gt. n_in_rec) call errquit('i2efr: uh?', 0)
         if (nleft_in_range .le. 0) goto 1000 ! end of range
c
         ncopy = min(nfree,nleft_in_range)
         do loop = 1, ncopy
            i(loop+nints) = labels(1,next_value+loop-1) + range(1,1)
            j(loop+nints) = labels(2,next_value+loop-1) + range(2,1)
            k(loop+nints) = labels(3,next_value+loop-1) + range(3,1)
            l(loop+nints) = labels(4,next_value+loop-1) + range(4,1)
            g(loop+nints) = values(next_value+loop-1)
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
      end if
c
      return
c
      end
      logical function int2e_buf_write()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
c
c     should test here for overflowing the file/memory cache
c
      call int2e_buf_cntr_pack()
      call int2e_buf_pack()
c
      if (cur_rec .le. max_cache_rec) then
         call dcopy(rec_len_in_dbl, buf, 1, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1)
      else
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
      integer i, nint_per_dbl, ind, v1, v2, mask
c
      nint_per_dbl = ma_sizeof(mt_dbl, 1, mt_int)
c
c     Copy count
c
      n_in_rec = n_in_buf
c     
c     Copy compressed values and labels 
c
      if (nint_per_dbl .eq. 2) then
         do i = 1, n_per_rec
            values(i) = precision * ibuf(i)
         enddo
         call util_unpack_8(4*n_per_rec, ibuf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         mask = ishft(1,32)-1  ! 2^32 - 1
         do i = 1, n_per_rec/2
            v1 = ishft(ibuf(i),-32)
            v2 = iand(ibuf(i),mask)
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
      integer i, nint_per_dbl, ind, v1, v2
c
      nint_per_dbl = ma_sizeof(mt_dbl, 1, mt_int)
c
c     Copy count
c
      n_in_buf = n_in_rec
c
c     Copy compressed values and labels 
c
      if (nint_per_dbl .eq. 2) then
         do i = 1, n_per_rec
            ibuf(i) = nint(values(i)*scale)
         enddo
         call util_pack_8(4*n_per_rec, ibuf(n_per_rec+1), labels)
      else if (nint_per_dbl .eq. 1) then
         ind = 1
         do i = 1, n_per_rec/2
            v1 = nint(values(ind  )*scale)
            v2 = nint(values(ind+1)*scale)
            ibuf(i) = ior(ishft(v1,32),v2)
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
      n = labels(1,cntr_ptr)
      labels(1,cntr_ptr) = ishft(n,-16)
      labels(2,cntr_ptr) = iand(ishft(n,-8),255)
      labels(3,cntr_ptr) = iand(n,255)
      values(cntr_ptr)   = 0.0d0
c
      end
      subroutine int2e_buf_cntr_unpack()
      implicit none
#include "cint2efile.fh"
c     
      integer i, j, k
c
      i = labels(1,cntr_ptr)
      j = labels(2,cntr_ptr)
      k = labels(3,cntr_ptr)
c
      nleft_in_range = ior(ishft(ior(ishft(i,8),j),8),k)
      if (i.gt.255 .or. j.gt.255 .or. k.gt.255 .or.
     $     nleft_in_range.gt.n_in_rec .or.
     $     values(cntr_ptr).ne.0.0d0) call errquit
     $     ('int2e_buf_cntr_unpack: invalid count', 
     $     labels(1,cntr_ptr))
c
      end
