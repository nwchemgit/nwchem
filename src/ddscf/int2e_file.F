      logical function int2e_file_open(
     $     filename,
     $     cachesize, filesize)
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      character*(*) filename
      integer cachesize, filesize
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
      otrace = util_print('int2efile',print_debug)
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
      rec_len_in_dbl = n_per_rec +
     $     ma_sizeof(mt_int, 2, mt_dbl) + 
     $     ma_sizeof(mt_byte, 8*n_per_rec, mt_dbl)
      call int2e_buf_clear()
c
      if (cachesize .gt. 0) then
         max_cache_rec = cachesize/rec_len_in_dbl
      else
         max_cache_rec = 0
      endif
      if (filesize .gt. 0) then
         max_file_rec = filesize/rec_len_in_dbl
      else
         max_file_rec = 100000000 ! Unlimited file storage.
      endif
c
      if (max_cache_rec .gt. 0) then
         if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $        'int2e', l_cache, k_cache)) then
            write(6,*) ' Failed to allocate cache '
            max_cache_rec = 0
         endif
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
      if (next_value .gt. (n_per_rec-1)) then
         int2e_set_bf_range = int2e_buf_write()
         if (.not. int2e_set_bf_range) return
      end if
c
      labels(1, next_value  ) = ilo + range_shift
      labels(2, next_value  ) = jlo + range_shift
      labels(3, next_value  ) = klo + range_shift
      labels(4, next_value  ) = llo + range_shift
      labels(1, next_value+1) = ihi + range_shift
      labels(2, next_value+1) = jhi + range_shift
      labels(3, next_value+1) = khi + range_shift
      labels(4, next_value+1) = lhi + range_shift
c
      next_value = next_value + 2
      n_in_rec   = next_value - 1
c
      if (otrace) then
         write(6,*) ' int2e_file: set range ',
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi
         call util_flush(6)
      end if
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
      if (next_value .gt. (n_in_rec-1)) then
         int2e_get_bf_range = int2e_buf_read()
         if (.not. int2e_get_bf_range) return
      end if
      ilo = labels(1, next_value  ) - range_shift
      jlo = labels(2, next_value  ) - range_shift
      klo = labels(3, next_value  ) - range_shift
      llo = labels(4, next_value  ) - range_shift
      ihi = labels(1, next_value+1) - range_shift
      jhi = labels(2, next_value+1) - range_shift
      khi = labels(3, next_value+1) - range_shift
      lhi = labels(4, next_value+1) - range_shift
      next_value = next_value + 2
c
      if (otrace) then
         write(6,*) ' int2e_file: get range ',
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi
         call util_flush(6)
      end if
c
      if (ilo.le.0 .or. jlo.le.0 .or. klo.le.0 .or. llo.le.0 .or.
     $     ihi.lt.ilo .or. jhi.lt.jlo .or. khi.lt.klo .or. lhi.lt.llo)
     $     call errquit('int2e_file_get_bf_range: invalid', 0)
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
      logical int2e_buf_write
      external int2e_buf_write
c
      integer ind
c
      int2e_file_write = .true.
c
      if (otrace) then
         write(6,*) ' int2e_file: write file ', nints
         call util_flush(6)
      end if
c
c     This invalidates all data beyond this point to reflect the
c     truncating effect of writing in the middle of a sequential file
c
      n_in_rec      = next_value - 1
      n_rec_in_file = cur_rec - 1
c
c     This copy can be made a lot more efficient
c     
      do ind = 1, nints
c
         if (n_in_rec .eq. n_per_rec) then
            int2e_file_write = int2e_buf_write()
            if (.not. int2e_file_write) return
         end if
c
         n_in_rec   = next_value
         next_value = next_value + 1
         labels(1,n_in_rec) = i(ind)
         labels(2,n_in_rec) = j(ind)
         labels(3,n_in_rec) = k(ind)
         labels(4,n_in_rec) = l(ind)
         values(n_in_rec)   = g(ind)
      end do
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
      integer loop
      logical int2e_buf_read
      external int2e_buf_read
c
c     This copy can be made a lot more efficient
c     
      int2e_file_read = .true.
      nints = 0
      do loop = 1, maxints
         if (next_value .gt. n_in_rec) then
            if (.not. int2e_buf_read()) goto 1000
         end if
         if (next_value .gt. n_in_rec) call errquit('i2efr: uh?', 0)
c
c     Negative labels are being used to store other info
c     (e.g.. blocking information) and indicate the end of this
c     block of integrals
c
         if (labels(1,next_value) .gt. range_shift) then
            goto 1000
         else
            nints = nints + 1
            i(nints) = labels(1,next_value)
            j(nints) = labels(2,next_value)
            k(nints) = labels(3,next_value)
            l(nints) = labels(4,next_value)
            g(nints) = values(next_value)
            next_value = next_value + 1
         end if
      end do
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
      call int2e_buf_pack()
c
      if (cur_rec .le. max_cache_rec) then
         call dcopy(rec_len_in_dbl, values, 1, 
     $        dbl_mb(k_cache+(cur_rec-1)*rec_len_in_dbl), 1)
      else
         call swrite(unitno, values, rec_len_in_dbl)
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
      if (cur_rec .gt. n_rec_in_file) then
         int2e_buf_read = .false.
         call int2e_buf_clear
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
     $        values, 1)
      else
         call sread(unitno, values, rec_len_in_dbl)
      endif
      call int2e_buf_unpack()
c
      next_value = 1
      cur_rec = cur_rec + 1
      if (n_in_rec .gt. 0) then
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
c     empty the integral buffer
c
      n_in_rec = 0
      next_value = 1
c
      end
      subroutine int2e_buf_pack()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
c
c     pack labels into 2 bytes each
c
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte) / 2
c
      if (npacked_per_int .eq. 2) then
         ind = 1
         do i = 1, n_per_rec
            ij = ior(ishft(labels(1,i),16),labels(2,i))
            kl = ior(ishft(labels(3,i),16),labels(4,i))
            packed_labels(ind  ) = ij
            packed_labels(ind+1) = kl
            ind = ind + 2
         enddo
      else if (npacked_per_int .eq. 4) then
         do i = 1, n_per_rec
            ij = ior(ishft(labels(1,i),16),labels(2,i))
            kl = ior(ishft(labels(3,i),16),labels(4,i))
            ijkl = ior(ishft(ij,32),kl)
            packed_labels(i) = ijkl
         enddo
      else
         call errquit('int2e_buf_pack: funny integer length',
     $        npacked_per_int)
      endif
c
      end
      subroutine int2e_buf_unpack()
      implicit none
#include "cint2efile.fh"
#include "mafdecls.fh"
c
c     unpack labels from 2 bytes each
c
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte) / 2
c
      if (npacked_per_int .eq. 2) then
         ind = 1
         do i = 1, n_per_rec
            ij = packed_labels(ind)
            kl = packed_labels(ind+1)
            labels(1,i) = ishft(ij,-16)
            labels(2,i) = iand(ij,65535)
            labels(3,i) = ishft(kl,-16)
            labels(4,i) = iand(kl,65535)
            ind = ind + 2
         enddo
      else if (npacked_per_int .eq. 4) then
         do i = 1, n_per_rec
            ijkl = packed_labels(i)
            labels(1,i) = ishft(ijkl,-48)
            labels(2,i) = iand(ishft(ijkl,-32),65535)
            labels(3,i) = iand(ishft(ijkl,-16),65535)
            labels(4,i) = iand(ijkl,65535)
         enddo
      else
         call errquit('int2e_buf_unpack: funny integer length',
     $        npacked_per_int)
      endif
c
      end
