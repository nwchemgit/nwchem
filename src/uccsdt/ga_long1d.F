      logical function ga_long1d_create(name, n, g_a)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      character*(*) name
      double precision n
      integer g_a
c
      integer b
      parameter (b = 128)
c
c     Create a 1-D global array in which the length is
c     specified as a double so as to exceed the 32-bit limit
c     (which limits a 1-D GA to just 16 GB).
c
c     Internally the array is stored as (128, n/128) so that
c     arrays up to 2 TB may be managed on 32 bit machines
c
c     The put/get/acc operations below have double precision indices
c
      ga_long1d_create = 
     $     ga_create(mt_dbl, b, nint(n/b)+1, name, b, 0, g_a)
c
      end
      subroutine ga_long1d_put(g_a, start, end, buf)
      implicit none
      integer g_a
      double precision start, end
      double precision buf(*)
c
c     A general operation may resolve into 3-pieces
c
c     An incomplete first column
c     A middle block of multiple complete columns
c     An incomplete final column
c
      integer b
      parameter (b = 128)
c
      integer ilo, ihi, ioff
      integer istart, iend, jstart, jend
c
c     Determine the coordinates of the start and end
c     in the actual 2-D array
c
      jstart = nint((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = nint((end-1d0)/dble(b)) + 1
      iend   = nint(end - dble(jend-1)*dble(b))
c
      ioff = 1
c
      if (istart.ne.1 .or. jstart.eq.jend) then
         ilo = istart
         ihi = b
         if (jstart .eq. jend) ihi = iend
         call ga_put(g_a, ilo, ihi, jstart, jstart, buf(ioff), 1)
         ioff = ioff + (ihi-ilo+1)
         if (jstart .eq. jend) return
         jstart = jstart + 1
      end if
c
      if (iend .ne. b) jend = jend - 1
      if (jstart .le. jend) then 
         call ga_put(g_a, 1, b, jstart, jend, buf(ioff), 1)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         ilo = 1
         if (jstart .eq. jend) ilo = istart
         ihi = iend
         call ga_put(g_a, ilo, ihi, jend, jend, buf(ioff), 1)
      end if
c
      end
      subroutine ga_long1d_get(g_a, start, end, buf)
      implicit none
      integer g_a
      double precision start, end
      double precision buf(*)
c
c     A general operation may resolve into 3-pieces
c
c     An incomplete first column
c     A middle block of multiple complete columns
c     An incomplete final column
c
      integer b
      parameter (b = 128)
c
      integer ilo, ihi, ioff
      integer istart, iend, jstart, jend
c
c     Determine the coordinates of the start and end
c     in the actual 2-D array
c
      jstart = nint((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = nint((end-1d0)/dble(b)) + 1
      iend   = nint(end - dble(jend-1)*dble(b))
c
      ioff = 1
c
      if (istart.ne.1 .or. jstart.eq.jend) then
         ilo = istart
         ihi = b
         if (jstart .eq. jend) ihi = iend
         call ga_get(g_a, ilo, ihi, jstart, jstart, buf(ioff), 1)
         ioff = ioff + (ihi-ilo+1)
         if (jstart .eq. jend) return ! Only one perhaps incomplete column
         jstart = jstart + 1
      end if
c
      if (iend .ne. b) jend = jend - 1
      if (jstart .le. jend) then 
         call ga_get(g_a, 1, b, jstart, jend, buf(ioff), 1)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         call ga_get(g_a, 1, iend, jend, jend, buf(ioff), 1)
      end if
c
      end
      subroutine ga_long1d_acc(g_a, start, end, buf, scale)
      implicit none
      integer g_a
      double precision start, end
      double precision buf(*)
      double precision scale
c
c     A general operation may resolve into 3-pieces
c
c     An incomplete first column
c     A middle block of multiple complete columns
c     An incomplete final column
c
      integer b
      parameter (b = 128)
c
      integer ilo, ihi, ioff
      integer istart, iend, jstart, jend
c
c     Determine the coordinates of the start and end
c     in the actual 2-D array
c
      jstart = nint((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = nint((end-1d0)/dble(b)) + 1
      iend   = nint(end - dble(jend-1)*dble(b))
c
      ioff = 1
c
      if (istart.ne.1 .or. jstart.eq.jend) then
         ilo = istart
         ihi = b
         if (jstart .eq. jend) ihi = iend
         call ga_acc(g_a, ilo, ihi, jstart, jstart, buf(ioff), 1, scale)
         ioff = ioff + (ihi-ilo+1)
         if (jstart .eq. jend) return
         jstart = jstart + 1
      end if
c
      if (iend .ne. b) jend = jend - 1
      if (jstart .le. jend) then 
         call ga_acc(g_a, 1, b, jstart, jend, buf(ioff), 1, scale)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         ilo = 1
         if (jstart .eq. jend) ilo = istart
         ihi = iend
         call ga_acc(g_a, ilo, ihi, jend, jend, buf(ioff), 1, scale)
      end if
c
      end
      subroutine ga_long1d_test()
      implicit none
c
      double precision rn, start, end
      integer n, lo, hi
      parameter (n = 10000)
      double precision buf(n)
      integer g_a
c
      logical ga_long1d_create
      double precision util_random
c
c     Create a 1-d array of length n
c
      rn = n
      if (.not. ga_long1d_create('test', rn, g_a))
     $     call errquit('ga_long1d_test: create?', n)
c
c     Fill it with 1...n
c
      do i = 1, n
         buf(i) = dble(i)
      end do
      start = 1
      end   = n
      call ga_long1d_put(g_a, start, end, buf)
c
c     Destroy local buffer
c
      call dfill(n, 0d0, buf, 1)
c
c     Test random sections
c
      start = util_random(12345) ! Initialize with seed
c
      do loop = 1, 100000
         start = int(util_random(0)*rn) + 1
         if (start .gt. rn) start = rn
         end   = start + nint(util_random(0)*512)
         if (end .gt. rn) end = rn
         call ga_long1d_get(g_a, start, end, 
         
      
