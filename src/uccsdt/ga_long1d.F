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
     $     ga_create(mt_dbl, b, nint(n/dble(b))+1, name, b, 0, g_a)
c
      end
      subroutine ga_long1d_put(g_a, start, finish, buf)
      implicit none
#include "errquit.fh"
      integer g_a
      double precision start, finish
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
c     Determine the coordinates of the start and finish
c     in the actual 2-D array
c
      jstart = int((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = int((finish-1d0)/dble(b)) + 1
      iend   = nint(finish - dble(jend-1)*dble(b))
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
         call ga_put(g_a, 1, b, jstart, jend, buf(ioff), b)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         call ga_put(g_a, 1, iend, jend, jend, buf(ioff), 1)
         ioff = ioff + iend
      end if
c
      if ((ioff-1) .ne. nint(finish-start+1d0)) call errquit
     $     ('ga_long1d: lost track of ioff', 0, UNKNOWN_ERR)
c
      end
      subroutine ga_long1d_get(g_a, start, finish, buf)
      implicit none
#include "errquit.fh"
      integer g_a
      double precision start, finish
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
c     Determine the coordinates of the start and finish
c     in the actual 2-D array
c
      jstart = int((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = int((finish-1d0)/dble(b)) + 1
      iend   = nint(finish - dble(jend-1)*dble(b))
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
         call ga_get(g_a, 1, b, jstart, jend, buf(ioff), b)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         call ga_get(g_a, 1, iend, jend, jend, buf(ioff), 1)
         ioff = ioff + iend
      end if
c
      if ((ioff-1) .ne. nint(finish-start+1d0)) call errquit
     $     ('ga_long1d: lost track of ioff', 0, UNKNOWN_ERR)
c
      end
      subroutine ga_long1d_acc(g_a, start, finish, buf, scale)
      implicit none
      integer g_a
      double precision start, finish
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
c     Determine the coordinates of the start and finish
c     in the actual 2-D array
c
      jstart = int((start-1d0)/dble(b)) + 1
      istart = nint(start - dble(jstart-1)*dble(b))
      jend   = int((finish-1d0)/dble(b)) + 1
      iend   = nint(finish - dble(jend-1)*dble(b))
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
         call ga_acc(g_a, 1, b, jstart, jend, buf(ioff), b, scale)
         ioff = ioff + b*(jend-jstart+1)
      end if
c
      if (iend .ne. b) then 
         jend = jend + 1        ! Since decremented above
         call ga_acc(g_a, 1, iend, jend, jend, buf(ioff), 1, scale)
         ioff = ioff + iend
      end if
c
      end
      subroutine ga_long1d_test()
      implicit none
#include "errquit.fh"
c
      double precision rn, start, finish
      integer n, i, loop
      parameter (n = 10000)
      double precision mirror(n), test(n), xx
      integer g_a
c
      logical ga_long1d_create
      double precision util_random
c
c     Create a 1-d array of length n
c
      rn = n
      if (.not. ga_long1d_create('test', rn, g_a))
     $     call errquit('ga_long1d_test: create?', n, GA_ERR)
c
c     Fill it with 1...n
c
      do i = 1, n
         mirror(i) = dble(i)
      end do
      start = 1
      finish   = n
      call ga_long1d_put(g_a, start, finish, mirror)
c
c     Test random sections for put/get
c
      start = util_random(12345) ! Initialize with seed
c
      do loop = 1, 100000
         start = int(util_random(0)*rn) + 1
         if (start .gt. rn) start = rn
         finish   = start + nint(util_random(0)*1024)
         if (finish .gt. rn) finish = rn
         call ga_long1d_get(g_a, start, finish, test(nint(start)))
         do i = nint(start), nint(finish)
            if (test(i).ne.mirror(i)) then
               write(6,*) ' bad ', i, test(i), mirror(i)
               call ga_long1d_get(g_a, dble(i), dble(i), xx)
               write(6,*) ' xx ', xx
               call errquit('bad',0, UNKNOWN_ERR)
            end if
         end do
         do i = nint(start), nint(finish)
            mirror(i) = mirror(i) + util_random(0)*0.001d0
         end do
         call ga_long1d_put(g_a, start, finish, mirror(nint(start)))
      end do
c
      write(6,*) ' PUT GET SEEM TO BE WORKING'
c
c     Now for acc/get
c
      do i = 1, n
         mirror(i) = dble(i)
      end do
      start = 1
      finish   = n
      call ga_long1d_put(g_a, start, finish, mirror)
      do loop = 1, 100000
         start = int(util_random(0)*rn) + 1
         if (start .gt. rn) start = rn
         finish   = start + nint(util_random(0)*1024)
         if (finish .gt. rn) finish = rn
         call ga_long1d_get(g_a, start, finish, test(nint(start)))
         do i = nint(start), nint(finish)
            if (abs(test(i)-mirror(i)) .gt. 1d-6) then
               write(6,*) ' bad ', i, test(i), mirror(i)
               call ga_long1d_get(g_a, dble(i), dble(i), xx)
               write(6,*) ' xx ', xx
               call errquit('bad',0, UNKNOWN_ERR)
            end if
         end do
         do i = nint(start), nint(finish)
            test(i) = util_random(0)*0.001d0
            mirror(i) = mirror(i) + test(i)*1.5d0
         end do
         call ga_long1d_acc(g_a, start, finish, test(nint(start)),1.5d0)
      end do
      write(6,*) ' ACC GET SEEM TO BE WORKING'
c
      end

         
      
c $Id$
