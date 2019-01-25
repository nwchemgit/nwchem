C
C     Computes how to partition n x n array for parallelization.
C     Returns the corresponding buffer size and start and end indices
C     for the current processor.
C
      subroutine calc_array_partitioning (n, m, istart, iend, bfr_size)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in)   :: n           ! size of array n x n

      
C     == Outputs ==
      integer, intent(out)  :: m           ! size of slice is m x n
      integer, intent(out)  :: istart      ! start index for this processor
      integer, intent(out)  :: iend        ! end index for this processor
      integer, intent(out)  :: bfr_size    ! = n x m


C     == Parameters ==
      character(*), parameter :: pname = "calc_array_partitioning: "

 
C     == Variables ==
      integer :: me, np, m0, ratio
      integer :: sratio, s0, endproc
      integer :: mindim
      integer :: min_slice_width

      
      me = ga_nodeid ()
      np = ga_nnodes ()
      
 
C
C     Determine how much work this current processor does, which is a
C     slice of size m x n.  Since the total matrix size is n x n, each
C     processor will get an n/np x n slice.  If n/np does not divide
C     evenly, the extra will be given to processor 0.
C
C     Find the start and end index for the column index.  That is, this
C     processor is working with an m x n slice, so we need to figure out
C     what part of the full matrix we are working with.  We use fortran
C     notation (starting at 1).  Proc 0 goes from 1 to m0, proc 1 goes
C     from m0 + 1 to m0 + 1 + m, etc.
C
C     However, if m is too small (i.e, smaller matrix and/or larger
C     number of processors) this is inefficient.  So instead we set a
C     minimum amount of work (min_slice_width) which means some procs
C     get nothing to do.
C
      ratio = n/np

      mindim = min(n, m)
      min_slice_width = min(mindim, 10) !if we have a larger array, use a min slice size of 10
      
      if (ratio .ge. min_slice_width) then
         m0 = ratio + mod (n, np)

         if (me.eq.0) then
            m = m0
            istart = 1
            iend = m0
         else
            m = ratio
            istart = m0 + 1 + (me-1)*m
            iend = istart + m - 1
         endif

      else 
         
         endproc = n / min_slice_width
         s0 = min_slice_width + mod (n, min_slice_width)

         if (me.eq.0) then      !do min x n work + remainder on this proc
            m = s0
            istart = 1
            iend = s0
         elseif ((me .gt. 0).and.(me .lt. endproc)) then !do min x n work on this proc, note < since procs start at 0
            m = min_slice_width
            istart = s0 + 1 + (me-1)*m
            iend = istart + m - 1
         else                   !do no work on this processor
            m = 0
            istart = 0
            iend = 0
         endif
         
      endif

      bfr_size = n*m

C     (print processor-dependent parameters)
c$$$      write (luout, "(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a)")
c$$$     $     "Processor ", me, ": ", m, " x ", n, ", [",
c$$$     $     istart, ":", iend, "]"

      end subroutine

c $Id$
