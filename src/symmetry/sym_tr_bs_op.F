      double precision function sym_trace_bas_op(r, ang_max, maxf,
     $     n, itype)
C$Id$
      implicit none
c
      integer ang_max, maxf, itype, n
      double precision r(maxf, maxf, 0:ang_max)
c
      integer i
      double precision trace
c
      trace = 0.0d0
      do i = 1, n
         trace = trace + r(i, i, itype)
      enddo
c
      sym_trace_bas_op = trace
c
      end
