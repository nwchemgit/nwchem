      subroutine print_integ_list(nints, i, j, k, l, g)
      integer nints
      integer i(*), j(*), k(*), l(*)
      double precision g(*)
c     $Id: print_integ.F,v 1.1 1995-12-13 01:25:33 d3g681 Exp $
      integer m
c
c     Debug printing of a list of labelled integrals
c
      do m = 1, nints
         write(6,1) i(m), j(m), k(m), l(m), g(m)
 1       format(4i6,f20.8)
      enddo
c
      end
