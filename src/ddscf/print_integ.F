      subroutine print_integ_list(nints, i, j, k, l, g)
      integer nints
      integer i(*), j(*), k(*), l(*)
      double precision g(*)
c     $Id: print_integ.F,v 1.2 1996-05-14 18:00:04 d3g681 Exp $
      integer m
c
c     Debug printing of a list of labelled integrals
c
      do m = 1, nints
         if (abs(g(m)) .gt. 1d-4)
     $        write(66,1) i(m), j(m), k(m), l(m), g(m)
 1       format(4i6,f20.6)
      enddo
c
      end
