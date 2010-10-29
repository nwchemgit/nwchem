      subroutine fock_eri_label_scale(neri, q4, i, j, k, l, eri)
c
c$Id$
c
      implicit none
      integer neri, i(neri), j(neri), k(neri), l(neri)
      double precision q4
      double precision eri(neri)
c
c     Scale labelled integrals by 
c
c     q4 / [(1 + D(i,j))(1 + D(k,l))(1 + D(ij,kl))]
c
c     where D(a,b) = delta(a,b) = 1 if a==b, 0 otherwise
c
      integer ind
      double precision g
c
      do ind = 1, neri
         g = eri(ind) * q4
         if (i(ind).eq.j(ind)) g = g * 0.5d0
         if (k(ind).eq.l(ind)) g = g * 0.5d0
         if (i(ind).eq.k(ind)) then
            if (j(ind).eq.l(ind)) g = g * 0.5d0
         endif
         eri(ind) = g
      enddo
c
      end




