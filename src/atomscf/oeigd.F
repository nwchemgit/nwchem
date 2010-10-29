      subroutine oeigd(fc, s, u, t, h)
C$Id$
      implicit none
      real *8 fc(*), s(*), u(*), t(*), h(*)
c.......................................................................
c     one-electron integrals. general.
c.......................................................................
c
#include "cguess.fh"
c
      integer l, i, j
      integer nbas1, k, nstep
      real *8 ufacl, expfac, zp, zq, zpq, term1, ppq, rpq
c
      ufacl = sqrt(8.0d0/3.1415926536d0)
      nstep = 0
      k = 0
      do l = 1 , nsym
         nbas1 = nbas(l)
         expfac = l + 0.5d0
         do i = 1 , nbas1
            do j = 1 , i
               k = k + 1
               zp = zeta(nstep+i)
               zq = zeta(nstep+j)
               zpq = 0.5d0*(zp+zq)
               term1 = sqrt(zpq)
               ppq = zp*zq/zpq
               rpq = sqrt(ppq/zpq)
               s(k) = rpq**expfac
               fc(k) = s(k)
               u(k) = ufacl*s(k)*term1
               t(k) = expfac*s(k)*ppq
               h(k) = t(k) - zn*u(k)
            enddo
         enddo
         nstep = nstep + nbas1
         ufacl = ufacl*2.0d0*l/(2.0d0*l+1.0d0)
      enddo
      return
      end
