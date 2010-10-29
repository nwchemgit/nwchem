      subroutine hamild(pcap, qcap, fc, fo, s, u, t, h, dos, dt,
     +                  c, smin, qmin, nbb)
C$Id$
      implicit none
c.......................................................................
c
c     construct fock matrices. this is a direct s c f procedure,
c     and the two-electron integrals are recalculated for every
c     iteration.
c.......................................................................
#include "cguess.fh"
c
      integer nbb
      real *8  c(*), smin(nbb,*), qmin(nbb,*)
      real *8  pcap(*), qcap(*), fc(*), fo(*), s(*)
      real *8  u(*), t(*), h(*), dos(*), dt(*)
c
      integer i, j, m, n
      integer nstep1, nstep2, nstep, nsh, nbas1, naddr, j1
      integer k, nosh1
      real *8 occucl, fact1, fact2, term1, term2
c
      call teigd (pcap, qcap, u, t, dt, dos)
c.......................................................................
c     compute smin and qmin
c.......................................................................
      nstep1 = 0
      nstep2 = 0
      nstep = 0
      do 50 i = 1 , nsym
         nsh = ncsh(i) + nosh(i)
         nbas1 = nbas(i)
         do j = 1 , nsh
            naddr = nstep2 + (j-1)*nbas1
            j1 = nstep + j
            do m = 1 , nbas1
               smin(m,j1) = 0.d0
               qmin(m,j1) = 0.d0
               do n = 1 , nbas1
                  k = max(m,n)*(max(m,n)-1)/2 + min(m,n) + nstep1
                  smin(m,j1) = smin(m,j1) + s(k)*c(n+naddr)
                  qmin(m,j1) = qmin(m,j1) + qcap(k)*c(n+naddr)
               enddo
            enddo
         enddo
         nstep = nstep + nsh
         nstep2 = nstep2 + nbas1**2
         nstep1 = nstep1 + n1(i)
 50   continue
c.......................................................................
c     compute fc and fo
c.......................................................................
      k = 1
      nstep = 0
      do i = 1 , nsym
         occucl = 4*i - 2
         nosh1 = nosh(i)
         nsh = ncsh(i) + nosh1
         fact1 = occucl/(occucl-nccup(i))
         fact2 = nccup(i)/(occucl-nccup(i))
         nbas1 = nbas(i)
         do m = 1 , nbas1
            do n = 1 , m
               term1 = 0.d0
               term2 = 0.d0
               do j = 1 , nsh
                  j1 = j + nstep
                  if (j.ne.nsh .or. nosh1.eq.0) then
                     term1 = term1 + smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  else
                     term2 = smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  end if
               enddo
               fo(k) = fact1*term1 + pcap(k) + h(k) - qcap(k)
               fc(k) = fact2*term2 + pcap(k) + h(k)
               k = k + 1
            enddo
         enddo
         nstep = nstep + nsh
      enddo
      return
      end
