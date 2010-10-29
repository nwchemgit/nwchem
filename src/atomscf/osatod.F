      function osatod(iat,ic,iiloc,iisch,nbb)
C$Id$
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "cguess.fh"
      integer nbb
      integer iat
      integer ic(6,nbb),
      integer iiloc(nbb,6)
      integer iisch(nbb,6)
      logical osatod
c
c...   function checks if this atom is the same as one for which info
c...   is in /cguess/ and from which we can use the density matrix
c...   for atomic startup
c...   if so osatod is true  and iiloc is updated
c...
c
      integer nbcn(5)
c..
      real *8 pi32, fac, ee, contn
      integer i, ii, iorb, ispdf, j, k
      integer mini, maxi, kk, if, is, kkzc, kh, isymax, nbasn
c
      data pi32/5.56832799683170d0/
c..
      osatod = .false.
      if (zan(iat).ne.zn) return
c
c...  gather  shell / symmetry info
c
      do i = 1 , 4
         nbcn(i) = 0
      enddo
c
c.. nbc  # shell's / symmetry
c.. iisch  contains index of shell
c.. iiloc  contains position of starting ao of shell in "real" world
c..  translate to 1 (s)
c
      do ii = 1 , nshell
         i = katom(ii)
         if (i.eq.iat) then
            mini = kmin(ii)
            maxi = kmax(ii)
            kk = ktype(ii)
            if (kk.eq.6) kk = 2
            do iorb = mini , maxi
               if (iorb.eq.1) then
                  nbcn(1) = nbcn(1) + 1
                  iisch(nbcn(1),1) = ii
                  iiloc(nbcn(1),1) = kloc(ii)
               else if (iorb.eq.2 .or. iorb.eq.5 .or. iorb.eq.11) then
c..  translate to 2 (p) 3(d) or  4(f)
                  ispdf = kk
                  nbcn(ispdf) = nbcn(ispdf) + 1
                  iisch(nbcn(ispdf),ispdf) = ii
                  iiloc(nbcn(ispdf),ispdf) = kloc(ii) + iorb - mini
               end if
            enddo
         end if
      enddo
c..     check nbcn
      do i = 1 , 4
         if (nbc(i).ne.nbcn(i)) return
      enddo
c..
c..     we gathered symmetry/shell info ; now check the real thing
c..
      kkzc = 0
      kh = 0
      isymax = 0
      do ispdf = 1 , 4
c..      nbas = total # primitives for this symmetry
         nbasn = 0
         if (nbc(ispdf).gt.0) isymax = ispdf
         do j = 1 , nbc(ispdf)
            ii = iisch(j,ispdf)
            is = kstart(ii)
            if = is + kng(ii) - 1
c..      ic = # number of primitives /contracted /symmetry
            if (ic(ispdf,j).ne.kng(ii)) return
            nbasn = nbasn + kng(ii)
c..      check the primitives / watch the subtle use of 2-dim cspd
            do k = is , if
               kkzc = kkzc + 1
               if (zeta(kkzc).ne.ex(k)) return
               contn = cspd(k,ispdf)
c...     get contraction coeff's as we are used to
               ee = 2*zeta(kkzc)
               fac = pi32/(ee*sqrt(ee))
               if (ispdf.eq.2) then
                  fac = 0.5d0*fac/ee
               else if (ispdf.eq.3) then
                  fac = 0.75d0*fac/(ee*ee)
               else if (ispdf.eq.4) then
                  fac = 1.875d0*fac/(ee**3)
               end if
               contn = contn*sqrt(fac)
               if (cont(kkzc).ne.contn) return
            enddo
         enddo
c...
         if (nbasn.ne.nbas(ispdf)) return
c...
      enddo
c..
c..      all checks out
c..
      osatod = .true.
c..
      return
      end
