      subroutine tracd(c, cc, nsqt, nsym, nbas, ncsh, nosh, nbc,
     +      nstrt, cont)
C$Id$
      implicit none
c.......................................................................
c
c     transform vectors from contr. to primitive basis functions.
c.......................................................................
      real *8 c(*),cc(*), cont(*)
      integer nsqt, nsym, nbas(*), ncsh(*), nosh(*), nbc(*), nstrt(*)
c
      integer i, l, n, k, j
      integer k1, k2, ii
      integer nstep, nstep1, nstep2, nbc1, nbasl, nsh1
c
      do i = 1 , nsqt
         c(i) = 0.0d0
      enddo
c.......................................................................
c
c     cc(i,j)  contains vectors over cont. functions.
c              i runs over functions, j runs over orbitals
c.......................................................................
      nstep = 0
      nstep1 = 0
      nstep2 = 0
      do l = 1 , nsym
         nbc1 = nbc(l)
         nbasl = nbas(l)
         nsh1 = ncsh(l) + nosh(l)
         ii = 1
         do n = 1 , nbc1
            k1 = nstrt(nstep+n)
            k2 = nstrt(nstep+n+1) - 1
            do k = k1 , k2
               do j = 1 , nsh1
                  c(nstep2+ii+(j-1)*nbasl) = cc(nstep1+n+(j-1)*nbc1)
     +               *cont(k)
               enddo
               ii = ii + 1
            enddo
         enddo
         nstep = nstep + nbc1
         nstep1 = nstep1 + nbc1**2
         nstep2 = nstep2 + nbasl**2
      enddo
      return
      end
