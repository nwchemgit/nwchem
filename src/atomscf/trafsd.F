      subroutine trafsd(nsym,nbas,ndim,a,nbc,cont,nstrt,ffc)
C$Id$
      implicit none
c.......................................................................
c     transform matrix a
c     from prim. to cont. basis.
c.......................................................................
      real *8 ffc(*), a(*), cont(*)
      integer nbas(*), nbc(*), nstrt(*)
      integer nsym, ndim
c
      real *8 sumc
      integer i, j, m, n, l
      integer mnstep, ijstep, k, ijbas, mnbas, nbc1, ij, jj
      integer m1, m2, n0, n2
c
      do i = 1 , ndim
         ffc(i) = a(i)
      enddo
      mnstep = 0
      ijstep = 0
      k = 1
      ijbas = 0
      mnbas = 0
      do  l = 1 , nsym
         nbc1 = nbc(l)
         do m = 1 , nbc1
            do n = 1 , m
               m1 = nstrt(mnbas+m)
               n0 = nstrt(mnbas+n)
               m2 = nstrt(mnbas+m+1) - 1
               n2 = nstrt(mnbas+n+1) - 1
               sumc = 0.0d0
               do i = m1 , m2
                  ij = ijstep + (i-ijbas-1)*(i-ijbas)/2
                  if (m.eq.n) then
                     n2 = i
                  end if
                  do j = n0 , n2
                     jj = j - ijbas
                     if (m.ne.n) then
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     else if (i.ne.j) then
                        sumc = sumc + 2.0d0*cont(i)*cont(j)*ffc(ij+jj)
                     else
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     end if
                  enddo
               enddo
               a(k) = sumc
               k = k + 1
            enddo
         enddo
         mnstep = mnstep + nbc(l)*(nbc(l)+1)/2
         ijstep = ijstep + nbas(l)*(nbas(l)+1)/2
         mnbas = mnbas + nbc(l)
         ijbas = ijbas + nbas(l)
      enddo
      return
      end
