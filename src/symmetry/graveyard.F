      subroutine sym_bas_munge_op(r, maxf, ang_max)
C$Id$
      implicit none
      integer maxf, ang_max
      double precision r(maxf, maxf, 0:ang_max)
c
      integer lmax              ! For basis normalization
      parameter (lmax=7)
      double precision dfac(0:lmax), scale
      integer i,j,k,ijk,l,mno,maxl
c
c     dfac(i) = sqrt((2*i-1)!!)
c
      dfac(0) = 1.0d0
      do i = 1,ang_max
         dfac(i) = dble(2*i-1)*dfac(i-1)
      enddo
      do i = 1, ang_max
         dfac(i)=sqrt(dfac(i))
      enddo
c
      do L = 0, ang_max
         ijk = 0
         maxl = ((L+1)*(L+2))/2
         do i = L,0,-1
            do j = L-i,0,-1
               k = L-i-j
               scale = dfac(i)*dfac(j)*dfac(k)
               ijk = ijk + 1
               do mno = 1, maxl
                  r(ijk,mno,L) = r(ijk,mno,L)*scale
               enddo
               do mno = 1, maxl
                  r(mno,ijk,L) = r(mno,ijk,L)/scale
               enddo
            enddo
         enddo
      enddo
c
      end
