


      subroutine pdfded(k,kdim,g_dens,dhelp,factor,dmult,lm,nbci,iiloc)
C$Id$
      implicit none
c.......................................................................
c
c     routine for distributing atomic densities to the molecular
c     density matrix. note that the number of orbitals differ in
c     the molecular and atomic case for d and f orbitals. transformat-
c     ion matrices are provided in dmult(*,*). to work out these tables
c     real atomic orbitals are needed. for f functions these are:
c              sqrt(1/60)*(2zzz-3zyy-3zxx)
c              sqrt(1/40)*(4zzy-yyy-xxy)
c              sqrt(1/40)*(4zzx-xxx-xyy)
c                    xyz
c              sqrt(1/4)*(xxz-yyz)
c              sqrt(1/24)*(3xxy-yyy)
c              sqrt(1/24)*(3xyy-xxx)
c     normalization of primitives is given by (xyz:xyz)=1, (xxy:xxy)=3
c     (xxx:xxx)=15.
c
c.......................................................................
c..
      integer g_dens
      integer k, kdim, lm, nbci
      real *8 dhelp(*),dmult(kdim,kdim)
      real *8 factor
      integer iiloc(nbci)
c
      integer l, na, m, nb
      integer ll,mm
      integer lmsave, nbrang
      real *8 delem,xx
c     
c..
c..
      do l = 1 , nbci
         lmsave = lm
         do na = 1 , kdim
            lm = lmsave
            do m = 1 , l
               lm = lm + 1
               delem = dhelp(lm)*factor
               nbrang = kdim
               if (m.eq.l) nbrang = na
               do nb = 1 , nbrang
                  xx = delem*dmult(na,nb)
                  ll = iiloc(l)+nb-1
                  mm = iiloc(m)+na-1
                  call ga_put(g_dens,ll,ll,mm,mm,xx,1)
                  call ga_put(g_dens,mm,mm,ll,ll,xx,1)
               enddo
            enddo
         enddo
      enddo
c..
      return
      end
