      subroutine get_col(g_a,n,i,vec)
C$Id$
      implicit none
c
c     Apparently this is a wrapper routine to get a single column out
c     of a Global Array. So logically VEC is intended as an array of
c     N-by-1, and hence the leading dimension fed to GA_GET should be N.
c
      integer g_a,n,i
      double precision vec(n,1)
      call ga_get(g_a,1,n,i,i,vec,n)
      return
      end


      subroutine put_col(g_a,n,i,vec)
      implicit none
c
c     Apparently this is a wrapper routine to put a single column into
c     a Global Array. So logically VEC is intended as an array of
c     N-by-1, and hence the leading dimension fed to GA_PUT should be N.
c
      integer g_a,n,i
      double precision vec(n,1)
      call ga_put(g_a,1,n,i,i,vec,n)
      return
      end
