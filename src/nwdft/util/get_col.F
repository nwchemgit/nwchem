      subroutine get_col(g_a,n,i,vec)
C$Id$
      implicit none
c
c     wrapper routine to get a single column out
c     of a Global Array. 
c
      integer g_a,n,i
      double precision vec(*)
      call ga_get(g_a,1,n,i,i,vec,n)
      return
      end


      subroutine put_col(g_a,n,i,vec)
      implicit none
c
c     wrapper routine to put a single column into
c     a Global Array.
c
      integer g_a,n,i
      double precision vec(*)
      call ga_put(g_a,1,n,i,i,vec,n)
      return
      end
