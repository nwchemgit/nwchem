      subroutine get_col(g_a,n,i,vec)
C$Id: get_col.F,v 1.2 1995-02-02 23:19:58 d3g681 Exp $
      implicit none
      integer g_a,n,i
      double precision vec(n)
      call ga_get(g_a,
     *     1,n,i,i,vec,1)
      return
      end
      subroutine put_col(g_a,n,i,vec)
      implicit none
      integer g_a,n,i
      double precision vec(n)
      call ga_put(g_a,
     *     1,n,i,i,vec,1)
      return
      end
