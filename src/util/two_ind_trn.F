      subroutine two_index_transf( g_a, g_lhs, g_rhs, g_tmp, g_b )
C     $Id: two_ind_trn.F,v 1.1 1995-12-13 01:37:42 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, g_lhs, g_rhs, g_tmp, g_b
      integer n, row, col, gtype
c     
c     2-index square matrix transform:   B = tr(LHS).A.RHS
c     
      call ga_inquire(g_a, gtype, row, col)
      if (row.ne.col) call errquit(
     $     'two_index_transf: invalid arguments',0)
      n = row
      call ga_inquire(g_lhs, gtype, row, col)
      if ((row.ne.col).or.(row.ne.n)) call errquit(
     $     'two_index_transf: invalid arguments',0)
      call ga_inquire(g_rhs, gtype, row, col)
      if ((row.ne.col).or.(row.ne.n)) call errquit(
     $     'two_index_transf: invalid arguments',0)
      call ga_inquire(g_tmp, gtype, row, col)
      if ((row.ne.col).or.(row.ne.n)) call errquit(
     $     'two_index_transf: invalid arguments',0)
      call ga_inquire(g_b, gtype, row, col)
      if ((row.ne.col).or.(row.ne.n)) call errquit(
     $     'two_index_transf: invalid arguments',0)

      call ga_dgemm('n', 'n', n, n, n, 1.d0, g_a, g_rhs, 0.d0, g_tmp)
      call ga_dgemm('t', 'n', n, n, n, 1.d0, g_lhs, g_tmp, 0.d0, g_b )
      
      return
      end
