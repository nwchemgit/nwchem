      subroutine two_index_transf( g_a, g_lhs, g_rhs, g_tmp, g_b )
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, g_lhs, g_rhs, g_tmp, g_b
c
      integer ncola, nrowa, ncolb, nrowb, ncoll, nrowl, ncolr, nrowr
      integer ncolt, nrowt
      integer gtype
c     
c     2-index square matrix transform:   B = tr(LHS).A.RHS
c
c     The tmp array should be dimensioned at least (nrowa,ncolb)
c     
      call ga_inquire(g_a,   gtype, nrowa, ncola)
      call ga_inquire(g_b,   gtype, nrowb, ncolb)
      call ga_inquire(g_lhs, gtype, nrowl, ncoll)
      call ga_inquire(g_rhs, gtype, nrowr, ncolr)
      call ga_inquire(g_tmp, gtype, nrowt, ncolt)
c
      if (nrowl.ne.nrowa .or. ncoll.ne.nrowb .or.
     $    nrowr.ne.ncola .or. ncolr.ne.ncolb .or.
     $    nrowt.lt.nrowa .or. ncolt.lt.ncolb)
     $     call errquit('two_index_transf: invalid arguments',0, GA_ERR)
c
      call ga_dgemm('n', 'n', nrowa, ncolb, ncola, 1.d0, g_a, g_rhs, 
     $     0.d0, g_tmp)
      call ga_dgemm('t', 'n', nrowb, ncolb, nrowa, 1.d0, g_lhs, g_tmp, 
     $     0.d0, g_b )
      
      return
      end
