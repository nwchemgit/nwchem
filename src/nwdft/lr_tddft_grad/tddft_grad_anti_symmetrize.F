      subroutine tddft_grad_anti_symmetrize(g_mat)
      implicit none
c
c     Take the matrix in g_mat and add the negative of its transpose
c     to it.
c
c     It seems that there is a routine ga_antisymmetrize already in 
c     the utils directory. Could that replace this one?
c
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer g_mat
      integer g_tmp, ndim, type
      integer ihi(2), ichnk(2)
c
      ichnk(1) = -1
      ichnk(2) = -1
c
      ndim = ga_ndim(g_mat)
      call nga_inquire(g_mat,type,ndim,ihi)
      if (.not.nga_create(mt_dbl,ndim,ihi,"tmp",ichnk,g_tmp)) call
     +  errquit("tddft_grad_anti_symmetrize: failed to create tmp",
     +          0, GA_ERR)
      call ga_transpose(g_mat,g_tmp)
      call ga_add(-0.5d0,g_tmp,0.5d0,g_mat,g_mat)
      if (.not.ga_destroy(g_tmp)) call errquit
     +  ("tddft_grad_anti_symmetrize: failed to destroy tmp",
     +   0, GA_ERR)
      end
      
c $Id$
