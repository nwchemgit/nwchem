C
C     Compute dS/dR (gradient of overlap matrix)
C
      subroutine rt_tddft_grads (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
#include "bas.fh"


C     == Inputs ==
      type(rt_params_t), intent(inout) :: params !< struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_grads: "

      
C     == Variables ==
      integer g_grads_x, g_grads_y, g_grads_z
      integer lbuf, l_buf, k_buf
      integer lscratch, l_scr, k_scr
      integer max1e, max2e, mscratch_1e, mscratch_2e
      integer idatom(2)
      
      
      
C
C     Allocation
C
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e)
      max2e = max(max2e,1296*100)
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)

      
      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "dS/dR x", 0, 0, g_grads_x))
     $     call errquit (pname//"created failed",0,GA_ERR)
      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "dS/dR y", 0, 0, g_grads_y))
     $     call errquit (pname//"created failed",0,GA_ERR)
      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "dS/dR z", 0, 0, g_grads_z))
     $     call errquit (pname//"created failed",0,GA_ERR)

      if (.not. ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit('grad:could not allocate buffer',lbuf,
     &     MA_ERR)

      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch', 
     $     l_scr, k_scr))call errquit('grad: scratch alloc failed',
     $     lscratch, MA_ERR)

      
      call ga_zero (g_grads_x)
      call ga_zero (g_grads_y)
      call ga_zero (g_grads_z)


c$$$      call intd_1eov(basis,ish1,basis,ish2,lscr,scr,
c$$$     &     lbuf,H,idatom)
      
      call intd_1eov (ao_bas_han, 1, ao_bas_han, 1,
     $     lscratch, dbl_mb(k_scr), lbuf, dbl_mb(k_buf), idatom)

      call halt ()
      

C
C     Clean up
C
      if (.not. ma_chop_stack (l_buf))
     $     call errquit (pname//"chop failed", 0, GA_ERR)
      
      if (.not. ga_destroy (g_grads_x))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      if (.not. ga_destroy (g_grads_y))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      if (.not. ga_destroy (g_grads_z))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      
      end subroutine
