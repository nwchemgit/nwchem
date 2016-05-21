C
C     Takes complex MO density matrix and stores the corresponding real
C     part of the movecs onto disk as file "fname".
C     
      subroutine rt_tddft_zdens_movecs (params, g_zfock_mo,
     $     g_zdens_ao, g_zdens_mo, fname)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zfock_mo 
      integer, intent(in)           :: g_zdens_ao
      integer, intent(in)           :: g_zdens_mo !complex MO density matrix (either CS or total, ie, alpha+beta for OS)
      character(len=*), intent(in)  :: fname
      

C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_zdens_tot_movecs: "


C     == Variables ==
      integer :: g_zevecs, g_zevecs_ao
      integer :: g_movecs_re, g_movecs_im
      integer :: g_ztmp1, g_ztmp2, g_ztmp3, g_ztmp4
      integer :: lvals, ivals
      integer :: lvals_re, ivals_re
      integer :: locc, iocc
      integer :: i
      integer :: me
      double complex :: zocc


C     == External ==
      logical, external :: movecs_write


      me = ga_nodeid()

      
C
C     Allocation
C
      
      if (params%ns_ao .ne. params%ns_mo)
     $     call errquit (pname//
     $     "alloc not working for nao != nmo",0,0)

      if (.not. ga_duplicate (g_zdens_mo, g_zevecs, "zevecs"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zdens_mo, g_zevecs_ao, "zevecs AO"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zevecs_ao, g_ztmp1, "ztmp1"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zevecs_ao, g_ztmp2, "ztmp2"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zevecs_ao, g_ztmp3, "ztmp3"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zevecs_ao, g_ztmp4, "ztmp4"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "movecs re", 0,0, g_movecs_re))
     $     call errquit (pname//"failed to creat g_movecs_re",0,GA_ERR)

      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "movecs re", 0,0, g_movecs_im))
     $     call errquit (pname//"failed to creat g_movecs_im",0,GA_ERR)


      
      if (.not.ma_push_get(mt_dcpl, params%ns_mo,'vals',lvals,ivals))
     $     call errquit(pname//'cannot allocate vals',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_mo,'vals_re',
     $     lvals_re,ivals_re))
     $     call errquit(pname//'cannot allocate vals_re',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_mo,'occ',
     $     locc,iocc))
     $     call errquit(pname//'cannot allocate vals_re',0, MA_ERR)
      

      
C
C     Diagonalize to find complex movecs
C
      call ga_zero (g_zevecs)
      call zmat_diag (g_zfock_mo, g_zevecs,
     $     dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy, .false. = dont check

      


C
C     Convert MO vectors to AO basis
C
      call ga_zero (g_zevecs_ao)

      if (params%ns_ao .ne. params%ns_mo)
     $     call errquit (pname//"not working for nao != nmo",0,0)
      
C     XXX fix dims
      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     z1, params%g_zX, g_zevecs, z0, g_zevecs_ao)

      
C
C     Dump real part pto disk.
C
c$$$      do i = 1, params%ns_ao
c$$$         if (i .le. 1) then
c$$$            dbl_mb(iocc + i - 1) = 2d0
c$$$         else
c$$$            dbl_mb(iocc + i - 1) = 0d0
c$$$         endif
c$$$
c$$$         dbl_mb(ivals_re + i - 1) = dble (dcpl_mb(ivals + i - 1))
c$$$      enddo
c$$$
c$$$      call convert_z2d (g_zevecs_ao, g_movecs_re, g_movecs_im)
c$$$      
c$$$      if (.not.movecs_write (params%rtdb, ao_bas_han, fname, "dft",
c$$$     $     title, params%ns_ao, ipol, params%ns_mo,
c$$$     $     dbl_mb(iocc), params%ns_ao, 
c$$$     $     dbl_mb(ivals_re), params%ns_ao, g_movecs_re))
c$$$     $     call errquit(pname//"movec_write failed", 0, DISK_ERR)



C
C     Check that we can rebuild the density matrix
C
c$$$      call ga_zero (g_ztmp1)
c$$$      call mat_hermconj (g_zevecs_ao, g_ztmp1)
c$$$
c$$$      
c$$$C     (diag matrix of occupations)
c$$$      call ga_zero (g_ztmp2)
c$$$      do i = 1, params%ns_ao
c$$$         if (i .le. 1) call ga_put (g_ztmp2, i, i, i, i, z2, 1)
c$$$      enddo
c$$$
c$$$C     D C^+
c$$$      call ga_zero (g_ztmp3)
c$$$      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
c$$$     $     z1, g_ztmp2, g_ztmp1, z0, g_ztmp3)  ! <- g_ztmp3 now has D C^+
c$$$
c$$$      
c$$$C     C D C^+
c$$$      call ga_zero (g_ztmp1)
c$$$      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
c$$$     $     z1, g_zevecs_ao, g_ztmp3, z0, g_ztmp1)  ! <- g_ztmp1 now has C D C^+

c$$$      call ga_print (g_ztmp1)


C     (make diag occup matrix)
      call ga_zero (g_ztmp3)
      do i = 1, params%ns_ao
         if (i .le. 1) then
            zocc = z2
         else
            zocc = z0
         endif
         
         call ga_fill_patch (g_ztmp3, 1, params%ns_ao,
     $        i, i, zocc)
      enddo
      
C     (project onto TD complex movecs)
C     Note: only for Hermitian matrices (.true. flag)
      call rt_tddft_project_check (params, g_ztmp3, g_zdens_mo, .true.)
      if (me.eq.0) write (luout,*) pname//"projection check passed"
      call ga_zero (g_ztmp1)
      call rt_tddft_project (params, "R", g_ztmp3,
     $     g_zdens_mo, g_ztmp1, .true.) ! g_ztmp1 now has projected occup mask
      
      call ga_zero (g_ztmp4)
      call ga_elem_multiply (g_ztmp1, g_zevecs_ao, g_ztmp4) ! g_ztmp4 now holds masked movecs in AO basis

      call ga_zero (g_ztmp1)
      call mat_hermconj (g_zevecs_ao, g_ztmp1)

      call ga_zero (g_ztmp2)
      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     z1, g_ztmp4, g_ztmp1, z0, g_ztmp2)  !g_ztmp2 now has AO dens mat


      if (.not. mat_is_same (g_ztmp2, g_zdens_ao, params%tol_zero)) then
         if (me.eq.0)
     $        write (luout, *) "zdens -> movecs conversion FAILED!"
      else
         if (me.eq.0)
     $        write (luout, *) "zdens -> movecs conversion passed"
      endif
      


C
C     Clean up
C      
      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ma_chop_stack (lvals))
     $   call errquit(pname//'cannot chop stack',0, MA_ERR)

      end subroutine
c $Id$
