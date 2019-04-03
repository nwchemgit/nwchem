c
c     Performs the matrix-vector products:  g_Ax = (A+B)g_x
c     LHS of Eq. 18 (Furche & Ahlrichs)
c     HvD 10/2007
c
      subroutine tddft_grad_solve_product(acc,g_x,g_Ax)
c
c     In the context of the linear systems solver cphf_solve2.
c     This means that g_x and g_Ax are stored a 2D matrices in the
c     format
c
c     g(1:nocca*nvirta+noccb*nvirtb,1:nroots)
c
c     This requires an extra transformation to get back to the 3D 
c     representation that is used in the rest of the code. The whole 
c     procedure is sketched out in subroutine tddft_grad_solve.
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "tddft_grad_solve_common.fh"
c
c     Input:
c
      integer g_x          ! guess vectors stored in 2D MO
      double precision acc ! a tolerance on the matrix-vector product
c
c     Output:
c
      integer g_Ax         ! the matrix vector products stored as 2D MO
c
c     Local:
c
      integer g_x_3d_mo(2)  ! handle for work space array
      integer g_x_3d_ao     ! handle for work space array
      integer g_Ax_3d_ao    ! handle for work space array
      integer g_Ax_3d_mo(2) ! handle for work space array
      integer g_ambt        ! handle for work space array
c
      integer idim(3)      ! dimensions
      integer ichnk(3)     ! chunking
      integer alo(3)       ! lower chunk limits on A
      integer ahi(3)       ! upper chunk limits on A
c     integer blo(3)       ! lower chunk limits on B
c     integer bhi(3)       ! upper chunk limits on B
c     integer clo(3)       ! lower chunk limits on C
c     integer chi(3)       ! upper chunk limits on C
      integer klo          ! dimensions for mat_to_vec
      integer col          ! dimensions for mat_to_vec
c
      integer ip           ! counter over spin components
      integer ir           ! counter over roots
c
      character*32 pname
      integer iwhich
c
      pname = "tddft_grad_solve_product: "
      iwhich = 0 ! call to tddft_nga_cont()
c
c     1. g_x 2D MO -> g_x_3d_mo
c
      call tddft_grad_create_r(tdg_ipol,tdg_naoc,tdg_nav,tdg_nroots,
     +                         g_x_3d_mo)
      do ip = 1, tdg_ipol
        do ir = 1, tdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = tdg_naoc(ip)
          alo(3) = 1
          ahi(3) = tdg_nav(ip)
          klo = (ip-1)*(tdg_naoc(1)*tdg_nav(1))+1
          col = ir
          call tddft_grad_vec_to_mat(g_x_3d_mo(ip),alo,ahi,g_x,klo,col)
        enddo
      enddo
c
c     2. g_x_3d_mo -> g_x_3d_ao
c
      idim(1) = tdg_nroots*tdg_ipol
      idim(2) = tdg_nao
      idim(3) = tdg_nao
      ichnk(1) = tdg_nroots*tdg_ipol
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'g_x_3d_ao',ichnk,g_x_3d_ao))
     +    call errquit(pname//'failed to create g_x_3d_ao',0,GA_ERR)
      call ga_zero(g_x_3d_ao)
      call tddft_grad_trans_mo2ao(tdg_ipol,tdg_nao,tdg_nfc,tdg_naoc,
     +     tdg_nocc,tdg_nav,tdg_nfv,tdg_nroots,
     +     1.0d0,0.0d0,"ib",tdg_g_mo,g_x_3d_mo,"ib",g_x_3d_ao)
c
c     3. compute (A+B)g_x currently we have to compute (A-B)g_x as well
c        although we do not want it.
c
      if (.not.nga_create(mt_dbl,3,idim,'vec (A+B)X',ichnk,g_Ax_3d_ao))
     +    call errquit(pname//'failed to create g_Ax_3d_ao',0,GA_ERR)
c
      if (.not.tdg_tda) then
        if (.not.nga_create(mt_dbl,3,idim,'vec (A-B)X',ichnk,g_ambt))
     +      call errquit(pname//'failed to create g_ambt',0,GA_ERR)
      endif
      call ga_zero(g_Ax_3d_ao)
      call tddft_nga_cont(tdg_rtdb,tdg_geom,tdg_bfao,g_x_3d_ao,
     +     g_Ax_3d_ao,g_ambt,tdg_nao,tdg_ipol,tdg_tol2e,tdg_tda,
     +     tdg_oskel,tdg_kfac,tdg_lhashf,tdg_triplet,tdg_nroots,iwhich)
      if (.not.tdg_tda) then
        if (.not.ga_destroy(g_ambt))
     +      call errquit(pname//'failed to destroy g_ambt',0,GA_ERR)
      endif
c
c     4. Destroy g_x_3d_ao
c
      if (.not.ga_destroy(g_x_3d_ao)) 
     +   call errquit(pname//'failed to destroy g_x_3d_ao',0,GA_ERR)
c
c     5. Transform g_Ax_3d_ao -> g_Ax_3d_mo
c
      call tddft_grad_create_r(tdg_ipol,tdg_naoc,tdg_nav,tdg_nroots,
     +     g_Ax_3d_mo)
      call ga_zero(g_Ax_3d_mo)
      call tddft_grad_trans_ao2mo(tdg_ipol,tdg_nao,tdg_nfc,tdg_naoc,
     +     tdg_nocc,tdg_nav,tdg_nfv,tdg_nroots,1.0d0,0.0d0,"ib",
     +     tdg_g_mo,g_Ax_3d_ao,g_Ax_3d_mo,"ib")
c
c     6. Destroy g_Ax_3d_ao
c
      if (.not.ga_destroy(g_Ax_3d_ao)) 
     +   call errquit(pname//'failed to destroy g_Ax_3d_ao',0,GA_ERR)
c
c     7. Add orbital energy difference term
c
      call tddft_grad_epsilon(g_x_3d_mo,g_Ax_3d_mo,tdg_ipol,tdg_nroots,
     +     tdg_nao,tdg_nocc,tdg_nfc,dbl_mb(tdg_e_mo))
c
c     8. destroy g_x_3d_mo
c
      call tddft_grad_destroy_r(tdg_ipol,g_x_3d_mo)
c
c     9. Transform g_Ax_3d_mo -> g_Ax 2D MO
c
      do ip = 1, tdg_ipol
        do ir = 1, tdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = tdg_naoc(ip)
          alo(3) = 1
          ahi(3) = tdg_nav(ip)
          klo = (ip-1)*(tdg_naoc(1)*tdg_nav(1))+1
          col = ir
          call tddft_grad_mat_to_vec(g_Ax_3d_mo(ip),alo,ahi,g_Ax,
     +                               klo,col,1.0d0,"=")
        enddo
      enddo
c
      call tddft_grad_destroy_r(tdg_ipol,g_Ax_3d_mo)
c
      end
c $Id$
