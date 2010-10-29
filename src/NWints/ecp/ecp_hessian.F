C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_hessian (
     &    xyz_A,zeta_A,coef_A,n_prim_A,n_cont_A,l_A,i_c_A,
     &    xyz_B,zeta_B,coef_B,n_prim_B,n_cont_B,l_B,i_c_B,
     &    xyz_C,zeta_C,coef_C,n_prim_C,n_colc_C,
     &    ind_z,ind_c,n_zeta_C,n_coef_C,
     &    l_C,i_cent_C,n_C,l_ecp_max,
     &    sphcart,csco,lcsco,
     &    ecp_hess,n_int,n_blk,n_cent,
     &    DryRun,scr,lscr,ibug)
*                                                                      *
*   Calculate ecp gradient integrals.                                  *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   xyz_A (inp) - coordinates of centre A                              *
*   zeta_A (inp) - exponents of primitive gaussians on centre A        *
*   coef_A (inp) - contraction coefficients on centre A                *
*   n_prim_A (inp) - number of primitive gaussians on centre A         *
*   n_cont_A (inp) - number of contracted functions on centre A        *
*   l_A (inp) - angular momentum of functions on centre A              *
*   i_c_A (inp) - index of centre A                                    *
*   xyz_B (inp) - coordinates of centre B                              *
*   zeta_B (inp) - exponents of primitive gaussians on centre B        *
*   coef_B (inp) - contraction coefficients on centre B                *
*   n_prim_B (inp) - number of primitive gaussians on centre B         *
*   n_cont_B (inp) - number of contracted functions on centre B        *
*   l_B (inp) - angular momentum of functions on centre B              *
*   i_c_B (inp) - index of centre B                                    *
*   xyz_C (inp) - coordinates of ECP centres C                         *
*   zeta_C (inp) - array of exponents of primitive gaussians on all    *
*                  centres C. These are stored in an array of single   *
*                  dimension, i.e. packed.                             *
*   coef_C (inp) - array of contraction coefficients on all centres C  *
*   n_prim_C (inp) - array of number of primitive gaussians for each   *
*                    power of r, l value and ECP centre. The highest   *
*                    l value is for the local part, thus the second    *
*                    dimension is l_ecp_max+2 (or 0:l_ecp_max+1)       *
*   n_colc_C (inp) - array of number of coefficients for each l value  *
*                    and ECP centre. This is n_prim_C summed over the  *
*                    first dimension.                                  *
*   ind_z (inp) - array of addresses of first exponent for each l      *
*                 value and ECP centre.                                *
*   ind_c (inp) - array of addresses of first coefficient for each l   *
*                 value and ECP centre.                                *
*   n_zeta_C (inp) - total number of ECP exponents.                    *
*   n_coef_C (inp) - total number of ECP coefficients.                 *
*   l_C (inp) - maximum angular momentum of projectors on centres C    *
*   i_cent_C - indices of ECP centres C                                *
*   n_C (inp) - number of ECP centres C                                *
*   l_ecp_max (inp) - maximum angular momentum of any projector on any *
*                 ECP centre                                           *
*   sphcart - 0 for cartesian integrals, 1 for spherical integrals     *
*   ecp_hess (out) - hessian integrals over ECPs                       *
*   n_int (inp) - number of ECP hessian integrals for a given centre   *
*                 and cartesian derivative pair. Should be equal to    *
*                 NCA*NCB*[(La+1)*(La+2)/2]*[(Lb+1)*(Lb+2)/2]          *
*   n_blk (inp) - 1 for scalar only, 3 for s-o only, 4 for both        *
*   n_cent (inp) - total number of centres in the molecule             *
*   DryRun (inp) - logical for dry run. If true, routine only returns  *
*                  maximum scratch space needed, if false, integrals   *
*                  are returned.                                       *
*   scr (scr) - scratch array for work space                           *
*   lscr (i/o) - length of scratch array. Value returned if DryRun is  *
*                true, used as dimension if false.                     *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
#include "ecp_consts.fh"
#include "util.fh"
#include "errquit.fh"
      integer i,j,k,l,m,
     &    n_prim_A,n_cont_A,l_A,i_c_A,
     &    n_prim_B,n_cont_B,l_B,i_c_B,
     &    n_zeta_C,n_coef_C,n_C,i_c_C,l_ecp_max,
     &    n_int,n_blk,n_cent,sphcart,lscr,ibug
      integer n_prim_C(0:4,-1:l_ecp_max,n_C*2),
     &    n_colc_C(-1:l_ecp_max,n_C*2),
     &    ind_z(-1:l_ecp_max,n_C*2),ind_c(-1:l_ecp_max,n_C*2),
     &    l_C(n_C),i_cent_C(n_C)
      integer i_xca,i_xxca,i_xcb,i_xxcb,i_lp,i_lm,i_grdp,i_grdm,
     &    n_cartp1_A,n_cartp1_B,n_cartp2_A,n_cartp2_B,
     &    n_cartm1_A,n_cartm1_B,n_cartm2_A,n_cartm2_B,
     &    n_allp1_A,n_allp1_B,n_allp2_A,n_allp2_B,
     &    n_allm1_A,n_allm1_B,n_allm2_A,n_allm2_B,
     &    n_intpa,n_intpb,n_intma,n_intmb,
     &    n_intpaa,n_intmaa,n_intpbb,n_intmbb,
     &    n_intppab,n_intpmab,n_intmpab,n_intmmab
      integer i_scr,i_s,i_zeta,i_coef,
     &    mem,max_type1,max_type2,l_scr,
     &    n_cart_a,n_cart_b,n_cont_ab,
     &    n_all_a,n_all_b,n_ab,n_abblk,n_ab9blk
      integer lcsco, i_loc, ibugt
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision 
     &    xyz_A(3),zeta_A(n_prim_A),coef_A(n_prim_A,n_cont_A),
     &    xyz_B(3),zeta_B(n_prim_B),coef_B(n_prim_B,n_cont_B),
     &    xyz_C(3,n_C),zeta_C(n_zeta_C),coef_C(n_zeta_C),
     &    scr(lscr),ecp_hess(n_int*n_blk,3,3,n_cent*(n_cent-1)/2+n_cent)
      double precision 
     &    X_AC,Y_AC,Z_AC,R_AC,X_BC,Y_BC,Z_BC,R_BC,
     &    tol
      double precision csco(lcsco)
      character*4 dbyd(3)
      character*16 printhead
!      data tol/1.0d-14/
      data tol/1.0d-16/
      data dbyd/'d/dx','d/dy','d/dz'/
      data printhead/'          U_ab  '/
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_hessian ...'
      if (debug_addresses) then
        write(LuOut,*)' lscr in ecp_hessian:',lscr
        write (LuOut,*) 'n_prim_A,n_cont_A,l_A',n_prim_A,n_cont_A,l_A
        write (LuOut,*) 'n_prim_B,n_cont_B,l_B',n_prim_B,n_cont_B,l_B
        write (LuOut,*) 'l_ecp_max,n_c,n_cent',l_ecp_max,n_c,n_cent
      end if
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_all_b = n_cart_b*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_cont_ab = n_cont_a*n_cont_b
      n_ab = n_all_a*n_all_b
      n_abblk = n_ab*n_blk
      n_ab9blk = n_abblk*9
      if (debug_addresses) write (LuOut,*) 
     &    'n_cart_a,n_cart_b,n_cont_ab,n_ab',
     &    n_cart_a,n_cart_b,n_cont_ab,n_ab
      if (n_int .lt. n_abblk) call errquit (
     &    'Mismatch of integral count in ecp_hessian',99, BASIS_ERR)
      if (sphcart .ne. 0) call errquit(
     &    'Do your own spherical transformation, lazy bum!',99,
     &       BASIS_ERR)
      if (DryRun) then
        max_type1 = 0
        max_type2 = 0
      else
        i_loc = n_cent * (n_cent-1) / 2 + n_cent
        call dfill (n_int*n_blk*9*i_loc,zero,ecp_hess,1)
      end if
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      n_cartp1_A = n_cart_A+l_A+2
      n_cartp1_B = n_cart_B+l_B+2
      n_cartp2_A = n_cartp1_A+l_A+3
      n_cartp2_B = n_cartp1_B+l_B+3
      n_cartm1_A = n_cart_A-l_A-1
      n_cartm1_B = n_cart_B-l_B-1
      n_cartm2_A = n_cartm1_A-l_A
      n_cartm2_B = n_cartm1_B-l_B
c     n_cartm2_A = n_cartm1_A-l_A-2
c     n_cartm2_B = n_cartm1_B-l_B-2
*
      n_allp1_A = n_cartp1_A*n_cont_A
      n_allp1_B = n_cartp1_B*n_cont_B
      n_allp2_A = n_cartp2_A*n_cont_A
      n_allp2_B = n_cartp2_B*n_cont_B
      n_allm1_A = n_cartm1_A*n_cont_A
      n_allm1_B = n_cartm1_B*n_cont_B
      n_allm2_A = n_cartm2_A*n_cont_A
      n_allm2_B = n_cartm2_B*n_cont_B
*
      n_intpa = n_allp1_A*n_all_B
      n_intpb = n_all_A*n_allp1_B
      n_intma = n_allm1_A*n_all_B
      n_intmb = n_all_A*n_allm1_B
      n_intpaa = n_allp2_A*n_all_B
      n_intmaa = n_allm2_A*n_all_B
      n_intpbb = n_all_A*n_allp2_B
      n_intmbb = n_all_A*n_allm2_B
      n_intppab = n_allp1_A*n_allp1_B
      n_intpmab = n_allp1_A*n_allm1_B
      n_intmpab = n_allm1_A*n_allp1_B
      n_intmmab = n_allm1_A*n_allm1_B
*
      i_xca = 1
      i_xxca = i_xca+n_prim_A*n_cont_A
      i_xcb = i_xxca+n_prim_A*n_cont_A
      i_xxcb = i_xcb+n_prim_B*n_cont_B
      i_lp = i_xxcb+n_prim_B*n_cont_B
      i_lm = i_lp+max(n_intpaa,n_intpbb,n_intppab,n_intmbb,
     &                n_intmaa,n_intpmab)*n_blk
      i_grdp = i_lm+max(n_ab,n_intmmab,n_intmpab)*n_blk
      i_grdm = i_grdp+max(n_intpa,n_intpb)*3*n_blk
      i_scr = i_grdm+max(n_intma,n_intmb)*3*n_blk
      if (debug_addresses) then
        write (LuOut,*) 'i_xca,i_xcb',i_xca,i_xcb
        write (LuOut,*) 'i_lp,i_lm',i_lp,i_lm
        write (LuOut,*) 'i_scr',i_scr
      end if
*
*   Set up coefficients multiplied by exponents
*
      if (.not.DryRun) then
        l_scr = lscr-i_scr+1
        if (l_scr .lt. 0) call errquit (
     &      'Insufficient scratch memory in ecp_hessian',99,
     &       MEM_ERR)
        k = i_xca-1
        m = n_prim_A*n_cont_A
        do j = 1,n_cont_a
          do i = 1,n_prim_A
            scr(k+i) = zeta_A(i)*coef_A(i,j)
            scr(k+m+i) = zeta_A(i)*scr(k+i)
          end do
          k = k+n_prim_A
        end do
        k = i_xcb-1
        m = n_prim_B*n_cont_B
        do j = 1,n_cont_B
          do i = 1,n_prim_B
            scr(k+i) = zeta_B(i)*coef_B(i,j)
            scr(k+m+i) = zeta_B(i)*scr(k+i)
          end do
          k = k+n_prim_B
        end do
      end if
*
*   Loop over ECP centres
*
      if (debug_arrays) call ecp_matpr (xyz_C,1,3,1,n_C,1,3,1,n_C,
     &    'coordinate array','E',81,5)
      if (debug_gen) write (LuOut,*) 'Number of ECP centers =',n_C
      do i = 1,n_C
        l = l_C(i)
        i_c_C = i_cent_C(i)
        if (debug_gen) write (LuOut,*) 'ECP center',i
        if (debug_gen) write (LuOut,*) '   Maximum angular momentum',l
*
*     Set up relative cartesian coordinates
*
        X_AC = xyz_C(1,i)-xyz_A(1)
        Y_AC = xyz_C(2,i)-xyz_A(2)
        Z_AC = xyz_C(3,i)-xyz_A(3)
        R_AC = sqrt(X_AC**2+Y_AC**2+Z_AC**2)
        X_BC = xyz_C(1,i)-xyz_B(1)
        Y_BC = xyz_C(2,i)-xyz_B(2)
        Z_BC = xyz_C(3,i)-xyz_B(3)
        R_BC = sqrt(X_BC**2+Y_BC**2+Z_BC**2)
        if (debug_arrays) then
          write (LuOut,'(3x,A,3F10.6)') 'Relative coords of center A:',
     &        X_AC,Y_AC,Z_AC
          write (LuOut,'(3x,A,3F10.6)') 'Relative coords of center B:',
     &        X_BC,Y_BC,Z_BC
          write (LuOut,'(3x,A,3F10.6)') 'Distance to center A and B:',
     &        R_AC,R_BC
        end if
*
*     Pointers for local part
*
        i_zeta = ind_z(-1,i)
        i_coef = ind_c(-1,i)
*
C        if (debug_arrays) then
C          write (LuOut,*) 'Parameters of ECP'
C          call ecp_matpi (n_prim_C,0,4,-1,l_ecp_max,0,2,0,l_C,
C     &        'n_prim_C',81,15)
C          call ecp_matpi (ind_C,-1,l_ecp_max,1,1,0,l_C,1,1,
C     &        'ind_C',81,15)
C          call ecp_matpi (n_colc_C,-1,l_ecp_max,1,1,0,l_C,1,1,
C     &        'n_colc_C',81,15)
C          call ecp_matpr (zeta_C,1,n_zeta_C,1,1,i_zeta,i_zeta+n-1,1,1,
C     &        'zeta_C','E',81,5)
C          call ecp_matpr (coef_C,1,n_zeta_C,1,1,i_zeta,i_zeta+n-1,1,1,
C     &        'coef_C','E',81,5)
C        end if
*
*       Contributions from second derivative wrt centre B
*       =================================================
*
        if (i_c_C .ne. i_c_B) then
          if (.not.DryRun) then
            call dfill (n_intpbb*n_blk,zero,scr(i_lp),1)
            call dfill (n_ab*n_blk,zero,scr(i_lm),1)
          end if
*
*         Calculate integrals for l_B+2 with coeffs * zeta**2
*
          if (i_c_C .eq. i_c_A) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+2,
     &          n_prim_B,n_cont_B,scr(i_xxcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpbb,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            call ecp_int1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+2,
     &          n_prim_B,n_cont_B,scr(i_xxcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpbb,n_blk,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+2,
     &          n_prim_B,n_cont_B,scr(i_xxcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpbb,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+2,
     &          n_prim_B,n_cont_B,scr(i_xxcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpbb,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Calculate integrals for l_B with coeffs * zeta
*
          if (i_c_C .eq. i_c_A) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            call ecp_int1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,n_blk,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Construct gradient for l_b+1 with coeffs * zeta.
*
          if (.not.DryRun) then
            call dfill (n_intpb*3*n_blk,zero,scr(i_grdp),1)
            call ecp_grad_b (l_B+1,n_cartp2_b,n_cartp1_b,n_cart_b,
     &          n_cont_B,scr(i_lp),scr(i_lm),scr(i_grdp),n_all_a*n_blk,
     &          ibug/10)
            call dfill (n_intmbb*n_blk,zero,scr(i_lp),1)
          end if
*
*         Calculate integrals for l_B-2
*
          if (l_B .gt. 1) then
            if (i_c_C .eq. i_c_A) then
              if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-2,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_c(i_zeta),coef_c(i_coef),0,4,
     &            tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmbb,.false.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
              call ecp_int1 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-2,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmbb,n_blk,.false.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
            else
              if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-2,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_C(i_zeta),coef_C(i_coef),0,4,
     &            tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmbb,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
              call ecp_int2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-2,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmbb,n_blk,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
            end if
          end if
*
*         Construct gradient for l_B-1
*
          if (.not.DryRun) then
            if (l_B .gt. 0) then
              call dfill (n_intmb*3*n_blk,zero,scr(i_grdm),1)
              call ecp_grad_b (l_B-1,n_cart_b,n_cartm1_b,n_cartm2_b,
     &            n_cont_B,scr(i_lm),scr(i_lp),scr(i_grdm),
     &            n_all_a*n_blk,ibug/10)
            end if
*
*         Construct Hessian and add to final array
*
            call dfill (n_ab9blk,zero,scr(i_scr),1)
            call ecp_grad_b (l_B,n_cartp1_b,n_cart_b,n_cartm1_b,
     &          n_cont_B,scr(i_grdp),scr(i_grdm),scr(i_scr),
     &          n_all_a*3*n_blk,ibug/10)
            i_s = i_scr
            do j = 1,3
              do k = 1,3
                i_loc = i_c_B * (i_c_B-1) / 2 + i_c_B
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,k,j,i_loc),1)
                i_loc = i_c_C * (i_c_C-1) / 2 + i_c_C
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,k,j,i_loc),1)
                if (i_c_B .lt. i_c_C) then
                   i_loc = i_c_C * (i_c_C-1) / 2 + i_c_B
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                else
                   i_loc = i_c_B * (i_c_B-1) / 2 + i_c_C
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                endif
                i_s = i_s+n_abblk
              end do
            end do
          end if
        end if
*
*       Contributions from second derivative wrt centre A
*       =================================================
*
        if (i_c_C .ne. i_c_A) then
          if (.not.DryRun) then
            call dfill (n_intpaa*n_blk,zero,scr(i_lp),1)
            call dfill (n_ab*n_blk,zero,scr(i_lm),1)
          end if
*
*         Calculate integrals for l_A+2 with coeffs * zeta**2
*
          if (i_c_C .eq. i_c_B) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+2,
     &          n_prim_A,n_cont_A,scr(i_xxca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpaa,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            call ecp_int1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+2,
     &          n_prim_A,n_cont_A,scr(i_xxca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpaa,n_blk,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+2,
     &          n_prim_A,n_cont_A,scr(i_xxca),zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpaa,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+2,
     &          n_prim_A,n_cont_A,scr(i_xxca),zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpaa,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Calculate integrals for l_A with coeffs * zeta
*
          if (i_c_C .eq. i_c_B) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            call ecp_int1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,n_blk,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_ab,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Construct gradient for l_A+1 with coeffs * zeta.
*
          if (.not.DryRun) then
            call dfill (n_intpa*3*n_blk,zero,scr(i_grdp),1)
            call ecp_grad_a (l_A+1,n_cartp2_a,n_cartp1_a,n_cart_a,
     &          n_cont_A*n_blk,scr(i_lp),scr(i_lm),scr(i_grdp),n_all_b,
     &          ibug/10)
            call dfill (n_intmaa*n_blk,zero,scr(i_lp),1)
          end if
*
*         Calculate integrals for l_A-2
*
          if (l_A .gt. 1) then
            if (i_c_C .eq. i_c_B) then
              if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-2,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_c(i_zeta),coef_c(i_coef),0,4,
     &            tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmaa,.true.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
              call ecp_int1 (mem,DryRun,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-2,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmaa,n_blk,.true.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
            else
              if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-2,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_C(i_zeta),coef_C(i_coef),0,4,
     &            tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmaa,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
              call ecp_int2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-2,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lp),n_intmaa,n_blk,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
            end if
          end if
*
*         Construct gradient for l_A-1
*
          if (.not.DryRun) then
            if (l_A .gt. 0) then
              call dfill (n_intma*3*n_blk,zero,scr(i_grdm),1)
              call ecp_grad_a (l_A-1,n_cart_a,n_cartm1_a,n_cartm2_a,
     &            n_cont_A*n_blk,scr(i_lm),scr(i_lp),scr(i_grdm),
     &            n_all_b,ibug/10)
            end if
*
*         Construct Hessian and add to final array
*
            call dfill (n_ab9blk,zero,scr(i_scr),1)
            call ecp_grad_a (l_A,n_cartp1_a,n_cart_a,n_cartm1_a,
     &          n_cont_A*3*n_blk,scr(i_grdp),scr(i_grdm),scr(i_scr),
     &          n_all_b,ibug/10)
            i_s = i_scr
            do j = 1,3
              do k = 1,3
                i_loc = i_c_A * (i_c_A-1) / 2 + i_c_A
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,k,j,i_loc),1)
                i_loc = i_c_C * (i_c_C-1) / 2 + i_c_C
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,k,j,i_loc),1)
                if (i_c_A .lt. i_c_C) then
                   i_loc = i_c_C * (i_c_C-1) / 2 + i_c_A
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                else  
                   i_loc = i_c_A * (i_c_A-1) / 2 + i_c_C
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                endif
                i_s = i_s+n_abblk
              end do
            end do
          end if
*
        end if
*
*     Mixed gradient terms
*     ====================
*
        if ((i_c_C .ne. i_c_A) .and. (i_c_C .ne. i_c_B)) then
*
*         Calculate integrals for l_A+1, l_B+1
*
          if (.not.DryRun) 
     &        call dfill (n_intppab*n_blk,zero,scr(i_lp),1)
          if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &        R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &        n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &        R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &        n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &        n_prim_C(0,-1,i),n_colc_C(-1,i),
     &        zeta_C(i_zeta),coef_C(i_coef),0,4,
     &        tol,0,scr(i_scr),l_scr,
     &        csco,lcsco,
     &        scr(i_lp),n_intppab,ibug/10)
          if (DryRun) max_type2 = max(mem,max_type2)
          call ecp_int2 (mem,DryRun,
     &        R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &        n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &        R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &        n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &        l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &        ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &        zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &        csco,lcsco,
     &        scr(i_lp),n_intppab,n_blk,ibug/10)
          if (DryRun) max_type2 = max(mem,max_type2)
*
*         Calculate integrals for l_A-1, l_B+1
*
          if (.not.DryRun) 
     &        call dfill (n_intmpab*n_blk,zero,scr(i_lm),1)
          if (l_A .gt. 0) then
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_intmpab,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lm),n_intmpab,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Construct gradient wrt A for l_A, l_B+1.
*
          if (.not.DryRun) then
            call dfill (n_intpb*3*n_blk,zero,scr(i_grdp),1)
            call ecp_grad_a (l_A,n_cartp1_A,n_cart_A,n_cartm1_A,
     &          n_cont_A*n_blk,scr(i_lp),scr(i_lm),scr(i_grdp),
     &          n_allp1_B,ibug/10)
          end if
*
*         Calculate integrals for l_A+1, l_B-1
*
          if (l_B .gt. 0) then
             if (.not.DryRun) then
               call dfill (n_intpmab*n_blk,zero,scr(i_lp),1)
               call dfill (n_intmmab*n_blk,zero,scr(i_lm),1)
             endif
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpmab,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpmab,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
*
*         Calculate integrals for l_A-1, l_B-1
*
            if (l_A .gt. 0) then
              if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_C(i_zeta),coef_C(i_coef),0,4,
     &            tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmmab,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
              call ecp_int2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),l_scr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmmab,n_blk,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
            end if
*
*         Construct gradient wrt A for l_A, l_B-1.
*
            if (.not.DryRun) then
              call dfill (n_intmb*3*n_blk,zero,scr(i_grdm),1)
              call ecp_grad_a (l_A,n_cartp1_A,n_cart_A,n_cartm1_A,
     &            n_cont_A*n_blk,scr(i_lp),scr(i_lm),scr(i_grdm),
     &            n_allm1_B,ibug/10)
            end if
          end if
*
*         Construct gradient wrt B
*
          if (.not.DryRun) then
            call dfill (n_ab9blk,zero,scr(i_scr),1)
            call ecp_grad_b (l_B,n_cartp1_b,n_cart_b,n_cartm1_b,
     &          n_cont_B,scr(i_grdp),scr(i_grdm),scr(i_scr),
     &          n_all_a*3*n_blk,ibug/10)
            i_s = i_scr
            do j = 1,3 ! centre B
              do k = 1,3 ! centre A
                if (i_c_A .eq. i_c_B) then
                   i_loc = i_c_A * (i_c_A-1) / 2 + i_c_B
                   call daxpy (n_abblk,one,scr(i_s),1,
     &                 ecp_hess(1,j,k,i_loc),1)
                endif
                i_loc = i_c_A * (i_c_A-1) / 2 + i_c_B
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,j,k,i_loc),1)
                i_loc = i_c_C * (i_c_C-1) / 2 + i_c_C
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,j,k,i_loc),1)
                call daxpy (n_abblk,one,scr(i_s),1,
     &              ecp_hess(1,k,j,i_loc),1)
                if (i_c_A .lt. i_c_C) then
                   i_loc = i_c_C * (i_c_C-1) / 2 + i_c_A
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                else
                   i_loc = i_c_A * (i_c_A-1) / 2 + i_c_C
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,j,k,i_loc),1)
                endif
                if (i_c_C .lt. i_c_B) then
                   i_loc = i_c_B * (i_c_B-1) / 2 + i_c_C
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,k,j,i_loc),1)
                else
                   i_loc = i_c_C * (i_c_C-1) / 2 + i_c_B
                   call daxpy (n_abblk,-one,scr(i_s),1,
     &                 ecp_hess(1,j,k,i_loc),1)
                endif
                i_s = i_s+n_abblk
              end do
            end do
          end if
*
        end if
      end do
      if (debug_arrays) then
        write (LuOut,*) 'Final ECP hessian'
        write (LuOut,*) '(Lower triangle of Centres)'
        do i = 1,n_cent
          do j = 1,i
            write (LuOut,*) ' Centres',i,j
            do k = 1,3
              printhead(6:9) = dbyd(k)
              do m = 1,3
                printhead(1:4) = dbyd(m)
                i_loc = i * (i-1) / 2 + j
                call ecp_matpr (ecp_hess(1,m,k,i_loc),
     &              1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &              printhead,'E',78,4)
              end do
            end do
          end do
        end do
      end if
      if (DryRun) lscr = i_scr-1
     &    +max(max_type1,max_type2,n_ab9blk)
      if (debug_gen) write (LuOut,*) 'Exiting ecp_hessian'
*
      return
      end
