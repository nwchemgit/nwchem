C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_gradient (
     &    xyz_A,zeta_A,coef_A,n_prim_A,n_cont_A,l_A,i_c_A,
     &    xyz_B,zeta_B,coef_B,n_prim_B,n_cont_B,l_B,i_c_B,
     &    xyz_C,zeta_C,coef_C,n_prim_C,n_colc_C,
     &    ind_z,ind_c,n_zeta_C,n_coef_C,
     &    l_C,i_cent_C,n_C,l_ecp_max,
     &    sphcart,csco,lcsco,
     &    ecp_grads,n_int,n_blk,n_cent,
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
*   i_cent_C (inp) - indices of ECP centres C                          *
*   n_C (inp) - number of ECP centres C                                *
*   l_ecp_max (inp) - maximum angular momentum of any projector on any *
*                 ECP centre                                           *
*   sphcart (inp) - 0 for cartesian, 1 for spherical integrals         *
*   ecp_grads (out) - gradient integrals over ECPs                     *
*   n_int (inp) - number of ECP gradient integrals for a given centre  *
*                 and cartesian derivative. Should be equal to         *
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
#include "util.fh"
#include "errquit.fh"
      integer i,j,k,l,
     &    n_prim_A,n_cont_A,l_A,i_c_A,
     &    n_prim_B,n_cont_B,l_B,i_c_B,
     &    n_zeta_C,n_coef_C,n_C,i_c_C,l_ecp_max,
     &    n_int,n_blk,n_cent,sphcart,lscr,ibug
      integer n_prim_C(0:4,-1:l_ecp_max,n_C*2),
     &    n_colc_C(-1:l_ecp_max,n_C*2),
     &    ind_z(-1:l_ecp_max,n_C*2),ind_c(-1:l_ecp_max,n_C*2),
     &    l_C(n_C),i_cent_C(n_C)
      integer i_scr,i_s,i_zeta,i_coef,i_xca,i_xcb,i_lm,i_lp,
     &    mem,max_type1,max_type2,memscr,
     &    n_cart_a,n_cart_b,n_cart_ab,n_cont_ab,
     &    n_all_a,n_all_b,n_ab,n_ab3,
     &    n_cartm_a,n_cartm_b,n_cartp_a,n_cartp_b,
     &    n_intpa,n_intpb,n_intma,n_intmb
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision 
     &    xyz_A(3),zeta_A(n_prim_A),coef_A(n_prim_A,n_cont_A),
     &    xyz_B(3),zeta_B(n_prim_B),coef_B(n_prim_B,n_cont_B),
     &    xyz_C(3,n_C),zeta_C(n_zeta_C),coef_C(n_zeta_C),
     &    scr(lscr),ecp_grads(n_int*n_blk,3,n_cent)
      double precision 
     &    X_AC,Y_AC,Z_AC,R_AC,X_BC,Y_BC,Z_BC,R_BC,
     &    tol,zero,one
      double precision csco(lcsco)
      character*4 dbyd(3)
      character*12 printhead
      parameter (zero = 0.0d00, one = 1.0d00)
!      data tol/1.0d-14/
      data tol/1.0d-16/
      data dbyd/'d/dx','d/dy','d/dz'/
      data printhead/'      U_ab  '/
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
      if (debug_gen) write (LuOut,*) 'ibug = ',ibug
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_gradient ...'
      if (debug_addresses) then
        write(LuOut,*)' lscr in ecp_gradient:',lscr
        write (LuOut,*) 'n_prim_A,n_cont_A,l_A',n_prim_A,n_cont_A,l_A
        write (LuOut,*) 'n_prim_B,n_cont_B,l_B',n_prim_B,n_cont_B,l_B
        write (LuOut,*) 'l_ecp_max,n_c',l_ecp_max,n_c
      end if
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = n_cont_a*n_cont_b
      n_all_b = n_cart_b*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_ab = n_cart_ab*n_cont_ab
      n_ab3 = n_ab*3
      if (debug_addresses) write (LuOut,*) 
     &    'n_cart_a,n_cart_b,n_cart_ab,n_cont_ab',
     &    n_cart_a,n_cart_b,n_cart_ab,n_cont_ab
      if (n_int .lt. n_ab) call errquit (
     &    'Mismatch of integral count in ecp_gradient',99, BASIS_ERR)
      if (sphcart .ne. 0) call errquit(
     &    'Do your own spherical transformation, lazy bum!',99,
     &       BASIS_ERR)

      max_type1 = 0     ! take care of compiler warnings
      max_type2 = 0

      if (.not.DryRun)
     &  call dfill (n_int*n_blk*3*n_cent,zero,ecp_grads,1)
*
*   Set up pointers to scratch space for coefficients multiplied by 
*   exponents and for integrals with shifted l values
*
      n_cartp_a = n_cart_a+l_A+2
      n_cartp_b = n_cart_b+l_B+2
      n_cartm_a = n_cart_a-l_A-1
      n_cartm_b = n_cart_b-l_B-1
      n_intpa = n_cartp_a*n_cart_b*n_cont_ab
      n_intma = n_cartm_a*n_cart_b*n_cont_ab
      n_intpb = n_cart_a*n_cartp_b*n_cont_ab
      n_intmb = n_cart_a*n_cartm_b*n_cont_ab
      i_xca = 1
      i_xcb = i_xca+n_prim_A*n_cont_A
      i_lp = i_xcb+n_prim_B*n_cont_B
      i_lm = i_lp+max(n_intpa,n_intpb)*n_blk
      i_scr = i_lm+max(n_intma,n_intmb)*n_blk
      if (debug_addresses) then
        write (LuOut,*) 'i_xca,i_xcb',i_xca,i_xcb
        write (LuOut,*) 'i_lp,i_lm',i_lp,i_lm
        write (LuOut,*) 'i_scr',i_scr
      end if
*
*   Set up coefficients multiplied by exponents
*
      if (.not.DryRun) then
        memscr = lscr-i_scr+1
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in ecp_gradient',99, MEM_ERR)
        k = i_xca-1
        do j = 1,n_cont_a
          do i = 1,n_prim_A
            scr(k+i) = zeta_A(i)*coef_A(i,j)
          end do
          k = k+n_prim_A
        end do
        k = i_xcb-1
        do j = 1,n_cont_B
          do i = 1,n_prim_B
            scr(k+i) = zeta_B(i)*coef_B(i,j)
          end do
          k = k+n_prim_B
        end do
      end if
*
*   Loop over ECP centres
*
      if (debug_arrays) call ecp_matpr (xyz_C,1,3,1,n_C,1,3,1,n_C,
     &    'coordinate array','E',78,4)
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
*
*       Derivatives wrt centre B
*       ========================
*
        if (i_c_B .ne. i_c_C) then
*
*         Increment angular momentum on centre B
*
          if (.not.DryRun) then
            call dfill (n_intpb*n_blk,zero,scr(i_lp),1)
            call dfill (n_intmb*n_blk,zero,scr(i_lm),1)
          end if
          if (i_c_A .eq. i_c_C) then
          if ((n_colc_C(-1,i) .gt. 0)) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpb,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            endif
            call ecp_int1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpb,n_blk,.false.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
          if ((n_colc_C(-1,i) .gt. 0)) then
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpb,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            endif
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B+1,
     &          n_prim_B,n_cont_B,scr(i_xcb),zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,
     &          n_prim_A,n_cont_A,coef_A,zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpb,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Decrement angular momentum on centre B
*
          if (l_B .gt. 0) then
            if (i_c_A .eq. i_c_C) then
          if ((n_colc_C(-1,i) .gt. 0)) then
              if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_c(i_zeta),coef_c(i_coef),0,4,
     &            tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmb,.false.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
              endif
              call ecp_int1 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmb,n_blk,.false.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
            else
          if ((n_colc_C(-1,i) .gt. 0)) then
              if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_C(i_zeta),coef_C(i_coef),0,4,
     &            tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmb,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
              endif
              call ecp_int2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B-1,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intmb,n_blk,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
            end if
          end if
*
*         Add gradient terms to output array
*
          if (.not.DryRun) then
            call dfill (n_ab3*n_blk,zero,scr(i_scr),1)
            call ecp_grad_b (l_B,n_cartp_b,n_cart_b,n_cartm_b,n_cont_B,
     &          scr(i_lp),scr(i_lm),scr(i_scr),n_all_a*n_blk,ibug/10)
            i_s = i_scr
            do j = 1,3
              call daxpy (n_ab*n_blk,one,scr(i_s),1,
     &            ecp_grads(1,j,i_c_B),1)
              call daxpy (n_ab*n_blk,-one,scr(i_s),1,
     &            ecp_grads(1,j,i_c_C),1)
              i_s = i_s+n_ab*n_blk
            end do
          end if
        end if
*
*       Derivatives wrt centre A
*       ========================
*
        if (i_c_A .ne. i_c_C) then
*
*         Increment angular momentum on centre A
*
          if (.not.DryRun) then
            call dfill (n_intpa*n_blk,zero,scr(i_lp),1)
            call dfill (n_intma*n_blk,zero,scr(i_lm),1)
          end if
          if (i_c_B .eq. i_c_C) then
          if ((n_colc_C(-1,i) .gt. 0)) then
            if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpa,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
            endif
            call ecp_int1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,
     &          coef_B,zeta_B,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpa,n_blk,.true.,ibug/10)
            if (DryRun) max_type1 = max(mem,max_type1)
          else
          if ((n_colc_C(-1,i) .gt. 0)) then
            if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,0,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpa,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
            endif
            call ecp_int2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,
     &          n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A+1,
     &          n_prim_A,n_cont_A,scr(i_xca),zeta_A,
     &          l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &          ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &          zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr(i_lp),n_intpa,n_blk,ibug/10)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
*
*         Decrement angular momentum on centre A
*
          if (l_A .gt. 0) then
            if (i_c_B .eq. i_c_C) then
          if ((n_colc_C(-1,i) .gt. 0)) then
              if (n_blk .ne. 3) call ecp_local1 (mem,DryRun,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_c(i_zeta),coef_c(i_coef),0,4,
     &            tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intma,.true.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
              endif
              call ecp_int1 (mem,DryRun,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intma,n_blk,.true.,ibug/10)
              if (DryRun) max_type1 = max(mem,max_type1)
            else
          if ((n_colc_C(-1,i) .gt. 0)) then
              if (n_blk .ne. 3) call ecp_local2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            n_prim_C(0,-1,i),n_colc_C(-1,i),
     &            zeta_C(i_zeta),coef_C(i_coef),0,4,
     &            tol,0,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intma,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
           endif
              call ecp_int2 (mem,DryRun,
     &            R_BC,X_BC,Y_BC,Z_BC,l_B,
     &            n_prim_B,n_cont_B,coef_B,zeta_B,
     &            R_AC,X_AC,Y_AC,Z_AC,l_A-1,
     &            n_prim_A,n_cont_A,coef_A,zeta_A,
     &            l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &            ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &            zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &            csco,lcsco,
     &            scr(i_lm),n_intma,n_blk,ibug/10)
              if (DryRun) max_type2 = max(mem,max_type2)
            end if
          end if
*
*         Add gradient terms to output array
*
          if (.not.DryRun) then
            call dfill (n_ab3*n_blk,zero,scr(i_scr),1)
            call ecp_grad_a (l_A,n_cartp_a,n_cart_a,n_cartm_a,
     &          n_cont_A*n_blk,scr(i_lp),scr(i_lm),scr(i_scr),n_all_b,
     &          ibug/10)
            i_s = i_scr
            do j = 1,3
              call daxpy (n_ab*n_blk,one,scr(i_s),1,
     &            ecp_grads(1,j,i_c_A),1)
              call daxpy (n_ab*n_blk,-one,scr(i_s),1,
     &            ecp_grads(1,j,i_c_C),1)
              i_s = i_s+n_ab*n_blk
            end do
          end if
        end if
      end do
      if (debug_arrays) then
        do i = 1,n_cent
          write (LuOut,*) ' Centre',i
          do j = 1,3
            printhead(1:4) = dbyd(j)
            call ecp_matpr (ecp_grads(1,j,i),1,n_all_b,1,n_all_a,
     &          1,n_all_b,1,n_all_a,printhead,'E',78,4)
          end do
        end do
      end if
      if (DryRun) lscr = i_scr-1
     &    +max(n_ab3*n_blk,max_type1,max_type2)
      if (debug_gen) write (LuOut,*) 'Exiting ecp_gradient'
*
      return
      end
