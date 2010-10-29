C $Id$
************************************************************************
*                                                                      *
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{ecp\_integral}
c:tex-this routine computes the scalar and spin-orbit ECP integrals. The 
c:tex-ECPs are defined by
c:tex-\begin{eqnarray*}
c:tex-U^{AREP} = U^{AREP}_L + \sum_{l=0}^{L-1} U^{AREP}_\ell 
c:tex-\sum_{m=-\ell}^\ell | \ell m \rangle \langle \ell m | \\
c:tex-U^{so} = \sum_{\ell=1}^L \frac1{\ell} \frac{2}{2\ell+1} 
C:tex-\Delta U^{so}_\ell 
c:tex-\sum_{m=-\ell}^\ell \sum_{m'=-\ell}^\ell | \ell m \rangle 
c:tex-\langle \ell m | -i\hat\ell | \ell m' \rangle \langle \ell m' |
c:tex-\end{eqnarray*}
c:tex-Note that the spin-orbit operator is over $-i\ell$ and that the
c:tex-scalar product with $\sigma$ not {\bf s} should be taken.
c:tex-Note also that the integrals come out in the order z,y,x. 
c:tex-
c:tex-
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine ecp_integral (
     &    xyz_A,zeta_A,coef_A,n_prim_A,n_cont_A,l_A,i_c_A,
     &    xyz_B,zeta_B,coef_B,n_prim_B,n_cont_B,l_B,i_c_B,
     &    xyz_C,zeta_C,coef_C,n_prim_C,n_colc_C,
     &    ind_z,ind_c,n_zeta_C,n_coef_C,
     &    l_C,i_cent_C,n_C,l_ecp_max,
     &    sphcart,csco,lcsco,
     &    ecp_ints,n_int,n_blk,
     &    DryRun,scr,lscr,ibug)
c:tex-\end{verbatim}
*                                                                      *
*   Calculate scalar and spin-orbit ecp integrals.                     *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   xyz_A (inp) - coordinates of centre A                              *
*   zeta_A (inp) - exponents of primitive gaussians on centre A        *
*   coef_A (inp) - contraction coefficients on centre A                *
*   n_prim_A (inp) - number of primitive gaussians on centre A         *
*   n_cont_A (inp) - number of contracted functions on centre A        *
*   l_A (inp) - angular momentum of functions on centre A              *
*   i_c_A (inp) index of centre A in centre list                       *
*   xyz_B (inp) - coordinates of centre B                              *
*   zeta_B (inp) - exponents of primitive gaussians on centre B        *
*   coef_B (inp) - contraction coefficients on centre B                *
*   n_prim_B (inp) - number of primitive gaussians on centre B         *
*   n_cont_B (inp) - number of contracted functions on centre B        *
*   l_B (inp) - angular momentum of functions on centre B              *
*   i_c_B (inp) index of centre B in centre list                       *
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
*                     ECP centre.                                      *
*   sphcart - 0 for cartesian integrals, 1 for spherical integrals     *
*   ecp_ints (out) - integrals over ECPs                               *
*   n_int (inp) - number of ECP integrals. Should be equal to          *
*                 NCA*NCB*[(La+1)*(La+2)/2]*[(Lb+1)*(Lb+2)/2]          *
*   n_blk (inp) - 1 for scalar only, 3 for s-o only, 4 for both        *
*                 s-o integrals come out in the order z,y,x.           *
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
      integer i,j,l,n,
     &    n_prim_A,n_cont_A,l_A,i_c_A,
     &    n_prim_B,n_cont_B,l_B,i_c_B,
     &    n_zeta_C,n_coef_C,n_C,i_c_C,l_ecp_max,
     &    n_int,n_blk,sphcart,lscr,ibug
      integer n_prim_C(0:4,-1:l_ecp_max,n_C*2),
     &    n_colc_C(-1:l_ecp_max,n_C*2),
     &    ind_z(-1:l_ecp_max,n_C*2),ind_c(-1:l_ecp_max,n_C*2),
     &    l_C(n_C),i_cent_C(n_C)
      integer i_scr,i_coef,i_zeta,memscr,mem,
     &    max_type0,max_type1,max_type2,n_intt,n_all_a,n_all_b,
     &    n_cart_a,n_cart_b,n_cart_ab,n_cont_ab
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision 
     &    xyz_A(3),zeta_A(n_prim_A),coef_A(n_prim_A,n_cont_A),
     &    xyz_B(3),zeta_B(n_prim_B),coef_B(n_prim_B,n_cont_B),
     &    xyz_C(3,n_C),zeta_C(n_zeta_C),coef_C(n_zeta_C),
     &    scr(lscr),ecp_ints(n_int*n_blk)
      double precision 
     &    X_AC,Y_AC,Z_AC,R_AC,X_BC,Y_BC,Z_BC,R_BC,
     &    tol
      double precision csco(lcsco)
      data tol/1.0d-16/
*
      if ((n_blk .ne. 1) .and. (n_blk .ne. 3) .and. (n_blk .ne. 4))
     &    call errquit('Illegal value of n_blk in ecp_integral',99,
     &       INT_ERR)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      if (debug_gen) write (LuOut,*) ibug
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_integral ...'
      if (debug_addresses) then
        write(LuOut,*)' lscr in ecp_integral:',lscr
        write (LuOut,*) 'n_prim_A,n_cont_A,l_A',n_prim_A,n_cont_A,l_A
        write (LuOut,*) 'n_prim_B,n_cont_B,l_B',n_prim_B,n_cont_B,l_B
        write (LuOut,*) 'l_ecp_max,n_c',l_ecp_max,n_c
      end if
      n_cart_a = (l_a+1)*(l_a+2)/2
      n_cart_b = (l_b+1)*(l_b+2)/2
      n_cart_ab = n_cart_a*n_cart_b
      n_cont_ab = n_cont_a*n_cont_b
      n_all_a = n_cart_a*n_cont_a
      n_all_b = n_cart_b*n_cont_b
      n_intt = n_cart_ab*n_cont_ab
      if (debug_addresses) write (LuOut,*) 
     &    'n_cart_a,n_cart_b,n_cart_ab,n_cont_ab',
     &    n_cart_a,n_cart_b,n_cart_ab,n_cont_ab
      if (sphcart .ne. 0) call errquit(
     &    'Do your own spherical transformation, lazy bum!',99,
     &       BASIS_ERR)
      if (n_int .lt. n_intt) then
        write (LuOut,*) 'n_int  =',n_int
        write (LuOut,*) 'n_intt =',n_intt
        call errquit (
     &      'Mismatch of integral count in ecp_integrals',99,
     &       BASIS_ERR)
      end if
      i_scr = n_intt*n_blk+1
      memscr = lscr-i_scr+1
      if (DryRun) then
        max_type0 = 0
        max_type1 = 0
        max_type2 = 0
      else
        call dcopy (n_intt*n_blk,zero,0,ecp_ints,1)
        if (debug_addresses) write (LuOut,*) 'memscr =',memscr
        if (memscr .lt. 0) call errquit (
     &      'Insufficient scratch memory in ecp_integral',99, MEM_ERR)
      end if
      if (debug_addresses) write (LuOut,*) 'i_scr',i_scr
      if (debug_arrays) then 
        call ecp_matpr(xyz_A,1,3,1,1,1,3,1,1,'Centre A coordinates',
     &      'E',78,4)
        call ecp_matpr(xyz_B,1,3,1,1,1,3,1,1,'Centre B coordinates',
     &      'E',78,4)
        call ecp_matpr(xyz_C,1,3,1,n_C,1,3,1,n_C,'ECP coordinate array',
     &      'E',78,4)
        call ecp_matpr(zeta_A,1,n_prim_A,1,1,1,n_prim_A,1,1,
     &      'A exponents','E',78,4)
        call ecp_matpr(coef_A,1,n_prim_A,1,n_cont_A,1,n_prim_A,
     &      1,n_cont_A,'A coefficients','E',78,4)
        call ecp_matpr(zeta_B,1,n_prim_B,1,1,1,n_prim_B,1,1,
     &      'B exponents','E',78,4)
        call ecp_matpr(coef_B,1,n_prim_B,1,n_cont_B,1,n_prim_B,
     &      1,n_cont_B,'B coefficients','E',78,4)
        call ecp_matpr(zeta_C,1,n_zeta_C,1,1,1,n_zeta_C,1,1,
     &      'ECP exponents','E',78,4)
        call ecp_matpr(coef_C,1,n_zeta_C,1,1,1,n_zeta_C,1,1,
     &      'ECP coefficients','E',78,4)
      end if
      if (debug_gen) write (LuOut,*) 'Number of ECP centers =',n_C
*
*   Loop over ECP centres
*
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
          write (LuOut,'(3x,A,3F10.6)') 'Distance to centers A and B:',
     &        R_AC,R_BC
        end if
*
*     Evaluate integrals
*
        if (debug_gen) write (LuOut,*) 'starting integrals'
C        if (debug_arrays) then
C          write (LuOut,*) 'Parameters of ECP'
C          call ecp_matpi (n_prim_C,0,4,-1,l_ecp_max,0,2,0,l_C,
C     &        'n_prim_C',81,15)
C          call ecp_matpi (ind_C,-1,l_ecp_max,1,2*n_C,0,l_C,1,2*n_C,
C     &        'ind_C',81,15)
C          call ecp_matpi (n_colc_C,-1,l_ecp_max,1,2*n_C,0,l_C,1,2*n_C,
C     &        'n_colc_C',81,15)
C          call ecp_matpr (zeta_C,1,n_zeta_C,1,1,i_zeta,i_zeta+n-1,1,1,
C     &        'zeta_C','E',81,5)
C          call ecp_matpr (coef_C,1,n_zeta_C,1,1,i_zeta,i_zeta+n-1,1,1,
C     &        'coef_C','E',81,5)
C        end if
        if (.not.DryRun) call dcopy (n_intt*n_blk,zero,0,scr,1)
        i_zeta = ind_z(-1,i)
        i_coef = ind_c(-1,i)
        if ((i_c_A .eq. i_c_C) .and. (i_c_B .eq. i_c_C)) then
*
*       One-centre case, A = B = C
*
          if (debug_gen) write (LuOut,*) 'One-centre case'
          if ((n_colc_C(-1,i) .gt. 0)) then
            if ((n_blk .ne. 3)) 
     &          call ecp_local0 (mem,DryRun,
     &          l_B,n_prim_B,n_cont_B,coef_B,zeta_B,n_cart_b,
     &          l_A,n_prim_A,n_cont_A,coef_A,zeta_A,n_cart_a,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr,ibug/10)
            if (debug_arrays .and. .not.DryRun) 
     &          call ecp_matpr (scr,1,n_all_b,1,n_all_a,1,n_all_b,
     &          1,n_all_a,'Local ECP integrals','E',78,4)
            if (DryRun) max_type0 = max(mem,max_type0)
          end if
          call ecp_int0 (mem,DryRun,
     &        l_B,n_prim_B,n_cont_B,coef_B,zeta_B,n_cart_b,
     &        l_A,n_prim_A,n_cont_A,coef_A,zeta_A,n_cart_a,
     &        l_C(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &        ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &        zeta_C,coef_C,0,4,tol,sphcart,
     &        scr(i_scr),memscr,csco,lcsco,
     &        scr,n_blk,ibug/10)
          if (DryRun) max_type0 = max(mem,max_type0)
        else if (i_c_A .eq. i_c_C) then
*
*       Two-centre case with A = C
*
          if (debug_gen) write (LuOut,*) 'Two-centre case A=C'
          if ((n_colc_C(-1,i) .gt. 0)) then
            if ((n_blk .ne. 3))
     &          call ecp_local1 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr,n_intt,.false.,ibug/10)
            if (debug_arrays .and. .not.DryRun) 
     &          call ecp_matpr (scr,1,n_all_b,1,n_all_a,1,n_all_b,
     &          1,n_all_a,'Local ECP integrals','E',78,4)
            if (DryRun) max_type1 = max(mem,max_type1)
          end if
          call ecp_int1 (mem,DryRun,
     &        R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &        R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &        l_C(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &        ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &        zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &        csco,lcsco,
     &        scr,n_intt,n_blk,.false.,ibug/10)
          if (DryRun) max_type1 = max(mem,max_type1)
        else if (i_c_B .eq. i_c_C) then
*
*       Two-centre case with B = C
*
          if (debug_gen) write (LuOut,*) 'Two-centre case B=C'
          if ((n_colc_C(-1,i) .gt. 0)) then
            if ((n_blk .ne. 3))
     &          call ecp_local1 (mem,DryRun,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_c(i_zeta),coef_c(i_coef),0,4,
     &          tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr,n_intt,.true.,ibug/10)
            if (debug_arrays .and..not.DryRun) 
     &          call ecp_matpr (scr,1,n_all_b,1,n_all_a,1,n_all_b,
     &          1,n_all_a,'Local ECP integrals','E',78,4)
            if (DryRun) max_type1 = max(mem,max_type1)
          end if
          call ecp_int1 (mem,DryRun,
     &        R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &        R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &        l_C(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &        ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &        zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &        csco,lcsco,
     &        scr,n_intt,n_blk,.true.,ibug/10)
          if (DryRun) max_type1 = max(mem,max_type1)
        else
*
*     Three-centre case
*
          if (debug_gen) write (LuOut,*) 'Three-centre case'
          if ((n_colc_C(-1,i) .gt. 0))then
            if ((n_blk .ne. 3))
     &          call ecp_local2 (mem,DryRun,
     &          R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &          R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &          n_prim_C(0,-1,i),n_colc_C(-1,i),
     &          zeta_C(i_zeta),coef_C(i_coef),0,4,
     &          tol,sphcart,scr(i_scr),memscr,
     &          csco,lcsco,
     &          scr,n_intt,ibug/10)
            if (debug_arrays .and..not.DryRun) 
     &          call ecp_matpr (scr,1,n_all_b,1,n_all_a,1,n_all_b,
     &          1,n_all_a,'Local ECP integrals','E',78,4)
            if (DryRun) max_type2 = max(mem,max_type2)
          end if
          call ecp_int2 (mem,DryRun,
     &        R_BC,X_BC,Y_BC,Z_BC,l_B,n_prim_B,n_cont_B,coef_B,zeta_B,
     &        R_AC,X_AC,Y_AC,Z_AC,l_A,n_prim_A,n_cont_A,coef_A,zeta_A,
     &        l_c(i),n_prim_C(0,-1,i),n_colc_C(-1,i),ind_z(-1,i),
     &        ind_c(-1,i),n_zeta_C,n_coef_C,l_ecp_max,n_C,
     &        zeta_C,coef_C,0,4,tol,sphcart,scr(i_scr),memscr,
     &        csco,lcsco,
     &        scr,n_intt,n_blk,ibug/10)
          if (DryRun) max_type2 = max(mem,max_type2)
        end if
*
*     Add integrals from this centre to the total.
*
        if (.not.DryRun) then
          call daxpy (n_intt*n_blk,one,scr,1,ecp_ints,1)
          if (debug_arrays) then
            n = 1
            do j = 1,n_blk
              call ecp_matpr (scr(n),1,n_all_b,1,n_all_a,1,n_all_b,
     &            1,n_all_a,'Non-local ECP integrals','E',78,4)
              n = n+n_intt
            end do
          end if
        end if
      end do
      if (DryRun) then
        lscr = max(max_type0,max_type1,max_type2)+n_intt*n_blk
        if (debug_addresses) then
          write (LuOut,*) 'max_type0',max_type0
          write (LuOut,*) 'max_type1',max_type1
          write (LuOut,*) 'max_type2',max_type2
          write (LuOut,*) 'n_intt',n_intt
          write (LuOut,*) 'lscr',lscr
        end if
      else if (debug_arrays) then
        n = 1
        do j = 1,n_blk
          call ecp_matpr (ecp_ints(n),1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'Final ECP integrals','E',78,4)
          n = n+n_intt
        end do
      end if
      if (debug_gen) write (LuOut,*) 'Exiting ecp_integral'
*
      return
      end
