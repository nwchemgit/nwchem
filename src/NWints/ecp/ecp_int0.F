C $Id: ecp_int0.F,v 1.4 1997-08-20 18:30:07 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine ecp_int0 (mem_max,DryRun,
     &    l_a,n_prim_a,n_cont_a,coef_a,zeta_a,n_na,
     &    l_b,n_prim_b,n_cont_b,coef_b,zeta_b,n_nb,
     &    l_ecp_max,n_prim_c,n_coef_c,ind_c,n_c_tot,coef_c,zeta_c,
     &    p_min,p_max,tol,sphcart,tmp,l_tmp,
     &    csco,lcsco,
     &    ecp_ints,n_blk,ibug)
*                                                                      *
*   Calculate Type 0 nonlocal integrals for a given ECP centre         *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   mem_max (out) - maximum scratch memory required                    *
*   DryRun (inp) - logical to only return memory if true               *
*   R_AC (inp) - distance between centres A and C                      *
*   X_AC,Y_AC,Z_AC (inp) - cartesian coordinates of centre C relative  *
*                          to centre A, X_AC = X_C - X_A, etc.         *
*   l_a (inp) - (maximum) angular momentum of functions on centre A    *
*   n_prim_a (inp) - number of primitive functions on centre A         *
*   n_cont_a (inp) - number of contracted functions on centre A        *
*   coef_a - contraction coefficients of basis functions on centre A   *
*   zeta_a - exponents of primitive basis functions on centre A        *
*   R_BC (inp) - distance between centres B and C                      *
*   X_BC,Y_BC,Z_BC (inp) - cartesian coordinates of centre C relative  *
*                          to centre B, X_BC = X_C - X_B, etc.         *
*   l_b (inp) - (maximum) angular momentum of functions on centre B    *
*   n_prim_b (inp) - number of primitive functions on centre B         *
*   n_cont_b (inp) - number of contracted functions on centre B        *
*   coef_b - contraction coefficients of basis functions on centre B   *
*   zeta_b - exponents of primitive basis functions on centre B        *
*   l_ecp_max (inp) - maximum angular momentum of ECP projector        *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_coef_c (inp) - array of number of coefficients/exponents for     *
*                    each projector                                    *
*   ind_c (inp) - address of first coefficient/exponent for each l     *
*   n_c_tot (inp) - total number of ECP coefficients/exponents         *
*   coef_c - ECP contraction coefficients                              *
*   zeta_c - ECP exponents                                             *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart (inp) - 1 for spherical basis, 0 for cartesian basis.      *
*   n_blk (inp) - 1 for scalar only, 3 for s-o only, 4 for both        *
*   tmp - work array                                                   *
*   l_tmp - length of tmp                                              *
*   ecp_ints - integrals over ECP                                      *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Notes:                                                             *
*   -----                                                              *
*                                                                      *
*   The integrals come out in the order cmpt_a, cont_a, cmpt_b, cont_b *
*      where cont = contracted functions, cmpt = cartesian components  *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
#include "ecp_consts.fh"
      integer l_a,n_prim_a,n_cont_a,n_na,l_b,n_prim_b,n_cont_b,n_nb,
     &    l_ecp_max,n_c_tot,l_tmp,p_min,p_max,sphcart,mem_max,ibug
      integer n_coef_c(-1:l_ecp_max),n_prim_c(p_min:p_max,-1:l_ecp_max),
     &    ind_c(-1:l_ecp_max)
      integer i,j,m,i_wrk,i_prm,i_tmp,i_Q_int,i_scr,i_free,i_ang,l_c,
     &    l_c_min,l_c_max,l_cmp,n_ab,nc_ab,ncab,n_l,n_c,n_blk,n_cont_c
      integer n_all_a,n_all_b,n_ang,n_x
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision coef_c(n_c_tot),zeta_c(n_c_tot),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(l_tmp),ecp_ints(n_na,n_cont_a,n_nb,n_cont_b,n_blk),
     &    tol,wa,wb,wc
      double precision csco(lcsco)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_int0 ...'
*
      if (debug_addresses) then
        write (LuOut,*) 'n_prim_a,n_prim_b,n_cont_a,n_cont_b',
     &      n_prim_a,n_prim_b,n_cont_a,n_cont_b
        write (LuOut,*) 'n_c_tot',n_c_tot
      end if
      n_ab = n_prim_a*n_prim_b
      nc_ab = n_prim_a*n_cont_b
      ncab = n_cont_a*n_cont_b
      n_all_a = n_cont_a*n_na
      n_all_b = n_cont_b*n_nb
      if (debug_addresses) then
        write (LuOut,*) 'n_ab,nc_ab,ncab',n_ab,nc_ab,ncab
        write (LuOut,*) 'mem_max',mem_max
      end if
      if (DryRun) mem_max = 0
*
*   Check angular momenta
*
      if (sphcart .eq. 0) then
        if (mod(l_a+l_b,2) .ne. 0) return
        l_c_min = mod(l_b,2)
        l_c_max = min(l_a,l_b,l_ecp_max)
      else
        if (l_a .ne. l_b) return
        if (l_a .gt. l_ecp_max) return
        l_c_min = l_a
        l_c_max = l_a
      end if
*
*     Loop over ECP angular momentum projectors
*
C      write (LuOut,*) 'Starting l_c loop'
      do l_c = l_c_min,l_c_max,2
        if (l_c .gt. 0) then
          n_cont_c = n_blk/4+1
          n_x = n_blk
        else if (n_blk .ne. 3) then
          n_cont_c = 1
          n_x = 1
        else
          go to 10
        end if
        n_c = n_coef_c(l_c)
        i_Q_int = 1
        i_tmp = i_Q_int+ncab*n_cont_c
        i_prm = i_tmp+n_ab*n_c
        i_wrk = i_prm+n_ab
        i_free = i_wrk+nc_ab
        if (debug_addresses) write (LuOut,*) 'i_Q_int,i_tmp,i_prm,',
     &      'i_wrk,i_free,mem_max',
     &      i_Q_int,i_tmp,i_prm,i_wrk,i_free,mem_max
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
        else
          if (i_free-1 .gt. l_tmp) call errquit(
     &        ' Insufficient memory in ecp_int0',99)
          i = ind_c(l_c)
          if (debug_addresses) write (LuOut,*) 'l_c,ind_c(l_c),n_c',
     &        l_c,ind_c(l_c),n_c
          call ecp_radint0 (p_min,p_max,
     &        n_prim_c(p_min,l_c),n_c,n_cont_c,zeta_c(i),coef_c(i),
     &        l_a,n_prim_a,n_cont_a,zeta_a,coef_a,
     &        l_b,n_prim_b,n_cont_b,zeta_b,coef_b,
     &        tmp(i_prm),tmp(i_tmp),tmp(i_wrk),tmp(i_Q_int),
     &        (ibug/10))
        end if
*
        n_l = 2*l_c+1
        if (sphcart .eq. 0) then
*
*       Cartesian integrals.
*
          wa = n_l
          wb = (pi+pi)/wa
          if (debug_addresses) write (LuOut,*) 'n_l,n_na',n_l,n_na
          i_ang = i_tmp+n_l*n_na
          n_ang = n_l**2
          i_scr = i_ang+n_ang*n_x
          i_free = i_scr+max(n_na*n_nb,n_ang)
          if (debug_addresses) write (LuOut,*) 'i_ang,i_free',i_ang,
     &        i_free
          if (DryRun) then
            mem_max = max(mem_max,i_free-1)
          else
            if (i_free-1 .gt. l_tmp) call errquit(
     &          ' Insufficient memory in ecp_int0',99)
*
*       Set up angular overlap/spin-orbit matrix
*
            call dfill (n_ang*n_x,zero,tmp(i_ang),1)
*
*         scalar potential
*
            i = i_ang
            if (n_x .ne. 3) then
              call dfill (n_l,wb,tmp(i),n_l+1)
              tmp(i_ang+2*l_c*(l_c+1)) = wb+wb
              i = i+n_ang
            end if
*
*         z spin-orbit potential 
*
            if (n_x .ne. 1) then
              j = i
              do m = -l_c,l_c
                j = j+l_c+l_c
                tmp(j) = m
                tmp(j) = tmp(j)*wb
              end do
*
*         y,x spin-orbit potentials
*
              i = i+n_ang
              j = i+n_ang-1
              do m = -l_c,l_c-1
                wc = (l_c-m)*(l_c+m+1)
                if ((m .gt. 0) .or. (m .lt. -1)) wc = wc*half
                wc = sqrt(wc)*wb
                i = i+1
                tmp(i) = wc
                i = i+l_c+l_c
                tmp(i) = -wc
                i = i+1
                j = j+l_c+l_c
                tmp(j) = -wc
              end do
            end if
*
            i = i_ang
            do j = 1,n_x
              if (debug_arrays) call ecp_matpr (tmp(i),-l_c,l_c,
     &            -l_c,l_c,-l_c,l_c,-l_c,l_c,'Spherical tensor overlap',
     &            'F',78,4)
              call ecp_cstrans (l_a,n_na,n_l,l_c,l_c,l_cmp,tmp(i_tmp),
     &            n_na,tmp(i),n_l,csco,lcsco,csco,-1,-1,1)
              if (debug_arrays) call ecp_matpr (tmp(i_tmp),1,n_na,-l_c,
     &            l_c,1,n_na,-l_c,l_c,'Transformed for l_a','F',78,4)
              call ecp_cstrans (l_b,n_nb,n_na,l_c,l_c,l_cmp,tmp(i_scr),
     &            n_na,tmp(i_tmp),n_na,csco,lcsco,csco,-1,1,1)
              if (debug_arrays) call ecp_matpr (tmp(i_scr),1,n_na,
     &            1,n_nb,1,n_na,1,n_nb,'Angular integrals','F',78,4)
              call ecp_angrad (n_na,n_cont_a,n_nb,n_cont_b,
     &            tmp(i_scr),tmp(i_Q_int),ecp_ints(1,1,1,1,j))
              if (debug_arrays) call ecp_matpr (ecp_ints(1,1,1,1,j),
     &            1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &            'Assembled integrals','E',78,4)
              if ((j .eq. 1) .and. (n_x .eq. 4)) 
     &            i_Q_int = i_Q_int+ncab
            end do
          end if
        else
          if (n_blk .gt. 1) call errquit (
     &        'Cannot do spherical ecp spin-orbit type 0 integrals',
     &        99)
*
*       Spherical integrals (assume here that sphericals are normalized)
*
          wa = 1
          m = i_Q_int
          do i = 1,n_cont_b
            do j = 1,n_l
              call daxpy (n_cont_a,wa,tmp(m),1,ecp_ints(j,1,j,i,1),n_l)
            end do
            m = m+n_cont_a
          end do
        end if
   10   continue
      end do
      if (debug_arrays) call ecp_matpr (ecp_ints,1,n_all_b,1,n_all_a,
     &    1,n_all_b,1,n_all_a,'ECP integrals','E',78,4)
      if (debug_gen) write (LuOut,*) 'Exiting ecp_int0'
*
      return
      end
