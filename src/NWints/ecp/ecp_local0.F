C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_local0 (mem_max,DryRun,
     &    l_a,n_prim_a,n_cont_a,coef_a,zeta_a,n_na,
     &    l_b,n_prim_b,n_cont_b,coef_b,zeta_b,n_nb,
     &    n_prim_c,n_coef_c,zeta_c,coef_c,p_min,p_max,
     &    tol,sphcart,tmp,l_tmp,
     &    csco,lcsco,
     &    ecp_ints,ibug)
*                                                                      *
*   Calculate Type 0 local integrals for a given ECP centre            *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   mem_max (out) - maximum scratch memory required                    *
*   DryRun (inp) - logical to only return memory if true               *
*   l_a (inp) - (maximum) angular momentum of functions on centre A    *
*   n_prim_a (inp) - number of primitive functions on centre A         *
*   n_cont_a (inp) - number of contracted functions on centre A        *
*   coef_a (inp) - centre A contraction coefficients                   *
*   zeta_a (inp) - centre A exponents                                  *
*   n_na (inp) - number of cartesian components for centre A           *
*   l_b (inp) - (maximum) angular momentum of functions on centre B    *
*   n_prim_b (inp) - number of primitive functions on centre B         *
*   n_cont_b (inp) - number of contracted functions on centre B        *
*   coef_b (inp) - centre B contraction coefficients                   *
*   zeta_b (inp) - centre B exponents                                  *
*   n_nb (inp) - number of cartesian components for centre B           *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_coef_c (inp) - number of coefficients/exponents in local potl.   *
*   zeta_c - ECP exponents                                             *
*   coef_c - ECP contraction coefficients                              *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart (inp) - 1 for spherical basis, 0 for cartesian basis.      *
*   tmp (scr) - work array                                             *
*   l_tmp (inp) - length of tmp                                        *
*   csco (inp) - cartesian/spherical transformation                    *
*   lcsco (inp) - length of cartesian/spherical transformation         *
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
#include "util.fh"
#include "errquit.fh"
      integer l_a,n_prim_a,n_cont_a,n_na,l_b,n_prim_b,n_cont_b,n_nb,
     &    n_coef_c,l_tmp,p_min,p_max,sphcart,mem_max,ibug
      integer n_prim_c(p_min:p_max)
      integer i,j,m,i_wrk,i_prm,i_tmp,i_Q_int,i_scr,i_free,i_ang,
     &    l_c,l_c_min,l_c_max,l_cmp,
     &    n_ab,nc_ab,ncab,n_l
      integer n_all_a,n_all_b,n_ang
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision zeta_c(n_coef_c),coef_c(n_coef_c),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(l_tmp),ecp_ints(n_na,n_cont_a,n_nb,n_cont_b),
     &    tol,wa,wb
      double precision csco(lcsco)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_local0 ...'
*
      if (debug_addresses) then
        write (LuOut,*) 'Scratch memory l_tmp = ',l_tmp
        write (LuOut,*) 'n_prim_a,n_prim_b,n_cont_a,n_cont_b',
     &      n_prim_a,n_prim_b,n_cont_a,n_cont_b
      end if
      n_ab = n_prim_a*n_prim_b
      nc_ab = n_prim_a*n_cont_b
      ncab = n_cont_a*n_cont_b
      n_all_a = n_cont_a*n_na
      n_all_b = n_cont_b*n_nb
      if (debug_addresses) then
        write (LuOut,*) 'n_ab,nc_ab,ncab',n_ab,nc_ab,ncab
      end if
      if (DryRun) mem_max = 0
*
*   Check angular momenta
*
      if (sphcart .eq. 0) then
        if (mod(l_a+l_b,2) .ne. 0) return
        l_c_min = mod(l_b,2)
        l_c_max = min(l_a,l_b)
      else
        if (l_a .ne. l_b) return
        l_c_min = l_a
        l_c_max = l_a
      end if
*
*     Set up pointers for radial integral evaluation
*
      i_Q_int = 1
      i_tmp = i_Q_int+ncab
      i_prm = i_tmp+n_ab
      i_wrk = i_prm+n_ab*n_coef_c
      i_free = i_wrk+nc_ab
      if (debug_addresses) write (LuOut,*) 
     &    'i_Q_int,i_tmp,i_prm,i_wrk,i_free',
     &    i_Q_int,i_tmp,i_prm,i_wrk,i_free
      if (DryRun) then
        mem_max = max(mem_max,i_free-1)
        if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
      else
        if (i_free-1 .gt. l_tmp) call errquit(
     &      ' Insufficient memory in ecp_local0',99, MEM_ERR)
        call ecp_radint0 (p_min,p_max,
     &      n_prim_c,n_coef_c,1,zeta_c,coef_c,
     &      l_a,n_prim_a,n_cont_a,zeta_a,coef_a,
     &      l_b,n_prim_b,n_cont_b,zeta_b,coef_b,
     &      tmp(i_prm),tmp(i_tmp),tmp(i_wrk),tmp(i_Q_int),
     &      (ibug/10))
      end if
*
      if (sphcart .eq. 0) then
*
*   Cartesian integrals.
*
*
*    Loop over l values
*
        if (debug_addresses) write (LuOut,*) 'l_c_min,l_c_max',
     &      l_c_min,l_c_max
        do l_c = l_c_min,l_c_max,2
          n_l = 2*l_c+1
          wa = n_l
          wa = pi/wa
          wb = wa+wa
          if (debug_addresses) write (LuOut,*) 'n_l,n_na',
     &        n_l,n_na
          i_ang = i_tmp+n_l*n_na
          n_ang = n_l**2
          i_scr = i_ang+n_ang
          i_free = i_scr+max(n_na*n_nb,n_ang)
          if (debug_addresses) write (LuOut,*) 'i_ang,n_ang,i_free',
     &        i_ang,n_ang,i_free
          if (DryRun) then
            mem_max = max(mem_max,i_free-1)
            if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
          else
            if (i_free-1 .gt. l_tmp) call errquit(
     &          ' Insufficient memory in ecp_local0',99, MEM_ERR)
*
*       Set up angular overlap matrix
*
            call dcopy (n_ang,zero,0,tmp(i_ang),1)
            call dfill (n_l,wb,tmp(i_ang),n_l+1)
            tmp(i_ang+2*l_c*(l_c+1)) = wb+wb
*
            if (debug_arrays) call ecp_matpr (tmp(i_ang),-l_c,l_c,
     &          -l_c,l_c,-l_c,l_c,-l_c,l_c,'Spherical tensor overlap',
     &          'F',78,4)
            call ecp_cstrans (l_a,n_na,n_l,l_c,l_c,l_cmp,tmp(i_tmp),
     &          n_na,tmp(i_ang),n_l,csco,lcsco,csco,-1,-1,1)
            if (debug_arrays) call ecp_matpr (tmp(i_tmp),1,n_na,-l_c,
     &          l_c,1,n_na,-l_c,l_c,'Transformed for l_a','F',78,4)
            call ecp_cstrans (l_b,n_nb,n_na,l_c,l_c,l_cmp,tmp(i_scr),
     &          n_na,tmp(i_tmp),n_na,csco,lcsco,csco,-1,1,1)
            if (debug_arrays) call ecp_matpr (tmp(i_scr),1,n_na,
     &          1,n_nb,1,n_na,1,n_nb,'Angular integrals','F',78,4)
            call ecp_angrad (n_na,n_cont_a,n_nb,n_cont_b,
     &          tmp(i_scr),tmp(i_Q_int),ecp_ints)
            if (debug_arrays) call ecp_matpr (ecp_ints,
     &          1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &          'Assembled integrals','E',78,4)
          end if
        end do
      else
*
*       Spherical integrals (assume here that sphericals are normalized)
*
        wa = 1
        m = i_Q_int
        n_l = 2*l_a+1
        do i = 1,n_cont_b
          do j = 1,n_l
            call daxpy (n_cont_a,wa,tmp(m),1,ecp_ints(j,1,j,i),n_l)
          end do
          m = m+n_cont_a
        end do
      end if
      if (debug_arrays) then
        call ecp_matpr (ecp_ints,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'ECP integrals','E',78,4)
      end if
      if (debug_gen) write (LuOut,*) 'Exiting ecp_local0'
*
      return
      end
