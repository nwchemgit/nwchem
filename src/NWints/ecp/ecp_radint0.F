C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_radint0 (p_min,p_max,n_prim_c,n_c,n_cont_c,
     &    zeta_c,coeff_c,l_a,n_prim_a,n_cont_a,zeta_a,coeff_a,
     &    l_b,n_prim_b,n_cont_b,zeta_b,coeff_b,prim,temp,wrk,radint,
     &    ibug)
*                                                                      *
*   Calculate Type 0 radial integrals for a given ECP centre, angular  *
*   projector and exponent p                                           *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   n_prim_c (inp) - array of numbers of primitive ECP functions for   *
*                    each power of r                                   *
*   n_c - total number of primitive ECP functions                      *
*   n_cont_c (inp) - number of ECP contractions                        *
*   zeta_c (inp) - exponents on centre c                               *
*   coeff_c (inp) - contraction coefficients on centre c               *
*   l_a (inp) - n quantum number of functions on centre a              *
*   n_prim_a (inp) - number of primitive functions on centre a         *
*   n_cont_a (inp) - number of contracted functions on centre a        *
*   zeta_a (inp) - exponents on centre a                               *
*   coeff_a (inp) - contraction coefficients on centre a               *
*   l_b (inp) - n quantum number of functions on centre b              *
*   n_prim_b (inp) - number of primitive functions on centre b         *
*   n_cont_b (inp) - number of contracted functions on centre b        *
*   zeta_b (inp) - exponents on centre b                               *
*   coeff_b (inp) - contraction coefficients on centre b               *
*   prim - (scr) - primitive integrals                                 *
*   temp - (scr) - integrals contracted over c                         *
*   wrk - (scr) - integrals contracted over b and c                    *
*   radint - (scr) - contracted radial integrals                       *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
#include "ecp_consts.fh"
      integer ibug,i_a,i_b,i_c,ind_c,l_a,n_prim_a,n_cont_a,l_b,
     &    n_prim_b,n_cont_b,n_ab,n_c,n_p,n_prim_c,n_cont_c,
     &    p,pow,p_min,p_max
      logical debug_gen,debug_addresses,debug_arrays
      double precision
     &    prefac,term,fac,
     &    coeff_c(n_c,n_cont_c),zeta_c(n_c),
     &    coeff_a(n_prim_a,n_cont_a),zeta_a(n_prim_a),
     &    coeff_b(n_prim_b,n_cont_b),zeta_b(n_prim_b),
     &    prim(n_prim_a,n_prim_b,n_c),temp(n_prim_a*n_prim_b),
     &    wrk(n_prim_a,n_cont_b),radint(n_cont_a,n_cont_b,n_cont_c)
*
      dimension n_prim_c(p_min:p_max)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
*   Loop over powers of r in ECP nonlocal term
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_radint0 ...'
      fac = zero
      ind_c = 0
      n_ab = n_prim_a*n_prim_b
      call dcopy(n_ab,zero,0,temp,1)
      do p = p_min,p_max
*
*     Calculate prefactor
*
        pow = l_a+l_b+p+1
        i_b = mod(p,2)
        if (i_b .eq. 1) then
          prefac = half
        else
          prefac = qpirt+qpirt
        end if
        do i_c = i_b+1,pow-2,2
          term = i_c
          prefac = prefac*term*half
        end do
C        write (LuOut,*) 'prefac',prefac
C        write (LuOut,*) 'pow',pow
*
*     Calculate primitive integrals
*
        n_p = n_prim_c(p)
        do i_c = 1,n_p
          do i_b = 1,n_prim_b
            do i_a = 1,n_prim_a
              term = one/sqrt(zeta_a(i_a)+zeta_b(i_b)+zeta_c(i_c+ind_c))
C              write (LuOut,*) i_a,i_b,i_c+ind_c,zeta_a(i_a),zeta_b(i_b),
C     &            zeta_c(i_c+ind_c)
              prim(i_a,i_b,i_c+ind_c) = prefac*term**pow
C              write (LuOut,*) i_a,i_b,i_c,term,prim(i_a,i_b,i_c+ind_c)
            end do
          end do
        end do
        ind_c = ind_c+n_p
      end do
*
*     Contract over ecp function
*
      do i_c = 1,n_cont_c
C      write (LuOut,*) 'contract'
C      write (LuOut,*) prefac
C      write (LuOut,*) (temp(i_a),i_a=1,n_ab)
C      write (LuOut,*) (coeff_c(i_a,i_c),i_a=1,n_c)
        call dgemv ('N',n_ab,n_c,one,prim,n_ab,coeff_c(1,i_c),1,zero,
     &      temp,1)
*
*   Contract over basis functions
*
        if (debug_arrays) call ecp_matpr (temp,1,n_prim_a,1,n_prim_b,
     &      1,n_prim_a,1,n_prim_b,'Primitive integrals','E',78,4)
        call dgemm ('N','N',n_prim_a,n_cont_b,n_prim_b,one,temp,
     &      n_prim_a,coeff_b,n_prim_b,zero,wrk,n_prim_a)
        if (debug_arrays) call ecp_matpr (wrk,1,n_prim_a,1,n_cont_b,
     &      1,n_prim_a,1,n_cont_b,'b-contracted integrals','E',78,4)
        call dgemm ('T','N',n_cont_a,n_cont_b,n_prim_a,one,coeff_a,
     &      n_prim_a,wrk,n_prim_a,zero,radint(1,1,i_c),n_cont_a)
        if (debug_arrays) call ecp_matpr (radint(1,1,i_c),1,n_cont_a,1,
     &      n_cont_b,1,n_cont_a,1,n_cont_b,'ECP radial integrals','E',
     &      78,4)
      end do
      if (debug_gen) write (LuOut,*) 'Exiting ecp_radint0'
*
      return
      end
