C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_locint1 (p_min,p_max,l_max,n_prim_a,n_prim_b,
     &    n_c,n_prim_c,coeff_c,bi,a_sqrt,x,prefactor,tol,
     &    n_ab,n_abc,n_int,temp,ind,Qabc,Qa,Q_int,ibug)
*                                                                      *
*   Calculate Type 1 radial integrals Q^k_m for a given ECP centre to  *
*   be used in the local part of the ECP                               *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   l_max (inp) -  maximum k value = l_b+l_a                           *
*   n_prim_a (inp) - number of primitive functions on centre a         *
*   n_prim_b (inp) - number of primitive functions on centre b         *
*   n_c (inp) - total number of ECP primitive functions                *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   coeff_c - contraction coefficients of potential on centre c        *
*   bi (inp) - values of 1/beta                                        *
*   a_sqrt (inp) - values of 1/sqrt(alpha)                             *
*   x (inp) - values of alpha/2*sqrt(beta)                             *
*   prefactor (inp) - exp[x^2-zeta_a*R_ac^2]                           *
*   tol (inp) - maximum relative error in bessel functions             *
*   n_ab (inp) - n_prim_a*n_prim_b                                     *
*   n_abc (inp) - n_prim_a*n_prim_b*n_c                                *
*   n_int - number of integral blocks to be generated                  *
*   temp - work array                                                  *
*   ind - index array                                                  *
*   Qabc - uncontracted Q integrals                                    *
*   Qa - Qabc integrals contracted over core potential, used in down-  *
*        ward recursion                                                *
*   Q_int - final Q integrals                                          *
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
      integer i,ibug,ind_l,ind_p,j,l,l_max,
     &    m,mm,m_hi,m_count,n_prim_a,n_prim_b,
     &    n_ab,n_abc,n_abp,n_c,n_int,n_p,
     &    p,p_max,p_min
      integer ind(n_abc),n_prim_c(p_min:p_max)
      logical debug_gen,debug_addresses,debug_arrays
      double precision bi(n_abc),a_sqrt(n_ab,n_c),x(n_ab,n_c),
     &    prefactor(n_ab,n_c),temp(n_abc,6),Qabc(n_ab,n_c),
     &    Qa(n_ab,0:l_max),Q_int(n_prim_a*n_prim_b,n_int),
     &    coeff_c(n_c),tol,wa,wb,wc,prefac
*
Cbug
C      debug_gen = .true.
C      debug_arrays = .true.
C      debug_addresses = .true.
Cstd
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
Cend
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_locint1 ...'
*
*   Set up initial values Q_L^{L+p} for L = l_min to l_max
*
      m_count = l_max/2
      if (debug_addresses) write (LuOut,*) 'm_count,l_max',
     &    m_count,l_max
      ind_l = 1
      do l = 0,l_max
        ind_p = 1
        if (debug_gen) write (LuOut,*) 'l =',l
        do p = p_min,p_max
          n_p = n_prim_c(p)
          n_abp = n_ab*n_p
          if (n_p .gt. 0) then
            if (p .eq. 4) then
              call ecp_t1_Q4int (n_abp,l,x(1,ind_p),a_sqrt(1,ind_p),
     &            Qabc(1,ind_p))
            else if (p .eq. 3) then
              call ecp_t1_Q3int (n_abp,l,x(1,ind_p),a_sqrt(1,ind_p),
     &            Qabc(1,ind_p),temp,ind,tol)
            else if (p .eq. 2) then
              call ecp_t1_Q2int (n_abp,l,x(1,ind_p),a_sqrt(1,ind_p),
     &            Qabc(1,ind_p))
            else if (p .eq. 1) then
              call ecp_t1_Q1int (n_abp,l,x(1,ind_p),a_sqrt(1,ind_p),
     &            Qabc(1,ind_p),temp,ind,tol)
            else if (p .eq. 0) then 
              call ecp_t1_Q0int (n_abp,l,x(1,ind_p),a_sqrt(1,ind_p),
     &            Qabc(1,ind_p),temp,ind,tol)
            else
              call errquit(
     &            'Illegal p value in routine ecp_locint1',99,
     &       UNKNOWN_ERR)
            end if
          end if
          ind_p = ind_p+n_p
        end do
        do i = 1,n_c
          do j = 1,n_ab
            Qabc(j,i) = Qabc(j,i)*prefactor(j,i)
          end do
        end do
        if (debug_arrays) call ecp_matpr (Qabc,1,n_ab,
     &      1,n_c,1,n_ab,1,n_c,'Primitive integrals',
     &      'E',78,4)
        call ecp_contract (n_ab,n_c,1,Qabc,coeff_c,Qa(1,l))
        if (debug_arrays) call ecp_matpr (Qa(1,l),1,n_prim_a,
     &      1,n_prim_b,1,n_prim_a,1,n_prim_b,'C-contracted integrals',
     &      'E',78,4)
        call dcopy (n_ab,Qa(1,l),1,Q_int(1,ind_l),1)
        ind_l = ind_l+l/2+1
      end do
*
*   Perform recursion to obtain values Q_{L-2n}^{L+p}
*
      m_hi = l_max-1
      ind_l = 1
      if (debug_addresses) write (LuOut,*) 'm_count',m_count
      mm = 2
      do m = 1,m_count
        if (debug_addresses) write (LuOut,*) 'm,mm,m_hi',m,mm,m_hi
        if (debug_arrays) call ecp_matpr (Qa(1,0),1,n_ab,0,
     &      l_max,1,n_ab,mm-1,l_max,'Qa before recursion','E',81,5)
        call ecp_down_m (1,m_hi,n_ab,bi,Qa(1,mm),Qa(1,mm-1),
     &      Qa(1,mm))
        if (debug_arrays) call ecp_matpr (Qa(1,0),1,n_ab,0,
     &      l_max,1,n_ab,mm,l_max,'Qa after recursion','E',81,5)
        ind_l = ind_l+m+m+1
        i = ind_l
        do l = mm,l_max
          if (debug_addresses) write (LuOut,*) 'm,l,ind_l,i',
     &        mm,l,ind_l,i
          call dcopy (n_ab,Qa(1,l),1,Q_int(1,i),1)
          i = i+l/2+1
        end do
        m_hi = m_hi-2
        mm = mm+2
      end do
*
*   Special code if beta = 0
*
      do j = 1,n_ab
        if (bi(j) .eq. zero) then
          do l = 2,n_int
            Q_int(j,l) = zero
          end do
*
          ind_p = 0
          do p = p_min,p_max
            n_p = n_prim_c(p)
            if (n_p .gt. 0) then
              if (p .eq. 0) then
                prefac = qpirt+qpirt
              else if (p .eq. 1) then
                prefac = half
              else if (p .eq. 2) then
                prefac = qpirt
              end if
              do i = ind_p+1,ind_p+n_p
                wa = p+1
                wa = wa*half
                wb = prefac*prefactor(j,i)*a_sqrt(j,i)**(p+1)
                wc = a_sqrt(j,i)**2
                do m = 1,m_count
                  l = (m+1)**2
                  wb = wb*wa*wc
                  Q_int(j,l) = Q_int(j,l)+wb*coeff_c(i)
                  wa = wa+one
                end do
              end do
            end if
            ind_p = ind_p+n_p
          end do
        end if
      end do
*
      if (debug_arrays) call ecp_matpr (Q_int(1,1),1,n_ab,1,n_int,
     &    1,n_ab,1,n_int,'ECP radial integrals','E',81,5)
      if (debug_gen) write (LuOut,*) 'Exiting ecp_locint1'
*
      return
      end
