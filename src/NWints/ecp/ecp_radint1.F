C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_radint1 (p_min,p_max,l_lo,l_min,l_max,m_count,
     &    m_skip,n_prim_a,n_cont_a,coeff_a,n_prim_b,n_cont_b,coeff_b,
     &    n_c,n_prim_c,n_cont_c,coeff_c,ai,a_sqrt,x,prefactor,tol,
     &    sphcart,n_ab,nc_ab,n_abc,n_int,temp,ind,work,Qabc,Qa,Q_int,
     &    transpose,ibug)
*                                                                      *
*   Calculate Type 1 radial integrals Q^k_m for a given ECP centre and *
*   angular projector                                                  *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   l_lo (inp) -  minimum k value for initial integrals Q^k_k          *
*   l_min (inp) -  minimum k value = l_b                               *
*   l_max (inp) -  maximum k value = l_b+l_a                           *
*   m_count (inp) - number of downward recursions required             *
*   m_skip (inp) - number of rows of m values to skip before storing   *
*                  integrals in output array                           *
*   n_prim_a (inp) - number of primitive functions on centre a         *
*   n_cont_a (inp) - number of contracted functions on centre a        *
*   coeff_a - contraction coefficients of basis functions on centre a  *
*   n_prim_b (inp) - number of primitive functions on centre b         *
*   n_cont_b (inp) - number of contracted functions on centre b        *
*   coeff_b - contraction coefficients of basis functions on centre b  *
*   n_c (inp) - total number of ECP primitive functions                *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_cont_c (inp) - number of ECP contractions                        *
*   coeff_c - contraction coefficients of potential on centre c        *
*   ai (inp) - values of 1/alpha                                       *
*   a_sqrt (inp) - values of sqrt(ai)                                  *
*   x (inp) - values of alpha/2*sqrt(beta)                             *
*   prefactor (inp) - exp[x^2-zeta_a*R_ac^2]                           *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart - 0 for cartesians, 1 for sphericals                       *
*   n_ab (inp) - n_prim_a*n_prim_b                                     *
*   nc_ab (inp) - n_prim_a*n_cont_b                                    *
*   n_abc (inp) - n_prim_a*n_prim_b*n_c                                *
*   n_int - number of integral blocks to be generated                  *
*   temp - work array                                                  *
*   work - work array                                                  *
*   Qabc - uncontracted Q integrals                                    *
*   Qa - Qabc integrals contracted over core potential and on-centre   *
*        basis functions.                                              *
*   Q_int - final fully contracted Q integrals                         *
*   transpose - logical flag for transposition of the final integrals  *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
#include "ecp_consts.fh"
#include "errquit.fh"
      integer i,j,ibug,inc_l,ind_l,ind_p,l,ll,lm,l_lo,l_min,l_max,
     &    m,mm,m_hi,m_lo,m_count,m_skip,n_prim_a,n_prim_b,n_cont_a,
     &    n_cont_b,n_cont_c,n_ab,n_abc,n_abp,n_c,nc_ab,ncab,n_int,n_p,
     &    p,p_max,p_min,sphcart
      integer ind(n_abc),n_prim_c(p_min:p_max)
      logical transpose,debug_gen,debug_addresses,debug_arrays
      double precision ai(n_abc),a_sqrt(n_abc),x(n_abc),
     &    prefactor(n_abc),temp(n_abc,6),work(n_ab),Qabc(n_abc),
     &    Qa(nc_ab,l_lo:l_max,n_cont_c),
     &    Q_int(n_cont_a*n_cont_b,n_int,n_cont_c),
     &    coeff_a(n_prim_a,n_cont_a),coeff_b(n_prim_b,n_cont_b),
     &    coeff_c(n_c,n_cont_c),tol
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_radint1 ...'
*
*   Set up initial values Q_L^{L+p} for L = l_min to l_max
*
      if (debug_addresses) write (LuOut,*) 'm_count,m_skip',
     &    m_count,m_skip
      ind_l = 1
      inc_l = 0
      if (debug_addresses) write (LuOut,*) 'l_lo,l_min,l_max',
     &    l_lo,l_min,l_max
      do l = l_lo,l_max
        call dcopy (n_ab,zero,0,work,1)
        ind_p = 1
        if (debug_gen) write (LuOut,*) 'l =',l
        do p = p_min,p_max
          n_p = n_prim_c(p)
          n_abp = n_ab*n_p
          if (debug_addresses) write (LuOut,'(3I5)') p,n_p,n_abp
          if (n_abp .gt. 0) then
            if (p .eq. 4) then
              call ecp_t1_Q4int (n_abp,l,x(ind_p),a_sqrt(ind_p),
     &            Qabc(ind_p))
            else if (p .eq. 3) then
              call ecp_t1_Q3int (n_abp,l,x(ind_p),a_sqrt(ind_p),
     &            Qabc(ind_p),temp,ind,tol)
            else if (p .eq. 2) then
              call ecp_t1_Q2int (n_abp,l,x(ind_p),a_sqrt(ind_p),
     &            Qabc(ind_p))
            else if (p .eq. 1) then
              call ecp_t1_Q1int (n_abp,l,x(ind_p),a_sqrt(ind_p),
     &            Qabc(ind_p),temp,ind,tol)
            else if (p .eq. 0) then 
              call ecp_t1_Q0int (n_abp,l,x(ind_p),a_sqrt(ind_p),
     &            Qabc(ind_p),temp,ind,tol)
            else
              call errquit(
     &            'Illegal p value in routine ecp_radint1',99,
     &       BASIS_ERR)
            end if
          end if
          ind_p = ind_p+n_abp
        end do
        do i = 1,n_abc
          Qabc(i) = Qabc(i)*prefactor(i)
        end do
	do j = 1,n_cont_c
          call ecp_contract (n_ab,n_c,1,Qabc,coeff_c(1,j),work)
          if (debug_arrays) call ecp_matpr (work,1,n_prim_a,1,n_prim_b,
     &        1,n_prim_a,1,n_prim_b,'Primitive integrals','E',78,4)
*
*     Contract over on-centre basis functions
*
          call dgemm ('N','N',n_prim_a,n_cont_b,n_prim_b,one,work,
     &        n_prim_a,coeff_b,n_prim_b,zero,Qa(1,l,j),n_prim_a)
          if (debug_arrays)  call ecp_matpr (Qa(1,l,j),1,n_prim_a,
     &        1,n_cont_b,1,n_prim_a,1,n_cont_b,'Qa(*,l)','E',78,4)
        end do
*
*    Copy integrals into output array if l_c .ge. l_a
*    with contraction over off-center basis functions
*
        if ((l .ge. l_min) .and. (m_skip .eq. 0)) then
          do j = 1,n_cont_c
            if (transpose) then
              call dgemm ('T','N',n_cont_b,n_cont_a,n_prim_a,one,
     &            Qa(1,l,j),n_prim_a,coeff_a,n_prim_a,zero,
     &            Q_int(1,ind_l,j),n_cont_b)
              if (debug_arrays) then
                write (LuOut,*) 'ind_l ',ind_l
                call ecp_matpr (Q_int(1,ind_l,j),1,n_cont_b,1,n_cont_a,
     &              1,n_cont_b,1,n_cont_a,'Q_int(*,ind_l)','E',78,4)
              end if
            else
              call dgemm ('T','N',n_cont_a,n_cont_b,n_prim_a,one,
     &            coeff_a,n_prim_a,Qa(1,l,j),n_prim_a,zero,
     &            Q_int(1,ind_l,j),n_cont_a)
              if (debug_arrays) then
                write (LuOut,*) 'ind_l ',ind_l
                call ecp_matpr (Q_int(1,ind_l,j),1,n_cont_a,1,n_cont_b,
     &              1,n_cont_a,1,n_cont_b,'Q_int(*,ind_l)','E',78,4)
              end if
            end if
          end do
          if (sphcart .eq. 0) then
            inc_l = min(inc_l,l/2)+1
          else
            inc_l = min(inc_l,m_count)+1
          end if
          ind_l = ind_l+inc_l
        end if
      end do
      if (debug_arrays) then
        do j = 1,n_cont_c
          call ecp_matpr (Qa(1,1,j),1,nc_ab,l_lo,l_max,
     &        1,nc_ab,l_lo,l_max,'Qa integrals','E',81,5)
        end do
      end if
*
*   Perform recursion to obtain values Q_{L-2n}^{L+p},
*   contract over off-centre basis functions
*
      ncab = n_cont_a*n_cont_b
      do j = 1,n_cont_c
        m_hi = l_max-1
        m_lo = max(l_lo,1)
        inc_l = 1
        ind_l = 1
        if (m_skip .gt. 0) ind_l = 0
        if (debug_addresses) write (LuOut,*) 'm_count',m_count
        mm = m_lo+1
        do m = 1,m_count
          if (debug_addresses) write (LuOut,*) 'm,mm',m,mm
          if (debug_addresses) write (LuOut,*) 'm_lo,m_hi',m_lo,m_hi
          if (debug_arrays) call ecp_matpr (Qa(1,l_lo,j),1,nc_ab,l_lo,
     &        l_max,1,nc_ab,mm-1,l_max,'Qa before recursion','E',81,5)
          call ecp_down_m (m_lo,m_hi,nc_ab,ai,Qa(1,mm,j),Qa(1,mm-1,j),
     &        Qa(1,mm,j))
          if (debug_arrays) call ecp_matpr (Qa(1,l_lo,j),1,nc_ab,l_lo,
     &        l_max,1,nc_ab,mm,l_max,'Qa after recursion','E',81,5)
          if (debug_addresses) write (LuOut,*) 'ind_l,inc_l',ind_l,inc_l
          ll = m-m_skip
          if (ll .ge. 0) then
            ind_l = ind_l+ll+1
            if ((sphcart .eq. 0) .and. (m+m_skip .gt. l_min)) 
     &          ind_l = ind_l+ll
            i = ind_l
            do l = mm,l_max
              lm = l-l_min
              if (debug_addresses) write (LuOut,*) 'm,l,ind_l,lm,i',
     &            mm,l,ind_l,lm,i
              if (transpose) then
                call dgemm ('T','N',n_cont_b,n_cont_a,n_prim_a,one,
     &              Qa(1,l,j),n_prim_a,coeff_a,n_prim_a,zero,
     &              Q_int(1,i,j),n_cont_b)
              else
                call dgemm ('T','N',n_cont_a,n_cont_b,n_prim_a,one,
     &              coeff_a,n_prim_a,Qa(1,l,j),n_prim_a,zero,
     &              Q_int(1,i,j),n_cont_a)
              end if
              if (sphcart .eq. 0) then
                i = i+min(l/2-m_skip,lm)+1
              else
                i = i+min(lm,m_count)+1
              end if
            end do
          end if
          m_lo = max(m_lo-1,1)
          m_hi = m_hi-2
          mm = mm+1
          if (m .ge. l_lo) mm = mm+1
        end do
        if (debug_arrays) call ecp_matpr (Q_int(1,1,j),1,ncab,1,n_int,
     &      1,ncab,1,n_int,'ECP radial integrals','E',81,5)
        if (debug_gen) write (LuOut,*) 'Exiting ecp_radint1'
      end do
*
      return
      end
