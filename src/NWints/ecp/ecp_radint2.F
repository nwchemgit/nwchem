C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_radint2 (p_min,p_max,
     &    l_a,n_prim_a,n_cont_a,n_int_a,coeff_a,ai,alpha,
     &    l_b,n_prim_b,n_cont_b,n_int_b,coeff_b,bi,beta,
     &    l_c,n_prim_c,n_cont_c,n_c,coeff_c,gamma,prefactor,tol,sphcart,
     &    n_ab,n_abc,n_Q,m_count,temp,Qabc,Qab,Q,Q_int,ibug)
*                                                                      *
*   Calculate Type 2 radial integrals for a given ECP centre, angular  *
*   projector and exponent p                                           *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   p_min - minimum power of r in ECP expansion                        *
*   p_max - maximum power of r in ECP expansion                        *
*   l_a (inp) - (maximum) angular momentum of functions on centre a    *
*   n_prim_a (inp) - number of primitive functions on centre a         *
*   n_cont_a (inp) - number of contracted functions on centre a        *
*   n_int_a - total number of m values for centre a                    *
*   coeff_a - contraction coefficients of basis functions on centre a  *
*   ai (inp) - values of 1/2*zeta_a*R_ac                               *
*   alpha (inp) - values of a/2*sqrt(c)                                *
*   l_b (inp) - (maximum) angular momentum of functions on centre b    *
*   n_prim_b (inp) - number of primitive functions on centre b         *
*   n_cont_b (inp) - number of contracted functions on centre b        *
*   n_int_b - total number of m values for centre b                    *
*   coeff_b - contraction coefficients of basis functions on centre b  *
*   bi (inp) - values of 1/2*zeta_b*R_bc                               *
*   beta (inp) - values of b/2*sqrt(c)                                 *
*   l_c (inp) - angular momentum of ECP projector on centre c          *
*   n_c (inp) - total number of ECP primitive functions                *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   coeff_c - contraction coefficients of potential on centre c        *
*   gamma (inp) - values of 1/sqrt(c)                                  *
*   prefactor (inp) - exp[(alpha+beta)^2-zeta_a*R_ac^2-zeta_b*R_bc^2]  *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart - 1 for spherical basis, 0 for cartesian basis.            *
*   n_ab (inp) - n_prim_a*n_prim_b                                     *
*   n_abc (inp) - n_prim_a*n_prim_b*n_c                                *
*   n_Q - dimension of intermediate array Q                            *
*   m_count (inp) - number of Q^2_{mm} to be generated                 *
*   temp - work array                                                  *
*   Qabc - uncontracted Q integrals used in upward recursion in k      *
*   Qab - Q integrals contracted over core potential with maximum m, n *
*         values                                                       *
*   Q - work array used in downward recursion                          *
*   Q_int - final fully contracted Q integrals                         *
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
      integer h,i,j,ibug,ind_a,ind_b,ind_p,ind_even,ind_odd,i_c,
     &    ind_max,ind_mid,ind_min,ind_hi,ind_mi,ind_lo,
     &    k,k_max,k_p,
     &    l_a,l_b,l_c,l_max,l_mid,l_min,
     &    l,lambda_a,lambda_b,la_min,lb_min,
     &    m,m_aeqb,m_count,m_max,m_min,m_min1,
     &    mc_a,mc_aeqb,mc_b,mc_even,mc_odd,
     &    mt_max,mt_mid,mt_min,
     &    n,n_c,n_cont_c,n_ab,n_abc,n_abp,n_Q,n_max,n_abm,
     &    n_prim_a,n_prim_b,n_cont_a,n_cont_b,n_int_a,n_int_b,
     &    p,p_max,p_min,sphcart
      integer n_prim_c(p_min:p_max)
      logical debug_gen,debug_addresses,debug_arrays
      double precision ai(n_abc),bi(n_abc),coeff_c(n_c),
     &    alpha(n_abc),beta(n_abc),gamma(n_abc),prefactor(n_abc),
     &    coeff_a(n_prim_a,n_cont_a),coeff_b(n_prim_b,n_cont_b),
     &    Qabc(n_abc,m_count,3),Qab(n_ab,0:l_a,0:l_b,n_cont_c),
     &    Q(n_ab,n_Q,n_cont_c),temp(n_abc,20),
     &    Q_int(n_cont_a*n_cont_b,n_int_a,n_int_b,n_cont_c),
     &    tol
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
*
      if (debug_gen) then
        write (LuOut,'(//A,/)') 'Entering ecp_radint2 ...'
        write (LuOut,*) 'ibug =',ibug
      end if
      if (debug_addresses) then
	write (LuOut,*) 'l_a,l_b,l_c',l_a,l_b,l_c
        write (LuOut,*) 'n_ab,n_abc,n_Q,m_count',
     &      n_ab,n_abc,n_Q,m_count
      end if
      k_max = l_a+l_b
      l_max = max(l_a,l_b)
      l_min = min(l_a,l_b)
      if (debug_addresses) write (LuOut,*) 'k_max,l_max,l_min',
     &    k_max,l_max,l_min
      mt_max = l_max*(l_max+1)/2
      mt_min = l_min*(l_min-1)/2
      mt_mid = mt_max
      if (l_a .eq. l_b) then
        m_aeqb = 0
      else
        mt_max = mt_max-l_max
        m_aeqb = 1
      end if
      if (debug_addresses) write (LuOut,*) 'mt_max,mt_mid,mt_min',
     &    mt_max,mt_mid,mt_min
*
      ind_max = 1
      ind_mid = ind_max+mt_max
      ind_min = ind_mid+mt_mid
      if (debug_addresses) write (LuOut,*) 'ind_max,ind_mid,ind_min',
     &    ind_max,ind_mid,ind_min
      if (l_a .ge. l_b) then
        j = 1
        h = 0
      else
        j = 0
        h = 1
      end if
      if (debug_addresses) write (LuOut,*) 'j,h',j,h
*
*   Zero out arrays for accumulation
*
      n_abm = n_ab*(mt_min+mt_mid+mt_max)
      call dcopy(n_abm,zero,0,Q,1)
      n_abm = (l_a+1)*(l_b+1)*n_ab*n_cont_c
      call dcopy(n_abm,zero,0,Qab,1)
*
*   Set up initial values for K = 0 and 1
*
      m_min = l_c
      m_max = m_min+m_count-1
      if (debug_addresses) write (LuOut,*) 'm_count,m_min,m_max',
     &    m_count,m_min,m_max
*
*   Loop over p values to generate primitive integrals
*
      ind_p = 1
      if (debug_addresses) write (LuOut,*) 'n_abc',n_abc
      do p = p_min,p_max
        n_abp = n_ab*n_prim_c(p)
        if (debug_addresses) write (LuOut,*) 'p =',p
        if (debug_addresses) write (LuOut,*) n_abp,n_ab,n_prim_c(p)
        if (p .eq. 4) then
          if (n_abp .gt. 0) call ecp_t2_init4 (n_abp,n_abc,m_min,m_max,
     &        j,h,tol,ai,bi,alpha(ind_p),beta(ind_p),gamma(ind_p),
     &        prefactor(ind_p),temp(1,2),temp(1,1),Qabc(ind_p,1,1))
        else if (p .eq. 3) then
          if (n_abp .gt. 0) call ecp_t2_init3 (n_abp,n_abc,m_min,m_max,
     &        j,h,tol,ai,bi,alpha(ind_p),beta(ind_p),gamma(ind_p),
     &        prefactor(ind_p),temp(1,2),temp(1,1),Qabc(ind_p,1,1))
        else if (p .eq. 2) then
          if (n_abp .gt. 0) call ecp_t2_init2 (n_abp,n_abc,m_min,m_max,
     &        h,tol,alpha(ind_p),beta(ind_p),gamma(ind_p),
     &        prefactor(ind_p),temp,temp(1,6),temp(1,7),Qabc(ind_p,1,1))
        else if (p .eq. 1) then
          if (n_abp .gt. 0) call ecp_t2_init1 (n_abp,n_abc,m_min,m_max,
     &        j,h,tol,ai,bi,alpha(ind_p),beta(ind_p),gamma(ind_p),
     &        prefactor(ind_p),temp(1,2),temp(1,1),Qabc(ind_p,1,1))
        else if (p .eq. 0) then
          if (n_abp .gt. 0) call ecp_t2_init0 (n_abp,n_abc,m_min,m_max,
     &        j,h,tol,ai,bi,alpha(ind_p),beta(ind_p),gamma(ind_p),
     &        prefactor(ind_p),temp(1,2),temp(1,1),Qabc(ind_p,1,1))
        else
          call errquit(
     &        'Illegal p value in routine ecp_radint2',
     &        p, BASIS_ERR)
        end if
        if (debug_addresses) write (LuOut,*) 'end p loop'
        ind_p = ind_p+n_abp
      end do
      if (debug_arrays) then
        call ecp_matpr (Qabc(1,1,1),1,n_abc,m_min,m_max,
     &      1,n_abc,m_min,m_max,'Qabc(i,m,1)','E',78,4)
        call ecp_matpr (Qabc(1,1,2),1,n_abc,m_min+1,m_max,
     &      1,n_abc,m_min+1,m_max,'Qabc(i,m,2)','E',78,4)
      end if
*
*   Contract over ECP functions for lambda_a, lambda_b = 0,0
*
      do i_c = 1,n_cont_c
        call ecp_contract (n_ab,n_c,1,Qabc,coeff_c,Qab(1,0,0,i_c))
        if (debug_arrays) call ecp_matpr (Qab(1,0,0,i_c),1,n_prim_a,
     &      1,n_prim_b,1,n_prim_a,1,n_prim_b,'Qab(i,0,0)','E',78,4)
        if (k_max .ge. 2) call ecp_contract (n_ab,n_c,l_max,Qabc(1,2,1),
     &      coeff_c,Q(1,ind_mid,i_c))
*
*   Contract over ECP functions for lambda_a, lambda_b = 0,1 or 1,0
*
        if (k_max .gt. 0) then
          call ecp_contract (n_ab,n_c,1,Qabc(1,1,2),coeff_c,
     &        Qab(1,j,h,i_c))
          if (debug_arrays) then
            if (j .gt. h) then
              call ecp_matpr (Qab(1,j,h,i_c),1,n_prim_a,1,n_prim_b,
     &            1,n_prim_a,1,n_prim_b,'Qab(i,1,0)','E',78,4)
            else
              call ecp_matpr (Qab(1,j,h,i_c),1,n_prim_a,1,n_prim_b,
     &            1,n_prim_a,1,n_prim_b,'Qab(i,0,1)','E',78,4)
            end if
          end if
          if (k_max .ge. 2) call ecp_contract (n_ab,n_c,l_max-m_aeqb,
     &        Qabc(1,2,2),coeff_c,Q(1,ind_max,i_c))
          if (debug_arrays) then
            if (mt_max .gt. 0) write (LuOut,'(1p5e15.7)') 
     &          (Q(1,ind_max+i,i_c),i=0,mt_max-1)
            if (mt_mid .gt. 0) write (LuOut,'(1p5e15.7)') 
     &          (Q(1,ind_mid+i,i_c),i=0,mt_mid-1)
            if (mt_min .gt. 0) write (LuOut,'(1p5e15.7)') 
     &          (Q(1,ind_min+i,i_c),i=0,mt_min-1)
          end if
        end if
      end do
*
*   Upward recursion in k to maximum. Recursion is performed for
*   m if l_a > l_b otherwise for n.
*
      if (debug_gen) write (LuOut,'(/A,/)') 'Starting upward recursion'
      ind_lo = 1
      ind_mi = 2
      ind_hi = 3
      mc_odd = l_max-m_aeqb
      mc_even = l_max-1
      ind_odd = ind_max+mc_odd
      mc_odd = mc_odd-1
      ind_even = ind_mid+l_max
*
      do k_p = 2,k_max
        l_mid = k_p/2
        if (debug_gen) write (LuOut,*) 'k_p,l_mid',k_p,l_mid
        if (mod(k_p,2) .eq. 0) then
*     
*     Even k recursion to obtain Q^k_{mm}
*
          m_max = m_max-1
          m_min = m_min+1
          ind_p = 1
          if (debug_addresses) write (LuOut,*) 'k,m_min,m_max',
     &        k_p,m_min,m_max
          do p = p_min,p_max
            n_abp = n_ab*n_prim_c(p)
            k = k_p+p-2
            if (debug_arrays) then
              write(LuOut,'(1p5e15.7)') (Qabc(1,i+2,ind_lo),i=0,
     &            m_max-m_min)
              write(LuOut,'(1p5e15.7)') (Qabc(1,i+1,ind_mi),i=0,
     &            m_max-m_min)
            end if
            call ecp_up_k (m_min,m_max,k,0,j,n_abp,n_abc,
     &          alpha(ind_p),beta(ind_p),gamma(ind_p),
     &          Qabc(ind_p,2,ind_lo),Qabc(ind_p,1,ind_mi),
     &          Qabc(ind_p,1,ind_hi))
            if (debug_arrays) write (LuOut,'(1p5e15.7)') 
     &          (Qabc(1,i+1,ind_hi),i=0,m_max-m_min)
            ind_p = ind_p+n_abp
          end do
*
*       Contract over ECP functions
*
          if (debug_addresses) then
            write (LuOut,*) 'ind_hi',ind_hi
            write (LuOut,*) 'ind_even,mc_even',ind_even,mc_even
          end if
          do i_c = 1,n_cont_c
            if (mc_even .gt. 0) call ecp_contract (n_ab,n_c,mc_even,
     &          Qabc(1,2,ind_hi),coeff_c,Q(1,ind_even,i_c))
            if ((l_mid .le. l_a) .and. (l_mid .le. l_b)) then
              call ecp_contract (n_ab,n_c,1,Qabc(1,1,ind_hi),
     &            coeff_c,Qab(1,l_mid,l_mid,i_c))
              if (debug_arrays) call ecp_matpr (Qab(1,l_mid,l_mid,i_c),
     &            1,n_prim_a,1,n_prim_b,1,n_prim_a,1,n_prim_b,
     &            'Qab(*,l_mid,l_mid)','E',78,4)
            end if
          end do
          ind_even = ind_even+mc_even
          mc_even = mc_even-1
        else
*
*     Odd k recursion to obtain Q^k_{m+1m} or Q^k_{mm+1}
*
          ind_p = 1
          do p = p_min,p_max
            n_abp = n_ab*n_prim_c(p)
            k = k_p+p-2
            call ecp_up_k (m_min,m_max,k,1,h,n_abp,n_abc,
     &          alpha(ind_p),beta(ind_p),gamma(ind_p),
     &          Qabc(ind_p,2,ind_lo),Qabc(ind_p,1,ind_mi),
     &          Qabc(ind_p,1,ind_hi))
            ind_p = ind_p+n_abp
          end do
*
*       Contract over ECP functions
*
          do i_c = 1,n_cont_c
            if (mc_odd .gt. 0) call ecp_contract (n_ab,n_c,mc_odd,
     &          Qabc(1,2,ind_hi),coeff_c,Q(1,ind_odd,i_c))
            if ((l_mid+j .le. l_a) .and. (l_mid+h .le. l_b)) then
              call ecp_contract (n_ab,n_c,1,Qabc(1,1,ind_hi),
     &            coeff_c,Qab(1,l_mid+j,l_mid+h,i_c))
              if (debug_arrays) call ecp_matpr (Qab(1,l_mid+j,l_mid+h,
     &            i_c),1,n_prim_a,1,n_prim_b,1,n_prim_a,1,n_prim_b,
     &            'Qab(*,l_mid+j,l_mid+h)','E',78,4)
            end if
          end do
          ind_odd = ind_odd+mc_odd
          mc_odd = mc_odd-1
        end if
        i = ind_lo
        ind_lo = ind_mi
        ind_mi = ind_hi
        ind_hi = i
      end do
*
*   Downward recursions in m,n
*
      m_min = l_c
      ind_odd = ind_max
      ind_even = ind_mid
      if (l_a .ge. l_b) then
        ind_a = ind_max
        ind_b = ind_min
      else
        ind_a = ind_min
        ind_b = ind_max
      end if
      mc_a = l_a-1
      mc_aeqb = l_a-m_aeqb
      mc_b = l_b-1
      mc_even = l_max
      if (debug_gen) write(LuOut,'(/A,/)') 'Starting downward recursion'
      if (debug_addresses) then
        write (LuOut,*) 'ind_max,ind_mid,ind_min',
     &      ind_max,ind_mid,ind_min
        write (LuOut,*) 'mt_max,mt_mid,mt_min',mt_max,mt_mid,mt_min
      end if
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_max+i),i=0,mt_max-1)
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_mid+i),i=0,mt_mid-1)
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_min+i),i=0,mt_min-1)
      do k_p = 1,k_max
        l_mid = k_p/2
        if (debug_gen) write (LuOut,*) 'k_p,l_mid',k_p,l_mid
*     
*     Recursion from Q^k_{mm} for even k
*
        if (mod(k_p,2) .eq. 0) then
          m_min = m_min+1
*
*       Shift angular momentum from centre b to centre a, by
*       downward recursion on m_b
*
          m_max = m_min+mc_a-1
          n_max = min(l_mid,l_a-l_mid)-1
          ind_hi = ind_even
          ind_lo = ind_a
          if (debug_addresses) 
     &        write (LuOut,*) 'k,m_min,m_max',k_p,m_min,m_max
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_max+i),i=0,mt_max-1)
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_mid+i),i=0,mt_mid-1)
C          write (LuOut,'(1p5e15.7)') (Q(1,ind_min+i),i=0,mt_min-1)
          do n = 0,n_max
            lambda_a = l_mid+n+1
            lambda_b = l_mid-n-1
            do i_c = 1,n_cont_c
              call ecp_down_m (m_min-n,m_max-2*n,n_ab,bi,
     &            Q(1,ind_hi,i_c),Q(1,ind_lo,i_c),Q(1,ind_lo,i_c))
              if ((lambda_a .le. l_a) .and. (lambda_b .le. l_b))
     &            call dcopy (n_ab,Q(1,ind_lo,i_c),1,
     &            Qab(1,lambda_a,lambda_b,i_c),1)
            end do
            ind_hi = ind_lo+1
            ind_lo = ind_hi-mc_aeqb-n
          end do
*
*       .. and shift angular momentum from centre a to centre b, by
*       downward recursion on m_a
*
          m_max = m_min+mc_b-1
          n_max = min(l_mid,l_b-l_mid)-1
          ind_hi = ind_even
          ind_lo = ind_b
          do n = 0,n_max
            lambda_a = l_mid-n-1
            lambda_b = l_mid+n+1
            do i_c = 1,n_cont_c
              call ecp_down_m (m_min-n,m_max-2*n,n_ab,ai,
     &            Q(1,ind_hi,i_c),Q(1,ind_lo,i_c),Q(1,ind_lo,i_c))
              if ((lambda_a .le. l_a) .and. (lambda_b .le. l_b))
     &            call dcopy (n_ab,Q(1,ind_lo,i_c),1,
     &            Qab(1,lambda_a,lambda_b,i_c),1)
            end do
            ind_hi = ind_lo+1
            ind_lo = ind_hi-mc_b-n
          end do
          ind_a = ind_a+mc_aeqb
          ind_b = ind_b+mc_b
          mc_a = mc_a-1
          mc_aeqb = mc_aeqb-1
          mc_b = mc_b-1
        else
*
*     Odd k : generate Q_mm+1 from Qm+1m or vice-versa, depending on
*     whether l_a > l_b or v.v.
*
          m_min1 = m_min+1
          if (l_a .ge. l_b) then
            m_max = m_min1+mc_b
            if (debug_addresses) then
              write (LuOut,*) 'ind_a+1,ind_even+1,ind_b',
     &            ind_a+1,ind_even+1,ind_b
              write (LuOut,*) 'm_min1,m_max',m_min1,m_max
            end if
            do i_c = 1,n_cont_c
              if (m_min1 .lt. m_max) call ecp_down_m (m_min1+1,m_max,
     &            n_ab,ai,Q(1,ind_a+1,i_c),Q(1,ind_even+1,i_c),
     &            Q(1,ind_b,i_c))
              if ((l_mid .le. l_a) .and. (l_mid+1 .le. l_b))
     &            call ecp_down_m (m_min1,m_min1,n_ab,ai,
     &            Q(1,ind_a,i_c),Q(1,ind_even,i_c),
     &            Qab(1,l_mid,l_mid+1,i_c))
            end do
          else
            m_max = m_min1+mc_a
            if (debug_addresses) then
              write (LuOut,*) 'ind_b+1,ind_even+1,ind_a',
     &            ind_b+1,ind_even+1,ind_a
              write (LuOut,*) 'm_min1,m_max',m_min1,m_max
            end if
            do i_c = 1,n_cont_c
              if (m_min1 .lt. m_max) call ecp_down_m (m_min1+1,m_max,
     &            n_ab,bi,Q(1,ind_b+1,i_c),Q(1,ind_even+1,i_c),
     &            Q(1,ind_a,i_c))
              if ((l_mid+1 .le. l_a) .and. (l_mid .le. l_b))
     &            call ecp_down_m (m_min1,m_min1,n_ab,bi,
     &            Q(1,ind_b,i_c),Q(1,ind_even,i_c),
     &            Qab(1,l_mid+1,l_mid,i_c))
            end do
          end if
          ind_even = ind_even+mc_even
          mc_even = mc_even-1
*     
*       Shift angular momentum from centre b to centre a, by
*       downward recursion on m_b
*
          m_max = m_min+mc_a-1
          n_max = min(l_mid,l_a-l_mid-1)-1
          ind_hi = ind_a
          ind_lo = ind_a-mc_aeqb
          if (debug_addresses) 
     &        write (LuOut,*) 'k,m_min,m_max',k_p,m_min,m_max
          do n = 0,n_max
            lambda_a = l_mid+n+2
            lambda_b = l_mid-n-1
            do i_c = 1,n_cont_c
              call ecp_down_m (m_min-n,m_max-2*n,n_ab,bi,
     &            Q(1,ind_hi,i_c),Q(1,ind_lo,i_c),Q(1,ind_lo,i_c))
              if ((lambda_a .le. l_a) .and. (lambda_b .le. l_b))
     &            call dcopy (n_ab,Q(1,ind_lo,i_c),1,
     &            Qab(1,lambda_a,lambda_b,i_c),1)
            end do
            ind_hi = ind_lo+1
            ind_lo = ind_lo-mc_aeqb-n
          end do
*
*       .. and shift angular momentum from centre a to centre b, by
*       downward recursion on m_a
*
          m_max = m_min+mc_b-1
          n_max = min(l_mid,l_b-l_mid-1)-1
          ind_hi = ind_b
          ind_lo = ind_b-mc_b
          do n = 0,n_max
            lambda_a = l_mid-n-1
            lambda_b = l_mid+n+2
            do i_c = 1,n_cont_c
              call ecp_down_m (m_min-n,m_max-2*n,n_ab,ai,
     &            Q(1,ind_hi,i_c),Q(1,ind_lo,i_c),Q(1,ind_lo,i_c))
              if ((lambda_a .le. l_a) .and. (lambda_b .le. l_b))
     &            call dcopy (n_ab,Q(1,ind_lo,i_c),1,
     &            Qab(1,lambda_a,lambda_b,i_c),1)
            end do
            ind_hi = ind_lo+1
            ind_lo = ind_lo-mc_b-n
          end do
        end if
      end do
*
*    At this point, the highest m,n values in each block have been
*    generated. Now use downward recursion again to get all members of
*    the block. This has to be done in a different order because the
*    final result is indexed by m and n, not by (m-n).
*
      do i_c = 1,n_cont_c
        ind_b = 0
        do lambda_b = 0,l_b
          if (sphcart .eq. 0) then
            lb_min = max(0,lambda_b-(lambda_b+l_c)/2)
          else
            lb_min = max(0,lambda_b-l_c)
          end if
*
*     First do recursion in b to generate top row for all a.
*
          do lambda_a = 0,l_a
            m = lambda_b+l_c+1
            do l = lambda_b-1,lb_min,-1
              m = m-2
              call ecp_down_m (m,m,n_ab,bi,
     &            Qab(1,lambda_a,l+1,i_c),Qab(1,lambda_a,l,i_c),
     &            Qab(1,lambda_a,l,i_c))
            end do
          end do
*
*     Now contract over b and do recursion over a to generate final 
*     integrals, contracting over a.
*
          do l = lambda_b,lb_min,-1
            ind_b = ind_b+1
            ind_a = 0
            do lambda_a = 0,l_a
              if (sphcart .eq. 0) then
                la_min = max(0,lambda_a-(lambda_a+l_c)/2)
              else
                la_min = max(0,lambda_a-l_c)
              end if
              if (debug_addresses) write (LuOut,*) 'lambda_a,l',
     &            lambda_a,l
              call dgemm ('N','N',n_prim_a,n_cont_b,n_prim_b,one,
     &            Qab(1,lambda_a,l,i_c),n_prim_a,coeff_b,n_prim_b,zero,
     &            Q(1,lambda_a+1,i_c),n_prim_a)
              ind_a = ind_a+1
              call dgemm ('T','N',n_cont_a,n_cont_b,n_prim_a,one,
     &            coeff_a,n_prim_a,Q(1,lambda_a+1,i_c),n_prim_a,zero,
     &            Q_int(1,ind_a,ind_b,i_c),n_cont_a)
              if (debug_gen) 
     &            write (LuOut,'(A,2I5)') 'ind_a,ind_b',ind_a,ind_b
              if (debug_arrays) call ecp_matpr(Q_int(1,ind_a,ind_b,i_c),
     &            1,n_cont_a,1,n_cont_b,1,n_cont_a,1,n_cont_b,
     &            'Q_int','E',78,4)
              m = lambda_a+l_c+1
              do k = lambda_a,la_min+1,-1
                m = m-2
                call ecp_down_m (m,m,n_prim_a*n_cont_b,
     &              ai,Q(1,k+1,i_c),Q(1,k,i_c),Q(1,k,i_c))
                ind_a = ind_a+1
                call dgemm ('T','N',n_cont_a,n_cont_b,n_prim_a,one,
     &              coeff_a,n_prim_a,Q(1,k,i_c),n_prim_a,zero,
     &              Q_int(1,ind_a,ind_b,i_c),n_cont_a)
                if (debug_gen) 
     &              write (LuOut,'(A,2I5)') 'ind_a,ind_b',ind_a,ind_b
                if (debug_arrays) call ecp_matpr (Q_int(1,ind_a,ind_b,
     &              i_c),1,n_cont_a,1,n_cont_b,1,n_cont_a,1,n_cont_b,
     &              'Q_int','E',78,4)
              end do
            end do
          end do
        end do
      end do
      if (debug_gen) write (LuOut,'(//A,/)') 'Exiting ecp_radint2'
*
      return
      end
