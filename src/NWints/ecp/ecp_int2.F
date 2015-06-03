C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_int2 (mem_max,DryRun,
     &    R_AC,X_AC,Y_AC,Z_AC,l_a,n_prim_a,n_cont_a,coef_a,zeta_a,
     &    R_BC,X_BC,Y_BC,Z_BC,l_b,n_prim_b,n_cont_b,coef_b,zeta_b,
     &    l_ecp,n_prim_c,n_coef_c,
     &    ind_z,ind_c,n_z_tot,n_c_tot,l_ecp_max,n_cent,
     &    zeta_c,coef_c,p_min,p_max,tol,sphcart,tmp,ltmp,
     &    csco,lcsco,
     &    ecp_ints,n_ints,n_blk,ibug)
*                                                                      *
*   Calculate Type 2 radial integrals for a given ECP centre, angular  *
*   projector and exponent p                                           *
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
*   coef_a (inp) - centre A contraction coefficients                   *
*   zeta_a (inp) - centre A exponents                                  *
*   R_BC (inp) - distance between centres B and C                      *
*   X_BC,Y_BC,Z_BC (inp) - cartesian coordinates of centre C relative  *
*                          to centre B, X_BC = X_C - X_B, etc.         *
*   l_b (inp) - (maximum) angular momentum of functions on centre B    *
*   n_prim_b (inp) - number of primitive functions on centre B         *
*   n_cont_b (inp) - number of contracted functions on centre B        *
*   coef_b (inp) - centre B contraction coefficients                   *
*   zeta_b (inp) - centre B exponents                                  *
*   l_ecp_max (inp) - maximum angular momentum of ECP projector        *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_coef_c (inp) - array of number of coefficients/exponents for     *
*                    each projector                                    *
*   ind_z (inp) - array of addresses of first exponent for each l      *
*   ind_c (inp) - array of addresses of first coefficient for each l   *
*   n_z_tot (inp) - total number of ECP exponents                      *
*   n_c_tot (inp) - total number of ECP coefficients                   *
*   n_cent (inp) - number of ECP centres                               *
*   coef_c (inp) - ECP contraction coefficients                        *
*   zeta_c (inp) - ECP exponents                                       *
*   p_min (inp) - minimum power of r in ECP expansion                  *
*   p_max (inp) - maximum power of r in ECP expansion                  *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart (inp) - 1 for spherical basis, 0 for cartesian basis.      *
*   tmp (scr) - work array                                             *
*   ltmp - number of elements in tmp                                  *
*   csco (inp) - cartesian/spherical transformation                    *
*   lcsco (inp) - length of cartesian/spherical transformation         *
*   ecp_ints (out) - integrals over ECP                                *
*   n_ints (inp) - number of ecp radial integrals                      *
*   n_blk (inp) - 1 for scalar only, 3 for s-o only, 4 for both        *
*   ibug (inp) - debug flag. 0 for no debug, 1 for address printing,   *
*          2 for array printing, 3 for both.                           *
*                                                                      *
*   Notes:                                                             *
*   -----                                                              *
*                                                                      *
*   The ECP centre is centre C.                                        *
*   The integrals come out in the order cont_a, cont_b, cmpt_a, cmpt_b *
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
      integer i,k,l,
     &    i_ai,ibug,i_ca,i_Ga,i_bi,i_cb,i_Gb,i_free,i_cont_c,i_co,
     &    i_a,i_b,i_c,i_d,i_k,i_t,i_x,i_y,i_z,i_int_a,i_int_b,i_coef,
     &    i_off,i_Q_int,i_Qh,i_Q,i_Qabc,i_gam,i_bet,i_alp,i_pre,i_tmp,
     &    i_Qa,i_Qb,i_Qab,i_ang_a,i_ang_b,i_ang,i_wa,i_wb,k_a,k_b,
     &    l_a,l_sa,l_b,l_sb,l_c,l_max,l_min,l_cut,l_ecp,l_ecp_max,lam_a,
     &    lam_b,la_min,lb_min,m_count,mt_max,mt_min,mem_max,
     &    n_na,n_all_a,n_nsa,n_ta,n_nb,n_all_b,n_nsb,n_tb,
     &    n_ab,n_abc,nc_ab,ncab,n_rad,n_Q,n_row,n_ang_a,n_ang_b,
     &    n_l,n_l_c,n_lam_a,n_lam_b,ltmp,n_ints,n_c_tot,n_z_tot,n_blk,
     &    n_prim_a,n_prim_b,n_cont_a,n_cont_b,n_cont_c,n_int_a,n_int_b,
     &    n_cent,n_co,n_co_tot,n_co_max,n_pass,n_x,nu_a,nu_b,
     &    p_max,p_min,sphcart
      integer n_coef_c(-1:l_ecp_max,n_cent,2),
     &    n_prim_c(p_min:p_max,-1:l_ecp_max,n_cent,2),
     &    ind_Z(-1:l_ecp_max,n_cent,2),ind_C(-1:l_ecp_max,n_cent,2)
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays,skip
      double precision zeta_c(n_z_tot),coef_c(n_c_tot),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(ltmp),ecp_ints(n_ints,n_blk),
     &    R_AC,X_AC,Y_AC,Z_AC,R_BC,X_BC,Y_BC,Z_BC,
     &    tol,fac,log_prefactor
      double precision csco(lcsco)
*
      logical ecp_skipint
      external ecp_skipint
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
*
      if (debug_gen) then
        write (LuOut,'(//A,/)') 'Entering ecp_int2 ...'
        write (LuOut,*) 'l_a,l_b',l_a,l_b
        write(LuOut,*) 'p_min,p_max', p_min,p_max,l_ecp_max
      end if
      if (DryRun) mem_max = 0
*
*   Allocate memory for ecp-independent quantities
*
      n_na = (l_a+1)*(l_a+2)/2
      n_all_a = n_na*(l_a+3)/3
      n_nb = (l_b+1)*(l_b+2)/2
      n_all_b = n_nb*(l_b+3)/3
      l_sa = l_ecp+l_a
      n_nsa = (l_sa+1)**2
      n_ta = (l_sa+1)*(l_sa+2)/2
      l_sb = l_ecp+l_b
      n_nsb = (l_sb+1)**2
      n_tb = (l_sb+1)*(l_sb+2)/2
*
      n_ab = n_prim_a*n_prim_b
      nc_ab = n_prim_a*n_cont_b
      ncab = n_cont_a*n_cont_b
      if (debug_addresses) then
        write (LuOut,*) 'n_prim_a,n_prim_b,n_ab,nc_ab,ncab',
     &      n_prim_a,n_prim_b,n_ab,nc_ab,ncab
      end if
*
      i_ca = 1
      i_cb = i_ca+n_na*n_all_a
      i_Ga = i_cb+n_nb*n_all_b
      i_Gb = i_Ga+n_nsa
      i_free = i_Gb+n_nsb
      i = i_free+max(3*l_sa+n_ta,3*l_sb+n_tb,3*l_a+n_na**2,
     &    3*l_b+n_nb**2)+3
      if (debug_addresses) then
        write (LuOut,*) 'n_na,n_all_a,l_sa,n_nsa,n_ta',
     &      n_na,n_all_a,l_sa,n_nsa,n_ta
        write (LuOut,*) 'n_nb,n_all_b,l_sb,n_nsb,n_tb',
     &      n_nb,n_all_b,l_sb,n_nsb,n_tb
        write (LuOut,*) 'i_ca,i_cb,i_Ga,i_Gb',i_ca,i_cb,i_Ga,i_Gb
      end if
      if (DryRun) then
        mem_max = max(mem_max,i-1)
        if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
      else
        if (i-1 .gt. ltmp) call errquit(
     &      ' Insufficient memory in ecp_int2',99, MEM_ERR)
*
*   Expand cartesian basis about ECP centre in spherical harmonics
*
        i_x = i_Ga
        i_y = i_x+l_a+1
        i_z = i_y+l_a+1
        i_t = i_z+l_a+1
        call ecp_cart_xpd (l_a,n_na,n_all_a,X_AC,Y_AC,Z_AC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_ca),1,
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr (tmp(i_ca),1,n_na,1,n_all_a,1,
     &      n_na,1,n_all_a,'Cartesian expansion for centre a','E',78,4)
        i_y = i_x+l_b+1
        i_z = i_y+l_b+1
        i_t = i_z+l_b+1
        call ecp_cart_xpd (l_b,n_nb,n_all_b,X_BC,Y_BC,Z_BC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_cb),1,
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr (tmp(i_cb),1,n_nb,1,n_all_b,1,
     &      n_nb,1,n_all_b,'Cartesian expansion for centre b','E',78,4)
*
*  Set up spherical tensors which multiply bessel functions
*
        i_x = i_free
        i_y = i_x+l_sa+1
        i_z = i_y+l_sa+1
        i_t = i_z+l_sa+1
        call ecp_sph_tens (l_sa,n_nsa,n_ta,R_AC,X_AC,Y_AC,Z_AC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_Ga),
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr (tmp(i_Ga),1,n_nsa,1,1,
     &      1,n_nsa,1,1,'Spherical tensors for centre a','E',78,4)
        i_y = i_x+l_sb+1
        i_z = i_y+l_sb+1
        i_t = i_z+l_sb+1
        call ecp_sph_tens (l_sb,n_nsb,n_tb,R_BC,X_BC,Y_BC,Z_BC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_Gb),
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr (tmp(i_Gb),1,n_nsb,1,1,
     &      1,n_nsb,1,1,'Spherical tensors for centre b','E',78,4)
      end if
*
*     Loop over ECP angular momentum projectors
*
      do l_c = 0,l_ecp
*
*     Define parameters for the different integral class cases
*
        call ecp_ind_var (l_c,n_blk,n_coef_c(l_c,1,1),n_coef_c(l_c,1,2),
     &      ind_c(l_c,1,1),ind_c(l_c,1,2),ind_z(l_c,1,1),ind_z(l_c,1,2),
     &      n_x,n_co_tot,n_co_max,i_off,n_pass,i_cont_c,n_cont_c,skip)
        if (skip) go to 10
*
*     Check magnitude of integrals
*
        if (.not.DryRun) then
          if (n_blk .eq. 1) then
            if (ecp_skipint (
     &          l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &          l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &          n_coef_c(l_c,1,1),zeta_c(ind_z(l_c,1,1)),
     &          coef_c(ind_c(l_c,1,1)))) go to 10
          else if (n_blk .eq. 3) then
            if (ecp_skipint (
     &          l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &          l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &          n_coef_c(l_c,1,2),zeta_c(ind_z(l_c,1,2)),
     &          coef_c(ind_c(l_c,1,2)))) go to 10
          else
            if (ecp_skipint (
     &          l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &          l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &          n_coef_c(l_c,1,1),zeta_c(ind_z(l_c,1,1)),
     &          coef_c(ind_c(l_c,1,1))) .and.
     &          ecp_skipint (
     &          l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &          l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &          n_coef_c(l_c,1,2),zeta_c(ind_z(l_c,1,2)),
     &          coef_c(ind_c(l_c,1,2)))) go to 10
          end if
        end if
*
*   Set up array dimensions for integrals
*
        m_count = l_a+l_b+1+abs(l_a-l_b)/2
        l_max = max(l_a,l_b)
        l_min = min(l_a,l_b)
        mt_max = l_max*(l_max+1)/2
        mt_min = l_min*(l_min-1)/2
        n_Q = 2*mt_max+mt_min
        if (l_a .ne. l_b) n_Q = n_Q-l_max
        n_Q = max(n_Q,l_max+1)
*
        l_min = min(l_c,l_a)
        l_cut = max(l_a-l_c,0)
        n_int_a = (l_min+1)*(2*l_cut+l_min+2)/2
        if (sphcart .eq. 0) then
          l_min = l_cut/2
          n_int_a = n_int_a+l_min*(l_cut-l_min)
        end if
*
        l_min = min(l_c,l_b)
        l_cut = max(l_b-l_c,0)
        n_int_b = (l_min+1)*(2*l_cut+l_min+2)/2
        if (sphcart .eq. 0) then
          l_min = l_cut/2
          n_int_b = n_int_b+l_min*(l_cut-l_min)
        end if
        n_rad = ncab*n_int_a*n_int_b
        n_abc = n_ab*n_co_max
        if (debug_addresses) then
          write (LuOut,*) 'n_rad,n_abc,n_int_a,n_int_b,n_Q',
     &        n_rad,n_abc,n_int_a,n_int_b,n_Q
        end if
*
*   Allocate scratch memory for integrals
*
        i_Q_int = i_free
        i_Qh = i_Q_int+n_rad*n_cont_c
        i_Qabc = i_free
        i_tmp = i_Qabc+3*m_count*n_abc
        i_Qab = max(i_tmp+20*n_abc,i_Qh)
        i_Q = i_Qab+n_ab*(l_a+1)*(l_b+1)*n_cont_c
        i_ai = i_Q+n_ab*n_Q*n_cont_c
        i_bi = i_ai+n_abc
        i_gam = i_bi+n_abc
        i_bet = i_gam+n_abc
        i_alp = i_bet+n_abc
        i_pre = i_alp+n_abc
        i_coef = i_pre+n_abc
        i_free = i_coef+n_co_tot
        if (debug_addresses) then
          write (LuOut,*) 'i_Q_int,i_Qh,i_Qabc,i_tmp,i_Qab,i_Q,i_ai,',
     &        'i_bi,i_gam,i_bet,i_alp,i_pre,i_free'
          write (LuOut,*) i_Q_int,i_Qh,i_Qabc,i_tmp,i_Qab,i_Q,i_ai,
     &        i_bi,i_gam,i_bet,i_alp,i_pre,i_free
        end if
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
        else
          if (i_free-1 .gt. ltmp) call errquit(
     &        ' Insufficient memory in ecp_int2',99, MEM_ERR)
          i_t = i_cont_c
          do i_y = 1,n_pass
            i_z = ind_z(l_c,1,i_t)
            i_co = ind_c(l_c,1,i_t)
            n_co = n_coef_c(l_c,1,i_t)
#if 1
      call ecp_int2_k(
     R     r_bc,r_ac,tmp(i_pre),
     G     tmp(i_gam),zeta_c,zeta_b,zeta_a,
     B     tmp(i_bet),tmp(i_alp),tmp(i_bi),tmp(i_ai),
     I     n_co,n_prim_b,n_prim_a,i_z)
#else
*
*   Set up argument values etc
*
            i = 0
            do i_c = 0,n_co-1
              do i_b = 1,n_prim_b
                do i_a = 1,n_prim_a
                  tmp(i_gam+i) = one/sqrt(zeta_c(i_z+i_c)+zeta_b(i_b)
     &                +zeta_a(i_a))
                  tmp(i_bet+i) = R_bc*zeta_b(i_b)*tmp(i_gam+i)
                  tmp(i_alp+i) = R_ac*zeta_a(i_a)*tmp(i_gam+i)
                  tmp(i_bi+i) = one/(two*R_bc*zeta_b(i_b))
                  tmp(i_ai+i) = one/(two*R_ac*zeta_a(i_a))
                  log_prefactor = (tmp(i_alp+i)+tmp(i_bet+i))**2
     &                -zeta_a(i_a)*R_ac**2
     &                -zeta_b(i_b)*R_bc**2
                  tmp(i_pre+i) = exp(log_prefactor)
#ifdef DEBUG
                  if (debug_arrays) then
                    write (LuOut,'(3i5,1p3e20.10)') i_a,i_b,i_c,
     &                  zeta_a(i_a),zeta_b(i_b),zeta_c(i_z+i_c)
                    write (LuOut,'(10x,i5,1p3e20.10)') i,tmp(i_alp+i),
     &                  tmp(i_bet+i),tmp(i_gam+i)
                    write (LuOut,'(15x,1p2e20.10)') tmp(i_ai+i),
     &                  tmp(i_bi+i)
                    write (LuOut,'(15x,1p2e20.10)') log_prefactor,
     &                  tmp(i_pre+i)
                  end if
#endif
                  i = i+1
                end do
              end do
            end do
#endif
*
            n_abc = n_ab*n_co
            call dcopy (n_co,coef_c(i_co),1,tmp(i_coef),1)
            if (n_cont_c .eq. 2) call dcopy 
     &          (n_co,coef_c(ind_c(l_c,1,2)),1,tmp(i_coef+n_co),1)
#ifdef DEBUG
            if (debug_addresses) write (LuOut,*) 
     &          'n_ab,n_abc,n_Q,m_count',
     &          n_ab,n_abc,n_Q,m_count
#endif 
*
*   Calculate radial integrals
*
            call ecp_radint2 (p_min,p_max,
     &          l_a,n_prim_a,n_cont_a,n_int_a,coef_a,tmp(i_ai),
     &          tmp(i_alp),
     &          l_b,n_prim_b,n_cont_b,n_int_b,coef_b,tmp(i_bi),
     &          tmp(i_bet),
     &          l_c,n_prim_c(p_min,l_c,1,i_t),n_cont_c,n_co,
     &          tmp(i_coef),tmp(i_gam),tmp(i_pre),tol,sphcart,
     &          n_ab,n_abc,n_Q,m_count,
     &          tmp(i_tmp),tmp(i_Qabc),tmp(i_Qab),tmp(i_Q),tmp(i_Q_int),
     &          ibug/10)
            i_t = i_t+1
          end do
        end if
*
*     Allocate memory for angular integrals
*
        n_l_c = 2*l_c+1
        l_min = min(l_a,l_c)
        n_ang_a = (l_min+1)*(l_a+1)**2-l_min*(l_min+1)*(2*l_min+1)/6
        l_min = min(l_b,l_c)
        n_ang_b = (l_min+1)*(l_b+1)**2-l_min*(l_min+1)*(2*l_min+1)/6
*
        i_ang = i_Qh
        i_ang_a = i_ang+n_ang_a*n_ang_b*n_x
        if (n_x .gt. 1) i_ang_a = max(i_ang_a,i_ang+n_ang_a*n_l_c)
        i_ang_b = i_ang_a+n_l_c*n_ang_a*n_x
        i_free = i_ang_b+n_l_c*n_ang_b
        if (debug_addresses) then
          write (LuOut,*) 'n_l_c,n_ang_a,n_ang_b',n_l_c,n_ang_a,n_ang_b
          write (LuOut,*) 'i_ang,i_ang_a,i_ang_b',i_ang,i_ang_a,i_ang_b
        end if
*
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
          if (debug_gen) write (LuOut,*) ' This is a dry run'
        else
          if (i_free-1 .gt. ltmp) call errquit(
     &        ' Insufficient memory in ecp_int2',99, MEM_ERR)
*
*       Set up angular coefficients and contract over components of
*       spherical tensors (sum over q and q').
*
          i = 0
          if (debug_addresses) write (LuOut,*) 'l_max',l_max
          do l = 0,l_max
            n_l = l+l+1
            do k = l+l_c,abs(l-l_c),-2
              i_k = k**2
              if (l .le. l_a)
     &            call ecp_angint (tmp(i_ang_a+i),l,k,l_c,tmp(i_Ga+i_k))
              if (l .le. l_b)
     &            call ecp_angint (tmp(i_ang_b+i),l,k,l_c,tmp(i_Gb+i_k))
              i = i+n_l_c*n_l
            end do
          end do
          if (debug_arrays) then
            call ecp_matpr (tmp(i_ang_a),-l_c,l_c,1,n_ang_a,
     &          -l_c,l_c,1,n_ang_a,'A Angular integrals','F',78,4)
            call ecp_matpr (tmp(i_ang_b),-l_c,l_c,1,n_ang_b,
     &          -l_c,l_c,1,n_ang_a,'B Angular integrals','F',78,4)
          end if
*
*       Perform sums over m,m'. For the scalar integrals, this requires
*       the insertion of a factor of (1+delta_m,0); for the spin-orbit
*       integrals a routine is called to sum over m.
*
          if (n_x .gt. 1) then
            call dcopy (n_ang_a*n_l_c,tmp(i_ang_a),1,tmp(i_ang),1)
            i = i_ang_a+(n_x-3)*n_l_c*n_ang_a
            call ecp_so_ang (l_c,tmp(i_ang),tmp(i),n_ang_a)
          end if
          if (n_x .ne. 3) call dscal (n_ang_a,two,tmp(i_ang_a+l_c),
     &          n_l_c)
*
*       Sum over m for scalar integrals
*
          fac = n_l_c*2
          fac = fac*pi
          if (n_x .ne. 3) then
            call dgemm ('T','N',n_ang_a,n_ang_b,n_l_c,fac,tmp(i_ang_a),
     &          n_l_c,tmp(i_ang_b),n_l_c,zero,tmp(i_ang),n_ang_a)
#ifdef DEBUG
            if (debug_arrays) then
              call ecp_matpr (tmp(i_ang_a),-l_c,l_c,1,n_ang_a,-l_c,l_c,
     &            1,n_ang_a,'Modified scalar A Angular integrals','F',
     &            78,4)
              call ecp_matpr (tmp(i_ang),1,n_ang_a,1,n_ang_b,1,n_ang_a,
     &            1,n_ang_b,'Combined scalar angular integrals','F',78,
     &            4)
            end if
#endif
          end if
*
*       Sum over m for spin-orbit integrals [factor 4\pi]
*
          if (n_x .ge. 3) then
            i_a = i_ang_a+(n_x-3)*n_ang_a*n_l_c
            i_b = i_ang+(n_x-3)*n_ang_a*n_ang_b
            do i = 1,3
              call dgemm ('T','N',n_ang_a,n_ang_b,n_l_c,fac,tmp(i_a),
     &            n_l_c,tmp(i_ang_b),n_l_c,zero,tmp(i_b),n_ang_a)
              if (debug_arrays) then
                call ecp_matpr (tmp(i_a),-l_c,l_c,1,n_ang_a,-l_c,l_c,1,
     &              n_ang_a,'Modified A SO Angular integrals','F',78,4)
                call ecp_matpr (tmp(i_b),1,n_ang_a,1,n_ang_b,1,n_ang_a,
     &              1,n_ang_b,'Combined SO angular integrals','F',78,4)
              end if
              i_a = i_a+n_ang_a*n_l_c
              i_b = i_b+n_ang_a*n_ang_b
            end do
          end if
        end if
*
*     Now loop over angular momenta of expanded functions and
*     perform contraction of angular intgrals with radial integrals 
*     and expansion coefficients
*
        i_wa = i_ang_a
        i_wb = i_wa+n_nb*n_ang_a
        i_free = i_wb+n_na*n_nb
        if (debug_addresses) write (LuOut,*) 'i_wa,i_wb',i_wa,i_wb
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
        else
          if (i_free-1 .gt. ltmp) call errquit(
     &        ' Insufficient memory in ecp_int2',99, MEM_ERR)
          n_row = ncab*n_int_a
          do i_x = 1,n_x
            i_b = i_cb
            i_d = i_Q_int
            do nu_b = 0,l_b
              i_Qb = i_d
              if (sphcart .eq. 0) then
                lb_min = mod(nu_b,2)
              else
                lb_min = nu_b
              end if
              do lam_b = nu_b,lb_min,-2
                l_min = min(lam_b-1,l_c)
                i_ang_b = i_ang+n_ang_a*((l_min+1)*lam_b**2
     &              -l_min*(l_min+1)*(2*l_min+1)/6)
                n_lam_b = 2*lam_b+1
                i_int_b = i_Qb
                do k_b = lam_b+l_c,abs(lam_b-l_c),-2
                  fac = 2*k_b+1
                  call dgemm ('N','T',n_ang_a,n_nb,n_lam_b,fac,
     &                tmp(i_ang_b),n_ang_a,tmp(i_b),n_nb,zero,
     &                tmp(i_wa),n_ang_a)
#ifdef DEBUG
                  if (debug_gen) write (LuOut,*) 'nu_b,lam_b,k_b',
     &                nu_b,lam_b,k_b
                  if (debug_arrays) call ecp_matpr (tmp(i_wa),
     &                1,n_ang_a,1,n_nb,1,n_ang_a,1,n_nb,
     &                'Angular integrals contracted with b expansion',
     &                'E',78,4)
#endif
                  i_ang_b = i_ang_b+n_lam_b*n_ang_a
*
                  i_a = i_ca
                  i_c = i_int_b
                  do nu_a = 0,l_a
                    i_Qa = i_c
                    if (sphcart .eq. 0) then
                      la_min = mod(nu_a,2)
                    else
                      la_min = nu_a
                    end if
                    do lam_a = nu_a,la_min,-2
                      l_min = min(lam_a-1,l_c)
                      i_ang_a = i_wa+((l_min+1)*lam_a**2
     &                    -l_min*(l_min+1)*(2*l_min+1)/6)
                      n_lam_a = 2*lam_a+1
                      i_int_a = i_Qa
                      do k_a = lam_a+l_c,abs(lam_a-l_c),-2
                        fac = 2*k_a+1
#ifdef DEBUG
                        if (debug_gen) write (LuOut,*) 'nu_a,lam_a,k_a',
     &                      nu_a,lam_a,k_a
#endif
                        call dgemm ('N','N',n_na,n_nb,n_lam_a,fac,
     &                      tmp(i_a),n_na,tmp(i_ang_a),n_ang_a,zero,
     &                      tmp(i_wb),n_na)
                        i_ang_a = i_ang_a+n_lam_a
                        call ecp_angrad (n_na,n_cont_a,n_nb,n_cont_b,
     &                      tmp(i_wb),tmp(i_int_a),ecp_ints(1,i_x))
                        if (debug_arrays) then
                          call ecp_matpr (tmp(i_wb),1,n_na,1,n_nb,
     &                        1,n_na,1,n_nb,
     &                        'Contracted angular integrals','E',78,4)
                          call ecp_matpr (tmp(i_int_a),1,n_cont_a,
     &                        1,n_cont_b,1,n_cont_a,1,n_cont_b,
     &                        'Radial integrals','E',78,4)
                          call ecp_matpr (ecp_ints(1,i_x),1,
     &                        n_cont_a*n_na,1,n_nb*n_cont_b,1,
     &                        n_cont_a*n_na,1,n_nb*n_cont_b,
     &                        'ECP integrals','E',78,4)
                        end if
                        i_int_a = i_int_a+ncab
                      end do
                      i_a = i_a+n_lam_a*n_na
                      i_Qa = i_Qa+ncab
                    end do
                    if (sphcart .eq. 0) then
                      i_c = i_c+(min(nu_a,(nu_a+l_c)/2)+1)*ncab
                    else
                      i_c = i_c+(min(nu_a,l_c)+1)*ncab
                    end if
                  end do
*
                  i_int_b = i_int_b+n_row
                end do
                i_b = i_b+n_lam_b*n_nb
                i_Qb = i_Qb+n_row
              end do
              if (sphcart .eq. 0) then
                i_d = i_d+(min(nu_b,(nu_b+l_c)/2)+1)*n_row
              else
                i_d = i_d+(min(nu_b,l_c)+1)*n_row
              end if
            end do
            i_ang = i_ang+n_ang_a*n_ang_b
            if ((n_x .eq. 4) .and. (i_x .eq. 1))
     &          i_Q_int = i_Q_int+ncab*n_int_a*n_int_b
          end do ! i_x
*
        end if
   10   continue
*
      end do
#ifdef DEBUG      
      do i = 1,n_blk
        if (debug_arrays) call ecp_matpr (ecp_ints(1,i),1,n_nb*n_cont_b,
     &      1,n_na*n_cont_a,1,n_nb*n_cont_b,1,n_na*n_cont_a,
     &      'ECP integrals','E',78,4)
      end do
      if (debug_gen) write (LuOut,*) 'Exiting ecp_int2'
#endif
*
      return
      end
