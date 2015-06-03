C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_local2 (mem_max,DryRun,
     &    R_AC,X_AC,Y_AC,Z_AC,l_a,n_prim_a,n_cont_a,coef_a,zeta_a,
     &    R_BC,X_BC,Y_BC,Z_BC,l_b,n_prim_b,n_cont_b,coef_b,zeta_b,
     &    n_prim_c,n_coef_c,zeta_c,coef_c,p_min,p_max,
     &    tol,sphcart,tmp,l_tmp,
     &    csco,lcsco,
     &    ecp_ints,n_int,ibug)
*                                                                      *
*   Calculate local ecp integrals for a given ECP centre for two       *
*   off-centre basis functions and add to integral array.              *
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
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_coef_c (inp) - number of coefficients/exponents for local potl.  *
*   zeta_c (inp) - ECP exponents                                       *
*   coef_c (inp) - ECP contraction coefficients                        *
*   p_min (inp) - minimum power of r in ECP expansion                  *
*   p_max (inp) - maximum power of r in ECP expansion                  *
*   tol (inp) - maximum relative error in bessel functions             *
*   sphcart (inp) - 0 for cartesian integrals, 1 for spherical         *
*   tmp (scr) - work array                                             *
*   l_tmp (inp) - length of tmp                                        *
*   ecp_ints (out) - integrals over ECP                                *
*   ibug (inp) - debug flag. 0 for no debug, 1 for address printing,   *
*           2 for array printing, 3 for both.                          *
*                                                                      *
*   Notes:                                                             *
*   -----                                                              *
*                                                                      *
*   The ECP centre is centre C. Centre B is assumed to coincide with   *
*   centre C.                                                          *
*   The integrals come out in the order cont_a, cont_b, cmpt_a, cmpt_b *
*      where cont = contracted functions, cmpt = cartesian components  *
*   The integrals are added to the array ecp_ints, i.e. ecp_ints is    *
*   incremented by the integrals from this routine.                    *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
#include "ecp_consts.fh"
#include "util.fh"
#include "errquit.fh"
      integer l_a,n_prim_a,n_cont_a,l_b,n_prim_b,n_cont_b,n_coef_c,
     &    l_tmp,n_int,p_min,p_max,sphcart,mem_max,ibug
      integer n_prim_c(p_min:p_max)
      integer i,ii,i_a,i_b,i_c,i_ab,i_ai,i_ca,i_cb,i_Ga,i_free,i_ind,
     &    i_inta,i_intp,i_m,i_n,i_p,i_q,i_r,i_s,i_t,i_x,i_y,i_z,
     &    i_Q_int,i_Qa,i_Qabc,i_Qp,i_gam,i_alp,i_pre,i_tmp,
     &    i_rpc,i_xpc,i_ypc,i_zpc,i_prod,
     &    j,k,l,l_c,l_c1,l_p,ll,ll1,ll_a,ll_b,m,n,n_n,n_all,
     &    n_na,n_all_a,n_bas_a,n_pbas_a,n_nb,n_all_b,n_bas_b,n_pbas_b,
     &    n_nc,n_ncs,n_all_c,n_ab,n_abc,n_abp,n_rad,
     &    xa,xb,ya,yb,xx,yy
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays
      double precision zeta_c(n_coef_c),coef_c(n_coef_c),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(l_tmp),ecp_ints(n_int),
     &    R_AC,X_AC,Y_AC,Z_AC,R_BC,X_BC,Y_BC,Z_BC,
     &    X_AP,Y_AP,Z_AP,X_BP,Y_BP,Z_BP,Rsq,
     &    tol,log_prefac,ab,fourpi
      double precision csco(lcsco)
*
      logical ecp_skipint
      external ecp_skipint
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
*
      if (DryRun) mem_max = 0
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_local2 ...'
      if (debug_gen) write (LuOut,*) 'ibug =',ibug
      if (debug_addresses) then
        write (LuOut,*) 'Scratch memory l_tmp = ',l_tmp
        write (LuOut,*) 'p_min,p_max',p_min,p_max
        write (LuOut,*) 'n_prim_c',(n_prim_c(i),i = p_min,p_max)
        write (LuOut,*) 'l_a,l_b',l_a,l_b
      end if
*
*     Check magnitude of integrals
*
      if (.not.DryRun) then
        if (ecp_skipint (
     &      l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &      l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &      n_coef_c,zeta_c,coef_c)) return
      end if
*
*   Allocate memory for ecp-independent quantities
*
      n_na = (l_a+1)*(l_a+2)/2
      n_all_a = n_na*(l_a+3)/3
      n_nb = (l_b+1)*(l_b+2)/2
      n_all_b = n_nb*(l_b+3)/3
      if (debug_addresses) 
     &    write (LuOut,*) 'n_na,n_all_a,n_nb,n_all_b',
     &    n_na,n_all_a,n_nb,n_all_b
*
      n_ab = n_prim_a*n_prim_b
      n_bas_a = n_na*n_cont_a
      n_bas_b = n_nb*n_cont_b
      n_pbas_a = n_na*n_prim_a
      n_pbas_b = n_nb*n_prim_b
      n_abp = n_bas_a*n_nb
      n_abc = n_ab*n_coef_c
      if (debug_addresses) then
        write (LuOut,*) 'n_ab,n_abp,n_abc',n_ab,n_abp,n_abc
        write (LuOut,*) 'n_bas_a,n_bas_b,n_pbas_a,n_pbas_b',
     &      n_bas_a,n_bas_b,n_pbas_a,n_pbas_b
      end if
*
      l_c = l_a+l_b
      l_c1 = l_c+1
      n_nc = (l_c1)*(l_c1+1)/2
      n_ncs = l_c1**2
      n_all_c = n_nc*(l_c1+2)/3
*
*
*   Set up array dimensions for integrals and pointers
*
      l = l_c/2
      n_rad = (l+1)*(l_c1-l)
      i_xpc = 1
      i_ypc = i_xpc+n_ab
      i_zpc = i_ypc+n_ab
      i_rpc = i_zpc+n_ab
      i_ab = i_rpc+n_ab
      i_prod = i_ab+n_ab
      j = i_prod+n_ab
*
      i_inta = i_ab
      i_Q_int = i_inta+n_bas_a*n_pbas_b
      i_Qa = i_Q_int+n_ab*n_rad
      i_Qabc = i_Qa+n_ab*l_c1
      i_ai = i_Qabc+n_abc
      i_gam = i_ai+n_abc
      i_alp = i_gam+n_abc
      i_pre = i_alp+n_abc
      i_tmp = i_pre+n_abc
      i_ind = i_tmp+n_abc*6
      i_free = i_ind+n_abc
*
      i_intp = i_Qa
      i_Qp = i_intp+n_pbas_a*n_nb
      i_Ga = i_Qp+n_all_c
      i_ca = i_Ga+n_ncs
      i_cb = i_ca+n_na*n_all_a
      i_x = max(i_ca+n_nc*n_all_c,i_cb+n_nb*n_all_b)
      i_y = i_x+l_c1
      i_z = i_y+l_c1
      i_t = i_z+l_c1
      i = i_t+n_nc**2
      if (debug_addresses) then
        write (LuOut,*) 'i_xpc,i_ypc,i_zpc,i_rpc',
     &      i_xpc,i_ypc,i_zpc,i_rpc
        write (LuOut,*) 'i_ab,i_prod,i_inta,i_intp',
     &      i_ab,i_prod,i_inta,i_intp
        write (LuOut,*) 'i_Q_int,i_Qa,i_Qabc',
     &      i_Q_int,i_Qa,i_Qabc
        write (LuOut,*) 'i_ai,i_gam,i_alp,i_pre,i_tmp,i_ind',
     &      i_ai,i_gam,i_alp,i_pre,i_tmp,i_ind
        write (LuOut,*) 'i_Qp,i_Ga,i_ca,i_cb',
     &      i_Qp,i_Ga,i_ca,i_cb
        write (LuOut,*) 'i_x,i_y,i_z,i_t',
     &      i_x,i_y,i_z,i_t
        write (LuOut,*) 'i_free,i',i_free,i
      end if
      if (DryRun) then
        mem_max = max(mem_max,i_free-1,i-1,j-1)
        if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
      else
        if (max(i_free,i,j)-1 .gt. l_tmp) call errquit(
     &      ' Insufficient memory in ecp_local2',99, MEM_ERR)
        fourpi = pi+pi
        fourpi = fourpi+fourpi
#if 1
      call ecp_local2_k(
     A       tmp(i_ab),tmp(i_xpc),tmp(i_ypc),tmp(i_zpc),
     Z       zeta_a,zeta_b,zeta_c,
     b       tmp(i_prod),tmp(i_gam),tmp(i_alp),tmp(i_ai),
     R       tmp(i_rpc),tmp(i_pre),
     X       x_ac,x_bc,y_ac,y_bc,z_ac,z_bc,
     I       n_prim_a,n_prim_b,n_coef_c,n_ab
     S       )

#else
*
*     Loop over primitives on A and B and do Gaussian product
*
        Rsq = (X_AC-X_BC)**2+(Y_AC-Y_BC)**2+(Z_AC-Z_BC)**2
        i = 0
        do i_b = 1,n_prim_b
          do i_a = 1,n_prim_a
            ab = zeta_b(i_b)+zeta_a(i_a)
            tmp(i_ab+i) = ab
            tmp(i_xpc+i) = (X_AC*zeta_a(i_a)+X_BC*zeta_b(i_b))/ab
            tmp(i_ypc+i) = (Y_AC*zeta_a(i_a)+Y_BC*zeta_b(i_b))/ab
            tmp(i_zpc+i) = (Z_AC*zeta_a(i_a)+Z_BC*zeta_b(i_b))/ab
            tmp(i_rpc+i) = sqrt(tmp(i_xpc+i)**2+tmp(i_ypc+i)**2
     &          +tmp(i_zpc+i)**2)
            tmp(i_prod+i) = Rsq*zeta_b(i_b)*zeta_a(i_a)/ab
            i = i+1
          end do
        end do
*
*     Set up argument values for radial integrals
*
        i = 0
        do i_c = 1,n_coef_c
          do j = 0,n_ab-1
            tmp(i_gam+i) = one/sqrt(zeta_c(i_c)+tmp(i_ab+j))
            tmp(i_alp+i) = tmp(i_rpc+j)*tmp(i_ab+j)*tmp(i_gam+i)
            if (tmp(i_rpc+j) .eq. zero) then
              tmp(i_ai+i) = zero
            else
              tmp(i_ai+i) = one/(two*tmp(i_rpc+j)*tmp(i_ab+j))
            end if
            log_prefac = tmp(i_alp+i)**2-tmp(i_ab+j)*tmp(i_rpc+j)**2
     &          -tmp(i_prod+j)
            tmp(i_pre+i) = zero
            if (log_prefac.gt.ln_thr_ecp) 
     &         tmp(i_pre+i) = exp(log_prefac)
            i = i+1
          end do
        end do
#endif
        if (debug_arrays) then
          call ecp_matpr(tmp(i_gam),1,n_ab,1,n_coef_c,
     &          1,n_ab,1,n_coef_c,'gam','E',78,4)
          call ecp_matpr(tmp(i_alp),1,n_ab,1,n_coef_c,
     &          1,n_ab,1,n_coef_c,'alp','E',78,4)
          call ecp_matpr(tmp(i_ai),1,n_ab,1,n_coef_c,
     &          1,n_ab,1,n_coef_c,'ai','E',78,4)
          call ecp_matpr(tmp(i_pre),1,n_ab,1,n_coef_c,
     &          1,n_ab,1,n_coef_c,'pre','E',78,4)
        end if
*
*      Evaluate radial integrals
*
!        call dfill(n_ab*n_rad,zero,tmp(i_Q_int),1)
        call dcopy(n_ab*n_rad,zero,0,tmp(i_Q_int),1)
        call ecp_locint1 (p_min,p_max,l_c,n_prim_a,n_prim_b,n_coef_c,
     &      n_prim_c,coef_c,tmp(i_ai),tmp(i_gam),tmp(i_alp),tmp(i_pre),
     &      tol,n_ab,n_abc,n_rad,tmp(i_tmp),tmp(i_ind),tmp(i_Qabc),
     &      tmp(i_Qa),tmp(i_Q_int),ibug/10)
*
*     Loop over products of exponents and sum over radial integrals
*     and expansion coefficients to get integrals at the Gaussian
*     product centre.
*
        i = 0
        i_n = i_inta
        do i_b = 1,n_prim_b
!          call dfill(n_nb*n_pbas_a,zero,tmp(i_intp),1)
          call dcopy(n_nb*n_pbas_a,zero,0,tmp(i_intp),1)
          i_m = i_intp
          do i_a = 1,n_prim_a
*
*       Set up spherical tensors which multiply bessel functions
*
            call ecp_sph_tens (l_c,n_ncs,n_nc,tmp(i_rpc+i),tmp(i_xpc+i),
     &          tmp(i_ypc+i),tmp(i_zpc+i),tmp(i_x),tmp(i_y),tmp(i_z),
     &          tmp(i_t),tmp(i_Ga),csco,lcsco)
            if (debug_arrays) call ecp_matpr(tmp(i_Ga),1,n_ncs,1,1,
     &          1,n_ncs,1,1,'Spherical tensors from GP','E',78,4)
*
*       Loop over angular momenta generated by expansion of Gaussian
*       product centre about ECP centre
*
            i_p = i_Qp
            do l_p = 0,l_c
*
*         Expand cartesian basis about ECP centre in spherical tensors
*
              n_n = (l_p+1)*(l_p+2)/2
              n_all = n_n*(l_p+3)/3
              call ecp_cart_xpd (l_p,n_n,n_all,tmp(i_xpc+i),
     &            tmp(i_ypc+i),tmp(i_zpc+i),tmp(i_x),tmp(i_y),tmp(i_z),
     &            tmp(i_t),tmp(i_ca),1,csco,lcsco)
              if (debug_arrays) call ecp_matpr(tmp(i_ca),1,n_n,1,n_all,
     &            1,n_n,1,n_all,'Cartesian expansion of GP','E',78,4)
*
*         Loop over sphericals in expansion; contract expansion of
*         cartesian with spherical tensor and contract result with radial
*         integrals. Result is cartesian integrals on GP centre.
*
              i_c = i_ca
              i_q = i_Q_int+i
              call dcopy(n_n,zero,0,tmp(i_p),1)
              do n = 0,l_p
                do l = n,0,-2
                  ll1 = l+l+1
                  i_s = i_Ga+l**2
                  call dgemv ('N',n_n,ll1,tmp(i_q),tmp(i_c),n_n,
     &                tmp(i_s),1,one,tmp(i_p),1)
#ifdef DEBUG
                  if (debug_arrays) call ecp_matpr(tmp(i_p),1,1,1,n_n,
     &                1,1,1,n_n,'intermediate sum','E',78,4)
#endif
                  i_c = i_c+n_n*ll1
                  i_q = i_q+n_ab
                end do
              end do
#ifdef DEBUG
              if (debug_arrays) call ecp_matpr(tmp(i_p),1,1,1,n_n,
     &            1,1,1,n_n,'Cartesian integrals on GP centre','E',78,4)
#endif
              i_p = i_p+n_n
            end do
*
*       Generate cartesian expansion of basis functions a and b around GP centre
*
            X_AP = X_AC-tmp(i_xpc+i)
            Y_AP = Y_AC-tmp(i_ypc+i)
            Z_AP = Z_AC-tmp(i_zpc+i)
            call ecp_cart_xpd (l_a,n_na,n_all_a,X_AP,Y_AP,Z_AP,tmp(i_x),
     &          tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_ca),0,csco,lcsco)
            if (debug_arrays) call ecp_matpr(tmp(i_ca),
     &          1,n_na,1,n_all_a,1,n_na,1,n_all_a,
     &          'Cartesian expansion of a around GP cent','E',78,4)
            X_BP = X_BC-tmp(i_xpc+i)
            Y_BP = Y_BC-tmp(i_ypc+i)
            Z_BP = Z_BC-tmp(i_zpc+i)
            call ecp_cart_xpd (l_b,n_nb,n_all_b,X_BP,Y_BP,Z_BP,tmp(i_x),
     &          tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_cb),0,csco,lcsco)
            if (debug_arrays) call ecp_matpr(tmp(i_cb),
     &          1,n_n,1,n_all_b,1,n_nb,1,n_all_b,
     &          'Cartesian expansion of b around GP cent','E',78,4)
*
*       Loop over expanded functions on centres a and b, contract
*       expansion coefficients with integrals
*
            i_s = i_cb
            do ll_b = 0,l_b
              i_r = i_ca
              do ll_a = 0,l_a
                ll = ll_a+ll_b
                i_q = i_Qp+ll*(ll+1)*(ll+2)/6
                j = i_s
                do xb = ll_b,0,-1
                  do yb = ll_b-xb,0,-1
                    k = i_r
                    do xa = ll_a,0,-1
                      xx = ll-xa-xb
                      ii = i_q+xx*(xx+1)/2
                      do ya = ll_a-xa,0,-1
                        yy = xx-ya-yb
                        m = ii+yy
                        call dger (n_na,n_nb,tmp(m),tmp(k),1,tmp(j),1,
     &                      tmp(i_m),n_pbas_a)
                        k = k+n_na
                      end do
                    end do
                    j = j+n_nb
                  end do
                end do
#ifdef DEBUG
                if (debug_arrays) call ecp_matpr(tmp(i_m),
     &              1,n_pbas_a,1,n_nb,1,n_na,1,n_nb,
     &              'intermediate integrals',
     &              'E',78,4)
#endif
                i_r = i_r+n_na*(ll_a+1)*(ll_a+2)/2
              end do
              i_s = i_s+n_nb*(ll_b+1)*(ll_b+2)/2
            end do
            if (debug_arrays) call ecp_matpr(tmp(i_m),
     &          1,n_pbas_a,1,n_nb,1,n_na,1,n_nb,
     &          'primitive integrals',
     &          'E',78,4)
*
            i = i+1
            i_m = i_m+n_na
          end do
*
*     Contract over basis functions on a.
*     The factor of 4*pi is introduced at this stage
*
          i_m = i_intp
          do j = 1,n_nb
            if (debug_arrays) call ecp_matpr(tmp(i_m),
     &          1,n_na,1,n_prim_a,1,n_na,1,n_prim_a,
     &          'array of primitive integrals',
     &          'E',78,4)
            call dgemm ('N','N',n_na,n_cont_a,n_prim_a,fourpi,tmp(i_m),
     &          n_na,coef_a,n_prim_a,zero,tmp(i_n),n_na)
            i_m = i_m+n_pbas_a
            i_n = i_n+n_bas_a
          end do
        end do
*
*     Finally, contract over basis functions on b
*
        if (debug_arrays) call ecp_matpr(tmp(i_inta),
     &      1,n_bas_a,1,n_pbas_b,1,n_bas_a,1,n_pbas_b,
     &      'Half-contracted integrals',
     &      'E',78,4)
        call dgemm ('N','N',n_abp,n_cont_b,n_prim_b,one,tmp(i_inta),
     &      n_abp,coef_b,n_prim_b,one,ecp_ints,n_abp)
*
        if (debug_arrays) then
          call ecp_matpr (ecp_ints,1,n_bas_a,1,n_bas_b,
     &        1,n_bas_b,1,n_bas_a,'ECP integrals','E',78,4)
        end if
*
      end if
      if (debug_gen) write (LuOut,*) 'Exiting ecp_local2'
*
      return
      end
