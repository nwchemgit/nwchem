C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_local1 (mem_max,DryRun,
     &    R_AC,X_AC,Y_AC,Z_AC,l_a,n_prim_a,n_cont_a,coef_a,zeta_a,
     &    R_BC,X_BC,Y_BC,Z_BC,l_b,n_prim_b,n_cont_b,coef_b,zeta_b,
     &    n_prim_c,n_coef_c,zeta_c,coef_c,p_min,p_max,
     &    tol,sphcart,tmp,l_tmp,
     &    csco,lcsco,
     &    ecp_ints,n_int,transpose,ibug)
*                                                                      *
*   Calculate Type 1 local integrals for a given ECP centre            *
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
*   n_int (inp) - number of ECP integrals                              *
*   transpose (inp) - true if centres A and B are to be transposed.    *
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
      integer i,i_a,i_b,i_c,i_ai,i_ang,i_ang_a,i_ca,i_Ga,i_free,i_ind,
     &    i_int_a,i_k,i_t,i_x,i_y,i_z,i_Q_int,i_Qa,
     &    i_Qabc,i_gam,i_alp,i_pre,i_tmp,i_wrk,k,l,l_c,l_c_min,
     &    l_c_max,l_sa,l_max,l_min,ll,l_lo,m_count,m_skip,n,n_na,n_nsa,
     &    n_ta,n_all_a,n_nb,n_all_b,n_abc,n_ab,nc_ab,ncab,n_rad,n_l,
     &    n_l_c,n_ang_a,n_ang_t,n_Qa
      integer lcsco
      logical DryRun,transpose,debug_gen,debug_addresses,debug_arrays
      double precision zeta_c(n_coef_c),coef_c(n_coef_c),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(l_tmp),ecp_ints(n_int),
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
      if (DryRun) mem_max = 0
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_local1 ...'
      if (debug_gen) write (LuOut,*) 'ibug =',ibug
      if (debug_addresses) then
        write (LuOut,*) 'Scratch memory l_tmp = ',l_tmp
        write (LuOut,*) p_min,p_max
        write (LuOut,*) (n_prim_c(i),i = p_min,p_max)
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
      l_sa = l_b+l_a
      n_nsa = (l_sa+1)**2
      n_ta = (l_sa+1)*(l_sa+2)/2
      n_nb = (l_b+1)*(l_b+2)/2
      if (debug_addresses) 
     &    write (LuOut,*) 'n_na,n_all_a,l_sa,n_nsa,n_ta,n_nb',
     &    n_na,n_all_a,l_sa,n_nsa,n_ta,n_nb
*
      n_ab = n_prim_a*n_prim_b
      nc_ab = n_prim_a*n_cont_b
      ncab = n_cont_a*n_cont_b
      n_abc = n_ab*n_coef_c
      if (debug_addresses) 
     &    write (LuOut,*) 'n_ab,nc_ab,ncab',n_ab,nc_ab,ncab
*
      i_ca = 1
      i_Ga = i_ca+n_na*n_all_a
      i_x = i_Ga+n_nsa
      i_y = i_x+l_sa+1
      i_z = i_y+l_sa+1
      i_t = i_z+l_sa+1
      i_free = i_t+n_ta
      if (debug_addresses) then
        write (LuOut,*) 'i_ca,i_Ga,i_x,i_y,i_z,i_t,i_free',
     &      i_ca,i_Ga,i_x,i_y,i_z,i_t,i_free
      end if
      if (DryRun) then
        mem_max = max(mem_max,i_free-1)
        if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
      else
        if (i_free-1 .gt. l_tmp) call errquit(
     &      ' Insufficient memory in ecp_local1',99, MEM_ERR)
*
*     Expand cartesian basis about ECP centre in spherical tensors
*
        if (debug_arrays)
     &      write (LuOut,*) 'X_AC,Y_AC,Z_AC',X_AC,Y_AC,Z_AC
        call ecp_cart_xpd (l_a,n_na,n_all_a,X_AC,Y_AC,Z_AC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_ca),1,
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr(tmp(i_ca),1,n_na,1,n_all_a,
     &      1,n_na,1,n_all_a,'Cartesian expansion','E',78,4)
*
*    Set up spherical tensors which multiply bessel functions
*
        call ecp_sph_tens (l_sa,n_nsa,n_ta,R_AC,X_AC,Y_AC,Z_AC,
     &      tmp(i_x),tmp(i_y),tmp(i_z),tmp(i_t),tmp(i_Ga),
     &      csco,lcsco)
        if (debug_arrays) call ecp_matpr(tmp(i_Ga),1,n_nsa,1,1,
     &      1,n_nsa,1,1,'Spherical tensors','E',78,4)
      end if
*
*
*   Set up argument values for radial integrals
*
      i_ai = i_free
      i_gam = i_ai+n_abc
      i_alp = i_gam+n_abc
      i_pre = i_alp+n_abc
      i_free = i_pre+n_abc
      if (debug_addresses) write (LuOut,*) 'i_ai,i_gam,i_alp,i_pre',
     &    i_ai,i_gam,i_alp,i_pre
      if (DryRun) then
        mem_max = max(mem_max,i_free-1)
        if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
      else
        if (i_free-1 .gt. l_tmp) call errquit(
     &      ' Insufficient memory in ecp_local1',99, MEM_ERR)
        i = 0
        do i_c = 1,n_coef_c
          do i_b = 1,n_prim_b
            do i_a = 1,n_prim_a
              tmp(i_gam+i) = one/sqrt(zeta_c(i_c)+zeta_b(i_b)
     &            +zeta_a(i_a))
              tmp(i_alp+i) = R_ac*zeta_a(i_a)*tmp(i_gam+i)
              tmp(i_ai+i) = one/(two*R_ac*zeta_a(i_a))
              log_prefactor = tmp(i_alp+i)**2-zeta_a(i_a)*R_ac**2
              tmp(i_pre+i) = exp(log_prefactor)
              i = i+1
            end do
          end do
        end do
      end if
*
*     Loop over on-centre angular momenta
*
      l_max = l_a+l_b
      l_min = l_b
      if (sphcart .eq. 0) then
        l_c_min = mod(l_b,2)
      else
        l_c_min = l_b
      end if
      l_c_max = l_b
      if (debug_addresses) write (LuOut,*) 'l_c_min,l_c_max',
     &    l_c_min,l_c_max
      do l_c = l_c_max,l_c_min,-2
*
*   Set up array dimensions for integrals
*
        ll = min(l_a,l_c)
        n_rad = (ll+1)*(l_a+1)-ll*(ll+1)/2
        if (sphcart .eq. 0) then
          m_count = (l_a+l_b)/2
          ll = max(l_a-l_c,0)
          l = ll/2
          n_rad = n_rad+l*(ll-l)
        else
          m_count = ll
        end if
        m_skip = (l_b-l_c)/2
        l_lo = l_min-m_skip
C        write (LuOut,*) 'm_count',m_count
        n_Qa = l_a+m_skip+1
*
*   Allocate scratch memory for integrals
*
        i_Q_int = i_free
        i_Qa = i_Q_int+ncab*n_rad
        i_Qabc = i_Qa+n_ab*n_Qa
        i_tmp = i_Qabc+n_abc
        i_ind = i_tmp+n_abc*6
        i_wrk = i_ind+n_abc
        i_free = i_wrk+n_ab
        if (debug_addresses) then
          write (LuOut,*) 'i_Q_int,i_Qa,i_Qabc,i_ai,i_gam,i_alp,i_pre,',
     &        'i_tmp,i_ind,i_wrk,i_free'
          write (LuOut,*) i_Q_int,i_Qa,i_Qabc,i_ai,i_gam,i_alp,i_pre,
     &        i_tmp,i_ind,i_wrk,i_free
        end if
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
        else
          if (i_free-1 .gt. l_tmp) call errquit(
     &        ' Insufficient memory in ecp_local1',99, MEM_ERR)
*
*   Calculate radial integrals
*
          call ecp_radint1 (p_min,p_max,
     &        l_lo,l_min,l_max,m_count,m_skip,
     &        n_prim_a,n_cont_a,coef_a,n_prim_b,n_cont_b,coef_b,
     &        n_coef_c,n_prim_c,1,coef_c,
     &        tmp(i_ai),tmp(i_gam),tmp(i_alp),tmp(i_pre),tol,sphcart,
     &        n_ab,nc_ab,n_abc,n_rad,tmp(i_tmp),tmp(i_ind),tmp(i_wrk),
     &        tmp(i_Qabc),tmp(i_Qa),tmp(i_Q_int),transpose,ibug/10)
        end if
*
*   Allocate memory for the contraction of radial and angular parts
*
        n_l_c = 2*l_c+1
        ll = min(l_a,l_c)
        n_ang_a = (ll+1)*(l_a+1)**2-ll*(ll+1)*(2*ll+1)/6
        n_ang_t = n_ang_a*n_l_c
        if (debug_addresses) write (LuOut,*) 'n_ang_a',n_ang_a
        i_ang = i_Qa
        i_ang_a = i_ang+max(n_ang_a*n_nb,n_ang_t)
        i_z = i_ang_a
        i_free = i_ang_a+max(n_ang_t,n_na*n_nb)
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
        else
          if (i_free-1 .gt. l_tmp) call errquit(
     &        ' Insufficient memory in ecp_local1',99, MEM_ERR)
*
*       Set up angular coefficients for centre A and contract over 
*       components of spherical tensors (sum over q).
*
          i = i_ang_a
          do l = 0,l_a
            n_l = 2*l+1
            do k = l+l_c,abs(l-l_c),-2
              i_k = k**2
C              write (LuOut,*) i,loc(tmp(i))
              call ecp_angint (tmp(i),l,k,l_c,tmp(i_Ga+i_k))
              i = i+n_l_c*n_l
            end do
          end do
          if (debug_arrays) call ecp_matpr (tmp(i_ang_a),-l_c,l_c,1,
     &        n_ang_a,-l_c,l_c,1,n_ang_a,'Angular integrals','E',78,4)
*
*       Perform sum over m. This involves transformation of the
*       cartesian functions to a spherical basis on centre B, followed 
*       by evaluation of an overlap integral, which  gives a factor of 
*       4\pi/(2\ell+1). This factor cancels with the factor from the 
*       projector, and the cartesian to spherical transformation 
*       transfers to the other angular integral, leaving a factor of 
*       2\pi(1+\delta_{m,0}).
*
          call dscal (n_ang_a,two,tmp(i_ang_a+l_c),n_l_c)
          fac = pi+pi
          call dscal (n_ang_t,fac,tmp(i_ang_a),1)
          if (debug_arrays) call ecp_matpr (tmp(i_ang_a),-l_c,l_c,1,
     &        n_ang_a,-l_c,l_c,1,n_ang_a,'Scaled angular integrals',
     &        'E',78,4)
          call ecp_cstrans (l_b,n_nb,n_ang_a,l_c,l_c,l,tmp(i_ang),
     &        n_nb,tmp(i_ang_a),n_l_c,csco,lcsco,csco,-1,-1,1)
          n_ang_t = n_ang_a*n_nb
          if (debug_arrays) call ecp_matpr (tmp(i_ang),1,n_nb,1,n_ang_a,
     &        1,n_nb,1,n_ang_a,'Ang ints summed over m','E',78,4)
*
*       Now loop over angular momenta of expanded function a and
*       perform contraction of angular intgrals with radial integrals 
*       and expansion coefficients
*
          i_a = i_ca
          i_c = i_Q_int
          do n = 0,l_a
            i_Qa = i_c
            if (sphcart .eq. 0) then
              l_lo = mod(n,2)
            else
              l_lo = n
            end if
            do l = n,l_lo,-2
              if (debug_gen) write (LuOut,*) 'n,l',n,l
              ll = min(l-1,l_c)
              i_ang_a = i_ang+n_nb*((ll+1)*l**2-ll*(ll+1)*(2*ll+1)/6)
              n_l = 2*l+1
              i_int_a = i_Qa
              do k = l+l_c,abs(l-l_c),-2
                if (debug_addresses) then
                  write (LuOut,*) 'ang coef',i_ang_a-i_ang+1
                  write (LuOut,*) 'cart exp',i_a-i_ca+1
                  write (LuOut,*) 'integral',i_int_a-i_Q_int+1
                end if
                fac = 2*k+1
                if (transpose) then
                  call dgemm ('N','T',n_nb,n_na,n_l,fac,
     &                tmp(i_ang_a),n_nb,tmp(i_a),n_na,zero,
     &                tmp(i_z),n_nb)
                  if (debug_arrays) call ecp_matpr (tmp(i_z),1,n_nb,
     &                1,n_na,1,n_nb,1,n_na,'Completed angular ints',
     &                'F',78,4)
                  call ecp_angrad (n_nb,n_cont_b,n_na,n_cont_a,
     &                tmp(i_z),tmp(i_int_a),ecp_ints)
                else
                  call dgemm ('N','T',n_na,n_nb,n_l,fac,
     &                tmp(i_a),n_na,tmp(i_ang_a),n_nb,zero,
     &                tmp(i_z),n_na)
                  if (debug_arrays) call ecp_matpr (tmp(i_z),1,n_na,
     &                1,n_nb,1,n_na,1,n_nb,'Completed angular ints',
     &                'F',78,4)
                  call ecp_angrad (n_na,n_cont_a,n_nb,n_cont_b,
     &                tmp(i_z),tmp(i_int_a),ecp_ints)
                end if
                i_ang_a = i_ang_a+n_l*n_nb
                i_int_a = i_int_a+ncab
              end do
              i_a = i_a+n_l*n_na
              i_Qa = i_Qa+ncab
            end do
            if (sphcart .eq. 0) then
              i_c = i_c+(min((n+l_c)/2,n)+1)*ncab
            else
              i_c = i_c+(min(n,l_c)+1)*ncab
            end if
          end do
        end if
      end do ! loop over l_c
*
      if (debug_arrays) then
        n_all_a = n_na*n_cont_a
        n_all_b = n_nb*n_cont_b
        call ecp_matpr (ecp_ints,1,n_all_b,1,n_all_a,
     &      1,n_all_b,1,n_all_a,'ECP integrals','E',78,4)
      end if
      if (debug_gen) write (LuOut,*) 'Exiting ecp_local1'
*
      return
      end

