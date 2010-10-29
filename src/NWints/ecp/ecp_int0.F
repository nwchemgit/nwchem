C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_int0 (mem_max,DryRun,
     &    l_a,n_prim_a,n_cont_a,coef_a,zeta_a,n_na,
     &    l_b,n_prim_b,n_cont_b,coef_b,zeta_b,n_nb,
     &    l_ecp,n_prim_c,n_coef_c,
     &    ind_z,ind_c,n_z_tot,n_c_tot,l_ecp_max,n_cent,
     &    zeta_c,coef_c,p_min,p_max,tol,sphcart,tmp,l_tmp,
     &    csco,lcsco,
     &    ecp_ints,n_blk,ibug)
*                                                                      *
*   Calculate Type 0 nonlocal integrals for a given ECP centre         *
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
*   l_ecp (inp) - maximum angular momentum of ECP projectors           *
*   n_prim_c (inp) - number of primitive functions for each power of r *
*                    in ECP expansion                                  *
*   n_coef_c (inp) - array of number of coefficients/exponents for     *
*                    each projector                                    *
*   ind_z (inp) - array of addresses of first exponent for each l      *
*   ind_c (inp) - array of addresses of first coefficient for each l   *
*   n_z_tot (inp) - total number of ECP exponents                      *
*   n_c_tot (inp) - total number of ECP coefficients                   *
*   l_ecp_max (inp) - maximum angular momentum of any ECP projector    *
*   n_cent (inp) - number of ECP centres                               *
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
*   n_blk (inp) - 1 for scalar only, 3 for s-o only, 4 for both        *
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
     &    l_ecp,l_ecp_max,n_z_tot,n_c_tot,n_cent,l_tmp,p_min,p_max,
     &    sphcart,mem_max,n_blk,ibug
      integer n_coef_c(-1:l_ecp_max,n_cent,2),
     &    n_prim_c(p_min:p_max,-1:l_ecp_max,n_cent,2),
     &    ind_z(-1:l_ecp_max,n_cent,2),ind_c(-1:l_ecp_max,n_cent,2)
      integer i,j,m,i_wrk,i_prm,i_tmp,i_Q_int,i_coef,i_scr,i_free,i_ang,
     &    i_ze,i_cont_c,i_off,i_co,l_c,l_c_min,l_c_max,l_cmp,
     &    n_ab,nc_ab,ncab,n_l,n_co,n_co_tot,n_co_max,n_cont_c
      integer n_all_a,n_all_b,n_ang,n_x,n_pass
      integer lcsco
      logical DryRun,debug_gen,debug_addresses,debug_arrays,skip
      double precision zeta_c(n_z_tot),coef_c(n_c_tot),
     &    coef_a(n_prim_a,n_cont_a),coef_b(n_prim_b,n_cont_b),
     &    zeta_a(n_prim_a),zeta_b(n_prim_b),
     &    tmp(l_tmp),ecp_ints(n_na,n_cont_a,n_nb,n_cont_b,n_blk),
     &    tol,wa,wb,wc
      double precision csco(lcsco)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = (mod(ibug,10)/2 .eq. 1) .and. .not.DryRun
*
      if (debug_gen) write (LuOut,'(//A,/)') 'Entering ecp_int0 ...'
*
      if (debug_addresses) then
        write (LuOut,*) 'n_prim_a,n_prim_b,n_cont_a,n_cont_b',
     &      n_prim_a,n_prim_b,n_cont_a,n_cont_b
        write (LuOut,*) 'n_cent,n_c_tot',n_cent,n_c_tot
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
        l_c_max = min(l_a,l_b,l_ecp)
      else
        if (l_a .ne. l_b) return
        if (l_a .gt. l_ecp) return
        l_c_min = l_a
        l_c_max = l_a
      end if
*
*   Loop over ECP angular momentum projectors (l values)
*
      if (debug_addresses) write (LuOut,*) 'l_c_min,l_c_max,l_ecp',
     &    l_c_min,l_c_max,l_ecp
      do l_c = l_c_min,l_c_max,2
*
*     Define parameters for the different integral class cases
*
        call ecp_ind_var (l_c,n_blk,n_coef_c(l_c,1,1),n_coef_c(l_c,1,2),
     &      ind_c(l_c,1,1),ind_c(l_c,1,2),ind_z(l_c,1,1),ind_z(l_c,1,2),
     &      n_x,n_co_tot,n_co_max,i_off,n_pass,i_cont_c,n_cont_c,skip)
        if (debug_addresses) write (LuOut,*) 'l_c,n_x,skip',l_c,n_x,skip
        if (skip) go to 10
*
*     Set up pointers for radial integral evaluation
*
        i_Q_int = 1
        i_tmp = i_Q_int+ncab*max(n_cont_c,n_pass)
        i_prm = i_tmp+n_ab
        i_wrk = i_prm+n_ab*n_co_max
        i_coef = i_wrk+nc_ab
        i_free = i_coef+n_co_tot
        if (debug_addresses) write (LuOut,*) 
     &      'i_Q_int,i_tmp,i_prm,i_wrk,i_coef,i_free',
     &      i_Q_int,i_tmp,i_prm,i_wrk,i_coef,i_free
        if (DryRun) then
          mem_max = max(mem_max,i_free-1)
          if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
        else
          if (i_free-1 .gt. l_tmp) call errquit(
     &        ' Insufficient memory in ecp_int0',99, MEM_ERR)
          i = i_cont_c
          j = i_Q_int
          do m = 1,n_pass
            i_ze = ind_z(l_c,1,i)
            i_co = ind_c(l_c,1,i)
            n_co = n_coef_c(l_c,1,i)
            call dcopy (n_co,coef_c(i_co),1,tmp(i_coef),1)
            if (n_cont_c .eq. 2) call dcopy 
     &          (n_co,coef_c(ind_c(l_c,1,2)),1,tmp(i_coef+n_co),1)
            if (debug_addresses) write (LuOut,*) 
     &          'l_c,i_ze,i_co,n_co',
     &          l_c,i_ze,i_co,n_co
            call ecp_radint0 (p_min,p_max,
     &          n_prim_c(p_min,l_c,1,i),n_co,n_cont_c,
     &          zeta_c(i_ze),tmp(i_coef),
     &          l_a,n_prim_a,n_cont_a,zeta_a,coef_a,
     &          l_b,n_prim_b,n_cont_b,zeta_b,coef_b,
     &          tmp(i_prm),tmp(i_tmp),tmp(i_wrk),tmp(j),
     &          (ibug/10))
            i = i+1
            j = j+ncab
          end do
        end if
*
        n_l = 2*l_c+1
        if (sphcart .eq. 0) then
*
*       Cartesian integrals.
*
          wa = n_l
          wa = pi/wa
          wb = wa+wa
          if (debug_addresses) write (LuOut,*) 'n_l,n_na,n_x',
     &        n_l,n_na,n_x
          i_ang = i_tmp+n_l*n_na
          n_ang = n_l**2
          i_scr = i_ang+n_ang*n_x
          i_free = i_scr+max(n_na*n_nb,n_ang)
          if (debug_addresses) write (LuOut,*) 'i_ang,i_free',i_ang,
     &        i_free
          if (DryRun) then
            mem_max = max(mem_max,i_free-1)
            if (debug_addresses) write (LuOut,*) 'mem_max',mem_max
          else
            if (i_free-1 .gt. l_tmp) call errquit(
     &          ' Insufficient memory in ecp_int0',99, MEM_ERR)
*
*       Set up angular overlap/spin-orbit matrix
*
            call dcopy (n_ang*n_x,zero,0,tmp(i_ang),1)
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
*         spin-orbit potentials 
*
            if ((n_x .ne. 1) .and. (l_c .gt. 0)) then
*
*           z spin-orbit potential 
*
              j = i
              do m = -l_c,l_c
                j = j+l_c+l_c
                tmp(j) = m
                tmp(j) = tmp(j)*wb
              end do
*
*           y,x spin-orbit potentials
*
              i = i+n_ang
              j = i+n_ang+l_c+l_c-1
              do m = -l_c,l_c-1
                wc = (l_c-m)*(l_c+m+1)
                wc = sqrt(wc)*wa
                if ((m .eq. 0) .or. (m .eq. -1)) wc = wc+wc
                if (m .ge. 0) wc = -wc
                if (m .eq. -1) then
                  i = i+n_l+1
                else
                  i = i+1
                  tmp(i) = -wc
                  i = i+l_c+l_c
                  tmp(i) = wc
                  i = i+1
                end if
                tmp(j) = wc
                j = j+l_c+l_c
                if ((m .lt. -1) .or. (m .gt. 0)) then
                  tmp(j+2) = wc
                end if
              end do
            end if
*
            i = i_ang
            do j = 1,n_x
              m = j+i_off
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
              if (debug_arrays) call ecp_matpr (tmp(i_Q_int),1,n_cont_a,
     &            1,n_cont_b,1,n_cont_a,1,n_cont_b,'Radial integrals',
     &            'F',78,4)
              call ecp_angrad (n_na,n_cont_a,n_nb,n_cont_b,
     &            tmp(i_scr),tmp(i_Q_int),ecp_ints(1,1,1,1,m))
              if (debug_arrays) call ecp_matpr (ecp_ints(1,1,1,1,m),
     &            1,n_all_b,1,n_all_a,1,n_all_b,1,n_all_a,
     &            'Assembled integrals','E',78,4)
              if ((j .eq. 1) .and. (n_x .eq. 4)) 
     &            i_Q_int = i_Q_int+ncab
              i = i+n_ang
            end do
          end if
        else
          if (n_blk .gt. 1) call errquit (
     &        'Cannot do spherical ecp spin-orbit type 0 integrals',
     &        99, INT_ERR)
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
      if (debug_arrays) then
        do j = 1,n_blk
          call ecp_matpr (ecp_ints(1,1,1,1,j),1,n_all_b,1,n_all_a,
     &        1,n_all_b,1,n_all_a,'ECP integrals','E',78,4)
        end do
      end if
      if (debug_gen) write (LuOut,*) 'Exiting ecp_int0'
*
      return
      end
