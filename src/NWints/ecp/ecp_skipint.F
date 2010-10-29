C $Id$
************************************************************************
*                                                                      *
      logical function ecp_skipint(
     &    l_a,n_prim_a,n_cont_a,coef_a,zeta_a,R_AC,
     &    l_b,n_prim_b,n_cont_b,coef_b,zeta_b,R_BC,
     &    n_coef_c,zeta_c,coef_c)
*                                                                      *
*   Do screening check for ECP integrals                               *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   l_a (inp) - (maximum) angular momentum of functions on centre A    *
*   n_prim_a (inp) - number of primitive functions on centre A         *
*   n_cont_a (inp) - number of contracted functions on centre A        *
*   coef_a (inp) - centre A contraction coefficients                   *
*   zeta_a (inp) - centre A exponents                                  *
*   R_AC (inp) - distance between centres A and C                      *
*   l_b (inp) - (maximum) angular momentum of functions on centre B    *
*   n_prim_b (inp) - number of primitive functions on centre B         *
*   n_cont_b (inp) - number of contracted functions on centre B        *
*   coef_b (inp) - centre B contraction coefficients                   *
*   zeta_b (inp) - centre B exponents                                  *
*   R_BC (inp) - distance between centres B and C                      *
*   n_coef_c (inp) - number of coefficients/exponents for ECP          *
*   zeta_c (inp) - ECP exponents                                       *
*   coef_c (inp) - ECP contraction coefficients                        *
*                                                                      *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer l_a,n_prim_a,n_cont_a
      integer l_b,n_prim_b,n_cont_b
      Integer n_coef_c
      integer i_a,i_b,i_c
      double precision zeta_c(n_coef_c)
      double precision coef_c(n_coef_c)
      double precision coef_a(n_prim_a,n_cont_a)
      double precision coef_b(n_prim_b,n_cont_b)
      double precision zeta_a(n_prim_a)
      double precision zeta_b(n_prim_b)
      double precision R_AB
      double precision R_AC
      double precision R_BC
      double precision log_prefac,pr,pre,pref,sum
      double precision a,b,c,ca,cb,aA,bB,ab,aAbB,abc,P,wa,wb
      logical first
*
      integer idamax
      external idamax
*
      wa = l_a
      wb = l_b
      first = .true.
      log_prefac = zero
*
      R_AB = abs(R_AC-R_BC)
      do i_b = 1,n_prim_b
        b = zeta_b(i_b)
        i_c = idamax(n_cont_b,coef_b(i_b,1),n_prim_b)
        if (abs(coef_b(i_b,i_c)) .ne. zero) then
          cb = log(abs(coef_b(i_b,i_c)))
        else
          cb = zero
        end if
        bB = b*R_BC 
        do i_a = 1,n_prim_a
          a = zeta_a(i_a)
          i_c = idamax(n_cont_a,coef_a(i_a,1),n_prim_a)
          if (abs(coef_a(i_a,i_c)) .ne. zero) then
            ca = log(abs(coef_a(i_a,i_c)))
          else
            ca = zero
          end if
          aA = a*R_AC 
          aAbB = aA+bB
          ab = zeta_a(i_a)+zeta_b(i_b)
          pref = zero
          sum = zero
          do i_c = 1,n_coef_c
            c = zeta_c(i_c)
            abc = ab+c
            P = aAbB/abc
            pre = 0.5d0*(log(pi)-log(abc))
     &          - (a*b*R_AB**2 + b*c*R_BC**2 + a*c*R_AC**2)/abc
            if (abs(R_AC-P) .gt. 1.0d-10) pre = pre+wa*log(abs(R_AC-P))
            if (abs(R_BC-P) .gt. 1.0d-10) pre = pre+wb*log(abs(R_BC-P))
            if (pref .eq. zero) then
              pref = pre
              sum = coef_c(i_c)
            else 
              pr = pre - pref
              if (pr .lt. zero) then
                if (pr .gt. ln_thr_ecp) 
     &              sum = sum+coef_c(i_c)*exp(pr)
              else
                if (-pr .gt. ln_thr_ecp) then
                  sum = sum*exp(-pr)
                else
                  sum = zero
                end if
                sum = sum+coef_c(i_c)
                pref = pre
              end if
            end if
          end do
          if (sum .ne. zero) pref = pref+log(abs(sum))
          pref = pref+ca+cb
          if (first) then
            log_prefac = pref
            first = .false.
          else
            log_prefac = max(log_prefac, pref)
          end if
        end do
      end do
*
      log_prefac = log_prefac + log(dble(n_prim_a*n_prim_b))
      ecp_skipint = log_prefac .lt. ln_thr_ecp
*
      return
      end


