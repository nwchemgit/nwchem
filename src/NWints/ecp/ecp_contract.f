C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_contract (n_ab,n_c,m_count,Qk,coeff_c,Q)
*                                                                      *
*   Contract ECP radial integrals over potential expansion             *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n_ab (inp) - product of numbers of functions on centres a and b    *
*   n_c (inp) - number of functions on centre c (ecp centre)           *
*   m_count (inp) - number of m values                                 *
*   Qk (inp) - uncontracted ECP integrals                              *
*   coeff_c (inp) - coefficients of functions on centre c              *
*   Q (out) - contracted ECP integrals                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i_ab,i_c,m,m_count,n_ab,n_c
      double precision Qk(n_ab,n_c,m_count),Q(n_ab,m_count),
     &    coeff_c(n_c)
*
      do m = 1,m_count
        do i_ab = 1,n_ab
          Q(i_ab,m) = 0.0d00
        end do
        do i_c = 1,n_c
C          write (6,'(1p4e20.10)') (Qk(i_ab,i_c,m),i_ab = 1,n_ab),
C     &        coeff_c(i_c)
          do i_ab = 1,n_ab
            Q(i_ab,m) = Q(i_ab,m)+Qk(i_ab,i_c,m)*coeff_c(i_c)
          end do
        end do
      end do
*
      return
      end
