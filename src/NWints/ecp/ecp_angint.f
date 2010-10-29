C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_angint (angint,lambda,k,l,G_kq)
*                                                                      *
*   Set up array of 3j products, contract over spherical tensors.      *
*   This is (in TeX) \sum_q [(1+\delta_{mu,0})/(1+\delta_{\rho,0})]    *
*      \times {\cal B}^{\lambda k \ell}_{\mu q m}  G_{kq}              *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer lambda,k,l,mu,q,m
      double precision angint(-l:l,-lambda:lambda),G_kq(-k:k),wa
*
C      write (6,*) 'ecp_angint'
C      write (6,*) loc(angint(-l,-lambda))
C      write (0,*) lambda,k,l
      do m = 1,l
        do mu = 1,lambda
          angint(m,mu) = 0.0d00
          angint(m,-mu) = 0.0d00
          angint(-m,mu) = 0.0d00
          angint(-m,-mu) = 0.0d00
          q = m+mu
          if (q .le. k) then
            call ecp_3j_prod (lambda,l,k,mu,m,wa)
C            write (6,'(1PE15.7)') wa
C            write (6,'(1P2E15.7)') G_kq(q),G_kq(-q)
            angint(m,mu) = angint(m,mu)+wa*G_kq(q)
            angint(-m,-mu) = angint(-m,-mu)-wa*G_kq(q)
            if (q .ne. 0) then
              angint(m,-mu) = angint(m,-mu)+wa*G_kq(-q)
              angint(-m,mu) = angint(-m,mu)+wa*G_kq(-q)
            end if
C            write (6,'(1P2E15.7)') angint(m,mu),angint(m,-mu)
C            write (6,'(1P2E15.7)') angint(-m,mu),angint(-m,-mu)
          end if
          q = abs(mu-m)
          if (q .le. k) then
            call ecp_3j_prod (lambda,l,k,mu,-m,wa)
C            write (6,'(1PE15.7)') wa
C            write (6,'(1P2E15.7)') G_kq(q),G_kq(-q)
            angint(m,mu) = angint(m,mu)+wa*G_kq(q)
            angint(-m,-mu) = angint(-m,-mu)+wa*G_kq(q)
            if (q .ne. 0) then
              if (mu .gt. m) wa = -wa
              angint(m,-mu) = angint(m,-mu)-wa*G_kq(-q)
              angint(-m,mu) = angint(-m,mu)+wa*G_kq(-q)
            end if
C            write (6,'(1P2E15.7)') angint(m,mu),angint(m,-mu)
C            write (6,'(1P2E15.7)') angint(-m,mu),angint(-m,-mu)
          end if
        end do
        angint(m,0) = 0.0d00
        angint(-m,0) = 0.0d00
        q = m
        if (q .le. k) then
          call ecp_3j_prod (lambda,l,k,0,m,wa)
C          write (6,'(1PE15.7)') wa
          wa = wa+wa
C          write (6,'(1P2E15.7)') G_kq(q),G_kq(-q)
          angint(m,0) = wa*G_kq(q)
          angint(-m,0) = wa*G_kq(-q)
C          write (6,'(1PE15.7)') angint(m,0)
C          write (6,'(1PE15.7)') angint(-m,0)
        end if
      end do
      do mu = 1,lambda
        angint(0,mu) = 0.0d00
        angint(0,-mu) = 0.0d00
        q = mu
        if (q .le. k) then
          call ecp_3j_prod (lambda,l,k,mu,0,wa)
C          write (6,'(1PE15.7)') wa
C          write (6,'(1P2E15.7)') G_kq(q),G_kq(-q)
          angint(0,mu) = wa*G_kq(q)
          angint(0,-mu) = wa*G_kq(-q)
C          write (6,'(1P2E15.7)') angint(0,mu),angint(0,-mu)
        end if
      end do
      call ecp_3j_prod (lambda,l,k,0,0,wa)
C      write (6,'(1PE15.7)') wa
C      write (6,'(1PE15.7)') G_kq(0)
      angint(0,0) = wa*G_kq(0)
C      write (6,'(1PE15.7)') angint(0,0)
*
      return
      end
