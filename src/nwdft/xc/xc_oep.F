      Subroutine read_vect(rho, vect1, nq, ipol, isp)
      implicit none
      integer n, nq, ipol, isp
      double precision rho(nq,ipol*(ipol+1)/2), vect1(nq)
      do n = 1,nq
        if (isp.eq.1) then
          vect1(n) = rho(n,2)
        else
          vect1(n) = rho(n,3)
        end if
      end do
      return
      end

      Subroutine read_delrho(delrho, delrho_sig, nq, ipol, isp)
      implicit none
      integer n, nq, ipol, isp
      double precision delrho(nq,3,ipol), delrho_sig(nq,3)
      do n = 1,nq
        delrho_sig(n,1) = delrho(n,1,isp)
        delrho_sig(n,2) = delrho(n,2,isp)
        delrho_sig(n,3) = delrho(n,3,isp)
      end do
      return
      end

      Subroutine part_rho_rel(tol_rho, nq, Amat, Cmat, rho, vect1,
     &                        ipol, occup_i, del_rho, delrho_sig,
     &                        isp, aux_ik)
      implicit none
      integer n, nq, ipol, isp, aux_ik
      double precision tol_rho, occup_i, compx, compy, compz,
     &                 Amat(nq,ipol), vect1(nq), Cmat(nq,3,ipol),
     &                 rho(nq,ipol*(ipol+1)/2),
     &                 del_rho(nq,3,ipol), delrho_sig(nq,3),
     &                 relation, rhosig
      do n = 1,nq
        rhosig = vect1(n)
        if (rhosig.gt.tol_rho) then
          relation = occup_i*rho(n,2)/rhosig
          compx = occup_i*del_rho(n,1,1)/rhosig -
     &                     rho(n,2)*delrho_sig(n,1)/(rhosig*rhosig)
          compy = occup_i*del_rho(n,2,1)/rhosig -
     &                     rho(n,2)*delrho_sig(n,2)/(rhosig*rhosig)
          compz = occup_i*del_rho(n,3,1)/rhosig -
     &                     rho(n,2)*delrho_sig(n,3)/(rhosig*rhosig)
        else
          if (aux_ik.eq.1) then
            relation = 1.0d00
          else
            relation = 0.d00
          end if
          compx = 0.0d00
          compy = 0.0d00
          compz = 0.0d00
        end if
        compx = compx*Cmat(n,1,1)
        compy = compy*Cmat(n,2,1)
        compz = compz*Cmat(n,3,1)
        Amat(n,1) = Amat(n,1)*relation + compx + compy + compz
        Amat(n,2) = Amat(n,2)*relation
        Cmat(n,1,1) = Cmat(n,1,1)*relation
        Cmat(n,2,1) = Cmat(n,2,1)*relation
        Cmat(n,3,1) = Cmat(n,3,1)*relation
        Cmat(n,1,2) = 0.0d00
        Cmat(n,2,2) = 0.0d00
        Cmat(n,3,2) = 0.0d00
      end do
      return
      end

      Subroutine total_Amat( nq, ipol, Amat, Amat_coul, qwght)
      implicit none
      integer n, ipol, nq
      double precision Amat(nq,ipol), Amat_coul(nq), qwght(nq)
      do n = 1,nq
        Amat(n,1) = Amat(n,1) - Amat_coul(n)*qwght(n)
        Amat(n,2) = qwght(n)
      end do
      return
      end
c $Id$
