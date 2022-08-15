      subroutine bse_get_eia(eig,eia,nocc,nvir,llpole,ulpole)
      implicit none
      integer,intent(in) :: nocc,nvir
      double precision,intent(in),dimension(nocc+nvir) :: eig
      double precision,intent(out),dimension(*) :: eia
      integer ipole,llpole,ulpole,imo,amo
c
      do ipole=llpole,ulpole
        imo = (ipole-1)/nvir + 1
        amo = ipole - (imo-1)*nvir + nocc
        eia(ipole) = eig(amo) - eig(imo)
      enddo
c      
      return
      end subroutine
