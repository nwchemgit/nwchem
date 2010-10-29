      subroutine fock_eri_scale( ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi, eri )
c     
c     $Id$
c     
      implicit none
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
c     
c     
c     
c     Pre-scale diagonal blocks of the integrals 
c     for permutational symmtery ( a la DEB )
c     
c     
      logical oij, okl, oikjl
      integer i, j, k, l
      integer jtop, ktop, ltop
c     
c     
c     
      oij = (ilo.eq.jlo)
      okl = (klo.eq.llo)
      oikjl = (ilo.eq.klo).and.(jlo.eq.llo)
c$$$  print*,'eri_scale: ',ilo,ihi,jlo,jhi,
c$$$  $                     klo,khi,llo,lhi,
c$$$  $                     eri(llo,klo,jlo,ilo),
c$$$  $                     eri(lhi,khi,jhi,ihi)
c     
c     
c     
      if (oij) then
         do i = ilo,ihi
            ktop = khi
            if (oikjl) ktop = i
            do k = klo,ktop
               ltop = lhi
               if (okl) ltop = k
               do l = llo,ltop
                  eri(l, k, i, i) = eri(l, k, i, i)*0.5d0
               end do
            end do
         end do
      end if
      if (okl) then
         do i = ilo,ihi
            jtop = jhi
            if (oij) jtop = i
            do j = jlo,jtop
               ktop = khi
               if (oikjl) ktop = i
               do k = klo,ktop
                  eri(k, k, j, i) = eri(k, k, j, i)*0.5d0
               end do
            end do
         end do
      end if
      if (oikjl) then
         do i = ilo,ihi
            jtop = jhi
            if (oij) jtop = i
            do j = jlo,jtop
               eri(j, i, j, i) = eri(j, i, j, i)*0.5d0
            end do
         end do
      end if

      return
      end




