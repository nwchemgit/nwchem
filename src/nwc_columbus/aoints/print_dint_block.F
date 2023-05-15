      subroutine print_dint_block(ilo, ihi, jlo, jhi,
     &  klo, khi, llo, lhi, geri, idatom)

      integer  ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision geri(llo:lhi, klo:khi, jlo:jhi, ilo:ihi, 3, 4)
      integer idatom(4)
      character*4 dcart(3)
      data dcart/'dx','dy','dz'/
      integer i, j, k, l, icart, iat

      do iat = 1, 4
       if(idatom(iat).ge.1) then
        do i=ilo, ihi
         do j=jlo, jhi
          do k=klo, khi
           do l=llo, lhi
            do icart = 1, 3
             write(6,2000) idatom(iat),dcart(icart),'[',i,',',j,';',
     &        k,',',l,'] = ', geri(l,k,j,i,icart,iat)
 2000 format(i2,a2,4(a1,i2),a,f20.10)
            enddo
           enddo
          enddo
         enddo
        enddo
       endif
      enddo
      return
      end
