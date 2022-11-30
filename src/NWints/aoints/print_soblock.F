      subroutine print_soblock(H, ilo, ihi, jlo, jhi, iat, icart,
     &  socomp, nat)

      integer  ilo, ihi, jlo, jhi, iat, icart, socomp, nat
      double precision H(jlo:jhi, ilo:ihi, 3, 3, nat)
      character*2 dcart(3)
      data dcart/'dx','dy','dz'/
      integer i, j

      do i=ilo, ihi
        do j=jlo, jhi
          write(6,2000) iat,dcart(icart),'<',i,'|Hso',socomp,'|',j,'>=',
     &      H(j,i,socomp,icart,iat)
 2000 format(i2,a2,3(a,i2),a,f20.10)
        enddo
      enddo
      return
      end
