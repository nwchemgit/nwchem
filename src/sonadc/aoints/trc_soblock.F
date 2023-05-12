      subroutine trc_soblock(H, ilo, ihi, jlo, jhi, iat, icart,
     &  socomp, nat, iabf, jabf, max_at_bf, frc, dblk)

      integer  ilo, ihi, jlo, jhi, iat, icart, socomp, nat, iabf, jabf,
     &  max_at_bf
      double precision H(jlo:jhi, ilo:ihi, 3, 3, nat)
      double precision frc(3,nat)
      double precision dblk(max_at_bf,max_at_bf)
      character*2 dcart(3)
      data dcart/'dx','dy','dz'/
      integer i, j, ip, jp
      double precision de

      de = 0d0
      do i=ilo, ihi
        ip=i-iabf+1
        do j=jlo, jhi
          jp=j-jabf+1
          de=de+H(j,i,socomp,icart,iat)*dblk(ip,jp)
*         write(6,2000) iat,dcart(icart),'<',i,'|Hso',socomp,'|',j,'>=',
*    &      H(j,i,socomp,icart,iat), ' den=',dblk(ip,jp)
 2000 format(i2,a2,3(a,i2),2(a,f20.10))
        enddo
      enddo
      frc(icart,iat)=frc(icart,iat)+de
*     write(*,'("forces",i2,9e13.5)')socomp,((frc(i,j),i=1,3),j=1,nat)
      return
      end
