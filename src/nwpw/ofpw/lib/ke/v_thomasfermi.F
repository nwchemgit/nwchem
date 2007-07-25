*
* $Id: v_thomasfermi.F,v 1.2 2007-07-25 01:29:57 bylaska Exp $
*

*  Computes the Thomas-Fermi potential and energy density

      subroutine v_thomasfermi(n2ft3d,ispin,dn,xcp,xce)
      implicit double precision(a-h, o-z)
      implicit integer (i-n)

      integer n2ft3d
      integer ispin
      real*8 dn(n2ft3d,2)
      real*8 xcp(n2ft3d,2)
      real*8 xce(n2ft3d)


      parameter (two3rd=2.0d0/3.0d0,five3rd=5.0d0/3.0d0)
      parameter (twotwothirds=1.587401052d0)
      parameter (CTF=2.871234d0)
      parameter (dncut=1.0d-30)


      call nwpw_timing_start(4)

      if (ispin.eq.1) then
         do k=1,n2ft3d
            xcp(k,1) = CTF*five3rd*(dn(k,1)**two3rd)
            xce(k)   = CTF*(dn(k,1)+dn(k,ispin))**two3rd
         end do

      else
         call dcopy(n2ft3d,0.0d0,0,xce,1)
         do ms=1,ispin
         do k=1,n2ft3d
            xcp(k,ms) = twotwothirds*CTF*five3rd*(dn(k,ms)**two3rd)
            xce(k) = xce(k) + twotwothirds*CTF*(dn(k,ms))**two3rd
         end do
         end do

      end if

      call nwpw_timing_end(4)

      return
      end
