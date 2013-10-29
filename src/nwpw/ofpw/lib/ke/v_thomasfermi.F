*
* $Id$
*

*     ***********************************************
*     *                                             *
*     *            v_thomasfermi                    *
*     *                                             *
*     ***********************************************

*  Computes the Thomas-Fermi potential and energy density

      subroutine v_thomasfermi(n2ft3d,ispin,dn,xcp,xce)
      implicit none
      integer n2ft3d
      integer ispin
      real*8 dn(n2ft3d,2)
      real*8 xcp(n2ft3d,2)
      real*8 xce(n2ft3d)


      integer ms,k
      real*8  n,nup,ndn,eta,f

      real*8 two3rd,five3rd,twotwothirds,CTF,ffac,dncut
      parameter (two3rd=2.0d0/3.0d0,five3rd=5.0d0/3.0d0)
      parameter (twotwothirds=1.587401052d0)
      parameter (CTF=2.871234d0)
      parameter (ffac=0.851207191959658d0)
      parameter (dncut=1.0d-30)


      call nwpw_timing_start(4)

      call dcopy(n2ft3d,0.0d0,0,xce,1)
      if (ispin.eq.1) then
         do k=1,n2ft3d
            xcp(k,1) = CTF*five3rd*(dn(k,1)+dn(k,ispin))**two3rd
            xce(k)   = xce(k) + CTF*(dn(k,1)+dn(k,ispin))**two3rd
         end do
      else
         do k=1,n2ft3d
            nup = dn(k,1)+0.5d0*dncut
            ndn = dn(k,2)+0.5d0*dncut
            n = nup+ndn
            eta = (nup-ndn)/n
            f = ffac*((1.0d0+eta)**five3rd+(1.0d0-eta)**five3rd-2.0d0)
            xce(k) = xce(k)+(1.0d0+f*(twotwothirds-1.0d0))*CTF*n**two3rd
            xcp(k,1) = twotwothirds*CTF*five3rd*nup**two3rd
            xcp(k,2) = twotwothirds*CTF*five3rd*ndn**two3rd
         end do
      end if

      call nwpw_timing_end(4)

      return
      end
