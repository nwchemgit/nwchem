C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_c2slower (l,Xp,Xm,ncp,ncm)
*                                                                      *
*   Multiply spherical harmonic transformation coefficients by r^2 to  *
*   generate coefficients for "contaminant" functions                  *
*                                                                      *
*   l (inp) - angular momentum of set to be generated                  *
*   Xp (out) - set of coefficients be generated.                       *
*   Xm (inp) - set of coefficients to be multiplied by r^2.            *
*   ncp (inp) - number of cartesians for l set, (l+1)*(l+2)/2          *
*   ncm (inp) - number of cartesians for l-2 set, l*(l-1)/2            *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,j,k,l,m,x,y,z,ncp,ncm
      double precision Xp(ncp,ncm),Xm(ncm,ncm)
*
      x = 0
      do i = l-2,0,-1
        k = l-i
        do j = k-2,0,-1
          x = x+1
          y = x+2*k-1
          z = y+2
          do m = 1,ncm
            Xp(x,m) = Xp(x,m)+Xm(x,m)
            Xp(y,m) = Xp(y,m)+Xm(x,m)
            Xp(z,m) = Xp(z,m)+Xm(x,m)
          end do
        end do
      end do
*
      return
      end
