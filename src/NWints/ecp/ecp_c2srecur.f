C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_c2srecur (l,Xp,Xo,Xm,ncp,nco,ncm)
*                                                                      *
*   Perform recursion to generate transformation coefficients from     *
*   cartesian monomials to unnormalized solid spherical harmonics      *
*                                                                      *
*   l (inp) - angular momentum of current set. Set to be generated     *
*             has angular momentum l+1                                 *
*   Xp (out) - X_(l+1)^m, set of u.s.s.h. to be generated.             *
*   Xo (inp) - X_l^m, current set of u.s.s.h.                          *
*   Xp (inp) - X_(l-1)^m, previous set of u.s.s.h.                     *
*   ncp (inp) - number of cartesians for l+1 set, (l+2)*(l+3)/2        *
*   nco (inp) - number of cartesians for l set, (l+1)*(l+2)/2          *
*   ncm (inp) - number of cartesians for l-1 set, l*(l+1)/2            *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,j,k,l,m,m1,x,y,z,ncp,nco,ncm
      double precision wa,wb
      double precision Xp(ncp,-l-1:l+1),Xo(nco,-l:l),Xm(ncm,-l+1:l-1)
*
      do m = 1,l
        m1 = m+1
        wa = l+m+1
        x = 0
        do i = l,0,-1
          k = l-i
          do j = k,0,-1
            x = x+1
            y = x+k+1
            z = y+1
            Xp(x,m1) = Xp(x,m1)+wa*Xo(x,m)
            Xp(x,-m1) = Xp(x,-m1)+wa*Xo(x,-m)
            Xp(y,m1) = Xp(y,m1)-wa*Xo(x,-m)
            Xp(y,-m1) = Xp(y,-m1)+wa*Xo(x,m)
            Xp(z,m) = Xp(z,m)+Xo(x,m)
            Xp(z,-m) = Xp(z,-m)+Xo(x,-m)
          end do
        end do
      end do
      wa = l+1
      wb = 2*l+1
      wb = wb/wa
      x = 0
      do i = l,0,-1
        k = l-i
        do j = k,0,-1
          x = x+1
          y = x+k+1
          z = y+1
          Xp(x,1) = Xp(x,1)+wa*Xo(x,0)
          Xp(y,-1) = Xp(y,-1)+wa*Xo(x,0)
          Xp(z,0) = Xp(z,0)+wb*Xo(x,0)
        end do
      end do
      wb = -l
      wb = wb/wa
      x = 0
      do i = l-1,0,-1
        k = l-i
        do j = k-1,0,-1
          x = x+1
          y = x+2*k+1
          z = y+2
          Xp(x,0) = Xp(x,0)+wb*Xm(x,0)
          Xp(y,0) = Xp(y,0)+wb*Xm(x,0)
          Xp(z,0) = Xp(z,0)+wb*Xm(x,0)
        end do
      end do
*
      return
      end
