C $Id: ecp_sph_tens.f,v 1.1 1996-09-30 19:29:29 d3e129 Exp $
************************************************************************
*                                                                      *
      subroutine ecp_sph_tens (l,n_n,n_t,R,X,Y,Z,xn,yn,zn,tmp,G_kq,
     &    csco,lcsco)
Cold     &    c_to_s,s_to_c,ldt,lstart)
*                                                                      *
*   Set up spherical tensors which arise from expansion of exponential *
*   about a new centre. Limit on angular momentum depends on projector *
*   on new centre.                                                     *
*                                                                      *
*   l (inp) - maximum angular momentum of spherical tensor
*   n_n (inp) - number of components of spherical tensor = (n+l+1)**2  *
*   n_t (inp) - dimension of work area = (n+l+1)*(n+l+2)/2             *
*   R (inp) - distance to new centre                                   *
*   X,Y,Z (inp) - relative cartesian coordinates of new centre         *
*   xn,yn,zn (scr) - work arrays to store powers of X,Y,Z              *
*   tmp (scr) - work array for transformation                          *
*   G_kq (out) - array of spherical tensors                            *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,j,k,l,m,ind_c,ind_s,n_n,n_t,mmc,mms
      integer lcsco
      double precision X,Y,Z,R,xr,yr,zr,xn(0:l),yn(0:l),zn(0:l),
     &    tmp(n_t),G_kq(n_n)
      double precision csco(lcsco)
*
*   Set up monomials in X, Y and Z
*
      xn(0) = 1
      yn(0) = 1
      zn(0) = 1
      xr = -X/R
      yr = -Y/R
      zr = -Z/R
      do i = 1,l
        xn(i) = xn(i-1)*xr
        yn(i) = yn(i-1)*yr
        zn(i) = zn(i-1)*zr
      end do
*
*   Loop over angular momenta on new centre
*
      ind_s = 1
      do m = 0,l
*
*     Loop over cartesian components on new centre
*
        ind_c = 0
        do i = m,0,-1
          do j = m-i,0,-1
            k = m-i-j
            ind_c = ind_c+1
            tmp(ind_c) = xn(i)*yn(j)*zn(k)
          end do
        end do
C        call ecp_matp(tmp,1,ind_c,'Cartesian tensor',81,5)
*
*     Transform cartesians to sphericals
*
        mmc = (m+1)*(m+2)/2
        mms = 2*m+1
Cnew
        call ecp_cstrans (m,mmc,1,m,m,k,tmp,mmc,G_kq(ind_s),mms,
     &      csco,lcsco,csco,1,-1,1)
C        call ecp_matp(G_kq(ind_s),1,mms,'Spherical tensor',81,5)
Cold
C        call ecp_cartsph (m,mmc,1,tmp,m,m,mms,G_kq(ind_s),1,-1,1)
Cend
        ind_s = ind_s+mms
      end do
*
      return
      end

            
