C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_sph_tens (l,n_n,n_t,R,X,Y,Z,xn,yn,zn,tmp,G_kq,
     &    csco,lcsco)
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
     &    tmp(n_t),G_kq(n_n),tooclose
      double precision csco(lcsco)
      data tooclose/1.0d-14/
*
*   Set up monomials in X, Y and Z
*
      if (R .lt. tooclose) then
        call dfill(n_n,0.0d00,G_kq,1)
        G_kq(1) = 1.0d00
      else
        xn(0) = 1.0d00
        yn(0) = 1.0d00
        zn(0) = 1.0d00
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
*
*     Transform cartesians to sphericals
*
          mmc = (m+1)*(m+2)/2
          mms = 2*m+1
          call ecp_cstrans (m,mmc,1,m,m,k,tmp,mmc,G_kq(ind_s),mms,
     &        csco,lcsco,csco,1,-1,1)
          ind_s = ind_s+mms
        end do
      end if
*
      return
      end

            
