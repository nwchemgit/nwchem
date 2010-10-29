C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_cstrans (l,ncart,ntr,l_hi,l_lo,nsph,cart,ldc,
     &    sph,lds,csco,lcsco,lstart,indicator,position,mode)
*                                                                      *
*   Transform from cartesian monomials x^i y^j z^k to real solid       *
*   spherical tensor functions r^n G_lm(theta,phi). These are          *
*   sqrt(2pi(1+delta_m,0)/2l+1) times the real spherical harmonics.    *
*                                                                      *
*   l (inp) - angular momentum of shell, l = i+j+k                     *
*   ncart (inp) - number of cartesian components = (l+1)*(l+2)/2       *
*   ntr (inp) - number of quantities to transform                      *
*   l_hi - maximum angular momentum of spherical quantities            *
*   l_lo - minimum angular momentum of spherical quantities            *
*   nsph (out) - number of spherical components                        *
*   cart (i/o) - array of cartesian quantities.                        *
*   ldc (inp) - leading dimension of cart                              *
*   sph (i/o) - array of spherical quantities                          *
*   lds (inp) - leading dimension of sph                               *
*   csco (inp) - transformation arrays                                 *
*   lcsco (inp) - length of transformation arrays                      *
*   lstart - pointers to beginning of transformation for given l       *
*   indicator (inp) - +1 = cartesian to spherical,                     *
*                     -1 = spherical to cartesian                      *
*   position (inp) - +1 = post-multiply, array is ntr*(ncart or nsph)  *
*                    -1 = premultiply, array is (ncart or nsph)*ntr    *
*   mode (inp) - +1 = transform with regular matrix,                   *
*                -1 = transform with inverse matrix                    *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer l,ncart,ntr,l_hi,l_lo,nsph,ldc,lds,lcsco,lstart(0:1,0:l),
     &    indicator,position,mode
      integer i,j,offset_start,offset_end
      double precision cart(ldc,ncart),sph(lds,ncart)
      double precision csco(lcsco),zero,one
      parameter (zero = 0.0d00, one = 1.0d00)
*
*   Test validity of input parameters
*
      if (abs(indicator) .ne. 1) call errquit (
     &      'Illegal value of indicator in ecp_cstrans',99, BASIS_ERR)
      if (abs(position) .ne. 1) call errquit (
     &      'Illegal value of position in ecp_cstrans',99, BASIS_ERR)
      if (abs(mode) .ne. 1) call errquit (
     &      'Illegal value of mode in ecp_cstrans',99, BASIS_ERR)
*
C      write (6,*) 'l_hi,l_lo',l_hi,l_lo
      offset_start = (l_hi+1)*(l_hi+2)/2
      offset_end = (l_lo-1)*(l_lo)/2
C      write (6,*) 'offset_start,offset_end',offset_start,offset_end
      nsph = offset_start-offset_end
*
      if (position .gt. 0) then
        if (lds .lt. ntr) call errquit (
     &      ' lds < ntr in ecp_cstrans',99, BASIS_ERR)
        if (ldc .lt. ntr) call errquit (
     &      ' ldc < ntr in ecp_cstrans',99, BASIS_ERR)
      else
        if (lds .lt. nsph) call errquit (
     &      ' lds < nsph in ecp_cstrans',99, BASIS_ERR)
        if (ldc .lt. ncart) call errquit (
     &      ' ldc < ncart in ecp_cstrans',99, BASIS_ERR)
      end if
*
      i = lstart(0,l)+ncart*(ncart-offset_start)
      j = lstart(1,l)+ncart-offset_start
C      write (6,*) lstart(0,l),lstart(1,l)
C      write (6,*) 'l,i,j',l,i,j
C      call ecp_matpr(csco(i),1,ncart,1,ncart,1,ncart,1,nsph,
C     &    'Cart to sph','F',81,5)
C      call ecp_matpr(csco(j),1,ncart,1,ncart,1,nsph,1,ncart,
C     &    'Sph to cart','F',81,5)
*
*     Cartesian to spherical
*
      if (indicator .gt. 0) then
        if (position .gt. 0) then
          if (mode .gt. 0) then
            call dgemm ('N','N',ntr,nsph,ncart,one,cart,ldc,
     &          csco(i),ncart,zero,sph,lds)
          else
            call dgemm ('N','T',ntr,nsph,ncart,one,cart,ntr,
     &          csco(j),ncart,zero,sph,lds)
          end if
        else
          if (mode .gt. 0) then
            call dgemm ('T','N',nsph,ntr,ncart,one,csco(i),
     &          ncart,cart,ldc,zero,sph,lds)
          else
            call dgemm ('N','N',nsph,ntr,ncart,one,csco(j),
     &          ncart,cart,ldc,zero,sph,lds)
          end if
        end if
*
*     Spherical to cartesian
*
      else
        if (position .gt. 0) then
          if (mode .gt. 0) then
            call dgemm ('N','N',ntr,ncart,nsph,one,sph,lds,
     &          csco(j),ncart,zero,cart,ldc)
          else
            call dgemm ('N','T',ntr,ncart,nsph,one,sph,lds,
     &          csco(i),ncart,zero,cart,ldc)
          end if
        else
          if (mode .gt. 0) then
            call dgemm ('T','N',ncart,ntr,nsph,one,csco(j),
     &          ncart,sph,lds,zero,cart,ldc)
          else
            call dgemm ('N','N',ncart,ntr,nsph,one,csco(i),
     &          ncart,sph,lds,zero,cart,ldc)
          end if
        end if
      end if
*
      return
      end
