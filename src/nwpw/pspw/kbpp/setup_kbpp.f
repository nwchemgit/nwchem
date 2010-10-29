*
* $Id$
*

*     *****************************************************
*     *                                                   *
*     *                setup_kbpp                         *
*     *                                                   *
*     *****************************************************

      subroutine setup_kbpp(nfft1,nfft2,nfft3,unita,unitg,G)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3),unitg(3,3)
      double precision G(nfft1/2+1,nfft2,nfft3,3)

*     **** local variables ***** 
      integer i,j,k
      integer k1,k2,k3
      integer nffth1,nffth2,nffth3
      double precision g1,g2,g3


      call get_unitg(unita,unitg)

      nffth3 = nfft3/2
      nffth2 = nfft2/2
      nffth1 = nfft1/2
*     **** less confusing algorithm ****
      do k3 = -nffth3+1, nffth3
         do k2 = -nffth2+1, nffth2
            do k1 = 0,nffth1
               g1 = k1*unitg(1,1) + k2*unitg(1,2) +k3*unitg(1,3)
               g2 = k1*unitg(2,1) + k2*unitg(2,2) +k3*unitg(2,3)
               g3 = k1*unitg(3,1) + k2*unitg(3,2) +k3*unitg(3,3)
               i=k1
               j=k2
               k=k3
c              if (i .lt. 0) i = i + nfft1
               if (j .lt. 0) j = j + nfft2
               if (k .lt. 0) k = k + nfft3

               G(i+1,j+1,k+1,1) = g1
               G(i+1,j+1,k+1,2) = g2
               G(i+1,j+1,k+1,3) = g3

            end do  
         end do
      end do


      return
      end




      subroutine get_unitg(unita,unitg)
      implicit none

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors                                *
*               in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             unita  --- primitive vectors in coordination space             *
*                                                                            *
*     Outputs:                                                               *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  dscal from BLAS                                              *
*                                                                            *
*     Last modification:  3/30/99  by Eric Bylaska                           *
*                                                                            *
******************************************************************************


*     ------------------
*     argument variables
*     ------------------
      double precision unita(3,3), unitg(3,3)

*     ---------------
*     local variables
*     ---------------
      double precision volume
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)


*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume=dabs(volume)

      return
      end
