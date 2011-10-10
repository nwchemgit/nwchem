*     ***************************************
*     *										*
*     *      fcoord_to_real					*
*     *										*
*     ***************************************
*
*     This routine converts from fractional coordinates to real
*     coordinates.
*
*     Entry - nion
*             a(3,3): lattice vectors
*             ion: fractional coordinates
*    Exit -
*             ion: real coordinates
*

      subroutine fcoord_to_real(nion,ion)
      implicit none
      integer nion
      real*8 ion(3,*)

c     **** local variables ****
      integer i,j
      real*8 tion(3)
      real*8 a(3,3)

*     *** external functions ****
      real*8   lattice_unita
      external lattice_unita
    
      do j=1,3
      do i=1,3
        a(i,j) = lattice_unita(i,j)
      end do
      end do

      do i=1,nion
         tion(1) = ion(1,i)
         tion(2) = ion(2,i)
         tion(3) = ion(3,i)

         ion(1,i) = a(1,1)*tion(1)
     >            + a(1,2)*tion(2)
     >            + a(1,3)*tion(3)
         ion(2,i) = a(2,1)*tion(1)
     >            + a(2,2)*tion(2)
     >            + a(2,3)*tion(3)
         ion(3,i) = a(3,1)*tion(1)
     >            + a(3,2)*tion(2)
     >            + a(3,3)*tion(3)
      end do

      return
      end


*     ***************************************
*     *										*
*     *      fcoord_to_frac 				*
*     *										*
*     ***************************************
*
*     This routine converts from real coordinates to fractional 
*     coordinates.
*
*     Entry - nion
*             a(3,3): lattice vectors
*             ion: real coordinates
*     Exit -
*             ion: fractional coordinates

      subroutine fcoord_to_frac(nion,ion)
      implicit none
      integer nion
      real*8 ion(3,*)

*     **** local variables ****
      integer i,j
      real*8 a(3,3),b(3,3),volume
      real*8 tion(3)

*     *** external functions ****
      real*8   lattice_unita
      external lattice_unita
    
      do j=1,3
      do i=1,3
        a(i,j) = lattice_unita(i,j)
      end do
      end do

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      volume = a(1,1)*b(1,1)
     >       + a(2,1)*b(2,1)
     >       + a(3,1)*b(3,1)
      
      volume = 1.0d0/volume
      call dscal(9,volume,b,1)


      do i=1,nion
         tion(1) = ion(1,i)
         tion(2) = ion(2,i)
         tion(3) = ion(3,i)
         ion(1,i) = b(1,1)*tion(1)
     >            + b(2,1)*tion(2)
     >            + b(3,1)*tion(3)

         ion(2,i) = b(1,2)*tion(1)
     >            + b(2,2)*tion(2)
     >            + b(3,2)*tion(3)

         ion(3,i) = b(1,3)*tion(1)
     >            + b(2,3)*tion(2)
     >            + b(3,3)*tion(3)

      end do

      return
      end 


c $Id$
