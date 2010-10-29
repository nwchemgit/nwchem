
*
* $Id$
*

*     ***************************
*     *				*
*     *   Ursenbach_Smoother	*
*     *				*
*     ***************************
*
*   This routine implements the density smoothing 
*  routine implemented by 
*
*         Ursenbach et al, J.Chem.Phys. vol 106, No. 7, page 2811.
*
*
      subroutine Ursenbach_Smoother(ngrid,f,fs,fsp)
      implicit none
      integer ngrid
      real*8 f(*)
      real*8 fs(*)
      real*8 fsp(*)

*     **** parameters ****
      real*8 tolerance
      parameter (tolerance=0.005d0)
      real*8 a1,a2,a3
      parameter (a1=40000.0d0)
      parameter (a2=(1.0d0/401.0d0) )
      parameter (a3=80000.0d0)

*     **** local variables ****
      integer i

      do i=1,ngrid
         if (f(i) .ge. tolerance) then
            fs(i)  = f(i)
            fsp(i) = 1.0d0

         else if (f(i) .lt. 0.0d0) then
            fs(i)  = a2*exp(f(i))
            fsp(i) = fs(i)

         else
            fs(i)  = a2*(a1*f(i)**2 + f(i) + 1.0d0) 
            fsp(i) = a2*(a3*f(i) + 1.0d0)
         end if


      end do

      return
      end
