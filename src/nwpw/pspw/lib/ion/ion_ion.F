*
* $Id: ion_ion.F,v 1.1 2003-04-07 22:05:51 windus Exp $
*


*     ********************************
*     *                              *
*     *         ion_ion_e            *
*     *                              *
*     ********************************
*
*   This function computes the (free-space) Coulomb energy between
*   ion cores.
*
*   Entry -
*       nion  ---- number of ions
*       katm  ---- element id number
*       r_ion ---- coordinates of ions
*

      real*8 function ion_ion_e()
      implicit none
      

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,e

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

      e = 0.0d0
      do j=1,ion_nion()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          e = e + psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r
        end do
      end do

      ion_ion_e = e
      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f            *
*     *                              *
*     ********************************
*
*   This routine computes the (free-space) Coulomb forces between
*   ion cores.
*
*   Entry -
*       nion  ---- number of ions
*       katm  ---- element id number
*       r_ion ---- coordinates of ions
*
*   Exit -
*       f_ion ---- force vectors
*

      subroutine ion_ion_f(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion().eq.1) return

      do j=1,ion_nion()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)

          v = psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r**3

          f_ion(1,i) = f_ion(1,i) + (x*v)
          f_ion(2,i) = f_ion(2,i) + (y*v)
          f_ion(3,i) = f_ion(3,i) + (z*v)

          f_ion(1,j) = f_ion(1,j) - (x*v)
          f_ion(2,j) = f_ion(2,j) - (y*v)
          f_ion(3,j) = f_ion(3,j) - (z*v)
        end do
      end do

      return
      end
