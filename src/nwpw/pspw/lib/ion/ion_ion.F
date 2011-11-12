*
* $Id$
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


*     ********************************
*     *                              *
*     *         ion_ion_e_qm         *
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

      real*8 function ion_ion_e_qm()
      implicit none
      

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,e

*     *** external functions ***
      integer  ion_nion_qm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion_qm,ion_katm
      external psp_zv,ion_rion

      e = 0.0d0
      do j=1,ion_nion_qm()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          e = e + psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r
        end do
      end do

      ion_ion_e_qm = e
      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f_qm         *
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

      subroutine ion_ion_f_qm(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion_qm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion_qm,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion_qm().eq.1) return

      do j=1,ion_nion_qm()
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


*     ********************************
*     *                              *
*     *         ion_ion_e_mm         *
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

      real*8 function ion_ion_e_mm()
      implicit none
      

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,e

*     *** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_nion_qm,ion_katm
      external psp_zv,ion_rion

      e = 0.0d0
      do j=ion_nion_qm()+1,ion_nion()
        do i=ion_nion_qm()+1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          e = e + psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r
        end do
      end do

      ion_ion_e_mm = e
      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f_mm         *
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

      subroutine ion_ion_f_mm(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion_mm().eq.1) return

      do j=ion_nion_qm()+1,ion_nion()
        do i=ion_nion_qm()+1,j-1
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


*     ********************************
*     *                              *
*     *         ion_ion_f_async      *
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

      subroutine ion_ion_f_async(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j,tnp,tid,dutask
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion().eq.1) return

      call Parallel_np(tnp)
      call Parallel_taskid(tid)

      dutask = 0
      do j=1,ion_nion()
        do i=1,j-1
         if (dutask.eq.tid) then
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
         end if
         dutask = mod((dutask+1),tnp)
        end do
      end do

      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f_qm_async   *
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

      subroutine ion_ion_f_qm_async(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j,tnp,tid,dutask
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion_qm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion_qm,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion_qm().eq.1) return

      call Parallel_np(tnp)
      call Parallel_taskid(tid)

      dutask=0
      do j=1,ion_nion_qm()
        do i=1,j-1
         if (dutask.eq.tid) then
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
         end if
         dutask = mod((dutask+1),tnp)
        end do
      end do

      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f_mm_async   *
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

      subroutine ion_ion_f_mm_async(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j,tnp,tid,dutask
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion_mm().eq.1) return

      call Parallel_np(tnp)
      call Parallel_taskid(tid)

      dutask = 0
      do j=ion_nion_qm()+1,ion_nion()
        do i=ion_nion_qm()+1,j-1
         if (dutask.eq.tid) then
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
         end if
         dutask = mod((dutask+1),tnp)
        end do
      end do

      return
      end



*     ********************************
*     *                              *
*     *         ion_ion_efg          *
*     *                              *
*     ********************************
*
*   This routine computes the (free-space) electric field tensor
*   ion cores.

      subroutine ion_ion_efg(efg)
      implicit none
      real*8 efg(3,3,*)

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,v,v2

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion().eq.1) return

      do j=1,ion_nion()
        do i=1,ion_nion()
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)

          v  = psp_zv(ion_katm(i))/r**3
          v2 = psp_zv(ion_katm(i))/r**5

          efg(1,1,j) = v - 3.0d0*x*x*v2
          efg(2,1,j) = -3.0d0*x*y*v2
          efg(3,1,j) = -3.0d0*x*z*v2

          efg(1,2,j) = -3.0d0*x*y*v2
          efg(2,2,j) = v - 3.0d0*y*y*v2
          efg(3,2,j) = -3.0d0*z*y*v2

          efg(1,3,j) = -3.0d0*x*z*v2
          efg(2,3,j) = -3.0d0*y*z*v2
          efg(3,3,j) = v - 3.0d0*z*z*v2
        end do
      end do

      return
      end



