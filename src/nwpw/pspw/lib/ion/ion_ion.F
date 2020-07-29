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
      integer i,j,nion,nion_q
      real*8 x,y,z,r,qi,qj
      integer qcoord_ptr

      !common /eiond_ejtmp/ e

      real*8 e1,e2
      common /pack_conjg_tcc_indot_sasb/ e1,e2

*     *** external functions ***
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external psp_zv,ion_rion
      real*8   nwpw_cosmo_screen_qc0,nwpw_cosmo_Eqq,nwpw_cosmo_EQionq
      external nwpw_cosmo_screen_qc0,nwpw_cosmo_Eqq,nwpw_cosmo_EQionq

      nion   = ion_nion()
      !nion_q = ion_nion_q()

!$OMP MASTER
      e1 = 0.0d0
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:e1)
      do j=1,nion
        qj = psp_zv(ion_katm(j))
        do i=1,j-1
          qi = psp_zv(ion_katm(i))
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          if (r.gt.1.0d-6) then
             e1 = e1 + qi*qj/r
          end if
        end do
      end do
!$OMP END DO

c!$OMP DO reduction(+:e2)
c      do i=nion+1,nion+nion_q
c         qi = nwpw_cosmo_screen_qc0(i-nion)
c         do j=1,nion
c           qj = psp_zv(ion_katm(j))
c           x = ion_rion(1,i) - ion_rion(1,j)
c           y = ion_rion(2,i) - ion_rion(2,j)
c           z = ion_rion(3,i) - ion_rion(3,j)
c           r = dsqrt(x*x + y*y + z*z)
c           if (r.gt.1.0d-6) then
c              e2 = e2 + qi*qj/r
c           end if
c        end do
c      end do
c!$OMP END DO

      ion_ion_e = e1 + nwpw_cosmo_EQionq() + nwpw_cosmo_Eqq()
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

!$OMP MASTER
      do j=1,ion_nion()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)

          if (r.gt.1.0d-6) then
             v = psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r**3

             f_ion(1,i) = f_ion(1,i) + (x*v)
             f_ion(2,i) = f_ion(2,i) + (y*v)
             f_ion(3,i) = f_ion(3,i) + (z*v)

             f_ion(1,j) = f_ion(1,j) - (x*v)
             f_ion(2,j) = f_ion(2,j) - (y*v)
             f_ion(3,j) = f_ion(3,j) - (z*v)
          end if
        end do
      end do
!$OMP END MASTER
!$OMP BARRIER

      call nwpw_cosmo_FQionq(f_ion)

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

!$OMP MASTER
      dutask = 0
      do j=1,ion_nion()
        do i=1,j-1
         if (dutask.eq.tid) then
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          if (r.gt.1.0d-6) then
             v = psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r**3

             f_ion(1,i) = f_ion(1,i) + (x*v)
             f_ion(2,i) = f_ion(2,i) + (y*v)
             f_ion(3,i) = f_ion(3,i) + (z*v)

             f_ion(1,j) = f_ion(1,j) - (x*v)
             f_ion(2,j) = f_ion(2,j) - (y*v)
             f_ion(3,j) = f_ion(3,j) - (z*v)
          end if
         end if
         dutask = mod((dutask+1),tnp)
        end do
      end do
!$OMP END MASTER
!$OMP BARRIER

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
!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER

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
!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER

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
          if (r.gt.1.0d-6) then
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
          end if
        end do
      end do

      return
      end



