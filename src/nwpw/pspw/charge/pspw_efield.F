c
c $Id$
c
  
***********************************************************
*                                                         *
*                   PSPW-Efield module                    *
*                                                         *
***********************************************************

*     **********************************
*     *	                               *
*     *          pspw_Efield_init      *
*     *                                *
*     **********************************
 
      subroutine pspw_Efield_init(rtdb,ispin,ne)
      implicit none
      integer rtdb
      integer ispin,ne(2)

#include "bafdecls.fh"
#include "btdb.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

      if (.not.btdb_get(rtdb,'nwpw:efield',mt_log,1,efield_on))
     >   efield_on = .false.

      if (efield_on) then
         if(.not.btdb_get(rtdb,'nwpw:efield_type',mt_int,1,efield_type))
     >      efield_type = 2

         if(.not.btdb_get(rtdb,'nwpw:efield_vector',
     >                    mt_dbl,3,efield)) then
            efield(1) = 0.0d0
            efield(2) = 0.0d0
            efield(3) = 0.0d0
         end if

         if(.not.btdb_get(rtdb,'nwpw:efield_center',
     >                    mt_dbl,3,efield_center)) then
            efield_center(1) = 0.0d0
            efield_center(2) = 0.0d0
            efield_center(3) = 0.0d0
         end if
         if (efield_type.eq.0) 
     >      call dipole_Efield_init(ispin,ne,efield,efield_center)

      end if

      return
      end



*     **********************************
*     *	                               *
*     *        pspw_Efield_end         *
*     *                                *
*     **********************************
      subroutine pspw_Efield_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

      if (efield_on) then
         if (efield_type.eq.0)  call dipole_Efield_end()
      end  if

      return
      end


*     **********************************
*     *                                *
*     *       pspw_Efield_Print        *
*     *                                *
*     **********************************

      subroutine pspw_Efield_Print(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      if (efield_on) then
         call Parallel_taskid(taskid)

         if (taskid.eq.MASTER) then
            if (efield_type.eq.0) write(unit,1160)
            if (efield_type.eq.1) write(unit,1170)
            if (efield_type.eq.2) write(unit,1180)
            write(unit,1190) efield
            write(unit,1191) efield_center
         end if
      end if

      return
 1160 format(/' periodic Electric field:')
 1170 format(/' APC Electric field:')
 1180 format(/' real space Electric field:')
 1190 format(5X, 'Electric Field (au) (',3F11.5,' )')
 1191 format(5X, '        Center (au) (',3F11.5,' )')
      end




*     **********************************
*     *	                               *
*     *         pspw_Efield_found      *
*     *                                *
*     **********************************
      logical function pspw_Efield_found()
      implicit none

#include "pspw_efield.fh"

      pspw_Efield_found = efield_on
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_Efield_type       *
*     *                                *
*     **********************************
      integer function pspw_Efield_type()
      implicit none

#include "pspw_efield.fh"

      pspw_Efield_type = efield_type
      return
      end


*     **********************************
*     *                                *
*     *    pspw_Efield_Energy_dipole   *
*     *                                *
*     **********************************
*    This routine calculates the energybetween the QM ions and efield
*  Note the ions are positive.
*
      real*8 function pspw_Efield_Energy_dipole(dipole)
      implicit none
      real*8 dipole(3)

#include "bafdecls.fh"
#include "pspw_efield.fh"

*     **** local variables ****
      integer ii,ia
      real*8  qii
      real*8  r,energy

*     **** external functions ****
      integer   ion_nion,ion_katm
      real*8    ion_rion,psp_zv
      external  ion_nion,ion_katm
      external  ion_rion,psp_zv

      energy = dipole(1)*efield(1)
     >       + dipole(2)*efield(2)
     >       + dipole(3)*efield(3)

      pspw_Efield_Energy_dipole = energy
      return
      end


*     **********************************
*     *                                *
*     *     pspw_Efield_Energy_ion     *
*     *                                *
*     **********************************
*    This routine calculates the energybetween the QM ions and efield
*  Note the ions are positive.
*
      real*8 function pspw_Efield_Energy_ion()
      implicit none

#include "bafdecls.fh"
#include "pspw_efield.fh"

*     **** local variables ****
      integer ii,ia
      real*8  qii
      real*8  r,energy

*     **** external functions ****
      integer   ion_nion,ion_katm
      real*8    ion_rion,psp_zv
      external  ion_nion,ion_katm
      external  ion_rion,psp_zv

!$OMP MASTER
      eetmp = 0.0d0
!$OMP END MASTER
!$OMP BARRIER

      if (efield_on) then
!$OMP DO reduction(+:eetmp)
         do ii=1,ion_nion()
            ia = ion_katm(ii)
            qii = psp_zv(ia)
            eetmp = eetmp
     >              + qii*efield(1)*(ion_rion(1,ii) - efield_center(1))
     >              + qii*efield(2)*(ion_rion(2,ii) - efield_center(2))
     >              + qii*efield(3)*(ion_rion(3,ii) - efield_center(3))
         end do
!$OMP END DO
      end if

      pspw_Efield_Energy_ion = eetmp
      return
      end

*     **********************************
*     *                                *
*     *         pspw_Efield_Fion       *
*     *                                *
*     **********************************
*    This routine calculates the forces between the QM ions and efield
*  Note the ions are positive.
*
*     Exit - fion: force on the QM ions
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_Efield_Fion(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "pspw_efield.fh"

*     **** local variables ****
      integer i
      real*8  qii

*     **** external functions ****
      integer  ion_katm,ion_nion
      external ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external psp_zv,ion_rion

      if (efield_on.and.(efield_type.eq.2)) then
!$OMP DO
         do i=1,ion_nion()
            qii = psp_zv(ion_katm(i))
            fion(1,i) = fion(1,i) - qii*efield(1)
            fion(2,i) = fion(2,i) - qii*efield(2)
            fion(3,i) = fion(3,i) - qii*efield(3)
         end do
!$OMP END DO
      end if

      return
      end




*     **********************************
*     *                                *
*     *     pspw_Efield_Generate_V     *
*     *                                *
*     **********************************
*
*   Applying Efield on to electron density.  
* Note the minus charge is because the charge is negative.
*
      subroutine pspw_Efield_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer k

      if (efield_on.and.(efield_type.eq.2)) then
!$OMP DO
         do k=1,n2ft3d
            Vqm(k) = Vqm(k) - efield(1)*(rgrid(1,k)-efield_center(1))
     >                      - efield(2)*(rgrid(2,k)-efield_center(2))
     >                      - efield(3)*(rgrid(3,k)-efield_center(3))
         end do
!$OMP END DO
      end if

      return
      end

