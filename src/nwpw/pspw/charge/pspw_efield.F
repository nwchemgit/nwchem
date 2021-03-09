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
 
      subroutine pspw_Efield_init()
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

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

*     **** local variables ****
      logical value

      return
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
*     *     pspw_Efield_Energy_ion     *
*     *                                *
*     **********************************

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

      energy = 0.0d0
      do ii=1,ion_nion()
         ia = ion_katm(ii)
         qii = psp_zv(ia)
         energy = -qii*efield(1)*(ion_rion(1,ii) - efield_center(1))
     >           - qii*efield(2)*(ion_rion(2,ii) - efield_center(2))
     >           - qii*efield(3)*(ion_rion(3,ii) - efield_center(3))
      end do

      pspw_Efield_Energy_ion = energy
      return
      end



*     **********************************
*     *                                *
*     *     pspw_Efield_Generate_V     *
*     *                                *
*     **********************************

      subroutine pspw_charge_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "bafdecls.fh"
#include "pspw_efield.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer k

!$OMP DO
      do k=1,n2ft3d
        Vqm(k) = Vqm(k) - efield(1)*(rgrid(1,k)-efield_center(1))
     >                  - efield(2)*(rgrid(2,k)-efield_center(2))
     >                  - efield(3)*(rgrid(3,k)-efield_center(3))
      end do
!$OMP END DO

      return
      end

