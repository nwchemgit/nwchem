*
* $Id: nose-hoover.F,v 1.5 2004-05-05 19:33:51 bylaska Exp $
*

*     ***************************
*     *				*
*     *	     Nose_Init          *
*     *				*
*     ***************************

      subroutine Nose_Init(eke0)
      implicit none
      real*8 eke0


#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "nose-hoover.fh"


*     **** boltzman constant ****
      double precision kb
      parameter (kb=3.16679d-6)


*     **** local variables ****
      logical value
      integer i,rtdb
      real*8 pi,am,fmass,tmpe
      real*8 Te
      real*8 Tr

*     **** external functions ****
      real*8   control_fake_mass
      real*8   control_Nose_Pe
      real*8   control_Nose_Te
      real*8   control_Nose_Pr
      real*8   control_Nose_Tr
      real*8   ion_amass
      integer  ion_nion,control_rtdb
      external control_fake_mass
      external control_Nose_Pe
      external control_Nose_Te
      external control_Nose_Pr
      external control_Nose_Tr
      external ion_amass
      external ion_nion,control_rtdb

*     ************************************
*     **** initialize the thermostats ****
*     ************************************
      eke0_init = eke0
      Pe = control_Nose_Pe()
      Te = control_Nose_Te()
      Pr = control_Nose_Pr()
      Tr = control_Nose_Tr()
      fmass = control_fake_mass()

      Xe0 = 0.0d0
      Xe1 = 0.0d0
      Xe2 = 0.0d0
      Xr0 = 0.0d0
      Xr1 = 0.0d0
      Xr2 = 0.0d0

*     **** Set Er0 ****
      if (ion_nion().gt.2) then
        Er0 = 0.5d0*(3.0d0*dble(ion_nion())-6.0d0)*kb*Tr
      else
        Er0 = 0.5d0*(1)*kb*Tr
      end if

      am = 0.0d0
      do i=1,ion_nion()
         am = am + ion_amass(i)
      end do

*     **** Set Ee0 - read from rtdb otherwise use current KE ****
      rtdb  = control_rtdb()
      value = rtdb_get(rtdb,'cpmd:eke0',mt_dbl,1,tmpe)
      write(*,*) "value:",value,tmpe
      if (value) then
        eke0_init = tmpe
      end if
      Ee0 = 4.0*kb*Te*fmass*dble(ion_nion())/am * eke0_init

*     **** Set Qe and Qr - read from rtdb otherwise set using periods ****
      pi = 4.0d0*datan(1.0d0)
      value = rtdb_get(rtdb,'cpmd:Qe',  mt_dbl,1,Qe)
     >   .and.rtdb_get(rtdb,'cpmd:Qr',  mt_dbl,1,Qr)
      write(*,*) "value:",value,Qe,Qr
      if (.not.value) then
        Qe = Pe*Pe*Ee0/(pi*pi)
        Qr = Pr*Pr*Er0/(pi*pi)
      else
        Pe = pi*dsqrt(Qe/Ee0)
        Pr = pi*dsqrt(Qr/Er0)
      end if

      return
      end


*     ***************************
*     *                         *
*     *      Nose_end           *
*     *                         *
*     ***************************

      subroutine Nose_end()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "nose-hoover.fh"


      !**** local variables ****
      logical value
      integer rtdb

      !**** external functions ****
      integer  control_rtdb
      external control_rtdb

      rtdb  = control_rtdb()
      value = rtdb_put(rtdb,'cpmd:eke0',mt_dbl,1,eke0_init)
     >   .and.rtdb_put(rtdb,'cpmd:Qe',  mt_dbl,1,Qe)
     >   .and.rtdb_put(rtdb,'cpmd:Qr',  mt_dbl,1,Qr)
      if (.not.value)
     >  call errquit(
     >       'Nose_End: error writing to rtdb', 0, RTDB_ERR)

      return
      end



*     ***************************
*     *				*
*     *	     Nose_reset_T       *
*     *				*
*     ***************************

      subroutine Nose_reset_T(Te_new,Tr_new)
      implicit none
      real*8 Te_new,Tr_new
     
#include "nose-hoover.fh"

*     **** boltzman constant ****
      double precision kb
      parameter (kb=3.16679d-6)

*     **** local variables ****
      integer i
      real*8 am,fmass

*     **** external functions ****
      real*8   control_fake_mass
      real*8   ion_amass
      integer  ion_nion
      external control_fake_mass
      external ion_amass
      external ion_nion

      fmass = control_fake_mass()

*     **** reSet Er0 ****
      if (ion_nion().gt.2) then
        Er0 = 0.5d0*(3.0d0*dble(ion_nion())-6.0d0)*kb*Tr_new
      else
        Er0 = 0.5d0*(1)*kb*Tr_new
      end if

      am = 0.0d0
      do i=1,ion_nion()
         am = am + ion_amass(i)
      end do

*     **** reSet Ee0 ****
      Ee0 = 4.0*kb*Te_new*fmass*dble(ion_nion())/am * eke0_init

      return
      end



*     ***************************
*     *				*
*     *	     Nose_Newton_Step   *
*     *				*
*     ***************************

      subroutine Nose_Newton_Step(eke,eki)
      implicit none
      real*8 eke,eki

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 FXe,FXr,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      FXe = 2.0d0*(eke-Ee0)
      FXr = 2.0d0*(eki-Er0)

      Xe2 = (0.5d0*dt*dt/Qe)*FXe
      Xr2 = (0.5d0*dt*dt/Qr)*FXr


      return
      end

*     ***************************
*     *                         *
*     *      Nose_Verlet_Step   *
*     *                         *
*     ***************************

      subroutine Nose_Verlet_Step(eke,eki)
      implicit none
      real*8 eke,eki

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 FXe,FXr,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      FXe = 2.0d0*(eke-Ee0)
      FXr = 2.0d0*(eki-Er0)

      Xe2 = 2.0d0*Xe1 - Xe0 + (dt*dt/Qe)*FXe
      Xr2 = 2.0d0*Xr1 - Xr0 + (dt*dt/Qr)*FXr

      return
      end

*     ***************************
*     *                         *
*     *      Nose_dXe           *
*     *                         *
*     ***************************

      real*8 function Nose_dXe()
      implicit none

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()
      Nose_dXe = (3.0d0*Xe1-4.0d0*Xe0+Xem)/(2.0d0*dt)
      return
      end

*     ***************************
*     *                         *
*     *      Nose_dXr           *
*     *                         *
*     ***************************

      real*8 function Nose_dXr()
      implicit none

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()
      Nose_dXr = (3.0d0*Xr1-4.0d0*Xr0+Xrm)/(2.0d0*dt)
      return
      end

*     ***************************
*     *                         *
*     *      Nose_sse           *
*     *                         *
*     ***************************

      real*8 function Nose_sse()
      implicit none

#include "nose-hoover.fh"


*     ***** local variables ****
      real*8 dXe,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      dXe         = (3.0d0*Xe1-4.0d0*Xe0+Xem)/(2.0d0*dt)
      Nose_sse = 1.0d0/(1.0d0+0.5d0*dXe*dt)
      return
      end

*     ***************************
*     *                         *
*     *      Nose_ssr           *
*     *                         *
*     ***************************

      real*8 function Nose_ssr()
      implicit none

#include "nose-hoover.fh"

*     ***** local variables ****
      real*8 dXr,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      dXr         = (3.0d0*Xr1-4.0d0*Xr0+Xrm)/(2.0d0*dt)
      Nose_ssr = 1.0d0/(1.0d0+0.5d0*dXr*dt)
      return
      end

*     ***************************
*     *                         *
*     *      Nose_shift         *
*     *                         *
*     ***************************

      subroutine Nose_shift()
      implicit none

#include "nose-hoover.fh"

      Xem = Xe0
      Xe0 = Xe1
      Xe1 = Xe2
      Xrm = Xr0
      Xr0 = Xr1
      Xr1 = Xr2
      return
      end

*     ***************************
*     *                         *
*     *      Nose_e_energy      *
*     *                         *
*     ***************************

      real*8 function Nose_e_energy()
      implicit none

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 dXe,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      dXe          = (3.0d0*Xe1-4.0d0*Xe0+Xem)/(2.0d0*dt)
      Nose_e_energy = 0.5d0*Qe*dXe**2 + 2.0d0*Ee0*Xe1
      return
      end

*     ***************************
*     *                         *
*     *      Nose_r_energy      *
*     *                         *
*     ***************************

      real*8 function Nose_r_energy()
      implicit none

#include "nose-hoover.fh"

*     **** local variables ****
      real*8 dXr,dt

*     **** external functions ****
      real*8   control_time_step
      external control_time_step

      dt = control_time_step()

      dXr          = (3.0d0*Xr1-4.0d0*Xr0+Xrm)/(2.0d0*dt)
      Nose_r_energy = 0.5d0*Qr*dXr**2 + 2.0d0*Er0*Xr1
      return
      end


*     *********************
*     *                   *
*     *      Nose_Qe      *
*     *                   *
*     *********************

      real*8 function Nose_Qe()
      implicit none

#include "nose-hoover.fh"

      Nose_Qe = Qe
      return
      end


*     *********************
*     *                   *
*     *      Nose_Pe      *
*     *                   *
*     *********************

      real*8 function Nose_Pe()
      implicit none

#include "nose-hoover.fh"

      Nose_Pe = Pe
      return
      end


*     *********************
*     *                   *
*     *      Nose_Qr      *
*     *                   *
*     *********************

      real*8 function Nose_Qr()
      implicit none

#include "nose-hoover.fh"

      Nose_Qr = Qr
      return
      end

*     *********************
*     *                   *
*     *      Nose_Pr      *
*     *                   *
*     *********************

      real*8 function Nose_Pr()
      implicit none

#include "nose-hoover.fh"

      Nose_Pr = Pr
      return
      end


*     *********************
*     *                   *
*     *      Nose_Ee0     *
*     *                   *
*     *********************

      real*8 function Nose_Ee0()
      implicit none

#include "nose-hoover.fh"

      Nose_Ee0 = Ee0
      return
      end

*     *********************
*     *                   *
*     *      Nose_Er0     *
*     *                   *
*     *********************

      real*8 function Nose_Er0()
      implicit none

#include "nose-hoover.fh"

      Nose_Er0 = Er0
      return
      end


