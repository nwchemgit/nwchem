c     $Id: cgsd_excited.F,v 1.2 2003-07-19 20:47:55 bylaska Exp $                       

*     ***************************
*     *				*
*     *		cgsd_excited    *
*     *				*
*     ***************************

      subroutine cgsd_excited()
      implicit none

#include "stdio.fh"


c#include "global.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"

      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled,value
      integer taskid
      integer minimizer
      integer i,j,ms,neall,NN
      integer nexcited(2)
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(30)
      integer rtdb

*     **** external functions ****
      logical     epsi_initialize,epsi_finalize,epsi_filefind
      integer     psi_ne_excited,control_excited_ne
      real*8      psi_virtual
      external    epsi_initialize,epsi_finalize,epsi_filefind
      external    psi_ne_excited,control_excited_ne
      external    psi_virtual
      


      call Parallel_taskid(taskid)
      neall = control_excited_ne(1)+control_excited_ne(2)

      if (neall.gt.0) then

      if (taskid.eq.MASTER) write(*,1301)

*     **** read in excited wavefunctions and initialize epsi ****
      if (.not.epsi_filefind()) then
        call epsi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if
      value = epsi_initialize()

      call psi_minimize_virtual()


      if (taskid.eq.MASTER) then
        write(*,1500)
        NN=psi_ne_excited(1)-psi_ne_excited(2)
        EV=27.2116
        do I=psi_ne_excited(1),(psi_ne_excited(2)+1),-1
          write(*,1510) psi_virtual(1,I),
     >                  psi_virtual(1,I)*EV
        end do

        do I=psi_ne_excited(2),1,-1
          write(*,1510)  psi_virtual(1,I),
     >                   psi_virtual(1,I)*EV,
     >                   psi_virtual(2,I),
     >                   psi_virtual(2,I)*EV
        end do
        write(*,*)
        write(*,*)

      end if

*     **** write out excited wavefunctions ***
      value = epsi_finalize(.true.)

      end if

      return
 1301 FORMAT(//'== Virtual Orbital Calculation =='/)
 1500 FORMAT(/' virtual orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
      end


