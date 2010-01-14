c     $Id: cgsd_excited.F,v 1.9 2008-09-11 21:26:47 bylaska Exp $                       

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
      integer nexcited(2),ispin
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(30)
      integer rtdb

*     **** external functions ****
      logical     control_DOS,control_Mulliken
      logical     epsi_initialize,epsi_finalize
      logical     control_check_number_virtuals
      integer     psi_ne_excited,control_excited_ne,control_rtdb
      real*8      psi_virtual
      external    control_DOS,control_Mulliken
      external    epsi_initialize,epsi_finalize
      external    control_check_number_virtuals
      external    psi_ne_excited,control_excited_ne,control_rtdb
      external    psi_virtual
      integer     control_symmetry,control_mapping1d,psi_ispin
      external    control_symmetry,control_mapping1d,psi_ispin
      character*4 psi_ab_irrep_name
      external    psi_ab_irrep_name


      call Parallel_taskid(taskid)
      neall = control_excited_ne(1)+control_excited_ne(2)

      if (neall.gt.0) then

      if (taskid.eq.MASTER) write(luout,1301)

*     **** read in excited wavefunctions and initialize epsi ****
      call control_ispin_set(psi_ispin())
      if (.not.control_check_number_virtuals()) then
        call epsi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if
      call psi_get_ne_excited(ispin,nexcited)
cccccccc  call Dnex_init(ispin,nexcited,control_mapping1d()) ccccccc

      value = epsi_initialize()
      call psi_gen_density_potentials(1)
      call psi_minimize_virtual()

      if (control_symmetry().gt.0) call psi_ab_gen_irrep_names(.true.)

      if (taskid.eq.MASTER) then
        write(luout,1500)
        NN=psi_ne_excited(1)-psi_ne_excited(2)
        EV=27.2116
        if (control_symmetry().eq.1) then

          do I=psi_ne_excited(1),(psi_ne_excited(2)+1),-1
            write(luout,1512) psi_virtual(1,I),
     >                    psi_virtual(1,I)*EV,
     >                    psi_ab_irrep_name(I)
          end do

          do I=psi_ne_excited(2),1,-1
            write(luout,1512)  psi_virtual(1,I),
     >                     psi_virtual(1,I)*EV,
     >                     psi_ab_irrep_name(I),
     >                     psi_virtual(2,I),
     >                     psi_virtual(2,I)*EV,
     >                     psi_ab_irrep_name(I+psi_ne_excited(1))
          end do
        else
          do I=psi_ne_excited(1),(psi_ne_excited(2)+1),-1
            write(luout,1510) psi_virtual(1,I),
     >                    psi_virtual(1,I)*EV
          end do

          do I=psi_ne_excited(2),1,-1
            write(luout,1510)  psi_virtual(1,I),
     >                     psi_virtual(1,I)*EV,
     >                     psi_virtual(2,I),
     >                     psi_virtual(2,I)*EV
          end do
        end if
        write(luout,*)
        write(luout,*)

      end if
      if (control_symmetry().gt.0) call psi_ab_kill_irrep_names()



*     **** calculate Density of States ***
      rtdb = control_rtdb()
      if (control_DOS()) then
         call epsi_DOS(rtdb)
      end if

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) then
         call epsi_Mulliken(rtdb)
         call kp_structure(rtdb)
      end if



*     **** write out excited wavefunctions ***
      value = epsi_finalize(.true.)

ccccccc  call Dnex_end() ccccccc
      end if

      return
 1301 FORMAT(//'== Virtual Orbital Calculation =='/)
 1500 FORMAT(/' virtual orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
      end


