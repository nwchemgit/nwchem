c     $Id$                       

*     ********************************
*     *		  	  	     *
*     *	     cgsd_oep_eigenvalues    *
*     *				     *
*     ********************************

      subroutine cgsd_oep_eigenvalues()
      implicit none

#include "stdio.fh"


c#include "global.fh"
c#include "bafdecls.fh"
c#include "btdb.fh"

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
      logical     control_oep
      integer     psi_ne_excited,control_excited_ne,psi_ne
      real*8      psi_virtual,psi_eigenvalue
      external    epsi_initialize,epsi_finalize,epsi_filefind
      external    control_oep
      external    psi_ne_excited,control_excited_ne,psi_ne
      external    psi_virtual,psi_eigenvalue
      


      call Parallel_taskid(taskid)

      if (control_oep()) then

      if (taskid.eq.MASTER) write(*,1301)

      value = epsi_initialize() !*** load in virtual orbitals ***

      call psi_gen_density_potentials(1) !**** generate oep ****
      call psi_1add_oep_to_vall()

      call psi_1to2() !** store psi1 orbitals temporarilly in psi2 ***
      call psi_2gen_hml()
      call psi_diagonalize_hml_assending()
      call psi_2rotate1()


      !*** minimize orbitals ****
      call psi_minimize_f_orb()
      call psi_minimize_virtual()


      if (taskid.eq.MASTER) then

        !**** write out KS virtual eigenvalues ****
        write(*,1500)
        NN=psi_ne_excited(1)-psi_ne_excited(2)
        EV=27.2116d0
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

        !**** write out KS filled eigenvalues ****
        write(*,1501)
        NN=psi_ne(1)-psi_ne(2)
        EV=27.2116
        do I=psi_ne(1),(psi_ne(2)+1),-1
          write(*,1510) psi_eigenvalue(1,I),
     >                  psi_eigenvalue(1,I)*EV
        end do

        do I=psi_ne(2),1,-1
          write(*,1510)  psi_eigenvalue(1,I),
     >                   psi_eigenvalue(1,I)*EV,
     >                   psi_eigenvalue(2,I),
     >                   psi_eigenvalue(2,I)*EV
        end do


        write(*,*)
        write(*,*)


      end if

*     **** do not write out excited wavefunctions ***
      value = epsi_finalize(.false.)

      call psi_2to1() !**** store put psi2 back to psi1 ****

      end if

      return
 1301 FORMAT(//'== Perturbative OEP Orbital Calculation =='/)
 1500 FORMAT(/' virtual orbital energies:')
 1501 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
      end


