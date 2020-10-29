c     $Id$                       

*     ***************************
*     *				*
*     *	      cgsd_excited      *
*     *				*
*     ***************************

      subroutine cgsd_excited()
      implicit none

#include "stdio.fh"

c#include "global.fh"
c#include "bafdecls.fh"
#include "btdb.fh"

      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled,value,gvirtual,bvirtual,confine,oprint,gci
      integer taskid
      integer minimizer,gfac
      integer i,j,ms,neall,NN
      integer nexcited(2),ispin
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(30)
      integer rtdb,civirt

*     **** external functions ****
      logical     control_DOS,control_Mulliken,control_MP2
      logical     epsi_initialize,epsi_finalize,control_2qintegrals
      logical     control_check_number_virtuals
      integer     psi_ne_excited,control_excited_ne,control_rtdb
      real*8      psi_virtual
      external    control_DOS,control_Mulliken,control_MP2
      external    epsi_initialize,epsi_finalize,control_2qintegrals
      external    control_check_number_virtuals
      external    psi_ne_excited,control_excited_ne,control_rtdb
      external    psi_virtual
      integer     control_symmetry,control_mapping1d,psi_ispin
      external    control_symmetry,control_mapping1d,psi_ispin
      character*4 psi_ab_irrep_name
      external    psi_ab_irrep_name
      logical     band_reformat_c_ewvfnc
      external    band_reformat_c_ewvfnc
      integer     epsi_get_version
      external    epsi_get_version
      logical     control_gradient_virtual,control_bound_virtual
      external    control_gradient_virtual,control_bound_virtual
      logical     control_confine_virtual
      external    control_confine_virtual
      logical     psi_2epsi_gradients
      external    psi_2epsi_gradients
      integer     control_gradient_virtual_fac
      external    control_gradient_virtual_fac
      integer     control_CI_virtual
      external    control_CI_virtual

      call Parallel_taskid(taskid)

      gvirtual = control_gradient_virtual()
      gfac     = control_gradient_virtual_fac()
      bvirtual = control_bound_virtual()
      civirt   = control_CI_virtual()
      gci      = (civirt.gt.0)
      confine  = control_confine_virtual()
      neall = control_excited_ne(1)+control_excited_ne(2)

      oprint = (taskid.eq.MASTER)

      if ((neall.gt.0).or.(gvirtual)) then

*        **** Using gradients of filled orbitals to define virtual space ****
         if (gvirtual) then
            value = psi_2epsi_gradients(gfac)

*        **** calculating regular virtual orbitals ****
         else
            if (oprint.and.(.not.bvirtual)) write(luout,1301)
            if (oprint.and.(bvirtual))      write(luout,1302)

*           **** read in excited wavefunctions and initialize epsi ****
            call control_ispin_set(psi_ispin())
            if (.not.control_check_number_virtuals()) then
              call epsi_new()
              newpsi = .true.
            else
              newpsi = .false.

*              **** convert from band format to pspw format ****
               if (epsi_get_version().eq.5) then
                 newpsi = .true.
                 value = btdb_parallel(.false.)
                 if (taskid.eq.MASTER) then
                   value= band_reformat_c_ewvfnc(1)
                 end if
                 value = btdb_parallel(.true.)
               end if

            end if
         end if

      call psi_get_ne_excited(ispin,nexcited)
ccccc      call Dnex_init(ispin,nexcited,control_mapping1d())

      value = epsi_initialize()

      if (confine)   call electron_create_confine()

      if (bvirtual) then
         call pspw_HFX_virtual_set(.false.)
         call psi_gen_density_potentials(3)
      else
         call psi_gen_density_potentials(1)
      end if

      if (gvirtual) then
         call psi_gen_hml_virtual(.true.)
      else if (gci) then
         !call psi_gen_hml_virtual(.true.)
         call psi_minimize_virtual_CI(civirt)
      else
         call psi_gen_hml_virtual(.true.)
         call psi_minimize_virtual()
      end if

      if (confine) call electron_destroy_confine()
      if (bvirtual) call pspw_HFX_virtual_set(.true.)
      if ((bvirtual).or.(confine.and.(.not.gvirtual))) then
         call psi_gen_density_potentials(1)
         call psi_gen_hml_virtual(.true.)
      end if
         

      if (control_symmetry().gt.0) call psi_ab_gen_irrep_names(.true.)

      if (oprint) then
        write(luout,1500)
        NN=psi_ne_excited(1)-psi_ne_excited(2)
        EV=27.2116d0
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

*     **** calculate MP2 Energy ***
      if (control_MP2()) then
         call psi_MP2_energy(rtdb)
      end if

      if (control_2qintegrals()) then
         call psi_2q_integrals(rtdb)
      end if


*     **** write out excited wavefunctions ***
      value = epsi_finalize(.true.)

cccc      call Dnex_end()
      end if

      return
 1301 FORMAT(//'== Virtual Orbital Calculation =='/)
 1302 FORMAT(//'== Bound Virtual Orbital Calculation =='/)
 1500 FORMAT(/' virtual orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
      end


