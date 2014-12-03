c
c     $Id$                       
c

*  ************************************************************
*  *                MPI cgminimize routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine paw_nominimize(E)
      implicit none
      real*8     E(*)


#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_energy_kin_atom.fh"
#include "paw_energy_vloc_atom.fh"
#include "paw_energy_ion_atom.fh"
#include "paw_energy_core_atom.fh"
#include "paw_energy_hartree_atom.fh"
#include "paw_xc.fh"


*     **** local variables ****
      logical value
      integer i,neall,npack1,k

      integer G1(2),R1(2)

      double precision E0,ehartree,eorbit,exc,pxc,eion
      double precision evloc_pw,evloc_atom,ehartr_pw,eke
      double precision ekin_atom
      double precision eion_atom
      double precision ecore_atom
      double precision ecore_ion_atom
      double precision ecore_self_atom
      double precision ehartree_atom
      double precision exc_atom,pxc_atom



*     **** external functions ****
      integer  paw_psi_ne,control_version
      real*8   paw_rho_error
      real*8   paw_dng_1ehartree
      real*8   paw_psi_1ke
      real*8   paw_psi_1vl
      real*8   paw_rho_1exc
      real*8   paw_rho_1pxc
      real*8   paw_psi_1eorbit
      real*8   paw_mult_energy_atom_comp !**no header file for paw_mult**
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult
   
      external paw_psi_ne,control_version
      external paw_rho_error
      external paw_dng_1ehartree
      external paw_psi_1ke
      external paw_psi_1vl
      external paw_rho_1exc
      external paw_rho_1pxc
      external paw_psi_1eorbit
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult


      call Pack_npack(1,npack1)
      neall = paw_psi_ne(1)+paw_psi_ne(2)

*     **** allocate G1, and R1 ****
      value = BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'R1',R1(2),R1(1))
      if (.not. value) 
     >  call errquit('paw_cgminimize:out of heap memory',0,MA_ERR)



*     ***** get the initial gradient and direction ****
      call paw_psi_1get_residual(dcpl_mb(R1(1)),dcpl_mb(G1(1)),E0)



*     **** plane-wave energies ****
      eorbit    = paw_psi_1eorbit()
      eke       = paw_psi_1ke()
      ehartr_pw = paw_dng_1ehartree()
      exc       = paw_rho_1exc()
      pxc       = paw_rho_1pxc()
      evloc_pw  = paw_psi_1vl()


*     **** atomic energies ***
      ekin_atom     = paw_energy_kin_atom()
      evloc_atom    = paw_energy_vloc_atom()

      ehartree_atom   = paw_energy_hartree_atom()
      eion_atom       = paw_energy_ion_atom()
      ecore_atom      = paw_energy_core_atom()
      ecore_ion_atom  = paw_energy_core_ion_atom()
      ecore_self_atom = paw_energy_core_self_atom()

      exc_atom        = paw_energy_xc_atom()


*     *** kinetic energy
      E(2) = eke
      E(3) = ekin_atom

*     *** coulomb contributions
      E(4) = ehartr_pw
      E(5) = eion_atom + ecore_atom + ehartree_atom +
     >       ecore_ion_atom + ecore_self_atom +
     >       paw_mult_energy_atom_self() +
     >       paw_mult_energy_atom_comp()

      E(6)=paw_mult_energy_atom_mult()

*     *** exch-correlation
      E(7) = exc
      E(8) = exc_atom

*     *** local pseudopot ***
      E(9)  = evloc_pw
      E(10) = evloc_atom

c      write(*,*) "eke            =",eke
c      write(*,*) "ekin_atom      =",ekin_atom
c      write(*,*)
c      write(*,*) "ehartr_pw      =",ehartr_pw
c
c      write(*,*)
c      write(*,*) "eion_atom      =",eion_atom
c      write(*,*) "ecore_atom     =",ecore_atom
c      write(*,*) "ehartree_atom  =",ehartree_atom
c      write(*,*) "ecore_ion_atom =",ecore_ion_atom
c      write(*,*) "ecore_self_atom=",ecore_self_atom
c      write(*,*) "mult_atom_self =",paw_mult_energy_atom_self()
c      write(*,*) "mult_atom_comp =",paw_mult_energy_atom_comp()
c      write(*,*)
c
c      write(*,*) "mult_mult      =",E(6)
c      write(*,*)
c      write(*,*) "exc            =",exc
c      write(*,*) "exc_atom       =",exc_atom
c      write(*,*)
c      write(*,*) "evloc_pw       =",evloc_pw
c      write(*,*) "evloc_atom     =",evloc_atom


*     *** total energy ***
      E(1) = 0.0d0
      do i=2,10
       E(1) = E(1) + E(i)
      end do

      E(11) = eorbit

      value = BA_free_heap(R1(2))
      value = value.and.BA_free_heap(G1(2))
      if (.not. value) 
     >  call errquit('paw_nominimize:error freeing heap',0, MA_ERR)

  
      return
      end
 

