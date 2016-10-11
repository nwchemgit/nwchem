c-----------------------------------------------------------------------
c
      subroutine wfn1_print_wave(wfn1_wave)
      implicit none
c
#include "wfn1_wfn.fh"
#include "global.fh"
#include "stdio.fh"
c
      type(wfn1_wfn), intent(in) :: wfn1_wave ! The wavefunction
c
      integer iproc ! the process rank
      integer nproc ! the number of process ranks
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
        write(LuOut,'(" The alpha natural orbitals")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_wave%g_na)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta natural orbitals")')
        write(LuOut,'(" -------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_wave%g_nb)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The alpha correlation functions")')
        write(LuOut,'(" -------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_wave%g_ca)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta correlation functions")')
        write(LuOut,'(" ------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_wave%g_cb)
c
      end subroutine wfn1_print_wave
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_deriv(wfn1_deriv)
      implicit none
c
#include "wfn1_drv.fh"
#include "global.fh"
#include "stdio.fh"
c
      type(wfn1_drv), intent(in) :: wfn1_deriv ! The derivative
c
      integer iproc ! the process rank
      integer nproc ! the number of process ranks
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
        write(LuOut,'(" The alpha natural orbital derivatives")')
        write(LuOut,'(" -------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dna)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta natural orbitals derivatives")')
        write(LuOut,'(" -------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dnb)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The alpha correlation function derivatives")')
        write(LuOut,'(" ------------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dca)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta correlation function derivatives")')
        write(LuOut,'(" -----------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dcb)
c
      end subroutine wfn1_print_deriv
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_energy(wfn1_energy)
      implicit none
#include "wfn1_en.fh"
#include "global.fh"
#include "stdio.fh"
      type(wfn1_en), intent(in) :: wfn1_energy
      integer iproc
      integer nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      if (iproc.eq.0) then
        write(LuOut,'(" Total energy      = ",f22.10)')wfn1_energy%E_tot
        write(LuOut,'(" 1-Electron energy = ",f22.10)')wfn1_energy%E_1el
        write(LuOut,'(" 2-Electron energy = ",f22.10)')wfn1_energy%E_2el
        write(LuOut,'(" Nuclear energy    = ",f22.10)')wfn1_energy%E_2el
      endif
c
      end subroutine wfn1_print_energy
c
c-----------------------------------------------------------------------
