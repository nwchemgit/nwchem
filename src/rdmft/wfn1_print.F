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
