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
        write(LuOut,'(" WFN1                     :")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,'(" Total WFN1 energy = ",f22.10)')wfn1_energy%E_tot
        write(LuOut,'(" 1-Electron energy = ",f22.10)')wfn1_energy%E_1el
        write(LuOut,'(" 2-Electron energy = ",f22.10)')wfn1_energy%E_2el
        write(LuOut,'(" Entropy energy    = ",f22.10)')wfn1_energy%E_ent
        write(LuOut,'(" Nuclear energy    = ",f22.10)')wfn1_energy%E_nuc
        write(LuOut,*)
        write(LuOut,'(" conventional Hartree-Fock:")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,'(" Total SCF energy  = ",f22.10)')
     &        wfn1_energy%E_tot-wfn1_energy%E_2el+wfn1_energy%E_2el_hf
        write(LuOut,'(" 2-Electron energy = ",f22.10)')
     &        wfn1_energy%E_2el_hf
      endif
c
      end subroutine wfn1_print_energy
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_error(wfn1_err)
      implicit none
#include "wfn1_error.fh"
#include "global.fh"
#include "stdio.fh"
      type(wfn1_error), intent(in) :: wfn1_err
      integer :: iproc
      integer :: nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      if (iproc.eq.0) then
        write(LuOut,'(" Tot:",f16.6," Na:",e16.6," Nb:",e16.6,
     &                              " Ca:",e16.6," Cb:",e16.6)')
     &    wfn1_err%total,wfn1_err%errna,wfn1_err%errnb,wfn1_err%errca,
     &    wfn1_err%errcb
      endif
c
      end subroutine wfn1_print_error
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_mat(dmat,nr,nc)
      implicit none
c
c     A simple serial matrix print routine
c
#include "global.fh"
#include "stdio.fh"
c
      integer, intent(in) :: nr ! the number of rows
      integer, intent(in) :: nc ! the number of columns
c
      double precision, intent(in) :: dmat(nr,nc)
c
      integer :: ir, ic, kk ! counters
      integer :: iproc, nproc
c
      integer, parameter :: mxc = 8 ! the maximum number of columns
                                    ! printed at a time
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      if (iproc.eq.0) then
c
        do ic = 1, nc, mxc
          write(LuOut,*)
          write(LuOut,'(5x,$)')
          do kk = ic, min(ic+mxc-1,nc)
            write(LuOut,'(i12,$)')kk
          enddo
          write(LuOut,*)
          write(LuOut,*)
          do ir = 1, nr
            write(LuOut,'(i4,$)')ir
            do kk = ic, min(ic+mxc-1,nc)
              write(LuOut,'(f12.6,$)')dmat(ir,kk)
            enddo
            write(LuOut,*)
          enddo
        enddo
c
      endif
c
      end subroutine wfn1_print_mat
c
c-----------------------------------------------------------------------
