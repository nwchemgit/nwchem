      subroutine xtb_print_param(calc, param)
      use tblite_xtb_calculator, only: xtb_calculator
      use tblite_scf_info, only : scf_info,atom_resolved,shell_resolved
      use tblite_param, only : param_record
      use tblite_coulomb_charge, only : coulomb_kernel

      implicit none

#include "stdio.fh"
#include "util.fh"
#include "util_params.fh"
#include "xtb.fh"

      type(xtb_calculator)           :: calc
      type(param_record)             :: param
      type(scf_info),allocatable     :: info

      character(len=255) temp, kernel

      integer i

      allocate(info)
      info = calc%coulomb%es3%variable_info()

      select case(info%charge)
      case (atom_resolved)
        temp = 'atom-resolved'
      case (shell_resolved)
        temp = 'shell-resolved'
      end select

      select case(param%charge%kernel)
      case(coulomb_kernel%effective)
        kernel = 'Klopman-Ohno'
      case(coulomb_kernel%dftb_gamma)
        kernel = 'DFTB gamma'
      end select


 5000 format(/,10X,'No. of atomic orbitals: ',I10,
     $       /,10X,'No. of shells         : ',I10,
     $       /,10X,'Electronic temp.      : ',F10.6,
     $       /,10X,'Accuracy              : ',F10.6,
     $       /,10X,'    Integral cutoff   : ',F10.6,
     $       /,10X,'    E. convergence    : ',F10.6,
     $       /,10X,'    W.F. convergence  : ',F10.6,
     $       /,10X,'Broyden damping       : ',F10.6,
     $       /,10X,'Coulomb kernel        : ',A,
     $       /,10X,'Coulomb 3rd order     : ',A,
     $       /,10X,'H0-scaling            : ',3(F10.6),
     $       /,10X,'s8                    : ',F10.6,
     $       /,10X,'a1                    : ',F10.6,
     $       /,10X,'a2                    : ',F10.6,
     $       /,10X,'s9                    : ',F10.6,
     $       /,10X,'kExp                  : ',2(F10.6),
     $       /,10X,'alpha                 : ',F10.6,
     $       /,10X,'a3                    : ',F10.6,
     $       /,10X,'a5                    : ',F10.6,
     $       /,10X,'cn-shift              : ',F10.6,
     $       /,10X,'cn-exp                : ',F10.6,
     $       /,10X,'max-rad               : ',F10.6)


      write(luout,5000) calc%bas%nao,calc%bas%nsh,kt/kb_au,acc,
     $      calc%bas%intcut, 1d-6*acc, 2d-5*acc,calc%mixer_damping,
     $      kernel(1:len_trim(kernel)),temp(1:len_trim(temp)),
     $      (param%hamiltonian%ksh(i,i),i=0,param%hamiltonian%lmax),
     $      param%dispersion%s8,param%dispersion%a1,param%dispersion%a2,
     $      param%dispersion%s9,param%repulsion%kexp,
     $      param%repulsion%klight,param%charge%gexp,
     $      param%multipole%dmp3,param%multipole%dmp5,
     $      param%multipole%shift,param%multipole%kexp,
     $      param%multipole%rmax

      call util_flush(luout)
      deallocate(info)

      end
