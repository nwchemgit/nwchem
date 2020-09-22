*
* $Id$
*

*     ***********************************
*     *             			*
*     *           psi_new		*
*     *             			*
*     ***********************************

      subroutine psi_new()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"


*    *** local variables ***
      logical success,oddelcfill
      real*8  icharge,tcharge,t
      integer mult,x,xx
      integer ispin,ne(2),nextra_orbs
      integer version,rtdb
      integer ngrid(3)
      real*8  unita(3,3)
      logical lcao_print,lcao_skip,adjust_skip
      integer lcao_iterations
      character*50 filename


*     ***** external functions ****
      character*50 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_ispin
      external control_ngrid,control_ispin
      integer  control_multiplicity
      external control_multiplicity
      real*8   control_totalcharge
      external control_totalcharge
      real*8   ion_TotalCharge_qm
      external ion_TotalCharge_qm
      logical  lcao_wvfnc,control_print,psi_filefind
      external lcao_wvfnc,control_print,psi_filefind
      integer  control_rtdb,control_version
      external control_rtdb,control_version
      integer  control_frac_occ_extra_orbitals
      external control_frac_occ_extra_orbitals
      logical  control_mult_fixed
      external control_mult_fixed

      

      filename = control_input_psi()
      version = control_version()
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)

      ispin   = control_ispin()
      mult    = control_multiplicity()
      tcharge = control_TotalCharge()
      icharge = ion_TotalCharge_qm()
      t = icharge - tcharge       !** total number of electrons **
      x = NINT(t)
      xx = x

      oddelcfill = .false.

      if (.not.control_mult_fixed()) then
*        **** reassign spin to agree with total number of electrons ****
         if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
            ispin = 2
            call control_ispin_set(ispin)
         end if

*        **** reassign mult to agree with total number of electrons ****
         if ((mod(x,2).ne.0) .and.(mod(mult,2).ne.0)) then !** odd number of electrons and mult odd **
            mult = mult - 1
            do while (mult.gt.(x+1))
               mult = mult - 2
            end do
            if (mult.lt.1) mult = 2
            call control_multiplicity_set(mult)
         end if
         if ((mod(x,2).eq.0) .and.(mod(mult,2).eq.0)) then !** even number of electrons and mult even **
            mult = mult - 1
            do while (mult.gt.(x+1))
               mult = mult - 2
            end do
            if (mult.lt.1) mult = 1
            call control_multiplicity_set(mult)
         end if
      else
         if ((mod(x,2).ne.0).and.(ispin.eq.1)) oddelcfill = .true.  !** odd number of electrons **
      end if


*     **** assign number of up and down electrons ****
      if (ispin.eq.1) then
         ne(1) = xx/2
         ne(2) = 0
      else
         x = mult-1
         ne(1) = (xx+x)/2
         ne(2) = (xx-x)/2
      end if
      !write(*,*) "ispin,ne=",ispin,ne
      !write(*,*) "mult=",mult


      nextra_orbs = control_frac_occ_extra_orbitals()
      if (ispin.eq.1) then
         ne(1) = ne(1) + nextra_orbs
      else
         ne(1) = ne(1) + nextra_orbs
         ne(2) = ne(2) + nextra_orbs
      end if

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,
     >         'nwpw:lcao_iterations',mt_int,1,lcao_iterations))
     >   lcao_iterations = 10

      if (.not.btdb_get(rtdb,
     >         'nwpw:lcao_print',mt_log,1,lcao_print))
     >   lcao_print = .false.
      lcao_print = lcao_print.or.control_print(print_high)
      if (.not.btdb_get(rtdb,
     >         'nwpw:lcao_skip',mt_log,1,lcao_skip))
     >   lcao_skip = .true.
      if (.not.btdb_get(rtdb,
     >         'nwpw:adjust_skip',mt_log,1,adjust_skip))
     >   adjust_skip = .false.

*     **** adjust occupations of existing wavefunction ****
      if (psi_filefind().and.(.not.adjust_skip)) then
         call wvfnc_adjust(filename,ispin,ne)

*     **** make wavefunction from scratch ****
      else
         success = .false.
         if (.not.lcao_skip) then
            success =  lcao_wvfnc(lcao_print,lcao_iterations,
     >                        filename,
     >                        version,
     >                        ngrid,
     >                        unita,
     >                        ispin,
     >                        ne)
         end if

         if (.not.success) then
            call wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  ispin,
     >                  ne,oddelcfill)
         end if

      end if

      return
      end


*     ***********************************
*     *             			*
*     *           epsi_new		*
*     *             			*
*     ***********************************

      subroutine epsi_new()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"

*    *** local variables ***
      logical success,adjust_skip
      integer ispin,ne(2)
      integer version,rtdb
      integer ngrid(3)
      real*8  unita(3,3)
      character*50 filename

*     ***** local functions ****
      character*50 control_input_epsi
      external     control_input_epsi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_ispin,control_excited_ne
      external control_ngrid,control_ispin,control_excited_ne
      integer  control_rtdb,control_version
      external control_rtdb,control_version
      logical  epsi_filefind
      external epsi_filefind

      filename = control_input_epsi()
      version = control_version()
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)

      ispin   = control_ispin()
      ne(1) = 0 
      ne(2) = 0 
      ne(1) = control_excited_ne(1)
      if (ispin.eq.2) ne(2) = control_excited_ne(2)

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,
     >         'nwpw:adjust_skip',mt_log,1,adjust_skip))
     >   adjust_skip = .false.
      if (epsi_filefind().and.(.not.adjust_skip)) then
         call wvfnc_adjust(filename,ispin,ne)
      else
         call wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  ispin,
     >                  ne,.false.)
      end if
      return
      end

