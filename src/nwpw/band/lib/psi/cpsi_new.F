*
* $Id$
*

*     ***********************************
*     *             			*
*     *           cpsi_new		*
*     *             			*
*     ***********************************

      subroutine cpsi_new()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"


*    *** local variables ***
      logical success,adjust_skip,oddelcfill
      real*8  icharge,tcharge,t
      integer mult,x,xx,nextra_orbs
      integer ispin,ne(2),nbrillioun
      integer version,rtdb
      integer ngrid(3)
      real*8  unita(3,3)
      logical lcao_print,lcao_skip,soc
      integer lcao_iterations
      character*50 filename


*     ***** local functions ****
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
      logical  control_spin_orbit,psi_filefind
      external control_spin_orbit,psi_filefind
      real*8   ion_TotalCharge_qm
      external ion_TotalCharge_qm
c     logical  lcao_wvfnc
c     external lcao_wvfnc
      integer  control_rtdb,brillioun_nbrillioun
      external control_rtdb,brillioun_nbrillioun
      integer  control_frac_occ_extra_orbitals
      external control_frac_occ_extra_orbitals
      logical  control_mult_fixed
      external control_mult_fixed

      filename = control_input_psi()
      version = 5   !version 5-band structure, version 3- pspw
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

      soc     = control_spin_orbit()
      ispin   = control_ispin()
      mult    = control_multiplicity()
      tcharge = control_totalcharge()
      icharge = ion_TotalCharge_qm()
      t = icharge - tcharge   !*** total number of electrons ***
      x = NINT(t)
      xx = x

      oddelcfill = .false.

      !**** spin orbit ****
      if (soc) then
         ispin=2
         ne(1)=t
         ne(2)=t


      !**** not spin orbit ****
      else

         if (.not.control_mult_fixed()) then
*           **** reassign spin to agree with total number of electrons ****
            if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
               ispin = 2
               call control_ispin_set(ispin)
            end if

*           **** reassign mult to agree with total number of electrons ****
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

*        **** assign number of up and down electrons ****
         if (ispin.eq.1) then
            ne(1) = xx/2
            ne(2) = 0
         else
            x = mult-1
            ne(1) = (xx+x)/2
            ne(2) = (xx-x)/2
         end if

      end if

      nbrillioun = brillioun_nbrillioun()

      nextra_orbs = control_frac_occ_extra_orbitals()
      if (ispin.eq.1) then
         ne(1) = ne(1) + nextra_orbs
      else
         ne(1) = ne(1) + nextra_orbs
         ne(2) = ne(2) + nextra_orbs
      end if

      rtdb = control_rtdb()
c      if (.not.btdb_get(rtdb,
c     >         'nwpw:lcao_iterations',mt_int,1,lcao_iterations))
c     >   lcao_iterations = 10
c      if (.not.btdb_get(rtdb,
c     >         'nwpw:lcao_print',mt_log,1,lcao_print))
c     >   lcao_print = .false.
c      if (.not.btdb_get(rtdb,
c     >         'nwpw:lcao_skip',mt_log,1,lcao_skip))
c     >   lcao_skip = .false.
      if (.not.btdb_get(rtdb,
     >         'nwpw:adjust_skip',mt_log,1,adjust_skip))
     >   adjust_skip = .false.


      if (psi_filefind().and.(.not.adjust_skip)) then
         call c_wvfnc_adjust(filename,ispin,ne,nbrillioun)
      else
         success = .false.
c        if (.not.lcao_skip .and. const)
c    >     success =  lcao_wvfnc(lcao_print,lcao_iterations,
c    >                           filename,
c    >                           version,
c    >                           ngrid,
c    >                           unita,
c    >                           ispin,
c    >                           ne)
         if (.not.success) then
            version = 5 
            call c_wvfnc_new(filename,
     >                       version,
     >                       ngrid,
     >                       unita,
     >                       soc,
     >                       ispin,
     >                       ne,
     >                       nbrillioun,oddelcfill)
         end if
      end if

      return
      end



*     ***********************************
*     *             			*
*     *           cpsi_new_ne		*
*     *             			*
*     ***********************************

      subroutine cpsi_new_ne(ispin,ne)
      implicit none 
      integer ispin,ne(2)

#include "bafdecls.fh"


*    *** local variables ***
      logical success,soc
      integer nbrillioun
      integer version
      integer ngrid(3)
      real*8  unita(3,3)
      character*50 filename


*     ***** local functions ****
      character*50 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      logical  control_spin_orbit
      external control_spin_orbit

      filename = control_input_psi()
      version = 5   !version 5-band structure, version 3- pspw
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
       
      soc=control_spin_orbit()
      nbrillioun = brillioun_nbrillioun()
      call c_wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  soc,
     >                  ispin,
     >                  ne,
     >                  nbrillioun,.false.)
      return
      end



*     ***********************************
*     *                                 *
*     *          ecpsi_new              *
*     *                                 *
*     ***********************************

      subroutine ecpsi_new()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"

*    *** local variables ***
      logical success,soc
      integer ispin,ne(2),nbrillioun
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
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      logical  control_spin_orbit
      external control_spin_orbit

      filename = control_input_epsi()
      version = 5   !version 5-band structure, version 3- pspw
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

      soc=control_spin_orbit()
      nbrillioun = brillioun_nbrillioun()
      call c_wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  soc,
     >                  ispin,
     >                  ne,
     >                  nbrillioun,.false.)

      return
      end


