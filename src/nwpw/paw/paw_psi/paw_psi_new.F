*
* $Id$
*

*     ***********************************
*     *                                 *
*     *           paw_psi_new           *
*     *                                 *
*     ***********************************

      subroutine paw_psi_new()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"


*    *** local variables ***
      logical success
      real*8  icharge,tcharge,t
      integer mult,x,xx
      integer ispin,ne(2)
      integer version,rtdb
      integer ngrid(3)
      real*8  unita(3,3)
      logical lcao_print,lcao_skip
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
      real*8   control_TotalCharge
      external control_TotalCharge
      real*8   ion_TotalCharge_qm
      external ion_TotalCharge_qm
      logical  lcao_wvfnc
      external lcao_wvfnc
      integer  control_rtdb,control_version
      external control_rtdb,control_version

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


*     **** reassign spin to agree with total number of electrons ****
      if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
         ispin = 2
         call control_ispin_set(ispin)
      end if

*     **** reassign mult to agree with total number of electrons ****
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

*     **** assign number of up and down electrons ****
      if (ispin.eq.1) then
         ne(1) = xx/2
         ne(2) = 0
      else
         x = mult-1
         ne(1) = (xx+x)/2
         ne(2) = (xx-x)/2
      end if

      call wvfnc_new(filename,
     >               version,
     >               ngrid,
     >               unita,
     >               ispin,
     >               ne)


      return
      end

