*
* $Id: psi_new.F,v 1.4 2001-12-05 01:51:20 bylaska Exp $
*

*     ***********************************
*     *             					*
*     *           psi_new				*
*     *             					*
*     ***********************************

      subroutine psi_new()
      implicit none 

#include "mafdecls.fh"


*    *** local variables ***
      logical success
      real*8  icharge,tcharge,t
      integer mult,x
      integer ispin,ne(2)
      integer version
      integer ngrid(3)
      real*8  unita(3,3)
      character*30 filename


*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_ispin
      external control_ngrid,control_ispin
      integer  control_multiplicity
      external control_multiplicity
      real*8   control_totalcharge
      external control_totalcharge
      real*8   ion_totalcharge
      external ion_totalcharge
      logical  lcao_wvfnc
      external lcao_wvfnc

      filename = control_input_psi()
      version = 3
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
      tcharge = control_totalcharge()
      icharge = ion_totalcharge()
      t = icharge - tcharge

      if (ispin.eq.1) then
         ne(1) = t/2
         ne(2) = 0
      else
         x = t
         if (mod(x,2).ne.0) mult = 2
         x = mult-1
         ne(1) = (t+x)/2
         ne(2) = (t-x)/2
      end if

      success =  lcao_wvfnc(.false.,10,
     >                      filename,
     >                      version,
     >                      ngrid,
     >                      unita,
     >                      ispin,
     >                      ne)

      if (.not.success)
     >   call wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  ispin,
     >                  ne)


      return
      end
