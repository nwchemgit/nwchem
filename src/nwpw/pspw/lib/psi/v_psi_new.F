*
* $Id: v_psi_new.F,v 1.1 2001-11-28 22:15:21 bylaska Exp $
*

*     ***********************************
*     *             					*
*     *           v_psi_new				*
*     *             					*
*     ***********************************

      subroutine v_psi_new(ispin,ne)
      implicit none 
      integer ispin,ne(2)

#include "mafdecls.fh"

*    *** local variables ***
      integer version
      integer ngrid(3)
      real*8  unita(3,3)
      character*30 filename

*     ***** local functions ****
      character*30 control_input_v_psi
      external     control_input_v_psi
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

      filename = control_input_v_psi()
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

      call v_wvfnc_new(filename,
     >               version,
     >               ngrid,
     >               unita,
     >               ispin,
     >               ne)


      return
      end
