*
* $Id: cpsi_new.F,v 1.3 2002-05-03 19:44:02 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           cpsi_new		*
*     *             			*
*     ***********************************

      subroutine cpsi_new()
      implicit none 

#include "mafdecls.fh"
#include "rtdb.fh"


*    *** local variables ***
      logical success
      real*8  icharge,tcharge,t
      integer mult,x
      integer ispin,ne(2),nbrillioun
      integer version,rtdb
      integer ngrid(3)
      real*8  unita(3,3)
      logical lcao_print,lcao_skip
      integer lcao_iterations
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
c     logical  lcao_wvfnc
c     external lcao_wvfnc
      integer  control_rtdb,brillioun_nbrillioun
      external control_rtdb,brillioun_nbrillioun

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
      nbrillioun = brillioun_nbrillioun()



      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,
     >         'nwpw:lcao_iterations',mt_int,1,lcao_iterations))
     >   lcao_iterations = 10
      if (.not.rtdb_get(rtdb,
     >         'nwpw:lcao_print',mt_log,1,lcao_print))
     >   lcao_print = .false.
      if (.not.rtdb_get(rtdb,
     >         'nwpw:lcao_skip',mt_log,1,lcao_skip))
     >   lcao_skip = .false.
    
      success = .false.
c     if (.not.lcao_skip)
c    >  success =  lcao_wvfnc(lcao_print,lcao_iterations,
c    >                        filename,
c    >                        version,
c    >                        ngrid,
c    >                        unita,
c    >                        ispin,
c    >                        ne)

      if (.not.success)
     >   call c_wvfnc_new(filename,
     >                  version,
     >                  ngrid,
     >                  unita,
     >                  ispin,
     >                  ne,
     >                  nbrillioun)


      return
      end
