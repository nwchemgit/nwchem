
c     ********************************** 
c     *                                *
c     *      ion_scaling_atoms_start   *
c     *                                *
c     ********************************** 
      subroutine ion_scaling_atoms_start(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "ion_scaling_atoms.fh"

      if(.not.btdb_get(rtdb,'nwpw:scaling_natms',
     >                 mt_int,1,scaling_natms)) then
         scaling_natms = -1
      end if

      if (scaling_natms.gt.0) then
         if (.not.BA_alloc_get(mt_int,scaling_natms,
     >           'ionscaling_atoms',scaling_atoms(2),scaling_atoms(1)))
     >      call errquit('scaling_atoms_init:out of memory',0,MA_ERR)

         if (.not.btdb_get(rtdb,'nwpw:scaling_atoms',
     >                mt_int,scaling_natms,int_mb(scaling_atoms(1))))
     >      call errquit('scaling_atoms_init',0,RTDB_ERR)
      end if

      return
      end

c     ********************************** 
c     *                                *
c     *      ion_scaling_atoms_end     *
c     *                                *
c     ********************************** 
      subroutine ion_scaling_atoms_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "ion_scaling_atoms.fh"

      if (scaling_natms.gt.0) then
         if (.not.BA_free_heap(scaling_atoms(2)))
     >      call errquit('scaling_atoms_init:freeing memory',0,MA_ERR)
      end if
      return
      end


c     ********************************** 
c     *                                *
c     *      ion_scaling_atoms_print   *
c     *                                *
c     ********************************** 
      subroutine ion_scaling_atoms_print(luout)
      implicit none
      integer luout

#include "bafdecls.fh"
#include "ion_scaling_atoms.fh"

      integer i

      if (scaling_natms.gt.0) then
        write(luout,10) (int_mb(scaling_atoms(1)+i-1),i=1,scaling_natms)
      end if 

      return
 10   FORMAT(5X, ' cooling/heating atoms: ',10I5)
      end


c     ********************************** 
c     *                                *
c     *      ion_scaling_atoms_natms   *
c     *                                *
c     ********************************** 
      integer function ion_scaling_atoms_natms()
      implicit none

#include "ion_scaling_atoms.fh"

      ion_scaling_atoms_natms = scaling_natms
      return
      end


c     ********************************** 
c     *                                *
c     *      ion_scaling_atoms         *
c     *                                *
c     ********************************** 
      subroutine ion_scaling_atoms(rti,nion,vion)
      implicit none
      real*8  rti
      integer nion
      real*8  vion(3,nion)

#include "bafdecls.fh"
#include "ion_scaling_atoms.fh"

*     **** local variables ****
      integer i,ii

      if (scaling_natms.gt.0) then
         do i=1,scaling_natms
            ii = int_mb(scaling_atoms(1)+i-1)
            vion(1,ii) = vion(1,ii)*rti
            vion(2,ii) = vion(2,ii)*rti
            vion(3,ii) = vion(3,ii)*rti
         end do
      else
         call dscal(3*nion,rti,vion,1)
      end if

      return
      end
