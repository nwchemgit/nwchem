*
* $Id$
*

*     ***********************************
*     *             			*
*     *           v_cpsi_new		*
*     *             			*
*     ***********************************

      subroutine v_cpsi_new(spin_orbit,ispin,ne,nbrillioun)
      implicit none 
      logical spin_orbit
      integer ispin,ne(2),nbrillioun

#include "bafdecls.fh"

*    *** local variables ***
      integer version
      integer ngrid(3)
      real*8  unita(3,3)
      character*50 filename

*     ***** local functions ****
      character*50 control_input_v_psi
      external     control_input_v_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      filename = control_input_v_psi()
      version = 5
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

      call c_v_wvfnc_new(filename,
     >               version,
     >               ngrid,
     >               unita,
     >               ispin,
     >               ne,nbrillioun,spin_orbit)

      return
      end


*     ***********************************
*     *             			*
*     *           v_cpsi_delete		*
*     *             			*
*     ***********************************

      subroutine v_cpsi_delete()
      implicit none 

#include "bafdecls.fh"

*     *** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)

      logical      found
      integer      l
      character*50  filename
      character*54  filename2
      character*255 full_filename,full_filename2

*     ***** local functions ****
      character*50 control_input_v_psi
      external     control_input_v_psi


      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         filename = control_input_v_psi()

         l = index(filename,' ') - 1
         filename2 = filename(1:l)//".old"

         call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         call util_file_name_noprefix(filename2,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)

         inquire(file=full_filename,exist=found)
         if (found) then
           write(*,109) filename2
           write(*,110) filename 
  109      format(' vmovecs filename moved to: ',A)
  110      format(' vmovecs deleted          : ',A)

           call util_file_copy(full_filename,full_filename2)
           call util_file_unlink(full_filename)
         end if
      end if

      return
      end
