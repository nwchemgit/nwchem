*
* $Id: psi_get_ne.F,v 1.1 2001-12-14 20:02:46 bylaska Exp $
*

#define TCGMSG
*     ***********************************
*     *             			*
*     *           psi_get_ne		*
*     *             			*
*     ***********************************

      subroutine psi_get_ne(ispin,ne)
      implicit none 
      integer ispin,ne(2)


#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer version,l
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)


*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(4,full_filename,l,'r',l)
         call iread(4,version,1)
         call iread(4,nfft,3)
         call dread(4,unita,9)
         call iread(4,ispin,1)
         call iread(4,ne,2)
         call closefile(4)
      end if

c     **** send header to all nodes ****
#ifdef MPI
      msglen = 1
      call MPI_BCAST(version,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 3
      call MPI_BCAST(nfft,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 9
      call MPI_BCAST(unita,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 1
      call MPI_BCAST(ispin,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 2
      call MPI_BCAST(ne,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
      msglen = 3
      call BRDCST(9+MSGINT,nfft,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,ispin,mitob(msglen),MASTER)
      msglen = 2
      call BRDCST(9+MSGINT,ne,mitob(msglen),MASTER)
#endif


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         if (taskid.eq.MASTER) then
            write(6,*) "Error reading psi"
         end if
      end if

      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
         if (taskid.eq.MASTER) then
            write(6,*) "Error reading psi"
         end if

      end if


      return
      end
