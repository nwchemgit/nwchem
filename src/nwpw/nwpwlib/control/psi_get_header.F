*
* $Id: psi_get_header.F,v 1.2 2001-12-12 02:42:26 bylaska Exp $
*


#define TCGMSG
*     ***********************************
*     *             					*
*     *           psi_get_header		*
*     *             					*
*     ***********************************

      subroutine psi_get_header(version,nfft,unita,ispin,ne)
      implicit none 
      integer version
      integer nfft(3)
      real*8  unita(3,3)
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
      integer l
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi

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



      return
      end
