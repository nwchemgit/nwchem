*
* $Id: cpsi_read.F,v 1.3 2002-02-03 01:54:12 bylaska Exp $
*

#define TCGMSG
*     ***************************
*     *          	   			*
*     *       cpsi_read			*
*     *             			*
*     ***************************
      subroutine cpsi_read(ispin,ne,nbrillioun,psi2)
      implicit none 
      integer ispin,ne(2),nbrillioun
      double complex psi2(*)

#include "mafdecls.fh"
#include "rtdb.fh"

#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer rtdb
      integer version,l
      integer nfft3d
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)
      logical brillioun_check,psi_nolattice
      integer n,nb,nbrill
      integer msglen,indx,indx1

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value

*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb
      external control_ngrid,control_rtdb

      call Parallel_taskid(taskid)
      call C3dB_nfft3d(1,nfft3d)


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('out of stack memory',0)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,ispin,1)
         call iread(5,ne,2)
         call iread(5,nbrill,1)
      end if

c     **** send header to all nodes ****
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
      msglen = 1
      call BRDCST(9+MSGINT,nbrill,mitob(msglen),MASTER)
#endif 


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        call errquit('Error reading wavefunctions', 0)
      end if

      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .false.

      if (.not.psi_nolattice) then
      if ( (dabs(unita(1,1)-control_unita(1,1)).gt.1.0d-6) .or.
     >     (dabs(unita(2,1)-control_unita(2,1)).gt.1.0d-6) .or.
     >     (dabs(unita(3,1)-control_unita(3,1)).gt.1.0d-6) .or.
     >     (dabs(unita(1,2)-control_unita(1,2)).gt.1.0d-6) .or.
     >     (dabs(unita(2,2)-control_unita(2,2)).gt.1.0d-6) .or.
     >     (dabs(unita(3,2)-control_unita(3,2)).gt.1.0d-6) .or.
     >     (dabs(unita(1,3)-control_unita(1,3)).gt.1.0d-6) .or.
     >     (dabs(unita(2,3)-control_unita(2,3)).gt.1.0d-6) .or.
     >     (dabs(unita(3,3)-control_unita(3,3)).gt.1.0d-6) ) then
            write(*,*) "unita(psi)      unita(lattice)"
            write(*,*) unita(1,1), control_unita(1,1)
            write(*,*) unita(2,1), control_unita(2,1)
            write(*,*) unita(3,1), control_unita(3,1)
            write(*,*) unita(1,2), control_unita(1,2)
            write(*,*) unita(2,2), control_unita(2,2)
            write(*,*) unita(3,2), control_unita(3,2)
            write(*,*) unita(1,3), control_unita(1,3)
            write(*,*) unita(2,3), control_unita(2,3)
            write(*,*) unita(3,3), control_unita(3,3)
        call errquit('Error reading wavefunctions', 1)
      end if
      end if

      if (.not.rtdb_get(rtdb,'nwpw:psi_brillioun_check',
     >                  mt_log,1,brillioun_check))
     >   brillioun_check = .true.

      if (brillioun_check) then
      if ( nbrill.ne.nbrillioun) then
        call errquit(
     >   'Error reading wavefunctions - wrong number of k-points', 1)
      end if
      end if


*     *************************
*     *** read in 3d blocks ***
*     *************************

*     * file may contain more data then necessary *
      if (nbrillioun.le.nbrill) then

        do nb=1,nbrillioun
        do n=1,(ne(1)+ne(2))
         indx = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2)) 
         call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)))
         call Cram_c_pack(nb,dcpl_mb(tmp2(1)))
         call Cram_c_Copy(nb,dcpl_mb(tmp2(1)),psi2(indx))
        end do
        end do

*     * file contains less data then necessary *
      else 

        do nb=1,nbrill
        do n=1,(ne(1)+ne(2))
         indx = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2)) 
         call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)))
         call Cram_c_pack(nb,dcpl_mb(tmp2(1)))
         call Cram_c_Copy(nb,dcpl_mb(tmp2(1)),psi2(indx))
        end do
        end do

*       **** use the last zone for filling in unknown zones ****
        do nb=nbrill+1,nbrillioun
        do n=1,(ne(1)+ne(2))
         indx1 = 1 + (n-1)*nfft3d + (nbrill-1)*nfft3d*(ne(1)+ne(2)) 
         indx  = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2)) 
         call Cram_c_Copy(nbrill,psi2(indx1),dcpl_mb(tmp2(1)))
         call Cram_c_unpack(nbrill,dcpl_mb(tmp2(1)))
         call Cram_c_pack(nb,dcpl_mb(tmp2(1)))
         call Cram_c_Copy(nb,dcpl_mb(tmp2(1)),psi2(indx))
        end do
        end do

      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))

*     end if

      return
      end
