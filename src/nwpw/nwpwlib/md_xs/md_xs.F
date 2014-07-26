c $Id: md_xs.F 25353 2014-03-17 14:12:25Z jhammond $

c     **********************************************
c     *                                            *
c     *                md_xs_analysis              *
c     *                                            *
c     **********************************************
      subroutine md_xs_analysis(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid,MASTER,np
      parameter (MASTER=0)

      logical oprint
      character*20 version


      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

#ifdef FEFF
      call feff_version(version)
      if (oprint) then
         write(luout,*) 
     >   "==== ","MD XS analysis: FEFF version ",version," ===="
      end if
      
#else
      if (oprint) then
         write(luout,*) 
     >     "FEFF library not inlucded in compile: ",
     >     "MD XS analysis not performed"
      end if
#endif
      return
      end



#ifdef FEFF
c     **********************************************
c     *                                            *
c     *                md_xs                       *
c     *                                            *
c     **********************************************
      subroutine md_xs(rtdb,nion,symbols,rion)
      implicit none
      integer rtdb
      integer nion
      character*2 symbols(*)
      real*8      rion(3,*)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid,MASTER,np
      parameter (MASTER=0)

      integer bigstr(2)

*     **** external functions ****


      call Parallel_taskid(taskid)
      call Parallel_np(np)


      return
      end


#endif

