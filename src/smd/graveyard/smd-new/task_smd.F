* $Log: not supported by cvs2svn $
* Revision 1.1  2008/04/18 17:40:32  marat
* Adding serial version of smd code
*
* Revision 1.8  2008/01/16 22:27:03  marat
* added nearly all input statements
*
* Revision 1.7  2008/01/16 06:09:58  marat
* added parsing of short range parameters
*
* Revision 1.6  2008/01/15 18:48:23  marat
* smd updates
*
* Revision 1.5  2008/01/12 00:44:40  marat
* updates to smd
*
* Revision 1.4  2008/01/04 00:00:33  marat
* Fixed uninitialised variables problem
*
* Revision 1.3  2008/01/03 18:52:09  marat
* restoring original files for now
*
* Revision 1.1  2007/12/28 01:11:09  bylaska
* ...EJB
*
* Revision 1.3  2007/12/27 23:22:54  bylaska
* ...EJB
*
* Revision 1.25  2007/11/17 22:45:48  bylaska
* ...EJB

      logical function task_smd(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

c      include 'p_input.inc'
c      include 'p_array.inc'
c
c      integer iseed,istep
c      integer ntype,nprnt,nstep,nequl
c      integer ncons,ntcons,consatm
c      integer nbond,ntbond,bondatm
c      integer nshel,ntshel,shelatm
c
c      real*8 ekin,ecoul,eshrt,ebond,eshel,ewald1
c      real*8 tstep,ivv,etime1,etime2
c
c      logical lveloc,lupdate
c
c      dimension ivv(mxatms,3)
c      dimension consatm(mxcons,2),bondatm(mxbond,2),shelatm(mxshel,2)
c      character*80 infile,outfile

      write(*,*) "Hello from task_smd"

c      call smd_init(rtdb)
c      call smd_end(rtdb)
c      goto 500
c      if (.not.rtdb_cget(rtdb,'smd:inputfile',1,infile))
c     >   infile='INPUT'
c
c      if (.not.rtdb_cget(rtdb,'smd:outputfile',1,outfile))
c     >   outfile='OUTPUT'
c
c      write(*,*) "inputfile=",infile
c      write(*,*) "outputfile=",outfile
c
c      etime1=0.0
c      etime2=0.0
c      call cpu_time(etime1)
c
c      call util_flush(6)
c      open(unit=output,file=outfile)
c
c      write(*,*) "processng input file"
c      call util_flush(6)
c
c      call md_inpt(infile,
c     $    iseed,tstep,nstep,nequl,nprnt,ntype,ncons,consatm,
c     $    nbond,bondatm,nshel,shelatm,lveloc)
c
c      write(*,*) "finished processing input file"
c      call util_flush(6)
c      call util_flush(output)
c
c      call md_init(iseed,ntype,ncons,ntcons,consatm,nbond,ntbond,
c     $    bondatm,nshel,ntshel,shelatm,lveloc,ewald1)
c
c      write(*,*) "finished initialization"
c      call util_flush(6)
c      call util_flush(output)
c
c      call cpu_time(etime2)
c      write(output,'(/,a,f20.3)')'Set-up CPU time : ',(etime2-etime1)
c
c      do istep=1,nstep
c
c       call verlt_test(tstep,ivv,lupdate)
c
c       write(*,*) "finished verlt_test"
c       call util_flush(6)
c       call util_flush(output)
c
c       if(lupdate)call list_verlt()
c
c       call md_frce(ntype,ecoul,eshrt,ebond,ntbond,eshel,ntshel,ewald1)
c
c       write(*,*) "finished md_frce"
c       call util_flush(6)
c       call util_flush(output)
c
c       if(ntcons.eq.0)call inte_leapf(tstep,ekin)
c       if(ntcons.gt.0)call inte_shake(tstep,ntcons,ekin)
c 
c       if(istep.le.nequl)call md_scle(ntshel)
c
c       call print_output(istep,nprnt,ekin,ecoul,eshrt,ebond,eshel)
c
c      enddo
c
c      call print_final()
c
c      call cpu_time(etime2)
c      write(output,'(/,a,f20.3)')'Total CPU time : ',(etime2-etime1)
c
c500   continue
      write(*,*) "Goodbye from task_smd"
      task_smd = .true.
      return
      END
c $Id$
