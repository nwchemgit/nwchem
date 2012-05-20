* $Log: not supported by cvs2svn $
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

      include 'p_input.inc'
      include 'p_array.inc'

      integer iseed,istep
      integer ntype,nprnt,nstep,nequl
      integer ncons,ntcons,consatm
      integer nbond,ntbond,bondatm
      integer nshel,ntshel,shelatm

      real*8 ekin,ecoul,eshrt,ebond,eshel,ewald1
      real*8 tstep,ivv,etime1,etime2

      logical lveloc,lupdate

      dimension ivv(mxatms,3)
      dimension consatm(mxcons,2),bondatm(mxbond,2),shelatm(mxshel,2)
      character*80 infile,outfile

      write(*,*) "Hello from task_smd"

      call smd_init(rtdb)
      call smd_end(rtdb)
      goto 500
      if (.not.rtdb_cget(rtdb,'smd:inputfile',1,infile))
     >   infile='INPUT'

      if (.not.rtdb_cget(rtdb,'smd:outputfile',1,outfile))
     >   outfile='OUTPUT'

      write(*,*) "inputfile=",infile
      write(*,*) "outputfile=",outfile

      etime1=0.0
      etime2=0.0
      call cpu_time(etime1)

      call util_flush(6)
      open(unit=output,file=outfile)

      write(*,*) "processng input file"
      call util_flush(6)

      call md_inpt(infile,
     $    iseed,tstep,nstep,nequl,nprnt,ntype,ncons,consatm,
     $    nbond,bondatm,nshel,shelatm,lveloc)

      write(*,*) "finished processing input file"
      call util_flush(6)
      call util_flush(output)

      call md_init(iseed,ntype,ncons,ntcons,consatm,nbond,ntbond,
     $    bondatm,nshel,ntshel,shelatm,lveloc,ewald1)

      write(*,*) "finished initialization"
      call util_flush(6)
      call util_flush(output)

      call cpu_time(etime2)
      write(output,'(/,a,f20.3)')'Set-up CPU time : ',(etime2-etime1)

      do istep=1,nstep

       call verlt_test(tstep,ivv,lupdate)

       write(*,*) "finished verlt_test"
       call util_flush(6)
       call util_flush(output)

       if(lupdate)call list_verlt()

       call md_frce(ntype,ecoul,eshrt,ebond,ntbond,eshel,ntshel,ewald1)

       write(*,*) "finished md_frce"
       call util_flush(6)
       call util_flush(output)

       if(ntcons.eq.0)call inte_leapf(tstep,ekin)
       if(ntcons.gt.0)call inte_shake(tstep,ntcons,ekin)
 
       if(istep.le.nequl)call md_scle(ntshel)

       call print_output(istep,nprnt,ekin,ecoul,eshrt,ebond,eshel)

      enddo

      call print_final()

      call cpu_time(etime2)
      write(output,'(/,a,f20.3)')'Total CPU time : ',(etime2-etime1)

500   continue
      write(*,*) "Goodbye from task_smd"
      task_smd = .true.
      return
      END
c $Id$
