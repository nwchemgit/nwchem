*
* $Id: mympi.h,v 1.2 1997-11-04 10:07:54 d3e129 Exp $
*
      integer me, nproc, ierr, status(MPI_STATUS_SIZE), istatus
      common /distvars/ me, nproc
      common /statvars/ ierr, istatus, status 
