C
C $Id: hess_anal.F,v 1.5 2000-07-26 23:39:24 windus Exp $
C
      logical function hess_anal(rtdb)
C
C  This function is the molecular calculation module for calculating
C  analytic second derivatives.  Specifically, the second derivative
C  of the energy with respect to the nuclear coordinates will be
C  calculated.  Everything must be either read from or written to
C  the RTDB or files.
c
c  This routine assumes that task_hessian is the calling routine.  It
c  assumes that checks for unsupported wavefunctions has been done.
C
C  The notation here will follow that of Johnson and Frisch, J.Chem.Phys.
C  100, 7429-7442 (1994) unless otherwise annotated.
C
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb

      logical task_energy
      logical bgj_get_scf_method
      external task_energy
      external bgj_get_scf_method
c
      hess_anal = .false.
c
c  Broadcast information to all of the nodes on read
c
      if (.not. rtdb_parallel(.true.)) call errquit
     *  ('hess_anal: unable to set rtdb_parallel to true',555)
c
c  Set up the print mechanism (push down the stack and loads the print
c  information for hessian)
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'hessian')
c
c  Check to make sure that the energy and wavefunction are calculated
c
      if(.not.task_energy(rtdb)) call errquit 
     *  ('hess_anal: energy failure', 555)
      write(6,*) 'after the energy call'
      call flush(6)
c
c     initialize hessian variables
c
      if (.not.bgj_get_scf_method(theory)) call errquit
     *  ('hess_anal: failure to get theory type', 555)
      write(6,*) 'before the hess_init call'
      call flush(6)
      call hess_init(rtdb)
      write(6,*) 'after the hess_init call'
      call flush(6)
c
c  Get the one-electron second derivative contributions
c
      call onedd_cont(rtdb)
      write(6,*) 'after the onedd_cont call'
      call flush(6)
c
c  Get the nuclear contribution to the hessian
c
      call nucdd_cont(rtdb)
      write(6,*) 'The Hessian after nucdd contributions:'
      call flush(6)
      call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
c
c  Get the two-electron second derivative contributions
c
      call twodd_cont(rtdb)
      write(6,*) 'after the twodd_cont call'
      call flush(6)
c
c  solve the CPSCF problem
c
      call hess_cphf(rtdb)
      write(6,*) 'after the hess_cphf call'
      call flush(6)
c
c     dump hessian file
c
      call hess_wrt()
c
c  clean up any hessian specific stuff
c
      call hess_tidy()

      hess_anal = .true.
      return 
      end
