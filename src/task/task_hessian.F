      logical function task_hessian(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
c     Generic NWChem interface to compute the analytic hessian.  
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
c     task:numerical (logical) - optional - if true use numerical 
c         differentiation. If absent or false use default selection.
c     
c     RTDB output parameters no for analytic hessian at the moment.
c     ----------------------
c     task:hessian(real array) - derivative w.r.t. geometry cart. coords.
c     task:status (logical)- T/F for success/failure
c     task:cputime (real)  - cpu time to execute the task
c     task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
c     If the method does not have analytic derivatives automatically call
c     the numerical derivative routine (not true in my case at the time).
c
      logical  scf,scf_hessian
      external scf,scf_hessian
c
      logical status
      double precision cpu, wall
      character*80 theory, prefix

c
      task_hessian = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:hessian: theory not specified',0)
c
c     Actually do the deed
c
      prefix = theory           ! Most common
c
      if (theory .eq. 'scf') then
         status = scf_hessian(rtdb)
      else
         call errquit('task_hessian: unknown theory',0)
      endif
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
c     if (status) then
c        key = ' '
c        write(key,'(a,'':hessian'')') prefix(1:inp_strlen(prefix))
c     endif
c
      task_hessian = status
c
      end

