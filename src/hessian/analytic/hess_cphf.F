C
C $Id: hess_cphf.F,v 1.1 1999-05-06 17:29:34 windus Exp $
C
      subroutine hess_cphf(rtdb)
c
c  Solves the CPHF equations to get the responses
c
      implicit none
c
#include "global.fh"
#include "hess_info.fh"
#include "util.fh"
c
      integer rtdb, i
      character*256 cphf_rhs, cphf_sol
      logical hf_job
c
      logical  cphf, file_write_ga, file_read_ga
      external cphf, file_write_ga, file_read_ga
c
      hf_job = theory .eq. 'HF'
c
c Put the rhs info on disk and then pick up the solutions.
c For now this is done with extra memory, eventually (if sticking with a
c CPHF call for each "atomic displacement") this could be taken to one
c solution matrix that gets used ncent*3 times and processed before making
c the next cphf call.
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
c
      do i = 1,ncent*3
       if(.not.file_write_ga(cphf_rhs,g_rhs(i)))
     1     call errquit('could not write cphf_rhs',i)
       if(.not.cphf(rtdb))call errquit('could not do cphf',i)
c
c  g_rhs_xc is now holding the solution
c
       if(.not.file_read_ga(cphf_sol,g_rhs_xc(i)))
     1    call errquit('could not read cphf_rhs',i)
       call ga_print(g_rhs_xc(i))
      enddo

c     call cphf_poliz(rtdb)
      end
