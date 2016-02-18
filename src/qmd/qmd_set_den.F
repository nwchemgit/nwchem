c
c     qmd_set_den
c
      subroutine qmd_set_den(nstates,cstate,zden)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
#include "qmd_common.fh"
c
      integer nstates            ! # of states
      integer cstate             ! initial state
      double complex zden(nstates,nstates)
c
      integer i, j
c
      character*30 pname
c
c     preliminaries
      pname='qmd_set_den: '
c
      do i=1,nstates
        do j=1,nstates
          zden(j,i)=(0.d0,0.d0)
        end do
      end do
      zden(cstate,cstate)=(1.d0,0.d0)
c
      return
      end
