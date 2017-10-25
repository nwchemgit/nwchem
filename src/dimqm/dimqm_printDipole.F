      subroutine dimqm_printDipole(rtdb, muind, icmplx, lpprint)
c
c
c
c
c
      implicit none      
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"

c   Input variables
      integer rtdb
      integer icmplx
      double precision muind(3,nDIM,icmplx)
c      double precision qind(nDIM, icmplx)
      logical lpprint
c
c   Local variables
      double precision diptot(3, icmplx)
      character*60 dd
      character*60 d

      if(ga_nodeid().ne.0) return
      diptot = SUM(muind, DIM = 2)
      dd =
     $ '============================================================='
      d =
     $ '-------------------------------------------------------------'
c   BLAS matrix-matrix multiply to add in charge term

      if(lpprint) then
       write(luout,'(/1x,a)') dd
       write(luout,'(1x,a)')
     $                   'Total induced dipole moment in DIM system :'
       write(luout,'(1x,a)') d
       if(icmplx == 1) then
        write(luout,'(2x,a,2x,3f16.8)') 'A.U.:', diptot(:,1)
        write(luout,'(2x,a,1x,3f16.8)') 'Debye:', diptot(:,1)*AU2DEBYE
       else
        write(luout,'(1x,a)') 'A.U.'
        write(luout,'(1x,a,3f16.8)') 'Real: ', diptot(:,1)
        write(luout,'(1x,a,3f16.8)') 'Imag: ', diptot(:,2)
        write(luout,'(1x,a)') d
        write(luout,'(1x,a)') 'Debye'
        write(luout,'(1x,a,3f16.8)') 'Real: ', diptot(:,1)*AU2DEBYE
        write(luout,'(1x,a,3f16.8)') 'Imag: ', diptot(:,2)*AU2DEBYE
       end if
       write(luout,'(1x,a/)') dd
      else
       write(luout,'(1x,a)') 
     $                   'Total induced dipole moment in DIM system :'
       if(icmplx == 1) then
        write(luout,'(3f14.8)') diptot(:,1)
       else
         write (luout,'(1x,a,3f16.8)') 'Real: ', diptot(:,1)
         write (luout,'(1x,a,3f16.8)') 'Imag: ', diptot(:,2)
       end if
      end if
      
      end subroutine dimqm_printDipole
