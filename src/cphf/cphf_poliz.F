
       subroutine cphf_poliz( rtdb, geom, basis, movecs )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
c
       integer rtdb
       integer geom
       integer basis
       character*255 movecs
       integer ibas,jbas
       integer g
       character*255 integ_type
       logical oskel
c
c
c
       write(6,900) movecs
 900   format(/,1x,'Inside CPHF_POLIZ',
     $        /,'Movecs file:',a80)
c...collect the dipole moment integrals from int_mpole:
c
c
       call int_dip_ga(ibas, jbas, g, integ_type, oskel)
c
c
c
       return
       end
