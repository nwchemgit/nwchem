      subroutine argos_diana_average()
c
c $Id$
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      character*3 fmt
      character*80 fil
      integer i
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
c     get the average coordinates
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
      call argos_diana_xaver(dbl_mb(ndx),dbl_mb(i_xadat))
c
      if(card(7:7).eq.' ') return
c
      if(me.eq.0) then
c
      read(card(7:80),'(a)') fil
c
      i=index(fil,'.')
      if(i.le.0) then
      fmt='pdb'
      else
      fmt=filcop(i+1:i+3)
      endif
c
      if(fmt.eq.'pdb') then
      call argos_diana_wrtpdb(lfncop,fil,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xadat),int_mb(i_wsel),dbl_mb(i_wdat),.false.,.false.,
     + int_mb(i_ndxw))
      elseif(fmt.eq.'xyz') then
      call argos_diana_wrtxyz(lfncop,fil,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xadat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.false.,int_mb(i_ndxw))
      endif
c
      write(*,'(/,a,a)')
     +  ' Averaged coordinates written to ',
     + fil(1:index(fil,' ')-1)
c
      endif
c
      return
      end
