      subroutine dft_scfcvg(rms,etold,etnew,Enuc,scfcon,
     &     ipol,igcon,iter,itrscf,idone,converged,rtdb)
C     $Id: dft_scfcvg.F,v 1.6 1996-10-01 21:31:05 d3h449 Exp $
      implicit none

      double precision rms(2)   ! [input]
      double precision etold    ! [input]
      double precision etnew    ! [output]
      double precision Enuc
      double precision scfcon
      integer igcon
      integer ipol
      integer iter
      integer itrscf
      integer idone
      integer rtdb
      logical converged

#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "rtdb.fh"
c     
      logical DENSITY,status,oprint_c,oprint_i
      integer me,nproc,iswitch
      double precision de,abde
c     
      idone=0
      converged = .false.
      me=ga_nodeid()
c     
c--   >  Check to see if convergence on density is required.
c     
      oprint_c=util_print('convergence',print_high)
      oprint_i=util_print('information',print_low)
      DENSITY = igcon.gt.0
c     
c--   >  Evaluate change in energy.
c     
      de = etnew-etold
      abde = dabs(de)
c     
c--   >  Check to see if energy is converged.
c     
      if((abde.lt.scfcon).and.(.not.DENSITY))then
        idone=1
        converged = .true.
        status = rtdb_put(
     &       rtdb, 'dft:converged', MT_LOG, 1, converged)
        return
      elseif((abde.lt.scfcon).and.(DENSITY))then
        idone=1
      endif
c     
c--   >  Check for density matrix convergence.
c     
      if( DENSITY )then
        iswitch=0
        if(dsqrt(rms(1)).le.10.0**(-igcon))iswitch=1
        idone=idone*iswitch
c        write(*,*)' iswitch, idone, rms(1): ',
c     &              iswitch, idone, rms(1)
        if(ipol.eq.2)then
          iswitch=0
          if(dsqrt(rms(2)).le.10.0**(-igcon))iswitch=1
          idone=idone*iswitch
c          write(*,*)' iswitch, idone, rms(2): ',
c     &                iswitch, idone, rms(2)
        endif
        if(idone.eq.1)then
          converged = .true.
          status = rtdb_put(
     &         rtdb, 'dft:converged', MT_LOG, 1, converged)
          return
        endif
      endif
c     
c--   >  Check iteration value.
c     
 1000 if(iter.lt.1)then
        return
      elseif(iter.eq.itrscf)then
        idone=1
        status = rtdb_put(
     &       rtdb, 'dft:converged', MT_LOG, 1, converged)
        return
      endif
c     
      return
c     
c     
      end
