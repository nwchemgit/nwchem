      subroutine dft_scfcvg(rms,etold,etnew,Enuc,e_conv,d_conv,
     &   g_conv,ipol,iter,iterations,idone,converged,rtdb)
C     $Id: dft_scfcvg.F,v 1.9 1997-02-27 02:52:50 gutowski Exp $
      implicit none

      double precision rms(2)   ! [input]
      double precision etold    ! [input]
      double precision etnew    ! [output]
      double precision Enuc
      double precision e_conv, d_conv, g_conv
      integer ipol
      integer iter
      integer iterations
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
      integer me,iswitch
      double precision de,abde
c     
      idone=0
      converged = .false.
      me=ga_nodeid()
c     
c--   >  Check to see if convergence on density is required.
c     
      oprint_c=util_print('convergence',print_default)
      oprint_i=util_print('information',print_low)
      DENSITY = d_conv.gt.0
c     
c--   >  Evaluate change in energy.
c     
      de = etnew-etold
      abde = dabs(de)
c     
c--   >  Check to see if energy is converged.
c     
      if((abde.lt.e_conv).and.(.not.DENSITY))then
        idone=1
        converged = .true.
        status = rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged)
        return
      elseif((abde.lt.e_conv).and.(DENSITY))then
        idone=1
      endif
c     
c--   >  Check for density matrix convergence.
c     
      if( DENSITY )then
        iswitch=0
        if(dsqrt(rms(1)).le.d_conv)iswitch=1
        idone=idone*iswitch
c        write(*,*)' iswitch, idone, rms(1): ',
c     &              iswitch, idone, rms(1)
        if(ipol.eq.2)then
          iswitch=0
          if(dsqrt(rms(2)).le.d_conv)iswitch=1
          idone=idone*iswitch
c          write(*,*)' iswitch, idone, rms(2): ',
c     &                iswitch, idone, rms(2)
        endif
        if(idone.eq.1)then
          converged = .true.
          status = rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged)
          return
        endif
      endif
c     
c--   >  Check iteration value.
c     
 1000 if(iter.lt.1)then
        return
      elseif(iter.eq.iterations)then
        idone=1
        status = rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged)
     &       
        return
      endif
c     
      return
c     
c     
      end
