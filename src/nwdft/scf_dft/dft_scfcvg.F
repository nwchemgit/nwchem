      subroutine dft_scfcvg(rms,etold,etnew,Enuc,scfcon,
     &igcon,iter,itrscf,idone,converged,rtdb)
C$Id: dft_scfcvg.F,v 1.4 1996-01-09 17:52:55 og845 Exp $
      implicit none

      double precision rms   ! [output]
      double precision etold ! [input]
      double precision etnew ! [output]
      double precision Enuc
      double precision scfcon
      integer igcon
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
c-->  Check to see if convergence on density is required.
c 
      oprint_c=util_print('convergence',print_high)
      oprint_i=util_print('information',print_low)
      DENSITY = igcon.gt.0
c
c-->  Evaluate change in energy.
c
      de = etnew-etold
      abde = dabs(de)
c
c-->  Check to see if energy increased.
c
      if( me.eq.0.and.oprint_c) then
       if(de.gt.0.D0) write(LuOut,9050)
      end if
c
c-->  Check to see if energy is converged.
c
      if((abde.lt.scfcon).and.(.not.DENSITY))then
        idone=1
        if( me.eq.0.and.oprint_i)then
         write(LuOut,9150)
         write(LuOut,9100)etnew,de
C         write(LuOut,9200)Enuc
         write(LuOut,9250)Enuc+etnew
        end if
        converged = .true.
        status = rtdb_put(
     &       rtdb, 'dft:converged', MT_LOG, 1, converged)
        return
      elseif((abde.lt.scfcon).and.(DENSITY))then
        idone=1
      endif
c
c-->  Check for density matrix convergence.
c
      if( DENSITY )then

          iswitch=0
      
          if(dsqrt(rms).le.10.0**(-igcon))iswitch=1

          idone=idone*iswitch

          if(idone.eq.1)then
            if( me.eq.0.and.oprint_i)then
             write(LuOut,9350)
             write(LuOut,9100)etnew,de
             write(LuOut,9400)dsqrt(rms)
C             write(LuOut,9200)Enuc
             write(LuOut,9250)Enuc+etnew
            end if
            converged = .true.
            status = rtdb_put(
     &           rtdb, 'dft:converged', MT_LOG, 1, converged)
            return
          endif

       endif
c      
c--    >  Check iteration value.
c      
 1000  if(iter.lt.1)then
         return
       elseif(iter.lt.itrscf)then
         if( me.eq.0.and.oprint_i )then
           write(LuOut,9100)etnew,de
           write(LuOut,9250)Enuc+etnew
           if(DENSITY)write(LuOut,9400)dsqrt(rms)
         end if
       elseif(iter.eq.itrscf)then
         idone=1
         if( me.eq.0.and.oprint_i)then
           write(LuOut,9300)
           write(LuOut,9100)etnew,de
           write(LuOut,9200)Enuc
           write(LuOut,9250)Enuc+etnew
         end if
         status = rtdb_put(
     &        rtdb, 'dft:converged', MT_LOG, 1, converged)
         return
       endif
c      
       return
c
9000  format(/,5x,'Iteration number',I4)
9050  format(/,'SCF energy increased this iteration.')
9100  format(/,15x,'Electronic energy:        ',f20.10,
     &       /,15x,'Delta E:                  ',e20.6)
9150  format(/,'SCF Convergence criterion on energy met.')
9200  format(/,15x,'Nuclear repulsion energy: ',f20.10)
9250  format(/,15x,'Total energy:             ',f20.10)
9300  format(/,5x,'Maximum number SCF cycles exceeded.')
9350  format(/,5x,'Requested convergence criteria on energy and',
     &         ' density met.')
9400  format(15x,'RMS change in density:    ',e20.6)
c
      end
