      subroutine schwarz_print(natoms, nshell)
C$Id$
      implicit none
#include "schwarz.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "msgids.fh"
      integer natoms, nshell
      double precision nat_pairs, nat_quart, nsh_pairs, nsh_quart
C
C     Special array just to hold totals for printout so that we
C     can continue accumulating in the base data structures.
C     
      Integer totals(nscreen)
c     
c     print out a summary of the performance of the screening
c    
c     Atom pair info is done by all processes
C     But number of Fock builds is not cummulative
c
      call icopy(nscreen, iscreen, 1, totals, 1)
      call ga_igop(msg_schpr, totals(1), nscreen-1, '+')
      if (ga_nodeid().eq.0) then
         nat_pairs = natoms*(natoms+1)/2
         nat_quart = nat_pairs*(nat_pairs+1)/2
         nsh_pairs = nshell*(nshell+1)/2
         nsh_quart = nsh_pairs*(nsh_pairs+1)/2
         write(6,5) totals(6),
c$$$     $        Totals(1), 100.d0*Totals(1)/(nat_pairs * Totals(6)),
c$$$     $        Totals(2), 100.d0*Totals(2)/(nat_quart * Totals(6)),
     $        Totals(3), 100.d0*Totals(3)/(nat_quart * Totals(6)),
     $        Totals(4), 100.d0*Totals(4)/(nsh_quart * Totals(6)),
     $        Totals(5), 100.d0*Totals(5)/(nsh_quart * Totals(6))
 5       format(/
     $        ' Integral screening statistics'/
     $        ' -----------------------------'//
     $        '       Direct Fock builds = ',i10/
c$$$     $        '      Atom pairs rejected = ',i10,' (',f6.2,'%)'/
c$$$     $        '   Atom quartets rejected = ',i10,' (',f6.2,'%)'/
     $        '   Atom quartets computed = ',i10,' (',f6.2,'%)'/
     $        '  Shell quartets rejected = ',i10,' (',f6.2,'%)'/
     $        '  Shell quartets computed = ',i10,' (',f6.2,'%)'/)
      endif
c     
      end
