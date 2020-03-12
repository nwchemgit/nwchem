C> \file minit.F
C> Memory management routines
C>
C> \ingroup selci
C> @{
C>
C> \brief Initialize memory administration
C>
C> This routine sets the maximum amount of memory the selected CI
C> code is allowed to use. Obviously this whole memory system could
C> be implemented in terms of the usual MA capabilities, but for 
C> some reason that has not been done.
C>
      integer function selci_minit(llword)
*
* $Id$
*
csw   1
c     implicit integer (a-z)
      implicit none
      integer llword !< [Input] The maximum number of words
c
#include "errquit.fh"
#include "selci_minit.fh"
c
c     initialise the memory management ... simple single direction
c     stack ... removed ability to grow from top down also
c
c     llword = value to initialise lword
c     lword = no. of real*8 words in stack
c     lmin = retains maximum amount of core ever used (info only)
c     ltop = (ltop+1) points to next free word on stack
c
      lmin = 0
      ltop=0
      lword = llword
      selci_minit = 0
      return
      end
c
      integer function selci_mleft(junk)
c
c     return the no. of REAL words free
c
      implicit none
#include "selci_minit.fh"      
      selci_mleft = lword - ltop
      return
      end
c
      integer function selci_mptr(nword)
c
c     return pointer into REAL stack for nword REALS
c
      implicit none
#include "errquit.fh"      
#include "selci_minit.fh"      
      if (nword.lt.0) call errquit('MALOC: nword < 0',nword, INPUT_ERR)
      mfr = lword - ltop
      if (nword .le. mfr) then
        selci_mptr = ltop + 1
        ltop = ltop + nword
        lmin = max0(ltop,lmin)
        return
      else
         selci_mptr=0
        write (6,20) nword,mfr
 20     format(///' core request:',i10,5x,'available:',i10/)
        call errquit('MPTR: need this much more memory', nword-mfr,
     &       MEM_ERR)
      endif
      end
c
      integer function selci_mfree(iad)
c
c     free REAL memory back to this pointer
c
      implicit none
#include "selci_minit.fh"      
#include "errquit.fh"      
      if (iad.lt.1 .or. iad.gt.(ltop+1))
     $     call errquit('MFREE: invalid argument',iad, INPUT_ERR)
      ltop = iad - 1
      selci_mfree = 0
      return
      end
c
      integer function selci_mstats(junk)
c
c     print out statistics about core usage
c
      implicit none
#include "selci_minit.fh"      
      write (6,30) lmin,lword-lmin
 30   format(/' minimal core high water =',i10,
     1     ' real numbers (spare =',i10,')')
      selci_mstats = 0
      return
c
      end
C>
C> @}
