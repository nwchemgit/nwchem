      integer function selci_minit(llword)
*
* $Id: minit.F,v 1.1 2003-10-17 22:32:27 carlfahl Exp $
*
csw   1
      implicit integer (a-z)
#include "errquit.fh"
      common /selci_corep/ lword, ltop, lmax, lmin
      save /selci_corep/
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
c
      entry selci_mleft(junk)
c
c     return the no. of REAL words free
c
      selci_mleft = lword - ltop
      return
c
      entry selci_mptr(nword)
c
c     return pointer into REAL stack for nword REALS
c
      if (nword.lt.0) call errquit('MALOC: nword < 0',nword, INPUT_ERR)
      mfr = lword - ltop
      if (nword .le. mfr) then
        selci_mptr = ltop + 1
        ltop = ltop + nword
        lmin = max0(ltop,lmin)
        return
      else
        write (6,20) nword,mfr
 20     format(///' core request:',i10,5x,'available:',i10/)
        call errquit('MPTR: need this much more memory', nword-mfr,
     &       MEM_ERR)
      endif
c
      entry selci_mfree(iad)
c
c     free REAL memory back to this pointer
c
      if (iad.lt.1 .or. iad.gt.(ltop+1))
     $     call errquit('MFREE: invalid argument',iad, INPUT_ERR)
      ltop = iad - 1
      selci_mfree = 0
      return
c
      entry selci_mstats(junk)
c
c     print out statistics about core usage
c
      write (6,30) lmin,lword-lmin
 30   format(/' minimal core high water =',i10,
     1     ' real numbers (spare =',i10,')')
      selci_mstats = 0
      return
c
      end
