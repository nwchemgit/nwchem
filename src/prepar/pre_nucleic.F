      logical function pre_nucleic(ffield,sgmnam)
c
c $Id: pre_nucleic.F,v 1.1 1999-06-16 23:56:59 d3j191 Exp $
c
      implicit none
c
      character*5 sgmnam
      character*80 ffield
c
      pre_nucleic=.true.
c
      if(sgmnam(1:3).eq.'A  '.or.sgmnam(1:3).eq.'+A '.or.
     + sgmnam(1:3).eq.'C  '.or.sgmnam(1:3).eq.'+C '.or.
     + sgmnam(1:3).eq.'G  '.or.sgmnam(1:3).eq.'+G '.or.
     + sgmnam(1:3).eq.'I  '.or.sgmnam(1:3).eq.'+I '.or.
     + sgmnam(1:3).eq.'T  '.or.sgmnam(1:3).eq.'+T '.or.
     + sgmnam(1:3).eq.'U  '.or.sgmnam(1:3).eq.'+U ') return
      if(ffield(1:5).eq.'amber') then
      if(sgmnam(1:3).eq.'DA '.or.sgmnam(1:3).eq.'RA '.or.
     + sgmnam(1:3).eq.'DC '.or.sgmnam(1:3).eq.'RC '.or.
     + sgmnam(1:3).eq.'DG '.or.sgmnam(1:3).eq.'RG '.or.
     + sgmnam(1:3).eq.'DT '.or.sgmnam(1:3).eq.'RU ') return
      endif
      if(ffield(1:6).eq.'charmm') then
      endif
c
      pre_nucleic=.false.
      return
      end
