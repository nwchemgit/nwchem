C $Id: restor.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_restor(fu30,lhead,jop,lsad,llopt)
C
C Read reactant, product, and saddle point information from fu30 file
C
      Implicit None
C
#include "drdyP.fh"
#include "stdio.fh"
C
*::passed:
      Integer fu30
      Logical lhead,lsad
      Integer llopt(40)
      Integer jop
*::local:
      Character*1 line(80)
      Integer fu6
      Integer iop, i, j
      Integer nend, nfreq, nbar
C
      nend   = 0
      nfreq  = 0
      nbar   = 0
      fu6    = luout
      lhead  = .false.
      jop    = 0
      lsad   = .false.
*
      read (fu30,3000,end=98) ctitle
C      write (fu6,3000) ctitle
      read (fu30,3005,end=99) llopt
C      write (fu6,3005) llopt
      lhead = .true.
      iop = 1
      do while (iop.le.4)
C Reactants and products
        read (fu30,3000,end=100) line
C        write (fu6,3000) line
        if (iop.eq.1) then
          read (fu30,*,end=100) sr,vr,vzero
C          write (fu6,*) sr,vr,vzero
        else if (iop.eq.2) then
          read (fu30,*,end=100) s,v
C          write (fu6,*) s,v
        else
          read (fu30,*,end=100) sp,vp
C          write (fu6,*) sp,vp
        endif
        if (icode(iop).eq.1) then
          read (fu30,3005,end=100) ifrgrp(1,iop)
C          write (fu6,3005) ifrgrp(1,iop)
        elseif (icode(iop).eq.2) then
          write (fu6,6000) iop
          call errquit('drdy_restor: fatal error',911)
        else
          nend = 3*nratom(iop)
          if (icode(iop).eq.3) then
            nfreq = nend - 5
          else
            nfreq = nend - 6
          endif
          read (fu30,3005,end=100) (ifrgrp(i,iop),i=1,nfreq)
C          write (fu6,3005) (ifrgrp(i,iop),i=1,nfreq)
          read (fu30,3000,end=100) line
C          write (fu6,*) line
          nbar = (nend*(nend+1))/2
          read (fu30,*,end=100) (hgrp(i,iop),i=1,nbar)
C          write (fu6,*) (hgrp(i,iop),i=1,nbar)
        endif
        jop = iop
        if ((iop.eq.1 .and. lgs(6).le.2) .or.
     &       iop.eq.2 .or.
     &      (iop.eq.3 .and. lgs(6).eq.1 .or. lgs(6).eq.3) .or.
     &       iop.eq.4) then
          iop = iop + 1
        else
          iop = iop + 2
        endif
      enddo
c saddle point
      if (icode(5).eq.3) then
        nfreq = n3 - 6
      else
        nfreq = n3 - 7
      endif
      read (fu30,3000,end=101) line
C      write (fu6,*) line
      read (fu30,*,end=101) s,vspsv
C      write (fu6,*) s,vspsv
C      read (fu30,3000,end=101) line
c      read (fu30,*,end=101) (xspsv(i),i=1,n3)
      read (fu30,3005,end=101) (ifrsp(i),i=1,nfreq)
      read (fu30,3000,end=101) line
      do i = 1, n3
        read (fu30,*,end=101) (fspsv(j,i),j=1,i)
      enddo
      do i = 2,n3
        do j = 1,i-1
          fspsv(i,j) = fspsv(j,i)
        enddo
      enddo
      lsad = .true.
      return
98    continue
      write (fu6,*) ' trouble in restor with title'
      return
99    continue
      write (fu6,*) ' trouble in restor with lopt'
      return
100   continue
      write (fu6,*) ' trouble in restor for iop=',iop
      return
101   continue
      write (fu6,*) ' trouble in restor for saddle point'
      Return
3000  Format (80a1)
3005  Format (1x,i2,39i3)
6000  Format (' For iop=',i5,', icode=2, Morse diatomic, was chosen',
     *      /, 'This option is not supported, reset icode to 3')
      End
