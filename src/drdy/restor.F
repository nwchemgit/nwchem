C $Id$
      Subroutine drdy_restor(lhead,jop,lsad,llopt)
C
C Read reactant, product, and saddle point information from fu30 file
C
      Implicit None
#include "errquit.fh"
C
#include "drdyP.fh"
#include "stdio.fh"
C
*::passed:
      Integer jop,llopt(40)
      Logical lhead,lsad
*::local:
      Integer i,ierr,ifreq,ii,iop,j,nbar,nend,nfreq,nfreq1
      Character*80 line
C
      lhead = .false.
      jop = 0
      lsad = .false.
*
      read (fu30,3000,end=98) ctitle
C      write(fu6,3000) ctitle
      read (fu30,3005,end=99) llopt
C      write(fu6,3005) llopt
      lhead = .true.
      iop = 1
      do while (iop.le.4)
C Reactants and products
         read (fu30,3000,end=100) line
C         write(fu6,3000) line
         if (iop.eq.1) then
            if (lgs(39).eq.0) then
               read (fu30,*,end=100) sr,vr,vzero
            else
               read (fu30,*,end=100) sr,vr,vzero,vrsp,vzerosp
            endif
C            write(fu6,*) sr,vr,vzero
         else if (iop.eq.2) then
            if (lgs(39).eq.0) then
               read (fu30,*,end=100) s,v
            else
               read (fu30,*,end=100) s,v,vspc
            endif
C            write(fu6,*) s,v
         else if (iop.eq.3) then
            if (lgs(39).eq.0) then
               read (fu30,*,end=100) s,vp
            else
               read (fu30,*,end=100) s,vp,vpsp
            endif
C            write(fu6,*) sp,vp
         else
            if (lgs(39).eq.0) then
               read (fu30,*,end=100) s,v
            else
               read (fu30,*,end=100) s,v,vspc
            endif
C            write(fu6,*) s,v
         endif
         if (icode(iop).eq.1) then
            read (fu30,3005,end=100) ifrgrp(1,iop)
C            write(fu6,3005) ifrgrp(1,iop)
         elseif (icode(iop).eq.2) then
            write(fu6,6000) iop
            call errquit('drdy_restor: fatal error',911, UNKNOWN_ERR)
         else
            nend = 3*nratom(iop)
            nfreq = nf(iop)
            read (fu30,3005,end=100) (ifrgrp(i,iop),i=1,nfreq)
C            write(fu6,3005) (ifrgrp(i,iop),i=1,nfreq)
            read (fu30,3000,end=100) line
C            write(fu6,*) line
            nbar = (nend*(nend+1))/2
            read (fu30,*,end=100) (hgrp(i,iop),i=1,nbar)
C            write(fu6,*) (hgrp(i,iop),i=1,nbar)
C Diagonalize hessian
            do i = 1,nbar
               hess(i) = hgrp(i,iop)
            enddo
            ierr = 0
            call drdy_fdiag(nend,nfreq,ierr,0)
            if (ierr.eq.0) then
C Store frequencies in freqr array
               if (iop.eq.1 .or. iop.eq.3) then
                  do ifreq = 1,nfreq
                     freqr(ifreq,iop) = scr2(ifreq)
                  enddo
               else
                  nfreq1 = nf(iop-1)
                  do ifreq = 1,nfreq
                     freqr(ifreq+nfreq1,iop-1) = scr2(ifreq)
                  enddo
               endif
            endif
         endif
         jop = iop
         if ((iop.eq.1 .and. lgs(6).le.2) .or. 
     *        iop.eq.2. or.
     *       (iop.eq.3 .and. (lgs(6).eq.1.or.lgs(6).eq.3)) .or.
     *        iop.eq.4) then
            iop = iop + 1
         else
            iop = iop + 2
         endif
      enddo
      nfreq = nf(5)
C      write(fu6,*) ' ***** nf=',nf
      if (lgs(1).ne.0) then
c saddle point
         read (fu30,3000,end=101) line
C         write(fu6,*) line
            if (lgs(39).eq.0) then
               read (fu30,*,end=101) s,vspsv
            else
               read (fu30,*,end=101) s,vspsv,vspspsv
            endif
C         write(fu6,*) s,vspsv
C         write(fu6,*) ' ***** nfreq=',nfreq
         read (fu30,3005,end=101) (ifrsp(i),i=1,nfreq)
C         write(fu6,3005) (ifrsp(i),i=1,nfreq)
         read (fu30,3000,end=101) line
C         write(fu6,*) line
         do i = 1, n3
            read (fu30,*,end=101) (fspsv(j,i),j=1,i)
C            write(fu6,*) (fspsv(j,i),j=1,i)
         enddo
         do i = 2,n3
            do j = 1,i-1
               fspsv(i,j) = fspsv(j,i)
            enddo
         enddo
C Diagonalize hessian
         ii = 0
         do i = 1,n3
            do j = 1,i
               ii = ii + 1
               hess(ii) = fspsv(j,i)
            enddo
         enddo
         ierr = 0
         call drdy_fdiag(n3,nfreq+1,ierr,0)
         if (ierr.eq.0) then
C Store frequencies in freqr array
            do i = 1,nfreq+1
               freqr(i,5) = scr2(i)
            enddo
         endif
C
      endif
      lsad = .true.
      return
98    continue
      write(fu6,*) ' trouble reading title line from fu30 file'
      return
99    continue
      write(fu6,*) ' trouble reading lopt line from fu30 file'
      return
100   continue
      write(fu6,*) ' trouble reading from fu30 file for iop=',iop
      return
101   continue
      write(fu6,*) ' trouble reading saddle point from fu30 file'
      return
3000  Format (a80)
3005  Format (1x,i2,39i3)
6000  Format (' For iop=',i5,', icode=2, Morse diatomic, was chosen',
     *   /, 'This option is not supported, reset icode to 3 and try',
     *   ' again')
      End
