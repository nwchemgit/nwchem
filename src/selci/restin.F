      subroutine selci_restin(iconf,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     &                  nroot,irange,nrange,iwpt,
     &                  noconf,ioconf,nintpo,nconmx,ncold)
      implicit real*8(a-h,o-z), integer(i-n)
      dimension ioconf(nintpo,nconmx),irange(nrange),
     &          ept(nroot),enew(nroot),ptnorm(nroot),
     &          ept_mp(nroot),ptnorm_mp(nroot)
c
c irest = unit for restart tape 'selci.r.?'
c         ? is 1 or 2
c istde = standard error
c istdo = standard output
c 
      parameter(irest = 69)
*mdc*if unix
      parameter(istde=0,istdo=6)
*mdc*endif
      logical doesit
#include "cstool.fh"
      me = nodeid()
      master = 0
      if (me.ne.master) goto 90000
c
c check existence of file
c
      ifile = 1
      doesit = .false.
      inquire(file=names(1),exist=doesit)
      if (doesit) goto 00001
      inquire(file=names(2),exist=doesit)
      if (doesit) then
         ifile = 2
      else
*         write(istde,'(a,2(1x,a9))')
*     &    ' no restart file exists ',names
*         write(istdo,'(a,2(1x,a9))')
*     &    ' no restart file exists ',names
         goto 90000
      endif
c      
00001 continue
      open(unit=irest,file=names(ifile),access='sequential',
     &     form='unformatted',status='old',
     &     err=90911)
c
      read(irest)iconf,incold,inroot,inrang,inocnf,inntpo,incnmx,iniwpt
      nloop = inocnf
      if (incold.ne.ncold)   goto 91911
      if (inroot.ne.nroot)   goto 91911
      if (inrang.ne.nrange)  goto 91911
      if (inntpo.ne.nintpo)  goto 91911
      if (iniwpt.ne.iwpt)  then
         iwpt = iniwpt
         write(6,*)' switching pt selection procedure '
c 1 is MP 0 is EN.  other is dead!!
         if (iwpt.lt.0 .or. iwpt.gt.3) call errquit('rest: iwpt',iwpt)
      endif
      read(irest)(irange(ii),ii=1,nrange)
      read(irest)(ept(ii),ii=1,nroot)
      read(irest)(ept_mp(ii),ii=1,nroot)
      read(irest)(enew(ii),ii=1,nroot)
      read(irest)(ptnorm(ii),ii=1,nroot)
      read(irest)(ptnorm_mp(ii),ii=1,nroot)
c
      if (inocnf.gt.nconmx) goto 91911
c
      do 00100 jj = 1,nloop
         read(irest)(ioconf(ii,jj),ii=1,nintpo)
c         write(6,00101)'node:',nodeid(),' selci_restin(1):ioconf ',
c     &    ii,nloop,(ioconf(ii,jj),ii=1,nintpo)
00100 continue
c00101 format(1x,a,i3,a,i3,i5,22(2x,z8))
      close(unit=irest,status='keep')
      noconf = inocnf
c
      write(6,123) iconf, ncold
 123  format(/' Restarting at configuration ',i6,' of ',i6/)
*      write(6,*)'node:',nodeid(),' selci_restin(2):',
*     & iconf,ncold,nroot,nrange,
*     & noconf,nintpo,nconmx
*      write(6,*)' irange: ',irange
*      write(6,*)' ept:    ',ept
*      write(6,*)' enew:   ',enew
*      write(6,*)' ptnorm: ',ptnorm
c
90000 continue
      call brdcst(50,iconf,mitob(1),master)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
90911 call errquit(' error opening restart file ',0)
91911 continue
c      write(istdo,*)'iconf,inroot,inrang,inocnf,inntpo,incnmx,incold'
c      write(istdo,*)iconf,inroot,inrang,inocnf,inntpo,incnmx,incold
c      write(istdo,*)'nroot,nrange,noconf,nintpo,nconmx,ncold'
c      write(istdo,*)nroot,nrange,noconf,nintpo,nconmx,ncold
c
c      write(istde,*)'iconf,inroot,inrang,inocnf,inntpo,incnmx,incold'
c      write(istde,*)iconf,inroot,inrang,inocnf,inntpo,incnmx,incold
c      write(istde,*)'nroot,nrange,noconf,nintpo,nconmx,ncold'
c      write(istde,*)nroot,nrange,noconf,nintpo,nconmx,ncold
c
c
c 
      write(istdo,*)' restart info is not correct. it will be ignored '
      write(istde,*)' restart info is not correct. it will be ignored '
      close(unit=irest,status='keep')
      goto 90000
      end
