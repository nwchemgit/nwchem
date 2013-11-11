C> \ingroup selci
C> @{
      subroutine selci_restin(iconf,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     &                  nroot,irange,nrange,iwpt,
     &                  noconf,ioconf,nintpo,nconmx,ncold)
*
* $Id$
*
      implicit real*8(a-h,o-z), integer(i-n)
#include "errquit.fh"
      dimension ioconf(nintpo,nconmx),irange(nrange),
     &          ept(nroot),enew(nroot),ptnorm(nroot),
     &          ept_mp(nroot),ptnorm_mp(nroot)
c
c     irest = unit for restart tape 'selci.r.?'
c            ? is 1 or 2
c     istde = standard error
c     istdo = standard output
c 
      parameter(irest = 69)
      parameter(istde=0,istdo=6)
      logical doesit
#include "cselcifiles.fh"
#include "global.fh"
      me = ga_nodeid()
      master = 0
      if (me.ne.master) goto 90000
c
c check existence of file
c
      ifile = 1
      doesit = .false.
      inquire(file=names(1),exist=doesit)
      if (doesit) goto 1
      inquire(file=names(2),exist=doesit)
      if (doesit) then
         ifile = 2
      else
         goto 90000
      endif
c      

 1    open(unit=irest,file=names(ifile),access='sequential',
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
         if (iwpt.lt.0 .or. iwpt.gt.3) call errquit('rest: iwpt',iwpt,
     &       UNKNOWN_ERR)
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
      do jj = 1,nloop
         read(irest)(ioconf(ii,jj),ii=1,nintpo)
      enddo
      close(unit=irest,status='keep')
      noconf = inocnf
c
      write(6,123) iconf, ncold
 123  format(/' Restarting at configuration ',i6,' of ',i6/)
c
90000 continue
      call ga_brdcst(50,iconf,mitob(1),master)
      return
c
90911 call errquit(' error opening restart file ',0, DISK_ERR)
91911 continue
c 
      write(istdo,*)' restart info is not correct. it will be ignored '
      write(istde,*)' restart info is not correct. it will be ignored '
      close(unit=irest,status='keep')
      goto 90000
      end
C> @}
