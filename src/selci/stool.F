      subroutine selci_stool(iconf,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     & nroot,irange,nrange,iwpt,
     & noconfi,ioconf,nintpo,nconmx,ncold)
      implicit real*8 (a-h,o-z), integer(i-n)
c
c     collect info from the other processes for restart
c
      dimension ioconf(nintpo,nconmx),irange(nrange),
     & ept(nroot),enew(nroot),ptnorm(nroot),
     & ept_mp(nroot),ptnorm_mp(nroot)
*mdc*if intel
*      character*9 zdate
*      character*8 ztime
*mdc*endif
c
c spare arrays needed
c
      dimension eptn(50), enewn(50), ptnormn(50), irangen(21)
      dimension eptn_mp(50), ptnormn_mp(50)
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
#include "cstool.fh"
      integer icall
      data icall/0/
c      data names /'selci.r.1','selci.r.2'/
c
      its = mtime()
      if (nroot.gt.50) call errquit('STOOL: hard dimension fail',nroot)
      me = nodeid()
      master = 0
c      
c
c check existence of file
c
      if (me.eq.master) then
         icall = icall + 1
         ifile  = 1
         idfile = 2
         if(mod(icall,2).eq.0) then
            ifile  = 2
            idfile = 1
         endif
      endif
c
c do not update vector length
c if a crash occurs then node 0 (master) will 
c hold all results of previous work.  
c
      noconf = noconfi
*mdc*if unix
      call flush(istdo)
*mdc*endif
      isync = 1
      if (me.ne.master) then
         call snd(21, (noconf-ncold), mitob(1), 0, isync)
         call snd(22, ioconf(1,ncold+1),mitob(nintpo*(noconf-ncold)),
     &    0,isync)
         call snd(23, ept, mdtob(nroot), 0, isync)
         call snd(24, enew, mdtob(nroot), 0, isync)
         call snd(25, ptnorm, mdtob(nroot), 0, isync)
         call snd(26, irange, mitob(21), 0, isync)
         call snd(27, ept_mp, mdtob(nroot), 0, isync)
         call snd(28, ptnorm_mp, mdtob(nroot), 0, isync)
         call dfill(nroot,0.0d0,ept,1)
         call dfill(nroot,0.0d0,ept_mp,1)
         call dfill(nroot,0.0d0,enew,1)
         call dfill(nroot,0.0d0,ptnorm,1)
         call dfill(nroot,0.0d0,ptnorm_mp,1)
         call ifill(21,0,irange,1)
         goto 90000
      else
         ipt = noconf+1
         do 00100 islave=1,nnodes()-1
            call rcv(21, new, mitob(1), len, islave, node, isync)
            noconf = noconf + new
            if (noconf.le.nconmx) then
               call rcv(22, ioconf(1,ipt), mitob(nintpo*new), len,
     $          islave, node, isync)
               ipt = ipt + new
            else
               call errquit(' stool: noconf exceeded nconmx ',
     &              noconf)
            endif
            call rcv(23, eptn, mdtob(nroot), len, islave, node,
     &       isync)
            call rcv(24, enewn, mdtob(nroot), len, islave, node,
     &       isync)
            call rcv(25, ptnormn, mdtob(nroot), len, islave, node,
     &       isync)
            call rcv(26, irangen, mitob(nrange), len, islave, node,
     &       isync)
            call rcv(27, eptn_mp, mdtob(nroot), len, islave, node,
     &       isync)
            call rcv(28, ptnormn_mp, mdtob(nroot), len, islave, node,
     &       isync)
            do 00200 i = 1,nroot
               ept(i) = ept(i) + eptn(i)
               ept_mp(i) = ept_mp(i) + eptn_mp(i)
               enew(i) = enew(i) + enewn(i)
               ptnorm(i) = ptnorm(i) + ptnormn(i)
               ptnorm_mp(i) = ptnorm_mp(i) + ptnormn_mp(i)
00200       continue
            do 00300 i = 1,nrange
               irange(i) = irange(i) + irangen(i)
00300       continue
00100    continue
      endif
*      write(6,*)' open unknown ',names(ifile)
      open(unit=irest,file=names(ifile),access='sequential',
     & form='unformatted',status='unknown',
     & err=90911)
      close(unit=irest,status='delete')
*      write(6,*)' open new ',names(ifile)
      open(unit=irest,file=names(ifile),access='sequential',
     & form='unformatted',status='unknown',
     & err=90911)
c special rewind just incase compiler is dead about open
      rewind irest
c
      write(irest)iconf,ncold,nroot,nrange,noconf,nintpo,nconmx,iwpt
      write(irest)(irange(ii),ii=1,nrange)
      write(irest)(ept(ii),ii=1,nroot)
      write(irest)(ept_mp(ii),ii=1,nroot)
      write(irest)(enew(ii),ii=1,nroot)
      write(irest)(ptnorm(ii),ii=1,nroot)
      write(irest)(ptnorm_mp(ii),ii=1,nroot)
*      write(6,*)' stool:',iconf,ncold,nroot,nrange,noconf,
*     & nintpo,nconmx,iwpt
c      write(6,*)' irange:    ',irange
*      write(6,*)' enew:      ',enew
*      write(6,*)' ept:       ',ept
*      write(6,*)' ptnorm:    ',ptnorm
*      write(6,*)' ept_mp:    ',ept_mp
*      write(6,*)' ptnorm_mp: ',ptnorm_mp
      nloop = noconf
      do 00400 jj = 1,nloop
         write(irest)(ioconf(ii,jj),ii=1,nintpo)
c         write(6,00401)'node:',nodeid(),' selci_stool(1):ioconf ',
c     &    ii,nloop,(ioconf(ii,jj),ii=1,nintpo)
00400 continue
c00401 format(1x,a,i3,a,i3,i5,22(2x,z8))
*mdc*if intel
*      call date(zdate)
*      call time(ztime)
*      write(istdo,'(a,i6,a,i6,1x,a,1x,a)')
*     & ' checkpointed at iconf = ',iconf,' | ',ncold,zdate,ztime
*mdc*else
      write(istdo,'(a,i6,a,i6)')
     & ' checkpointed at iconf = ',iconf,' of ',ncold
*mdc*endif
c
c close restart tape
c
      close(unit=irest,status='keep')
*mdc*if unix
      call flush(istdo)
*mdc*endif
c
c delete other restart tape
c
**      write(6,*)' open other and delete ',names(idfile)
      open(unit=irest,file=names(idfile),access='sequential',
     & form='unformatted',status='unknown',
     & err=90912)
      close(unit=irest,status='delete')
      itf = mtime()
*      write(istdo,*)' time to checkpoint: ',((itf-its)/100),' seconds'
90000 continue
      call synch(88)
c
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
90911 call errquit('stool: error opening restart file',0)
90912 call errquit('stool: error opening restart file',0)
c
      end
