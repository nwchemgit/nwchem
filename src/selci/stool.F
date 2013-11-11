C> \ingroup selci
C> @{
      subroutine selci_stool(iconf,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     &     nroot,irange,nrange,iwpt,
     &     noconfi,ioconf,nintpo,nconmx,ncold)
*
* $Id$
*
      implicit real*8 (a-h,o-z), integer(i-n)
#include "errquit.fh"
c     
c     collect info from the other processes for restart
c     
      dimension ioconf(nintpo,nconmx),irange(nrange),
     &     ept(nroot),enew(nroot),ptnorm(nroot),
     &     ept_mp(nroot),ptnorm_mp(nroot)
c     
c     spare arrays needed
c     
      dimension eptn(50), enewn(50), ptnormn(50), irangen(21)
      dimension eptn_mp(50), ptnormn_mp(50)
      integer pmap(0:1023)      ! 1024 max process count
c     
c     irest = unit for restart tape 'selci.r.?'
c     ? is 1 or 2
c     istde = standard error
c     istdo = standard output
c     
      parameter(irest = 69)
      parameter(istde=0,istdo=6)
#include "cselcifiles.fh"
#include "global.fh"
      integer icall
      data icall/0/
c     
      its = mtime()
      if (nroot.gt.50) call errquit('STOOL: hard dimension fail',
     $     nroot, INPUT_ERR)
      me = ga_nodeid()
      master = 0
      if (ga_nnodes() .gt. 1024) call errquit('ptidy: too many procs',
     $     0, INPUT_ERR)
      call ga_list_nodeid(pmap, ga_nnodes()) ! Map to msg passing layer
*     write(6,*) ga_nodeid(),' processor map ', 
*     $     (pmap(i),i=0,ga_nnodes()-1)
      node0 = pmap(0)
c     
c     check existence of file
c     The following two lines are to take care of compiler warnings.
c     
      ifile  = 1
      idfile = 2

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
c     do not update vector length
c     if a crash occurs then node 0 (master) will 
c     hold all results of previous work.  
c     
      noconf = noconfi
      isync = 1
      if (me.ne.master) then
         call snd(21, (noconf-ncold), mitob(1), node0, isync)
         call snd(22, ioconf(1,ncold+1),mitob(nintpo*(noconf-ncold)),
     &        node0,isync)
         call snd(23, ept, mdtob(nroot), node0, isync)
         call snd(24, enew, mdtob(nroot), node0, isync)
         call snd(25, ptnorm, mdtob(nroot), node0, isync)
         call snd(26, irange, mitob(21), node0, isync)
         call snd(27, ept_mp, mdtob(nroot), node0, isync)
         call snd(28, ptnorm_mp, mdtob(nroot), node0, isync)
         call dfill(nroot,0.0d0,ept,1)
         call dfill(nroot,0.0d0,ept_mp,1)
         call dfill(nroot,0.0d0,enew,1)
         call dfill(nroot,0.0d0,ptnorm,1)
         call dfill(nroot,0.0d0,ptnorm_mp,1)
         call ifill(21,0,irange,1)
      else
         ipt = noconf+1
         do iislave=1,ga_nnodes()-1
            islave = pmap(iislave)
            call rcv(21, new, mitob(1), len, islave, node, isync)
            noconf = noconf + new
            if (noconf.le.nconmx) then
               call rcv(22, ioconf(1,ipt), mitob(nintpo*new), len,
     $              islave, node, isync)
               ipt = ipt + new
            else
               call errquit(' stool: noconf exceeded nconmx ',
     &              noconf, INPUT_ERR)
            endif
            call rcv(23, eptn, mdtob(nroot), len, islave, node,
     &           isync)
            call rcv(24, enewn, mdtob(nroot), len, islave, node,
     &           isync)
            call rcv(25, ptnormn, mdtob(nroot), len, islave, node,
     &           isync)
            call rcv(26, irangen, mitob(nrange), len, islave, node,
     &           isync)
            call rcv(27, eptn_mp, mdtob(nroot), len, islave, node,
     &           isync)
            call rcv(28, ptnormn_mp, mdtob(nroot), len, islave, node,
     &           isync)
            do i = 1,nroot
               ept(i) = ept(i) + eptn(i)
               ept_mp(i) = ept_mp(i) + eptn_mp(i)
               enew(i) = enew(i) + enewn(i)
               ptnorm(i) = ptnorm(i) + ptnormn(i)
               ptnorm_mp(i) = ptnorm_mp(i) + ptnormn_mp(i)
            enddo
            do i = 1,nrange
               irange(i) = irange(i) + irangen(i)
            enddo
         enddo
         open(unit=irest,file=names(ifile),access='sequential',
     &        form='unformatted',status='unknown',
     &        err=90911)
         close(unit=irest,status='delete')
         open(unit=irest,file=names(ifile),access='sequential',
     &        form='unformatted',status='unknown',
     &        err=90911)
         rewind irest
c     
         write(irest)iconf,ncold,nroot,nrange,noconf,nintpo,nconmx,
     $        iwpt
         write(irest)(irange(ii),ii=1,nrange)
         write(irest)(ept(ii),ii=1,nroot)
         write(irest)(ept_mp(ii),ii=1,nroot)
         write(irest)(enew(ii),ii=1,nroot)
         write(irest)(ptnorm(ii),ii=1,nroot)
         write(irest)(ptnorm_mp(ii),ii=1,nroot)
         nloop = noconf
         do jj = 1,nloop
            write(irest)(ioconf(ii,jj),ii=1,nintpo)
         enddo
         write(istdo,'(a,i6,a,i6)')
     &        ' checkpointed at iconf = ',iconf,' of ',ncold
c     
c     close restart tape
c     
         close(unit=irest,status='keep')
         if(me.eq.0)call util_flush(istdo)
c     
c     delete other restart tape
c     
         open(unit=irest,file=names(idfile),access='sequential',
     &        form='unformatted',status='unknown',
     &        err=90912)
         close(unit=irest,status='delete')
         itf = mtime()
      endif
c
      call ga_sync()
c     
      return
c     
90911 call errquit('stool: error opening restart file',0, DISK_ERR)
90912 call errquit('stool: error opening restart file',0, DISK_ERR)
c     
      end
C> @}
