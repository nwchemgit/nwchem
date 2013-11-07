C> \ingroup selci
C> @{
C>
      subroutine selci_hpp(rtdb, hpp,nbig,npvec,energy)
c
c     modified from sigma.F
c
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "util.fh"
#include "cbuff.fh"
#include "cselcifiles.fh"
      integer rtdb,nbig,npvec,iflham,iflhdg,nci
      double precision hd(npvec),hpp(nbig,nbig), hpptemp(npvec,npvec),
     &     temp(npvec,npvec), htemp(npvec,npvec), c0(npvec,npvec),
     &     energy(*)
      parameter(iflham = 1, iflhdg = 2)
c     
      call dfill(nbig*nbig,0.0d0,hpp,1)
      call dfill(npvec*npvec,0.0d0,hpptemp,1)
      call dfill(npvec*npvec,0.0d0,temp,1)
      call dfill(npvec*npvec,0.0d0,htemp,1)
      call dfill(npvec*npvec,0.0d0,c0,1)
      do i = 1, npvec
         do j = 1, npvec
            if(i.eq.npvec-j+1)c0(i,j)=1.0d0
         enddo
      enddo
c     
      write(*,*)"transformation matrix c0"
      call output(c0,1,npvec,1,npvec,npvec,npvec,1)
c     
c     get iflham from rtdb?
c     
c      if(.not.rtdb_get(rtdb,"selci:iflham",mt_int,1,iflham))
c     &     call errquit("don't know iflham,selci_hpp",0)
c      if(.not.rtdb_get(rtdb,"selci:iflhdg",mt_int,1,iflhdg))
c     &     call errquit("don't know iflhdg,selci_hpp",0)
      if(ga_nodeid().eq.0)then
         open(iflham,status='old',form='unformatted',file=hamil)
         open(iflhdg,status='old',form='unformatted',file=hamdg)
         read (iflham) nci
         if(nci.ne.npvec)
     &        call errquit("nci doesn't match npvec", nci, INPUT_ERR)
      endif
c     
      call selci_sread(iflhdg,hd,nci)
      do i = 1, npvec
         htemp(i,i)=hd(i)
      enddo
c     
 30   call selci_sread(iflham,rbuff,lenbuf)
      iptr = 0
      ipti = 0
c     
 40   i = ibuff(ipti+1)
c      write(6,*) ' sigma: i=',i
      if (i.eq.iendr) goto 30   
      if (i.eq.iendf) goto 1000 
      i = i + imagic            
      nj = ibuff(ipti+2)        
      ipti = ipti + 2
c
      do jj = 1,nj
         j = ibuff(ipti+jj) + imagic
         hij = rbuff(iptr+jj)
         if(i.gt.npvec.or.j.gt.npvec) 
     &        call errquit("p space dimension mismatch?", npvec,
     &       INPUT_ERR)
         htemp(i,j)=hij
         htemp(j,i)=htemp(i,j)
      enddo
      ipti = ipti + nj
      iptr = iptr + nj
      goto 40
c     
 1000 continue
c
      close(iflham,status='keep')
      close(iflhdg,status='keep')
c     
      write(*,*)"hamiltonian from iflham"
      call output(htemp,1,npvec,1,npvec,npvec,npvec,1)
c     
      call dgemm("n","n",npvec,npvec,npvec,1.0d0,htemp,npvec,c0,npvec,
     $     0.0d0,temp,npvec)
      write(*,*)"temp=htemp*c0"
      call output(temp,1,npvec,1,npvec,npvec,npvec,1)
      write(*,*)"energy=",energy(1)
      energy(1)=energy(1)+0.1
      do i=1,npvec
         call daxpy(npvec,-energy(1),c0(1,i),1,temp(1,i),1)
      enddo
      call dgemm("t","n",npvec,npvec,npvec,1.0d0,c0,npvec,temp,npvec,
     $     0.0d0,hpptemp,npvec)
      write(*,*)"hamiltonian after transformation"
      call output(hpptemp,1,npvec,1,npvec,npvec,npvec,1)
      do i=1,npvec
         call dcopy(npvec,hpptemp(1,i),1,hpp(1,i),1)
      enddo
      write(*,*)"hamiltonian after transformation, hpp"
      call output(hpp,1,npvec,1,npvec,nbig,nbig,1)
c     
c     what do I do here?
c
c
c      do 1010 iroot = 1, nroot
c         call ga_dgop (99, si(1, iroot), nci, '+')
c 1010 continue
c     
      end
C>
C> @}
c $Id$
