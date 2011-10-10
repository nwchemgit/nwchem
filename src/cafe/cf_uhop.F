      subroutine cf_uhop(lpbc,lpbcs,lseq,isl,isgm,isga,isq3,xs,nsaloc,
     + chg,lda,rda,uda)
c
      implicit none
c
#include "cf_common.fh"
#include "bitops.fh"
#include "util.fh"
#include "mafdecls.fh"
c
      logical lpbc,lpbcs
      real*8 xs(msa,3),chg(mqt,mqp,mset)
      integer nsaloc
c
      integer isl(msa,mis2)
      integer lda(16,*),isgm(msa),isga(msa),isq3(msa),lseq(mseq)
      real*8 rda(11,*),uda(4,*)
c
      real*8 q,rwx(3),rwy(3),rw,rwi,uhop(4),facthp,cutoff,uh(4)
      integer isa,jsa,iss,jss,ndxi
      integer ilda,idon,idong,idons,idonp,iacc,iaccg,iaccs,iaccp
      integer idown,iaown,idfr,idto,iafr,iato
      integer jsprev,i
      logical lskeep
c
      integer ndon,nacc
c
      facthp=1.0d0
      cutoff=1.0d3
      jsprev=0
      lskeep=.false.
c
c      write(*,5555) stime
c 5555 format('TIME ',f12.6)
      do 1 ilda=1,nlda
c
      uda(1,ilda)=0.0d0
      uda(2,ilda)=0.0d0
      uda(3,ilda)=0.0d0
      uda(4,ilda)=0.0d0
c
      uhop(1)=0.0d0
      uhop(2)=0.0d0
      uhop(3)=0.0d0
      uhop(4)=0.0d0
c
      idon=lda(2,ilda)
      idong=lda(1,ilda)
      idons=lda(3,ilda)
      idonp=lda(4,ilda)
      iacc=lda(7,ilda)
      iaccg=lda(6,ilda)
      iaccs=lda(8,ilda)
      iaccp=lda(9,ilda)
c
      idown=lda(11,ilda)
      idfr=lda(12,ilda)
      idto=lda(13,ilda)
      iaown=lda(14,ilda)
      iafr=lda(15,ilda)
      iato=lda(16,ilda)
c
      ndon=idto-idfr+1
      nacc=iato-iafr+1
c
c      write(*,'(i5,a,5i5)') me,' LDA d ',ilda,nlda,idown,idfr,idto
c
c     get donor and acceptor atom coordinates
c
      call sp_getda(idown,idfr,idto,isl,isga,isq3,xs,nsaloc)
c      write(*,'(i5,a,5i5)') me,' LDA a ',ilda,nlda,iaown,iafr,iato
      call sp_getda(iaown,iafr,iato,isl,isga,isq3,xs,nsaloc+ndon)
c
      do 2001 isa=1,ndon
c      write(*,3001) isa,isga(nsaloc+isa),isq3(nsaloc+isa),
c     + (xs(nsaloc+isa,jsa),jsa=1,3)
c 3001 format(' donor',3i5,3f12.6)
      isgm(nsaloc+isa)=idons
 2001 continue
      do 2002 isa=1,nacc
c      write(*,3002) isa,isga(nsaloc+ndon+isa),isq3(nsaloc+ndon+isa),
c     + (xs(nsaloc+ndon+isa,jsa),jsa=1,3)
c 3002 format(' acceptor',3i5,3f12.6)
      isgm(nsaloc+ndon+isa)=iaccs
 2002 continue
c
c      write(*,'(a,8i5)') 'DA ',idon,idong,idons,idonp,
c     + iacc,iaccg,iaccs,iaccp
c
c     loop over the atoms in the donor and acceptor residues
c
cx      do 2 isa=1,nsaloc
      do 2 isa=nsaloc+1,nsaloc+ndon+nacc+1
c
      if(isgm(isa).eq.idons.or.isgm(isa).eq.iaccs) then
c
c     loop over all atoms not in either donor or acceptor residues
c
      do 3 jsa=1,nsaloc
c
      iss=isgm(isa)
      jss=isgm(jsa)
c
c      write(*,'(4i5)') isa,jsa,iss,jss
c
      if(jsprev.ne.jss) then
      if(lskeep) then
      uda(1,ilda)=uda(1,ilda)+facthp*uhop(1)
      uda(2,ilda)=uda(2,ilda)+facthp*uhop(2)
      uda(3,ilda)=uda(3,ilda)+facthp*uhop(3)
      uda(4,ilda)=uda(4,ilda)+facthp*uhop(4)
      endif
      uhop(1)=0.0d0
      uhop(2)=0.0d0
      uhop(3)=0.0d0
      uhop(4)=0.0d0
      jsprev=jss
      lskeep=.false.
      endif
c
      if(isgm(jsa).ne.idons.and.isgm(jsa).ne.iaccs) then
c
c     find the periodic boundary offset from the midpoint
c
c      write(*,3323) ilda,isa,xs(isa,1),xs(isa,2),xs(isa,3),
c     + jsa,xs(jsa,1),xs(jsa,2),xs(jsa,3)
c 3323 format(2i5,3f12.6,/,5x,i5,3f12.6)
c
      if(isgm(isa).eq.idons) then
      rwy(1)=xs(jsa,1)-rda(6,ilda)
      rwy(2)=xs(jsa,2)-rda(7,ilda)
      rwy(3)=xs(jsa,3)-rda(8,ilda)
      else
      rwy(1)=xs(jsa,1)-rda(9,ilda)
      rwy(2)=xs(jsa,2)-rda(10,ilda)
      rwy(3)=xs(jsa,3)-rda(11,ilda)
      endif
c
c      write(*,3335) rwy
      if(lpbc.or.lpbcs) call cf_pbc(0,rwy,1,rwx,1,0,1,1)
c
c      write(*,3333) (xs(jsa,i),i=1,3)
c 3333 format('xs(jsa,1:3)',3f12.6)
c      write(*,3334) (rda(i,ilda),i=6,11)
c 3334 format(3x,'rda(6:8)',3f12.6,/,'  rda(9:11)',3f12.6)
c      write(*,3335) rwx
c 3335 format(8x,'rwx',3f12.6)
c
c     evaluate the distance between atoms isa and jsa
c
      if(isga(isa).eq.iaccg) then
      rwx(1)=-rwx(1)+xs(jsa,1)-rda(1,ilda)
      rwx(2)=-rwx(2)+xs(jsa,2)-rda(2,ilda)
      rwx(3)=-rwx(3)+xs(jsa,3)-rda(3,ilda)
c      write(*,3337) (rda(i,ilda),i=1,3)
c 3337 format(3x,'rda(1:3)',3f12.6)
      else
      rwx(1)=-rwx(1)+xs(jsa,1)-xs(isa,1)
      rwx(2)=-rwx(2)+xs(jsa,2)-xs(isa,2)
      rwx(3)=-rwx(3)+xs(jsa,3)-xs(isa,3)
c      write(*,3336) (xs(isa,i),i=1,3)
c 3336 format('xs(isa,1:3)',3f12.6)
      endif 
      rw=sqrt(rwx(1)**2+rwx(2)**2+rwx(3)**2)
      rwi=1.0d0/rw
c
      if(rw.le.cutoff) lskeep=.true.
c
c      write(*,'(3i5,f12.6)') ilda,isa,jsa,rw
c
c     donor energies
c
      uh(1)=0.0d0
      uh(2)=0.0d0
      uh(3)=0.0d0
      uh(4)=0.0d0
c
      if(iss.eq.idons) then
c
c     before hopping
c
      q=chg(isq3(isa),3,lseq(iss))*chg(isq3(jsa),3,lseq(jss))
      uh(1)=q*rwi
      uhop(1)=uhop(1)+q*rwi
c      write(*,6665) stime,1,isga(isa),chg(isq3(isa),3,lseq(iss)),
c     + isga(jsa),chg(isq3(jsa),3,lseq(jss)),uh(1)
c 6665 format(f12.6,'q',2i5,f12.6,i5,f12.6,e12.5)
c
c     after hopping - skip if donor site
c
      if(isga(isa).ne.idong) then
      ndxi=isa
      if(isga(isa).eq.idong) ndxi=isa+lda(5,ilda)
      if(isga(isa).eq.idong+lda(5,ilda)) ndxi=isa-lda(5,ilda)
      q=chg(isq3(ndxi),3,idonp)*chg(isq3(jsa),3,lseq(jss))
      uh(3)=q*rwi
      uhop(3)=uhop(3)+q*rwi
c      write(*,6665) stime,2,isga(isa),chg(isq3(ndxi),3,idonp),
c     + isga(jsa),chg(isq3(jsa),3,lseq(jss)),uh(3)
      endif
      endif
c
c     acceptor energies
c
      if(iss.eq.iaccs) then
c
c     before hopping - skip if acceptor site
c
      if(isga(isa).ne.iaccg) then
      q=chg(isq3(isa),3,lseq(iss))*chg(isq3(jsa),3,lseq(jss))
      uh(2)=q*rwi
      uhop(2)=uhop(2)+q*rwi
c      write(*,6665) stime,3,isga(isa),chg(isq3(isa),3,lseq(iss)),
c     + isga(jsa),chg(isq3(jsa),3,lseq(jss)),uh(2)
      endif
c
c     after hopping
c
      ndxi=isa
      if(isga(isa).eq.iaccg) ndxi=isa+lda(10,ilda)
      if(isga(isa).eq.iaccg+lda(10,ilda)) ndxi=isa-lda(10,ilda)
      q=chg(isq3(ndxi),3,iaccp)*chg(isq3(jsa),3,lseq(isgm(jsa)))
      uh(4)=q*rwi
      uhop(4)=uhop(4)+q*rwi
c      write(*,6665) stime,4,isga(isa),chg(isq3(ndxi),3,iaccp),
c     + isga(jsa),chg(isq3(jsa),3,lseq(isgm(jsa))),uh(4)
      endif
c
c      write(*,6666) isga(isa),iss,isga(jsa),jss,lskeep,rw,uh
c 6666 format(4i5,4x,l1,f12.6,4e12.5)
c
      endif
c
    3 continue
c
      endif
c
    2 continue
c
      if(lskeep) then
      uda(1,ilda)=uda(1,ilda)+facthp*uhop(1)
      uda(2,ilda)=uda(2,ilda)+facthp*uhop(2)
      uda(3,ilda)=uda(3,ilda)+facthp*uhop(3)
      uda(4,ilda)=uda(4,ilda)+facthp*uhop(4)
      endif
c
c      write(*,'(a,8f12.3)') 'uhop ',(uda(ndxi,ilda),ndxi=1,4)
c
c
    1 continue
c
      return
      end
c $Id$
