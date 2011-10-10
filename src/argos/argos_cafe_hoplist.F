      subroutine argos_cafe_hoplist(issgm,nsaloc,nhopl,lda,rda)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer issgm(msa),nsaloc,nhopl(np)
      integer lda(16,*)
      real*8 rda(11,*)
c
      integer i,j,ioff
c
      do 1 i=1,np
      nhopl(i)=0
    1 continue
      nhopl(me+1)=nlda
c      write(*,'(i5,a,i5)') me,' nlda= ',nlda
c
      call ga_igop(mcf_76,nhopl,np,'+')
c
      do 2 i=2,np
      nhopl(i)=nhopl(i)+nhopl(i-1)
    2 continue
c      write(*,'(i5,a,10i5)') me,' nhopl ',(nhopl(i),i=1,np)
c
c
      ioff=0
      if(me.gt.0) then
      ioff=nhopl(me)
      do 3 i=1,nlda
      do 4 j=1,16
      lda(j,i+ioff)=lda(j,i)
    4 continue
      do 5 j=1,8
      rda(j,i+ioff)=rda(j,i)
    5 continue
    3 continue
      do 6 i=1,nhopl(me)
      do 7 j=1,16
      lda(j,i)=0
    7 continue
      do 8 j=1,8
      rda(j,i)=0.0d0
    8 continue
    6 continue
      endif
      if(me.lt.np-1) then
      do 9 i=nhopl(me+1)+1,nhopl(np)
      do 10 j=1,16
      lda(j,i)=0
   10 continue
      do 11 j=1,8
      rda(j,i)=0.0d0
   11 continue
    9 continue
      endif
c
      nldat=nhopl(np)
      call ga_igop(mcf_77,lda,16*nldat,'+')
c
c      do 322 i=1,nldat
c      write(*,'(18i5)') me,i,(lda(j,i),j=1,16)
c  322 continue
c
      do 201 i=1,nldat
      do 202 j=1,nsaloc
      if(lda(3,i).eq.issgm(j)) then
      if(lda(12,i).eq.0) lda(12,i)=j
      lda(13,i)=j
      endif
      if(lda(8,i).eq.issgm(j)) then
      if(lda(15,i).eq.0) lda(15,i)=j
      lda(16,i)=j
      endif
  202 continue
  201 continue
c
      if(me.gt.0) then
      do 221 j=1,nldat
      do 222 i=1,11
      lda(i,j)=0
  222 continue
      lda(14,j)=0
  221 continue
      endif
c
      call ga_igop(mcf_77,lda,16*nldat,'+')
      call ga_dgop(mcf_78,rda,11*nldat,'+')
c
c      do 2222 i=1,nlda
c      write(*,'(18i5)') me,i,(lda(j,i),j=1,16)
c 2222 continue
c
      nlda=0
      do 12 i=1,nldat
      do 13 j=1,nlda
      if(lda(1,i).eq.lda(1,j).and.lda(6,i).eq.lda(6,j)) goto 12
   13 continue
      nlda=nlda+1
      do 14 j=1,16
      lda(j,nlda)=lda(j,i)
   14 continue
      do 15 j=1,8
      rda(j,nlda)=rda(j,i)
   15 continue
   12 continue
      nldat=nlda
c
c
      if(me.eq.0) then
      if(util_print('qhop',print_high)) then
      write(lfnhop,110) nldat
  110 format(/,'Number of donor-acceptor pairs is ',i5,/,
     + '---------- Donor --------- --------- Acceptor ----',
     + ' --- Donor ---  ---Acceptor---',/,
     + ' glob  loc  sgm prot  off glob  loc  sgm prot  off',
     + '  own from   to  own from   to',/)
      do 212 i=1,nldat
      write(lfnhop,111) (lda(j,i),j=1,16),(rda(j,i),j=1,3)
  111 format(16i5,3f12.6)
  212 continue
      endif
      endif
c
      return
      end
c $Id$
