      logical function pre_nonbon(iatm,natm,lbnd,mbnd,ibnd,nbnd,
     + lang,mang,iang,nang,l3rd,m3rd,i3rd,n3rd,lnon,mnon,nnon)
c
c $Id$
c
      implicit none
c
      integer iatm,natm,mbnd,ibnd,nbnd,mang,iang,nang
      integer m3rd,i3rd,n3rd,mnon,nnon
      integer lbnd(4,mbnd),lang(5,mang),l3rd(2,m3rd),lnon(2,mnon)
c
      integer i,j,k
c
      do 1 i=iatm,natm-1
      do 2 j=i+1,natm
      do 3 k=ibnd,nbnd
      if(lbnd(1,k).eq.i.and.lbnd(2,k).eq.j) goto 2
    3 continue
      do 4 k=iang,nang
      if(lang(1,k).eq.i.and.lang(3,k).eq.j) goto 2
    4 continue
      do 5 k=i3rd,n3rd
      if(l3rd(1,k).eq.i.and.l3rd(2,k).eq.j) goto 2
    5 continue
      do 6 k=1,nnon
      if(lnon(1,k).eq.i.and.lnon(2,k).eq.j) goto 2
    6 continue
      nnon=nnon+1
      if(nnon.gt.mnon) call md_abort('increase mnon',9999)
      lnon(1,nnon)=i
      lnon(2,nnon)=j
    2 continue
    1 continue
c
      pre_nonbon=.true.
      return
      end
