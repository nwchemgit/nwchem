      logical function pre_atype(lfnout,lfnpar,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + latm,catm,matm,lbnd,mbnd,nbnd,
     + jlo,ilo,ihi,jhi,ltyp,mtyp,lring,aring,mring,
     + nring3,nring4,nring5,nring6,
     + latmt,matmt,natmt,lbndt,mbndt,nbndt)
c
c $Id: pre_atype.F,v 1.15 2001-11-06 06:15:20 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
c     in  : lfnpar            = atom types file number
c           filtyp            = atom types file name
c           latm(2,matm)      = atomic number
c                3            = number of bonds
c                4            = center type
c           catm(1,matm)      = atom name
c                2            = atom name (corrected)
c           matm              = dimension atom list
c           natm              = length atom list
c           lbnd(1:2,mbnd)    = bond indices
c           mbnd              = dimension bond list
c           nbnd              = length bond list
c           lring(1:6,mring)  = ring indices
c           aring(mring)      = logical true if aromatic ring
c
c     out : catm(3,matm)      = atom type
c           ltyp(1,mtyp)      = saturation
c                2            = aliphatic ring
c                3            = aromatic ring
c                4            = number of hydrogen neigbors
c                5            = number of neighbors
c                5+i          = index of i-th neighbor
c
      integer maty
      parameter (maty=1000)
c
      integer lfnpar,lfnout
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
      character*255 par_s,par_x,par_q,par_u,par_t,par_c
      character*255 filnam
      integer matm,matmt,mtyp,natmt
      integer latm(5,matm),latmt(5,matmt)
      character*6 catm(3,matm)
      integer mbnd,nbnd,mbndt,nbndt
      integer lbnd(2,mbnd),lbndt(2,mbndt)
      integer jlo,ilo,ihi,jhi
      integer ltyp(15,jlo:jlo+mtyp-1)
      integer mring,nring3,nring4,nring5,nring6
      integer lring(6,mring)
      logical aring(mring)
c
      integer i,j,k,kk,l,ll,m,ntype,lt,ld,nadd
      character*4 atype(maty)
      integer itype(20,maty)
      character*100 card
      integer iatnum,nhydr,nonh,irng,ifr,ito
      logical skipa(5),skipt(5),skipaa(5),skiptt(5)
c
c     setup typ array with latm data
c     ------------------------------
c
      do 101 i=1,mtyp
      do 102 j=1,5
      ltyp(j,i)=0
  102 continue
  101 continue
c
      do 1 i=jlo,jhi
      ltyp(1,i)=latm(2,i)
    1 continue
c
c     find index for each neighbor
c     ----------------------------
c
      do 2 i=1,nbnd
      if(lbnd(1,i).ge.jlo.and.lbnd(1,i).le.jhi) then
      ltyp(5,lbnd(1,i))=ltyp(5,lbnd(1,i))+1
      ltyp(5+ltyp(5,lbnd(1,i)),lbnd(1,i))=lbnd(2,i)
      if(latm(2,lbnd(2,i)).eq.1) ltyp(4,lbnd(1,i))=ltyp(4,lbnd(1,i))+1
      endif
      if(lbnd(2,i).ge.jlo.and.lbnd(2,i).le.jhi) then
      ltyp(5,lbnd(2,i))=ltyp(5,lbnd(2,i))+1
      ltyp(5+ltyp(5,lbnd(2,i)),lbnd(2,i))=lbnd(1,i)
      if(latm(2,lbnd(1,i)).eq.1) ltyp(4,lbnd(2,i))=ltyp(4,lbnd(2,i))+1
      endif
    2 continue
c
      ntype=jhi
      ifr=jlo
      ito=jhi
c 
      do 103 i=ifr,ito
      if(latm(1,i).gt.0) then
      ntype=ntype+1
      ltyp(5,i)=ltyp(5,i)+1
      ltyp(5+ltyp(5,i),i)=ntype
      ltyp(1,ntype)=latmt(2,latm(1,i))
      latm(2,ntype)=latmt(2,latm(1,i))
      nadd=0
      do 104 j=1,nbndt
      if(lbndt(1,j).eq.latm(1,i)) then
      nadd=nadd+1
      ltyp(5,ntype)=ltyp(5,ntype)+1
      ltyp(5+ltyp(5,ntype),ntype)=ntype+nadd
      ltyp(1,ntype+nadd)=latmt(2,lbndt(2,j))
      latm(2,ntype+nadd)=latmt(2,lbndt(2,j))
      if(latmt(2,lbndt(2,j)).eq.1) ltyp(4,ntype)=ltyp(4,ntype)+1
      elseif(lbndt(2,j).eq.latm(1,i)) then
      nadd=nadd+1
      ltyp(5,ntype)=ltyp(5,ntype)+1
      ltyp(5+ltyp(5,ntype),ntype)=ntype+nadd
      ltyp(1,ntype+nadd)=latmt(2,lbndt(1,j))
      latm(2,ntype+nadd)=latmt(2,lbndt(1,j))
      if(latmt(2,lbndt(1,j)).eq.1) ltyp(4,ntype)=ltyp(4,ntype)+1
      endif
  104 continue
      ntype=ntype+nadd
      if(ltyp(1,i).eq.6.and.ltyp(5,i).eq.3.and.latm(4,i).eq.0) then
      latm(4,i)=1
      if(latm(5,i).eq.4) then
      do 107 j=1,3
      do 108 k=1,3
      if(ltyp(1,ltyp(5+j,i)).eq.7.and.ltyp(1,ltyp(5+k,i)).eq.8)
     + latm(5,i)=2
  108 continue
  107 continue
      endif
      endif
      if(ltyp(1,i).eq.7.and.ltyp(5,i).eq.3.and.latm(4,i).eq.0) then
      do 105 j=1,3
      if(ltyp(1,ltyp(5+j,i)).eq.6.and.ltyp(5,ltyp(5+j,i)).eq.3) then
      do 106 k=1,3
      if(ltyp(1,ltyp(5+k,ltyp(5+j,i))).eq.8) then
      latm(4,i)=1
      if(latm(5,i).eq.3) latm(5,i)=1
      endif
  106 continue
      endif
  105 continue
      endif
      endif
  103 continue
c
c     order each neighbor list
c     ------------------------
c
      do 3 i=jlo,jhi
      do 4 j=1,ltyp(5,i)-1
      do 5 k=j+1,ltyp(5,i)
      if(ltyp(5+j,i).gt.ltyp(5+k,i)) then
      l=ltyp(5+j,i)
      ltyp(5+j,i)=ltyp(5+k,i)
      ltyp(5+k,i)=l
      endif
    5 continue
    4 continue
    3 continue
c
c     saturation
c     ----------
c
c     C sp2     : 3
c     C in C=O  : 2
c     N sp2     : 3
c     O in =O   : 2
c
      do 6 i=jlo,jhi
      if(latm(2,i).eq.6.and.ltyp(5,i).eq.3) then
      ltyp(1,i)=3
      do 7 j=1,ltyp(5,i)
      if(latm(2,ltyp(5+j,i)).eq.8.and.ltyp(5,ltyp(5+j,i)).eq.1) then
      ltyp(1,i)=2
      endif
    7 continue
      endif
      if(latm(2,i).eq.7.and.ltyp(5,i).eq.3) then
      ltyp(1,i)=3
      endif
      if(latm(2,i).eq.8.and.ltyp(5,i).eq.1) then
      ltyp(1,i)=2
      endif
    6 continue
c
c     C=C in aromatic 6-ring with 2 sp2 N : 2
c
      do 8 i=nring5+1,nring6
      if(aring(i)) then
      k=0
      do 9 j=1,6
      if(latm(2,i).eq.7) k=k+1
    9 continue
      if(k.eq.2) then
      do 10 j=1,nbnd
      k=0
      do 11 l=1,6
      if(latm(2,lring(l,i)).eq.6) then
      if(lbnd(1,j).eq.lring(l,i).or.lbnd(2,j).eq.lring(l,i)) k=k+1
      endif
   11 continue
      if(k.eq.2) then
      ltyp(1,lbnd(1,j))=2
      ltyp(1,lbnd(2,j))=2
      goto 8
      endif
   10 continue
      endif
      endif
    8 continue
c
c     ring types
c     ----------
c
      do 40 i=1,nring3
      do 41 j=1,3
      if(aring(i)) then
      ltyp(3,lring(j,i))=10*ltyp(3,lring(j,i))+3
      else
      ltyp(2,lring(j,i))=10*ltyp(2,lring(j,i))+3
      endif
   41 continue
   40 continue
      do 42 i=nring3+1,nring4
      do 43 j=1,4
      if(aring(i)) then
      ltyp(3,lring(j,i))=10*ltyp(3,lring(j,i))+4
      else
      ltyp(2,lring(j,i))=10*ltyp(2,lring(j,i))+4
      endif
   43 continue
   42 continue
      do 12 i=nring4+1,nring5
      do 13 j=1,5
      if(aring(i)) then
      ltyp(3,lring(j,i))=10*ltyp(3,lring(j,i))+5
      else
      ltyp(2,lring(j,i))=10*ltyp(2,lring(j,i))+5
      endif
   13 continue
   12 continue
      do 14 i=nring5+1,nring6
      do 15 j=1,6
      if(aring(i)) then
      ltyp(3,lring(j,i))=10*ltyp(3,lring(j,i))+6
      else
      ltyp(2,lring(j,i))=10*ltyp(2,lring(j,i))+6
      endif
   15 continue
   14 continue
c
      do 39 i=jlo,jhi
      if(ltyp(1,i).eq.3) then
      if(ltyp(3,i).eq.0) then
      ltyp(1,i)=2
      else
      ltyp(1,i)=0
      endif
      endif
   39 continue
c
c     read the atom type data
c     -----------------------
c
      ntype=0
      lt=index(par_s,' ')-1
      ld=index(dir_s,' ')-1
      l=ld+lt
      filnam=dir_s(1:ld)//par_s(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=901)
   18 continue
      read(lfnpar,1000,end=17,err=9999) card
      if(card(1:10).ne.'Atom types') goto 18
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
 2000 format(' Atom type definitions',t40,a)
      endif
   16 continue
      read(lfnpar,1000,end=17,err=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 16
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
 1001 format(a4,i7,i3,2i5,i3,i7,i3,3i7)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
 1002 format(27x,i7,i3,3i7)
      goto 16
   17 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
c
  901 continue
      lt=index(par_x,' ')-1
      ld=index(dir_x,' ')-1
      l=ld+lt
      filnam=dir_x(1:ld)//par_x(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=902)
  118 continue
      read(lfnpar,1000,end=117,err=9999) card
      if(card(1:10).ne.'Atom types') goto 118
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
      endif
  116 continue
      read(lfnpar,1000,end=117,err=9999) card
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 116
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
      goto 116
  117 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
c
  902 continue
      lt=index(par_q,' ')-1
      ld=index(dir_q,' ')-1
      l=ld+lt
      filnam=dir_q(1:ld)//par_q(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=903)
  218 continue
      read(lfnpar,1000,end=217,err=9999) card
      if(card(1:10).ne.'Atom types') goto 218
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
      endif
  216 continue
      read(lfnpar,1000,end=217,err=9999) card
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 216
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
      goto 216
  217 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
c
  903 continue
      lt=index(par_u,' ')-1
      ld=index(dir_u,' ')-1
      l=ld+lt
      filnam=dir_u(1:ld)//par_u(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=904)
  318 continue
      read(lfnpar,1000,end=317,err=9999) card
      if(card(1:10).ne.'Atom types') goto 318
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
      endif
  316 continue
      read(lfnpar,1000,end=317,err=9999) card
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 316
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
      goto 316
  317 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
c
  904 continue
      lt=index(par_t,' ')-1
      ld=index(dir_t,' ')-1
      l=ld+lt
      filnam=dir_t(1:ld)//par_t(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=905)
  418 continue
      read(lfnpar,1000,end=417,err=9999) card
      if(card(1:10).ne.'Atom types') goto 418
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
      endif
  416 continue
      read(lfnpar,1000,end=417,err=9999) card
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 416
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
      goto 416
  417 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
c
  905 continue
      lt=index(par_c,' ')-1
      ld=index(dir_c,' ')-1
      l=ld+lt
      filnam=dir_c(1:ld)//par_c(1:lt)
      open(unit=lfnpar,file=filnam(1:l),status='old',form='formatted',
     + err=906)
  518 continue
      read(lfnpar,1000,end=517,err=9999) card
      if(card(1:10).ne.'Atom types') goto 518
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:l)
      endif
  516 continue
      read(lfnpar,1000,end=517,err=9999) card
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 516
      ntype=ntype+1
      if(ntype.gt.maty) call errquit('increase maty',9999)
      read(card,1001) atype(ntype),(itype(i,ntype),i=1,10)
      read(lfnpar,1002) (itype(i,ntype),i=11,15)
      read(lfnpar,1002) (itype(i,ntype),i=16,20)
      goto 516
  517 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:l),' closed' 
      endif
  906 continue
      if(util_print('files',print_medium)) then
      write(lfnout,2001)
 2001 format(' ')
      endif
c
c     determine the atom types
c     ------------------------
c
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(i5,a4,5i5)')
     +  (i,catm(2,i),(ltyp(j,i),j=1,5),i=ilo,ihi)
      endif
c
      do 38 i=ilo,ihi
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(a,a4,6i5)') 'Atom ',catm(2,i),latm(2,i),
     + (ltyp(j,i),j=1,5)
      endif
      if(latm(2,i).gt.0) then
      do 19 j=1,ntype
      iatnum=itype(1,j)
      irng=0
      if(iatnum.ge.60000) then
      irng=6
      iatnum=iatnum-60000
      endif
      if(iatnum.ge.50000) then
      irng=5
      iatnum=iatnum-50000
      endif
      nonh=0
      do 20 k=1,4
      if(iatnum.ge.1000)  then
      nonh=nonh+1
      iatnum=iatnum-1000
      endif
   20 continue
      nhydr=0
      do 21 k=1,4
      if(iatnum.ge.200)  then
      nhydr=nhydr+1
      iatnum=iatnum-200
      endif
      if(nhydr.eq.4) nhydr=5
   21 continue
c
c     itype(1,*) match atom number
c           2    match explicit saturation
c           3    match aliphatic ring
c           4    match aromatic ring
c           5    match explicit number of neighbors
c
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(a,a4,6i5,1x,6l1)') 'Type ',atype(j)(1:4),iatnum,
     + itype(2,j),itype(3,j),itype(4,j),nhydr,itype(5,j),
     + iatnum.eq.latm(2,i),
     + (itype(2,j).eq.0.or.itype(2,j).eq.ltyp(1,i)),
     + (itype(3,j).eq.0.or.itype(3,j).eq.ltyp(2,i).or.
     +   (itype(3,j).eq.1.and.ltyp(2,i).gt.0)),
     + (itype(4,j).eq.0.or.itype(4,j).eq.ltyp(3,i).or.
     +   (itype(4,j).eq.1.and.ltyp(3,i).gt.0)),
     + (nhydr.eq.0.or.(nhydr.eq.5.and.ltyp(4,i).eq.0).or.
     +   nhydr.eq.ltyp(4,i)),
     + (itype(5,j).eq.0.or.itype(5,j).eq.ltyp(5,i))
      endif
c
      if(iatnum.eq.latm(2,i).and.
     + (nhydr.eq.0.or.(nhydr.eq.5.and.ltyp(4,i).eq.0).or.
     +   nhydr.eq.ltyp(4,i)).and.
     + (itype(2,j).eq.0.or.itype(2,j).eq.ltyp(1,i)).and.
     + (itype(3,j).eq.0.or.itype(3,j).eq.ltyp(2,i).or.
     +   (itype(3,j).eq.1.and.ltyp(2,i).gt.0)).and.
     + (itype(4,j).eq.0.or.itype(4,j).eq.ltyp(3,i).or.
     +   (itype(4,j).eq.1.and.ltyp(3,i).gt.0)).and.
     + (itype(5,j).eq.0.or.itype(5,j).eq.ltyp(5,i))) then
c
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(10x,a)') 'try'
      endif
c
c     initialize skip vectors for neighbor atoms : skipa(1:neighbors)
c     ----------------------- for atom type cond : skipt(1:3)
c
      do 22 k=1,5
      skipa(k)=.false.
   22 continue
      do 23 l=1,3
      skipt(l)=.false.
      if(itype(1+l*5,j).eq.0.and.itype(2+l*5,j).eq.0) skipt(l)=.true.
   23 continue
c
c     loop over neighbor conditions for the current atom type : l (1:3)
c     -------------------------------------------------------
c
      do 24 l=1,3
      if(.not.skipt(l)) then
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(5x,a,i5)') 'neighbor condition ',itype(1+l*5,j)
      endif
c
c     loop over neighbor atoms for the current atom : k (1:neighbors)
c     ---------------------------------------------
c
      do 25 k=1,ltyp(5,i)
      if(.not.skipa(k).and..not.skipt(l)) then
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(5x,a,i5)') 'neighbor atom      ',
     + latm(2,ltyp(5+k,i))
      endif
      iatnum=itype(1+l*5,j)
      irng=0
      if(iatnum.ge.60000) then
      irng=6
      iatnum=iatnum-60000
      endif
      if(iatnum.ge.50000) then
      irng=5
      iatnum=iatnum-50000
      endif
      nonh=0
      do 26 m=1,4
      if(iatnum.ge.1000)  then
      nonh=nonh+1
      iatnum=iatnum-1000
      endif
   26 continue
      nhydr=0
      do 27 m=1,4
      if(iatnum.ge.200)  then
      nhydr=nhydr+1
      iatnum=iatnum-200
      endif
      if(nhydr.eq.4) nhydr=5
   27 continue
c
c     check atomic number of neighbor k
c     ---------------------------------
c
      if(iatnum.gt.0.and.iatnum.ne.latm(2,ltyp(5+k,i))) goto 25
c
c     check bonded hydrogens to neighbor k
c     ------------------------------------
c
      if(nhydr.gt.0) then
      if(nhydr.eq.5) then
      if(ltyp(4,ltyp(5+k,i)).gt.0) goto 25
      else
      if(nhydr.ne.ltyp(4,ltyp(5+k,i))) goto 25
      endif
      endif
c
c     check bonded nonh-hydrogens to neighbor k
c     -----------------------------------------
c
      if(nonh.gt.0) then
      if(nonh.ne.ltyp(5,ltyp(5+k,i))-ltyp(4,ltyp(5+k,i))) goto 25
      endif
c
c     check number of neighbors
c     -------------------------
c
      if(itype(2+l*5,j).gt.0.and.
     + itype(2+l*5,j).ne.ltyp(5,ltyp(5+k,i))) goto 25
c
c     check if neighbor in ring
c     -------------------------
c
      if(irng.gt.0) then
      if(irng.ne.ltyp(3,ltyp(5+k,i))) goto 25
      endif
c
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(5x,a,2i5)') 'neighbor accepted ',
     + itype(1+l*5,j),latm(2,ltyp(5+k,i))
      endif
c
c     set skip vectors
c     ----------------
c
      do 28 kk=1,ltyp(5,ltyp(5+k,i))
      skipaa(kk)=ltyp(5+kk,ltyp(5+k,i)).eq.i
   28 continue
      do 29 ll=1,3
      skiptt(ll)=itype(2+l*5+ll,j).eq.0
   29 continue
c
c     loop over the neighbor of neighbor conditions
c     ---------------------------------------------
c
      do 30 ll=1,3
      if(.not.skiptt(ll)) then
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(10x,a,i5)') 'neighbor condition ',
     + itype(2+l*5+ll,j)
      endif
c
c     loop over neighbor of neighbor atoms
c     ------------------------------------
c
      do 31 kk=1,ltyp(5,ltyp(5+k,i))
      if(.not.skipaa(kk).and..not.skiptt(ll)) then
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(10x,a,i5)') 'neighbor atom      ',
     + latm(2,ltyp(5+kk,ltyp(5+k,i)))
      endif
c
      iatnum=itype(2+l*5+ll,j)
      irng=0
      if(iatnum.ge.60000) then
      irng=6
      iatnum=iatnum-60000
      endif
      if(iatnum.ge.50000) then
      irng=5
      iatnum=iatnum-50000
      endif
      nonh=0
      do 32 m=1,4
      if(iatnum.ge.1000)  then
      nonh=nonh+1
      iatnum=iatnum-1000
      endif
   32 continue
      nhydr=0
      do 33 m=1,4
      if(iatnum.ge.200) then
      nhydr=nhydr+1
      iatnum=iatnum-200
      endif
      if(nhydr.eq.4) nhydr=5
   33 continue
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(30x,a,3i5)') 'condition ia nn nh ',
     + iatnum,nonh,nhydr
      write(lfnout,'(30x,a,3i5)') 'found     ia nn nh ',
     + latm(2,ltyp(5+kk,ltyp(5+k,i))),
     + ltyp(5,ltyp(5+kk,ltyp(5+k,i)))-ltyp(4,ltyp(5+kk,ltyp(5+k,i))),
     + ltyp(4,ltyp(5+kk,ltyp(5+k,i)))
      endif
c
c     check atomic number of neighbor of neighbor k
c     ---------------------------------------------
c
      if(iatnum.gt.0.and.iatnum.ne.latm(2,ltyp(5+kk,ltyp(5+k,i))))
     +  goto 31
c
c     check bonded hydrogens to neighbor k
c     ------------------------------------
c
      if(nhydr.gt.0) then
      if(nhydr.eq.5) then
      if(ltyp(4,ltyp(5+kk,ltyp(5+k,i))).gt.0) goto 31
      else
      if(nhydr.ne.ltyp(4,ltyp(5+kk,ltyp(5+k,i)))) goto 31
      endif
      endif
c
c     check bonded nonh-hydrogens to neighbor k
c     -----------------------------------------
c
      if(nonh.gt.0) then
      if(nonh.ne.ltyp(5,ltyp(5+kk,ltyp(5+k,i)))-
     + ltyp(4,ltyp(5+kk,ltyp(5+k,i)))) goto 31
      endif
c
c     check if neighbor of neighbor in ring
c     -------------------------------------
c
      if(irng.gt.0) then
      if(irng.ne.ltyp(3,ltyp(5+kk,ltyp(5+k,i)))) goto 31
      endif
c
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(10x,a,4i5)') 'neighbor accepted ',
     + itype(2+l*5+ll,j),latm(2,ltyp(5+kk,ltyp(5+k,i))),kk,ll
      endif
c
c     neighbor of neighbor condition satisfied
c     ----------------------------------------
c
      skipaa(kk)=.true.
      skiptt(ll)=.true.
      endif
   31 continue
      endif
   30 continue
c
c     test if all neighbor of neighbor conditions are satisfied
c     ---------------------------------------------------------
c
      do 34 ll=1,3
      if(.not.skiptt(ll)) goto 25
   34 continue
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(20x,a)') 'neighbors of neighbor accepted'
      endif
c
c     neighbor condition satisfied
c     ----------------------------
c
      skipa(k)=.true.
      skipt(l)=.true.
      endif
   25 continue
      endif
   24 continue
c
c     test if all neighbor conditions are satisfied
c     ---------------------------------------------
c
      do 35 l=1,3
      if(.not.skipt(l)) goto 19
   35 continue
      if(util_print('atomtypes',print_debug)) then
      write(lfnout,'(30x,a)') 'neighbor accepted'
      endif
c
      catm(3,i)(1:4)=atype(j)
      if(util_print('atomtypes',print_debug)) then
       write(lfnout,'(a,a,a,a)') ' Accepted for atom ',catm(1,i),
     + ' type ',catm(3,i)
      endif
c
      endif
   19 continue
      endif
   38 continue
c
      pre_atype=.true.
      return
c
 9999 pre_atype=.false.
      return
      end
