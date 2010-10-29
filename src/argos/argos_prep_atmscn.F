      logical function argos_prep_atmscn(lfnsgm,filsgm,atmnam,num,link)
c
c $Id$
c
      implicit none
c
      integer lfnsgm,num,link
      character*255 filsgm
      character*4 atmnam
c
      character*80 card
      integer i,j,length,natm,ilnk,npar
      character*6 name
      real*8 version
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='old',err=9999)
c
      num=0
      link=0
c
    1 continue
      read(lfnsgm,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
      read(card,1001,err=9999) version
 1001 format(f12.6)
      read(lfnsgm,1002,err=9999) natm,npar
 1002 format(i5,25x,i5)
      do 2 i=1,npar
      read(lfnsgm,1000,err=9999,end=9999) card
    2 continue
      do 3 i=1,natm
      read(lfnsgm,1003) name,ilnk
 1003 format(5x,a6,10x,i5)
      if(name(1:4).eq.atmnam) then
      link=ilnk
      num=i
      close(unit=lfnsgm)
      argos_prep_atmscn=.true.
      return
      endif
      do 4 j=1,npar
      read(lfnsgm,1000,err=9999,end=9999) card
    4 continue
    3 continue
      close(unit=lfnsgm)
 9999 continue
      argos_prep_atmscn=.false.
      return
      end
      character*2 function argos_prep_atnam(number)
c
c $Id$
c
      implicit none
      integer number
c
      character*2 aname(0:105)
c
      data aname / '  ',
     + ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     + 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     + 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br',' R','Rb','Sr',' Y','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr','Rf','Ha' /
c
      argos_prep_atnam=aname(number)
c
      return
      end
