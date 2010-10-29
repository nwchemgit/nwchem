      logical function pre_multip(lseq,mseq,nseq,latm,matm,natm,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp)
c
c
c $Id$
c
      implicit none
c
      integer mseq,nseq
      integer lseq(6,mseq)
      integer matm,natm
      integer latm(11,matm)
      integer mang,nang
      integer lang(5,mang)
      integer mdih,ndih
      integer ldih(6,mdih)
      integer mimp,nimp
      integer limp(6,mimp)
c
      integer i
c
c     check number of segments involved in angles
c     -------------------------------------------
c
      do 76 i=1,nang
      if(lseq(1,latm(5,lang(1,i))).ne.lseq(1,latm(5,lang(2,i))).and.
     + lseq(1,latm(5,lang(1,i))).ne.lseq(1,latm(5,lang(3,i))).and.
     + lseq(1,latm(5,lang(2,i))).ne.lseq(1,latm(5,lang(3,i)))) then
      lang(5,i)=-1
      else
      lang(5,i)=0
      endif
   76 continue
c
c     check number of segments involved in dihedrals
c     ----------------------------------------------
c
      do 77 i=1,ndih
      ldih(6,i)=0
      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(2,i))).and.
     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(3,i)))) then
      ldih(6,i)=-1
      else
      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(2,i))).and.
     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
      ldih(6,i)=-1
      else
      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
     + lseq(1,latm(5,ldih(3,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
      ldih(6,i)=-1
      else
      if(lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
     + lseq(1,latm(5,ldih(3,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
      ldih(6,i)=-1
      endif
      endif
      endif
      endif
   77 continue
c
c     check number of segments involved in improper dihedrals
c     -------------------------------------------------------
c
      do 78 i=1,nimp
      limp(6,i)=0
      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(2,i))).and.
     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(3,i))).and.
     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(3,i)))) then
      limp(6,i)=-1
      else
      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(2,i))).and.
     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(4,i))).and.
     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(4,i)))) then
      limp(6,i)=-1
      else
      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(3,i))).and.
     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(4,i))).and.
     + lseq(1,latm(5,limp(3,i))).ne.lseq(1,latm(5,limp(4,i)))) then
      limp(6,i)=-1
      else
      if(lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(3,i))).and.
     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(4,i))).and.
     + lseq(1,latm(5,limp(3,i))).ne.lseq(1,latm(5,limp(4,i)))) then
      limp(6,i)=-1
      endif
      endif
      endif
      endif
   78 continue
c
      pre_multip=.true.
      return
      end
