      logical function argos_prep_ring(lato,mato,nato,lang,mang,nang,
     + ldih,mdih,ndih)
c
      implicit none
c
      integer mato,nato,mang,nang,mdih,ndih
      integer lato(5,mato),lang(3,mang),ldih(4,mdih)
c
      integer i,j,k
c
      argos_prep_ring=.true.
c
      do 1 i=1,nato
      lato(4,i)=0
    1 continue
c
      do 2 i=1,nang
      do 3 j=1,ndih
      if((lang(1,i).eq.ldih(1,j).and.lang(3,i).eq.ldih(4,j)).or.
     + (lang(1,i).eq.ldih(4,j).and.lang(3,i).eq.ldih(1,j))) then
      if(lato(3,lang(2,i)).ne.3) goto 3
      lato(4,lang(2,i))=lato(4,lang(2,i))+1
      do 4 k=1,4
      if(lato(3,ldih(k,j)).ne.3) goto 3
      lato(4,ldih(k,j))=lato(4,ldih(k,j))+1
    4 continue
      endif
    3 continue
    2 continue
c
      do 6 i=1,ndih-1
      do 7 j=i+1,ndih
      if((ldih(1,i).eq.ldih(1,j).and.ldih(4,i).eq.ldih(4,j)).or.
     + (ldih(1,i).eq.ldih(4,j).and.ldih(4,i).eq.ldih(1,j))) then
      if(lato(3,ldih(2,i)).ne.3) goto 7
      if(lato(3,ldih(3,i)).ne.3) goto 7
      lato(4,ldih(2,i))=lato(4,ldih(2,i))+2
      lato(4,ldih(3,i))=lato(4,ldih(3,i))+2
      do 8 k=1,4
      if(lato(3,ldih(k,j)).ne.3) goto 7
      lato(4,ldih(k,j))=lato(4,ldih(k,j))+2
    8 continue
      endif
    7 continue
    6 continue
c
      return
      end
c $Id$
