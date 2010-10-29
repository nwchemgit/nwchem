      logical function argos_prep_dih(lang,mang,nang,ldih,mdih,ndih)
c
c $Id$
c
      implicit none
c
      integer mang,mdih
      integer nang,ndih
      integer lang(3,mang),ldih(4,mdih)
      integer i,j,k,ltemp
c
c
c     determine the list of dihedrals
c     -------------------------------
c
      argos_prep_dih=.false.
c
      do 1 i=1,nang-1
      do 2 j=i+1,nang
      if(lang(2,i).eq.lang(1,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) return
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(3,j)
      else
      if(lang(2,i).eq.lang(3,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) return
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(1,j)
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(3,j)) then
      if(ndih.ge.mdih) return
      ndih=ndih+1
      ldih(1,ndih)=lang(1,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(1,j)) then
      if(ndih.ge.mdih) return
      ndih=ndih+1
      ldih(1,ndih)=lang(3,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      endif
      endif
      endif
      endif
    2 continue
    1 continue
c
      do 3 i=1,ndih
      if(ldih(2,i).gt.ldih(3,i)) then
      ltemp=ldih(2,i)
      ldih(2,i)=ldih(3,i)
      ldih(3,i)=ltemp
      ltemp=ldih(1,i)
      ldih(1,i)=ldih(4,i)
      ldih(4,i)=ltemp
      endif
    3 continue
c
      do 4 i=1,ndih-1
      do 5 j=i+1,ndih
      if(ldih(2,i).gt.ldih(2,j).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).gt.ldih(3,j)).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).eq.ldih(3,j).and.
     +  ldih(1,i).gt.ldih(1,j)).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).eq.ldih(3,j).and.
     +  ldih(1,i).eq.ldih(1,j).and.ldih(4,i).gt.ldih(4,j))) then
      do 6 k=1,4
      ltemp=ldih(k,i)
      ldih(k,i)=ldih(k,j)
      ldih(k,j)=ltemp
    6 continue
      endif
    5 continue
    4 continue
c
      argos_prep_dih=.true.
      return
      end
