      integer function argos_cafe_numcon(idb,ndb,ia)
c
      implicit none
c
      integer ia,ndb
      integer idb(ndb,4)
c
      integer i,number
c
      number=0
c
      do 1 i=1,ndb
      if(idb(i,3).eq.1.and.(idb(i,1).eq.ia.or.idb(i,2).eq.ia))
     + number=number+1
    1 continue
c
      argos_cafe_numcon=number
c
      return
      end
c $Id$
