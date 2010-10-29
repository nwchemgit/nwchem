c
c $Id$
c

      SUBROUTINE list_excld(ntype)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_elst.inc'

      integer i,j,iexcl,nelist,exclpair,nexcl
      integer ntype,eatm

      dimension iexcl(mxtype),exclpair(mxtype,mxtype)

      do i=1,ntype
       iexcl(i)=0
      enddo

      do i=1,ntype-1

       iexcl(i)=0

       do j=i+1,ntype

        if(typmol(i).eq.typmol(j))then

         iexcl(i)=iexcl(i)+1
         if(iexcl(i).gt.mxtype)then
          write(output,"(/,1x,'mxtype exceeded')")
          stop
         endif
         exclpair(i,iexcl(i))=j-i
        
        endif

       enddo

      enddo

      nelist=0

      do i=1,natms

       epoint(i)=nelist+1
       eatm=atmtype(i)
       nexcl=iexcl(eatm)

       if(nexcl.gt.0)then

        do j=1,nexcl

         nelist=nelist+1
         if(nelist.gt.(mxatms*mxelist))then
          write(output,"(/,1x,'mxatms*mxelist exceeded')")
          stop
         endif
         if((nelist-epoint(i)+1).gt.mxnlist)then
          write(output,"(/,1x,'mxnlist exceeded')")
          stop
         endif

         elist(nelist)=i+exclpair(eatm,j)

        enddo

       endif

      enddo

      epoint(natms)=nelist+1

      return

      END
