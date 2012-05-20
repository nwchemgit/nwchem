      SUBROUTINE list_verlt()

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_cuto.inc'
      include 'cm_latt.inc'
      include 'cm_vlst.inc'
      include 'cm_elst.inc'

      integer i,j,k
      integer eatm,nlist

      real*8 rij,rijsq

      dimension rij(mxatms,3)

      nlist=0

      do i=1,natms-1

       k=0
       point(i)=nlist+1
       if(epoint(i).ne.epoint(i+1))eatm=epoint(i)

       do j=i+1,natms

        k=k+1
        rij(k,1)=ccc(i,1)-ccc(j,1)
        rij(k,2)=ccc(i,2)-ccc(j,2)
        rij(k,3)=ccc(i,3)-ccc(j,3)

       enddo
c
       call tool_rebox(k,mxatms,latt,rlatt,rij)

c

       k=0

       do j=i+1,natms

        k=k+1

        if((epoint(i).ne.epoint(i+1)).and.(elist(eatm).eq.j))then

         eatm=min(eatm+1,(epoint(i+1)-1))

        else

         rijsq=rij(k,1)*rij(k,1)+rij(k,2)*rij(k,2)+rij(k,3)*rij(k,3)


         if(rijsq.lt.vcutsq)then

          nlist=nlist+1

          if(nlist.gt.(mxatms*mxvlist))then
           write(output,"(/,1x,'mxatms*mxvlist exceeded')")
           stop
          endif
          if((nlist-point(i)+1).gt.mxnlist)then
           write(output,"(/,1x,'mxnlist exceeded')")
           stop
          endif

          list(nlist)=j

         endif

        endif

       enddo

      enddo

      point(natms)=nlist+1



      return

      END
c $Id$
