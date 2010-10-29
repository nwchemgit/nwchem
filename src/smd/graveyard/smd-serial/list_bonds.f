c
c $Id$
c

      SUBROUTINE list_bonds(nbond,ntbond,bondatm)

      implicit none

      include 'p_array.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_bond.inc'

      integer i,j,k
      integer nbond,ntbond,bondatm 

      dimension bondatm(mxbond,2)

      ntbond=0

      do i=1,natms

       do j=1,nbond

        if(atmtype(i).eq.bondatm(j,1))then

         ntbond=ntbond+1
         if(ntbond.gt.mxbond2)then
          write(output,*)'Increase mxbond2 to ',ntbond
          stop
         endif

         k=bondatm(j,2)-bondatm(j,1)
         if(typmol(atmtype(i)).ne.typmol(atmtype(i+k)))then
          write(output,*)'Bonded atoms must be from the same molecule'
          stop
         endif

         bonlist(ntbond,1)=i
         bonlist(ntbond,2)=i+k
         bonlist(ntbond,3)=j

        endif

       enddo

      enddo
 
      return

      END
