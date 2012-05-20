      SUBROUTINE list_const(ncons,ntcons,consatm)

      implicit none

      include 'p_array.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_cons.inc'

      integer i,j,k
      integer ncons,ntcons,consatm

      dimension consatm(mxcons,2)

      ntcons=0

      do i=1,natms

       do j=1,ncons

        if(atmtype(i).eq.consatm(j,1))then

         ntcons=ntcons+1
         if(ntcons.gt.mxcons2)then
          write(output,*)'Increase mxcons2 to ',ntcons
          stop
         endif

         k=consatm(j,2)-consatm(j,1)
         conlist(ntcons,1)=i
         conlist(ntcons,2)=i+k
         conlist(ntcons,3)=j

        endif

       enddo

      enddo
 
      write(99,*) ntcons
      do i=1,ntcons
       write(99,*) conlist(i,1),conlist(i,2)

      end do
      return

      END
c $Id$
