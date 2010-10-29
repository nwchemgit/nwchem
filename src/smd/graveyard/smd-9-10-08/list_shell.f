c
c $Id$
c

      SUBROUTINE list_shell(nshel,ntshel,shelatm)

      implicit none

      include 'p_array.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_shel.inc'

      integer i,j,k
      integer nshel,ntshel,shelatm 

      dimension shelatm(mxshel,2)

      ntshel=0

      do i=1,natms

       do j=1,nshel

        if(atmtype(i).eq.shelatm(j,1))then

         ntshel=ntshel+1
         if(ntshel.gt.mxshel2)then
          write(output,*)'Increase mxshel2 to ',ntshel
          stop
         endif

         k=shelatm(j,2)-shelatm(j,1)
         if(typmol(atmtype(i)).ne.typmol(atmtype(i+k)))then
          write(output,*)'Shell and core must be from the same molecule'
          stop
         endif

         shllist(ntshel,1)=i
         shllist(ntshel,2)=i+k
         shllist(ntshel,3)=j

        endif

       enddo

      enddo
 
      return

      END
