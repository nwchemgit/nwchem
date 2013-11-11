C> \ingroup selci
C> @{
      subroutine selci_wrtcon(ifllog,iocc,norbs)
*
* $Id$
*
      dimension iocc(norbs),list(255),locc(255)
c
c     print out the orbital occupation in iocc
c     0 = uocc, 1 = socc, 3 = docc
c
      igot = 0
      do 10 i = 1,norbs
         if (iocc(i).gt.0) then
            igot = igot + 1
            list(igot) = i
            if (iocc(i) .eq. 1) then
               locc(igot) = 1
            else
               locc(igot) = 2
            endif
         endif
 10   continue
c
      if (norbs.ge.100) then
         write(ifllog,1) (list(i),locc(i),i=1,igot)
 1       format(1x,9(i4,'(',i1,')':))
      else
         write(ifllog,2) (list(i),locc(i),i=1,igot)
 2       format(1x,11(i3,'(',i1,')':))
      endif
c
      end
C> @}

            
      
