      subroutine selci_check(node)
*
* $Id: selci_check.f,v 1.1 2003-04-07 20:55:54 windus Exp $
*
      logical opn
      character*80 filename
      do i = 1, 99
         if (i.ne.5 .and. i.ne.6) then
            inquire(i,opened=opn,name=filename)
            if (opn) then
               write(6,*) ' node ', node, ' unit ', i,
     $              ' name ', filename
            endif
         endif
      enddo
      end

            
