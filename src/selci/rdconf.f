C> \ingroup selci
C> @{
      subroutine selci_rdconf(iflcon,ioconf,indxci,nintpo,noconf)
*
* $Id$
*
      integer ioconf(nintpo*noconf),indxci(noconf+1)
c
c     read occupations and index vector from ciconf file
c
      read (iflcon) ioconf
      read (iflcon) indxci
c
      end
C> @}
