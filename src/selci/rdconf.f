      subroutine selci_rdconf(iflcon,ioconf,indxci,nintpo,noconf)
*
* $Id: rdconf.f,v 1.2 1997-10-31 23:42:22 d3e129 Exp $
*
      integer ioconf(nintpo*noconf),indxci(noconf+1)
c
c     read occupations and index vector from ciconf file
c
      read (iflcon) ioconf
      read (iflcon) indxci
c
      end
