      subroutine selci_prtcon(ifllog, norbs, ioconf, nintpo, nbitpi)
*
* $Id: prtcon.f,v 1.2 1997-10-31 23:42:20 d3e129 Exp $
*
      dimension ioconf(nintpo),iocc(255)
c
      call selci_upkcon(norbs, iocc, ioconf, nintpo, nbitpi)
      call selci_wrtcon(ifllog, iocc, norbs)
c
      end
