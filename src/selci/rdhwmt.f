C> \ingroup selci
C> @{
      subroutine selci_rdhwmt(iflwmt,mmulti, nsmax, nf, nfmax, nfmax2)
*
* $Id$
*
      dimension nf(0:32)
c
      call ifill(33,0,nf,1)
      read (iflwmt,*) mmulti, nsmax
      read (iflwmt,*) (nf(i),i=mod(nsmax,2),nsmax,2)
      nfmax = nf(nsmax)
      nfmax2 = nf(nsmax-2)
c
      end
C> @}
