C> \file mkindb.f
C> Make indbar array
C>
C> \ingroup selci
C> @{
C>
C> \brief Make indbar array and lists of socc/docc orbitals
C>
      subroutine selci_mkindb(norbs, iocc, indbar, listd, lists, ns, nd)
      implicit none
*
* $Id$
*
      integer norbs       !< [Input] The number of orbitals
      integer iocc(norbs) !< [Input] The occupation of each orbital
                          !< - 0: empty
                          !< - 1: singly occupied
                          !< - 3: doubly occupied
      integer indbar(norbs) !< [Output] One more than the number of
                            !< singly occupied orbitals of lower orbital
                            !< number
      integer lists(norbs)  !< [Output] Ordered list of singly occupied
                            !< orbitals
      integer listd(norbs)  !< [Output] Ordered list of doubly occupied
                            !< orbitals
      integer ns            !< [Output] The number of single occupied
                            !< orbitals
      integer nd            !< [Output] The number of doubly occupied
                            !< orbitals
c
      integer i
c
c     make indbar array and lists of socc/docc orbitals
c     indbar(i) = position of orbital i if placed in the ordered
c                 set of singly occupied orbitals
c
      ns = 0
      nd = 0
      do 10 i = 1,norbs
         indbar(i) = ns + 1
         if (iocc(i).eq.1) then
            ns = ns + 1
            lists(ns) = i
         else if (iocc(i).eq.3) then
            nd = nd + 1
            listd(nd) = i
         endif
 10   continue
c
      end
C>
C> @}
