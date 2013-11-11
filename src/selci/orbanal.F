C> \ingroup selci
C> @{
      subroutine selci_orbanal(noconf, ioconf, nintpo, nbitpi, norbs,
     $     map)
      implicit none
*
* $Id$
*
      integer noconf, nintpo, nbitpi, norbs
      integer ioconf(nintpo, noconf)
      integer map(norbs)
c
c     Count the number of times each orbital is referenced in
c     an orbital configuration ... an attempt to partition the
c     reference space into computationally useful weakly and
c     strongly occupied orbital subspaces.
c
      integer iocc(255), npres(255)
      integer iconf, i, j, itmp
c
      call ifill(norbs, 0, npres, 1)
c
      do iconf = 1, noconf
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         do i = 1, norbs
            if (iocc(i) .eq. 1) then
               npres(i) = npres(i) + 1
            else if (iocc(i) .eq. 3) then
               npres(i) = npres(i) + 2
            endif
         enddo
      enddo
c
      write(6,1) 
      write(6,2) (dble(npres(i))/noconf, i=1,norbs)
 1    format(/' Average occupancy of each orbital ')
 2    format(15f5.2)
c
c     Now sort the orbitals by their occupation
c
      do i = 1, norbs
         map(i) = i
         do j = 1, i-1
            if (npres(map(j)).lt.npres(map(i))) then
              itmp = map(i)
              map(i) = map(j)
              map(j) = itmp
           endif
        enddo
      enddo
c
      write(6,*)
      write(6,*) ' Orbitals sorted by occupation number '
      write(6,3) (map(i),i=1,norbs)
 3    format(19i4)
      write(6,*)
c
      end
      subroutine makeicon4(noconf, ioconf, nintpo, nbitpi, norbs,
     $     map, icon4)
      implicit none
      integer noconf
      integer nintpo
      integer nbitpi
      integer ioconf(nintpo,noconf)
      integer norbs
      integer map(norbs)
      integer icon4(4,noconf)
c
c     In order to expedite figuring which configurations interact
c     compute the occupation of each reference in 4 arbitrary 
c     orbital subspaces.  The map is used so we can play with 
c     the ordering of orbitals.
c
      integer iconf, iocc(255),k
c
      call ifill(4*noconf, 0, icon4, 1)
c
      do iconf = 1, noconf
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         call doicon4(norbs,iocc,map,icon4(1,iconf))
         write(6,1) iconf,(icon4(k,iconf),k=1,4)
 1       format(1x,i5,5x,4i4)
      enddo
c
      end
      subroutine doicon4(norbs,iocc,map,icon4)
      implicit none
      integer norbs
      integer iocc(norbs)
      integer map(norbs)
      integer icon4(4)
c
      integer ilo, ihi, i, n4, k, itest
c
      n4 = norbs/4
      if (n4.eq.0) n4 = 1
c
      icon4(1)=0
      icon4(2)=0
      icon4(3)=0
      icon4(4)=0
c
      k = 1
      do ilo = 1, norbs, n4
         ihi = min(ilo+n4-1,norbs)
         do i = ilo, ihi
            itest = iocc(map(i))
            if (itest .eq. 1) then
               icon4(k) = icon4(k) + 1
            else if (itest .eq. 3) then
               icon4(k) = icon4(k) + 2
            endif
         enddo
         k = k + 1
      enddo
c
      end
C> @}
