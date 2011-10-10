      subroutine hnd_prop_dens_make (rtdb,geom,basis,nbf,nmo,nclosed,
     &                               nopen,nvirt,scftyp,
     &                               vectors,occ,g_dens,ndens)
c
c     Calculate the density matrices from the movecs 
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer     rtdb          ! [input] database handle
      integer     geom          ! [input] geometry handle
      integer     basis         ! [input] handles to basis
      integer     nbf        ! [input] number of basis functions
      integer     nmo        ! [input] number of MOs
      integer     nclosed(2) ! [input] number of closed shells
      integer     nopen(2)   ! [input] number of open shells
      integer     nvirt(2)   ! [input] number of virtuals
      character*3 scftyp     ! [input] wave function type (RHF or UHF)
      integer     vectors(2) ! [input] handles for MO vectors
      integer     g_dens(3)  ! [output] density handles
      integer     ndens      ! [output] number of active density handles
      double precision     occ(nbf,2) ! [input] array with occupation numbers
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c    
      integer i, iopen
c
c     Create density matrix handles
c
      ndens = 1
      if (scftyp.eq.'UHF') ndens = 3
c
      do i = 1, ndens
        g_dens(i) = ga_create_atom_blocked(geom,basis,'density matrix')
        call ga_zero(g_dens(i))
      enddo
c
c     Create density matrix/ces
c
      if (scftyp.eq.'RHF') then
        if (nclosed(1).gt.0) call ga_matmul_patch('n','t',2.0D0,0.0d0,
     $          vectors(1),1,nbf,1,nclosed(1),vectors(1),1,nclosed(1),
     $          1,nbf,g_dens(1),1,nbf,1,nbf)
        if (nopen(1) .gt. 0) then
           do iopen = 1, nopen(1)
             call ga_matmul_patch(
     $             'n','t',occ(nclosed(1)+iopen,1),1.d0,
     $       vectors(1),1,nbf,nclosed(1)+iopen,nclosed(1)+iopen,
     $       vectors(1),nclosed(1)+iopen,nclosed(1)+iopen,1,nbf,
     $       g_dens(1),1,nbf,1,nbf)
           enddo
        endif
        call ga_symmetrize(g_dens(1))
      else  ! UHF
        do iopen = 1, nopen(1)
           call ga_matmul_patch('n', 't', occ(iopen,1), 1.0d0,
     $                          vectors(1),1,nbf,iopen,iopen,
     $                          vectors(1),iopen,iopen,1,nbf,
     $                          g_dens(1),1,nbf,1,nbf)
        enddo
        do iopen = 1, nopen(2)
           call ga_matmul_patch('n', 't', occ(iopen,2), 1.0d0,
     $                          vectors(2),1,nbf,iopen,iopen,
     $                          vectors(2),iopen,iopen,1,nbf,
     $                          g_dens(2),1,nbf,1,nbf)
        enddo
c
c       For UHF dens(1) and dens(2) are alpha and beta densities
c       Also create a total density
c
        call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2), g_dens(3) )
        call ga_symmetrize(g_dens(1))
        call ga_symmetrize(g_dens(2))
        call ga_symmetrize(g_dens(3))
      endif
c
      end
c $Id$
