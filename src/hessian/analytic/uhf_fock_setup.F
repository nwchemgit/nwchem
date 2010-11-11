c
      subroutine uhf_fock_setup(g_dens, g_fock, jfac, kfac, ndens,
c
     &   nfock, nbf, oprint)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
c  This routine takes the multidimensional density and the fock matrices
c  that are set up for alpha and beta and expand them to be alpha and
c  beta/coulomb and exchange (2 -> 4) so that the correct summation can
c  be performed after the shell_fock_build call.
c
      integer g_dens       ! [input/output] density
      integer g_fock       ! [input/output] Fock matrices
      integer ndens        ! [input/output] number of density matrices
      integer nfock        ! [input/output] number of Fock matrices
      integer nbf          ! [input] number of basis functions
      double precision jfac(ndens)  ! [input/output] Coulomb scale factors
      double precision kfac(ndens)  ! [input/output] Exchange scale factors
      logical oprint       ! [input] Print output
c
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer g_dens2, g_fock2
      integer i, noffset, noffset2, noffset3
      double precision zero
      parameter (zero = 0.0D0)
c
c First set up jfac and kfac correctly.  This assumes that there is
c already enough space in jfac and kfac to expand from 2 to 4.
c
      noffset = ndens/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      do i = 1, noffset
        jfac(noffset3+i) = zero
        jfac(noffset2+i) = jfac(noffset+i)
        jfac(noffset+i)  = zero
        kfac(noffset3+i) = kfac(noffset+i)
        kfac(noffset2+i) = zero
        kfac(noffset+i)  = kfac(i)
        kfac(i)          = zero
      enddo
c
c Next expand the density.  We create a new ga and substitute it
c for the original.
c
      dims(1) = ndens * 2
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'Density',chunk,
     &     g_dens2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_dens2',555, GA_ERR)
c
c Copy "backwards" so that we don't loose any data
c
      alo(1) = noffset + 1 ! Beta
      ahi(1) = ndens
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = noffset*3 + 1
      bhi(1) = noffset*4
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      alo(1) = 1       ! Alpha
      ahi(1) = noffset
      blo(1) = noffset + 1
      bhi(1) = ndens
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
c
c  Now get rid of original density
c
      if (.not. ga_destroy(g_dens)) call errquit
     *   ('uhf_fock_setup: failed to free g_dens', ndens, GA_ERR)
c
c  Assign g_dens to be the new matrix and update ndens
c
      g_dens = g_dens2
      ndens = ndens * 2
c
c Now expand the fock matrices doing what we did above for
c the densities.
c
      dims(1) = nfock * 2
      chunk(1) = dims(1)
      if (.not. nga_create (MT_DBL, 3, dims, 'Fock matrices',
     &    chunk, g_fock2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_fock2',555, GA_ERR)
c
      noffset = nfock/2
      alo(1) = noffset + 1  ! Beta
      ahi(1) = nfock
      blo(1) = noffset*3 + 1
      bhi(1) = noffset*4
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      alo(1) = 1       ! Alpha
      ahi(1) = noffset
      blo(1) = noffset + 1
      bhi(1) = nfock
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
c
      if (.not. ga_destroy(g_fock)) call errquit
     *   ('uhf_fock_setup: failed to free g_fock', nfock, GA_ERR)
c
      g_fock = g_fock2
      nfock = nfock * 2
c
      if (oprint) then
        if (ga_nodeid().eq.0) write(LuOut,*)
     &     'Density and Fock matrices after uhf_fock_setup'
        do i=1,ndens
           ahi(1)=i
           alo(1)=i
           call nga_print_patch(g_dens,alo,ahi,0)
        enddo
        call ga_print(g_fock)
      endif
c
      return
      end
