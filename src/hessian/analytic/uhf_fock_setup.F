c
c     $Id$
c
c****f* uhf/uhf_fock_setup
c
c     NAME
c
c       uhf_fock_setup -- set up pre-factors for Coulomb and exchange
c       terms, density matrices and Fock matrices for the UHF 
c       shell_fock_build call
c
c     FUNCTION
c
c       This routine takes the multi-dimensional density and Fock 
c       matrices corresponding to alpha and beta spin channel 
c       quantities and expands them to alpha Coulomb and exchange
c       matrices and beta Coulomb and exchange matrices.
c
c       After the call to shell_fock_build the sister routine 
c       uhf_fock_finish will combine the resulting matrices to 
c       construct alpha Fock matrices and beta Fock matrices.
c
c       The mathematics used in uhf_fock_finish imposes an interesting
c       constraint on how we copy the input Fock matrices here in
c       uhf_fock_setup. Looking for example at the alpha-Fock matrices,
c       these are constructed as
c
c         alpha-Fock = alpha-Coulomb + beta-Coulomb - alpha-Exchange
c
c       This means that alpha-Coulomb and alpha-exchange should be
c       initialized with the beta input Fock matrices (beta-Fock0 for
c       short), likewise beta-Coulomb and beta-exchange should be
c       initialized with alpha-Fock0. To verify this simply substitute
c       those quantities in the above expression to get:
c
c         alpha-Fock = beta-Fock0 + alpha-Fock0 - beta-Fock0
c                    = alpha-Fock0
c
c       to see that this way the input alpha-Fock matrices contribute
c       to the output alpha-Fock matrices. There is no problem with
c       respect to the Coulomb and exchange contributions that 
c       shell_fock_build has added, because they appear in the right
c       places.
c
c       Schematically the operation of this routine can be presented
c       as (ndens and nfock are the corresponding values on input):
c
c          Input                        Output
c          =====                        ======
c
c          -------                      -------
c            ^                            ^
c            |                          ndens/2  alpha-density matrices
c          ndens/2  alpha-density         v
c            |      matrices            -------
c            v                            ^
c          -------                ===>  ndens/2  alpha-density matrices
c            ^                            v
c            |                          -------
c          ndens/2  beta-density          ^
c            |      matrices            ndens/2  beta-density matrices
c            v                            v
c          -------                      -------
c                                         ^
c                                       ndens/2  beta-density matrices
c                                         v
c                                       -------
c
c
c                                       -------
c                                         ^
c          -------                      nfock/2  alpha-Coulomb matrices
c            ^                            v
c            |                          -------
c          nfock/2  alpha-Fock            ^
c            |      matrices         >  nfock/2  alpha-exchange matrices
c            v                    \ /     v
c          -------                 X    -------
c            ^                    / \     ^
c            |                       >  nfock/2  beta-Coulomb matrices
c          nfock/2  beta-Fock             v
c            |      matrices            -------
c            v                            ^
c          -------                      nfock/2  beta-exchange matrices
c                                         v
c                                       -------
c
c       The pre-factors for the Coulomb terms (jfac) and the exchange
c       terms (kfac) undergo a corresponding copying operation as well.
c
c     SYNOPSIS
c
      subroutine uhf_fock_setup(g_dens, g_fock, jfac, kfac,
     &                          ndens, nfock, nbf, oprint)
c
c     SOURCE
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
c     INPUTS
c
      integer g_dens       ! [input/output] density
      integer g_fock       ! [input/output] Fock matrices
      integer ndens        ! [input/output] number of density matrices
      integer nfock        ! [input/output] number of Fock matrices
      integer nbf          ! [input] number of basis functions
      double precision jfac(ndens) ! [input/output] Coulomb scale factors
      double precision kfac(ndens) ! [input/output] Exchange scale factors
      logical oprint       ! [input] Print output
c
c     SOURCE
c
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer g_dens2, g_fock2
      integer i, noffset, noffset2, noffset3
      double precision zero
      parameter (zero = 0.0D0)
c
c     First set up jfac and kfac correctly.  This assumes that there is
c     already enough space in jfac and kfac to expand from 2 to 4.
c
      noffset = ndens/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      do i = 1, noffset
        jfac(noffset3+i) = zero
        jfac(noffset2+i) = jfac(noffset+i)
        jfac(noffset+i)  = zero
c       jfac(i)          = already set correctly
        kfac(noffset3+i) = kfac(noffset+i)
        kfac(noffset2+i) = zero
        kfac(noffset+i)  = kfac(i)
        kfac(i)          = zero
      enddo
c
c     Next expand the density.  We create a new ga and substitute it
c     for the original.
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
c     Copy "backwards" so that we don't loose any data
c     Copy g_dens{a(1:ndens/2),b(1:ndens/2)}
c     --> g_dens2{a(1:ndens/2),a(1:ndens/2),b(1:ndens/2),b(1:ndens/2)}
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
c     Now get rid of original density
c
      if (.not. ga_destroy(g_dens)) call errquit
     *   ('uhf_fock_setup: failed to free g_dens', ndens, GA_ERR)
c
c     Assign g_dens to be the new matrix and update ndens
c
      g_dens = g_dens2
      ndens = ndens * 2
c
c     Now expand the fock matrices doing what we did above for
c     the densities.
c     Copy g_fock{a(1:nfock/2),b(1:nfock/2)}
c     --> g_fock2{b(1:nfock/2),b(1:nfock/2),a(1:nfock/2),a(1:nfock/2)}
c
      dims(1) = nfock * 2
      chunk(1) = dims(1)
      if (.not. nga_create (MT_DBL, 3, dims, 'Fock matrices',
     &    chunk, g_fock2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_fock2',555, GA_ERR)
c
      noffset = nfock/2
      alo(1) = 1             ! Alpha input Fock matrices
      ahi(1) = noffset
      blo(1) = noffset*3 + 1 ! Beta exchange matrices
      bhi(1) = noffset*4
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = noffset*2 + 1 ! Beta Coulomb matrices
      bhi(1) = noffset*3
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      alo(1) = noffset + 1   ! Beta input Fock matrices
      ahi(1) = nfock
      blo(1) = noffset + 1   ! Alpha exchange matrices
      bhi(1) = nfock
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = 1             ! Alpha Coulomb matrices
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
        if (ga_nodeid().eq.0) then
          write(LuOut,*)'Density and Fock matrices after uhf_fock_setup'
          call util_flush(LuOut)
        endif
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
c
c******
