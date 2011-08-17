c
c     $Id$
c
c****f* uhf/uhf_fock_finish
c
c     NAME
c
c       uhf_fock_finish -- take the Coulomb and exchange matrices from
c       shell_fock_build and construct alpha and beta Fock matrices
c       from them
c
c     FUNCTION
c
c       This routine compresses the alpha and beta density matrices
c       back into the usual UHF form. Likewise the Coulomb and exchange
c       pre-factors (jfac and kfac respectively) are compressed.
c
c       In addition this routine uses the Coulomb and exchange, alpha
c       and beta Fock matrices and combines these to form alpha and
c       beta Fock matrices. The expressions used in this transformation
c       are:
c       * alpha-Fock = alpha-Coulomb + beta-Coulomb - alpha-exchange
c       * beta-Fock  = alpha-Coulomb + beta-Coulomb - beta-exchange
c       Schematically the operation can be presented as (ndens and 
c       nfock are the corresponding values on output):
c
c         Input                               Output
c         =====                               ======
c
c         -------                             -------
c           ^                                   ^
c         ndens/2  alpha-density matrices       |
c           v                                 ndens/2  alpha-density
c         -------                               |      matrices
c           ^                                   v
c         ndens/2  alpha-density matrices     -------
c           v                                   ^
c         -------                         ==>   |
c           ^                                 ndens/2  beta-density
c         ndens/2  beta-density matrices        |      matrices
c           v                                   v
c         -------                             -------
c           ^
c         ndens/2  beta-density matrices
c           v
c         -------
c
c
c         -------                              -------
c           ^                                    ^
c         nfock/2  alpha-Coulomb matrices        |
c           v                                  nfock/2  alpha-Fock
c         -------                                |      matrices
c           ^                                    v
c         nfock/2  alpha-exchange matrices     -------
c           v                                    ^
c         -------                          ==>   |
c           ^                                  nfock/2  beta-Fock
c         nfock/2  beta-Coulomb matrices         |      matrices
c           v                                    v
c         -------                              -------
c           ^
c         nfock/2  beta-exchange matrices
c           v
c         -------
c
c
c     SYNOPSIS
c
      subroutine uhf_fock_finish(g_dens, g_fock, jfac, kfac,
     &                           ndens, nfock, nbf, oprint)
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
c  This routine takes the results of shell_fock_build and forms the
c  complete alpha and beta Fock matrices (4 -> 2).
c
c     INPUTS
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
c     SOURCE
c
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer g_dens2, g_fock2
      integer i, noffset, noffset2, noffset3
      double precision zero, one, mone
      parameter (zero =  0.0D0)
      parameter (one  =  1.0D0)
      parameter (mone = -1.0D0)
c
c     Put jfac and kfac back to their original form.
c
      ndens = ndens/2
      noffset = ndens/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      do i = 1, noffset
        jfac(noffset+i)  = jfac(noffset2+i)
        kfac(i)          = kfac(noffset+i)
        kfac(noffset+i)  = kfac(noffset3+i)
      enddo
c
c     Put the density back to its original.  We create a new ga and 
c     substitute it for the original.
c
      dims(1) = ndens
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'Density',chunk,
     &     g_dens2)) call errquit
     &     ('uhf_fock_finish: could not allocate g_dens2',555, GA_ERR)
c
c     Copy data
c
      alo(1) = noffset*3 + 1  ! Beta
      ahi(1) = noffset*4
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = noffset + 1
      bhi(1) = ndens
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      alo(1) = 1   ! Alpha
      ahi(1) = noffset
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
c
c     Now get rid of original density
c
      if (.not. ga_destroy(g_dens)) call errquit
     *   ('uhf_fock_finish: failed to free g_dens', ndens, GA_ERR)
c
c     Assign g_dens to be the new matrix and update ndens
c
      g_dens = g_dens2
c
c     Now get the Fock matrices set up.
c
      nfock = nfock/2
      dims(1) = nfock
      chunk(1) = dims(1)
      if (.not. nga_create (MT_DBL, 3, dims, 'Fock matrices',
     &    chunk, g_fock2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_fock2',555, GA_ERR)
c
      noffset = nfock/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      alo(1) = 1  ! Coulomb part (add alpha and beta)
      ahi(1) = noffset
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      clo(1) = 1
      chi(1) = noffset
      clo(2) = 1
      chi(2) = nbf
      clo(3) = 1
      chi(3) = nbf
      call nga_add_patch(one,g_fock,alo,ahi,
     &                    one,g_fock,blo,bhi,
     &                    g_fock2,clo,chi)
      alo(1) = noffset + 1
      ahi(1) = nfock
      call nga_copy_patch('N',g_fock2,clo,chi,
     &    g_fock2,alo,ahi)
c
c     Add exchange of same spin
c
      alo(1) = noffset + 1   ! Alpha
      ahi(1) = nfock
      blo(1) = 1
      bhi(1) = noffset
      call nga_add_patch(mone,g_fock,alo,ahi,
     &                   one,g_fock2,blo,bhi,
     &                    g_fock2,blo,bhi)
c
      alo(1) = noffset*3 + 1 ! Beta
      ahi(1) = noffset*4
      blo(1) = noffset + 1
      bhi(1) = nfock
      call nga_add_patch(mone,g_fock,alo,ahi,
     &                   one,g_fock2,blo,bhi,
     &                    g_fock2,blo,bhi)
c
c     Now get rid of g_fock and assign the new one (g_fock2) to g_fock
c
      if (.not. ga_destroy(g_fock)) call errquit
     &   ('uhf_fock_setup: failed to free g_fock', nfock, GA_ERR)
c
      g_fock = g_fock2
c
      if (oprint) then
        if (ga_nodeid().eq.0) then
           write(LuOut,*)
     &        'Density and Fock matrices after uhf_fock_finish'
           call util_flush(LuOut)
        endif
        do i=1,ndens
           ahi(1)=i
           alo(1)=i
           call nga_print_patch(g_dens,alo,ahi,0)
        enddo
        do i=1,nfock
           ahi(1)=i
           alo(1)=i
           call nga_print_patch(g_fock,alo,ahi,0)
        enddo
      endif
c
      return
      end
c
c******
