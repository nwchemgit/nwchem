      subroutine tddft_epsilon(ivector,g_trials,g_prod_s,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,evl,lecut,ecut,diff_max,
     2  lowin, owstart, owend, lewin, ewinl, ewinh, vshift)
c
c $Id$
c
c Add (epsilon_a-epsilon_i) contribution to the product vector.
c This subroutine is explicitly parallelized.
c
c ! Caution ! Never destroy ivector
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer ivector          ! Vector to be processed
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_s(2)      ! GA handle for trial vectors
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nbf_ao           ! Number of AO basis functions
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer iproc            ! Processor number
      integer ilo,ihi,jlo,jhi  ! Distribution of GA data
      double precision evl(nbf_ao,2)   ! Eigenvalues
      integer vshift
      logical lecut
      double precision ecut
      logical lewin
      double precision ewinl, ewinh
      logical lowin
      integer owstart(2)
      integer owend(2)
      double precision diff_max
c
      integer i,j,k,l
      double precision t,p
      double precision evl_diff  ! eigenvalue difference
c
c --------------------------------------------------------
c Add (epsilon_a-epsilon_i)*trial vector diagonal elements
c --------------------------------------------------------
c
      iproc=ga_nodeid()
      call ga_sync()
c
      do i=1,ipol
       call ga_distribution(g_trials(i),iproc,ilo,ihi,jlo,jhi)
       if ((jlo.le.ivector).and.(jhi.ge.ivector)) then
        do l=ilo,ihi
          if (abs(vshift).gt.0) then
             k=nocc(i)+1+vshift
          else
             k=mod(l-1,nmo(i)-nfv(i)-nocc(i))+nocc(i)+1
          end if
          j=(l-k+nocc(i))/(nmo(i)-nfv(i)-nocc(i))+nfc(i)+1
          call ga_get(g_trials(i),l,l,ivector,ivector,t,1)
          call ga_get(g_prod_s(i),l,l,ivector,ivector,p,1)
          evl_diff = evl(nmo(i),i)-evl(1,i)  ! span the eigenvalues
          if (lecut) then
            if (evl(j,i).le.ecut) then
              evl_diff = evl(k,i)-evl(j,i)
            end if
          else if (lowin) then
            if (j.ge.owstart(i) .and. j.le.owend(i)) then
              evl_diff = evl(k,i)-evl(j,i)
            end if
          else if (lewin) then
            if (evl(j,i).le.ewinh .and. evl(j,i).ge.ewinl) then
              evl_diff = evl(k,i)-evl(j,i)
            end if
          else
            evl_diff = evl(k,i)-evl(j,i)
          endif
          p=p+(evl_diff)*t
          call ga_put(g_prod_s(i),l,l,ivector,ivector,p,1)
        enddo
       endif ! jlo, jhi check
      enddo
      call ga_sync()
c
c ------
c Return
c ------
c
      return
      end
