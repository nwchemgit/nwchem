      subroutine tddft_epsilon(ivector,g_trials,g_prod_s,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,evl)
c
c $Id: tddft_epsilon.F,v 1.5 2003-10-17 22:58:29 carlfahl Exp $
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
c
      integer i,j,k,l
      double precision t,p
c
c --------------------------------------------------------
c Add (epsilon_a-epsilon_i)*trial vector diagonal elements
c --------------------------------------------------------
c
      iproc=ga_nodeid()
      call ga_sync()
      do i=1,ipol
        call ga_distribution(g_trials(i),iproc,ilo,ihi,jlo,jhi)
        if ((jlo.gt.ivector).or.(jhi.lt.ivector))
     1    call errquit('tddft_epsilon: an internal error detected',0,
     &       GA_ERR)
        do l=ilo,ihi
          k=mod(l-1,nmo(i)-nfv(i)-nocc(i))+nocc(i)+1
          j=(l-k+nocc(i))/(nmo(i)-nfv(i)-nocc(i))+nfc(i)+1
          call ga_get(g_trials(i),l,l,ivector,ivector,t,1)
          call ga_get(g_prod_s(i),l,l,ivector,ivector,p,1)
          p=p+(evl(k,i)-evl(j,i))*t
          call ga_put(g_prod_s(i),l,l,ivector,ivector,p,1)
        enddo
      enddo
      call ga_sync()
c
c ------
c Return
c ------
c
      return
      end
