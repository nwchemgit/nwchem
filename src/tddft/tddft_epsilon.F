      subroutine tddft_epsilon(ivector,g_trials,g_prod_s,
     1  ipol,nbf_ao,nocc,nmo,evl)
c
c $Id: tddft_epsilon.F,v 1.1 2002-03-22 02:10:36 sohirata Exp $
c
c Add (epsilon_a-epsilon_i) contribution to the product vector.
c
c ! Caution ! Never destroy ivector
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
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
      double precision evl(nbf_ao,2)   ! Eigenvalues
c
      integer i,j,k,l
      double precision t,p
c
c --------------------------------------------------------
c Add (epsilon_a-epsilon_i)*trial vector diagonal elements
c --------------------------------------------------------
c
      do i=1,ipol
        l=0
        do j=1,nocc(i)
          do k=nocc(i)+1,nmo(i)
            l=l+1
            call ga_get(g_trials(i),l,l,ivector,ivector,t,1)
            call ga_get(g_prod_s(i),l,l,ivector,ivector,p,1)
            p=p+(evl(k,i)-evl(j,i))*t
            call ga_put(g_prod_s(i),l,l,ivector,ivector,p,1)
          enddo
        enddo
      enddo
c
c ------
c Return
c ------
c
      return
      end
