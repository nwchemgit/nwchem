      subroutine tddft_subspace(ipol,nocc,nmo,g_trials,
     1  g_prod_p,g_prod_m,ntrials,tda,g_apb,g_amb)
c
c $Id: tddft_subspace.F,v 1.1 2002-03-22 02:10:37 sohirata Exp $
c
c Form subspace representation of (A+B) and (A-B) matrices
c using trial and product vectors.  The elements of (A+B)
c and (A-B) matrices are simply:
c
c (A+B)_pq = (A+B) product vector (p) * trial vector (q)
c (A-B)_pq = (A-B) product vector (p) * trial vector (q)
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
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer g_apb            ! GA handle for subspace (A+B) matrix
      integer g_amb            ! GA handle for subspace (A-B) matrix
c
      integer nov(2)           ! Number of occupied-virtual pairs
      double precision a
      integer i,m,n
c
c --------------------
c Allocate work arrays
c --------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nocc(i))*nocc(i)
      enddo
      if (ipol.eq.1) nov(2)=0
c
c -----------------
c Form dot products
c -----------------
c
c Orthonormality of trial vectors (debug)
c
      if (util_print('subspace',print_debug)) then
        do m=1,ntrials
          do n=1,ntrials
            a=0.0d0
            do i=1,ipol
              a=a+ga_ddot_patch(g_trials(i),'N',1,nov(i),m,m,
     1                          g_trials(i),'N',1,nov(i),n,n)
            enddo
            call ga_put(g_apb,m,m,n,n,a,1)
          enddo
        enddo
        call ga_print(g_apb)
      endif
c
c Form A+B and A-B
c
      do m=1,ntrials
        do n=1,ntrials
          a=0.0d0
          do i=1,ipol
            a=a+ga_ddot_patch(g_trials(i),'N',1,nov(i),m,m,
     1                        g_prod_p(i),'N',1,nov(i),n,n)
          enddo
          call ga_put(g_apb,m,m,n,n,a,1)
          if (.not.tda) then
            a=0.0d0
            do i=1,ipol
              a=a+ga_ddot_patch(g_trials(i),'N',1,nov(i),m,m,
     1                          g_prod_m(i),'N',1,nov(i),n,n)
            enddo
            call ga_put(g_amb,m,m,n,n,a,1)
          endif
        enddo
      enddo
c
c --------------
c Print matrices
c --------------
c
      if (util_print('subspace',print_high)) then
        call ga_print(g_apb)
        if (.not.tda) call ga_print(g_amb)
      endif
c           
c ------
c Return
c ------
c
      return
      end
