      subroutine tddft_diagon(ntrials,tda,g_apb,g_apbvec,
     1  apbval,g_amb,g_ambvec,ambval)
c
c $Id: tddft_diagon.F,v 1.1 2002-03-22 02:10:36 sohirata Exp $
c
c Diagonalize the subspace representation of (A-B) and
c (A-B)^(1/2)(A+B)(A-B)^(1/2) matrices and return the
c eigenvalues and eigenvectors.  (A+B) (g_apb) and (A-B)
c (g_amb) will not be destroyed.
c
c When (not tda) g_ambvec contains subspace representation
c of (X-Y) and g_apbvec contains subspace representation 
c of (X+Y) upon return, after bi-orthogonalization is 
c imposed.
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
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer g_apb            ! GA handle for subspace (A+B) matrix
      integer g_apbvec         ! GA handle for (X+Y)
      integer g_amb            ! GA handle for subspace (A-B) matrix
      integer g_ambvec         ! GA handle for (X-Y)
      double precision apbval(ntrials) ! (A-B)^(1/2)(A+B)(A-B)^(1/2) eigenvalues
      double precision ambval(ntrials) ! (A-B) eigenvalues
c
      integer g_work1          ! GA handle for a work array
      integer g_work2          ! GA handle for a work array
      integer i,n,m
      double precision a
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c ------------------
c Tamm-Dancoff AZ=wZ
c ------------------
c
      if (tda) then
        call ga_diag_std(g_apb,g_apbvec,apbval)
        if (util_print('diagonization',print_high)) then
          if (nodezero) call ma_print(apbval,ntrials,1,'Eigenvalues')
          call ga_print(g_apbvec)
        endif
c
c --------------------------------------------------------
c Full linear response (A-B)^(1/2)(A+B)(A-B)^(1/2) Z=w^2 Z
c --------------------------------------------------------
c
      else
c
c Allocate work arrays for (A-B)^(1/2) and (A-B)^(-1/2)
c
        if (.not.ga_create(mt_dbl,ntrials,ntrials,'(A+B)^(1/2)',
     1    -1,-1,g_work1)) call errquit
     2    ('tddft_diagon: failed to create g_work1',0)
        if (.not.ga_create(mt_dbl,ntrials,ntrials,'(A-B)^(-1/2)',
     1    -1,-1,g_work2)) call errquit
     2    ('tddft_diagon: failed to create g_work2',0)
c
c Diagonalize (A-B) (even when it is already diagonal)
c
        call ga_diag_std(g_amb,g_ambvec,ambval)
        if (util_print('diagonization',print_high)) then
          if (nodezero) call ma_print(ambval,ntrials,1,'Eigenvalues')
          call ga_print(g_ambvec)
        endif
c
c Make (A-B)^(1/2) and (A-B)^(-1/2)
c
        call ga_zero(g_work1)
        call ga_zero(g_work2)
        do i=1,ntrials
          if (ambval(i).lt.0.0d0) then
            call errquit
     1      ('tddft_diagon: A-B is not positive-definite; try TDA',0)
          else
            call ga_put(g_work1,i,i,i,i,dsqrt(ambval(i)),1)
            call ga_put(g_work2,i,i,i,i,1.0d0/dsqrt(ambval(i)),1)
          endif
        enddo
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_ambvec,
     1    g_work1,0.0d0,g_apbvec)
        call ga_dgemm('N','T',ntrials,ntrials,ntrials,1.0d0,g_apbvec,
     1    g_ambvec,0.0d0,g_work1)
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_ambvec,
     1    g_work2,0.0d0,g_apbvec)
        call ga_dgemm('N','T',ntrials,ntrials,ntrials,1.0d0,g_apbvec,
     1    g_ambvec,0.0d0,g_work2)
        if (util_print('diagonization',print_high)) then
          call ga_print(g_work1)
          call ga_print(g_work2)
        endif
c
c Make (A-B)^(1/2)(A+B)(A-B)^(1/2)
c
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_work1,
     1    g_apb,0.0d0,g_apbvec)
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_apbvec,
     1    g_work1,0.0d0,g_ambvec)
c
c Diagonalize (A-B)^(1/2)(A+B)(A-B)^(1/2)
c
        call ga_diag_std(g_ambvec,g_apbvec,apbval)
        do i=1,ntrials
          if (apbval(i).lt.0.0d0) then
            call errquit
     1      ('tddft_diagon: negative excitation energy',0)
          else
            apbval(i)=dsqrt(apbval(i))
          endif
        enddo
        if (util_print('diagonization',print_high)) then
          if (nodezero) call ma_print(apbval,ntrials,1,'Eigenvalues')
          call ga_print(g_apbvec)
        endif
c
c Make (X-Y) = w(A-B)^(-1/2) Z
c
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_work2,
     1    g_apbvec,0.0d0,g_ambvec)
        do i=1,ntrials
          call ga_scale_patch(g_ambvec,1,ntrials,i,i,apbval(i))
        enddo
c
c Make (X+Y) = (A-B)^(1/2) Z
c
        call ga_dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,g_work1,
     1    g_apbvec,0.0d0,g_work2)
        call ga_copy(g_work2,g_apbvec)
c
c Bi-orthogonalize
c
        do n=1,ntrials
          do m=1,ntrials
            if (n.ne.m) then
              a=ga_ddot_patch(g_apbvec,'N',1,ntrials,n,n,
     1                        g_ambvec,'N',1,ntrials,m,m)
              call ga_add_patch(1.0d0,g_apbvec,1,ntrials,n,n,
     1                          -a,   g_ambvec,1,ntrials,m,m,
     2                                g_apbvec,1,ntrials,n,n)
            endif
          enddo
        enddo
c
c Normalize
c
        do n=1,ntrials
          a=ga_ddot_patch(g_apbvec,'N',1,ntrials,n,n,
     1                    g_ambvec,'N',1,ntrials,n,n)
          if (a.eq.0.0d0) then
            call errquit('tddft_diagon: failed to bi-orthonormalize',0)
          else
            call ga_scale_patch(g_apbvec,1,ntrials,n,n,1.0d0/dsqrt(a))
            call ga_scale_patch(g_ambvec,1,ntrials,n,n,1.0d0/dsqrt(a))
          endif
        enddo
c
c Destroy work GA
c
        if (.not.ga_destroy(g_work1)) call errquit
     1    ('tddft_diagon: failed to destroy g_work1',0)
        if (.not.ga_destroy(g_work2)) call errquit
     1    ('tddft_diagon: failed to destroy g_work2',0)
      endif
c           
c ------
c Return
c ------
c
      return
      end
