      subroutine tddft_diagon(ntrials,tda,apb,apbvec,
     1  apbval,amb,ambvec,ambval)
c
c $Id$
c
c Diagonalize the subspace representation of (A-B) and
c (A-B)^(1/2)(A+B)(A-B)^(1/2) matrices and return the
c eigenvalues and eigenvectors.  (A+B) (apb) and (A-B)
c (amb) will not be destroyed.
c
c When (not tda) g_ambvec contains subspace representation
c of (X-Y) and g_apbvec contains subspace representation 
c of (X+Y) upon return, after bi-orthogonalization is 
c imposed.
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
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      double precision apb(ntrials,ntrials)     ! Subspace (A+B) matrix
      double precision apbvec(ntrials,ntrials)  ! (X+Y)
      double precision amb(ntrials,ntrials)     ! Subspace (A-B) matrix
      double precision ambvec(ntrials,ntrials)  ! (X-Y)
      double precision apbval(ntrials) ! (A-B)^(1/2)(A+B)(A-B)^(1/2) eigenvalues
      double precision ambval(ntrials) ! (A-B) eigenvalues
c
      integer l_work1,k_work1  ! MA handle for a work array
      integer l_work2,k_work2  ! MA handle for a work array
      integer l_work3,k_work3  ! MA handle for a work array
      integer l_work4,k_work4  ! MA handle for a work array
      integer info
      integer i,j,n,m
      double precision a
      double precision dummy1(ntrials),dummy2
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
      if (.not.ma_push_get(mt_dbl,ntrials*ntrials,'Work 1',
     1  l_work1,k_work1)) call errquit
     2  ('tddft_diagon: failed to allocate work1',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,ntrials*ntrials,'Work 2',
     1  l_work2,k_work2)) call errquit
     2  ('tddft_diagon: failed to allocate work2',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,10*ntrials,'Work 3',
     1  l_work3,k_work3)) call errquit
     2  ('tddft_diagon: failed to allocate work3',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,ntrials*ntrials,'Work 4',
     1  l_work4,k_work4)) call errquit
     2  ('tddft_diagon: failed to allocate work4',0, MA_ERR)
c
c ------------------
c Tamm-Dancoff AZ=wZ
c ------------------
c
      if (tda) then
        do i=1,ntrials
          do j=1,ntrials
            dbl_mb(k_work4+(i-1)*ntrials+j-1)=apb(j,i)
          enddo
        enddo
        call util_dgeev('N','V',ntrials,dbl_mb(k_work4),ntrials,
     1    apbval,dummy1,dummy2,1,apbvec,ntrials,
     2    dbl_mb(k_work3),10*ntrials,info)
        call tddft_schmidt(ntrials,apbvec)
        call util_nwsort(ntrials,apbval,apbvec,'A')
        if (util_print('diagonization',print_debug)) then
          if (nodezero) then
            call ma_print(apbval,ntrials,1,'Eigenvalues')
            call ma_print(apbvec,ntrials,ntrials,'Eigenvectors') 
          endif
        endif
c
c --------------------------------------------------------
c Full linear response (A-B)^(1/2)(A+B)(A-B)^(1/2) Z=w^2 Z
c --------------------------------------------------------
c
      else
c
c Diagonalize (A-B) (even when it is already diagonal)
c
        do i=1,ntrials
          do j=1,ntrials
            dbl_mb(k_work4+(i-1)*ntrials+j-1)=amb(j,i)
          enddo
        enddo
        call util_dgeev('N','V',ntrials,dbl_mb(k_work4),ntrials,
     1    ambval,dummy1,dummy2,1,ambvec,ntrials,
     2    dbl_mb(k_work3),10*ntrials,info)
        call tddft_schmidt(ntrials,ambvec)
        call util_nwsort(ntrials,ambval,ambvec,'A')
        if (util_print('diagonization',print_debug)) then
          if (nodezero) then
            call ma_print(ambval,ntrials,1,'Eigenvalues')
            call ma_print(ambvec,ntrials,ntrials,'Eigenvectors') 
          endif
        endif
c
c Make (A-B)^(1/2) and (A-B)^(-1/2)
c
        do i=1,ntrials*ntrials
          dbl_mb(k_work1+i-1)=0.0d0
          dbl_mb(k_work2+i-1)=0.0d0
        enddo
        do i=1,ntrials
          if (ambval(i).lt.0.0d0) then
            call errquit
     1      ('tddft_diagon: A-B is not positive-definite; try TDA',0,
     &       INPUT_ERR)
          else
            dbl_mb(k_work1+(i-1)*ntrials+i-1)=dsqrt(ambval(i))
            dbl_mb(k_work2+(i-1)*ntrials+i-1)=1.0d0/dsqrt(ambval(i))
          endif
        enddo
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,ambvec,
     1    ntrials,dbl_mb(k_work1),ntrials,0.0d0,apbvec,ntrials)
        call dgemm('N','T',ntrials,ntrials,ntrials,1.0d0,apbvec,
     1    ntrials,ambvec,ntrials,0.0d0,dbl_mb(k_work1),ntrials)
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,ambvec,
     1    ntrials,dbl_mb(k_work2),ntrials,0.0d0,apbvec,ntrials)
        call dgemm('N','T',ntrials,ntrials,ntrials,1.0d0,apbvec,
     1    ntrials,ambvec,ntrials,0.0d0,dbl_mb(k_work2),ntrials)
        if (util_print('diagonization',print_debug)) then
          call ma_print(dbl_mb(k_work1),ntrials,ntrials,'(A-B)^(1/2)')
          call ma_print(dbl_mb(k_work2),ntrials,ntrials,'(A-B)^(-1/2)')
        endif
c
c Make (A-B)^(1/2)(A+B)(A-B)^(1/2)
c
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,
     1    dbl_mb(k_work1),ntrials,apb,ntrials,0.0d0,apbvec,ntrials)
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,
     1    apbvec,ntrials,dbl_mb(k_work1),ntrials,0.0d0,ambvec,ntrials)
c
c Diagonalize (A-B)^(1/2)(A+B)(A-B)^(1/2)
c
        do i=1,ntrials
          do j=1,ntrials
            dbl_mb(k_work4+(i-1)*ntrials+j-1)=ambvec(j,i)
          enddo
        enddo
        call util_dgeev('N','V',ntrials,dbl_mb(k_work4),ntrials,
     1    apbval,dummy1,dummy2,1,apbvec,ntrials,
     2    dbl_mb(k_work3),10*ntrials,info)
        call tddft_schmidt(ntrials,apbvec)
        call util_nwsort(ntrials,apbval,apbvec,'A')
        do i=1,ntrials
          if (apbval(i).lt.0.0d0) then
            call errquit
     1      ('tddft_diagon: negative excitation energy',0, CALC_ERR)
          else
            apbval(i)=dsqrt(apbval(i))
          endif
        enddo
        if (util_print('diagonization',print_debug)) then
          if (nodezero) then
            call ma_print(apbval,ntrials,1,'Eigenvalues')
            call ma_print(apbvec,ntrials,ntrials,'Eigenvectors') 
          endif
        endif
c
c Make (X-Y) = w(A-B)^(-1/2) Z
c
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,
     1    dbl_mb(k_work2),ntrials,apbvec,ntrials,0.0d0,ambvec,ntrials)
        do i=1,ntrials
          do j=1,ntrials
            ambvec(j,i)=ambvec(j,i)*apbval(i)
          enddo
        enddo
c
c Make (X+Y) = (A-B)^(1/2) Z
c
        call dgemm('N','N',ntrials,ntrials,ntrials,1.0d0,
     1    dbl_mb(k_work1),ntrials,apbvec,ntrials,0.0d0,dbl_mb(k_work2),
     2    ntrials)
        do i=1,ntrials
          do j=1,ntrials
            apbvec(j,i)=dbl_mb(k_work2+(i-1)*ntrials+j-1)
          enddo
        enddo
c
c Bi-orthogonalize
c
        do n=1,ntrials
          do m=1,ntrials
            if (n.ne.m) then
              a=0.0d0
              do j=1,ntrials
                a=a+apbvec(j,n)*ambvec(j,m)
              enddo
              do j=1,ntrials
                apbvec(j,n)=apbvec(j,n)-ambvec(j,m)*a
              enddo
            endif
          enddo
        enddo
c
c Normalize
c
        do n=1,ntrials
          a=0.0d0
          do j=1,ntrials
            a=a+apbvec(j,n)*ambvec(j,n)
          enddo
          if (a.ne.0.d0) then
            do j=1,ntrials
              apbvec(j,n)=apbvec(j,n)/dsqrt(a)
              ambvec(j,n)=ambvec(j,n)/dsqrt(a)
            enddo
          endif ! a .ne. 0.0
        enddo
      endif
c           
c ------
c Return
c ------
c
      if (.not.ma_pop_stack(l_work4)) call errquit
     1  ('tddft_diagon: failed to deallocate work4',0, MA_ERR)
      if (.not.ma_pop_stack(l_work3)) call errquit
     1  ('tddft_diagon: failed to deallocate work3',0, MA_ERR)
      if (.not.ma_pop_stack(l_work2)) call errquit
     1  ('tddft_diagon: failed to deallocate work2',0, MA_ERR)
      if (.not.ma_pop_stack(l_work1)) call errquit
     1  ('tddft_diagon: failed to deallocate work1',0, MA_ERR)
      return
      end
