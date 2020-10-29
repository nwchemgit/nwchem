
*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   CMatrix library			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 5/19/06					*
*     *								*
*     ***********************************************************
c
c
c
c

*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm1_rot2     *
*     *                                 *
*     ***********************************

      subroutine CMatrix_zgemm1_rot2(m,n,k,
     >                  alpha,
     >                  A,lda,ma,na,
     >                  B,ldb,mb,nb,
     >                  beta,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  Bcol,Bwork,work1,work2)
      implicit none
      integer m,n,k
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)

      integer ldb,mb(*),nb(*)
c      real*8  B(ldb,*)
      complex*16  B(*)

      complex*16  beta

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  Bcol(*),Bwork(*)
      complex*16  work1(*),work2(*)

#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** local variables ****
      logical jeven
      integer i,j,j1,ii,jj,ne0
      integer iwrk,jcur,ierr,bshift,iwrk2
      integer request1(4),request2(4)
      integer bshift2(0:np_j)

      bshift     = 1
      bshift2(0) = 1
      do jj=1,np_j
         bshift      = bshift + na(jj)*nb(taskid_j+1)
         bshift2(jj) = bshift
      end do

*     *** collect B into columns ***
      ne0 = 0
      do i=1,np_i
         ne0 = ne0+mb(i)
      end do
      j1 = 0
      do jj=1,taskid_j
         j1 = j1 + nb(jj)
      end do
   
      bshift = 0
      iwrk2 = nb(taskid_j+1)
      ii    = 0
      do jcur=0,np_j-1
         iwrk  = na(jcur+1)
!$OMP DO
         do j=1,iwrk2
         do i=1,iwrk
            Bwork(bshift+i+(j-1)*iwrk) = B(ii+i+(j+j1-1)*ne0)
         end do
         end do
!$OMP END DO
         bshift = bshift + iwrk*iwrk2
         ii     = ii     + iwrk
      end do


*     *** C = beta*C ***
c      call dscal(ldc*nc(taskid_j+1),beta,C,1)
!$OMP DO
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j) = beta*C(i,j)
         end do
      end do
!$OMP END DO

      call Cmatrix_start_rot(1,
     >                       taskid_j,np_j,comm_j,
     >                       A,work1,lda,na,
     >                       request1)
      jcur = taskid_j
      iwrk = na(jcur+1)

      if ((mc(taskid_i+1).gt.0).and.
     >    (nc(taskid_j+1).gt.0).and.
     >    (iwrk.gt.0))
     >     call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                A, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

      jeven = .true.
      do j=2,np_j-1
         if (jeven) then
            jeven = .false.
            jcur = mod(jcur-1+np_j,np_j)
            iwrk = na(jcur+1)
            call Cmatrix_end_rot(request1)
            call Cmatrix_start_rot(j,
     >                             taskid_j,np_j,comm_j,
     >                             A,work2,lda,na,
     >                             request2)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                work1, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

         else
            jeven = .true.
            jcur = mod(jcur-1+np_j,np_j)
            iwrk = na(jcur+1)
            call Cmatrix_end_rot(request2)
            call Cmatrix_start_rot(j,
     >                             taskid_j,np_j,comm_j,
     >                             A,work1,lda,na,
     >                             request1)
            if ((mc(taskid_i+1).gt.0).and.
     >          (nc(taskid_j+1).gt.0).and.
     >          (iwrk.gt.0))
     >          call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),
     >                iwrk,
     >                alpha,
     >                work2, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

         end if 
      end do
      if (jeven) then
         jcur = mod(jcur-1+np_j,np_j)
         iwrk = na(jcur+1)
         call Cmatrix_end_rot(request1)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >             alpha,
     >             work1, ma(taskid_i+1),
     >             Bwork(bshift2(jcur)), iwrk,
     >             dcmplx(1.0d0,0.0d0),
     >             C, ldc)

      else
         jcur = mod(jcur-1+np_j,np_j)
         iwrk = na(jcur+1)
         call Cmatrix_end_rot(request2)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >             alpha,
     >             work2, ma(taskid_i+1),
     >             Bwork(bshift2(jcur)), iwrk,
     >             dcmplx(1.0d0,0.0d0),
     >             C, ldc)
      end if 
      return
      end





*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm1_rot      *
*     *                                 *
*     ***********************************

      subroutine CMatrix_zgemm1_rot(m,n,k,
     >                  alpha,
     >                  A,lda,ma,na,
     >                  B,ldb,mb,nb,
     >                  beta,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  Bcol,Bwork,work1,work2)
      implicit none
      integer m,n,k
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)

      integer ldb,mb(*),nb(*)
c      real*8  B(ldb,*)
      complex*16  B(*)

      complex*16  beta

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  Bcol(*),Bwork(*)
      complex*16  work1(*),work2(*)

#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** local variables ****
      logical jeven
      integer i,j,ii,jj,mbmax,nbmax,ne0
      integer iwrk,jcur,ierr,bshift,iwrk2
      integer request1(4),request2(4)
      integer bshift2(0:np_j)

      bshift     = 1
      bshift2(0) = 1
      do jj=1,np_j
         bshift      = bshift + na(jj)*nb(taskid_j+1)
         bshift2(jj) = bshift
      end do

*     *** collect B into columns ***
      ne0 = 0
      mbmax = 0
      do i=1,np_i
         if (mb(i).gt.mbmax) mbmax = mb(i)
         ne0 = ne0+mb(i)
      end do
      nbmax = 0
      do j=1,np_j
         if (nb(j).gt.nbmax) nbmax = nb(j)
      end do
      if (np_i.gt.1) then


         do j=1,nb(taskid_j+1)


          !*** Allgather is flaky on my mac laptop ***
          !call dcopy(ne0,0.0,0,Bcol(1+(j-1)*ne0),1)
          call Parallel_shared_vector_zero(.true.,ne0,Bcol(1+(j-1)*ne0))
          bshift = 0
          do ii=1,taskid_i
             bshift = bshift + mb(ii)
          end do
          do i=1,mb(taskid_i+1)
            Bcol(bshift+i+(j-1)*ne0) = B(i+(j-1)*mb(taskid_i+1))
          end do
          call C3dB_Vector_SumAll(2*ne0,Bcol(1+(j-1)*ne0))
         end do

         bshift = 0
         iwrk2 = nb(taskid_j+1)
         ii    = 0
         do jcur=0,np_j-1
            iwrk  = na(jcur+1)
!$OMP DO
            do j=1,iwrk2
            do i=1,iwrk
               Bwork(bshift+i+(j-1)*iwrk) = Bcol(ii+i+(j-1)*ne0)

            end do
            end do
!$OMP END DO
            bshift = bshift + iwrk*iwrk2
            ii     = ii     + iwrk
         end do
      else
         bshift = 0
         iwrk2 = nb(taskid_j+1)
         ii    = 0
         do jcur=0,np_j-1
            iwrk  = na(jcur+1)
!$OMP DO
            do j=1,iwrk2
            do i=1,iwrk
               Bwork(bshift+i+(j-1)*iwrk) = B(ii+i+(j-1)*ne0)
            end do
            end do
!$OMP END DO
            bshift = bshift + iwrk*iwrk2
            ii     = ii     + iwrk
         end do
      end if



*     *** C = beta*C ***
c      call dscal(ldc*nc(taskid_j+1),beta,C,1)
!$OMP DO
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j) = beta*C(i,j)
         end do
      end do
!$OMP END DO

      call Cmatrix_start_rot(1,
     >                       taskid_j,np_j,comm_j,
     >                       A,work1,lda,na,
     >                       request1)
      jcur = taskid_j
      iwrk = na(jcur+1)

      if ((mc(taskid_i+1).gt.0).and.
     >    (nc(taskid_j+1).gt.0).and.
     >    (iwrk.gt.0))
     >     call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                A, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

      jeven = .true.
      do j=2,np_j-1
         if (jeven) then
            jeven = .false.
            jcur = mod(jcur-1+np_j,np_j)
            iwrk = na(jcur+1)
            call Cmatrix_end_rot(request1)
            call Cmatrix_start_rot(j,
     >                             taskid_j,np_j,comm_j,
     >                             A,work2,lda,na,
     >                             request2)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                work1, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

         else
            jeven = .true.
            jcur = mod(jcur-1+np_j,np_j)
            iwrk = na(jcur+1)
            call Cmatrix_end_rot(request2)
            call Cmatrix_start_rot(j,
     >                             taskid_j,np_j,comm_j,
     >                             A,work1,lda,na,
     >                             request1)
            if ((mc(taskid_i+1).gt.0).and.
     >          (nc(taskid_j+1).gt.0).and.
     >          (iwrk.gt.0))
     >          call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),
     >                iwrk,
     >                alpha,
     >                work2, ma(taskid_i+1),
     >                Bwork(bshift2(jcur)), iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)

         end if 
      end do
      if (jeven) then
         jcur = mod(jcur-1+np_j,np_j)
         iwrk = na(jcur+1)
         call Cmatrix_end_rot(request1)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >             alpha,
     >             work1, ma(taskid_i+1),
     >             Bwork(bshift2(jcur)), iwrk,
     >             dcmplx(1.0d0,0.0d0),
     >             C, ldc)

      else
         jcur = mod(jcur-1+np_j,np_j)
         iwrk = na(jcur+1)
         call Cmatrix_end_rot(request2)
         if ((mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0).and.
     >       (iwrk.gt.0))
     >       call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >             alpha,
     >             work2, ma(taskid_i+1),
     >             Bwork(bshift2(jcur)), iwrk,
     >             dcmplx(1.0d0,0.0d0),
     >             C, ldc)
      end if 
      return
      end


      subroutine CMatrix_Bwork_copy(m,n,A,lda,ii,B,ldb)
      implicit none
      integer n,m
      integer lda,ii,ldb
      complex*16 A(LDA,*)
      complex*16 B(LDB,*)

*     *** local variables ***
      integer i,j

      do j=1,n
         do i=1,m
            B(i,j) = A(ii+i,j)
         end do
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *         CMatrix_start_rot       *
*     *                                 *
*     ***********************************

      subroutine CMatrix_start_rot(j,
     >                             taskid_j,np_j,comm_j,
     >                             A,W,lda,na,
     >                             request)
      implicit none
      integer j
      integer taskid_j,np_j,comm_j
      complex*16 A(*),W(*)
      integer lda,na(*)
      integer request(*)

#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer proc_to,proc_from,msgtype,amsglen,wmsglen,mpierr

      proc_to   = mod(taskid_j+j,np_j)
      proc_from = mod(taskid_j-j+np_j,np_j)
      msgtype   = j
      amsglen = lda*na(taskid_j+1)
      wmsglen = lda*na(proc_from+1)
!$OMP MASTER
#ifdef MPI4
            if (wmsglen.gt.0) then
               stupid_msglen = wmsglen
               stupid_type   = msgtype
               stupid_taskid = proc_from
               call MPI_IRECV(W,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_j,
     >                    stupid_request,stupid_ierr)
               request(1) = stupid_request
               request(3) = 1
            else
               request(3) = 0
            end if

            if (amsglen.gt.0) then
               stupid_msglen = amsglen
               stupid_type   = msgtype
               stupid_taskid = proc_to
               call MPI_ISEND(A,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_j,
     >                     stupid_request,stupid_ierr)
               request(2) = stupid_request
               request(4) = 1
            else
               request(4) = 0
            end if
#else
            if (wmsglen.gt.0) then
               call MPI_IRECV(W,wmsglen,MPI_DOUBLE_COMPLEX,
     >                    proc_from,
     >                    msgtype,comm_j,
     >                    request(1),mpierr)
               request(3) = 1
            else
               request(3) = 0
            end if
            if (amsglen.gt.0) then
               call MPI_ISEND(A,amsglen,MPI_DOUBLE_COMPLEX,
     >                     proc_to,
     >                     msgtype,comm_j,
     >                     request(2),mpierr)
               request(4) = 1
            else
               request(4) = 0
            end if
#endif

      if ((request(3).eq.1).and.(request(4).eq.1)) then
         request(3) = 1
      else if (request(3).eq.1) then
         request(3) = 2
      else if (request(4).eq.1) then
         request(3) = 3
      else 
         request(3) = 4
      end if
!$OMP END MASTER

      return
      end

*     ***********************************
*     *                                 *
*     *         CMatrix_end_rot         *
*     *                                 *
*     ***********************************

      subroutine CMatrix_end_rot(request)
      implicit none
      integer request(*)

*     **** wait for completion of mp_send, also do a sync ****
      if (request(3).eq.1) then
         call Parallel_mpiWaitAll(2,request)
      else if (request(3).eq.2) then
         call Parallel_mpiWaitAll(1,request)
      else if (request(3).eq.3) then
         call Parallel_mpiWaitAll(1,request(2))
      endif

      return
      end





*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm1          *
*     *                                 *
*     ***********************************

      subroutine CMatrix_zgemm1(m,n,k,nblock,
     >                  alpha,
     >                  A,lda,ma,na,
     >                  B,ldb,mb,nb,
     >                  beta,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer m,n,k,nblock
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)

      integer ldb,mb(*),nb(*)
      complex*16  B(ldb,*)

      complex*16  beta

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  work1(*),work2(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift


      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j) = beta*C(i,j)
         end do
      end do

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.k)
         iwrk = min(nblock, mb(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)

!$OMP MASTER
*        **** pack current iwrk columns of A into work1 ***
         if (taskid_j.eq.jcur) then
            call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
         end if

*        **** pack current iwrk rows of B into work2 ***
         if (taskid_i.eq.icur) then
            call zlacpy("G", iwrk,nb(taskid_j+1),
     >                   B(ii+1,1), ldb,
     >                   work2,  iwrk)
         end if

#ifdef MPI4
c        **** broadcast work1  within my row ***
         stupid_msglen = iwrk*ma(taskid_i+1)
         stupid_taskid = jcur
         call MPI_Bcast(work1,stupid_msglen,stupid_complex,
     >                  stupid_taskid,stupid_comm_j,stupid_ierr)

c        **** broadcast work2  within my column ***
         stupid_msglen = iwrk*nb(taskid_j+1)
         stupid_taskid = icur
         call MPI_Bcast(work2,stupid_msglen,stupid_complex,
     >                  stupid_taskid,stupid_comm_i,stupid_ierr)
#else
c        **** broadcast work1  within my row ***
         call MPI_Bcast(work1,iwrk*ma(taskid_i+1),MPI_DOUBLE_COMPLEX,
     >                  jcur,comm_j,ierr)

c        **** broadcast work2  within my column ***
         call MPI_Bcast(work2,iwrk*nb(taskid_j+1),MPI_DOUBLE_COMPLEX,
     >                  icur,comm_i,ierr)
#endif

!$OMP END MASTER
!$OMP BARRIER


         if ((iwrk.gt.0)          .and.
     >       (mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0))
     >     call zgemm_omp('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                work1, ma(taskid_i+1),
     >                work2, iwrk,
     >                dcmplx(1.0d0,0.0d0),
     >                C, ldc)


         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mb(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do

      return
      end


*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm2          *
*     *                                 *
*     ***********************************
      subroutine CMatrix_zgemm2(m,n,k,nblock,
     >                  alpha,
     >                  A,lda,ma,na,
     >                  B,ldb,mb,nb,
     >                  beta,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer m,n,k,nblock
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)

      integer ldb,mb(*),nb(*)
      complex*16  B(ldb,*)

      complex*16  beta

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  work1(*),work2(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift

!$OMP DO
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j) = beta*C(i,j)
         end do
      end do
!$OMP END DO

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)

*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_complex,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_COMPLEX,jcur,comm_j,ierr)
#endif

!$OMP END MASTER
!$OMP BARRIER


c            if ((iwrk.gt.0)          .and.
c     >          (nb(taskid_j+1).gt.0).and.
c     >          (ma(taskid_i+1).gt.0))
            if ((iwrk.gt.0)          .and.
     >          (nb(taskid_j+1).gt.0))
     >        call zgemm_omp('C','N',iwrk,nb(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, ldb,
     >                   dcmplx(0.0d0,0.0d0),
     >                   work2, iwrk)

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            !call dcopy(2*nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
            call Parallel_shared_vector_zero(.true.,
     >                                       2*nc(taskid_j+1)*iwrk,
     >                                       work2)
         end if


c        **** summ to node that holds current rows of C ****
!$OMP MASTER
#ifdef MPI4
         stupid_msglen = nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_complex,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_COMPLEX,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call daxpy_omp(2*nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do


      return
      end

*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm2c_omp     *
*     *                                 *
*     ***********************************
      subroutine CMatrix_zgemm2c_omp(m,n,k,nblock,
     >                  A,B,lda,ma,na,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  tid,nthr,
     >                  work1,work2,
     >                  privatetmp)
      implicit none
      integer m,n,k,nblock

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)
      complex*16  B(lda,*)

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      integer tid,nthr

      complex*16 work1(*),work2(*)
      complex*16 privatetmp(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      integer npack2,n2

      n2     = na(taskid_j+1)
      npack2 = ma(taskid_i+1)

      call Parallel_shared_vector_zero(.false.,
     >                                 2*mc(taskid_i+1)*nc(taskid_j+1),
     >                                 C)

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)

         call Parallel_shared_vector_zero(.false.,
     >                                    2*nc(taskid_j+1)*iwrk,work2)
*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1, ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_complex,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_COMPLEX,jcur,comm_j,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER

            if ((iwrk.gt.0).and.(na(taskid_j+1).gt.0)) then
              call zgemm2c_omp_group(1,tid,nthr,iwrk,n2,npack2,
     >                       work1,B,work2,privatetmp)
            end if
         end if

c        **** summ to node that holds current rows of C ****
!$OMP BARRIER
!$OMP MASTER
#ifdef MPI4
         stupid_msglen = nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_complex,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_COMPLEX,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call DAXPY_OMP(2*nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do


      return
      end


*     ***********************************
*     *                                 *
*     *         CMatrix_dgemm2c_omp     *
*     *                                 *
*     ***********************************
      subroutine CMatrix_dgemm2c_omp(m,n,k,nblock,
     >                  A,B,lda,ma,na,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  tid,nthr,
     >                  work1,work2,
     >                  privatetmp)
      implicit none
      integer m,n,k,nblock

      integer lda,ma(*),na(*)
      real*8  A(lda,*)
      real*8  B(lda,*)

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      integer tid,nthr

      real*8  work1(*),work2(*)
      real*8  privatetmp(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif



*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      integer npack2,n2

      n2     = na(taskid_j+1)
      npack2 = ma(taskid_i+1)

      call Parallel_shared_vector_zero(.false.,
     >                                 mc(taskid_i+1)*nc(taskid_j+1),
     >                                 C)

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)


         call Parallel_shared_vector_zero(.false.,
     >                                    nc(taskid_j+1)*iwrk,work2)
*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1, ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_double,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_PRECISION,jcur,comm_j,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


            if ((iwrk.gt.0).and.(na(taskid_j+1).gt.0)) then
              call dgemm2c0_omp_group(1,tid,nthr,iwrk,n2,npack2,
     >                       work1,B,work2,privatetmp)
            end if
         end if

c        **** summ to node that holds current rows of C ****
!$OMP BARRIER
!$OMP MASTER
#ifdef MPI4
         stupid_msglen = nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_double,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call DAXPY_OMP(nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do


      return
      end




*     ***********************************
*     *                                 *
*     *         CMatrix_zgemm3          *
*     *                                 *
*     ***********************************

      subroutine CMatrix_zgemm3(m,n,k,nblock,
     >                  alpha,
     >                  A,lda,ma,na,
     >                  B,ldb,mb,nb,
     >                  beta,
     >                  C,ldc,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer m,n,k,nblock
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)

      integer ldb,mb(*),nb(*)
      complex*16  B(ldb,*)

      complex*16  beta

      integer ldc,mc(*),nc(*)
      complex*16  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  work1(*),work2(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      real*8  dum

!$OMP DO
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j) = beta*C(i,j)
         end do
      end do
!$OMP END DO

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
      do while (kk.lt.n)
         iwrk = min(nblock, mb(icur+1)-ii)
         iwrk = min(iwrk,   nc(jcur+1)-jj)


!$OMP MASTER
         if (taskid_i.eq.icur) then
            call zlacpy("G", iwrk,nb(taskid_j+1),
     >                   B(ii+1,1), ldb,
     >                   work2,     iwrk)
         end if

#ifdef MPI4
        stupid_msglen = iwrk*nb(taskid_j+1)
        stupid_taskid = icur
        call MPI_Bcast(work2,stupid_msglen,stupid_complex,
     >                  stupid_taskid,stupid_comm_i,stupid_ierr)
#else
        call MPI_Bcast(work2,iwrk*nb(taskid_j+1),MPI_DOUBLE_COMPLEX,
     >                  icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER

         if ((iwrk.gt.0)          .and.
     >       (na(taskid_j+1).gt.0).and.
     >       (mc(taskid_i+1).gt.0))
     >      call zgemm_omp('N','C',mc(taskid_i+1),iwrk,na(taskid_j+1),
     >              alpha,
     >              A, lda,
     >              work2, iwrk,
     >              dcmplx(0.0d0,0.0d0),
     >              work1, mc(taskid_i+1))

#ifdef MPI4
        stupid_msglen = mc(taskid_i+1)*iwrk
        stupid_taskid = jcur
        call MPI_Reduce(work1,work2,stupid_msglen,stupid_complex,
     >                  stupid_sum,stupid_taskid,
     >                  stupid_comm_j,stupid_ierr)
#else
        call MPI_Reduce(work1,work2,mc(taskid_i+1)*iwrk,
     >                   MPI_DOUBLE_COMPLEX,MPI_SUM,jcur,comm_j,ierr)
#endif


         if (taskid_j.eq.jcur) then
            shift = 1
            do j=jj,(jj+iwrk-1)
               call daxpy_omp(2*mc(taskid_i+1),
     >                    1.0d0,
     >                    work2(shift),1,
     >                    C(1,j+1),1)
               shift = shift + mc(taskid_i+1)
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.nc(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mb(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do


      return
      end





*     ***********************************
*     *                                 *
*     *         CMatrix_tqliq           *
*     *                                 *
*     ***********************************

      subroutine CMatrix_tqliq(n,eig,tu,
     >                  Q,ldq,mq,nq,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer n

      integer ldq,mq(*),nq(*)
      complex*16  Q(ldq,*)
      real*8  eig(*),tu(*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      complex*16  work1(*),work2(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      integer MAXITER
      parameter (MAXITER = 100)
      real*8  tole
      parameter (tole=1.0d-15)

      logical notdone
      integer i,j,l,m,iter
      integer ii,jj0,jj1,jcur0,jcur1,ierr,istat
      real*8  b,c,f,g,p,r,s


      do l=1,n-1
         iter = 0

         do m=l,n-1
         if (dabs(tu(m)).lt.tole) go to 2
         end do
         m = n
  2      continue
         if (m.eq.l) then
            notdone = .false.
         else
            notdone = .true.
         end if
         do while ((iter.lt.MAXITER).and.(notdone))
            g = (eig(l+1)-eig(l))/(2.0d0*tu(l))
            r = dsqrt(g**2+1.0d0)
ccccneed to fixccc            g = eig(m)-eig(l)+tu(l)/(g+dsign(r,g))
            s = 1.0d0
            c = 1.0d0
            p = 0.0d0
            do i = m-1,l,-1
               f = s*tu(i)
               b = c*tu(i)
               if (dabs(f).ge.dabs(g)) then
                  c = g/f
                  r = dsqrt(c**2+1.0d0)
                  tu(i+1) = f*r
                  s = 1/r
                  c = c*s
               else
                  s = f/g
                  r = dsqrt(s**2+1.0d0)
                  tu(i+1) = g*r
                  c = 1/r
                  s = s*c
               end if
               g = eig(i+1)-p
               r = (eig(i)-g)*s + 2.0d0*c*b
               p = s*r
               eig(i+1) = g+p
               g = c*r-b


*              **** update eigenvectors ****
               jcur0 = 0
               jj0   = 1
               do j=1,i-1
                 jj0 = jj0 + 1
                 if (jj0.gt.nq(jcur0+1)) then
                    jcur0 = jcur0 + 1
                    jj0   = 1
                 end if
               end do
               jcur1 = jcur0
               jj1   = jj0 + 1
               if (jj1.gt.nq(jcur1+1)) then
                  jcur1 = jcur1 + 1
                  jj1   = 1
               end if

               if (jcur0.eq.taskid_j)
     >             call zcopy(mq(taskid_i+1),Q(1,jj0),1,work1,1)
               if (jcur1.eq.taskid_j)
     >             call zcopy(mq(taskid_i+1),Q(1,jj1),1,work2,1)

#ifdef MPI4
               stupid_msglen = mq(taskid_i+1)
               stupid_taskid = jcur0
               stupid_type   = jcur1
               call MPI_Bcast(work1,stupid_msglen,stupid_complex,
     >                  stupid_taskid,stupid_comm_j,stupid_ierr)
               call MPI_Bcast(work2,stupid_msglen,stupid_complex,
     >                  stupid_type,stupid_comm_j,stupid_ierr)
#else
               call MPI_Bcast(work1,mq(taskid_i+1),MPI_DOUBLE_COMPLEX,
     >                  jcur0,comm_j,ierr)
               call MPI_Bcast(work2,mq(taskid_i+1),MPI_DOUBLE_COMPLEX,
     >                  jcur1,comm_j,ierr)
#endif

               if (jcur0.eq.taskid_j) then
                  do ii=1,mq(taskid_i+1)
                    Q(ii,jj0) = c*Q(ii,jj0) - s*work2(ii)
                  end do
               end if

               if (jcur1.eq.taskid_j) then
                  do ii=1,mq(taskid_i+1)
                    Q(ii,jj1) = c*Q(ii,jj1) + s*work1(ii)
                  end do
               end if

            end do
            eig(l) = eig(l) - p
            tu(l)  = g
            tu(m)  = 0.0d0


            do m=l,n-1
            if (dabs(tu(m)).lt.tole) go to 3
            end do
            m = n
  3         continue
            if (m.eq.l) then
               notdone = .false.
            else
               notdone = .true.
            end if

            iter = iter + 1
         end do

      end do

      return
      end



*     ***********************************
*     *                                 *
*     *         CMatrix_houseq          *
*     *                                 *
*     ***********************************
      subroutine CMatrix_houseq(jcol,
     >                  n,
     >                  A,V,Q,lda,ma,na,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer jcol,n

      integer lda,ma(*),na(*)
      complex*16  A(lda,*),V(lda,*),Q(lda,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  work1(*),work2(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      integer ii0,icur0,ii1,icur1,ii2,icur2
      integer jj0,jcur0,jj1,jcur1
      real*8  beta,mu0,mu,v20,v2

      !call dcopy(2*ma(taskid_i+1)*na(taskid_j+1),0.0d0,0,V,1)
      call Parallel_shared_vector_zero(.true.,
     >                                 2*ma(taskid_i+1)*na(taskid_j+1),
     >                                 V)

      jcur0 = 0
      jj0   = 1
      do j=1,jcol-1
        jj0 = jj0 + 1
        if (jj0.gt.na(jcur0+1)) then
           jcur0 = jcur0 + 1
           jj0 = 1
        end if
      end do
      jcur1 = jcur0
      jj1   = jj0 + 1
      if (jj1.gt.na(jcur1+1)) then
           jcur1 = jcur1 + 1
           jj1 = 1
      end if

      icur0 = 0
      ii0   = 1
      do i=1,jcol-1
        ii0 = ii0 + 1
        if (ii0.gt.ma(icur0+1)) then
           icur0 = icur0 + 1
           ii0 = 1
        end if
      end do
      icur1 = icur0
      ii1   = ii0 + 1
      if (ii1.gt.ma(icur1+1)) then
           icur1 = icur1 + 1
           ii1 = 1
      end if
      icur2 = icur1
      ii2   = ii1 + 1
      if (ii2.gt.ma(icur2+1)) then
           icur2 = icur2 + 1
           ii2 = 1
      end if

      if (jcur0.eq.taskid_j) then

         icur = icur1
         ii   = ii1
         do i=jcol+1,n
            if (icur.eq.taskid_i) V(ii,jj0) = A(ii,jj0)
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do


         mu0 = 0.0d0
         icur = icur1
         ii   = ii1
         do i=jcol+1,n
            if (icur.eq.taskid_i) 
     >         mu0 = mu0 + dconjg(V(ii,jj0))*V(ii,jj0)
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do
#ifdef MPI4
         stupid_msglen = 1
         call MPI_AllReduce(mu0,mu,stupid_msglen,
     >                      stupid_double,stupid_sum,
     >                      stupid_comm_i,stupid_ierr)
#else
         call MPI_AllReduce(mu0,mu,1,
     >                      MPI_DOUBLE_PRECISION,MPI_SUM,comm_i,ierr)
#endif
         mu = dsqrt(mu)


         if (mu.ne.0.0d0) then
cccc need to fix           if (icur1.eq.taskid_i)
cccc need to fix   >        beta = V(ii1,jj0) + dsign(mu,V(ii1,jj0))
#ifdef MPI4
           stupid_msglen = 1
           stupid_taskid = icur1
           call MPI_Bcast(beta,stupid_msglen,stupid_double,
     >                  stupid_taskid,stupid_comm_i,stupid_ierr)
#else
           call MPI_Bcast(beta,1,MPI_DOUBLE_PRECISION,
     >                  icur1,comm_i,ierr)
#endif

           icur = icur2
           ii   = ii2
           do i=jcol+2,n
              if (icur.eq.taskid_i) V(ii,jj0) = V(ii,jj0)/beta
              ii = ii + 1
              if (ii.gt.ma(icur+1)) then
                 icur = icur + 1
                 ii = 1
              end if
           end do
         end if
         if (icur1.eq.taskid_i)  V(ii1,jj0) = dcmplx(1.0d0,0.0d0)
         if (icur0.eq.taskid_i)  V(ii0,jj0) = dcmplx(0.0d0,0.0d0)

         v20 = 0.0d0
         icur = icur0
         ii   = ii0
         do i=jcol,n
            if (icur.eq.taskid_i) 
     >         v20 = v20 + dcmplx(V(ii,jj0))*V(ii,jj0)
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do
#ifdef MPI4
         stupid_msglen = 1
         call MPI_AllReduce(v20,v2,stupid_msglen,
     >                      stupid_double,stupid_sum,
     >                      stupid_comm_i,stupid_ierr)
#else
         call MPI_AllReduce(v20,v2,1,
     >                      MPI_DOUBLE_PRECISION,MPI_SUM,comm_i,ierr)
#endif
         v2 = 2.0d0/v2
      end if
#ifdef MPI4
      stupid_msglen = 1
      stupid_taskid = jcur0
      call MPI_Bcast(v2,stupid_msglen,stupid_double,
     >               stupid_taskid,stupid_comm_j,stupid_ierr)
#else
      call MPI_Bcast(v2,1,MPI_DOUBLE_PRECISION,
     >               jcur0,comm_j,ierr)
#endif

      call CMatrix_eye(n,n,dcmplx(1.0d0,0.0d0),
     >                 Q,lda,ma,na,taskid_i,taskid_j)
      call CMatrix_zgemm3(n,n,n,64,
     >             dcmplx(-v2,0.0d0),
     >             V,ma(taskid_i+1), ma,na,
     >             V,ma(taskid_i+1), ma,na,
     >             dcmplx(1.0d0,0.0d0),
     >             Q,ma(taskid_i+1), ma,na,
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             work1,work2)


      return
      end


*     ***********************************
*     *                                 *
*     *         CMatrix_eigsrtq         *
*     *                                 *
*     ***********************************
      subroutine CMatrix_eigsrtq(n,eig,
     >                  Q,ldq,mq,nq,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer n

      integer ldq,mq(*),nq(*)
      real*8     eig(*)
      complex*16 Q(ldq,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      complex*16 work1(*),work2(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical notdone
      integer i,j,k,l,m,iter
      integer ii,jj0,jj1,jcur0,jcur1,ierr,istat
      real*8  b,c,f,g,p,r,s


      do i=1,n-1
         k = i
         p = eig(i)
         do j=i+1,n
            if (eig(j).ge.p) then
               k = j
               p = eig(j)
            end if
         end do
         if (k.ne.i) then
            eig(k) = eig(i)
            eig(i) = p

            jcur0 = 0
            jj0   = 1
            do j=1,i-1
               jj0 = jj0 + 1
               if (jj0.gt.nq(jcur0+1)) then
                  jcur0 = jcur0 + 1
                  jj0 = 1
               end if
            end do
            jcur1 = 0
            jj1   = 1
            do j=1,k-1
               jj1 = jj1 + 1
               if (jj1.gt.nq(jcur1+1)) then
                  jcur1 = jcur1 + 1
                  jj1 = 1
               end if
            end do


            if (jcur0.eq.taskid_j)
     >         call zcopy(mq(taskid_i+1),Q(1,jj0),1,work1,1)
            if (jcur1.eq.taskid_j)
     >         call zcopy(mq(taskid_i+1),Q(1,jj1),1,work2,1)

#ifdef MPI4
            stupid_msglen = mq(taskid_i+1)
            stupid_taskid = jcur0
            stupid_type   = jcur1
            call MPI_Bcast(work1,stupid_msglen,stupid_complex,
     >                     stupid_taskid,stupid_comm_j,stupid_ierr)
            call MPI_Bcast(work2,stupid_msglen,stupid_complex,
     >                     stupid_type,stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,mq(taskid_i+1),MPI_DOUBLE_COMPLEX,
     >                     jcur0,comm_j,ierr)
            call MPI_Bcast(work2,mq(taskid_i+1),MPI_DOUBLE_COMPLEX,
     >                     jcur1,comm_j,ierr)
#endif

            if (jcur0.eq.taskid_j)
     >         call zcopy(mq(taskid_i+1),work2,1,Q(1,jj0),1)
            if (jcur1.eq.taskid_j)
     >         call zcopy(mq(taskid_i+1),work1,1,Q(1,jj1),1)

         end if

      end do

      return
      end


*     ***********************************
*     *                                 *
*     *         CMatrix_getdiags        *
*     *                                 *
*     ***********************************
      subroutine CMatrix_getdiags(n,eig,tu,
     >                  A,lda,ma,na,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i,comm_j,
     >                  work1)
      implicit none
      integer n

      integer lda,ma(*),na(*)
      real*8     eig(*),tu(*)
      complex*16 A(lda,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      real*8 work1(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif



*     **** local variables ****
      integer i,j,ii,jj,is,ie,js,je
      integer icur,jcur,ierr

*     **************************
*     **** gather diagonals ****
*     **************************
      !call dcopy(n,0.0d0,0,work1,1)
      !call dcopy(n,0.0d0,0,eig,1)
      call Parallel_shared_vector_zero(.false.,n,work1)
      call Parallel_shared_vector_zero(.true.,n,eig)
      js = 1
      do jcur = 0,taskid_j-1
        js = js + na(jcur+1)
      end do
      jcur = taskid_j
      je   = js-1 + na(jcur+1)
      jj   = 1
      do j=js,je

         icur=0
         ii = 1
         do i=1,j-1
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do
         work1(j) = dble(A(ii,jj))

#ifdef MPI4
         stupid_msglen = 1
         stupid_taskid = icur
         call MPI_Bcast(work1(j),stupid_msglen,stupid_double,
     >                  stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Bcast(work1(j),1,MPI_DOUBLE_PRECISION,
     >                  icur,comm_i,ierr)
#endif
         jj = jj + 1
         if (jj.gt.na(jcur+1)) then
            jcur = jcur + 1
            jj = 1
         end if
      end do

#ifdef MPI4
      stupid_msglen = n
      call MPI_AllReduce(work1,eig,stupid_msglen,
     >                   stupid_double,stupid_sum,
     >                   stupid_comm_j,stupid_ierr)
#else
      call MPI_AllReduce(work1,eig,n,
     >                   MPI_DOUBLE_PRECISION,MPI_SUM,comm_j,ierr)
#endif


*     ******************************
*     **** gather off-diagonals ****
*     ******************************
      !call dcopy(n,0.0d0,0,work1,1)
      !call dcopy(n,0.0d0,0,tu,1)
      call Parallel_shared_vector_zero(.false.,n,work1)
      call Parallel_shared_vector_zero(.true.,n,tu)
      is = 1
      do icur = 0,taskid_i-1
        is = is + ma(icur+1)
      end do
      icur = taskid_i
      ie   = is-1 + ma(icur+1)
      if (ie.ge.n) ie=ie-1
      ii   = 1
      do i=is,ie

         jcur=0
         jj = 1
         do j=1,i
            jj = jj + 1
            if (jj.gt.na(jcur+1)) then
               jcur = jcur + 1
               jj = 1
            end if
         end do
         work1(i) = dble(A(ii,jj))
#ifdef MPI4
         stupid_msglen = 1
         stupid_taskid = jcur
         call MPI_Bcast(work1(i),stupid_msglen,stupid_double,
     >                  stupid_taskid,stupid_comm_j,stupid_ierr)
#else
         call MPI_Bcast(work1(i),1,MPI_DOUBLE_PRECISION,
     >                  jcur,comm_j,ierr)
#endif
         ii = ii + 1
         if (ii.gt.ma(icur+1)) then
            icur = icur + 1
            ii = 1
         end if
      end do
#ifdef MPI4
      stupid_msglen = n-1
      call MPI_AllReduce(work1,tu,stupid_msglen,
     >                   stupid_double,stupid_sum,
     >                   stupid_comm_i,stupid_ierr)
#else
      call MPI_AllReduce(work1,tu,n-1,
     >                   MPI_DOUBLE_PRECISION,MPI_SUM,comm_i,ierr)
#endif
      return
      end




*     ***********************************
*     *                                 *
*     *         CMatrix_MaxAll          *
*     *                                 *
*     ***********************************
      subroutine CMatrix_MaxAll(sum)
c     implicit none
      complex*16  sum

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      complex*16  sumall

*     **** external functions ****

      call Parallel_np(np)
      if (np.gt.1) then
        msglen = 1
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_complex,
     >                      stupid_max,stupid_world,stupid_ierr)
#else
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_COMPLEX,
     >                      MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif
        sum = sumall
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *         CMatrix_SumAll          *
*     *                                 *
*     ***********************************
      subroutine CMatrix_SumAll(sum)
c     implicit none
      complex*16  sum

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      complex*16 sumall

*     **** external functions ****

      call Parallel_np(np)
      if (np.gt.1) then
        msglen = 1
#ifdef MPI4
        stupid_msglen = msglen
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_complex,
     >                      stupid_sum,stupid_world,stupid_ierr)
#else
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_COMPLEX,
     >                      MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        sum = sumall
      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         CMatrix_mm_transpose     *
*     *                                 *
*     ***********************************

      subroutine CMatrix_mm_transpose(n,A,B,ldq,mq,nq)
      implicit none
      integer n
      integer ldq,mq(*),nq(*)
      complex*16  A(ldq,*)
      complex*16  B(ldq,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      logical value
      integer taskid
      integer i,j
      integer ii,jj,rr,ss
      integer icur,jcur,rcur,scur
      integer psend,precv,msglen,msgtype,mpierr,status(2)

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij


      call Parallel_taskid(taskid)
      msglen  = 1
      msgtype = 1

*    **** allocate memory ****
      value = BA_push_get(mt_int,MPI_STATUS_SIZE,
     >                     'status',status(2),status(1))
      if (.not. value)
     > call errquit(' CMatrix_m_transpose:out of stack',0,MA_ERR)


      jj   = 1
      jcur = 0
      rr   = 1
      rcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         ss   = 1
         scur = 0
         do i=1,n

          
            psend = Parallel2d_convert_taskid_ij(icur,jcur)
            precv = Parallel2d_convert_taskid_ij(rcur,scur)

            if (psend.eq.precv) then
               if (psend.eq.taskid) B(rr,ss) = A(ii,jj)
            else
#ifdef MPI4
               !**** send ****
               if (psend.eq.taskid) then
                  stupid_msglen = 1
                  stupid_type   = msgtype
                  stupid_taskid = precv
                  call MPI_SEND(A(ii,jj),
     >                          stupid_msglen,stupid_complex,
     >                          stupid_taskid,
     >                          stupid_type,stupid_world,stupid_ierr)
               end if

               !**** recv ****
               if (precv.eq.taskid) then
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = psend
                  call MPI_RECV(B(rr,ss),
     >                       stupid_msglen,stupid_complex,
     >                       stupid_taskid,
     >                       stupid_type,stupid_world,
     >                       int_mb(status(1)),stupid_ierr)
               end if

#else
               if (psend.eq.taskid)
     >         call MPI_SEND(A(ii,jj),1,
     >                       MPI_DOUBLE_COMPLEX,precv,
     >                       msgtype,MPI_COMM_WORLD,mpierr)

               !**** recv ****
               if (precv.eq.taskid)
     >         call MPI_RECV(B(rr,ss),msglen,
     >                       MPI_DOUBLE_COMPLEX,psend,
     >                       msgtype,MPI_COMM_WORLD,
     >                       int_mb(status(1)),mpierr)
#endif
            end if
           

            ii = ii + 1
            if (ii.gt.mq(icur+1)) then
              icur = icur + 1
              ii   = 1
            end if

            ss = ss + 1
            if (ss.gt.nq(scur+1)) then
              scur = scur + 1
              ss   = 1
            end if

         end do

         jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
         end if

         rr = rr + 1
         if (rr.gt.mq(rcur+1)) then
           rcur = rcur + 1
           rr   = 1
         end if

      end do

      !*** deallocate memory ***
      value = BA_pop_stack(status(2))
      if (.not. value)
     > call errquit(' CMatrix_m_transpose:popping stack',0,MA_ERR)

      return
      end






*     ***********************************
*     *                                 *
*     *       CMatrix_combo_zgemm2      *
*     *                                 *
*     ***********************************
      subroutine CMatrix_combo_zgemm2(m,n,k,nblock,
     >                  alpha,
     >                  A,B,lda,ma,na,
     >                  beta,
     >                  C,ldc,ldc2,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer m,n,k,nblock
      complex*16  alpha

      integer lda,ma(*),na(*)
      complex*16  A(lda,*)
      complex*16  B(lda,*)

      complex*16  beta

      integer ldc,ldc2,mc(*),nc(*)
      complex*16  C(ldc,ldc2,3)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      complex*16  work1(*),work2(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift,shft,shft2,shft3

!$OMP DO
      do kk=1,3
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j,kk) = beta*C(i,j,kk)
         end do
      end do
      end do
!$OMP END DO

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)


*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

            shft  = iwrk*ma(taskid_i+1)
            shft2 = iwrk*nc(taskid_j+1)
            shft3 = shft2+shft2

!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))

               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   B(1,jj+1), lda,
     >                   work1(1+shft),  ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = 2*iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_complex,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,2*iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_COMPLEX,jcur,comm_j,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


            if ((iwrk.gt.0)          .and.
     >          (na(taskid_j+1).gt.0)) then

              call zgemm_omp('C','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   A, lda,
     >                   dcmplx(0.0d0,0.0d0),
     >                   work2, iwrk)
              call zgemm_omp('C','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, lda,
     >                   dcmplx(0.0d0,0.0d0),
     >                   work2(1+shft2), iwrk)
              call zgemm_omp('C','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1(1+shft), ma(taskid_i+1),
     >                   B, lda,
     >                   dcmplx(0.0d0,0.0d0),
     >                   work2(1+shft3), iwrk)
            end if

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            !call dcopy(3*nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
            call Parallel_shared_vector_zero(.true.,
     >                                 3*nc(taskid_j+1)*iwrk,work2)
         end if


!$OMP MASTER
c        **** summ to node that holds current rows of C ****
#ifdef MPI4
         stupid_msglen = 3*nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_complex,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,3*nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_COMPLEX,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
              call daxpy_omp(2*nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1,1),mc(taskid_i+1))
              call daxpy_omp(2*nc(taskid_j+1),1.0d0,work1(shift+shft2),
     >                                    iwrk,
     >                                    C(i+1,1,2),mc(taskid_i+1))
              call daxpy_omp(2*nc(taskid_j+1),1.0d0,work1(shift+shft3),
     >                                    iwrk,
     >                                    C(i+1,1,3),mc(taskid_i+1))
              shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do

      return
      end



*     ***********************************
*     *                                 *
*     *   CMatrix_combo_dgemm2c_omp     *
*     *                                 *
*     ***********************************
      subroutine CMatrix_combo_dgemm2c_omp(m,n,k,nblock,
     >                  A,B,lda,ma,na,
     >                  C,ldc,ldc2,mc,nc,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  tid,nthr,
     >                  work1,work2,
     >                  privatetmp)
      implicit none
      integer m,n,k,nblock

      integer lda,ma(*),na(*)
      complex*16 A(lda,*)
      complex*16 B(lda,*)

      integer ldc,ldc2,mc(*),nc(*)
      complex*16  C(ldc,ldc2,3)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      integer tid,nthr

      complex*16 work1(*),work2(*)
      complex*16 privatetmp(*)


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj,ik
      integer kk,iwrk,icur,jcur,ierr,shift,shft,shft2,shft3
      integer nthr3,nthr3l,ftid1,ftid2,ftid3,npack2,n2

      if (nthr.eq.2) then
         nthr3  = 2
         nthr3l = 2
         ftid1 = 0
         ftid2 = 0
         ftid3 = 0
      else
         nthr3  = max(1,(nthr - mod(nthr,3))/3)
         nthr3l = max(1,nthr-2*nthr3)
         ftid1 = 0
         ftid2 = min(nthr-1,nthr3)
         ftid3 = min(nthr-1,2*nthr3)
      end if
      n2     = na(taskid_j+1)
      npack2 = ma(taskid_i+1)

      call Parallel_shared_vector_zero(.false.,
     >                           2*3*nc(taskid_j+1)*mc(taskid_i+1),
     >                           C)

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)

         shft  = iwrk*ma(taskid_i+1)
         shft2 = iwrk*nc(taskid_j+1)
         shft3 = shft2+shft2

*        **** set iwrk*nc(taskid_j+1) submatrix ==0 ****
         call Parallel_shared_vector_zero(.false.,
     >                                    2*3*nc(taskid_j+1)*iwrk,
     >                                    work2)

*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then


!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))

               call zlacpy("G", ma(taskid_i+1),iwrk,
     >                   B(1,jj+1), lda,
     >                   work1(1+shft),  ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = 2*iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_complex,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,2*iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_COMPLEX,jcur,comm_j,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


            if ((iwrk.gt.0)          .and.
     >          (na(taskid_j+1).gt.0)) then

               if ((tid.ge.ftid1).and.(tid.lt.(ftid1+nthr3))) then
                  call zgemm2c_omp_group(1,tid-ftid1,nthr3,
     >                                   iwrk,n2,npack2,
     >                                   work1,A,work2,
     >                                   privatetmp)
               end if

               if ((tid.ge.ftid2).and.(tid.lt.(ftid2+nthr3))) then
                  call zgemm2c_omp_group(2,tid-ftid2,nthr3,
     >                                   iwrk,n2,npack2,
     >                                   work1,B,work2(1+shft2),
     >                                   privatetmp)
               end if

               if ((tid.ge.ftid3).and.(tid.lt.(ftid3+nthr3l))) then
                  call zgemm2c_omp_group(3,tid-ftid3,nthr3l,
     >                                   iwrk,n2,npack2,
     >                                   work1(1+shft),B,work2(1+shft3),
     >                                   privatetmp)
               end if
            end if

         end if

!$OMP BARRIER
!$OMP MASTER
c        **** summ to node that holds current rows of C ****
#ifdef MPI4
         stupid_msglen = 3*nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_complex,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,3*nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_COMPLEX,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call daxpy_omp(2*nc(taskid_j+1),1.0d0,
     >                           work1(shift),iwrk,
     >                           C(i+1,1,1),mc(taskid_i+1))
               call daxpy_omp(2*nc(taskid_j+1),1.0d0,
     >                        work1(shift+shft2),iwrk,
     >                        C(i+1,1,2),mc(taskid_i+1))
               call daxpy_omp(2*nc(taskid_j+1),1.0d0,
     >                        work1(shift+shft3),iwrk,
     >                       C(i+1,1,3),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do

      return
      end

