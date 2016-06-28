
*
* $Id$
*


*     ***********************************************************
*     *								*
*     *   		   DMatrix library			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 5/19/06					*
*     *								*
*     ***********************************************************




*     ***********************************
*     *                                 *
*     *         DMatrix_dgemm1_rot2     *
*     *                                 *
*     ***********************************

      subroutine DMatrix_dgemm1_rot2(m,n,k,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
c      real*8  B(ldb,*)
      real*8  B(*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  Bcol(*),Bwork(*)
      real*8  work1(*),work2(*)


      call errquit("TCGMSG DMatrix_dgemm1_rot2 not implemented",0,0)
      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_dgemm1_rot      *
*     *                                 *
*     ***********************************

      subroutine DMatrix_dgemm1_rot(m,n,k,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
c      real*8  B(ldb,*)
      real*8  B(*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  Bcol(*),Bwork(*)
      real*8  work1(*),work2(*)


      call errquit("TCGMSG DMatrix_dgemm1_rot not implemented",0,0)
      return
      end


*     ***********************************
*     *                                 *
*     *         DMatrix_dgemm1          *
*     *                                 *
*     ***********************************
c
      subroutine DMatrix_dgemm1(m,n,k,nblock,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
      real*8  B(ldb,*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

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


*        **** pack current iwrk columns of A into work1 ***
         if (taskid_j.eq.jcur) then
            call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
         end if

*        **** pack current iwrk rows of B into work2 ***
         if (taskid_i.eq.icur) then
            call dlacpy("G", iwrk,nb(taskid_j+1),
     >                   B(ii+1,1), ldb,
     >                   work2,  iwrk)
         end if

c        **** broadcast work1  within my row ***
         if (np_j.gt.1)
     >    call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1, mdtob(iwrk*ma(taskid_i+1)),
     >                         jcur)

c        **** broadcast work2  within my column ***
         if (np_i.gt.1)
     >    call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
     >                         work2, mdtob(iwrk*nb(taskid_j+1)),
     >                         icur)

         if ((iwrk.gt.0)          .and.
     >       (mc(taskid_i+1).gt.0).and.
     >       (nc(taskid_j+1).gt.0))
     >     call dgemm('N','N',mc(taskid_i+1),nc(taskid_j+1),iwrk,
     >                alpha,
     >                work1, ma(taskid_i+1),
     >                work2, iwrk,
     >                1.0d0,
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




      subroutine DMatrix_dgemm2(m,n,k,nblock,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
      real*8  B(ldb,*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

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
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)


*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
         if (np_j.gt.1)
     >    call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1,mdtob(iwrk*ma(taskid_i+1)),
     >                         jcur)

c            if ((iwrk.gt.0)          .and.
c     >          (nb(taskid_j+1).gt.0).and.
c     >          (ma(taskid_i+1).gt.0))
            if ((iwrk.gt.0)          .and.
     >          (nb(taskid_j+1).gt.0))
     >        call dgemm('T','N',iwrk,nb(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, ldb,
     >                   0.0d0,
     >                   work2, iwrk)

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            call dcopy(nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
         end if


c        **** summ to node that holds current rows of C ****
        if (np_i.gt.1)
     >   call GA_PGROUP_DGOP(comm_i, 9+MSGDBL,
     >                       work2,nc(taskid_j+1)*iwrk,'+')


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call daxpy(nc(taskid_j+1),1.0d0,work2(shift),iwrk,
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



      subroutine DMatrix_dgamm2(m,n,k,nblock,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
      real*8  B(ldb,*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

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
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)


*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
         if (np_j.gt.1)
     >    call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1,mdtob(iwrk*ma(taskid_i+1)),
     >                         jcur)

c            if ((iwrk.gt.0)          .and.
c     >          (nb(taskid_j+1).gt.0).and.
c     >          (ma(taskid_i+1).gt.0))
            if ((iwrk.gt.0)          .and.
     >          (nb(taskid_j+1).gt.0))
     >        call dgamm_omp('T','N',iwrk,nb(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, ldb,
     >                   0.0d0,
     >                   work2, iwrk)

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            call dcopy(nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
         end if


c        **** summ to node that holds current rows of C ****
        if (np_i.gt.1)
     >   call GA_PGROUP_DGOP(comm_i, 9+MSGDBL,
     >                       work2,nc(taskid_j+1)*iwrk,'+')


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call daxpy(nc(taskid_j+1),1.0d0,work2(shift),iwrk,
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




      subroutine DMatrix_dgemm3(m,n,k,nblock,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)

      integer ldb,mb(*),nb(*)
      real*8  B(ldb,*)

      real*8  beta

      integer ldc,mc(*),nc(*)
      real*8  C(ldc,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      real*8  dum


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
      do while (kk.lt.n)
         iwrk = min(nblock, mb(icur+1)-ii)
         iwrk = min(iwrk,   nc(jcur+1)-jj)


         if (taskid_i.eq.icur) then
            call dlacpy("G", iwrk,nb(taskid_j+1),
     >                   B(ii+1,1), ldb,
     >                   work2,     iwrk)
         end if

         if (np_i.gt.1)
     >      call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
     >                         work2, mdtob(iwrk*nb(taskid_j+1)),
     >                         icur)

         if ((iwrk.gt.0)          .and.
     >       (na(taskid_j+1).gt.0).and.
     >       (mc(taskid_i+1).gt.0))
     >      call dgemm('N','T',mc(taskid_i+1),iwrk,na(taskid_j+1),
     >              alpha,
     >              A, lda,
     >              work2, iwrk,
     >              0.0d0,
     >              work1, mc(taskid_i+1))

         if (np_j.gt.1)
     >    call GA_PGROUP_DGOP(comm_j, 9+MSGDBL,
     >                       work1,mc(taskid_i+1)*iwrk,'+')


         if (taskid_j.eq.jcur) then
            shift = 1
            do j=jj,(jj+iwrk-1)
               call daxpy(mc(taskid_i+1),
     >                    1.0d0,
     >                    work1(shift),1,
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




      subroutine DMatrix_tqliq(n,eig,tu,
     >                  Q,ldq,mq,nq,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer n

      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*),eig(*),tu(*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      real*8  work1(*),work2(*)


#include "tcgmsg.fh"
#include "msgtypesf.h"

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
            g = eig(m)-eig(l)+tu(l)/(g+dsign(r,g))
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
     >             call dcopy(mq(taskid_i+1),Q(1,jj0),1,work1,1)
               if (jcur1.eq.taskid_j)
     >             call dcopy(mq(taskid_i+1),Q(1,jj1),1,work2,1)

               if (np_j.gt.1) then
                 call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                           work1, mdtob(mq(taskid_i+1)),
     >                           jcur0)
                 call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                           work2, mdtob(mq(taskid_i+1)),
     >                           jcur1)
               end if


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


      subroutine DMatrix_houseq(jcol,
     >                  n,
     >                  A,V,Q,lda,ma,na,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer jcol,n

      integer lda,ma(*),na(*)
      real*8  A(lda,*),V(lda,*),Q(lda,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift
      integer ii0,icur0,ii1,icur1,ii2,icur2
      integer jj0,jcur0,jj1,jcur1
      real*8  beta,mu0,mu,v20,v2

      call dcopy(ma(taskid_i+1)*na(taskid_j+1),0.0d0,0,V,1)

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
            if (icur.eq.taskid_i) mu0 = mu0 + V(ii,jj0)**2
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do
         mu = mu0
         if (np_i.gt.1) call GA_PGROUP_DGOP(comm_i,9+MSGDBL,mu,1,'+')
         mu = dsqrt(mu)


         if (mu.ne.0.0d0) then
           if (icur1.eq.taskid_i)
     >        beta = V(ii1,jj0) + dsign(mu,V(ii1,jj0))
              if (np_i.gt.1) 
     >        call GA_PGROUP_BRDCST(comm_i,9+MSGDBL,beta,mdtob(1),icur1)

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
         if (icur1.eq.taskid_i)  V(ii1,jj0) = 1.0d0
         if (icur0.eq.taskid_i)  V(ii0,jj0) = 0.0d0

         v20 = 0.0d0
         icur = icur0
         ii   = ii0
         do i=jcol,n
            if (icur.eq.taskid_i) v20 = v20 + V(ii,jj0)**2
            ii = ii + 1
            if (ii.gt.ma(icur+1)) then
               icur = icur + 1
               ii = 1
            end if
         end do
         v2 = v20
         if (np_i.gt.1) call GA_PGROUP_DGOP(comm_i,9+MSGDBL,v2,1,'+')

         v2 = 2.0d0/v2
      end if
      if (np_j.gt.1)
     >   call GA_PGROUP_BRDCST(comm_j,9+MSGDBL,v2,mdtob(1),jcur0)


      call DMatrix_eye(n,n,1.0d0,Q,lda,ma,na,taskid_i,taskid_j)
      call DMatrix_dgemm3(n,n,n,64,
     >             (-v2),
     >             V,ma(taskid_i+1), ma,na,
     >             V,ma(taskid_i+1), ma,na,
     >             1.0d0,
     >             Q,ma(taskid_i+1), ma,na,
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             work1,work2)


      return
      end


      subroutine DMatrix_eigsrtq(n,eig,
     >                  Q,ldq,mq,nq,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i, comm_j,
     >                  work1,work2)
      implicit none
      integer n

      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*),eig(*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      real*8  work1(*),work2(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

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
     >         call dcopy(mq(taskid_i+1),Q(1,jj0),1,work1,1)
            if (jcur1.eq.taskid_j)
     >         call dcopy(mq(taskid_i+1),Q(1,jj1),1,work2,1)

              if (np_j.gt.1) then
              call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1, mdtob(mq(taskid_i+1)),
     >                         jcur0)
              call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work2, mdtob(mq(taskid_i+1)),
     >                         jcur1)
              end if


            if (jcur0.eq.taskid_j)
     >         call dcopy(mq(taskid_i+1),work2,1,Q(1,jj0),1)
            if (jcur1.eq.taskid_j)
     >         call dcopy(mq(taskid_i+1),work1,1,Q(1,jj1),1)

         end if

      end do

      return
      end


      subroutine DMatrix_getdiags(n,eig,tu,
     >                  A,lda,ma,na,
     >                  taskid_i,taskid_j,
     >                  np_i,np_j,
     >                  comm_i,comm_j,
     >                  work1)
      implicit none
      integer n

      integer lda,ma(*),na(*)
      real*8  A(lda,*),eig(*),tu(*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j
      real*8  work1(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"


*     **** local variables ****
      integer i,j,ii,jj,is,ie,js,je
      integer icur,jcur

*     **************************
*     **** gather diagonals ****
*     **************************
      call dcopy(n,0.0d0,0,work1,1)
      call dcopy(n,0.0d0,0,eig,1)
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
         work1(j) = A(ii,jj)

         if (np_i.gt.1)
     >      call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
     >                           work1(j), mdtob(1),
     >                           icur)
         jj = jj + 1
         if (jj.gt.na(jcur+1)) then
            jcur = jcur + 1
            jj = 1
         end if
      end do
      if (np_j.gt.1) call GA_PGROUP_DGOP(comm_j,9+MSGDBL,work1,n,'+')
      call dcopy(n,work1,1,eig,1)


*     ******************************
*     **** gather off-diagonals ****
*     ******************************
      call dcopy(n,0.0d0,0,work1,1)
      call dcopy(n,0.0d0,0,tu,1)
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
         work1(i) = A(ii,jj)
         if (np_j.gt.1)
     >    call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1(i), mdtob(1),
     >                         jcur)
         ii = ii + 1
         if (ii.gt.ma(icur+1)) then
            icur = icur + 1
            ii = 1
         end if
      end do
      if (np_i.gt.1) call GA_PGROUP_DGOP(comm_i,9+MSGDBL,work1,n-1,'+')
      call dcopy(n-1,work1,1,tu,1)

      return
      end




*     ***********************************
*     *                                 *
*     *         DMatrix_MaxAll          *
*     *                                 *
*     ***********************************
      subroutine DMatrix_MaxAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer msglen,mpierr,np

*     **** external functions ****

      call Parallel_np(np)
      if (np.gt.1) then
        call GA_DGOP(9+MSGDBL,sum,1,'max')
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_SumAll          *
*     *                                 *
*     ***********************************
      subroutine DMatrix_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer mpierr,np

*     **** external functions ****

      call Parallel_np(np)
      if (np.gt.1) then
        call GA_DGOP(9+MSGDBL,sum,1,'+')
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *         DMatrix_mm_transpose    *
*     *                                 *
*     ***********************************

      subroutine DMatrix_mm_transpose(n,A,B,ldq,mq,nq)
      implicit none
      integer n
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)
      real*8  B(ldq,*)

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer taskid
      integer i,j
      integer ii,jj,rr,ss
      integer icur,jcur,rcur,scur
      integer psend,precv,msglen,msgtype,mpierr
      integer rcv_len,rcv_proc

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij


      call Parallel_taskid(taskid)
      msglen  = 1
      msgtype = 1


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
               !**** send ****
               if (psend.eq.taskid)
     >         call SND(9+MSGDBL,
     >                  A(ii,jj),mdtob(msglen),
     >                  precv,
     >                  1)


               !**** recv ****
               if (precv.eq.taskid)
     >         call RCV(9+MSGDBL,
     >                  B(rr,ss),mdtob(msglen),rcv_len,
     >                  psend,
     >                  rcv_proc,1)

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

      return
      end



*     ***********************************
*     *                                 *
*     *       DMatrix_combo_dgemm2      *
*     *                                 *
*     ***********************************
      subroutine DMatrix_combo_dgemm2(m,n,k,nblock,
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
      real*8  alpha

      integer lda,ma(*),na(*)
      real*8  A(lda,*)
      real*8  B(lda,*)

      real*8  beta

      integer ldc,ldc2,mc(*),nc(*)
      real*8  C(ldc,ldc2,3)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)


#include "tcgmsg.fh"
#include "msgtypesf.h"


*     **** local variables ****
      logical docalc1,docalc2
      integer i,j,ii,jj
      integer kk,iwrk,icur,jcur,ierr,shift,shft,shft2,shft3

      do kk=1,3
      do j=1,nc(taskid_j+1)
         do i=1,mc(taskid_i+1)
            C(i,j,kk) = beta*C(i,j,kk)
         end do
      end do
      end do

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

*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))

               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   B(1,jj+1), lda,
     >                   work1(1+shft),     ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
            if (np_j.gt.1)
     >        call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1,mdtob(2*iwrk*ma(taskid_i+1)),
     >                         jcur)

            if ((iwrk.gt.0)          .and.
     >          (na(taskid_j+1).gt.0)) then

              call dgemm('T','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   A, lda,
     >                   0.0d0,
     >                   work2, iwrk)
              call dgemm('T','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, lda,
     >                   0.0d0,
     >                   work2(1+shft2), iwrk)
              call dgemm('T','N',iwrk,na(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1(1+shft), ma(taskid_i+1),
     >                   B, lda,
     >                   0.0d0,
     >                   work2(1+shft3), iwrk)
            end if

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            call dcopy(3*nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
         end if

c        **** summ to node that holds current rows of C ****
        if (np_i.gt.1)
     >   call GA_PGROUP_DGOP(comm_i, 9+MSGDBL,
     >                       work2,3*nc(taskid_j+1)*iwrk,'+')


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call daxpy(nc(taskid_j+1),1.0d0,work2(shift),iwrk,
     >                                    C(i+1,1,1),mc(taskid_i+1))
               call daxpy(nc(taskid_j+1),1.0d0,work2(shift+shft2),iwrk,
     >                                    C(i+1,1,2),mc(taskid_i+1))
               call daxpy(nc(taskid_j+1),1.0d0,work2(shift+shft3),iwrk,
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

