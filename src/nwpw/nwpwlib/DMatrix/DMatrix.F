
*
* $Id: DMatrix.F,v 1.2 2006-09-20 19:18:48 bylaska Exp $
*

*     ***********************************************************
*     *								*
*     *   		   DMatrix library			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 5/19/06					*
*     *								*
*     ***********************************************************
c
c
c
c

c     ****************************************
c     *                                      *
c     *        DMatrix_block_mapping         *
c     *                                      *
c     ****************************************

c     Defines a SUMMA block mapping of an n1xn2 matrix.
c
      subroutine DMatrix_block_mapping(n1,n2,
     >                                 np_i,np_j,
     >                                 block_i,block_j)
      implicit none
      integer n1,n2
      integer np_i,np_j
      integer block_i(0:np_i-1),block_j(0:np_j-1)

*     **** local variables ****
      integer i,j,k
      
*     *** determine block mapping ****
      call icopy(np_i,0,1,block_i,1)
      call icopy(np_j,0,1,block_j,1)

      i = 0
      do k=1,n1
        block_i(i) = block_i(i) + 1
        i = mod(i+1,np_i)
      end do
      j = 0
      do k=1,n2
        block_j(j) = block_i(j) + 1
        j = mod(j+1,np_j)
      end do

      return
      end
 


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
         call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1, mdtob(iwrk*ma(taskid_i+1)),
     >                         jcur)

c        **** broadcast work2  within my column ***
         call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
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


*        **** pack current iwrk columns of A into work1 ***
         if (taskid_j.eq.jcur) then
            call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
         end if

c        **** broadcast work1  within my row ***
         call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1,mdtob(iwrk*ma(taskid_i+1)),
     >                         jcur)

         if ((iwrk.gt.0)          .and.
     >       (nb(taskid_j+1).gt.0).and.
     >       (ma(taskid_i+1).gt.0))
     >     call dgemm('T','N',iwrk,nb(taskid_j+1),ma(taskid_i+1),
     >                alpha,
     >                work1, ma(taskid_i+1),
     >                B, ldb,
     >                0.0d0,
     >                work2, iwrk)

c        **** summ to node that holds current rows of C ****
         call GA_PGROUP_DGOP(comm_i, 9+MSGDBL,
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

         call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
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


         call GA_PGROUP_DGOP(comm_j, 9+MSGDBL,
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


      subroutine DMatrix_tredq(n,
     >                 A,Q,lda,ma,na,
     >                 taskid_i,taskid_j,
     >                 np_i,np_j,
     >                 comm_i, comm_j,
     >                 work1,work2,
     >                 Q1,V1)
      implicit none
      integer n


      integer lda,ma(*),na(*)
      real*8  A(lda,*),Q(lda,*)

      integer taskid_i,taskid_j
      integer np_i,np_j
      integer comm_i,comm_j

      real*8  work1(*),work2(*)
      real*8  Q1(*),V1(*)

*     **** local variables ****
      integer k


      call DMatrix_eye(n,n,1.0d0,Q,lda,ma,na,taskid_i,taskid_j)
      do k = 1,n-2

         call DMatrix_houseq(k,n,A,V1,Q1,lda,ma,na,
     >            taskid_i,taskid_j,np_i,np_j,comm_i,comm_j,
     >            work1,work2)

*        **** A = Q1*A*Q1 ****
         call DMatrix_dgemm1(n,n,n,64,
     >             1.0d0,
     >             Q1,lda, ma,na,
     >             A, lda, ma,na,
     >             0.0d0,
     >             V1,lda, ma,na,
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             work1,work2)
         call DMatrix_dgemm1(n,n,n,64,
     >             1.0d0,
     >             V1,lda, ma,na,
     >             Q1,lda, ma,na,
     >             0.0d0,
     >             A, lda, ma,na,
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             work1,work2)

*        **** Q = Q*Q1 ****
        call dcopy(ma(taskid_i+1)*na(taskid_j+1),Q,1,V1,1)
         call DMatrix_dgemm1(n,n,n,64,
     >             1.0d0,
     >             V1,lda, ma,na,
     >             Q1,lda, ma,na,
     >             0.0d0,
     >             Q,lda, ma,na,
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             work1,work2)

      end do


      return
      end




      subroutine DMatrix_eye(m,n,
     >               alpha,
     >               Q,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer m,n
      real*8  alpha

      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)

      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,m

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
               if (i.eq.j) then
                  Q(ii,jj) = alpha
               else
                  Q(ii,jj) = 0.0d0
               end if
            end if

            ii = ii + 1
            if (ii.gt.mq(icur+1)) then
              icur = icur + 1
              ii   = 1
            end if
         end do

         jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
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


               call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1, mdtob(mq(taskid_i+1)),
     >                         jcur0)

               call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work2, mdtob(mq(taskid_i+1)),
     >                         jcur1)


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
         call GA_PGROUP_DGOP(comm_i,9+MSGDBL,mu,1,'+')
         mu = dsqrt(mu)


         if (mu.ne.0.0d0) then
           if (icur1.eq.taskid_i)
     >        beta = V(ii1,jj0) + dsign(mu,V(ii1,jj0))
              call GA_PGROUP_BRDCST(comm_i,9+MSGDBL,beta,mdtob(1),icur1)

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
         call GA_PGROUP_DGOP(comm_i,9+MSGDBL,v2,1,'+')


         v2 = 2.0d0/v2
      end if
      call GA_PGROUP_BRDCST(comm_j,9+MSGDBL,v2,mdtob(1),jcur0)




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

              call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1, mdtob(mq(taskid_i+1)),
     >                         jcur0)
              call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work2, mdtob(mq(taskid_i+1)),
     >                         jcur1)


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

         call GA_PGROUP_BRDCST(comm_i, 9+MSGDBL,
     >                         work1(j), mdtob(1),
     >                         icur)
         jj = jj + 1
         if (jj.gt.na(jcur+1)) then
            jcur = jcur + 1
            jj = 1
         end if
      end do
      call GA_PGROUP_DGOP(comm_j,9+MSGDBL,work1,n,'+')
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
         call GA_PGROUP_BRDCST(comm_j, 9+MSGDBL,
     >                         work1(i), mdtob(1),
     >                         jcur)
         ii = ii + 1
         if (ii.gt.ma(icur+1)) then
            icur = icur + 1
            ii = 1
         end if
      end do
      call GA_PGROUP_DGOP(comm_i,9+MSGDBL,work1,n-1,'+')
      call dcopy(n-1,work1,1,tu,1)

      return
      end


