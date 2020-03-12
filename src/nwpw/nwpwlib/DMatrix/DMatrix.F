
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
 




*     ***********************************
*     *                                 *
*     *         DMatrix_tredq           *
*     *                                 *
*     ***********************************
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
      real*8  r2


      call DMatrix_eye(n,n,1.0d0,Q,lda,ma,na,taskid_i,taskid_j)
      do k = 1,n-2
          write(*,*) "k=",k

c         call DMatrix_houseq_V(k,n,A,V1,Q1,lda,ma,na,r2,
c     >            taskid_i,taskid_j,np_i,np_j,comm_i,comm_j,
c     >            work1,work2)

         call nwpw_timing_start(26)
         call DMatrix_houseq(k,n,A,V1,Q1,lda,ma,na,
     >            taskid_i,taskid_j,np_i,np_j,comm_i,comm_j,
     >            work1,work2)
         call nwpw_timing_end(26)
          write(*,*) "out houseq",k

         call nwpw_timing_start(27)
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
         call nwpw_timing_end(27)
          write(*,*) "out dgemms",k

      end do


      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_eye             *
*     *                                 *
*     ***********************************

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

!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER

      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_setdiag         *
*     *                                 *
*     ***********************************
      subroutine DMatrix_setdiag(m,n,
     >               alpha,
     >               Q,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer m,n
      real*8  alpha(*)

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
                  Q(ii,jj) = alpha(i)
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



*     ***********************************
*     *                                 *
*     *         DMatrix_s22             *
*     *                                 *
*     ***********************************
      subroutine DMatrix_s22(n,
     >               dte,
     >               Q,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  dte
      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

!$OMP SINGLE
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
               if (i.eq.j) then
                  Q(ii,jj) = (1.0d0-Q(ii,jj))*0.5d0/dte
               else
                  Q(ii,jj) = -Q(ii,jj)*0.5d0/dte
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
!$OMP END SINGLE NOWAIT

      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_s21             *
*     *                                 *
*     ***********************************
      subroutine DMatrix_s21(n,
     >               dte,
     >               Q,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  dte
      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

!$OMP SINGLE
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
               if (i.eq.j) then
                  Q(ii,jj) = (1.0d0-Q(ii,jj))*0.5d0
               else
                  Q(ii,jj) = -Q(ii,jj)*0.5d0
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
!$OMP END SINGLE NOWAIT

      return
      end


*     ***********************************
*     *                                 *
*     *         DMatrix_s11             *
*     *                                 *
*     ***********************************

      subroutine DMatrix_s11(n,
     >               dte,
     >               Q,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  dte
      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

!$OMP SINGLE
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  Q(ii,jj) = -Q(ii,jj)*0.5d0*dte
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
!$OMP END SINGLE NOWAIT

      return
      end

*     ***********************************
*     *                                 *
*     *     DMatrix_Kiril_BTransform    *
*     *                                 *
*     ***********************************
      subroutine DMatrix_Kiril_BTransform(n,
     >               s12,s21,
     >               ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      integer ldq,mq(*),nq(*)
      real*8  s12(ldq,*),s21(ldq,*)
      integer taskid_i,taskid_j

*     **** local variables ****
      integer i,j,ii,jj,im,jm,i1,j1,k,l,index
      integer icur,jcur
      real*8  tmp1

      do k=1,n
         do l=1,(k-1)
            index = (k-1) + l*(k-1) !****index = (i-1) + (j-1)*n****
            im = mod(index,n)
            jm = (index-im)/n
            i1 = im + 1
            j1 = jm + 1

            jj   = 1
            jcur = 0
            do j=1,j1-1
                jj = jj + 1
               if (jj.gt.nq(jcur+1)) then
                 jcur = jcur + 1
                 jj   = 1
               end if
            end do
            ii   = 1
            icur = 0
            do i=1,i1-1
               ii = ii +1
               if (ii.gt.mq(icur+1)) then
                  icur = icur + 1
                  ii   = 1
               end if
            end do
            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
               tmp1 = 0.5d0*(s12(ii,jj)+s21(ii,jj))
               s12(ii,jj) = tmp1
               s21(ii,jj) = tmp1
            end if
         end do
      end do
      return
      end




*     ***********************************
*     *                                 *
*     *      DMatrix_SASBtimesVtrans    *
*     *                                 *
*     ***********************************

      subroutine DMatrix_SASBtimesVtrans(n,
     >               SA,SB,
     >               Vt,A,B,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  SA(*),SB(*)
      integer ldq,mq(*),nq(*)
      real*8  Vt(ldq,*), A(ldq,*), B(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur


!$OMP SINGLE
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  A(ii,jj) = Vt(ii,jj)*SA(i)
                  B(ii,jj) = Vt(ii,jj)*SB(i)
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
!$OMP END SINGLE

      return
      end



*     ***********************************
*     *                                 *
*     *         DMatrix_trace           *
*     *                                 *
*     ***********************************
      double precision function DMatrix_trace(n,
     >                           Q,ldq,mq,nq,
     >                           taskid_i,taskid_j)
      implicit none
      integer n
      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur
      double precision sum

      sum = 0.0d0
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
               if (i.eq.j) then
                  sum = sum + Q(ii,jj) 
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

      call DMatrix_SumAll(sum)
      DMatrix_trace = sum
      return
      end




*     ***********************************
*     *                                 *
*     *         DMatrix_sqr_trace       *
*     *                                 *
*     ***********************************
      double precision function DMatrix_sqr_trace(n,
     >                           Q,ldq,mq,nq,
     >                           taskid_i,taskid_j)
      implicit none
      integer n
      integer ldq,mq(*),nq(*)
      real*8  Q(ldq,*)
      integer taskid_i,taskid_j

*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur
      double precision sum

      sum = 0.0d0
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  sum = sum + Q(ii,jj)**2
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

      call DMatrix_SumAll(sum)
      DMatrix_sqr_trace = sum
      return
      end






*     ***********************************
*     *                                 *
*     *         DMatrix_m_ijtoiijjp     *
*     *                                 *
*     ***********************************
      subroutine DMatrix_m_ijtoiijjp(i1,j1,ii,jj,p,ldq,mq,nq)
      implicit none
      integer i1,j1,ii,jj,p
      integer ldq,mq(*),nq(*)

*     **** local variables ****
      integer i,j,icur,jcur

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      jj   = 1
      jcur = 0
      do j=1,j1-1
          jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
         end if
      end do

      ii   = 1
      icur = 0
      do i=1,i1-1
         ii = ii +1
         if (ii.gt.mq(icur+1)) then
            icur = icur + 1
            ii   = 1
         end if
      end do
      p = Parallel2d_convert_taskid_ij(icur,jcur)
      return
      end


*     ***********************************
*     *                                 *
*     *         DMatrix_m_get_value     *
*     *                                 *
*     ***********************************
      real*8 function DMatrix_m_get_value(i1,j1,A,ldq,mq,nq)
      implicit none
      integer i1,j1
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)

*     **** local variables ****
      logical value
      integer taskid
      integer i,j
      integer ii,jj
      integer icur,jcur
      integer psend
      real*8 w

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)


      jj   = 1
      jcur = 0
      do j=1,j1-1
          jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
         end if
      end do

      ii   = 1
      icur = 0
      do i=1,i1-1
         ii = ii +1
         if (ii.gt.mq(icur+1)) then
            icur = icur + 1
            ii   = 1
         end if
      end do

      w = 0.0d0
      psend = Parallel2d_convert_taskid_ij(icur,jcur)
      if (psend.eq.taskid) then
          w = A(ii,jj)
      end if
      call Parallel_Brdcst_value(psend,w)

      Dmatrix_m_get_value = w
      return
      end


*     ***********************************
*     *                                 *
*     *       DMatrix_m_set_value       *
*     *                                 *
*     ***********************************
      subroutine DMatrix_m_set_value(w,i1,j1,A,ldq,mq,nq)
      implicit none
      real*8  w
      integer i1,j1
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)

*     **** local variables ****
      logical value
      integer taskid
      integer i,j
      integer ii,jj
      integer icur,jcur
      integer precv

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)

      jj   = 1
      jcur = 0
      do j=1,j1-1
          jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
         end if
      end do

      ii   = 1
      icur = 0
      do i=1,i1-1
         ii = ii +1
         if (ii.gt.mq(icur+1)) then
            icur = icur + 1
            ii   = 1
         end if
      end do

      precv = Parallel2d_convert_taskid_ij(icur,jcur)
      if (precv.eq.taskid) then
          A(ii,jj) = w
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *       DMatrix_m_add_value       *
*     *                                 *
*     ***********************************
      subroutine DMatrix_m_add_value(w,i1,j1,A,ldq,mq,nq)
      implicit none
      real*8  w
      integer i1,j1
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)

*     **** local variables ****
      logical value
      integer taskid
      integer i,j
      integer ii,jj
      integer icur,jcur
      integer precv

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)

      jj   = 1
      jcur = 0
      do j=1,j1-1
          jj = jj + 1
         if (jj.gt.nq(jcur+1)) then
           jcur = jcur + 1
           jj   = 1
         end if
      end do

      ii   = 1
      icur = 0
      do i=1,i1-1
         ii = ii +1
         if (ii.gt.mq(icur+1)) then
            icur = icur + 1
            ii   = 1
         end if
      end do

      precv = Parallel2d_convert_taskid_ij(icur,jcur)
      if (precv.eq.taskid) then
          A(ii,jj) = A(ii,jj) + w
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *      DMatrix_HmltimesSA         *
*     *                                 *
*     ***********************************

      subroutine DMatrix_HmltimesSA(n,SA,
     >               A,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  SA(*)
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  A(ii,jj) = A(ii,jj)*SA(j)
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




*     ***********************************
*     *                                 *
*     *      DMatrix_HmldivideSA        *
*     *                                 *
*     ***********************************

      subroutine DMatrix_HmldivideSA(n,SA,
     >               A,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  SA(*)
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  A(ii,jj) = A(ii,jj)/SA(j)
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





*     ***********************************
*     *                                 *
*     *      DMatrix_HmldivideDplusD    *
*     *                                 *
*     ***********************************

      subroutine DMatrix_HmldivideDplusD(n,D,
     >               A,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  D(*)
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

!$OMP MASTER
      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  A(ii,jj) = A(ii,jj)/(D(i)+D(j))
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
!$OMP END MASTER
!$OMP BARRIER

      return
      end





*     ***********************************
*     *                                 *
*     *      DMatrix_Hmlfweightscale    *
*     *                                 *
*     ***********************************
*
*     A(i,j) = A(i,j) * (2*fw(i)/(fw(i)+fw(j)))

      subroutine DMatrix_Hmlfweightscale(n,fw,
     >               A,ldq,mq,nq,
     >               taskid_i,taskid_j)
      implicit none
      integer n
      real*8  fw(*)
      integer ldq,mq(*),nq(*)
      real*8  A(ldq,*)
      integer taskid_i,taskid_j


*     **** local variables ****
      integer i,j,ii,jj
      integer icur,jcur

      jj   = 1
      jcur = 0
      do j=1,n
         ii   = 1
         icur = 0
         do i=1,n

            if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j).and.
     >          ((fw(i)+fw(j)).gt.1.0d-9)) then
                  A(ii,jj) = A(ii,jj)* (2.0d0*fw(i)/(fw(i)+fw(j)))
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

