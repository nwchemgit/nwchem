
c     ****************************************
c     *                                      *
c     *         eDneall_ffm_zeroMatrix       *
c     *                                      *
c     ****************************************
      subroutine eDneall_ffm_zeroMatrix(matrix,m,n)
      implicit none
      integer m,n
      real*8 matrix(m,n)
      integer k,j

!$OMP DO schedule(static)
      do k=1,n
!$OMP SIMD
      do j=1,m
        matrix(j,k) = 0.0
      end do
!$OMP END SIMD
      end do
!$OMP END DO
      end

c     ****************************************
c     *                                      *
c     *   eDneall_ffm_sym_Multiply_reduce    *
c     *                                      *
c     ****************************************
      subroutine eDneall_ffm_sym_Multiply_reduce(A1,A2,
     >                         nida,nidb,ne,hml,tmp,thrhml,
     >                        shift,shift2,n)
      implicit none
      complex*16 A1(*),A2(*)
      integer nida,nidb,ispin,ne(2)
      real*8 hml(*),tmp(*)
      real*8 thrhml(*)

      integer n,shift,shift2
      integer tid

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(15)    
      tid  = Parallel_threadid()
         call epack_ccm_sym_dot_reduce(nida,nidb,n,
     >                     A1(shift),
     >                     A2(shift),
     >                     hml(shift2), tmp,
     >                     thrhml((tid)*ne(1)*ne(1)*8+ shift2))
      call nwpw_timing_end(15)    
      return
      end


*     ***********************************
*     *                                 *
*     *      epack_ccm_sym_dot_reduce   *       
*     *                                 *
*     ***********************************
      subroutine epack_ccm_sym_dot_reduce(nida,nidb,n,A,B,matrix,tmp,
     >                                                    thrmatrix)
      implicit none
      integer    nida,nidb,n
      real*8 A(*)
      real*8 B(*)
      real*8     matrix(n,n)
      real*8     thrmatrix(n,n)
      real*8     tmp(*)

*     **** local variables ****
      integer j,k,nk
      integer np,npack,npack2
      integer tid,nthr
      integer offsetk,bk,bkc
      integer offsetThread

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(2)

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      npack  = (nida+nidb)
      npack2 = 2*npack

      !compute offsets in A and B
      nk = npack2 - 2*nida
      bk = floor(REAL(nk)/REAL(nthr))
      offsetk = (tid)*bk
      if(tid==nthr-1) bk = nk-offsetk

      call DGEMM('T','N',n,n,bk,
     >             2.0d0,
     >             A(1+ 2*nida + offsetk),npack2,
     >             B(1+ 2*nida + offsetk),npack2,
     >             0.0d0,
     >             thrmatrix, n)
      if(tid.eq.0) then
      call DGEMM('T','N',n,n,2*nida,
     >             1.0d0,
     >             A,npack2,
     >             B,npack2,
     >             1.0d0,
     >             thrmatrix, n)
      end if

      !Do the reduction
!$OMP CRITICAL
      do k=1,n
!$OMP SIMD
      do j=1,n
        matrix(j,k) = matrix(j,k) + thrmatrix(j,k) 
      end do
!$OMP END SIMD
      end do
!$OMP END CRITICAL

      call nwpw_timing_end(2)
      return
      end










c     ****************************************
c     *                                      *
c     *        eDneall_m_scale_s22           *
c     *                                      *
c     ****************************************

      subroutine eDneall_m_scale_s22(mb,ispin,ne,dte,s22)
      implicit none
      integer mb,ispin,ne(2)
      real*8 dte
      real*8 s22(*)
        

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift2 = (ms-1)*ishift2

!$OMP DO private(k) schedule(static)
        do k=1,ne(ms)
           indx = k + (k-1)*ne(ms) + shift2
           s22(indx) = (1.0d0 - s22(indx))*0.5d0/dte

           do j=k+1,ne(ms)
              indx  = j + (k-1)*ne(ms) + shift2
              indxt = k + (j-1)*ne(ms) + shift2

              s22(indx)  = -s22(indx)*0.5d0/dte
              s22(indxt) = s22(indx)
           end do
        end do
!$OMP END DO

 30     continue
      end do

      return
      end



c     ****************************************
c     *                                      *
c     *       eDneall_m_scale_s21            *
c     *                                      *
c     ****************************************

      subroutine eDneall_m_scale_s21(mb,ispin,ne,dte,s21)
      implicit none
      integer mb,ispin,ne(2)
      real*8 dte
      real*8 s21(*)

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt
      integer tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift2 = (ms-1)*ishift2

        do k=tid+1,ne(ms),nthr
           indx = k + (k-1)*ne(ms) + shift2
           s21(indx) = (1.0d0 - s21(indx))*0.5d0

           do j=k+1,ne(ms)
              indx  = j + (k-1)*ne(ms) + shift2
              indxt = k + (j-1)*ne(ms) + shift2

              s21(indx)  = -s21(indx)*0.5d0
              s21(indxt) = s21(indx)
           end do
        end do

 30     continue
      end do
      return
      end


c     ****************************************
c     *                                      *
c     *        eDneall_m_scale_s11           *
c     *                                      *
c     ****************************************

      subroutine eDneall_m_scale_s11(mb,ispin,ne,dte,s11)
      implicit none
      integer mb,ispin,ne(2)
      real*8 dte
      real*8 s11(*)

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt
      integer tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift2 = (ms-1)*ishift2

        do k=tid+1,ne(ms),nthr
           indx = k + (k-1)*ne(ms) + shift2
           s11(indx) = -s11(indx)*0.5d0*dte

           do j=k+1,ne(ms)
              indx  = j + (k-1)*ne(ms) + shift2
              indxt = k + (j-1)*ne(ms) + shift2

              s11(indx)  = -s11(indx)*0.5d0*dte
              s11(indxt) = s11(indx)
           end do
        end do

 30     continue
      end do
      return
      end


