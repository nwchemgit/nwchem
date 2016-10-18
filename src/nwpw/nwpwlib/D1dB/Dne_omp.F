
c     ****************************************
c     *                                      *
c     *         eDneall_ffm_zeroMatrix       *
c     *                                      *
c     ****************************************
      subroutine eDneall_ffm_zeroMatrix(matrix,m,n)
      implicit none
      integer m,n
      real*8 matrix(m*n)
      integer k,j

!$OMP DO schedule(static)
      do k=1,m*n
        matrix(k) = 0.0
      end do
!$OMP END DO
      end

c     ****************************************
c     *                                      *
c     *   eDneall_ffm_sym_Multiply_reduce    *
c     *                                      *
c     ****************************************

      subroutine eDneall_ffm_sym_Multiply_reduce_concur(A1,A2,
     >                         nida,nidb,ne,hml,thrhml,
     >           shift,shift2,n,ftid,nwthr,reduce_lock)
      USE omp_lib
      implicit none
      complex*16 A1(*),A2(*)
      integer nida,nidb,ispin,ne(2)
      real*8 hml(*)
      real*8 thrhml(*)

      real*8 ratio
      integer ftid,nwthr,nthrc,nthrr,tidc,tidr
      integer n,shift,shift2,bm,offsetm
      integer nk,bk,offsetk,minchunk,nchunkc,nchunkr
      integer itid,tid,nthr,chunkc,ne1,chunkr
      integer innernthr

      INTEGER(kind=omp_nest_lock_kind) reduce_lock
*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads,Parallel_maxthreads
      external Parallel_threadid,Parallel_nthreads,Parallel_maxthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      nk = 2*nidb

      if( (tid .ge. ftid) .and. (tid .lt. ftid+nwthr)) then 
        nthrr = 1
        nthrc = nwthr / nthrr
        nwthr = nthrc*nthrr

        bk = max(min(nk,1),floor(REAL(nk)/REAL(nthrc)))
        nthrc = floor(REAL(nk) / REAL(bk))
        nwthr = nthrc*nthrr
        
        if(tid .lt. ftid+nwthr) then
          tidc = mod((tid-ftid),nthrc)
          tidr = (tid-ftid) / nthrc
        

          nchunkc = nthrc
          chunkc=tidc

          !compute offsets in A and B
          bk = floor(REAL(nk)/REAL(nchunkc))
          offsetk = (chunkc)*bk
          if(chunkc==nchunkc-1) bk = nk-offsetk
          
          ne1 = ne(1)

          nchunkr = nthrr
          chunkr = tidr
                


          bm = floor(REAL(n)/REAL(nchunkr))
          offsetm = (chunkr)*bm
          if(chunkr==nchunkr-1) bm = n-offsetm

             call epack_ccm_sym_dot_reduce_concur(nida,nidb,n,
     >                               A1(shift),
     >                               A2(shift),
     >                               hml(shift2), 
     >                               thrhml((tid)*ne1*ne1*8+ shift2),
     >                      nk,bk,offsetk,bm,offsetm,reduce_lock,ftid)
        end if
      end if
      return
      end


      subroutine epack_ccm_sym_dot_reduce_concur(nida,nidb,n,A,B,matrix,
     >             thrmatrix,nk,bk,offsetk,bm,offsetm,reduce_lock,ftid)
      USE omp_lib
      implicit none
      integer    nida,nidb,n,ftid
      real*8 A(*)
      real*8 B(*)
      real*8     matrix(n,n)
      real*8     thrmatrix(n,n)

*     **** local variables ****
      integer j,k,nk
      integer np,npack,npack2
      integer tid
      integer offsetk,bk,bm,bkc,offsetm
      integer offsetThread
      integer bkk,ibkk,kk

      INTEGER(kind=omp_nest_lock_kind) reduce_lock


*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads



      tid  = Parallel_threadid()

      bkk = bk
      npack  = (nida+nidb)
      npack2 = 2*npack


#if 1

      do kk=0,bk-1,bkk
        ibkk = min(bkk,bk-kk)

        if(kk.eq.0)then
        call DGEMM('T','N',bm,n,ibkk,
     >               2.0d0,
     >               A(1+ 2*nida + offsetk +kk + offsetm*npack2),npack2,
     >               B(1+ 2*nida + offsetk +kk ),npack2,
     >               0.0d0,
     >               thrmatrix(1+offsetm,1), n)
        else
        call DGEMM('T','N',bm,n,ibkk,
     >               2.0d0,
     >               A(1+ 2*nida + offsetk +kk + offsetm*npack2),npack2,
     >               B(1+ 2*nida + offsetk +kk ),npack2,
     >               1.0d0,
     >               thrmatrix(1+offsetm,1), n)
        endif
      end do

      if(tid.eq.ftid) then
      !write(*,*) tid,"I'm doing the first thing", 2*nida
      call DGEMM('T','N',bm,n,2*nida,
     >             1.0d0,
     >             A(1 + offsetm*npack2),npack2,
     >             B,npack2,
     >             1.0d0,
     >             thrmatrix(1+offsetm,1), n)
      end if

      if(offsetm .eq. 0 .and. bm .eq. n) then
        call omp_set_nest_lock(reduce_lock)
        call daxpy(n*n,(1.0d0),thrmatrix,1,matrix,1)
        call omp_unset_nest_lock(reduce_lock)
      else
        call omp_set_nest_lock(reduce_lock)
        do k=1,n
!$OMP SIMD
          do j=1+offsetm,offsetm+bm
            matrix(j,k) = matrix(j,k) + thrmatrix(j,k) 
          end do
!$OMP END SIMD
        end do
        call omp_unset_nest_lock(reduce_lock)
      end if
#endif

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
!$OMP DO private(k)
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

!$OMP DO private(k)
        do k=1,ne(ms)
           indx = k + (k-1)*ne(ms) + shift2
           s21(indx) = (1.0d0 - s21(indx))*0.5d0

           do j=k+1,ne(ms)
              indx  = j + (k-1)*ne(ms) + shift2
              indxt = k + (j-1)*ne(ms) + shift2

              s21(indx)  = -s21(indx)*0.5d0
              s21(indxt) = s21(indx)
           end do
        end do
!$OMP END DO

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

!$OMP DO private(k)
        do k=1,ne(ms)
           indx = k + (k-1)*ne(ms) + shift2
           s11(indx) = -s11(indx)*0.5d0*dte

           do j=k+1,ne(ms)
              indx  = j + (k-1)*ne(ms) + shift2
              indxt = k + (j-1)*ne(ms) + shift2

              s11(indx)  = -s11(indx)*0.5d0*dte
              s11(indxt) = s11(indx)
           end do
        end do
!$OMP END DO

 30     continue
      end do
      return
      end






      subroutine eDneall_mmm_Multiply_concur(mb,ispin,ne,A,B,alpha,C,
     >                                  beta,tmp,ms,shift2,n,ftid,nwthr)
      implicit none
      integer mb,ispin,ne(2)
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta
      real*8 tmp(*)

*     **** local variables ****
      integer MASTER,taskid,np,tid,nthr
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
      integer mstart,mend,nstart,nend,i,j,pindx
      integer ishiftA,ishiftB,ishiftC
      integer ftid,nwthr
*     **** matrix_blocking common block ****
      integer mblock(2),nblock(2),algorithm(2)
      common /matrix_blocking/ mblock,nblock,algorithm

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  Parallel_index_1dblock
      external Parallel_index_1dblock
      external dlafill

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if( (tid .ge. ftid) .and. (tid .lt. ftid+nwthr)) then 
        tid = tid - ftid
        

         !*** completely serial ****
         if (algorithm(ms).lt.0) then
            if (tid.eq.MASTER)
     >         call DGEMM('N','N',n,n,n,
     >                alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
          else
            pindx = tid + taskid*nwthr
            i = mod(pindx,mblock(ms))
            j = (pindx-i)/mblock(ms)
            mstart = Parallel_index_1dblock(n,mblock(ms),i)
            mend   = Parallel_index_1dblock(n,mblock(ms),i+1)
            nstart = Parallel_index_1dblock(n,nblock(ms),j)
            nend   = Parallel_index_1dblock(n,nblock(ms),j+1)
            ishiftA = shift2 + mstart 
            ishiftB = shift2 + nstart*n
            ishiftC = shift2 + mstart + nstart*n

            !*** just threaded ****
            if (algorithm(ms).lt.1) then
               call DGEMM('N','N',mend-mstart,nend-nstart,n,
     >                alpha,
     >                A(ishiftA), n,
     >                B(ishiftB), n,
     >                beta,
     >                C(ishiftC), n)
            !*** threads and cpus ****
            else
               call dlacpy('G',(mend-mstart),(nend-nstart),
     >                     C(ishiftC),n,tmp(ishiftC),n)
               call DGEMM('N','N',mend-mstart,nend-nstart,n,
     >                alpha,
     >                A(ishiftA), n,
     >                B(ishiftB), n,
     >                beta,
     >                tmp(ishiftC), n)
            end if
          end if

      end if
      return
      end





      subroutine eDneall_mmm_Multiply_concur_red(mb,ispin,ne,A,B,alpha,
     >                                  C,beta,tmp,ms,shift2,ftid,nwthr)
      implicit none
      integer mb,ispin,ne(2)
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta
      real*8 tmp(*)
      integer ftid,nwthr

*     **** local variables ****
      integer MASTER,taskid,np,tid,nthr
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
      integer mstart,mend,nstart,nend,i,j,pindx
      integer ishiftA,ishiftB,ishiftC
*     **** matrix_blocking common block ****
      integer mblock(2),nblock(2),algorithm(2)
      common /matrix_blocking/ mblock,nblock,algorithm

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  Parallel_index_1dblock
      external Parallel_index_1dblock
      external dlafill

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (mb.eq.0) then
         if ((algorithm(1).lt.1).and.(algorithm(2).lt.1)) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1)+ne(2)*ne(2),C)

         else if (algorithm(1).lt.1) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1),C)
            call Parallel_Vector_SumAll2(
     >                  ne(2)*ne(2),tmp(ne(1)*ne(1)+1),C(ne(1)*ne(1)+1))

         else if (algorithm(2).lt.1) then
            call Parallel_Vector_SumAll2(
     >                   ne(1)*ne(1),tmp,C)
            call Parallel_Brdcst_values(MASTER,
     >                   ne(2)*ne(2),C(ne(1)*ne(1)+1))

         else 
            call Parallel_Vector_SumAll2(
     >                   ne(1)*ne(1)+ne(2)*ne(2),tmp,C)
         end if
      else
         if (algorithm(mb).lt.1) then
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),C)
         else
            call Parallel_Vector_SumAll2(ne(mb)*ne(mb),tmp,C)
         end if
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *        eDneall_fmf_Multiply          *
c     *                                      *
c     ****************************************
          
      subroutine eDneall_fmf_Multiply(mb,ispin,ne,Ain,npack1,
     >                                hml,alpha,
     >                                Aout,beta)
      implicit none
      integer    mb,ispin,ne(2)
      complex*16 Ain(*)
      integer    npack1
      real*8     hml(*)
      real*8     alpha
      complex*16 Aout(*)
      real*8     beta
        

*     **** local variables ****
      integer tid,nthr,mstart,mend
      integer ms,ms1,ms2,n,shift,shift2,shift3,ishift2,ishift3

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  Parallel_index_1dblock
      external Parallel_index_1dblock

      call nwpw_timing_start(16)    
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      mstart = Parallel_index_1dblock(npack1,nthr,tid)
      mend   = Parallel_index_1dblock(npack1,nthr,tid+1)

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
         n     = ne(ms)
         if (n.le.0) go to 30
         shift  = 1 + (ms-1)*ne(1)*npack1
         shift2 = 1 + (ms-1)*ishift2

        
         !if npack / nthr >> ne divide along npack
        !if () then
        !else
         call DGEMM('N','N',2*(mend-mstart),n,n,
     >             (alpha),
     >             Ain(shift+mstart), 2*npack1,
     >             hml(shift2),    n,
     >             (beta),
     >             Aout(shift+mstart),2*npack1)
         !end if
   30    continue
      end do

      call nwpw_timing_end(16)
      return
      end




c     ****************************************
c     *                                      *
c     *        eDneall_mm_Expand             *
c     *                                      *
c     ****************************************

      subroutine eDneall_mm_Expand(mb,ne,A,A0)
      implicit none
      integer mb,ne(2)
      real*8 A(*),A0(*)


*     **** local variables ****
      integer shift2,nn

      shift2 = 1
      if (mb.eq.0) then
         nn     = ne(1)*ne(1) + ne(2)*ne(2)
         shift2 = 1
      else if (mb.eq.1) then
         nn     = ne(1)*ne(1)
         shift2 = 1
      else if (mb.eq.2) then
         nn     = ne(2)*ne(2)
         shift2 = 1+ne(1)*ne(1)
      end if


      call dcopy(nn,A,1,A0(shift2),1)
      return
      end


c     ****************************************
c     *                                      *
c     *        eDneall_m_dmax                *
c     *                                      *
c     ****************************************

      double precision function eDneall_m_dmax(mb,ispin,ne,A)
      implicit none
      integer mb,ispin,ne(2)
      real*8 A(*)


*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2
      double precision adiff1, adiff2

      integer  idamax
      external idamax

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      adiff1 = 0.0d0
      adiff2 = 0.0d0
      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift2 = 1 + (ms-1)*ishift2

        adiff1 = adiff2
        adiff2 = A(shift2-1+idamax(ne(ms)*ne(ms),A(shift2),1))
        adiff2 = dabs(adiff2)
        if (adiff2.gt.adiff1) adiff1 = adiff2
 30     continue
      end do

      eDneall_m_dmax = adiff1
      return
      end
