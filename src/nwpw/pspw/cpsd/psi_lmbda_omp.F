*     ********************************************
*     *                                          *
*     *                psi_lmbda_omp             *
*     *                                          *
*     ********************************************
!#define DUMP

      subroutine psi_lmbda_omp(ispin,ne,nemaxq,nida,nidb,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr,thrlmbda,taskid,adiff,notgram)
      USE omp_lib
      implicit none
      integer ispin,ne(2),nemaxq,nida,nidb
      complex*16 psi1(nida+nidb,nemaxq)
      complex*16 psi2(nida+nidb,nemaxq)
      real*8     dte
      real*8     lmbda(*)
      real*8     thrlmbda(*)
      real*8     tmp(*)
      integer    ierr

      integer MASTER,tid,nthr
      parameter (MASTER=0)

*     **** local variables ****
      logical notgram
      integer taskid
      integer ms,mb,mb1,mb2,it
      integer nn,i
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      real*8  adiff

*     **** matrix_blocking common block ****
      integer mblock(2),nblock(2),algorithm(2)
      common /matrix_blocking/ mblock,nblock,algorithm

*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd
      real*8  convg,convg2
      parameter (itrlmd=120, convg=1.0d-15,convg2=1.0e-12)

      integer shift,shift2,ishift2
      integer n,npack1

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads,Parallel_totalthreads
      external Parallel_threadid,Parallel_nthreads,Parallel_totalthreads

      integer np

      real*8   Dneall_m_dmax
      external Dneall_m_dmax
      double precision   eDneall_m_dmax
      external eDneall_m_dmax

      integer wtid,wnthr

      INTEGER(kind=omp_nest_lock_kind) reduce_lock1
      INTEGER(kind=omp_nest_lock_kind) reduce_lock2
      INTEGER(kind=omp_nest_lock_kind) reduce_lock3
      common / reduce_ffm / reduce_lock1,reduce_lock2,reduce_lock3


#ifdef DUMP
      integer out_unit
      character(len=20) :: filename
      logical exist
      integer j
#endif


      call nwpw_timing_start(3)

      call Parallel_np(np)
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::




      npack1 = nida+nidb

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

#ifdef DUMP
      out_unit = tid

   10 format(2E14.5)
          if (tid < 10) then
          write (filename, "(A8,I1)") "results.", tid
          else
          write (filename, "(A8,I2)") "results.", tid
          end if
          inquire(file=trim(filename), exist=exist)
          if (exist) then
            open(unit=out_unit,file=trim(filename),
     >      status="old", position="append", action="write")
          else
            open(unit=out_unit,file=trim(filename),
     >      action="write",status="new")
          end if
#endif



      DO 640 ms=1,ispin
        notgram=.true.
        IF(ne(ms).le.0) GO TO 640

        call Dneall_m_size(ms,nn)
#if 0
        s11  = 0*nn + 1
        s21  = 1*nn + 1
        s22  = 2*nn + 1
        s12  = 3*nn + 1
        sa0  = 4*nn + 1
        sa1  = 5*nn + 1
        st1  = 6*nn + 1
        st2  = 7*nn + 1
#else
        s11  = 0*nn + 1
        s21  = 1*nn + 1
        s22  = 2*nn + 1
        s12  = 3*nn + 1
        st2  = 4*nn + 1
        sa1  = 5*nn + 1
        sa0  = 6*nn + 1
        st1  = 7*nn + 1
#endif

c      if (tid.eq.0) then
c          write(*,*) 'Reduce FFMs'
c      end if

      if (ms.eq.0) then
         mb1 = 1
         mb2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         mb1 = ms
         mb2 = ms
         ishift2 = 0
      end if



      do mb=mb1,mb2
         shift  = 1 + (mb-1)*ne(1)*npack1
         shift2 = 1 + (mb-1)*ishift2
         n     = ne(mb)
         if (n.le.0) go to 30


!zero out the three matrices in one operation
!        call eDneall_ffm_zeroMatrix(tmp(s11),n,3*n)
        call eDneall_ffm_zeroMatrix(tmp(s11),n,n)
        call eDneall_ffm_zeroMatrix(tmp(s21),n,n)
        call eDneall_ffm_zeroMatrix(tmp(s22),n,n)

      call nwpw_timing_start_thr(15)    
      call nwpw_timing_start_thr(2)
        wnthr = max(1,(nthr - mod(nthr,3))/3)


        call eDneall_ffm_sym_Multiply_reduce_concur(psi2,psi2,
     >                                nida,nidb,ne,
     >                                tmp(s22),thrlmbda,
     >                       shift,shift2,n,0,wnthr,mb,reduce_lock1)

        call eDneall_ffm_sym_Multiply_reduce_concur(psi2,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s21),thrlmbda,
     >          shift,shift2,n,min(nthr-1,wnthr),wnthr,mb,reduce_lock2)

        call eDneall_ffm_sym_Multiply_reduce_concur(psi1,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s11),thrlmbda,shift,shift2,n,
     >         min(nthr-1,2*wnthr),max(1,nthr-2*wnthr),mb,reduce_lock3)

      call nwpw_timing_end_thr(2)
      call nwpw_timing_end_thr(15)

      !do one mpi reduce for all three ffms
      if (np.gt.1) then
!$OMP BARRIER
      call nwpw_timing_start_thr(15)    
!$OMP MASTER
        !do only one big reduce and then shift data accordingly
        !call Parallel_Vector_SumAll_master(3*n*n,tmp(s11))
        call Parallel_Vector_SumAll_master_tmp(3*nn,tmp(s11),tmp(st2))
        if(shift2 .ne. 0) then
          tmp(st2:st2+nn-shift2)=tmp(st2+shift2:st2+nn)
          tmp(sa1:sa1+nn-shift2)=tmp(sa1+shift2:sa1+nn)
          tmp(sa0:sa0+nn-shift2)=tmp(sa0+shift2:sa0+nn)
        end if
!$OMP END MASTER
      call nwpw_timing_end_thr(15)    
      end if

!$OMP BARRIER
  30     continue
        end do
#ifdef DUMP
!$OMP BARRIER
          write(out_unit,*) 's11:'
          write(out_unit,10) (tmp(s11+i),i=0,nn-1)
          write(out_unit,*) 's21:'
          write(out_unit,10) (tmp(s21+i),i=0,nn-1)
          write(out_unit,*) 's22:'
          write(out_unit,10) (tmp(s22+i),i=0,nn-1)
!$OMP BARRIER
#endif


#if 1
c        call Dne_ffm_combo_sym_Multiply(ms,psi1,psi2,npack1,tmp(s11))
*       ***** scale the overlap matrices ****
        call Dneall_m_scale_s22(ms,dte,tmp(s22))
        call Dneall_m_scale_s21(ms,dte,tmp(s21))
        !call dcopy(nn,tmp(s21),1,tmp(s12),1)
        call Parallel_shared_vector_copy(.true.,nn,tmp(s21),tmp(s12))
        call Dneall_m_scale_s11(ms,dte,tmp(s11))

        !call dcopy(nn,tmp(s22),1,tmp(sa0),1)
        call Parallel_shared_vector_copy(.true.,nn,tmp(s22),tmp(sa0))

        do it=1,itrlmd
          !call dcopy(nn,tmp(s22),1,tmp(sa1),1)
          call Parallel_shared_vector_copy(.true.,nn,tmp(s22),tmp(sa1))

          call Dneall_mmm_Multiply(ms,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >                              tmp(sa1),1.0d0)
          !call dcopy(nn,tmp(sa1),1,tmp(st1),1)
          call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(st1))
          call daxpy_omp(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)

          adiff = Dneall_m_dmax(ms,tmp(st1))
          if(adiff.lt.convg) GO TO 630
          if (adiff.gt.1.0d10) go to 620
          !call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
          call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(sa0))
        end do
        if (adiff.lt.convg2) GO TO 630

  620   continue
        ierr=10 
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          WRITE(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          WRITE(6,*) '        +Try using a smaller time step'
          WRITE(6,*) '        +Gram-Schmidt being performed, spin:',ms
        end if
c        call Dneall_f_ortho(ms,psi2,npack1)
        call Dneall_f_GramSchmidt(ms,psi2,npack1)
        notgram = .false.

  630   continue

*       :::::::::::::::::  correction due to the constraint  :::::::::::::::::
        if (notgram)
     >     call Dneall_fmf_Multiply(ms,
     >                              psi1,npack1,
     >                              tmp(sa1), dte,
     >                              psi2,1.0d0)

        call Dneall_mm_Expand(ms,tmp(sa1),lmbda)

#else




#if 0

!&OMP SECTIONS
!&OMP SECTION
        call eDneall_m_scale_s22(ms,ispin,ne,dte,tmp(s22))
!&OMP END SECTION
!&OMP SECTION
        call eDneall_m_scale_s21(ms,ispin,ne,dte,tmp(s21))
!&OMP END SECTION
!&OMP SECTION
        call eDneall_m_scale_s11(ms,ispin,ne,dte,tmp(s11))
!&OMP END SECTION
!&OMP END SECTIONS

!&OMP SECTIONS
!&OMP SECTION
        call dcopy(nn,tmp(s22),1,tmp(sa0),1)
!&OMP END SECTION
!&OMP SECTION
        call dcopy(nn,tmp(s21),1,tmp(s12),1)
!&OMP END SECTION
!&OMP END SECTIONS

#else

*       ***** scale the overlap matrices ****
        call Dneall_m_scale_s22(ms,dte,tmp(s22))
        call Dneall_m_scale_s21(ms,dte,tmp(s21))
        !call dcopy(nn,tmp(s21),1,tmp(s12),1)
        call Parallel_shared_vector_copy(.true.,nn,tmp(s21),tmp(s12))
        call Dneall_m_scale_s11(ms,dte,tmp(s11))

        !call dcopy(nn,tmp(s22),1,tmp(sa0),1)
        call Parallel_shared_vector_copy(.true.,nn,tmp(s22),tmp(sa0))
#endif

#if 0
        do it=1,itrlmd

#if 1
!$OMP SINGLE
          ! sa1 = s22
          call dcopy(nn,tmp(s22),1,tmp(sa1),1)
!$OMP END SINGLE
#else
          call Parallel_shared_vector_copy(.true.,nn,tmp(s22),tmp(sa1))
#endif

!!!!!      call nwpw_timing_start(49)    
      wnthr = 1!(nthr - mod(nthr,3))/3

      do mb=mb1,mb2
          n     = ne(mb)
          if (n.gt.0) then
          shift2 = 1 + (mb-1)*ishift2


          if (algorithm(mb).gt.1) then
                  call eDneall_ffm_zeroMatrix(thrlmbda(0*nn+1),n,n)
                  call eDneall_ffm_zeroMatrix(thrlmbda(1*nn+1),n,n)
                  call eDneall_ffm_zeroMatrix(thrlmbda(2*nn+1),n,n)
                  call eDneall_ffm_zeroMatrix(thrlmbda(3*nn+1),n,n)
!$OMP BARRIER
          end if

          call eDneall_mmm_Multiply_concur(ms,ispin,ne,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0,thrlmbda((0)*nn+ 1),
     >                                              mb,shift2,n,0,wnthr)
          call eDneall_mmm_Multiply_concur(ms,ispin,ne,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                   thrlmbda((3)*nn+ 1),0.0d0,thrlmbda((1)*nn+ 1),
     >                             mb,shift2,n,min(nthr-1,wnthr),wnthr)
          call eDneall_mmm_Multiply_concur(ms,ispin,ne,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0,thrlmbda((2)*nn+ 1),
     >                           mb,shift2,n,min(2*wnthr,nthr-1),wnthr)
!$OMP BARRIER
        end if
      end do

          !do the mpi reduce here
          call eDneall_mmm_Multiply_concur_red(ms,ispin,ne,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0,thrlmbda((0)*nn+ 1),
     >                                              mb,shift2,n,0,wnthr)
          call eDneall_mmm_Multiply_concur_red(ms,ispin,ne,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                   thrlmbda((3)*nn+ 1),0.0d0,thrlmbda((1)*nn+ 1),
     >                              mb,shift2,n,min(nthr-1,wnthr),wnthr)
          call eDneall_mmm_Multiply_concur_red(ms,ispin,ne,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0,thrlmbda((2)*nn+ 1),
     >                            mb,shift2,n,min(nthr-1,2*wnthr),wnthr)

!$OMP SINGLE
      call daxpy(nn,(1.0d0),thrlmbda(3*nn+1),1,tmp(sa1),1)
!$OMP END SINGLE 

      wnthr = 1
      do mb=mb1,mb2
          n     = ne(mb)
          if (n.gt.0) then
          shift2 = 1 + (mb-1)*ishift2


          if (algorithm(mb).gt.1) then
!!$OMP SINGLE
!                 call dcopy(n*n,0.0d0,0,tmp(st2),1)
!!$OMP END SINGLE
                 call eDneall_ffm_zeroMatrix(tmp(st2),n,n)
          end if

          call eDneall_mmm_Multiply_concur(ms,ispin,ne,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >            tmp(sa1),1.0d0,tmp(st2),mb,shift2,n, 0,wnthr)
!$OMP BARRIER
        end if
      end do
      !mpi reduce ?
      call eDneall_mmm_Multiply_concur_red(ms,ispin,ne,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >            tmp(sa1),1.0d0,tmp(st2),mb,shift2, 0,wnthr)
!$OMP BARRIER



!!!!!      call nwpw_timing_end(49)    


!!$OMP CRITICAL
!          write(out_unit,*) 'st1:'
!          write(out_unit,10) (tmp(sa1+i),i=0,nn-1)
!          write(out_unit,*) 'sa0:'
!          write(out_unit,10) (tmp(sa0+i),i=0,nn-1)
!!$OMP END CRITICAL
!!$OMP BARRIER

!!$OMP SINGLE
!!          !st1 = sa1
!!          call dcopy(nn,tmp(sa1),1,tmp(st1),1)
!!          !  st1 -= sa0
!!          call daxpy(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)
!!          adiff = eDneall_m_dmax(ms,ispin,ne,tmp(st1))
!!$OMP END SINGLE 

          call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(st1))
          call daxpy_omp(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)
          adiff = Dneall_m_dmax(ms,tmp(st1))


          if(adiff.lt.convg) GO TO 630
          if (adiff.gt.1.0d10) go to 620

!!!$OMP SINGLE
!!          !sa0 = sa1
!!          call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
!!          !sa1 = s22
!!          call dcopy(nn,tmp(s22),1,tmp(sa1),1)
!!!$OMP END SINGLE
        call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(sa0))

        end do
        if (adiff.lt.convg2) GO TO 630


#else
        do it=1,itrlmd
          call Parallel_shared_vector_copy(.true.,nn,tmp(s22),tmp(sa1))

          call Dneall_mmm_Multiply(ms,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >                              tmp(sa1),1.0d0)
          !call dcopy(nn,tmp(sa1),1,tmp(st1),1)
          call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(st1))
          call daxpy_omp(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)

          adiff = Dneall_m_dmax(ms,tmp(st1))
          if(adiff.lt.convg) GO TO 630
          if (adiff.gt.1.0d10) go to 620
          !call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
          call Parallel_shared_vector_copy(.true.,nn,tmp(sa1),tmp(sa0))
        end do
        if (adiff.lt.convg2) GO TO 630
#endif




  620   continue
       
        ierr=10 
        if (taskid.eq.MASTER) then
!$OMP MASTER
          write(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          write(6,*) '        +Try using a smaller time step'
          write(6,*) '        +Gram-Schmidt being performed, spin:',ms
!$OMP END MASTER
        end if
        call Dneall_f_GramSchmidt(ms,psi2,npack1)
        notgram = .false.
  630   continue

*       :::::::::::::::::  correction due to the constraint  :::::::::::::::::
        if (notgram)
     >     call Dneall_fmf_Multiply(ms,
     >                              psi1,npack1,
     >                              tmp(sa1), dte,
     >                              psi2,1.0d0)

        call Dneall_mm_Expand(ms,tmp(sa1),lmbda)
#endif

  640 continue


c*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      call nwpw_timing_end(3)

#ifdef DUMP
          close (out_unit)
#endif
      return
      end



