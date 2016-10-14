*     ********************************************
*     *                                          *
*     *                psi_lmbda_omp             *
*     *                                          *
*     ********************************************

      subroutine psi_lmbda_omp(ispin,ne,nemaxq,nida,nidb,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr,thrlmbda)
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
      integer mb,ms,ms1,ms2,it
      integer nn,i
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      real*8  adiff

*     **** matrix_blocking common block ****
      integer mblock(2),nblock(2),algorithm(2)
      common /matrix_blocking/ mblock,nblock,algorithm

*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd
      real*8  convg
      parameter (itrlmd=120, convg=1.0d-15)

      integer shift,shift2,ishift2
      integer n,npack1

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads,Parallel_totalthreads
      external Parallel_threadid,Parallel_nthreads,Parallel_totalthreads

      integer np

      real*8   eDneall_m_dmax
      external eDneall_m_dmax

      integer wtid,wnthr

      INTEGER(kind=omp_nest_lock_kind) reduce_lock1
      INTEGER(kind=omp_nest_lock_kind) reduce_lock2
      INTEGER(kind=omp_nest_lock_kind) reduce_lock3
      common / reduce_ffm / reduce_lock1,reduce_lock2,reduce_lock3
      call nwpw_timing_start(3)

      call Parallel_np(np)
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::




      npack1 = nida+nidb

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      DO 640 mb=1,ispin
        notgram=.true.
        IF(ne(mb).le.0) GO TO 640

        nn = ne(mb)*ne(mb)
        s11  = 0*nn + 1
        s21  = 1*nn + 1
        s22  = 2*nn + 1
        s12  = 3*nn + 1
        st2  = 4*nn + 1
        sa1  = 5*nn + 1
        sa0  = 6*nn + 1
        st1  = 7*nn + 1

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
         shift  = 1 + (ms-1)*ne(1)*npack1
         shift2 = 1 + (ms-1)*ishift2
         n     = ne(ms)
         if (n.le.0) go to 30


!zero out the three matrices in one operation
        call eDneall_ffm_zeroMatrix(tmp(s11),n,3*n)


      call nwpw_timing_start_thr(15)    
      call nwpw_timing_start_thr(2)
        wnthr = max(1,(nthr - mod(nthr,3))/3)
        call eDneall_ffm_sym_Multiply_reduce_concur(psi2,psi2,
     >                                nida,nidb,ne,
     >                                tmp(s22),thrlmbda,
     >                       shift,shift2,n,0,wnthr,reduce_lock1)
        call eDneall_ffm_sym_Multiply_reduce_concur(psi2,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s21),thrlmbda,
     >             shift,shift2,n,min(nthr-1,wnthr),wnthr,reduce_lock2)
        call eDneall_ffm_sym_Multiply_reduce_concur(psi1,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s11),thrlmbda,shift,shift2,n,
     >            min(nthr-1,2*wnthr),max(1,nthr-2*wnthr),reduce_lock3)
      call nwpw_timing_end_thr(2)
      call nwpw_timing_end_thr(15)

      if (np.gt.1) then
!$OMP BARRIER
      call nwpw_timing_start_thr(15)    
!$OMP MASTER
        !do only one big reduce and then shift data accordingly
        call Parallel_Vector_SumAll_master(3*n*n,tmp(s11))
        if(shift2 .ne. 0) then
          tmp(st2:st2+n*n-shift2)=tmp(st2+shift2:st2+n*n)
          tmp(sa1:sa1+n*n-shift2)=tmp(sa1+shift2:sa1+n*n)
          tmp(sa0:sa0+n*n-shift2)=tmp(sa0+shift2:sa0+n*n)
        end if
!$OMP END MASTER
      call nwpw_timing_end_thr(15)    
      end if

!$OMP BARRIER
  30     continue
        end do


!&OMP SECTIONS
!&OMP SECTION
        call eDneall_m_scale_s22(mb,ispin,ne,dte,tmp(s22))
!&OMP END SECTION
!&OMP SECTION
        call eDneall_m_scale_s21(mb,ispin,ne,dte,tmp(s21))
!&OMP END SECTION
!&OMP SECTION
        call eDneall_m_scale_s11(mb,ispin,ne,dte,tmp(s11))
!&OMP END SECTION

!&OMP END SECTIONS

!$OMP SINGLE
        call dcopy(nn,tmp(s22),1,tmp(sa0),1)
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
        call dcopy(nn,tmp(s21),1,tmp(s12),1)
!$OMP END SINGLE


        do it=1,itrlmd
!$OMP SINGLE
          ! sa1 = s22
          call dcopy(nn,tmp(s22),1,tmp(sa1),1)
!$OMP END SINGLE
      call nwpw_timing_start(49)    
          wnthr = 1!(nthr - mod(nthr,3))/3

      do ms=ms1,ms2
          n     = ne(ms)
          if (n.gt.0) then
          shift2 = 1 + (ms-1)*ishift2


          if (algorithm(ms).gt.1) then
!$OMP SINGLE
                 call dcopy(n*n,0.0d0,0,thrlmbda((0)*nn+ 1),1)
!$OMP END SINGLE NOWAIT

!$OMP SINGLE
                 call dcopy(n*n,0.0d0,0,thrlmbda((1)*nn+ 1),1)
!$OMP END SINGLE NOWAIT

!$OMP SINGLE
                 call dcopy(n*n,0.0d0,0,thrlmbda((2)*nn+ 1),1)
!$OMP END SINGLE NOWAIT

!$OMP BARRIER
          end if

          call eDneall_mmm_Multiply_concur(mb,ispin,ne,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0,thrlmbda((0)*nn+ 1),
     >                                              ms,shift2,n,0,wnthr)
          call eDneall_mmm_Multiply_concur(mb,ispin,ne,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                   thrlmbda((3)*nn+ 1),0.0d0,thrlmbda((1)*nn+ 1),
     >                             ms,shift2,n,min(nthr-1,wnthr),wnthr)
          call eDneall_mmm_Multiply_concur(mb,ispin,ne,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0,thrlmbda((2)*nn+ 1),
     >                           ms,shift2,n,min(2*wnthr,nthr-1),wnthr)


!$OMP BARRIER
        end if
      end do

          !do the mpi reduce here

          call eDneall_mmm_Multiply_concur_red(mb,ispin,ne,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0,thrlmbda((0)*nn+ 1),
     >                                              ms,shift2,n,0,wnthr)
          call eDneall_mmm_Multiply_concur_red(mb,ispin,ne,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                   thrlmbda((3)*nn+ 1),0.0d0,thrlmbda((1)*nn+ 1),
     >                              ms,shift2,n,min(nthr-1,wnthr),wnthr)
          call eDneall_mmm_Multiply_concur_red(mb,ispin,ne,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0,thrlmbda((2)*nn+ 1),
     >                            ms,shift2,n,min(nthr-1,2*wnthr),wnthr)

!$OMP SINGLE
      call daxpy(nn,(1.0d0),thrlmbda(3*nn+1),1,tmp(sa1),1)
!$OMP END SINGLE 

      wnthr = 1
      do ms=ms1,ms2
          n     = ne(ms)
          if (n.gt.0) then
          shift2 = 1 + (ms-1)*ishift2


          if (algorithm(ms).gt.1) then
!$OMP SINGLE
                 call dcopy(n*n,0.0d0,0,tmp(st2),1)
!$OMP END SINGLE
          end if


          
          call eDneall_mmm_Multiply_concur(mb,ispin,ne,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >            tmp(sa1),1.0d0,tmp(st2),ms,shift2,n, 0,wnthr)
!$OMP BARRIER
        end if
      end do
      !mpi reduce ?
      call eDneall_mmm_Multiply_concur_red(mb,ispin,ne,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >            tmp(sa1),1.0d0,tmp(st2),ms,shift2, 0,wnthr)

!$OMP BARRIER
      call nwpw_timing_end(49)    

!$OMP SINGLE
          !st1 = sa1
          call dcopy(nn,tmp(sa1),1,tmp(st1),1)
          !  st1 -= sa0
          call daxpy(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)
          adiff = eDneall_m_dmax(mb,ispin,ne,tmp(st1))
!$OMP END SINGLE 

          if(adiff.lt.convg) then
                GO TO 630
          end if
          if (adiff.gt.1.0d10) then
                go to 620
          end if

!$OMP SINGLE
          !sa0 = sa1
          call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
          !sa1 = s22
          call dcopy(nn,tmp(s22),1,tmp(sa1),1)
!$OMP END SINGLE

        end do

  620   continue
       
        ierr=10 
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
!$OMP MASTER
          write(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          write(6,*) '        +Try using a smaller time step'
          write(6,*) '        +Gram-Schmidt being performed, spin:',ms
!$OMP END MASTER
        end if
c        call Dneall_f_GramSchmidt(ms,psi2,npack1)
        notgram = .false.
  630   continue

*       :::::::::::::::::  correction due to the constraint  :::::::::::::::::
        if (notgram) then

          call eDneall_fmf_Multiply(mb,ispin,ne,
     >                              psi1,nida+nidb,
     >                              tmp(sa1), dte,
     >                              psi2,1.0d0)
        end if

!$OMP SINGLE
        call eDneall_mm_Expand(mb,ne,tmp(sa1),lmbda)
!$OMP END SINGLE

  640 continue

c*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      call nwpw_timing_end(3)

      return
      end



