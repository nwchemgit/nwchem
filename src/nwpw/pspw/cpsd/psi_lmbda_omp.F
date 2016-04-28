*     ********************************************
*     *                                          *
*     *                psi_lmbda_omp             *
*     *                                          *
*     ********************************************

      subroutine psi_lmbda_omp(ispin,ne,nemaxq,nida,nidb,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr,thrlmbda)
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


      call nwpw_timing_start(3)

      call Parallel_np(np)
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::

      !call omp_debug_ids()

      npack1 = nida + nidb

!$OMP PARALLEL default(none) firstprivate(npack1) private(i,n,tid,nthr)
!$OMP& private(s22,s12,st2,sa1,sa0,st1,ms1,ms2,ishift2) 
!$OMP& private(shift,shift2,mb,notgram,nn,s11,s21,dummy)
!$OMP& firstprivate(ispin,nida,nidb,np) 
!$OMP& shared(ne,psi1,psi2,dte) 

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


        call eDneall_ffm_sym_Multiply_reduce(psi2,psi2,
     >                                nida,nidb,ne,
     >                                tmp(s22),tmp(sa0),thrlmbda,
     >                                 shift,shift2,n)

        call eDneall_ffm_sym_Multiply_reduce(psi2,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s21),tmp(sa1),thrlmbda,
     >                                 shift,shift2,n)

        call eDneall_ffm_sym_Multiply_reduce(psi1,psi1,
     >                                nida,nidb,ne,
     >                                tmp(s11),tmp(st1),thrlmbda,
     >                                 shift,shift2,n)

!$OMP BARRIER
      call nwpw_timing_start(15)    
      if (np.gt.1) then
        !do only one big reduce and then shift data accordingly
        call Parallel_Vector_SumAll(3*n*n,tmp(s11),tmp(st2))

        if(shift2 .ne. 0) then
!$OMP SINGLE
        tmp(st2:st2+n*n-shift2)=tmp(st2+shift2:st2+n*n)
        tmp(sa1:sa1+n*n-shift2)=tmp(sa1+shift2:sa1+n*n)
        tmp(sa0:sa0+n*n-shift2)=tmp(sa0+shift2:sa0+n*n)
!$OMP END SINGLE
        end if
      end if
      call nwpw_timing_end(15)    

  30     continue
        end do


*       ***** scale the overlap matrices ****
        call eDneall_m_scale_s22(mb,ispin,ne,dte,tmp(s22))
        call eDneall_m_scale_s21(mb,ispin,ne,dte,tmp(s21))
        call eDneall_m_scale_s11(mb,ispin,ne,dte,tmp(s11))

  640 continue
!$OMP END PARALLEL 


      DO 650 mb=1,ispin
        notgram=.true.
        IF(ne(mb).le.0) GO TO 650

        nn = ne(mb)*ne(mb)
        s11  = 0*nn + 1
        s21  = 1*nn + 1
        s22  = 2*nn + 1
        s12  = 3*nn + 1
        st2  = 4*nn + 1
        sa1  = 5*nn + 1
        sa0  = 6*nn + 1
        st1  = 7*nn + 1

        call dcopy(nn,tmp(s21),1,tmp(s12),1)
        call dcopy(nn,tmp(s22),1,tmp(sa0),1)

        do it=1,itrlmd
          call dcopy(nn,tmp(s22),1,tmp(sa1),1)
          call Dneall_mmm_Multiply(mb,ispin,ne,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0,tmp(st2))
          call Dneall_mmm_Multiply(mb,ispin,ne,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                              tmp(sa1),1.0d0,tmp(st2))
          call Dneall_mmm_Multiply(mb,ispin,ne,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0,tmp(st2))
          call Dneall_mmm_Multiply(mb,ispin,ne,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >                              tmp(sa1),1.0d0,tmp(st2))

          call dcopy(nn,tmp(sa1),1,tmp(st1),1)
          call daxpy(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)

          adiff = Dneall_m_dmax(mb,ispin,ne,tmp(st1))
          if(adiff.lt.convg) GO TO 630
          if (adiff.gt.1.0d10) go to 620

          call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
        end do

  620   continue
        ierr=10 
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          write(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          write(6,*) '        +Try using a smaller time step'
          write(6,*) '        +Gram-Schmidt being performed, spin:',ms
        end if
        call Dneall_f_GramSchmidt(ms,psi2,npack1)
        notgram = .false.

  630   continue

*       :::::::::::::::::  correction due to the constraint  :::::::::::::::::
        if (notgram)
     >     call Dneall_fmf_Multiply(mb,ispin,ne,
     >                              psi1,npack1,
     >                              tmp(sa1), dte,
     >                              psi2,1.0d0)

        call Dneall_mm_Expand(mb,ne,tmp(sa1),lmbda)
  650 continue

c*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      call nwpw_timing_end(3)

      return
      end



