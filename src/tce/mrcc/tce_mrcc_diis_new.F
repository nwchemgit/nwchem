      subroutine tce_mrcc_diis_new(transpose,iter,t1,t2,t3,t4,
     1  d_r1m,d_t1m,k_t1_offsetm,size_t1m,
     2  d_r2m,d_t2m,k_t2_offsetm,size_t2m,
     3  d_r3m,d_t3m,k_t3_offsetm,size_t3m,
     4  d_r4m,d_t4m,k_t4_offsetm,size_t4m,
     5  rtdb,iref)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
#include "rtdb.fh"

      integer iref
      logical transpose
      logical t1, t2, t3, t4
      integer d_r1m(maxref)
      integer d_t1m(maxref)
      integer k_t1_offsetm(maxref)
      integer d_r1_diis
      integer d_t1_diis
      integer size_t1m(maxref)
      integer d_r2m(maxref)
      integer d_t2m(maxref)
      integer k_t2_offsetm(maxref)
      integer d_r2_diis
      integer d_t2_diis
      integer size_t2m(maxref)
      integer d_r3m(maxref)
      integer d_t3m(maxref)
      integer k_t3_offsetm(maxref)
      integer d_r3_diis
      integer d_t3_diis
      integer size_t3m(maxref)
      integer d_r4m(maxref)
      integer d_t4m(maxref)
      integer k_t4_offsetm(maxref)
      integer d_r4_diis
      integer d_t4_diis
      integer size_t4m(maxref)
      logical nodezero
      character*255 filename
      integer i,j
      double precision a(maxdiism+1,maxdiism+1)
      double precision b(maxdiism+1)
      integer iwork(maxdiism+1)
      integer info
      double precision ddotfile
      external ddotfile
      double precision c(400)
c      integer nref
c      logical mrccdebug,lstatus
      logical lstatus
      integer rtdb,g_addot
      character*3 namechunk
      character*3 namechunk2
c
c     ======================
c     Initialize parameters
c     ======================
c
c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c     1 mrccdebug = .false.
c       if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
c     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)
c       if (.not.rtdb_get(rtdb,'mrcc:diistype', mt_int, 1, idiis))
c     1  idiis=0

       nodezero = (ga_nodeid() .eq. 0)

c
c     ===================================
c     Trivial case: Jacobi if DIIS is off
c     ===================================
c
      if ((diis.le.0).or.lconverged) then

        do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

         g_movecs(1) = g_movecsm(iref,1)
         g_movecs(2) = g_movecsm(iref,2)

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
c         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)
         k_evl_sorted = k_evl_sortedm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         mrccshift = epsilon - dbl_mb(k_heff+iref-1+(iref-1)*nref) ! Global denominator shift
         if(lconverged) mrccshift = 0.0d0 ! a posteriori corr
c         if(model.eq.'succsd')mrccshift = 0.0d0

        if (transpose) then
          if (t1) call tce_jacobi_l1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref))
          if (t2) call tce_jacobi_l2(d_r2m(iref),d_t2m(iref),
     2 k_t2_offsetm(iref))
          if (t3) call tce_jacobi_l3(d_r3m(iref),d_t3m(iref),
     3 k_t3_offsetm(iref))
          if (t4) call tce_jacobi_l4(d_r4m(iref),d_t4m(iref),
     4 k_t4_offsetm(iref))
        else
         if(.not.lusesub) then
          if (t1) call tce_mrcc_jacobi_t1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref),iter,iref)
          if (t2) call tce_mrcc_jacobi_t2(d_r2m(iref),d_t2m(iref),
     2 k_t2_offsetm(iref),iter,iref)
          if (t3) call tce_jacobi_t3(d_r3m(iref),d_t3m(iref),
     3 k_t3_offsetm(iref),iref)
          if (t4) call tce_jacobi_t4(d_r4m(iref),d_t4m(iref),
     4 k_t4_offsetm(iref))
          else
       if (t1) call tce_jacobi_sub_t1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref),iter,iref)
       if (t2) call tce_jacobi_sub_t2(d_r2m(iref),d_t2m(iref),
     1 k_t2_offsetm(iref),iter,iref)
          endif
        endif
        endif
        enddo

        return
      endif
c
c SG
c
      if(lusesub) then
      do iref=1,nref
        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then
          t1_pointerm(1) = t1_pointerm(iref)
        endif
      enddo
      endif
c
c     ==============
c     DIIS or Jacobi
c     ==============
c
      IF (t1_pointerm(1).le.diis) THEN

        do iref=1,nref

      if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

         g_movecs(1) = g_movecsm(iref,1)
         g_movecs(2) = g_movecsm(iref,2)

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
c         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)
         k_evl_sorted = k_evl_sortedm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb


       mrccshift = epsilon - dbl_mb(k_heff+iref-1+(iref-1)*nref) ! Global denominator shift
       if(lconverged) mrccshift = 0.0d0 ! a posteriori corr
c       if(model.eq.'succsd')mrccshift = 0.0d0

        if (transpose) then
          if (t1) call tce_jacobi_l1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref))
          if (t2) call tce_jacobi_l2(d_r2m(iref),d_t2m(iref),
     2 k_t2_offsetm(iref))
          if (t3) call tce_jacobi_l3(d_r3m(iref),d_t3m(iref),
     3 k_t3_offsetm(iref))
          if (t4) call tce_jacobi_l4(d_r4m(iref),d_t4m(iref),
     4 k_t4_offsetm(iref))
        else
         if(.not.lusesub) then
          if (t1) call tce_mrcc_jacobi_t1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref),iter,iref)
          if (t2) call tce_mrcc_jacobi_t2(d_r2m(iref),d_t2m(iref),
     2 k_t2_offsetm(iref),iter,iref)
          if (t3) call tce_jacobi_t3(d_r3m(iref),d_t3m(iref),
     3 k_t3_offsetm(iref),iref)
          if (t4) call tce_jacobi_t4(d_r4m(iref),d_t4m(iref),
     4 k_t4_offsetm(iref))
          else
       if (t1) call tce_jacobi_sub_t1(d_r1m(iref),d_t1m(iref),
     1 k_t1_offsetm(iref),iter,iref)
       if (t2) call tce_jacobi_sub_t2(d_r2m(iref),d_t2m(iref),
     1 k_t2_offsetm(iref),iter,iref)
          endif
        endif

        endif ! SG
        enddo

      END IF

c        write(6,*)'TEST 1x'

c
c     Make Matrix A (at the end of DIIS cycle)
c
      IF(t1_pointerm(1).eq.(diis)) THEN 
      if(nodezero.and.mrccdebug) then
      write(6,*)'MICROCYCLE DIIS UPDATE:',
     &          iter,r1_pointerm(iref)
      end if



       do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

        do i=1,diis

         if(i.gt.1) then
          call ga_add(1.0d0,t1_diism(i,iref),-1.0d0,t1_diism(i-1,iref),
     1 r1_diism(i-1,iref))
         endif

         if(i.gt.1) then
          call ga_add(1.0d0,t2_diism(i,iref),-1.0d0,t2_diism(i-1,iref),
     1 r2_diism(i-1,iref))
         endif
         
         if(t3) then
         if(i.gt.1) then
          call ga_add(1.0d0,t3_diism(i,iref),-1.0d0,t3_diism(i-1,iref),
     1 r3_diism(i-1,iref))
         endif
         endif

         if(i.eq.1) then
          call ga_add(1.0d0,d_t1m(iref),-1.0d0,t1_diism(diis,iref),
     1 r1_diism(diis,iref))
          call ga_add(1.0d0,d_t2m(iref),-1.0d0,t2_diism(diis,iref),
     1 r2_diism(diis,iref))
         if(t3) then
          call ga_add(1.0d0,d_t3m(iref),-1.0d0,t3_diism(diis,iref),
     1 r3_diism(diis,iref))
         endif
         endif

        enddo
        endif
       enddo


c ============================
c         idiis = 0
c ============================

        if(idiis.eq.0) then
c        do iref=1,nref

c        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

        do i = 1, diis
          do j = i, diis
            a(i,j) = 0.0d0
            do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

c          if(lusesub) then

c           if(nodezero)write(6,"('Weight:',F16.12,i5)")
c     1 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1),iref

c      if(nodezero .and. (iter .lt. 15)) write(*,*) "I have no dcmu"
            if (t1) a(i,j) = a(i,j)
     1      + ga_ddot(r1_diism(i,iref),r1_diism(j,iref)) !*
c     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
c     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

            if (t2) a(i,j) = a(i,j)
     1      + ga_ddot(r2_diism(i,iref),r2_diism(j,iref)) !*
c     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
c     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

            if (t3) a(i,j) = a(i,j)
     1      + ga_ddot(r3_diism(i,iref),r3_diism(j,iref)) !*
c     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
c     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

            if (t4) a(i,j) = a(i,j)
     1      + ga_ddot(r4_diism(i,iref),r4_diism(j,iref))

            a(j,i) = a(i,j)

        endif

            enddo
          enddo
        enddo

        do i = 1,diis
          a(i,diis+1) = -1.0d0
          a(diis+1,i) = -1.0d0
        enddo
        a(diis+1,diis+1) = 0.0d0
        do i=1,maxdiism+1
         do j=diis+2,maxdiism+1
          a(i,j) = 0.0d0
          a(j,i) = 0.0d0
         enddo
        enddo
c        if (nodezero) then
c          call ma_print(a,1,diis,'A matrix')
c          call util_flush(LuOut)
c        endif

        do i = 1,maxdiism+1
          b(i) = 0.0d0
        enddo
        b(diis+1) = -1.0d0
        if (nodezero.and.mrccdebug) then
          call ma_print(b,maxdiism+1,1,'B vector')
          call util_flush(LuOut)
        endif
c
c     Solve AX = B
c
        if(nodezero) then
c        write(*,*) "I am here diis1"
        call dgesv(diis+1,1,a,maxdiism+1,iwork,b,maxdiism+1,info)
        if (info .ne. 0)
     1    call errquit('tce_diis: LU decomposition failed',0,
     2    UNKNOWN_ERR)
        endif
        call ga_brdcst(65330,b,
     +     MA_sizeof(MT_DBL,1,MT_BYTE)*(maxdiism+1)*1,0)
        call ga_sync()
        if(nodezero) then
        call ma_print(b,1,diis,'X matrix')
        call util_flush(LuOut)
        endif
c
c     Make New t's
c
       do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

        if (t1) call mrcc_daxfile(diis,b,t1_diism,d_t1m(iref),
     1 size_t1m(iref),iref)
        if (t2) call mrcc_daxfile(diis,b,t2_diism,d_t2m(iref),
     2 size_t2m(iref),iref)
        if (t3) call mrcc_daxfile(diis,b,t3_diism,d_t3m(iref),
     3 size_t3m(iref),iref)
        if (t4) call mrcc_daxfile(diis,b,t4_diism,d_t4m(iref),
     4 size_t4m(iref),iref)

        endif

        enddo
        else ! idiis=0
c ===================
c     idiis = 1
c ===================

        do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

c        do i = 1, diis
c          do j = i, diis
c            a(i,j) = 0.0d0

cc            do iref=1,nref

cc        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

cc          if(lusesub) then

cc           if(nodezero)write(6,"('Weight:',F16.12,i5)")
cc     1 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1),iref

c            if (t1) a(i,j) = a(i,j)
c     1      + ga_ddot(r1_diism(i,iref),r1_diism(j,iref))
cc     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
cc     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

c            if (t2) a(i,j) = a(i,j)
c     1      + ga_ddot(r2_diism(i,iref),r2_diism(j,iref))
cc     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
cc     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

c            if (t3) a(i,j) = a(i,j)
c     1      + ga_ddot(r3_diism(i,iref),r3_diism(j,iref))
cc     2 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
cc     3 dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

c            if (t4) a(i,j) = a(i,j)
c     1      + ga_ddot(r4_diism(i,iref),r4_diism(j,iref))
c
c            a(j,i) = a(i,j)

cc        endif

c            enddo
c          enddo
cc       enddo

        do i = 1, diis
          do j = i, diis
            a(i,j) = 0.0d0
            if (t1) a(i,j) = a(i,j)
     1      + ddotfile(r1_diism(i,iref),r1_diism(j,iref),size_t1m(iref))
            if (t2) a(i,j) = a(i,j)
     1      + ddotfile(r2_diism(i,iref),r2_diism(j,iref),size_t2m(iref))
            if (t3) a(i,j) = a(i,j)
     1      + ddotfile(r3_diism(i,iref),r3_diism(j,iref),size_t3m(iref))
            if (t4) a(i,j) = a(i,j)
     1      + ddotfile(r4_diism(i,iref),r4_diism(j,iref),size_t4m(iref))
            a(j,i) = a(i,j)
          enddo
        enddo


        do i = 1,diis
          a(i,diis+1) = -1.0d0
          a(diis+1,i) = -1.0d0
        enddo
        a(diis+1,diis+1) = 0.0d0
        do i=1,maxdiism+1
         do j=diis+2,maxdiism+1
          a(i,j) = 0.0d0
          a(j,i) = 0.0d0
         enddo
        enddo

c        if (nodezero) then
c          call ma_print(a,1,diis,'A matrix')
c          call util_flush(LuOut)
c        endif

        do i = 1,maxdiism+1
          b(i) = 0.0d0
        enddo
        b(diis+1) = -1.0d0
        if (nodezero.and.mrccdebug) then
          call ma_print(b,maxdiism+1,1,'B vector')
          call util_flush(LuOut)
        endif
c
c     Solve AX = B
c
        call dgesv(diis+1,1,a,maxdiism+1,iwork,b,maxdiism+1,info)
        if (info .ne. 0)
     1    call errquit('tce_diis: LU decomposition failed',0,
     2    UNKNOWN_ERR)

c        if (nodezero) then
c          call ma_print(b,1,diis,'X matrix')
c          call util_flush(LuOut)
c        endif
c        do i=1,diis
c         if(abs(b(i)).gt.1000) then
c           write(6,"('Large X: ',F16.1,2i5)")b(i),ga_nodeid(),iref
c           if(b(i).gt.0.0d0) then
c             b(i) = 0.0d0
c           else
c             b(i) = 0.0d0
c           endif
c         endif
c        enddo
c
c     Make New t's
c
c       do iref=1,nref

c        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

        if (t1) call mrcc_daxfile(diis,b,t1_diism,
     1 d_t1m(iref),size_t1m(iref),iref)
        if (t2) call mrcc_daxfile(diis,b,t2_diism,d_t2m(iref),
     2 size_t2m(iref),iref)
        if (t3) call mrcc_daxfile(diis,b,t3_diism,d_t3m(iref),
     3 size_t3m(iref),iref)
        if (t4) call mrcc_daxfile(diis,b,t4_diism,d_t4m(iref),
     4 size_t4m(iref),iref)

        endif
        enddo

        endif ! idiis

        if(lusesub) then
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
          call ga_sync()
        endif

c       endif

c       enddo

c      if(lusesub)lstatus = ga_destroy(g_addot)

      END IF

 9010 format(1x,A,' file name    = ',A)
 9020 format(1x,A,' file handles = ',20i5)

       return
       end

c $Id$
