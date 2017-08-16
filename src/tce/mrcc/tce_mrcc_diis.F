      subroutine tce_mrcc_diis(transpose,iter,t1,t2,t3,t4,
     1  d_r1,d_t1,k_t1_offset,size_t1,
     2  d_r2,d_t2,k_t2_offset,size_t2,
     3  d_r3,d_t3,k_t3_offset,size_t3,
     4  d_r4,d_t4,k_t4_offset,size_t4,
     5  rtdb,iref)
c     5  rtdb,iref,nref)
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

      logical transpose
      logical t1, t2, t3, t4
      integer d_r1
      integer d_t1
      integer k_t1_offset
      integer d_r1_diis
      integer d_t1_diis
      integer size_t1
      integer d_r2
      integer d_t2
      integer k_t2_offset
      integer d_r2_diis
      integer d_t2_diis
      integer size_t2
      integer d_r3
      integer d_t3
      integer k_t3_offset
      integer d_r3_diis
      integer d_t3_diis
      integer size_t3
      integer d_r4
      integer d_t4
      integer k_t4_offset
      integer d_r4_diis
      integer d_t4_diis
      integer size_t4
      logical nodezero
      character*255 filename
      integer i,j
      double precision a(maxdiism+1,maxdiism+1)
      double precision b(maxdiism+1)
      integer iwork(maxdiism+1)
      integer info
      double precision ddotfile
      external ddotfile

c      integer nref,iref
      integer iref
c      logical mrccdebug
      integer rtdb
      character*3 namechunk
      character*3 namechunk2
c
c     ======================
c     Initialize parameters
c     ======================
c
c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c      mrccdebug = .false.
c       if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
c     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)

       nodezero = (ga_nodeid() .eq. 0)
c
c     =====================
c     Jacobi if DIIS is off
c     =====================
c
      if (nodezero.and.mrccdebug) then
        write(LuOut,8000) 'createfile   ', cpusecs(1), cpusecs(51)
        write(LuOut,8000) 'deletefile   ', cpusecs(2), cpusecs(52)
        write(LuOut,8000) 'put_block    ', cpusecs(3), cpusecs(53)
        write(LuOut,8000) 'add_block    ', cpusecs(4), cpusecs(54)
        write(LuOut,8000) 'get_block    ', cpusecs(5), cpusecs(55)
        write(LuOut,8000) 'reconcilefile', cpusecs(6), cpusecs(56)
      endif
 8000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
c
c     ===================================
c     Trivial case: Jacobi if DIIS is off
c     ===================================
c
      if ((diis.le.0).or.lconverged) then
        if (transpose) then
          if (t1) call tce_jacobi_l1(d_r1,d_t1,k_t1_offset)
          if (t2) call tce_jacobi_l2(d_r2,d_t2,k_t2_offset)
          if (t3) call tce_jacobi_l3(d_r3,d_t3,k_t3_offset)
          if (t4) call tce_jacobi_l4(d_r4,d_t4,k_t4_offset)
        else
          if(.not.lusesub) then
          if (t1) call tce_mrcc_jacobi_t1(d_r1,d_t1,k_t1_offset,iter,
     +            iref)
          if (t2) call tce_mrcc_jacobi_t2(d_r2,d_t2,k_t2_offset,iter,
     +            iref)
          if (t3) call tce_jacobi_t3(d_r3,d_t3,k_t3_offset,iref)
          if (t4) call tce_jacobi_t4(d_r4,d_t4,k_t4_offset)
          else
       if (t1) call tce_jacobi_sub_t1(d_r1,d_t1,k_t1_offset,iter,iref)
       if (t2) call tce_jacobi_sub_t2(d_r2,d_t2,k_t2_offset,iter,iref)
          endif
        endif
        return
      endif
c
c     =====================
c         Copy T files
c     =====================
c
      if (t1) then
      t1_pointerm(iref) = t1_pointerm(iref) + 1
      if (t1_pointerm(iref) .gt. diis) t1_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t1_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.t1_diis_existm(t1_pointerm(iref),iref)) then
        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")t1_pointerm(iref)
        call tce_filename('t1.d'//namechunk//namechunk2,filename)
        if(lusesub) then
        call pgcreatefile(filename,d_t1_diis,size_t1,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
        call createfile(filename,d_t1_diis,size_t1)
        endif
c        call gatoeaf(d_t1_diis)
c        call gatodra(d_t1_diis)
        t1_diism(t1_pointerm(iref),iref) = d_t1_diis
        t1_diis_existm(t1_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t1_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t1_diis = t1_diism(t1_pointerm(iref),iref)
      endif
      if(.not.lusesub) then
      call copyfile(d_t1,d_t1_diis,size_t1)
      else
      call tce_mrcc_gacp(d_t1,d_t1_diis,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()),size_t1)
      endif
      endif

      if (t2) then
      t2_pointerm(iref) = t2_pointerm(iref) + 1
      if (t2_pointerm(iref) .gt. diis) t2_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t2_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.t2_diis_existm(t2_pointerm(iref),iref)) then
        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")t2_pointerm(iref)
        call tce_filename('t2.d'//namechunk//namechunk2,filename)
        if(lusesub) then
        call pgcreatefile(filename,d_t2_diis,size_t2,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
        call createfile(filename,d_t2_diis,size_t2)
        endif
c        call gatoeaf(d_t2_diis)
c        call gatodra(d_t2_diis)
        t2_diism(t2_pointerm(iref),iref) = d_t2_diis
        t2_diis_existm(t2_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t2_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t2_diis = t2_diism(t2_pointerm(iref),iref)
      endif
      if(.not.lusesub) then
      call copyfile(d_t2,d_t2_diis,size_t2)
      else
      call tce_mrcc_gacp(d_t2,d_t2_diis,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()),size_t2)
      endif
      endif
c
c     =====================
c         Copy r files
c     =====================
c
      if (t1) then
      r1_pointerm(iref) = r1_pointerm(iref) + 1
      if (r1_pointerm(iref) .gt. diis) r1_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*) 'diis_exist',(r1_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.r1_diis_existm(r1_pointerm(iref),iref)) then

        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")r1_pointerm(iref)
        call tce_filename("r1.d"//namechunk//namechunk2,filename)

        if(lusesub) then
        call pgcreatefile(filename,d_r1_diis,size_t1,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
        call createfile(filename,d_r1_diis,size_t1)
        endif
c        call gatoeaf(d_r1_diis)
c        call gatodra(d_r1_diis)

        r1_diism(r1_pointerm(iref),iref) = d_r1_diis
        r1_diis_existm(r1_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r1_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_r1_diis = r1_diism(r1_pointerm(iref),iref)
      endif
      if(.not.lusesub) then
      call copyfile(d_r1,d_r1_diis,size_t1)
      else
      call tce_mrcc_gacp(d_r1,d_r1_diis,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()),size_t1)
      endif
      endif

      if (t2) then
      r2_pointerm(iref) = r2_pointerm(iref) + 1
      if (r2_pointerm(iref) .gt. diis) r2_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r2_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.r2_diis_existm(r2_pointerm(iref),iref)) then
        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")r2_pointerm(iref)
        call tce_filename("r2.d"//namechunk//namechunk2,filename)
        if(lusesub) then
        call pgcreatefile(filename,d_r2_diis,size_t2,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
        call createfile(filename,d_r2_diis,size_t2)
        endif
c        call gatoeaf(d_r2_diis)
c        call gatodra(d_r2_diis)
        r2_diism(r2_pointerm(iref),iref) = d_r2_diis
        r2_diis_existm(r2_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r2_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_r2_diis = r2_diism(r2_pointerm(iref),iref)
      endif
      if(.not.lusesub) then
      call copyfile(d_r2,d_r2_diis,size_t2)
      else
      call tce_mrcc_gacp(d_r2,d_r2_diis,
     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()),size_t2)
      endif
      endif

c
c     ==============
c     DIIS or Jacobi
c     ==============
c
      IF (r2_pointerm(iref) .le. diis) THEN
        if (transpose) then
          if (t1) call tce_jacobi_l1(d_r1,d_t1,k_t1_offset)
          if (t2) call tce_jacobi_l2(d_r2,d_t2,k_t2_offset)
          if (t3) call tce_jacobi_l3(d_r3,d_t3,k_t3_offset)
          if (t4) call tce_jacobi_l4(d_r4,d_t4,k_t4_offset)
        else
          if(.not.lusesub) then
          if (t1) call tce_mrcc_jacobi_t1(d_r1,d_t1,k_t1_offset,iter,
     +            iref)
          if (t2) call tce_mrcc_jacobi_t2(d_r2,d_t2,k_t2_offset,iter,
     +            iref)
          if (t3) call tce_jacobi_t3(d_r3,d_t3,k_t3_offset,iref)
          if (t4) call tce_jacobi_t4(d_r4,d_t4,k_t4_offset)
          else
       if (t1) call tce_jacobi_sub_t1(d_r1,d_t1,k_t1_offset,iter,iref)
       if (t2) call tce_jacobi_sub_t2(d_r2,d_t2,k_t2_offset,iter,iref)
          endif
        endif
      END IF
       
c
c     Make Matrix A (at the end of DIIS cycle)
c
      IF(r1_pointerm(iref).eq.diis) THEN 

      if(nodezero.and.mrccdebug) then
      write(6,*)'MICROCYCLE DIIS UPDATE:',
     &          iter,r1_pointerm(iref)
      end if

cjb MAKE TRUE RESIDUES
c        write(6,*)'YES1',ga_nodeid(),iref

        do i=1,diis
         
         if(i.gt.1) then
          call ga_add(1.0d0,t1_diism(i,iref),-1.0d0,t1_diism(i-1,iref),
     1 r1_diism(i,iref))
         endif

         if(i.gt.1) then
          call ga_add(1.0d0,t2_diism(i,iref),-1.0d0,t2_diism(i-1,iref),
     1 r2_diism(i,iref))
         endif

        enddo

c        write(6,*)'YES2',ga_nodeid(),iref

c        if(lusesub) call ga_pgroup_sync(
c     1 int_mb(k_innodes+ga_nnodes()+ga_nodeid()))

cjb
        do i = 1, diis
          do j = i, diis
            a(i,j) = 0.0d0
            if (t1) a(i,j) = a(i,j)
     1      + ddotfile(r1_diism(i,iref),r1_diism(j,iref),size_t1)
            if (t2) a(i,j) = a(i,j)
     1      + ddotfile(r2_diism(i,iref),r2_diism(j,iref),size_t2)
            if (t3) a(i,j) = a(i,j)
     1      + ddotfile(r3_diism(i,iref),r3_diism(j,iref),size_t3)
            if (t4) a(i,j) = a(i,j)
     1      + ddotfile(r4_diism(i,iref),r4_diism(j,iref),size_t4)
            a(j,i) = a(i,j)
          enddo
        enddo
c        write(6,*)'YES3',ga_nodeid(),iref

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
        if (nodezero.and.mrccdebug) then
          call ma_print(a,maxdiism+1,maxdiism+1,'A matrix')
          call util_flush(LuOut)
        endif
c
c     Make Vector B
c
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
        if (nodezero.and.mrccdebug) then
          call ma_print(b,maxdiism+1,1,'X matrix')
          call util_flush(LuOut)
        endif
c        write(6,*)'YES4',ga_nodeid(),iref

c
c     Make New t's
c
       if(.not.lusesub) then
        if (t1) call mrcc_daxfile(diis,b,t1_diism,d_t1,size_t1,iref)
        if (t2) call mrcc_daxfile(diis,b,t2_diism,d_t2,size_t2,iref)
        if (t3) call mrcc_daxfile(diis,b,t3_diism,d_t3,size_t3,iref)
        if (t4) call mrcc_daxfile(diis,b,t4_diism,d_t4,size_t4,iref)
       else
        if (t1) call tce_mrcc_daxp(diis,b,t1_diism,d_t1,size_t1,iref)
        if (t2) call tce_mrcc_daxp(diis,b,t2_diism,d_t2,size_t2,iref)
        if (t3) call tce_mrcc_daxp(diis,b,t3_diism,d_t3,size_t3,iref)
        if (t4) call tce_mrcc_daxp(diis,b,t4_diism,d_t4,size_t4,iref)
       endif
c        write(6,*)'YES5',ga_nodeid(),iref

        if(lusesub) then
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
        else
          call ga_sync()
        endif
      END IF



 9010 format(1x,A,' file name    = ',A)
 9020 format(1x,A,' file handles = ',20i5)

       return
       end
c
c ===================================================
c
      subroutine amp_stabilization(vl,vr,nref)
      double precision vl(nref,nref),vr(nref,nref)
      double precision tol,x,xsum
      integer aux1(nref)
      integer i,j,k,idiff,jdiff,jtot
c
      tol=1.0d-9
c vr part
      do j=1,nref  ! column
       do i=1,nref
        aux1(i)=0
       enddo
       idiff=1
       do i=1,nref
        if(aux1(i).eq.0) then
         x=dabs(vr(i,j))
         if(x.lt.tol) then
          aux1(i)=idiff
          go to 11
         end if
         do k=i,nref
          if(dabs(dabs(vr(k,j))-dabs(x)).le.tol) aux1(k)=idiff
         enddo
 11      continue
         idiff=idiff+1
        end if
       enddo
       idiff=idiff-1
       do jdiff=1,idiff
        jtot=0
        xsum=0.0d0
        do i=1,nref
         if(aux1(i).eq.jdiff) then
          jtot=jtot+1
          xsum=xsum+dabs(vr(i,j))
         end if
        enddo
c        xsum=xsum/dreal(jtot)
        xsum=xsum/dble(jtot)
        do i=1,nref
         if(aux1(i).eq.jdiff) then
           if(vr(i,j).lt.0.0d0) then
            vr(i,j)=-xsum
           else
            vr(i,j)=xsum
           end if
         end if
        enddo
       enddo
      enddo
c vl part
      do j=1,nref  ! column
       do i=1,nref
        aux1(i)=0
       enddo
       idiff=1
       do i=1,nref
        if(aux1(i).eq.0) then
         x=dabs(vl(i,j))
         if(x.lt.tol) then
          aux1(i)=idiff
          go to 12
         end if
         do k=i,nref
          if(dabs(dabs(vl(k,j))-dabs(x)).le.tol) aux1(k)=idiff
         enddo
 12      continue
         idiff=idiff+1
        end if
       enddo
       idiff=idiff-1
       do jdiff=1,idiff
        jtot=0
        xsum=0.0d0
        do i=1,nref
         if(aux1(i).eq.jdiff) then
          jtot=jtot+1
          xsum=xsum+dabs(vl(i,j))
         end if
        enddo
c        xsum=xsum/dreal(jtot)
        xsum=xsum/dble(jtot)
        do i=1,nref
         if(aux1(i).eq.jdiff) then
           if(vl(i,j).lt.0.0d0) then
            vl(i,j)=-xsum
           else
            vl(i,j)=xsum
           end if
         end if
        enddo
       enddo
      enddo
c
      return
      end
c
c ===================================================
c

c $Id$
