      subroutine tce_mrcc_diis_save(transpose,iter,t1,t2,t3,t4,
     1  d_r1m,d_t1m,k_t1_offsetm,size_t1m,
     2  d_r2m,d_t2m,k_t2_offsetm,size_t2m,
     3  d_r3m,d_t3m,k_t3_offsetm,size_t3m,
     4  d_r4m,d_t4m,k_t4_offsetm,size_t4m,
     5  rtdb)
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
c     1 mrccdebug = .false.
c       if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
c     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)

       nodezero = (ga_nodeid() .eq. 0)

       if ((diis.le.0).or.lconverged) then

       return

       endif

c
c     =====================
c         Copy T files
c     =====================
c
      do iref=1,nref

      if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then ! SG

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
      if(.not.lusesub) then
        call createfile(filename,d_t1_diis,size_t1m(iref))
      else
        call pgcreatefile(filename,d_t1_diis,size_t1m(iref),mypgid)
      endif
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
      if(lusesub) then
      call tce_mrcc_gacp(d_t1m(iref),d_t1_diis,mypgid,size_t1m(iref))
      else
      call copyfile(d_t1m(iref),d_t1_diis,size_t1m(iref))
      endif
c      if((t1_pointerm(iref).eq.1).and.(t1_diis_existm(diis,iref))) then
c          call ga_add(1.0d0,t1_diism(t1_pointerm(iref),iref),
c     1 -1.0d0,t1_diism(diis,iref),r1_diism(t1_pointerm(iref),iref))
c      endif
      endif
c
c T2
c
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
      if(.not.lusesub) then
        call createfile(filename,d_t2_diis,size_t2m(iref))
      else
        call pgcreatefile(filename,d_t2_diis,size_t2m(iref),mypgid)
      endif
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
      if(lusesub) then
      call tce_mrcc_gacp(d_t2m(iref),d_t2_diis,mypgid,size_t2m(iref))
      else
      call copyfile(d_t2m(iref),d_t2_diis,size_t2m(iref))
      endif
c      if((t2_pointerm(iref).eq.1).and.(t2_diis_existm(diis,iref))) then
c          call ga_add(1.0d0,t2_diism(t2_pointerm(iref),iref),
c     1 -1.0d0,t2_diism(diis,iref),r2_diism(t2_pointerm(iref),iref))
c      endif
      endif
c
c T3
c
      if (t3) then
      t3_pointerm(iref) = t3_pointerm(iref) + 1
      if (t3_pointerm(iref) .gt. diis) t3_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t3_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.t3_diis_existm(t3_pointerm(iref),iref)) then
        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")t3_pointerm(iref)
        call tce_filename('t3.d'//namechunk//namechunk2,filename)
        call createfile(filename,d_t3_diis,size_t3m(iref))
        call gatoeaf(d_t3_diis)
        call gatodra(d_t3_diis)
        t3_diism(t3_pointerm(iref),iref) = d_t3_diis
        t3_diis_existm(t3_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t3_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t3_diis = t3_diism(t3_pointerm(iref),iref)
      endif
      call copyfile(d_t3m(iref),d_t3_diis,size_t3m(iref))
      endif

      endif ! SG

      enddo !iref

c
c     =====================
c         Copy r files
c     =====================
c

      do iref=1,nref

      if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then ! SG

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
      if(.not.lusesub) then
        call createfile(filename,d_r1_diis,size_t1m(iref))
      else
        call pgcreatefile(filename,d_r1_diis,size_t1m(iref),mypgid)
      endif
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
c      call tce_r1_divide(d_r1,k_t1_offset)
      if(iter.eq.1) then
      if(lusesub) then
      call tce_mrcc_gacp(d_r1m(iref),d_r1_diis,mypgid,size_t1m(iref))
      else
      call copyfile(d_r1m(iref),d_r1_diis,size_t1m(iref))
      endif
      endif
      endif
c
c R2
c
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
      if(.not.lusesub) then
        call createfile(filename,d_r2_diis,size_t2m(iref))
      else
        call pgcreatefile(filename,d_r2_diis,size_t2m(iref),mypgid)
      endif
        call gatoeaf(d_r2_diis)
        call gatodra(d_r2_diis)
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
c      call tce_r2_divide(d_r2,k_t2_offset)
      if(iter.eq.1) then
      if(lusesub) then
      call tce_mrcc_gacp(d_r2m(iref),d_r2_diis,mypgid,size_t2m(iref))
      else
      call copyfile(d_r2m(iref),d_r2_diis,size_t2m(iref))
      endif
      endif
      endif
c
c R3
c
      if (t3) then
      r3_pointerm(iref) = r3_pointerm(iref) + 1
      if (r3_pointerm(iref) .gt. diis) r3_pointerm(iref) = 1
      if (nodezero.and.mrccdebug) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r3_diis_existm(i,iref),i=1,diis)
        call util_flush(LuOut)
      endif
      if (.not.r3_diis_existm(r3_pointerm(iref),iref)) then
        write(namechunk,"(I3.3)")iref
        write(namechunk2,"(I3.3)")r3_pointerm(iref)
        call tce_filename("r3.d"//namechunk//namechunk2,filename)
        call createfile(filename,d_r3_diis,size_t3m(iref))
        call gatoeaf(d_r3_diis)
        call gatodra(d_r3_diis)
        r3_diism(r3_pointerm(iref),iref) = d_r3_diis
        r3_diis_existm(r3_pointerm(iref),iref) = .true.
        if (nodezero.and.mrccdebug) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r3_diism(i,iref),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_r3_diis = r3_diism(r3_pointerm(iref),iref)
      endif
c      call tce_r3_divide(d_r3,k_t3_offset)
      call copyfile(d_r3m(iref),d_r3_diis,size_t3m(iref))
      endif

      endif ! SG
      enddo !iref

 9010 format(1x,A,' file name    = ',A)
 9020 format(1x,A,' file handles = ',20i5)


       return
       end

c $Id$
