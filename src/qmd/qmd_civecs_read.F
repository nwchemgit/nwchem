c
c     qmd_civecs_read
c
      subroutine qmd_civecs_read(rtdb,g_cix1,g_cix2,g_ciy1,g_ciy2,tda,
     $                       nocc,nvir,nstates,ipol)
c
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "msgids.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
      integer g_cix1(2)
      integer g_cix2(2)
      integer g_ciy1(2)
      integer g_ciy2(2)
      integer nocc(2)
      integer nvir(2)
      integer nstates
      integer ipol
      logical tda
c
      integer ionum1
      integer ionum2
      integer ioerr
      integer istate
      integer imo, amo
      integer ispin
      integer i,n
      integer lo(3)
      integer ld(2)
      double precision scr
c
      logical tda_ver
      integer ipol_ver
      integer nroots_ver
      integer nocc_ver(2)
      integer nmo_ver(2)
      integer nfc_ver(2)
      integer nfv_ver(2)
      integer nov_ver(2)
c
      integer l_scr1,k_scr1
      integer l_scr2,k_scr2
      integer iscr
c
      double precision scr1,scr2
c
      integer naoc(2)
      integer nav(2)
      integer nov(2)
c
      integer len_fn_civecs           ! length of civecs file name
      character*255 fn_civecs         ! file name for TDDFT vectors
      character*255 fn_civecs_current ! file name for TDDFT current vectors
      character*255 fn_civecs_old     ! file name for TDDFT old vectors
c
      integer inntsize, ok
c
      character*30 pname
c
c     preliminaries
      pname='qmd_civecs_read: '
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
c     CI Vectors files
      if (rtdb_cget(rtdb,'tddft:civecs',1,fn_civecs)) then
        len_fn_civecs = inp_strlen(fn_civecs)
        if (ipol.eq.1) then
          fn_civecs_current = fn_civecs(1:len_fn_civecs)//"_singlet"
        else
          fn_civecs_current = fn_civecs(1:len_fn_civecs)
        end if
        len_fn_civecs = inp_strlen(fn_civecs_current)
        fn_civecs_old = fn_civecs_current(1:len_fn_civecs)//"_previous"
      else
        call errquit(pname//'no saved tddft results?',0,CALC_ERR)
      endif
c      
      if (ga_nodeid().eq.0) then
       if (.not.util_get_io_unit(ionum1))
     $    call errquit(pname//'failed to get unit',1,DISK_ERR)
       open(unit=ionum1,file=fn_civecs_old,form='unformatted',
     $      status='old')
c
       if (.not.util_get_io_unit(ionum2))
     $    call errquit(pname//'failed to get unit',2,DISK_ERR)
       open(unit=ionum2,file=fn_civecs_current,form='unformatted',
     $      status='old')
c
c      Read headers
c
       rewind(ionum1) ! civecs file 1
       read(ionum1,err=90000) tda_ver    ! tda flag
       read(ionum1,err=90000) ipol_ver   ! open or closed shell
       read(ionum1,err=90000) nroots_ver  ! number of roots
       read(ionum1,err=90000) nocc_ver(1), nocc_ver(2) !occupied 1,2 
       read(ionum1,err=90000) nmo_ver(1), nmo_ver(2) ! mos 1,2
       read(ionum1,err=90000) nfc_ver(1), nfc_ver(2) ! frozen core 1,2
       read(ionum1,err=90000) nfv_ver(1), nfv_ver(2) ! frozen vir 1,2
       read(ionum1,err=90000) nov_ver(1), nov_ver(2) ! mo pairs 1,2
       read(ionum1,err=90000)
c
       rewind(ionum2) ! civecs file 2
       read(ionum2,err=90000) tda_ver
       read(ionum2,err=90000) ipol_ver
       read(ionum2,err=90000) nroots_ver
       read(ionum2,err=90000) nocc_ver(1), nocc_ver(2)
       read(ionum2,err=90000) nmo_ver(1), nmo_ver(2)
       read(ionum2,err=90000) nfc_ver(1), nfc_ver(2)
       read(ionum2,err=90000) nfv_ver(1), nfv_ver(2)
       read(ionum2,err=90000) nov_ver(1), nov_ver(2)
       read(ionum2,err=90000)
c
       do ispin = 1,ipol
        naoc(ispin)= nocc_ver(ispin)-nfc_ver(ispin)
        nav(ispin) = nmo_ver(ispin)-nocc_ver(ispin)-nfv_ver(ispin)
        nov(ispin) = naoc(ispin)*nav(ispin)
        if (nocc(ispin).ne.nocc_ver(ispin)) 
     &    call errquit(trim(pname)//"nocc problem",0,0)
       enddo
       if (ipol.eq.1) then 
         naoc(2) = 0
         nav(2)  = 0
         nov(2)  = 0
       endif 
c
c      Read CI vectors into GAs
c
       do istate=1,nstates
c
c       energy
        read(ionum1,err=90000) ! energy line
        read(ionum2,err=90000) ! s2 line
c
c       <S^2>
        read(ionum1,err=90000) ! energy line
        read(ionum2,err=90000) ! s2 line
c
        do ispin=1,ipol
c
         ld(1)=naoc(ispin)
         ld(2)=nav(ispin)
c
c        Allocate memory
         if (.not.ma_push_get(mt_dbl,nov(ispin),"slice 1",
     &    l_scr1, k_scr1))call errquit(trim(pname)//"slice 1",0,0)
         if (.not.ma_push_get(mt_dbl,nov(ispin),"slice 2",
     &    l_scr2, k_scr2))call errquit(trim(pname)//"slice 2",0,0)
c
c        Read in chunks from both civec files
         call sread(ionum1,dbl_mb(k_scr1),nov(ispin)) !read file 1
         call sread(ionum2,dbl_mb(k_scr2),nov(ispin)) !read file 2
c
         iscr = 0
         do imo=1,naoc(ispin)
          do amo=1,nav(ispin)
           iscr = iscr + 1
c
           lo(1)=imo
           lo(2)=amo
           lo(3)=istate
c
           scr1 = dbl_mb(k_scr1+iscr-1)
           call nga_put(g_cix1(ispin),lo,lo,scr1,ld)
c
           scr2 = dbl_mb(k_scr2+iscr-1)
           call nga_put(g_cix2(ispin),lo,lo,scr2,ld)
c
          end do !amo
         end do !imo
c
         if (.not.tda) then
c
c          Read in chunks from both civec files
           call sread(ionum1,dbl_mb(k_scr1),nov(ispin))! read file 1
           call sread(ionum2,dbl_mb(k_scr2),nov(ispin))! read file 2
c
           iscr = 0
           do imo=1,nocc(ispin)
             do amo=1,nvir(ispin)
              iscr = iscr + 1
c
              lo(1)=imo
              lo(2)=amo
              lo(3)=istate
c
              scr1 = dbl_mb(k_scr1+iscr-1)
              call nga_put(g_ciy1(ispin),lo,lo,scr1,ld)
c
              scr2 = dbl_mb(k_scr2+iscr-1)
              call nga_put(g_ciy2(ispin),lo,lo,scr2,ld)
c
             enddo !amo
           enddo! imo
c
         endif ! .not.tda
c
c        Deallocate memory
         if (.not.ma_pop_stack(l_scr2))
     &    call errquit(trim(pname)//"failed to pop stack: slice 2",0,0)
         if (.not.ma_pop_stack(l_scr1))
     &    call errquit(trim(pname)//"failed to pop stack: slice 1",0,0)
c
        enddo !ispin
c
       enddo !istate
c
       close(ionum1)
       close(ionum2)
       ok=1
c
      end if !ga_nodeid
c
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      return
c
90000 continue
      call errquit(pname//'error reading ci vectors',0,DISK_ERR)
      return
c
      return
      end
