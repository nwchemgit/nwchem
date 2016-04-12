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
       rewind(ionum1)
c      TDA flag, ipol, #roots, #occupied, #MOs, #frozen core, 
c       #frozen virtuals, #MO pairs, blank
       read(ionum1,err=90000) tda_ver
       read(ionum1,err=90000) ipol_ver
       read(ionum1,err=90000) nroots_ver
       read(ionum1,err=90000) nocc_ver(1), nocc_ver(2)
       read(ionum1,err=90000) nmo_ver(1), nmo_ver(2)
       read(ionum1,err=90000) nfc_ver(1), nfc_ver(2)
       read(ionum1,err=90000) nfv_ver(1), nfv_ver(2)
       read(ionum1,err=90000) nov_ver(1), nov_ver(2)
       read(ionum1,err=90000)
c
       rewind(ionum2)
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
c      Read CI vectors into GAs
c
       do istate=1,nstates
        lo(3)=istate
c       energy
        read(ionum1,err=90000) ! energy line
        read(ionum2,err=90000) ! s2 line
c       <S^2>
        read(ionum1,err=90000) ! energy line
        read(ionum2,err=90000) ! s2 line
        do ispin=1,ipol
         ld(1)=nocc(ispin)
         ld(2)=nvir(ispin)
         do imo=1,nocc(ispin)
          lo(1)=imo
          do amo=1,nvir(ispin)
           lo(2)=amo
           read(ionum1,err=90000) scr
           call nga_put(g_cix1(ispin),lo,lo,scr,ld)
           read(ionum2,err=90000) scr
           call nga_put(g_cix2(ispin),lo,lo,scr,ld)
          end do !amo
         end do !imo
        end do !ispin
c
        if (.not.tda) then
         do ispin=1,ipol
          ld(1)=nocc(ispin)
          ld(2)=nvir(ispin)
          do imo=1,nocc(ispin)
           lo(1)=imo
           do amo=1,nvir(ispin)
            lo(2)=amo
            read(ionum1,err=90000) scr
            call nga_put(g_ciy1(ispin),lo,lo,scr,ld)
            read(ionum2,err=90000) scr
            call nga_put(g_ciy2(ispin),lo,lo,scr,ld)
           end do !amo
          end do !imo
         end do !ispin
        end if !tda
       end do !istate
       close(ionum1)
       close(ionum2)
       ok=1
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
