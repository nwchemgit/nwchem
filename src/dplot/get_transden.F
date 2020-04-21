c
c        construct the transition density
         subroutine get_transden(spin,iroot,ipol,filename,basis,
     &        g_movecs, g_tdens)
c
         implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
c
         character*8 spin
         integer iroot
         integer ipol
         character*256 filename
         integer basis            ! AO basis set handle
         integer g_movecs(ipol)   ! MO vectors
         integer g_tdens(ipol)    ! Transition density matrix
c
         integer l_trials, k_trials
c
         integer i,m,n
         logical tda
         integer nroots
         integer nocc(2)
         integer nmo(2)
         integer nfc(2)
         integer nfv(2)
         integer nov(2)
         integer icntr,itmom
         double precision r,cntr(3),tmom(20)
         integer nbf_ao
c
         logical nodezero
c
c        CI vectors are written out as X+Y and X-Y
         integer g_x(2),g_y(2) 
         integer g_temp(2) ! scratch space
c
         integer inntsize,ddblsize,logisize,ok
c
         character*32 pname
         pname = 'get_transden: '
c
c        preliminaries
         inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
         ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
         logisize=MA_sizeof(MT_LOG,1,MT_BYTE)
         ok = 0
         call ga_sync()
         nodezero=(ga_nodeid().eq.0)
         if (.not.bas_numbf(basis,nbf_ao))
     &      call errquit(pname//'bas_numbf failed',0,0)
c
c        initialization
         do i=1,ipol
	   call ga_zero(g_tdens(i))
         end do ! ipol
         do icntr=1,3
           cntr(icntr)=0.0d0
         enddo
         do itmom=1,20
           tmom(itmom)=0.0d0
         enddo
c
c        Read header information from civecs file
         if (nodezero) then      
           open(unit=69,file=filename,form='unformatted',
     &          status='unknown',err=1000)
            rewind(69)
            read(69,err=1001) tda
            read(69,err=1001) ipol
            read(69,err=1001) nroots
            read(69,err=1001) nocc(1),nocc(2)
            read(69,err=1001) nmo(1),nmo(2)
            read(69,err=1001) nfc(1),nfc(2)
            read(69,err=1001) nfv(1),nfv(2)
            read(69,err=1001) nov(1),nov(2)
            read(69,err=1001)
c
            if (ipol.eq.1) nocc(2)=0
            if (ipol.eq.1) nmo(2)=0
            if (ipol.eq.1) nfc(2)=0
            if (ipol.eq.1) nfv(2)=0
c
            do i=1,ipol
             nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
            end do ! ipol
            if (ipol.eq.1) nov(2)=0
           close(unit=69,status='keep',err=1002) ! file
           ok = 1
         end if ! nodezero
c
c        broadcast status and variables to other nodes
         call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, tda, logisize, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT,ipol, inntsize, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT,nroots, inntsize, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, nocc, inntsize*2, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, nmo, inntsize*2, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, nfc, inntsize*2, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, nfv, inntsize*2, 0)
         call ga_brdcst(Msg_Vec_NMO+MSGINT, nov, inntsize*2, 0)
         call ga_sync()
c
c        Set up X, Y vectors and transition density matrices
         do i=1,ipol
            if (.not.ga_create(mt_dbl,nov(i),nroots,'X vector',
     1        -1,-1,g_x(i))) call errquit
     2       (pname//'failed to create g_x',0, GA_ERR)
             call ga_zero(g_x(i))
             if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'temp',
     1        -1,-1,g_temp(i))) call errquit
     2       (pname//'failed to create g_temp',0, GA_ERR)
             call ga_zero(g_temp(i))
c
             if (.not.tda) then
              if (.not.ga_create(mt_dbl,nov(i),nroots,'Y vector',
     1          -1,-1,g_y(i))) call errquit
     2          (pname//'failed to create g_y',0, GA_ERR)
              call ga_zero(g_y(i))
             end if ! .not. tda
         end do  ! ipol
c
c        Read remainder of the civecs file
         if (nodezero) then      
            open(unit=69,file=filename,form='unformatted',
     &                status='unknown',err=1000)
            rewind(69)
            read(69,err=1001) tda
            read(69,err=1001) ipol
            read(69,err=1001) nroots
            read(69,err=1001) nocc(1),nocc(2)
            read(69,err=1001) nmo(1),nmo(2)
            read(69,err=1001) nfc(1),nfc(2)
            read(69,err=1001) nfv(1),nfv(2)
            read(69,err=1001) nov(1),nov(2)
            read(69,err=1001)
         end if ! nodezero
c
         do n = 1,nroots
c
            if (nodezero) then
               read(69) r  ! energy of root
               read(69) r  ! s2
            end if
c
            do i=1,ipol
c
c            Allocate memory
             if (.not.ma_push_get(mt_dbl,nov(i),"slice",l_trials,
     &          k_trials))
     &          call errquit(trim(pname)//"failed to alloc slice",0,0)
c
             if (.not.tda) then
               if (nodezero) then
                call sread(69,dbl_mb(k_trials),nov(i))
                call ga_put(g_x(i),1,nov(i),n,n,dbl_mb(k_trials),nov(i))
                call sread(69,dbl_mb(k_trials),nov(i))
                call ga_put(g_y(i),1,nov(i),n,n,dbl_mb(k_trials),nov(i))
               end if ! nodezero
             else 
               if (nodezero) then
                call sread(69,dbl_mb(k_trials),nov(i))
                call ga_put(g_x(i),1,nov(i),n,n,dbl_mb(k_trials),nov(i))
               end if ! nodezero
             end if  !tda
c
c            Deallocate memory
             if (.not.ma_pop_stack(l_trials))
     &          call errquit(trim(pname)//"failed to pop stack",0,0)
c
            end do ! ipol
         end do ! nroots
c
         if (nodezero) close(unit=69,status='keep',err=1002) ! file
c
         call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
         call ga_sync()
c
c        Separate out X and Y vectors if the calculations are not TDA
c        Nothing needs to be done for TDA as Y is zero
         if (.not.tda) then
           do i=1,ipol
              call ga_add(0.5d0,g_x(i), 0.5d0,g_y(i),g_x(i))
              call ga_add(1.0d0,g_x(i),-1.0d0,g_y(i),g_y(i))
           enddo
         end if
c
c        calculate X component of the transition density matrix
         call tddft_transfm(iroot,g_x,g_movecs,nbf_ao,nocc,nmo,
     &             nfc,nfv,ipol,g_temp,1)  ! x-transition density
         do i = 1,ipol
          call multipole_density(basis,cntr,3,g_temp(i),tmom,20)  ! transition moments
          call ga_copy(g_temp(i),g_tdens(i))
         end do
c
c        calculate Y component of the transition density matrix
         if (.not.tda) then
           do i = 1,ipol
                call ga_zero(g_temp(i))
           end do
           call tddft_transfm(iroot,g_y,g_movecs,nbf_ao,nocc,nmo,
     &             nfc,nfv,ipol,g_temp,1)  ! y-transition density
c
c          accumulate the Y component of the transition density matrix
           do i = 1,ipol
              call multipole_density(basis,cntr,3,g_temp(i),tmom,20)  ! transition moments
              call ga_add(1.d0,g_tdens(i),1.d0,g_temp(i),g_tdens(i))
           end do
         end if  ! tda
c
         if (ipol.eq.1) then
          do i=1,20
            tmom(i)=tmom(i)*dsqrt(2.0d0)
          enddo
         end if 
c
         if (ga_nodeid().eq.0) then
                write(luout,*) " *** tmom(2)***: ", tmom(2)
                write(luout,*) " *** tmom(3)***: ", tmom(3)
                write(luout,*) " *** tmom(4)***: ", tmom(4)
         end if
c
c        symmetrize the transition density matrix
         do i = 1,ipol
             call ga_symmetrize(g_tdens(i))
         enddo
c
c        calculate total, spin density or individual components
         If (ipol.eq.2) Then
            If (Spin.eq.'TOTAL') Then
              Call GA_dAdd(1.d0,g_tdens(1),1.d0,g_tdens(2),g_tdens(1))
            Else If (Spin.eq.'SPINDENS') Then
              Call GA_dAdd(1.d0,g_tdens(1),-1.d0,g_tdens(2),g_tdens(1))
            Else If (Spin.eq.'ALPHA') Then
            Else If (Spin.eq.'BETA') Then
               Call GA_Copy(g_tdens(2),g_tdens(1))
            End If
         End If ! ipol check
c
c        cleanup
         do i=1,ipol
            if (.not.ga_destroy(g_x(i))) call errquit
     2       (pname//'failed to destroy g_x',0, GA_ERR)
            if (.not.ga_destroy(g_temp(i))) call errquit
     2       (pname//'failed to destroy g_temp',0, GA_ERR)
            if (.not.tda) then
              if (.not.ga_destroy(g_y(i))) call errquit
     2          (pname//'failed to destroy g_y',0, GA_ERR)
            end if ! tda
         enddo ! ipol
c
      return
c
 1000 call errquit(pname//'failed to open file',0,-1)
 1001 call errquit(pname//'failed to read file',0,-1)
 1002 call errquit(pname//'failed to close file',0,-1)
c
      end
c $Id$
