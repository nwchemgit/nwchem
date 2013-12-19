c
c        get the unrelaxed difference density
         subroutine get_diffden(spin,iroot,ipol,filename,basis,
     &        g_movecs, g_ddens)
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
         character*80 filename
         integer basis            ! AO basis set handle
         integer g_movecs(ipol)   ! MO vectors
         integer g_ddens(ipol)    ! Unrelaxed difference density
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
         integer g_x(2)  ! X+Y
         integer g_y(2)  ! X-Y
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
	   call ga_zero(g_ddens(i))
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
           open(unit=69,file=filename,form='formatted',
     &          status='unknown',err=1000)
            rewind(69)
            read(69,*,err=1001) tda
            read(69,*,err=1001) ipol
            read(69,*,err=1001) nroots
            read(69,*,err=1001) nocc(1),nocc(2)
            read(69,*,err=1001) nmo(1),nmo(2)
            read(69,*,err=1001) nfc(1),nfc(2)
            read(69,*,err=1001) nfv(1),nfv(2)
            read(69,*,err=1001) nov(1),nov(2)
            read(69,*,err=1001)
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
c        Set up memory for X(X+Y), Y(X-Y) vectors
         do i=1,ipol
            if (.not.ga_create(mt_dbl,nov(i),nroots,'X vector',
     1        -1,-1,g_x(i))) call errquit
     2       (pname//'failed to create g_x',0, GA_ERR)
            call ga_zero(g_x(i))
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
          open(unit=69,file=filename,form='formatted',
     &                status='unknown',err=1000)
            rewind(69)
            read(69,*,err=1001) tda
            read(69,*,err=1001) ipol
            read(69,*,err=1001) nroots
            read(69,*,err=1001) nocc(1),nocc(2)
            read(69,*,err=1001) nmo(1),nmo(2)
            read(69,*,err=1001) nfc(1),nfc(2)
            read(69,*,err=1001) nfv(1),nfv(2)
            read(69,*,err=1001) nov(1),nov(2)
            read(69,*,err=1001)
c
            do n = 1,nroots
             if (tda) then
               read(69,*) r  ! energy of root
               read(69,*) r  ! s2_save(n)
               do i=1,ipol
                 do m=1,nov(i)
                   read(69,*) r
                   call ga_put(g_x(i),m,m,n,n,r,1)
                 end do
               end do ! ipol
             else   ! full tddft
               read(69,*) r  ! energy of root
               read(69,*) r  ! s2_save(n)
               do i=1,ipol
                 do m=1,nov(i)
                   read(69,*) r
                   call ga_put(g_x(i),m,m,n,n,r,1)
                 end do
               end do ! ipol
c
               do i=1,ipol
                 do m=1,nov(i)
                   read(69,*) r
                   call ga_put(g_y(i),m,m,n,n,r,1)
                 end do
               end do ! ipol
             end if  ! tda
            end do ! nroots
          close(unit=69,status='keep',err=1002) ! file
          ok = 1
         end if  ! nodezero
c
         call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
         call ga_sync()
c
c        calculate unrelaxed difference density matrix
         call calc_diffden(tda,iroot,g_x,g_y,
     &       g_movecs,nbf_ao,nocc,nmo,nfc,nfv,ipol,g_ddens)
c
c        symmetrize the difference density matrix
         do i = 1,ipol
             call ga_symmetrize(g_ddens(i))
         enddo
c
c        calculate total, spin density or individual components
         If (ipol.eq.2) Then
            If (Spin.eq.'TOTAL') Then
              Call GA_dAdd(1.d0,g_ddens(1),1.d0,g_ddens(2),g_ddens(1))
            Else If (Spin.eq.'SPINDENS') Then
              Call GA_dAdd(1.d0,g_ddens(1),-1.d0,g_ddens(2),g_ddens(1))
            Else If (Spin.eq.'ALPHA') Then
            Else If (Spin.eq.'BETA') Then
               Call GA_Copy(g_ddens(2),g_ddens(1))
            End If
         End If ! ipol check
c
c        cleanup
         do i=1,ipol
            if (.not.ga_destroy(g_x(i))) call errquit
     2       (pname//'failed to destroy g_x',0, GA_ERR)
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
c
c     calculate the unrelaxed difference density 
      subroutine calc_diffden(tda,iroot,g_xpy,g_xmy,
     &    g_movecs,nbf_ao,nocc,nmo,nfc,nfv,ipol,g_tao)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
c
      double precision thresh  ! Threshold for printing
      parameter (thresh=5.0d-2)
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
      double precision mijtol  ! Transition moment tolerance
      parameter (mijtol=1.d-04)
      double precision au2debye  ! Convert au to debye
      parameter (au2debye=2.541766d0)
c
      logical tda
      integer iroot          
      integer g_xpy(2)           ! X+Y
      integer g_xmy(2)           ! X-Y
      integer g_movecs(2)        ! mo vectors
      integer nbf_ao             ! Number of AO basis functions
      integer nocc(2)            ! Number of occupied orbitals
      integer nmo(2)             ! Number of orbitals
      integer nfc(2),nfv(2)      ! Number of frozen cores/virtuals
      integer ipol               ! =1 (restricted), =2 (unrestricted)
      integer g_tao(2)           ! GA handle for AO density matrix
c
      integer i
      integer g_xpy1(2)          ! X+Y for 1 root
      integer g_xmy1(2)          ! X-Y for 1 root
      integer g_tmo(2) 
      integer g_temp             ! temp array
      integer nov(2)
      character*32 pname
      logical nodezero
c
      double precision zero, one, two, three, half, third
      parameter (zero = 0d0, one=1d0, two=2d0, three=3d0,
     &   half=one/two, third=one/three)
c
c     Preliminaries
      nodezero=(ga_nodeid().eq.0)
      pname="calc_diffden: "
c
c     Number of occupied virtual pairs
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
c
c     Set up the memory
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nmo(i)-nocc(i)-nfv(i),
     &      nocc(i)-nfc(i),'X+Y for 1 root',-1,-1,g_xpy1(i))) 
     &   call errquit(pname//'failed to create g_xpy1',0, GA_ERR)
        call ga_zero(g_xpy1(i))
        if (.not.tda) then
         if (.not.ga_create(mt_dbl,nmo(i)-nocc(i)-nfv(i),
     &      nocc(i)-nfc(i),'X-Y for 1 root',-1,-1,g_xmy1(i)))
     &             call errquit
     &    (pname//'failed to create g_xmy1',0, GA_ERR)
          call ga_zero(g_xmy1(i))
        endif
        if (.not.ga_create(mt_dbl,nmo(i),nmo(i),
     &      'T in the MO basis',-1,-1,g_tmo(i))) call errquit
     &    (pname//'failed to create g_tmo',0, GA_ERR)
        call ga_zero(g_tmo(i))
      enddo
c
c     Copy one vector
      do i=1,ipol
        call ga_copy_patch('n',g_xpy(i),1,nov(i),iroot,iroot,
     &    g_xpy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i))
        if (.not.tda) then
          call ga_copy_patch('n',g_xmy(i),1,nov(i),iroot,iroot,
     &       g_xmy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i))
        endif
        call ga_zero(g_tmo(i))
      enddo ! ipol
c
c     Assemble the unrelaxed difference density in MO space
      do i=1,ipol
c
c         Compute the occupied-occupied block
          call ga_matmul_patch('t','n',-0.5d0,0.0d0,
     &         g_xpy1(i),1,nocc(i)-nfc(i),1,nmo(i)-nocc(i)-nfv(i),
     &         g_xpy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i),
     &         g_tmo(i),nfc(i)+1,nocc(i),nfc(i)+1,nocc(i))
          if (.not.tda) then
            call ga_matmul_patch('t','n',-0.5d0,1.0d0,
     &           g_xmy1(i),1,nocc(i)-nfc(i),1,nmo(i)-nocc(i)-nfv(i),
     &           g_xmy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i),
     &           g_tmo(i),nfc(i)+1,nocc(i),nfc(i)+1,nocc(i))
          endif
c
c         Compute the virtual-virtual block
          call ga_matmul_patch('n','t',0.5d0,0.0d0,
     &         g_xpy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i),
     &         g_xpy1(i),1,nocc(i)-nfc(i),1,nmo(i)-nocc(i)-nfv(i),
     &         g_tmo(i),nocc(i)+1,nmo(i)-nfv(i),
     &                nocc(i)+1,nmo(i)-nfv(i))
          if (.not.tda) then
            call ga_matmul_patch('n','t',0.5d0,1.0d0,
     &           g_xmy1(i),1,nmo(i)-nocc(i)-nfv(i),1,nocc(i)-nfc(i),
     &           g_xmy1(i),1,nocc(i)-nfc(i),1,nmo(i)-nocc(i)-nfv(i),
     &           g_tmo(i),nocc(i)+1,nmo(i)-nfv(i),
     &                  nocc(i)+1,nmo(i)-nfv(i))
          endif
      enddo ! ipol
c
c     Scratch space
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'temp',-1,-1,g_temp))
     &    call errquit(pname//'failed to create g_temp',0,GA_ERR)
c
c     Compute the unrelaxed difference density in the AO basis
      do i=1,ipol
          call ga_zero(g_temp)
          call ga_matmul_patch('n','t',1.0d0/(3-ipol),0.0d0,
     &         g_tmo(i),1,nmo(i),1,nmo(i),
     &         g_movecs(i),1,nmo(i),1,nbf_ao,
     &         g_temp,1,nmo(i),1,nbf_ao)
          call ga_matmul_patch('n','n',1.0d0,0.0d0,
     &         g_movecs(i),1,nbf_ao,1,nmo(i),
     &         g_temp,1,nmo(i),1,nbf_ao,
     &         g_tao(i),1,nbf_ao,1,nbf_ao)
      enddo
c
c     Delete memory
      do i=1,ipol
        if (.not.ga_destroy(g_xpy1(i))) call errquit
     &    (pname//'failed to destroy g_xpy1',0, GA_ERR)
        if (.not.tda) then
         if (.not.ga_destroy(g_xmy1(i))) call errquit
     &    (pname//'failed to destroy g_xmy1',0, GA_ERR)
        end if ! tda
        if (.not.ga_destroy(g_tmo(i))) call errquit
     &    (pname//'failed to destroy g_tmo',0, GA_ERR)
      enddo
      if (.not.ga_destroy(g_temp)) call errquit
     &    (pname//'failed to destroy g_temp',0, GA_ERR)
      if (nodezero) call util_flush(LuOut)
c
      return
      end
c $Id: $
