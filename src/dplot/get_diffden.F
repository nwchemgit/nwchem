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
         character*256 filename
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
         integer naoc(2)     ! the number of active occupied orbitals
         integer nav(2)      ! the number of active virtual orbitals
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
         pname = 'get_diffden: '
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
c        Get the number of active occupied and virtuals
         do i = 1,ipol
            naoc(i)= (nocc(i)-nfc(i))       
            nav(i) = (nmo(i)-nocc(i)-nfv(i))
         enddo
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
c
            do n = 1,nroots
             if (tda) then
               read(69) r  ! energy of root
               read(69) r  ! s2_save(n)
               do i=1,ipol
                 do m=1,nov(i)
                   read(69) r
                   call ga_put(g_x(i),m,m,n,n,r,1)
                 end do
               end do ! ipol
             else   ! full tddft
               read(69) r  ! energy of root
               read(69) r  ! s2_save(n)
               do i=1,ipol
                 do m=1,nov(i)
                   read(69) r
                   call ga_put(g_x(i),m,m,n,n,r,1)
                 end do
               end do ! ipol
c
               do i=1,ipol
                 do m=1,nov(i)
                   read(69) r
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
c        calculate unrelaxed difference density matrix T
         call compute_T(tda,ipol,nbf_ao,nmo,naoc,nav,iroot,g_movecs,
     &          g_x,g_y,g_ddens)
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
c     Compute the unrelaxed difference density matrix T
c
      subroutine compute_T(tda,ipol,nbf_ao,nmo,naoc,nav,iroot,g_movecs,
     &     g_xpy,g_xmy,g_tao)
c
c     Eq. 19: Furche & Ahlrich JCP (2002)
c     Tijs = -1/2 Sum_a [(X+Y)ias(X+Y)jas+(X-Y)ias(X-Y)jas] (occ-occ part)
c     Tabs =  1/2 Sum_i [(X+Y)ias(X+Y)ibs+(X-Y)ias(X-Y)ibs] (vir-vir part)
c     Tias = Tais = 0  (occ-vir = zero for unrelaxed difference density)
c
      implicit none
c
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
c
c     Input:
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nbf_ao    ! number of ao basis functions
      integer nmo(2)    ! number of molecular orbitals
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer iroot     ! the root of interest
      integer g_movecs(2)  ! molecular orbitals
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
c     Output:
      integer g_tao(2)   ! global array handle for T
c
c     Local:
      integer g_tmo(2)  ! global array handle for T
      integer alo(3)    ! lower chunk limits on A
      integer ahi(3)    ! upper chunk limits on A
      integer blo(3)    ! lower chunk limits on B
      integer bhi(3)    ! upper chunk limits on B
      integer clo(3)    ! lower chunk limits on C
      integer chi(3)    ! upper chunk limits on C
      integer ip        ! counter on polarization
      logical oroot
      integer idim(3)   ! dimension array for the nga calls
      integer ichunk(3) ! chunk size array for the nga calls
      integer i         ! general index
      integer g_temp    ! temporary array
c
      character*32 pname
      pname = 'compute_T: '
c
c     Create T
      do i = 1, ipol
        idim(1)   = 1  ! for one root
        idim(2)   = naoc(i)+nav(i)
        idim(3)   = naoc(i)+nav(i)
        ichunk(1) =  1
        ichunk(2) = -1
        ichunk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'T',ichunk,g_tmo(i)))
     +    call errquit(pname//'failed to create g_tmo',0,GA_ERR)
      enddo
c
c     Loop over spins
      do ip = 1, ipol
        call ga_zero(g_tmo(ip))
c
c       Do the occupied-occupied block
        alo(1) = iroot
        ahi(1) = iroot
        alo(2) = 1
        ahi(2) = naoc(ip)
        alo(3) = 1
        ahi(3) = nav(ip)
        blo(1) = iroot
        bhi(1) = iroot
        blo(2) = 1
        bhi(2) = nav(ip)
        blo(3) = 1
        bhi(3) = naoc(ip)
        clo(1) = iroot
        chi(1) = iroot
        clo(2) = 1
        chi(2) = naoc(ip)
        clo(3) = 1
        chi(3) = naoc(ip)
c
        if (tda) then
c
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)jas = Xjas
c (X-Y)ias = Xias , (X-Y)jas = Xjas
c Tij = - sum_a Xias*Xjas 
c
            call nga_matmul_patch('n','t',-0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
            call nga_matmul_patch('n','t',-0.5d0,1.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
          else
c (X+Y)ias(X+Y)jas 
            call nga_matmul_patch('n','t',-0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
c (X-Y)ias(X-Y)jas
            call nga_matmul_patch('n','t',-0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
        endif ! tda
c
c       Do the virtual-virtual block
        alo(1) = iroot
        ahi(1) = iroot
        alo(2) = 1
        ahi(2) = nav(ip)
        alo(3) = 1
        ahi(3) = naoc(ip)
        blo(1) = iroot
        bhi(1) = iroot
        blo(2) = 1
        bhi(2) = naoc(ip)
        blo(3) = 1
        bhi(3) = nav(ip)
        clo(1) = iroot
        chi(1) = iroot
        clo(2) = naoc(ip)+1
        chi(2) = naoc(ip)+nav(ip)
        clo(3) = naoc(ip)+1
        chi(3) = naoc(ip)+nav(ip)
c
        if (tda) then
c
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)ibs = Xibs
c (X-Y)ias = Xias , (X-Y)ibs = Xibs
c Tab = sum_i Xias*Xibs 
c
            call nga_matmul_patch('t','n',0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
            call nga_matmul_patch('t','n',0.5d0,1.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
          else
c (X+Y)ias(X+Y)ibs 
            call nga_matmul_patch('t','n',0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
c (X-Y)ias(X-Y)ibs 
            call nga_matmul_patch('t','n',0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tmo(ip),clo,chi)
        endif  ! tda
      enddo ! ip = 1, ipol
c
c     Compute T in the AO basis 
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Temp',-1,-1,g_temp))
     &    call errquit(pname//'failed to create g_temp',0,GA_ERR)
      do i=1,ipol
          call ga_matmul_patch('n','t',1.0d0/(3-ipol),0.0d0,
     &         g_tmo(i),1,nmo(i),1,nmo(i),
     &         g_movecs(i),1,nmo(i),1,nbf_ao,
     &         g_temp,1,nmo(i),1,nbf_ao)
          call ga_matmul_patch('n','n',1.0d0,0.0d0,
     &         g_movecs(i),1,nbf_ao,1,nmo(i),
     &         g_temp,1,nmo(i),1,nbf_ao,
     &         g_tao(i),1,nbf_ao,1,nbf_ao)
      enddo ! ipol
c
c     Clear T,...
      do i = 1,ipol
        if (.not.ga_destroy(g_tmo(i)))
     +    call errquit(pname//'failed to destroy g_tmo',0,GA_ERR)
      enddo
      if (.not.ga_destroy(g_temp)) 
     &    call errquit(pname//'failed to destroy g_temp',0, GA_ERR)
c
      end
c $Id$
