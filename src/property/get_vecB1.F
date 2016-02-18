      subroutine get_vecB1_opt1(
     &                     g_vecB1,    ! out:
     &                     g_rhs,      ! in : g_rhs vector
     &                     vectors,    ! in : MO vectors
     &                     nbf,        ! in : nr. basis functions
     &                     nmo,        ! in : nr. MOs
     &                     npol,       ! in : nr polarizations
     &                     nocc,       ! in : nr. occupied MOs
     &                     nvirt,      ! in : nr. virtual  MOs
     &                     debug)      ! in : = .true. allow debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         
c        To be used in fiao_b1_movecs.F, it uses less memory

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer npol ! nr of polarizations =1 (  restricted) 
                   !                     =2 (unrestricted) calc    
      integer g_vecB1(npol)
      integer vectors(npol),
     &        g_rhs,g_u
      integer ifld,nbf,nmo,disp
      logical debug
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        plo(3), phi(3),
     &        qlo(3), qhi(3)
c ------- Added for unrestricted calc ----- START
      integer ndir,    ! nr directions (x,y,z)
     &        ispin,
     &        nocc(npol),  ! store nr occupations     
     &        nvirt(npol)
c ------- Added for unrestricted calc ----- END
      ndir=3 ! = nr directions (x,y,z)
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)] valid only for closed shells
c     For unrestricted calc. or open shell systems
c     D1=D1_A + D1_B
c     D1_A=[(C1*C0+) - (C0*C1+)]_A
c     D1_B=[(C1*C0+) - (C0*C1+)]_B
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf  
       ahi(2) = nocc(ispin)
       ahi(3) = 1
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('giao_b1: nga_create failed g_u',0,GA_ERR)
       call ga_zero(g_vecB1(ispin)) ! initialize 
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1    
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 1
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nmo 
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nocc(ispin)
       dlo(1) = 1
       dhi(1) = nbf
       dlo(2) = 1
       dhi(2) = nocc(ispin)
       dlo(3) = 1
       dhi(3) = 1
       disp=nocc(1)*nvirt(1)*(ispin-1)
       plo(1) = disp+1
       phi(1) = disp+nocc(ispin)*nvirt(ispin)
       qlo(1) = nocc(ispin)+1
       qhi(1) = nmo
       qlo(2) = 1
       qhi(2) = nocc(ispin)
       qlo(3) = 1
       qhi(3) = 1
       do ifld=1,ndir
        clo(3) = ifld
        chi(3) = ifld
c     Make C1
c ======== Including g_u_ov (g_rhs --> g_u) ==== START
        plo(2) = ifld
        phi(2) = ifld
        call ga_zero(g_u)
        call nga_copy_patch('n',g_rhs,plo,phi,
     &                          g_u  ,qlo,qhi) 
c ======== Including g_u_ov (g_rhs --> g_u) ==== END
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                vectors(ispin),blo,bhi,  
     &                g_u           ,alo,ahi,
     &                g_vecB1(ispin),clo,chi)
c        This patch of g_vecB1 now has the perturbed MO
c        coefficients. let's print them for debug purposes:
         if (debug) then
           write (luout,*) 'fiao_b1: perturbed C, direction ',ifld
           call nga_print_patch(g_vecB1(ispin),clo,chi,1)
         end if
       enddo ! end-loop-ifld
      if (.not.ga_destroy(g_u)) call 
     &    errquit('fiao_b1: ga_destroy failed g_d1',0,GA_ERR)
      enddo ! end-loop-ispin
      return
      end

      subroutine get_vecB1_opt2(
     &                     g_vecB1,    ! out:
     &                     g_rhs,      ! in : g_rhs vector (occ-virt of g_u)
     &                     g_u_oo,     ! in : occ-occ of g_u
     &                     vectors,    ! in : MO vectors
     &                     nbf,        ! in : nr. basis functions
     &                     nmo,        ! in : nr. MOs
     &                     npol,       ! in : nr polarizations
     &                     nocc,       ! in : nr. occupied MOs
     &                     nvirt,      ! in : nr. virtual  MOs
     &                     debug)      ! in : = .true. allow debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         
c        To be used in giao_b1_movecs.F, it uses less memory

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer npol ! nr of polarizations =1 (  restricted) 
                   !                     =2 (unrestricted) calc    
      integer g_vecB1(npol)
      integer vectors(npol),
     &        g_rhs,g_u,g_u_oo(npol)
      integer ifld,nbf,nmo,disp
      logical debug
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        plo(3), phi(3),
     &        qlo(3), qhi(3),
     &        rlo(3), rhi(3),
     &        slo(3), shi(3)
c ------- Added for unrestricted calc ----- START
      integer ndir,    ! nr directions (x,y,z)
     &        ispin,
     &        nocc(npol),  ! nr. occupied MOs    
     &        nvirt(npol)  ! nr. virtual  MOs
c ------- Added for unrestricted calc ----- END
      ndir=3 ! = nr directions (x,y,z)
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)] valid only for closed shells
c     For unrestricted calc. or open shell systems
c     D1=D1_A + D1_B
c     D1_A=[(C1*C0+) - (C0*C1+)]_A
c     D1_B=[(C1*C0+) - (C0*C1+)]_B
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nmo
       ahi(2) = nocc(ispin)
       ahi(3) = 1
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('giao_b1: nga_create failed g_u',0,GA_ERR)
       call ga_zero(g_vecB1(ispin)) ! initialize 
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 1
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nmo 
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nocc(ispin)
       dlo(1) = 1
       dhi(1) = nbf
       dlo(2) = 1
       dhi(2) = nocc(ispin)
       dlo(3) = 1
       dhi(3) = 1
       disp=nocc(1)*nvirt(1)*(ispin-1)
       plo(1) = disp+1
       phi(1) = disp+nocc(ispin)*nvirt(ispin)
       qlo(1) = nocc(ispin)+1
       qhi(1) = nmo
       qlo(2) = 1
       qhi(2) = nocc(ispin)
       qlo(3) = 1
       qhi(3) = 1
       rlo(1) = 1
       rhi(1) = nocc(ispin)
       rlo(2) = 1
       rhi(2) = nocc(ispin)
       rlo(3) = 1
       rhi(3) = 1
       slo(1) = 1
       shi(1) = nocc(ispin)
       slo(2) = 1
       shi(2) = nocc(ispin)
       do ifld=1,ndir
        clo(3) = ifld
        chi(3) = ifld
c     Make C1
c ======== Including g_u_ov (g_rhs --> g_u) ==== START
        plo(2) = ifld
        phi(2) = ifld
        call ga_zero(g_u)
        call nga_copy_patch('n',g_rhs,plo,phi,
     &                          g_u  ,qlo,qhi)
c ======== Including g_u_ov (g_rhs  --> g_u) ==== END
c ======== Including g_u_oo (g_u_oo --> g_u) ==== START
        slo(3)= ifld
        shi(3)= ifld
        call nga_copy_patch('n',g_u_oo(ispin),slo,shi,
     &                          g_u          ,rlo,rhi)
c ======== Including g_u_oo (g_u_oo --> g_u) ==== END
         call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                vectors(ispin),blo,bhi,  
     &                g_u           ,alo,ahi,
     &                g_vecB1(ispin),clo,chi)
c        This patch of g_vecB1 now has the perturbed MO
c        coefficients. let's print them for debug purposes:
         if (debug) then
           write (luout,*) 'fiao_b1: perturbed C, direction ',ifld
           call nga_print_patch(g_vecB1(ispin),clo,chi,1)
         end if
       enddo ! end-loop-ifld
      if (.not.ga_destroy(g_u)) call 
     &    errquit('fiao_b1: ga_destroy failed g_u',0,GA_ERR)
      enddo ! end-loop-ispin
      return
      end

      subroutine update_rhs_angmomAO(
     &               g_rhs,  !in/out:
     &               vectors,!in : MO vectors
     &               nocc,   !in : nr. occup   MOs
     &               nvirt,  !in : nr. virtual MOs   
     &               npol,   !in : nr. polarizations
     &               nbf,    !in : nr. basis functions
     &               nmo,    !in : nr. MOs
     &               basis,  !in : basis handle
     &               debug)  !in : logical var for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "prop.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "bgj.fh"
#include "case.fh"
      integer basis
      integer npol
      integer vectors(npol)
      integer nocc(npol),nvirt(npol),
     &        iocc,disp,ispin,
     &        ndir,shift,
     &        alo(3),ahi(3),
     &        blo(3),bhi(3)
      integer nbf,nmo
      integer g_rhs,g_s10,g_s10_1
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      logical oskel,debug
      integer nat
      parameter (nat=1)
      external giao_aotomo
      oskel=.false.
      ndir=3 ! nr. directions (x,y,z)
c     Get S10 in GA and transform to MO set (virt,occ)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10_1)) 
     &    call 
     &    errquit('fiao_b1: nga_create failed g_s01',0,GA_ERR)
c     Get H10 in GA, using g_s10 array
      call ga_zero(g_s10_1)
      ahi(3) = npol*ndir
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10)) 
     &    call 
     &    errquit('fiao_b1: nga_create failed g_s01',0,GA_ERR)
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,
     &                   g_s10_1,'angmom', ! out: g_s10
     &                   origin,nat,oskel)
      call ga_scale (g_s10_1,0.5d0)
      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '--------angmom-0.5 ------ START'
        call ga_print(g_s10_1)
       if (ga_nodeid().eq.0)
     &  write(*,*) '--------angmom-0.5 ------ END'
      endif ! end-if-debug

c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs
c --------- g_s10_1 --> g_s10 --------- START
       blo(1) = 1
       bhi(1) = nbf ! nmo fix lindep 05-02-12
       blo(2) = 1
       bhi(2) = nbf ! nmo fix lindep 05-02-12
       blo(3) = 1
       bhi(3) = ndir
      do ispin=1,npol  
       disp=ndir*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf ! nmo fix lindep 05-02-12
       alo(2) = 1
       ahi(2) = nbf ! nmo fix lindep 05-02-12
       alo(3) = disp+1
       ahi(3) = disp+ndir  
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin
c --------- g_s10_1 --> g_s10 --------- END
      if (debug) then
       if (ga_nodeid().eq.0)
     & write(*,*) '--------angmom-0.5x2 ------ START'
       call ga_print(g_s10)
       if (ga_nodeid().eq.0)
     & write(*,*) '--------angmom-0.5x2------ END'
       if (ga_nodeid().eq.0)
     & write(*,*) '--------vectors ------ START'
       call ga_print(vectors(1))
       if (ga_nodeid().eq.0)
     & write(*,*) '--------vectors ------ END'
       if (ga_nodeid().eq.0) then
        write(*,10) nocc(1),nocc(2),
     &              nvirt(1)  ,nvirt(2),
     &              npol,ndir,nbf
 10    format('nclos=(',i3,',',i3,') ',
     &        'nvirt=(',i3,',',i3,') ',
     &        ' (npol,ndir,nbf)=(',i3',',i3,',',i3,')')
       endif
      endif ! end-if-debug
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,ndir,nbf) 
      if (debug) then
       if (ga_nodeid().eq.0)
     & write(*,*) '-----angmom-s10-aotomo ------ START'
       call ga_print(g_s10)
       if (ga_nodeid().eq.0)
     & write(*,*) '-----angmom-s10-aotomo ------ END'
      endif ! end-if-debug
      do ispin=1,npol
       shift=ndir*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+ndir
       disp=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = ndir
       call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs,blo,bhi)
      enddo ! end-loop-ispin
c
c     Cleanup g_s10 as we do not need it right now
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('fiao_b1: ga_destroy failed g_s10',0,GA_ERR)
      if (.not.ga_destroy(g_s10_1)) call 
     &    errquit('fiao_b1: ga_destroy failed g_s10',0,GA_ERR)
      return
      end
c $Id$
