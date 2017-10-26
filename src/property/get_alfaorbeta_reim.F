      subroutine get_alfaorbeta_reim(
     &                AorBre,      ! in/out: alpha or beta real part
     &                AorBim,      ! in/out: alpha or beta im   part
     &                g_vecE1,     ! in : 1st-order pert vec RE
     &                g_vecE1_im,  ! in : 1st-order pert vec IM
     &                g_dipEM,     ! in : dipole electric or magnetic
     &                g_vectors,   ! in : MO vectors
     &                idir,        ! in : = 1,2,3=x,y,z directions
     &                iresp,       ! in : = 1,2,3
     &                coeffre,     ! in : coeff for real part
     &                coeffim,     ! in : coeff for imag part
     &                caseAO,      ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &                nbf,         ! in : nr. basis functions
     &                nocc,        ! in : nr. occupied alpha or beta
     &                lifetime,    ! in : logical var for damping
     &                debug,       ! in : logical var for debugging
     &                g_dipel_i,
     &                omega,
     &                g_temp)      ! in : scratch GA array
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
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "util.fh"
      double precision AorBre, ! OUTPUT
     &                 AorBim  ! OUTPUT 
      integer idir,iresp,ndir,
     &        nbf,nocc,
     &        caseAO,ind1,ind2
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3), 
     &        clo(3),chi(3)
      integer g_temp,     ! scratch ga arrays (input)
     &        g_vecE1,    ! IN
     &        g_vecE1_im, ! IN
     &        g_vectors,  ! IN
     &        g_dipEM     ! IN : = g_dipel or g_dipmag
      logical lifetime,debug
      double precision sum,sumre,sumim,
     &                 coeffre,coeffim,trace
      external trace,get_C1MCtrace
c  DIM/QM stuff
      integer g_dipel_i
      double precision omega
      logical lclfld, use_dimqm, ldimqm

      use_dimqm=util_module_avail("dimqm")

      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlclfld(lclfld)
      else
         ldimqm   =  .false.
         lclfld   =  .false.
      endif 
c  DIM/QM stuff END

      if (debug) then
        if (ga_nodeid().eq.0) then
           write(*,2) idir,iresp,caseAO,nbf,
     &                nocc
 2         format('(idir,iresp,caseAO,nbf,nocc)=(',
     &            i3,',',i3,',',i3,',',i3,',',i3,')')
          endif
         if (ga_nodeid().eq.0) then
          write(*,10) idir,iresp
 10       format('---- g_vecE1-re(',i3,',',i3,')-------- START')
         endif
          call ga_print(g_vecE1)
         if (ga_nodeid().eq.0) then
          write(*,11) idir,iresp
 11       format('---- g_vecE1-re(',i3,',',i3,')-------- END')
         endif
         if (lifetime) then
         if (ga_nodeid().eq.0) then
          write(*,12) idir,iresp
 12       format('---- g_vecE1-im(',i3,',',i3,')-------- START')
         endif
          call ga_print(g_vecE1_im)
         if (ga_nodeid().eq.0) then
          write(*,13) idir,iresp
 13       format('---- g_vecE1-im(',i3,',',i3,')-------- END')
         endif
         endif ! end-if-lifetime
c         if (ga_nodeid().eq.0) then
c          write(*,14) idir,iresp
c 14       format('---- g_dipEM(',i3,',',i3,')-------- START')
c         endif
c          call ga_print(g_dipEM)
c         if (ga_nodeid().eq.0) then
c          write(*,15) idir,iresp
c 15       format('---- g_dipEM(',i3,',',i3,')-------- END')
c         endif
c         if (ga_nodeid().eq.0) then
c          write(*,16) idir,iresp
c 16       format('---- g_vectors(',i3,',',i3,')-------- START')
c         endif
c          call ga_print(g_vectors)
c         if (ga_nodeid().eq.0) then
c          write(*,17) idir,iresp
c 17       format('---- g_vectors(',i3,',',i3,')-------- END')
c         endif
      endif ! end-if-debug
 
       call get_C1MCtrace(
     &               sumre,       ! out: trace(transp(vecE1  )*g_temp)
     &               sumim,       ! out: trace(transp(vecE1im)*g_temp)
     &               lifetime,    ! in : =T => returns sumim
     &               g_vecE1,     ! in : 1st-order pert vec RE
     &               g_vecE1_im,  ! in : 1st-order pert vec IM
     &               g_dipEM,     ! in : dipole electric or magnetic
     &               g_vectors,   ! in : MO vectors
     &               idir,        ! in : = 1,2,3=x,y,z directions
     &               iresp,       ! in : = 1,2,3
     &               caseAO,      ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &               nbf,         ! in : nr. basis functions
     &               nocc,        ! in : nr. occupied alpha or beta
     &               debug,       ! in : logical var for debugging
     &               g_temp)      ! in : scratch GA array -> out

c           the factor of two is for the orbital occupations,
c           assuming that ispin is never equal to two     
       AorBre=AorBre+coeffre*sumre
       if (lifetime) 
     & AorBim=AorBim+coeffim*sumim

c  jbecca START
c  In the case of DIM/QM complex local fields, a complex dipole operator
c  now exists. Call the trace routine again with the imaginary dipole
c  portion to get the cross terms to add back into the polarizability.

      if (ldimqm.and.lclfld.and.(omega.gt.1d-8).and.lifetime) then
       call get_C1MCtrace(
     &               sumre,       ! out: trace(transp(vecE1  )*g_temp)
     &               sumim,       ! out: trace(transp(vecE1im)*g_temp)
     &               lifetime,    ! in : =T => returns sumim
     &               g_vecE1,     ! in : 1st-order pert vec RE
     &               g_vecE1_im,  ! in : 1st-order pert vec IM
     &               g_dipel_i,     ! in : dipole electric or magnetic
     &               g_vectors,   ! in : MO vectors
     &               idir,        ! in : = 1,2,3=x,y,z directions
     &               iresp,       ! in : = 1,2,3
     &               caseAO,      ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &               nbf,         ! in : nr. basis functions
     &               nocc,        ! in : nr. occupied alpha or beta
     &               debug,       ! in : logical var for debugging
     &               g_temp)      ! in : scratch GA array -> out

         AorBre=AorBre-coeffim*sumim
         AorBim=AorBim+coeffre*sumre
      endif
c  jbecca END

         
      return
      end

      subroutine get_C1MC(
     &                g_work,      ! out: C(E) M C
     &                g_vecE1,     ! in : 1st-order pert vec RE
     &                g_dipEM,     ! in : dipole electric or magnetic
     &                g_vectors,   ! in : MO vectors
     &                idir,        ! in : = 1,2,3=x,y,z directions
     &                iresp,       ! in : = 1,2,3
     &                caseAO,      ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &                nbf,         ! in : nr. basis functions
     &                nocc,        ! in : nr. occupied alpha or beta
     &                debug,       ! in : logical var for debugging
     &                g_temp)      ! in : scratch GA array
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
#include "rtdb.fh"
#include "stdio.fh"
      integer g_work  ! = C(E) M C
      integer idir,iresp,
     &        nbf,nocc,
     &        caseAO,ind1,ind2
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3), 
     &        clo(3),chi(3)
      integer g_temp,     ! IN: scratch ga arrays (input)
     &        g_vecE1,    ! IN
     &        g_vectors,  ! IN
     &        g_dipEM     ! IN : = g_dipel or g_dipmag
      logical debug
c Note.- (ind1,ind2)=(iresp,1    ) for caseAO=1 (g_dipEM ne g_smat0)
c        (ind1,ind2)=(1    ,iresp) for caseAO=2 (g_dipEM eq g_smat0)
c        (ind1,ind2)=(1    ,idir ) for caseAO=3 (g_dipEM eq g_smat0) in aor_r1_beta_anl
c        (ind1,ind2)=(idir ,1    ) for caseAO=4 (g_dipEM eq g_sket1) in aor_r1_beta_anl
            if      (caseAO .eq. 1) then
             ind1=iresp
             ind2=1
            else if (caseAO .eq. 2) then
             ind1=1
             ind2=iresp
            else if (caseAO .eq. 3) then
             ind1=1
             ind2=idir
            else if (caseAO .eq. 4) then
             ind1=idir
             ind2=1
            else
             call errquit
     &       ('get_C1MC: caseAO ne 1,2,3 or 4',
     &       0,GA_ERR)
            endif
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = ind1 ! pick direction iresp for g_dipEM
            ahi(3) = ind1
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            blo(3) = ind2 
            bhi(3) = ind2     
            clo(1) = 1
            chi(1) = nbf
            clo(2) = 1
            chi(2) = nocc
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,18) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18           format('FA-1::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
             endif
            endif ! end-if-debug

            call ga_zero(g_temp)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_dipEM  ,alo,ahi,
     &         g_vectors,blo,bhi,
     &         g_temp   ,clo,chi)
            if (debug) write (luout,*)
     &         'alfa: h(E) C(0) intermediate complete'
            alo(1) = 1
            ahi(1) = nocc
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = idir
            ahi(3) = idir
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            clo(1) = 1
            chi(1) = nocc
            clo(2) = 1
            chi(2) = nocc
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,19) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 19           format('FA-2::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
             endif
            endif ! end-if-debug
            call ga_zero(g_work)
            if (debug) then
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_vecE1 --------- START'
              call ga_print(g_vecE1)
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_vecE1 --------- END'
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_temp --------- START'
              call ga_print(g_temp)
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_temp --------- END'
            endif ! end-if-debug
            call nga_matmul_patch('t','n',1d0,0d0,
     &                  g_vecE1,alo,ahi,
     &                  g_temp ,blo,bhi,
     &                  g_work ,clo,chi)
            if (debug) then
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_work-inside----- START'
              call ga_print(g_work)
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_work-inside------- END'
            endif ! end-if-debug
      return
      end

      subroutine get_C1MCtrace(
     &                sumre,       ! out: trace(transp(vecE1  )*g_temp)
     &                sumim,       ! out: trace(transp(vecE1im)*g_temp)
     &                lifetime,    ! in : =T => returns sumim
     &                g_vecE1,     ! in : 1st-order pert vec RE
     &                g_vecE1_im,  ! in : 1st-order pert vec IM
     &                g_dipEM,     ! in : dipole electric or magnetic
     &                g_vectors,   ! in : MO vectors
     &                idir,        ! in : = 1,2,3=x,y,z directions
     &                iresp,       ! in : = 1,2,3
     &                caseAO,      ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &                nbf,         ! in : nr. basis functions
     &                nocc,        ! in : nr. occupied alpha or beta
     &                debug,       ! in : logical var for debugging
     &                g_temp)      ! in : scratch GA array
c Note.- g_temp= g_dipEM * g_vectors
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
#include "rtdb.fh"
#include "stdio.fh"
      integer idir,iresp,
     &        nbf,nocc,
     &        caseAO,ind1,ind2
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3), 
     &        clo(3),chi(3)
      integer g_temp,     ! IN: scratch ga arrays (input)
     &        g_vecE1,    ! IN
     &        g_vecE1_im, ! IN
     &        g_vectors,  ! IN
     &        g_dipEM     ! IN : = g_dipel or g_dipmag
      double precision trace,sumre,sumim
      external trace
      logical lifetime,debug
c Note.- (ind1,ind2)=(iresp,1    ) for caseAO=1 (g_dipEM ne g_smat0)
c        (ind1,ind2)=(1    ,iresp) for caseAO=2 (g_dipEM eq g_smat0)
c        (ind1,ind2)=(1    ,idir ) for caseAO=3 (g_dipEM eq g_smat0) in aor_r1_beta_anl
c        (ind1,ind2)=(idir ,1    ) for caseAO=4 (g_dipEM eq g_sket1) in aor_r1_beta_anl
            if      (caseAO .eq. 1) then
             ind1=iresp
             ind2=1
            else if (caseAO .eq. 2) then
             ind1=1
             ind2=iresp
            else if (caseAO .eq. 3) then
             ind1=1
             ind2=idir
            else if (caseAO .eq. 4) then
             ind1=idir
             ind2=1
            else
             call errquit
     &       ('get_C1MC: caseAO ne 1,2,3 or 4',
     &       0,GA_ERR)
            endif
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = ind1 ! pick direction iresp for g_dipEM
            ahi(3) = ind1
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            blo(3) = ind2 
            bhi(3) = ind2     
            clo(1) = 1
            chi(1) = nbf
            clo(2) = 1
            chi(2) = nocc
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,18) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18           format('FA-1::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
             endif
            endif

            call ga_zero(g_temp)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_dipEM  ,alo,ahi,
     &         g_vectors,blo,bhi,
     &         g_temp   ,clo,chi)
            if (debug) write (luout,*)
     &         'alfa: h(E) C(0) intermediate complete'
            alo(1) = 1
            ahi(1) = nocc
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = idir
            ahi(3) = idir
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            clo(1) = 1
            chi(1) = nocc
            clo(2) = 1
            chi(2) = nocc
            if (debug) then
             if (ga_nodeid().eq.0) then
             write(*,19) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 19          format('FA-2::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
             endif
            endif ! end-if-debug
            if (debug) then
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_vecE1 --------- START'
              call ga_print(g_vecE1)
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_vecE1 --------- END'
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_temp --------- START'
              call ga_print(g_temp)
             if (ga_nodeid().eq.0)
     &        write(*,*) '----g_temp --------- END'
            endif ! end-if-debug
c jbecca START
c        This should just be a sanity check
            sumre=0d0
c jbecca END
            sumre=trace(      ! out: trace of transp(A)* B
     &               g_vecE1, ! in : 
     &               g_temp,  ! in : 
     &               nocc,    ! in :  
     &               nbf,idir)! in : 
            sumim=0.0d0
            if (lifetime) then
            sumim=trace(        ! out: trace of transp(A)* B
     &               g_vecE1_im,! in : 
     &               g_temp,    ! in : 
     &               nocc,      ! in :  
     &               nbf,idir)  ! in : 
            endif ! end-if-lifetime
      return
      end

      double precision function 
     &                     trace( ! out: trace of transp(A)* B
     &                     g_A,   ! in : GA matrix A
     &                     g_B,   ! in : GA matrix B
     &                     nocc,  ! 
     &                     nbf,   !
     &                     idir)  !
c      Purpose: Calculate trace(transpose(A)*B) without the need
c               of doing a matrix multiplication A*B
c      Note1.- If we want  trace(A*B) swap: (1,2) for (2,1) in (aho,ahi)
c        alo(2) = 1
c        ahi(2) = nbf
c        alo(1) = i
c        ahi(1) = i
c      Note2.- In nga_ddot_patch() the op1,op2='n','t' do not work
c              or maybe it works if the resulting patch is not a vector
c              if the resulting patch was a matrix then it could work
c              transposing the resulting patch.
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
#include "geom.fh"
#include "stdio.fh"
       integer g_A,g_B
       double precision ac_trace
       integer idir,nbf,nocc,i
       integer alo(3),ahi(3),
     &         blo(3),bhi(3)
        trace=0.0d0
        alo(1) = 1
        ahi(1) = nbf
        alo(3) = idir
        ahi(3) = idir
        blo(1) = 1
        bhi(1) = nbf
        blo(3) = 1
        bhi(3) = 1
       do i=1,nocc
        alo(2) = i
        ahi(2) = i
        blo(2) = i
        bhi(2) = i
        ac_trace=nga_ddot_patch(g_A,'n',alo,ahi,
     &                          g_B,'n',blo,bhi) 
        trace=trace+ac_trace
       enddo ! end-loop-i
      return
      end
c $Id$
