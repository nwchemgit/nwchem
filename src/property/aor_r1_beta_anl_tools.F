      subroutine debug_update_g_anl_quad1(
     &            g_quadel,
     &            g_anl, 
     &            g_vecE1,
     &            g_vectors,
     &  ispin)
      
c Purpose: Debugging input to update_g_anl_quad
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
        
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
       integer ispin
       integer g_quadel,g_anl, 
     &         g_vecE1(2,2),g_vectors(2)
           if (ga_nodeid().eq.0)
     &     write(*,*) '---- g_quadel-------- START'
           call ga_print(g_quadel)
           if (ga_nodeid().eq.0)
     &     write(*,*) '---- g_quadel-------- END'

           if (ga_nodeid().eq.0) then
            write(*,17) ispin
 17         format('---- g_anl-lquad-bef(',i3,')-------- START')
           endif
           call ga_print(g_anl)
           if (ga_nodeid().eq.0) then
            write(*,18) ispin
 18         format('---- g_anl-lquad-bef(',i3,')-------- END')
           endif

           if (ga_nodeid().eq.0) then
            write(*,13) ispin
 13         format('---- g_vecE1(',i3,')-------- START')
           endif
           call ga_print(g_vecE1(ispin,1))
           if (ga_nodeid().eq.0) then
            write(*,14) ispin
 14         format('---- g_vecE1(',i3,')-------- END')
           endif

           if (ga_nodeid().eq.0) then
            write(*,15) ispin
 15         format('---- g_vectors(',i3,')-------- START')
           endif
           call ga_print(g_vectors(ispin))
           if (ga_nodeid().eq.0) then
            write(*,16) ispin
 16         format('---- g_vectors(',i3,')-------- END')
           endif
      return
      end

      subroutine update_g_anl_quad(
     &            g_anl,    ! in/out:
     &            g_quadel, ! in    : quadrupole AO integral
     &            g_vecE1,  ! in    : 1st-order elect. pert. MO vector
     &            g_vectors,! in    : MO vector
     &            scaling,  ! in    : scaling factor
     &            g_temp,   ! in    : scratch GA array
     &            g_work,   ! in    : scratch GA array
     &            idir,     ! in    : =1,2,3=x,y,z
     &            nocc,     ! in    : nr. occupied MOs      
     &            nbf,      ! in    : nr. basis functions  
     &            debug,    ! in    : debugging flag  
     &            lstatic)  ! in    : static    flag
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
        
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "prop.fh"
      integer npol
      integer g_anl,   ! in/out
     &        g_quadel,  
     &        g_vecE1,g_vectors,
     &        g_temp,g_work   
      logical lstatic,debug
      integer nbf,nocc,iresp,idir
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3), 
     &        clo(3),chi(3)      
      integer ind_dkl(2,3)
      data ind_dkl/ 2, 3,  ! dkl=123
     &              3, 1,  ! dkl=231
     &              1, 2 / ! dkl=312
      integer ind_p(2,2),p,ind1,ind2,k,l
      data ind_p/1, 2,   ! p12   p=1
     &           2, 1 /  ! p21   p=2
      integer LCred(2)
      data LCred/1,-1/
      integer qindex(3,3)
      double precision scl,scaling,one,two,four
      parameter(one=1.0d0,two=2.0d0,four=4.0d0)
      external get_C1MC
c     define translation table for quadrupole indices in
c     packed storage
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      qindex(1,1) = 1
      qindex(1,2) = 2
      qindex(2,1) = 2
      qindex(1,3) = 3
      qindex(3,1) = 3
      qindex(2,2) = 4
      qindex(2,3) = 5
      qindex(3,2) = 5
      qindex(3,3) = 6             
      do p=1,2 
        ind1=ind_p(1,p)
        ind2=ind_p(2,p)
c Allowing: (dkl)=(123,132,231,213,312,321)  d=idir
        k=ind_dkl(ind1,idir)  
        l=ind_dkl(ind2,idir)   
        scl=-LCred(p)*scaling               
        iresp = qindex(l,idir)  

        call get_C1MC(g_work,    ! out: C(E) M C
     &                g_vecE1,   ! in : 1st-order pert vec RE
     &                g_quadel,  ! in : dipole electric or magnetic
     &                g_vectors, ! in : MO vectors
     &                k,         ! in : = 1,2,3=x,y,z directions
     &                iresp,     ! in : = 1,2,3
     &                1,         ! in : indices in (alo,ahi)(3) (blo,bhi)(3)
     &                nbf,       ! in : nr. basis functions
     &                nocc,      ! in : nr. occupied alpha or beta
     &                debug,     ! in : logical var for debugging
     &                g_temp)    ! in : scratch GA array -> out

        clo(1) = 1
        chi(1) = nocc
        clo(2) = 1
        chi(2) = nocc           
        call nga_add_patch(scl,g_work,clo,chi,
     &                     one,g_anl ,clo,chi,
     &                         g_anl ,clo,chi)  
      enddo ! end-loop-p  
      return
      end

      subroutine update_g_anl_elmag(
     &              g_anl,  ! in/out:
     &              g_vecE1,! in    : C(E) pert MO vect.
     &              g_smatx,! in    : 0th/1st overlap deriv.
     &              g_vec,  ! in    : MO vect or C(B) MO vect
     &              scaling,! in    : scaling factor
     &              idir,   ! in    : direction (x,y or z)
     &              caseAO, ! in    : 3 or 4
     &              nbf,    ! in    : nr. basis functions
     &              nocc,   ! in    : nr. occupied MOs
     &              lstatic,! in    : flag for static calc.
     &              debug,  ! in    : flag for debugging
     &              g_work, ! in    : scratch GA array
     &              g_temp) ! in    : scratch GA array
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
        
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "bas.fh"
#include "stdio.fh"
       integer clo(3),chi(3)
       integer caseAO,nbf,nocc,idir
       logical lstatic,debug
       integer g_smatx,g_vecE1,g_vec,g_anl,
     &         g_temp,g_work ! scratch GA arrays
       double precision scaling,one,two,four
       parameter (one=1d0,two=2d0,four=4.0d0)
       external get_C1MC
       if (caseAO.ne.3 .and. caseAO.ne.4) then
        call errquit('update_g_anl_elmag: caseAO ne 3 or 4',
     &               1,INPUT_ERR)
       endif

         call get_C1MC(
     &           g_work,  ! out: C(E) M C
     &           g_vecE1, ! in : 1st-order pert vec RE
     &           g_smatx, ! in : dipole electric or magnetic
     &           g_vec,   ! in : MO vectors
     &           idir,    ! in : = 1,2,3=x,y,z directions
     &           0,       ! in : for caseAO=3,4 is dummy
     &           caseAO,  ! in : caseAO:(ind1,ind2)=(idir,1)
     &           nbf,     ! in : nr. basis functions
     &           nocc,    ! in : nr. occupied alpha or beta
     &           debug,   ! in : logical var for debugging
     &           g_temp)  ! in : scratch GA array   

       if (debug) then 
        if (ga_nodeid().eq.0) then
         write(*,5) 1
 5       format('---- g_work(',i3,')-------- START')
        endif
        call ga_print(g_work)
        if (ga_nodeid().eq.0) then
         write(*,6) 1
 6       format('---- g_work(',i3,')-------- END')
        endif
       endif ! end-if-debug
       clo(1) = 1
       chi(1) = nocc
       clo(2) = 1
       chi(2) = nocc
       call nga_add_patch(scaling,g_work,clo,chi,
     &                        one,g_anl ,clo,chi,
     &                            g_anl ,clo,chi)
      return
      end

      subroutine get_trace_ganl(g_anl,  ! in/out:
     &                          g_tran, ! in:
     &                          idir,   ! in: 1,2,3=x,y,z
     &                          oprint, ! in: logical var
     &                          nocc)   ! in: nr. occupied MOs
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

c Note (JA: by whom?) .- This routine is for testing only
c ====>  it will be better to use the idea of 
c        trace() defined in aor_r1_tensor.F
c        test: symmetrize the g_anl matrix before the LMO trafo:
      
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer nocc,i,idir,
     &        l_diag,k_diag
      integer g_tran, ! in
     &        g_anl,  ! in/out
     &        g_tmpanl! scratch GA array
      logical oprint
      character*(256) cstemp
      double precision sum
      double precision ga_trace_diag,trace,ac_trace
      external ga_trace_diag

      
      if (.not. ma_push_get(mt_dbl,nocc,'diag',l_diag,k_diag))
     &  call errquit('error alloc MA diag', 0, MA_ERR)
      
      
      if (oprint)
     &  write (luout,*) 'Message from beta_anl: Symmetrizing X'
      call ga_symmetrize(g_anl)
      
      write (cstemp,'(a)') 'aor_beta: tmpanl'
      if(.not.ga_create(MT_DBL,nocc,nocc,trim(cstemp),
     &  -1,-1,g_tmpanl))
     &  call errquit (trim(cstemp),0,GA_ERR)           
      call ga_zero(g_tmpanl)
      
      call ga_dgemm('t', 'n',nocc,nocc,nocc, 
     $  1.0d0,g_tran,g_anl,0.0d0,
     &  g_tmpanl)               ! out: g_tmpanl
      
      call ga_zero(g_anl)
      
c     ========= Added by FA to compute trace ====== START
      trace=0.0d0
      do i=1,nocc
        ac_trace=ga_ddot_patch(g_tmpanl,'n',i,i   ,1,nocc,
     &    g_tran  ,'n',1,nocc,i,i) 
        trace=trace+ac_trace
        
c$$$        if (ga_nodeid().eq.0) then
c$$$          write(*,2) i,trace
c$$$ 2        format('trace(',i3,')=',f15.8)
c$$$        endif
        
      enddo                     ! end-loop-i
c     ========= Added by FA to compute trace ====== END
      
      call ga_dgemm('n', 'n',nocc,nocc,nocc, 
     $  1.0d0,g_tmpanl,g_tran,0.0d0,
     &  g_anl)                  ! out: g_anl
      if (.not.ga_destroy(g_tmpanl)) call errquit
     &  ('aor_beta: ga_destroy failed g_tmpanl',0,GA_ERR)
      
      if (oprint) write (luout,
     &  '(/t12,a,t26,a/t11,6(''-''),t22,12(''-''))')
     &  'LMO #','contrib.'
      
      call ga_get_diagonal(g_anl,dbl_mb(k_diag))
      
      sum = 0.0d0
      do i=1,nocc
        
        sum = sum + dbl_mb(k_diag+i-1)
        
c$$$        if (ga_nodeid().eq.0) then
c$$$          write(*,1) i,sum
c$$$ 1        format('sum(',i3,')=',f15.8)
c$$$        endif
        
        if (oprint) write (luout,'(t11,i6,t22,f12.4)')
     &    i,dbl_mb(k_diag+i-1)
      enddo                     ! end-loop-i
      
c$$$      if (ga_nodeid().eq.0) then
c$$$        write(*,23) 1,sum,ga_trace_diag(g_anl),trace
c$$$ 23     format('g_anl-trace-lmo(',i3,')=(',
c$$$     &    f15.8,',',f15.8,',',f15.8,')')
c$$$      endif
      
      if (oprint)
     &  write (luout,'(1x,a,2i1,a,f12.4)')
     &  'Component ',idir,idir,': Sum = ',sum
      
      if (oprint) write (luout,'(1x,40(''-''))')
      
c$$$      sum = ga_trace_diag(g_anl)
c$$$      if (oprint) write (luout,*) 'sum from ga_trace: ',sum
      
      if (.not. ma_pop_stack(l_diag))
     &  call errquit('error deloc MA diag',0, MA_ERR)
      
      return
      end

      subroutine get_tracelessQtensor(g_quadel, ! in/out: quadrupole tensor
     &                                nbf,      ! in: nr. basis functions
     &                                g_work)   ! in: scratch GA array
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
        
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer g_quadel,g_work
      integer nbf,iresp
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3)
      integer qindex(3,3)
      double precision tenm8,one,two,three,
     &                 zero,half,third,four
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3d0,
     &           zero=0d0, half=one/two,
     &           third=one/three, four=4.0d0)
c     define translation table for quadrupole indices in
c     packed storage
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      qindex(1,1) = 1
      qindex(1,2) = 2
      qindex(2,1) = 2
      qindex(1,3) = 3
      qindex(3,1) = 3
      qindex(2,2) = 4
      qindex(2,3) = 5
      qindex(3,2) = 5
      qindex(3,3) = 6   
c     a) from trace -> g_work
      call ga_zero(g_work)  ! use for trace
      do iresp = 1,3          
        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = qindex(iresp,iresp) 
        ahi(3) = qindex(iresp,iresp)
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf           
        call nga_add_patch(one,g_quadel,alo,ahi,
     &                     one,g_work  ,blo,bhi,
     &                         g_work  ,blo,bhi)
       enddo ! end-loop-iresp         
c      b) scale quadel by 3
       call ga_scale(g_quadel,three)         
c      c) subtract trace from diagonal        
       do iresp = 1,3           
        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = qindex(iresp,iresp) 
        ahi(3) = qindex(iresp,iresp)
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf          
        call nga_add_patch(one,g_quadel,alo,ahi,
     &                    -one,g_work  ,blo,bhi,
     &                         g_quadel,alo,ahi)
       enddo ! end-loop-iresp         
c      d) divide the result by two, then by three
c         because of the factor by which the result enters
c         the Buckingham-Dunn tensor        
       call ga_scale(g_quadel,half*third)        
      return
      end

      subroutine write_vects(
     &              rtdb,geom,basis, ! in: handles
     &              g_vecE1,         ! in: C^(1,E)
     &              g_vecB1,         ! in: C^(1,B)
     &              g_tran,          ! in:
     &              g_vectmp,        ! in: scratch GA array
     &              npol,nocc,       ! in: nr. polariz,occ MOs
     &              nocct,nvirt,     ! in: nr. occ,virt MOs
     &              nopen,nmot,      ! in: nr. open shells,nmot=nocc*nvirt
     &              froct,epst,      ! in:
     &              nbf,             ! in: nr. basis functions
     &              lmo,             ! in: logical flag
     &              debug)           ! in: logical flag for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
        
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "prop.fh"
      integer rtdb,geom,basis
      integer g_vecE1,g_vecB1,g_tran,
     &        g_vectmp(2)
      integer idir,nbf,npol,nocc
      logical lmo,debug
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3),
     &        clo(3),chi(3)
      character*(256) cstemp
      integer nocct(npol),nvirt(npol),
     &        nopen(npol),nmot(npol)
      double precision froct(nbf,2),epst(nbf,2)
c        ---------------------------------------------------
c        part of the analysis is storing the perturbed
c        MOs. As elsewhere, we assume a closed shell system.
c        We write a lot of data here, e-field and b-field
c        perturbed canonical MOs and, if applicable, LMOs
c       ----------------------------------------------------
        do idir = 1,3
          alo(1) = 1
          ahi(1) = nbf
          alo(2) = 1
          ahi(2) = nbf
          alo(3) = idir         
          ahi(3) = idir
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nbf
          call ga_zero(g_vectmp(1))
          call nga_copy_patch('n',g_vecE1    ,alo,ahi,
     &                            g_vectmp(1),blo,bhi)
          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,24) idir
 24         format('---- g_vectmp-vecE1(',i3,')-------- START')
           endif
           call ga_print(g_vectmp(1))
           if (ga_nodeid().eq.0) then
            write(*,25) idir
 25         format('---- g_vectmp-vecE1(',i3,')-------- END')
           endif
          endif ! end-if-debug
          write(cstemp,'(a,i1,a)') 'cmo-efield',idir,'.movecs'
          call hnd_vec_write(rtdb,geom,basis,
     &                       nbf,nocct,nopen,nvirt,'scf',
     &                       g_vectmp,froct,epst,nmot,cstemp)
          call ga_zero(g_vectmp(1))
          call nga_copy_patch('n',g_vecB1    ,alo,ahi,
     &                            g_vectmp(1),blo,bhi)
          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,26) idir
 26         format('---- g_vectmp-vecB1(',i3,')-------- START')
           endif
           call ga_print(g_vectmp(1))
           if (ga_nodeid().eq.0) then
            write(*,27) idir
 27         format('---- g_vectmp-vecB1(',i3,')-------- END')
           endif
          endif ! end-if-debug
          write(cstemp,'(a,i1,a)') 'cmo-bfield',idir,'.movecs'
          call hnd_vec_write(rtdb,geom,basis,
     &                       nbf,nocct,nopen,nvirt,'scf',
     &                       g_vectmp,froct,epst,nmot,cstemp)
          if (lmo) then
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nocc
            alo(3) = idir         
            ahi(3) = idir
            blo(1) = 1
            bhi(1) = nocc
            blo(2) = 1
            bhi(2) = nocc
            clo(1) = 1
            chi(1) = nbf
            clo(2) = 1
            chi(2) = nocc
            call ga_zero(g_vectmp(1))
            call nga_matmul_patch('n','n',1d0,0d0,
     &                            g_vecE1    ,alo,ahi,
     &                            g_tran     ,blo,bhi,
     &                            g_vectmp(1),clo,chi)
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,28) idir
 28           format('---- g_vectmp-vecE1-lmo(',i3,')---- START')
             endif
             call ga_print(g_vectmp(1))
             if (ga_nodeid().eq.0) then
              write(*,29) idir
 29           format('---- g_vectmp-vecE1-lmo(',i3,')---- END')
             endif
            endif ! end-if-debug
            write(cstemp,'(a,i1,a)') 'lmo-efield',idir,'.movecs'
            call hnd_vec_write(rtdb,geom,basis,
     &                         nbf,nocct,nopen,nvirt,'scf',
     &                         g_vectmp,froct,epst,nmot,cstemp)
            call ga_zero(g_vectmp(1))
            call nga_matmul_patch('n','n',1d0,0d0,
     &                            g_vecB1    ,alo,ahi,
     &                            g_tran     ,blo,bhi,
     &                            g_vectmp(1),clo,chi)
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,30) idir
 30           format('---- g_vectmp-vecB1-lmo(',i3,')---- START')
             endif
             call ga_print(g_vectmp(1))
             if (ga_nodeid().eq.0) then
              write(*,31) idir
 31           format('---- g_vectmp-vecB1-lmo(',i3,')---- END')
             endif
            endif ! end-if-debug
            write(cstemp,'(a,i1,a)') 'lmo-bfield',idir,'.movecs'
            call hnd_vec_write(rtdb,geom,basis,
     &                         nbf,nocct,nopen,nvirt,'scf',
     &                         g_vectmp,froct,epst,nmot, cstemp)
          end if ! lmo
        enddo ! idir
      return
      end
c $Id$
