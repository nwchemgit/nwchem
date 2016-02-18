      subroutine rohf_hessv_2e2_opt(basis,  ! in : basis handle
     &                              geom,   ! in : geom  handle
     &                              nbf,    ! in : nr. basis functions
     &                              nmo,    ! in : nr. MOs vecs
     &                              nclosed,! in : nr. occupied MOs 
     &                              nopen,  ! in : nr. open shells (unpaired e's)
     $                              g_movec,! in : MO vec coeffs 
     &                              oskel,  ! in :
     &                              noskew, ! in : symm density ?
     &                              acc,    ! in : accuracy Fock construction   
     &                              g_x_re, ! in : RE g_x 
     &                              g_x_im, ! in : IM g_x
     &                              g_ax_re,! in/ou : RE g_ax Hessian product
     &                              g_ax_im,! in/ou : IM g_ax Hessian product
     &                              lifetime)
c
c   Purpose: Optimization of rohf_hessv_2e2()
c   Author : Fredy W. Aquino, Northwestern University
c   Date   : 03-28-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cscfps.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "case.fh"
c     
c     Return the ROHF orbital 2e-Hessian vector product, g_ax = A * g_x
      integer basis, geom       ! basis & geom handle
      integer nbf, nclosed, nopen ! Basis size and occupation
      integer nmo               ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel
      integer g_x_re ,g_x_im    ! Argument
      integer g_ax_re,g_ax_im   ! Hessian product
      double precision acc      ! Accuracy of "Fock" construction
    
      integer nvir, voff
      integer g_tmp,g_tmp1 ! scratch GA arrays

      integer ivec,nvec,gtype,vlen

      integer g_dens(2),g_fock(2) ! 2 stands for RE-IM blocks
      double precision tol2e
      logical odebug
      logical oprint
      logical noskew

      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one,zero,mone,four,
     &                 half,mhalf,two,mtwo
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      integer ncomp,cnt,nset,nmul,nblock,npol,ipm
      logical lifetime,debug
      data npol /1/ ! for unrestricted calculations
      external get_dens_reorim_1,shell_fock_build2

      if (nopen.ne.0) call errquit
     $     ('rohf_h2e2: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x_re,gtype,vlen,nvec)
c 
c      debug=.true. ! printouts allowed
      debug=.false. ! no printouts
c ======= set (nmul,nset,nblock) values ===== START
      ncomp=1 ! one single component
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=1 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=2 ! for RE-IM
      endif

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) lifetime,nset,nblock,npol,nmul,tol2e
 2001   format('(lifetime,nset,nblock,npol,nmul,tol2e)=(',
     &         L1,',',i3,',',i3,',',
     &         i3,',',i3,',',f15.8,')')
       endif
      endif ! end-if-debug
c ======= set (nmul,nset,nblock) values ===== END
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir  = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
c ======= create (g_dens,g_fock) ========== START
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      do ipm = 1,nblock ! =1 or 2  for ncomp=1
c ... allocate g_dens=[g_dens_re g_dens_im]
      dims(1)  = npol*nvec
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('uhf_hessv_2e2_opt: could not allocate g_dens',
     &                  555,GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
      dims(1)  = nmul*npol*nvec ! if npol=2 nmul=2 to store J+K integrals
                                ! this is done in shell_fock_build2()
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('uhf_hessv_2e2_opt: could not allocate g_fock',
     &                  555,GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm
c ======= create (g_dens,g_fock) ========== END
c ======= get dens for shell_fock_build ========= START
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug

           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    1,         ! in   : =1 1st block RE
     &                    g_x_re,    ! in   : 
     &                    g_movec,   ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nclosed,   ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    1,         ! in   : nr. polarizations
     &                    noskew,    ! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts

        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug

           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    2,         ! in   : =2 2nd block IM
     &                    g_x_im,    ! in   : IM
     &                    g_movec,   ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nclosed,   ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    1,         ! in   : nr. polarizations
     &                    noskew,    ! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts

        endif ! end-if-lifetime
c ======= get dens for shell_fock_build ========= END

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       noskew, ! in : =.false. for nonsymm dens
     &                       tol2e,  ! in :
     &                       debug)  ! in : = .true. -> debugging printouts

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-0(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-0(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0------ START' 
         call ga_print(g_ax_re)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0------ END'
       if (lifetime) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0------ START' 
         call ga_print(g_ax_im)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0------ END'
       endif ! end-if-lifetime
      endif ! end-if-debug
    
      if (oskel) 
     $   call errquit('rohf_h2e2: no way',0, UNKNOWN_ERR)

      g_tmp = ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp')
      g_tmp1= ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp1')
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nclosed  
      do cnt=1,nset ! loop in RE,IM
       do ivec = 1, nvec
         alo(1) = ivec
         ahi(1) = ivec

           if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,117) cnt,ivec
 117         format('XX:(cnt,ivec)=(',i3,',',i3,')')
            endif
           endif ! end-if-debug
c     
c     P      =  4(ij|kl) - (ik|jl) - (il|kj)
c     ij,kl
c     
c     K      =  (ik|jl) + (il|kj)
c     ij,kl
c     
c     cv         cv          pv   cp
c     Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c     pv          cv           cp   pv
c     Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c     
c     cp          cv           cp   pv
c     Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c     
c     Add the Fock matrices together overwriting the density
c     matrices to form the results above 
c     
c     Closed-Virtual bit
c     
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- START'
         call ga_print(g_fock(cnt)) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- END'
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- START'
         call ga_print(g_movec) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- END'
        endif ! end-if-debug

          call nga_matmul_patch('n','n',four,zero,
     &                          g_fock(cnt),alo,ahi,
     &                          g_movec    ,blo,bhi,
     &                          g_tmp      ,blo,bhi)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- START'
         call ga_print(g_tmp) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- END'
        endif ! end-if-debug

          call ga_zero(g_tmp1)
          call ga_matmul_patch('t','n',one,zero,
     $                         g_movec,voff,nmo ,1,nbf,     ! MO coefficients
     $                         g_tmp  ,1   ,nbf ,1,nclosed, ! result from step 1 
     $                         g_tmp1 ,1   ,nvir,1,nclosed) ! vir-occ Fock matrix

         if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,3701) cnt,ivec
3701     format('----- CvnFnnCno(',i3,',',i3,')------START')
          endif
           call ga_print(g_tmp1)
          if (ga_nodeid().eq.0) then
           write(*,3702) cnt,ivec
3702     format('----- CvnFnnCno(',i3,',',i3,')------END')
          endif
         endif ! end-if-debug

          if      (cnt.eq.1) then

          if (debug) then
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_ax-re-BEF-------- START'
           call ga_print(g_ax_re) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_ax-re-BEF-------- END'
          endif ! end-if-debug

c Note.- The operation below does:
c        g_ax_re= g_ax_re + 4 [4 C^T F C]  --> I am not sure if this is right.
          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
     $                       g_ax_re,1,ivec,four,'+')

          if (debug) then
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-re-AFT-------- START'
           call ga_print(g_ax_re) 
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-re-AFT-------- END'
          endif ! end-if-debug

          else if (cnt.eq.2) then

           if (debug) then
            if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-BEF-------- START'
            call ga_print(g_ax_im) 
            if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-BEF-------- END'
           endif ! end-if-debug

          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
     $                       g_ax_im,1,ivec,four,'+')

            if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-im-AFT-------- START'
             call ga_print(g_ax_im) 
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-im-AFT-------- END'
            endif ! end-if-debug

          endif ! end-if-cnt   

       enddo ! end-loop-ivec
      enddo ! end-loop-cnt

            if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-re-AFT-end---- START'
             call ga_print(g_ax_re) 
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-re-AFT-end---- END'
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-im-AFT-end---- START'
             call ga_print(g_ax_im) 
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_ax-im-AFT-end---- END'
            endif ! end-if-debug

        if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'FA-stop-test'
        endif ! end-if-debug

      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    
        if (.not.ga_destroy(g_tmp)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
        if (.not.ga_destroy(g_tmp1)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
      end

      subroutine rohf_hessv_2e2_opt_cmplx(
     &                          g_z,
     &                          g_Az1,
     &                          nsub,
     &                          basis,  ! in : basis handle
     &                          geom,   ! in : geom  handle
     &                          nbf,    ! in : nr. basis functions
     &                          nmo,    ! in : nr. MOs vecs
     &                          nclosed,! in : nr. occupied MOs 
     &                          nopen,  ! in : nr. open shells (unpaired e's)
     $                          g_movec,! in : MO vec coeffs 
     &                          oskel,  ! in :
     &                          noskew, ! in : symm density ?
     &                          acc,    ! in : accuracy Fock construction   
     &                          lifetime)
c
c   Purpose: Adaptation from rohf_hessv_2e2() for complex case.
c   Author : Fredy W. Aquino, Northwestern University
c   Date   : 03-28-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cscfps.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "case.fh"
c     
c     Return the ROHF orbital 2e-Hessian vector product, g_ax = A * g_x
      integer g_z,
     &        g_Az1, ! contains history of g_Az
     &        nsub,  ! points to block to be updated in g_Az
     &        g_xreim
      integer basis, geom       ! basis & geom handle
      integer nbf, nclosed, nopen ! Basis size and occupation
      integer nmo               ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel
      double precision acc      ! Accuracy of "Fock" construction
    
      integer nvir, voff
      integer g_tmp,g_tmp1 ! scratch GA arrays

      integer ivec,nvec,gtype,vlen

      integer g_dens(2),g_fock(2) ! 2 stands for RE-IM blocks
      double precision tol2e
      logical odebug
      logical oprint
      logical noskew

      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one,zero,mone,four,
     &                 half,mhalf,two,mtwo
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      integer ncomp,cnt,nset,nmul,nblock,npol,ipm
      logical lifetime,debug
      data npol /1/ ! for unrestricted calculations
      external get_dens_reorim_1,shell_fock_build2,
     &         update_gz_reorim1

      if (nopen.ne.0) call errquit
     $     ('rohf_h2e2: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
      call ga_inquire(g_z,gtype,vlen,nvec) ! out: nvec,vlen
       if (.not. ga_create(MT_DBL,vlen,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
c 
       debug=.false. ! no printouts
c ======= set (nmul,nset,nblock) values ===== START
      ncomp=1 ! one single component
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=1 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=2 ! for RE-IM
      endif

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) lifetime,nset,nblock,npol,nmul,tol2e
 2001   format('(lifetime,nset,nblock,npol,nmul,tol2e)=(',
     &         L1,',',i3,',',i3,',',
     &         i3,',',i3,',',f15.8,')')
       endif
      endif ! end-if-debug
c ======= set (nmul,nset,nblock) values ===== END
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir  = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
c ======= create (g_dens,g_fock) ========== START
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(2) = -1
      chunk(3) = -1
      do ipm = 1,nblock ! =1 or 2  for ncomp=1
c ... allocate g_dens=[g_dens_re g_dens_im]
      dims(1)  = npol*nvec
      chunk(1) = dims(1)
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('uhf_hessv_2e2_opt: could not allocate g_dens',
     &                  555,GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
      dims(1)  = nmul*npol*nvec ! if npol=2 nmul=2 to store J+K integrals
                                ! this is done in shell_fock_build2()
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('uhf_hessv_2e2_opt: could not allocate g_fock',
     &                  555,GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm
c ======= create (g_dens,g_fock) ========== END
c ======= get dens for shell_fock_build ========= START
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         call ga_zero(g_xreim)
         call getreorim(g_xreim,! out : real or im arr
     &                  g_z,    ! in  : = complx(g_xre,g_xim)
     &                  nvir,   ! in  : nr. virtual  MOs
     &                  nclosed,! in  : nr. occupied MOs
     &                  1)      ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END
           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    1,         ! in   : =1 1st block RE
     &                    g_xreim,   ! in   : 
     &                    g_movec,   ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nclosed,   ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    1,         ! in   : nr. polarizations
     &                    noskew,    ! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts
        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         call ga_zero(g_xreim)
         call getreorim(g_xreim,! out : real or im arr
     &                  g_z,    ! in  : = complx(g_xre,g_xim)
     &                  nvir,   ! in  : nr. virtual  MOs
     &                  nclosed,! in  : nr. occupied MOs
     &                  2)      ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END
           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    2,         ! in   : =2 2nd block IM
     &                    g_xreim,   ! in   : IM
     &                    g_movec,   ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nclosed,   ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    1,         ! in   : nr. polarizations
     &                    noskew,    ! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts
        endif ! end-if-lifetime
c ======= get dens for shell_fock_build ========= END
       if (.not.ga_destroy(g_xreim)) call errquit(
     &     'rohf_hessv2: ga_destroy failed g_xreim',0,GA_ERR)  

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       noskew, ! in : =.false. for nonsymm dens
     &                       tol2e,  ! in :
     &                       debug)  ! in : = .true. -> debugging printouts

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-0(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-0(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az1-0------ START' 
         call ga_print(g_Az1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az1-0------ END'
      endif ! end-if-debug
    
      if (oskel) 
     $   call errquit('rohf_h2e2: no way',0, UNKNOWN_ERR)

      g_tmp = ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp')
      g_tmp1= ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp1')
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nclosed  
      do cnt=1,nset ! loop in RE,IM
       do ivec = 1, nvec
         alo(1) = ivec
         ahi(1) = ivec

           if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,117) cnt,ivec
 117         format('XX:(cnt,ivec)=(',i3,',',i3,')')
            endif
           endif ! end-if-debug
c     
c     P      =  4(ij|kl) - (ik|jl) - (il|kj)
c     ij,kl
c     
c     K      =  (ik|jl) + (il|kj)
c     ij,kl
c     
c     cv         cv          pv   cp
c     Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c     pv          cv           cp   pv
c     Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c     
c     cp          cv           cp   pv
c     Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c     
c     Add the Fock matrices together overwriting the density
c     matrices to form the results above 
c     
c     Closed-Virtual bit
c     
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- START'
         call ga_print(g_fock(cnt)) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- END'
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- START'
         call ga_print(g_movec) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- END'
        endif ! end-if-debug

          call nga_matmul_patch('n','n',four,zero,
     &                          g_fock(cnt),alo,ahi,
     &                          g_movec    ,blo,bhi,
     &                          g_tmp      ,blo,bhi)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- START'
         call ga_print(g_tmp) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- END'
        endif ! end-if-debug

          call ga_zero(g_tmp1)
          call ga_matmul_patch('t','n',one,zero,
     $                         g_movec,voff,nmo ,1,nbf,     ! MO coefficients
     $                         g_tmp  ,1   ,nbf ,1,nclosed, ! result from step 1 
     $                         g_tmp1 ,1   ,nvir,1,nclosed) ! vir-occ Fock matrix

         if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,3701) cnt,ivec
3701     format('----- CvnFnnCno(',i3,',',i3,')------START')
          endif
           call ga_print(g_tmp1)
          if (ga_nodeid().eq.0) then
           write(*,3702) cnt,ivec
3702     format('----- CvnFnnCno(',i3,',',i3,')------END')
          endif
         endif ! end-if-debug

          if      (cnt.eq.1) then

          if (debug) then
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_Az1-re-BEF-------- START'
           call ga_print(g_Az1) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_Az1-re-BEF-------- END'
          endif ! end-if-debug

c Note.- The operation below does:
c        g_ax_re= g_ax_re + 4 [4 C^T F C]  --> I am not sure if this is right.
          call update_gz_reorim1(g_Az1,  ! out: = complx(g_xre,g_xim)
     &                           g_tmp1, ! in : real      arr 
     &                           1,      ! in  : =1 -> re =2 -> im
     &                           nsub,   ! in : index to sub-block in g_z
     &                           1,      ! in : = 1 or 2 index for component
     &                           vlen,   ! in : = nocc*nvir
     &                           four,   ! in  : scaling factor
     &                           nvir,
     &                           nclosed,
     &                           ivec)

          if (debug) then
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_Az1-re-AFT-------- START'
           call ga_print(g_Az1) 
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_Az1-re-AFT-------- END'
          endif ! end-if-debug

          else if (cnt.eq.2) then

           if (debug) then
            if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_Az1-im-BEF-------- START'
            call ga_print(g_Az1) 
            if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_Az1-im-BEF-------- END'
           endif ! end-if-debug

          call update_gz_reorim1(g_Az1,  ! out: = complx(g_xre,g_xim)
     &                           g_tmp1, ! in : real      arr 
     &                           2,      ! in  : =1 -> re =2 -> im
     &                           nsub,   ! in : index to sub-block in g_z
     &                           1,      ! in : = 1 or 2 index for component
     &                           vlen,   ! in : = nocc*nvir
     &                           four,   ! in  : scaling factor
     &                           nvir,
     &                           nclosed,
     &                           ivec)

            if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_Az1-im-AFT-------- START'
             call ga_print(g_Az1) 
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_Az1-im-AFT-------- END'
            endif ! end-if-debug

          endif ! end-if-cnt   

       enddo ! end-loop-ivec
      enddo ! end-loop-cnt

            if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_Az1-AFT-end---- START'
             call ga_print(g_Az1) 
             if (ga_nodeid().eq.0)
     &       write(*,*) '--------- g_Az1-AFT-end---- END'
            endif ! end-if-debug

        if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'FA-stop-test'
        endif ! end-if-debug

      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    
        if (.not.ga_destroy(g_tmp)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
        if (.not.ga_destroy(g_tmp1)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
      end
c $Id$
