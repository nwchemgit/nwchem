c ++++++++++++++++++ FA: uhf_hessv2_damp +++++++++++++++++++ START
      subroutine uhf_hessv3(acc,     ! in: accuracy of products
     &                      g_x,     ! in: A-matrix elements for density matrix REAL
     &                      g_ax,    ! in: Perturbed Fock operator REAL
     &                      g_x_im,  ! in: A-matrix elements for density matrix IMAG
     &                      g_ax_im, ! in: Perturbed Fock operator IMAG
     &                      omega,   ! in:
     &                      limag,   ! in: =.true. includes imaginary part
     &                      lifetime,! in: =.true. includes damping 
     &                      gamwidth,! in: 
     &                      ncomp)   ! in: nr. components
c
c Author: Fredy W. Aquino, Northwestern University (Oct 2012)
c Date  : 03-15-12
c
c Modified from uhf_hessv2 and following rohf_hessv3 (located in rohf_hessv3.F)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"

      logical limag,lifetime
      double precision omega,gamwidth
      integer ipm,ncomp

      double precision acc      ! [input] required accuracy of products
      integer g_x(ncomp),       ! [input] handle to input vectors
     &        g_x_im(ncomp)     ! [input] handle to input vectors
      integer g_ax(ncomp),      ! [input] handle to output products
     &        g_ax_im(ncomp)    ! [input] handle to output products

      integer gtype, vlen, nvec, nvecp, ivec
      double precision dnrm,wls,wlsim,coeffw
      double precision omg(ncomp),gam(ncomp)
      integer ilo(2), ihi(2)
      logical oprint, olprint
      external uhf_hessv_2e2_opt
      logical debug


      oprint = util_print("hessv",print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
c     Multiply a set of vectors by the level-shifted UHF hessian.
c     
c     Check dimensions    
      if(.not.cuhf_init_flag)
     $     call errquit('uhf_hessv2-dmp: UHF internal block invalid',
     &                  0,INPUT_ERR)

      do ipm=1,ncomp
       call ga_inquire(g_x(ipm),gtype,vlen,nvecp)
       if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gx-re',
     &                0, GA_ERR)
       call ga_inquire(g_ax(ipm),gtype,vlen,nvec)
       if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gax-re',
     &                0, GA_ERR)
       if (lifetime) then
        call ga_inquire(g_x_im(ipm),gtype,vlen,nvecp)
        if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gx-im',
     &                0, GA_ERR)
        call ga_inquire(g_ax_im(ipm),gtype,vlen,nvec)
        if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gax-im',
     &                0, GA_ERR)
       endif ! end-if-lifetime
      enddo ! end-loop-ipm

      if (nvecp .ne. nvec) 
     $     call errquit('uhf_hessv2-dmp: invalid nvecp',
     &                  nvecp, INPUT_ERR)

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x," g_x = ",i4,f24.8)')
     +      'uhf_hessv3',ivec,dnrm
          endif
        enddo
      endif

      do ipm = 1,ncomp
        call ga_zero(g_ax(ipm))
        if (lifetime) call ga_zero(g_ax_im(ipm))
      end do    

      coeffw=2.0d0 ! u-dft
      omg(1)=-omega
      omg(2)= omega
      gam(1)=-gamwidth
      gam(2)= gamwidth
      if (.not.lifetime) then
c       no damping: initialize Ax with terms proportional omega
         do ipm=1,ncomp
          wls   = lshift + coeffw * omg(ipm)
          call ga_dadd(wls  ,g_x(ipm),
     &                 0.0d0,g_ax(ipm),
     &                       g_ax(ipm))
         enddo ! end-lopp-ipm
      else                    ! lifetime
c       take care of damping here: Re and Im are coupled by gamwidth
         do ipm=1,ncomp
          wls   = lshift + coeffw * omg(ipm)
          wlsim = -coeffw * gam(ipm) 
          call ga_dadd(wls  ,g_x(ipm),
     &                 wlsim,g_x_im(ipm),
     &                       g_ax(ipm))
          wls   =  coeffw * gam(ipm) 
          wlsim = lshift + coeffw * omg(ipm)
          call ga_dadd(wls  ,g_x(ipm),
     &                 wlsim,g_x_im(ipm),
     &                       g_ax_im(ipm))  
         enddo ! end-lopp-ipm
      endif                   ! .not.lifetime
c ============== debug g_ax ==================== START
      debug=.false. ! .true.

      if (debug) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0-c(',ipm,')------ START' 
         call ga_print(g_ax(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0-c(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0-c(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0-c(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug
      call ga_sync()

c       next: add (e_a - e_i) times A (also called U) matrix to Ax
        do ipm=1,ncomp
         call uhf_hessv_1e(acc, g_x(ipm)   , g_ax(ipm)   , nvec) ! update g_ax 
        enddo ! end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call uhf_hessv_1e(acc, g_x_im(ipm), g_ax_im(ipm), nvec) ! update g_ax_im 
        enddo ! end-loop-ipm
        endif ! end-if-lifetime

      if (debug) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0-d(',ipm,')------ START' 
         call ga_print(g_ax(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0-d(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0-d(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0-d(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug
c ============== debug g_ax ==================== END

      if (pflg .gt. 1)then 
        if (ncomp.gt.1) then

            call uhf_hessv_2e3(acc, 
     &                         g_x, 
     &                         g_x_im,
     &                         g_ax,
     &                         g_ax_im, 
     &                         nvec,
     &                         limag,
     &                         lifetime)

        else                    ! call static 2e code

         debug=.false.
         if (debug) then
           do ipm=1,1
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_ax_re-0-x(',ipm,')------ START' 
            call ga_print(g_ax(ipm))
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_ax_re-0-x(',ipm,')------ END'
           enddo ! end-loop-ipm
          if (lifetime) then
           do ipm=1,1
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_ax_im-0-x(',ipm,')------ START' 
            call ga_print(g_ax_im(ipm))
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_ax_im-0-x(',ipm,')------ END'
           enddo ! end-loop-ipm
          endif ! end-if-lifetime
         endif ! end-if-debug     

            call uhf_hessv_2e2_opt(
     &                         acc, 
     &                         g_x(1), 
     &                         g_x_im(1),
     &                         g_ax(1),
     &                         g_ax_im(1), 
     &                         nvec,
     &                         lifetime)
        endif ! end-if-ncomp
      endif ! end-if-pflg

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_ax,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": g_ax = ",i4,f24.8)')
     +      'uhf_hessv3',ivec,dnrm
          endif
        enddo
      endif
      end

c ============= optimized alg ==================== START
      subroutine uhf_hessv_2e2_opt(
     &                         acc, 
     &                         g_x_re, 
     &                         g_x_im,
     &                         g_ax_re,
     &                         g_ax_im, 
     &                         nvec,
     &                         lifetime)
c
c Author: Fredy W. Aquino, Northwestern University (Oct 2012)
c Date  : 03-15-12
c         -> modifying uhf_hessv_2e2() for optimization
c Note.- Mimic uhf_hessv_2e2() but including two sets RE and IM.
c        (g_x_re,g_x_im) (g_ax_re,g_ax_im)

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
c     
      double precision acc      ! [input] required accuracy of products
      integer g_x_re ,g_x_im    ! [input] handle to input vectors
      integer g_ax_re,g_ax_im   ! [input] handle to output products  
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istart2(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,g_tmp1,g_tmp2
      integer g_dens(2),g_fock(2) ! 2 stands for RE-IM blocks
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim_1
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,ncomp,nblock,nmul,istartx,
     &        nocc_max,nvir_max
      logical lifetime,debug
      integer ip,im,ipm,cnt
      data npol /2/ ! for unrestricted calculations

      debug=.false. ! no printouts

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
      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'FA-noskew_uhf=',noskew_uhf 
      endif ! end-if-debug

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)

        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug

           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    1,         ! in   : =1 1st block RE
     &                    g_x_re,    ! in   : 
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    istartx,   ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
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
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    istartx,   ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts
        endif ! end-if-lifetime
      enddo ! end-loop-iset

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .false.,! in : =.false. for nonsymm dens
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

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 
   
c     start loop over components of perturbing field
      do cnt=1,nset 
       if      (cnt.eq.1) then
        g_ax=g_ax_re
       else if (cnt.eq.2) then
        g_ax=g_ax_im
       endif
       do iset = 1, npol
         if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)

        do ivec = 1, nvec

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,ipm
 117        format('XX:(cnt,iset,ivec,ipm)=(',
     &             i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug

           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_k        ,blo,bhi) 
c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-1 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_ax,     ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

              if (debug) then 
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-1 update_ax_fock'   
              endif ! end-if-debug
    
            endif ! end-if-same-spin-contrib
c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

             if (debug) then
              if (ga_nodeid().eq.0) then
               write(*,109) nocc_diff,nvir_diff 
 109           format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
              endif           
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-2 update_ax_fock'
             endif ! end-if-debug

              call update_ax_fock(
     &               g_ax,       ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug

           endif ! end-if-diff-spin-contrib
        enddo ! end-loop-ivec
        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
        enddo ! end-loop-iset (spin A,B)
          if      (cnt.eq.1) then
           g_ax_re=g_ax
          else if (cnt.eq.2) then
           g_ax_im=g_ax
          endif
       enddo ! end-loop-cnt

      if (debug) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-1------ START' 
         call ga_print(g_ax_re)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-1------ END'
       if (lifetime) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-1------ START' 
         call ga_print(g_ax_im)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-1------ END'
       endif ! end-if-lifetime
      endif ! end-if-debug

       if (.not. ga_destroy(g_j)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
       if (.not. ga_destroy(g_k)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
       do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'uhf_hv2e: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('uhf_hv2e: ga_destroy failed g_fock',0,GA_ERR)
       enddo ! end-loop-ipm    

      end

        subroutine get_dens_reorim_1(
     &                    g_dens,     ! out  : perturbed density matrix
     &                    cnt,        ! in/ou: counter of g_dens, =1 or 2
     &                    g_x,        ! in   : 
     &                    g_vecs,     ! in   : MO coefficients
     &                    nbf,        ! in   : nr. basis functions
     &                    nmo,        ! in   : nr. MOs
     &                    istart,     ! in   : shift nocc-nvirt block
     &                    nocc,       ! in   : nr. occupied MOs
     &                    nvir,       ! in   : nr. virtual  MOs 
     &                    nvec,       ! in   : nr. directions (x,y,z)
     &                    ipol,       ! in   : nr. polarizations
     &                    noskew_uhf, ! in : symmetry of density matrix
     &                    debug)      ! in   : =.true. -> debugging printouts  
c
c Author: Fredy W. Aquino, Northwestern University (Oct 2012)
c Date  : 03-15-12
c
c Note.- This routine works only for ncomp=1
c        It is being used when: g_x=g_x_re,g_x_im
c        g_dens acummulates the transformations C g_x C^t
c        g_dens(1) contains only transformations C g_x_re C^t
c        g_dens(2) contains only transformations C g_x_im C^t

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

      integer g_x               ! Argument: g_x_re or g_x_im
      integer g_dens(*)         ! size= 1 RE only or 2 RE+IM
      integer nbf,nmo,nocc,nvir
      integer g_vecs            ! MO coefficients
      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ivec,nvec,ipm,
     &        ipol, ! =1 for Alpha =2 for Beta
     &        shift,cnt
      character*(255) cstemp
      integer g_tmp1,g_tmp2,g_d
      integer istart,iend,vbase
      logical noskew_uhf,debug
      double precision one, zero, mone, 
     &                 four, half, mhalf, two, mtwo
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)

      if (.not. ga_create(MT_DBL,nbf,nbf,'gdens_reorim_1: g_d',
     $                    0, 0, g_d)) 
     &       call errquit('gdens_reorim_1: g_d', 0, GA_ERR)
      if (.not. ga_create(MT_DBL,nvir,nocc,'gdens_reorim_1: tmp1',
     $                    0, 0, g_tmp1)) 
     &       call errquit('gdens_reorim_1: tmp1', 0, GA_ERR)
      if (.not. ga_create(MT_DBL,nvir,nbf,'gdens_reorim_1: tmp2',
     $                    0, 0, g_tmp2)) 
     &       call errquit('gdens_reorim_1: tmp2', 0, GA_ERR)
      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      shift=(ipol-1)*nvec
      iend = istart + nocc*nvir - 1   
      vbase=nocc+1 ! First virtual    
      if (noskew_uhf) then
       do ivec=1,nvec
       call ga_zero(g_d)
       alo(1) = shift+ivec
       ahi(1) = alo(1)
       call ga_copy_patch('n', ! Reshape vector into matrix Xbj
     $                    g_x   ,istart,iend,ivec,ivec,
     $                    g_tmp1,1     ,nvir,1   ,nocc)
       call ga_matmul_patch('n', 't',one,zero,
     $                      g_tmp1,1,nvir,1,nocc,
     $                      g_vecs,1,nocc,1,nbf,
     $                      g_tmp2,1,nvir,1,nbf)
       call ga_matmul_patch('n', 'n', one, zero,
     $                      g_vecs,1,nbf ,vbase,nmo,
     $                      g_tmp2,1,nvir,1    ,nbf,
     $                      g_d   ,1,nbf ,1    ,nbf)
       call ga_symmetrize(g_d)
       call nga_copy_patch('N',g_d        ,blo,bhi,
     &                         g_dens(cnt),alo,ahi)
       enddo ! end-loop-ivecs
      else ! enter-if-noskew=.false.
       do ivec=1,nvec
       call ga_zero(g_d)
       alo(1) = shift+ivec
       ahi(1) = alo(1)
       call ga_copy_patch('n', ! Reshape vector into matrix Xbj
     $                    g_x   ,istart,iend,ivec,ivec,
     $                    g_tmp1,1     ,nvir,1   ,nocc)
       call ga_matmul_patch('n', 't',one,zero,
     $                      g_tmp1,1,nvir,1,nocc,
     $                      g_vecs,1,nocc,1,nbf,
     $                      g_tmp2,1,nvir,1,nbf)
       call ga_matmul_patch('n', 'n', one, zero,
     $                      g_vecs,1,nbf ,vbase,nmo,
     $                      g_tmp2,1,nvir,1    ,nbf,
     $                      g_d   ,1,nbf ,1    ,nbf)
       call ga_antisymmetrize(g_d)
       call nga_copy_patch('N',g_d        ,blo,bhi,
     &                         g_dens(cnt),alo,ahi)
       enddo ! end-loop-ivecs
      endif ! end-if-noskew_uhf
         if (.not.ga_destroy(g_tmp2)) call errquit(
     &         'gdens_reorim_1:ga_destroy failed g_tmp2',0,GA_ERR)
         if (.not.ga_destroy(g_tmp1)) call errquit(
     &         'gdens_reorim_1:ga_destroy failed g_tmp1',0,GA_ERR)
         if (.not.ga_destroy(g_d)) call errquit(
     &         'gdens_reorim_1:ga_destroy failed g_d',0,GA_ERR)
      return
      end

      subroutine uhf_hessv_2e3(acc, 
     &                         g_x_re, 
     &                         g_x_im,
     &                         g_ax_re,
     &                         g_ax_im, 
     &                         nvec,
     &                         limag,
     &                         lifetime)
c
c     Author : Fredy W. Aquino, Northwestern Unniversity (Oct 2012)
c     Date   : 03-24-12
c         -> modifying uhf_hessv_2e2()
c Note.- This routine works ONLY for ncomp=2
c        but should work for two cases:
c        lifetime=F  --> g_x_re,g_ax_re
c        lifetime=T  --> g_x_re,g_ax_re g_x_im,g_ax_im

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
c     
      logical limag
      double precision acc      ! [input] required accuracy of products
      integer g_x_re(2),        ! [input] handle to input vectors
     &        g_x_im(2)
      integer g_ax_re(2),       ! [input] handle to output products
     &        g_ax_im(2)
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istartx,istart2(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,g_tmp1,g_tmp2
      integer g_dens(4),g_fock(4)
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         CalcPerturbedTDPmat1_opt,
     &         get_undosymm_fock,update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,nblock,nmul
      logical lifetime,debug
      double precision coef(2,2)
      character*(255) cstemp
      integer ip,im,ipm,ncomp,nocc_max,nvir_max
      integer g_h1mat(2),           ! scratch GA array   
     &        g_pmats(2),g_pmata(2),! scratch GA array
     &        cnt,ind,indx(2,2)
      data npol /2/ ! for unrestricted calculations
      data indx /1,2, ! indx(1,1),indx(1,2)
     &           3,4/ ! indx(2,1),indx(2,2)
c     DIM/QM
      integer g_dens_tot(2), g_dim, clo(3), chi(3)
      integer g_dens_tot_i(2)
      logical ldimqm
c     DIM/QM END

      debug=.false. ! no printouts
      ncomp=2 ! using two components
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=2 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=4 ! for RE-IM
      endif

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) lifetime,nset,nblock,npol,nmul
 2001   format('(lifetime,nset,nblock,npol,nmul)=(',
     &         L1,',',i3,',',i3,',',i3,',',i3,')')
       endif
      endif ! end-if-debug

c ============= FA-defining scratch GA array g_h1mat() === START
      dims(1)  = nbf
      dims(2)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      do ipm = 1,ncomp
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmata(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmata(ipm))
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_h1mat(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo ! end-loop-ipm
c ============= FA-defining scratch GA array g_h1mat() === END
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1

      do ipm = 1,nblock ! =2 or 4 for ncomp=2
c ... allocate g_dens=[g_dens_re g_dens_im]
      dims(1)  = npol*nvec
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
      dims(1)  = nmul*npol*nvec ! if npol=2 nmul=2 to store J+K integrals
                                ! this is done in shell_fock_build2()
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm

      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'BEF get_dens_reorim-RE'
       endif ! end-if-debug

           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    1,      ! in   : =1 1st block RE
     &                    g_x_re, ! in   : 
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    istartx,! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array   

        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug

           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    2,      ! in   : =2 2nd block IM
     &                    g_x_im, ! in   : IM
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    istartx,! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array   

        endif ! end-if-lifetime
      enddo ! end-loop-iset

c jbecca START: adding in dim/qm functions here to handle UHF freq.
c              dependent case. 
c     DIM/QM jbecca
c     Calculate the DIM induced dipoles from the TOTAL density
c                                                -----
c     as written here, this has no extension to the damped (lifetime)
c     case. This is okay because this routine is only called for
c     FD no damping as of now
      if (.not.rtdb_get(bgj_get_rtdb_handle(), 'dimqm:lrsp', mt_log,
     $         1, ldimqm)) ldimqm = .false. 
      if (ldimqm) then 
         dims(1)     = nvec
         chunk(1)    = dims(1)
         do ipm = 1, ncomp
         if (.not. nga_create(MT_DBL, 3, dims, 'totdens', chunk,
     $         g_dens_tot(ipm))) call errquit('uhf_hessv3_2e3: DIM dens
     $                  create', 0, GA_ERR)

         alo(1)   = 1
         ahi(1)   = nvec

c     Copy alpha into total density array

         call nga_copy_patch('n', g_dens(ipm),    alo, ahi,
     $                        g_dens_tot(ipm),    alo, ahi)

         alo(1)   =  nvec+1
         ahi(1)   =  2*nvec
         clo(1)   =  1
         chi(1)   =  nvec
         clo(2)   =  1
         chi(2)   =  nbf
         clo(3)   =  1
         chi(3)   =  nbf

c     Add beta density to g_dens_tot
c     divided by 4 to match closed shell static routines. 
c     TODO: make sure that this is also being scaled in a similar fashion
c           in the closed shell FD cases

         call nga_add_patch(0.25d0, g_dens(ipm),     alo,  ahi,
     $                      0.25d0, g_dens_tot(ipm), clo,  chi,
     $                              g_dens_tot(ipm), clo,  chi)
         enddo       !ipm

c         call dimqm_addop_uhf(g_x_re, g_x_im, g_ax_re, g_ax_im, ncomp,
c     $                     limag, lifetime, g_dens_tot, g_dens_tot_i)
c
c         do ipm = 1, ncomp
c         if(.not. ga_destroy(g_dens_tot(ipm)))
c     $      call errquit('uhf_hessv3_2e3: DIM destroy dens', 0, GA_ERR)
c         enddo

c  Reset dimensions to previous values

         dims(1)  =  nmul*npol*nvec
         chunk(1) =  dims(1)
         alo(1)   =  0
         ahi(1)   =  0
      endif       ! ldimqm
      
c jbecca END

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.

      oskel_local = oskel .and. (nvec.eq.1)
      tol2e_local = min(max(acc,itol_floor),itol_ceil)

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .true., ! in : = .true for symm dens
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

      call get_undosymm_fock(
     &                g_fock,  ! in/ou: fock matrix
     &                nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &                nvec,    ! in   : nr. directions (x,y,z)
     &                nbf,     ! in   : nr. basis functions
     &                npol,    ! in   : nr. polarizations
     &                nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &                g_pmats, ! in   : scratch GA array
     &                limag)   ! in   : =.true. imaginary comp. exists

c ------- Remove GA arrays:
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_pmats(ipm))) call errquit(
     &     'uhf_hessv_2e3: ga_destroy failed g_pmats',0,GA_ERR)       
       if (.not.ga_destroy(g_pmata(ipm))) call errquit(
     &      'uhf_hessv_2e3: ga_destroy failed g_pmata',0,GA_ERR)
       if (.not.ga_destroy(g_h1mat(ipm)))  call errquit(
     &      'uhf_hessv_2e3: ga_destroy failed g_h1mat',0,GA_ERR)
      enddo ! end-loop-ipm

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0(',ipm,')------ START' 
         call ga_print(g_ax_re(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 
   
c     start loop over components of perturbing field
      do cnt=1,nset 
       do iset = 1, npol
        if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
        else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
        do ivec = 1, nvec
         do ipm=1,ncomp ! loop over Fock matrix components +/- here 
          ind=indx(ipm,cnt)

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,ipm,ind
 117        format('XX:(cnt,iset,ivec,ipm,ind)=(',
     &             i3,',',i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug

           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_k        ,blo,bhi) 

          if      (cnt.eq.1) then
           g_ax=g_ax_re(ipm)
          else if (cnt.eq.2) then
           g_ax=g_ax_im(ipm)
          endif
c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

            if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-1 update_ax_fock'
            endif ! end-if-debug

              call update_ax_fock(
     &               g_ax,     ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

             if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'AFT-1 update_ax_fock'      
              endif ! end-if-same-spin-contrib
             endif ! end-if-debug

c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

              if (debug) then
               if (ga_nodeid().eq.0) then
                write(*,109) nocc_diff,nvir_diff 
 109            format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
               endif           
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-2 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_ax,       ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug

           endif ! end-if-diff-spin-contrib

          if      (cnt.eq.1) then
           g_ax_re(ipm)=g_ax
          else if (cnt.eq.2) then
           g_ax_im(ipm)=g_ax
          endif

         enddo ! end-loop-ipm
        enddo ! end-loop-ivec

        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        enddo ! end-loop-iset (spin A,B)
       enddo ! end-loop-cnt
c  jbecca START
      if (ldimqm) then
         call dimqm_addop_uhf(g_x_re, g_x_im, g_ax_re, g_ax_im, ncomp,
     $                     limag, lifetime, g_dens_tot, g_dens_tot_i)

         do ipm = 1, ncomp
         if(.not. ga_destroy(g_dens_tot(ipm)))
     $      call errquit('uhf_hessv3_2e3: DIM destroy dens', 0, GA_ERR)
         enddo
      endif
c  jbecca END

      if (debug) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-1(',ipm,')------ START' 
         call ga_print(g_ax_re(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-1(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-1(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-1(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug

       if (.not. ga_destroy(g_j)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
       if (.not. ga_destroy(g_k)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
       do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
       enddo ! end-loop-ipm    

      end
c ============= optimized alg ==================== END

      subroutine get_undosymm_fock(
     &                  g_fock2, ! in/ou: fock matrix
     &                  nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &                  nvec,    ! in   : nr. directions (x,y,z)
     &                  nbf,     ! in   : nr. basis functions
     &                  npol,    ! in   : nr. polarizations
     &                  nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &                  g_p,     ! in   : scratch GA array
     &                  limag)   ! in   : =.true. -> imaginary component exist 
c           
c     Purpose: undo Fock symmetrization
c     Author : Fredy W. Aquino, Northwestern Unniversity (Oct 2012)
c     Date   : 03-24-12
c
c     If nset=1 => g_fock2(ipm) ipm=1,2
c     If nset=2 => g_fock2(ipm) ipm=1,2,3,4

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
      double precision pre_factor,
     &                 coeff(2,2)
      integer g_fock2(4), ! dim(g_fock2)= 2 or 4
     &        g_p(2)
      integer ivec,nvec,ipm,nbf,nset,npol,nmul,
     &        ind,ind1,ind2,cnt,indx(2,2),
     &        alo(3),ahi(3), 
     &        blo(3),bhi(3)  
      logical limag 
      data indx / 1,2, ! indx(1,1),indx(1,2)
     &            3,4/ ! indx(2,1),indx(2,2)
      
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
c ... jochen: next, we need to undo the symmetrization in order to
c     recover F(+) and F(-) separately. THIS HAS NOT YET BEEN ADAPTED
c     TO THE ROHF CASE!!!
c     apparently the density matrices are not needed any further
c ----- Construct coeffs for P(S),P(A) ------- START
      pre_factor = 1.0d0
      coeff(1,1)= pre_factor
      coeff(1,2)= pre_factor
      coeff(2,1)= pre_factor
      coeff(2,2)=-pre_factor
      if (limag) then
      coeff(1,1)= pre_factor
      coeff(1,2)= pre_factor
      coeff(2,1)=-pre_factor
      coeff(2,2)= pre_factor
      endif ! end-if-limag
c ----- Construct coeffs for P(S),P(A) ------- END
      do cnt=1,nset
       ind1=indx(1,cnt)
       ind2=indx(2,cnt)
       do ivec = 1,nvec*nmul*npol
        alo(1) = ivec
        ahi(1) = ivec 
c       use g_pmats for temp storage of the fock matrices
        call nga_copy_patch('N',g_fock2(ind1),alo,ahi,g_p(1),blo,bhi)
        call nga_copy_patch('N',g_fock2(ind2),alo,ahi,g_p(2),blo,bhi)
c
c ... jochen: the following should not make a diff. for PURE DFT
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        do ipm=1,2
          ind=indx(ipm,cnt)
          call nga_add_patch(coeff(ipm,1),g_p(1)      ,blo,bhi,
     &                       coeff(ipm,2),g_p(2)      ,blo,bhi,
     &                                    g_fock2(ind),alo,ahi)
        enddo ! end-loop-ipm
       enddo ! end-loop-ivec
      enddo ! end-loop-cnt
      return
      end

      subroutine update_ax_fock(
     &               g_ax,  ! in/ou: 
     &               g_fck, ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,! in   : MO vecs
     &               ivec,  ! in   : ivec-th MO
     &               vbase, ! in   : virtual base index for g_vecs
     &               nbf,   ! in   : nr. basis functions
     &               nmo,   ! in   : nr. MOs
     &               nocc,  ! in   : nr. occupied MOs
     &               nvir,  ! in   : nr. virtual  MOs
     &               istart,! in   : istart in g_ax(ipm)
     &               iend,  ! in   : iend   in g_ax(ipm)
     &               g_t1,  ! in   : scratch GA array
     &               g_t2,  ! in   : scratch GA array
     &               debug) ! in   : =.true. shows debugging printouts
c
c     Purpose: update g_ax with MO fock matrices
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
        
        double precision zero,one,four
        parameter (zero=0.0d0,one=1.0d0,four=4.0d0)
        integer g_ax,g_fck,g_vecs,
     &          g_t1,g_t2,
     &          nbf,nmo,nocc,nvir,vbase,
     &          ivec,istart,iend
        logical debug

        call ga_zero(g_t1)
        call ga_zero(g_t2)
        if (debug) then
         if (ga_nodeid().eq.0) then
         write(*,1) istart,iend,ivec,nvir,nocc
 1       format('In update_ax_fock:(istart,iend,ivec,nvir,nocc)=(',
     &           i3,',',i3,',',i3,',',i3,',',i3,')')
         endif
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- START'
         call ga_print(g_fck) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- END'
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- START'
         call ga_print(g_vecs) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- END'
        endif ! end-if-debug

        call ga_matmul_patch('n', 'n', one, zero,
     $                       g_fck ,1,nbf,1,nbf,
     $                       g_vecs,1,nbf,1,nocc,
     $                       g_t1  ,1,nbf,1,nocc)
        if (debug) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- FnnCno -------- START'
         call ga_print(g_t1) 
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- FnnCno -------- END'
        endif ! end-if-debug

        if (debug) then
         if (ga_nodeid().eq.0) then
           write(*,192) vbase,nmo,nbf,nocc,nvir
 192       format('update_ax_fock:(vbase,nmo,nbf,nocc,nvir)=(',
     &              i4,',',i4,',',i4,',',i4,',',i4,')')
         endif
        endif

        call ga_matmul_patch('t', 'n', one, zero,
     $                       g_vecs,vbase,nmo ,1,nbf,
     $                       g_t1  ,1    ,nbf ,1,nocc,
     $                       g_t2  ,1    ,nvir,1,nocc)

        if (debug) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- CvnFnnCno -------- START'
         call ga_print(g_t2) 
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- CvnFnnCno -------- END'
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- g_ax-BEF-------- START'
         call ga_print(g_ax) 
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- g_ax-BEF-------- END'
        endif ! end-if-debug
         call ga_dadd_patch(four,g_t2,1     ,nvir,1   ,nocc,
     $                      one ,g_ax,istart,iend,ivec,ivec,
     $                           g_ax,istart,iend,ivec,ivec)
        if (debug) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- g_ax-AFT-------- START'
         call ga_print(g_ax) 
         if (ga_nodeid().eq.0)
     &   write(*,*) '--------- g_ax-AFT-------- END'
        endif ! end-if-debug
      return
      end
c ++++++++++++++++++ FA: uhf_hessv2_damp +++++++++++++++++++ END
c 000000000000000000000000000000000000000000000000000000000000
c ++++++++++++++++++ udft calc using g_Az1 +++++++++++++ START
c 000000000000000000000000000000000000000000000000000000000000
c FA-04-25-12
      subroutine uhf_hessv3_cmplx(
     &                      acc,     ! in: accuracy of products
     &                      g_z,     ! in : z
     &                      g_Az1,
     &                      nsub,
     &                      omega,   ! in:
     &                      limag,   ! in: =.true. includes imaginary part
     &                      lifetime,! in: =.true. includes damping 
     &                      gamwidth,! in: 
     &                      ncomp)   ! in: nr. components
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c     Note.- Modified from uhf_hessv2 and following rohf_hessv3 (located in rohf_hessv3.F)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"

      logical limag,lifetime
      double precision omega,gamwidth
      integer ipm,ncomp

      double precision acc      ! [input] required accuracy of products

      integer g_x(ncomp),       ! [input] handle to input vectors
     &        g_x_im(ncomp)     ! [input] handle to input vectors
      integer g_ax(ncomp),      ! [input] handle to output products
     &        g_ax_im(ncomp)    ! [input] handle to output products

      integer g_z(ncomp),
     &        g_Az1,    ! history of g_Az (output)
     &        g_xreim,  ! scratch GA arrays
     &        g_Axreim, ! scratch GA arrays
     &        nsub,n,n1,maxsub,
     &        m1,m2,p1,p2,nreim,iset,npol,
     &        shift1(2),nocc1(2),
     &        nvir1(2),pretty
      double complex wls_cmplx
      double precision omg(ncomp),gam(ncomp)
      external conv2complex4,getreorim,getreorim1,
     &         getreorim1_u1,getreorim_u1,
     &         conv2complex4_u1

      integer gtype, vlen, nvec, nvecp, ivec
      double precision dnrm,wls,wlsim,coeffw
      integer ilo(2), ihi(2),alo(2),ahi(2)
      logical oprint, olprint
      external uhf_hessv_2e2_opt
      logical debug
      debug=.false.
      oprint = util_print("hessv",print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
c     Multiply a set of vectors by the level-shifted UHF hessian.
c     
c     Check dimensions    
      if(.not.cuhf_init_flag)
     $     call errquit('uhf_hessv2-dmp: UHF internal block invalid',
     &                  0,INPUT_ERR)

      goto 10 ! skip checking
      do ipm=1,ncomp
       call ga_inquire(g_x(ipm),gtype,vlen,nvecp)
       if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gx-re',
     &                0, GA_ERR)
       call ga_inquire(g_ax(ipm),gtype,vlen,nvec)
       if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gax-re',
     &                0, GA_ERR)
       if (lifetime) then
        call ga_inquire(g_x_im(ipm),gtype,vlen,nvecp)
        if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gx-im',
     &                0, GA_ERR)
        call ga_inquire(g_ax_im(ipm),gtype,vlen,nvec)
        if (vlen.ne.cuhf_vlen)
     $   call errquit('uhf_hessv2-dmp: invalid vector length-gax-im',
     &                0, GA_ERR)
       endif ! end-if-lifetime
      enddo ! end-loop-ipm

      if (nvecp .ne. nvec) 
     $     call errquit('uhf_hessv2-dmp: invalid nvecp',
     &                  nvecp, INPUT_ERR)

 10   continue
        call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
c ----- Clear sub-block ---- START
        call ga_inquire(g_Az1,gtype,n1,maxsub)
        alo(1)=1
        ahi(1)=n1
        alo(2)=nsub+1
        ahi(2)=nsub+nvec
        call nga_zero_patch(g_Az1,alo,ahi)
c ----- Clear sub-block ---- END

        coeffw=2.0d0 ! u-dft
        omg(1)=-omega
        omg(2)= omega
        gam(1)=-gamwidth
        gam(2)= gamwidth
      if (.not.lifetime) then
        write(*,*) 'uhf_hessv3_cmplx:STOP-not implemented'
        stop
c       no damping: initialize Ax with terms proportional omega
c        wls = lshift - coeffw * omega
c        call ga_dadd( wls, g_x(1), 0.d0, g_ax(1), g_ax(1) )
c        if (ncomp.gt.1) then
c        wls = lshift + coeffw * omega
c        call ga_dadd( wls, g_x(2), 0.d0, g_ax(2), g_ax(2) )
c        endif
      else                    ! lifetime
c       take care of damping here: Re and Im are coupled by gamwidth
c 00000000000000000000000000000000000 START
c        take care of damping here: Re and Im are coupled by gamwidth
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          wls   = lshift + coeffw * omg(ipm)
          wlsim = -coeffw * gam(ipm)
          wls_cmplx=dcmplx(wls,-wlsim)
          call ga_copy_patch('n',g_z(ipm),1 ,n ,1 ,nvec,
     &                           g_Az1   ,m1,m2,p1,p2)
          call ga_scale_patch(g_Az1,m1,m2,p1,p2,wls_cmplx)
          m1=m1+n
          m2=m2+n
         enddo ! end-lopp-ipm

c 00000000000000000000000000000000000 END
      endif                   ! .not.lifetime

      if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az1-c------ START' 
            call ga_print(g_Az1)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az1-c------ END'
      endif ! end-if-debug

c     next: add (e_a - e_i) times A (also called U) matrix to Ax
c 0000000000000000000000000000000000000000000000000000000000
c 0000000000 Adding 1e contrib to g_Az1 00000000000000 START
c 0000000000000000000000000000000000000000000000000000000000
       call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)

       shift1(1)=0
       shift1(2)=nalpha*(nmo-nalpha)*ncomp
       nocc1(1)=nalpha
       nocc1(2)=nbeta
       nvir1(1)=nmo-nocc1(1)
       nvir1(2)=nmo-nocc1(2)
       npol=2
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('uhf_hessv3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'uhf_hessv3_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('uhf_hessv3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)  
      do nreim=1,2 ! loop in RE,IM
        do ipm=1,ncomp
         call getreorim_u1(
     &                  g_xreim, ! out : real or im arr
     &                  g_z(ipm),! in  : = complx(g_xre,g_xim)
     &                  npol,    ! in  : nr. polarizations
     &                  nvir1,   ! in  : nr. virtual  MOs
     &                  nocc1,   ! in  : nr. occupied MOs
     &                  nreim)   ! in  : =1 -> re =2 -> im
         call getreorim1_u1(
     &                  g_Axreim,! out : real or im arr
     &                  g_Az1,   ! in  : = complx(g_xre,g_xim)
     &                  nsub,    ! in  : subblock index
     &                  ipm,     ! in  : = 1,2 to access slctd component
     &                  npol,    ! in  : nr. polarizations
     &                  nvir1,   ! in  : nr. virtual  MOs
     &                  nocc1,   ! in  : nr. occupied MOs
     &                  nreim)   ! in  : =1 -> re =2 -> im
         call uhf_hessv_1e(acc, 
     &                     g_xreim, 
     &                     g_Axreim, 
     &                     nvec) 
c ++++++++ update g_Az +++++++++ START
        call conv2complex4_u1(
     &                g_Az1,    ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                nsub,     ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                nreim)    ! in : =1 -> re =2 -> im
c ++++++++ update g_Az +++++++++ END
        enddo ! end-loop-ipm
      enddo ! end-loop-nreim
        if (.not. ga_destroy(g_xreim))  call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
        if (.not. ga_destroy(g_Axreim)) call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
c 0000000000000000000000000000000000000000000000000000000000
c 0000000000 Adding 1e contrib to g_Az1 00000000000000 END
c 0000000000000000000000000000000000000000000000000000000000
c ============== debug g_ax ==================== START
      if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az1-d------ START' 
            call ga_print(g_Az1)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az1-d------ END'
      endif ! end-if-debug
c ============== debug g_ax ==================== END

      if (pflg .gt. 1)then 
        if (ncomp.gt.1) then
            call uhf_hessv_2e3_cmplx(
     &                  acc, 
     &                  g_z,
     &                  g_Az1,   ! in: (n1,maxsub) history of Az matrix (large matrix)
     &                  nsub,    ! in: point to (n1,nvec) block to be updated in g_Az1
     &                  nvec,
     &                  limag,
     &                  lifetime)
        else                    ! call static 2e code

      if (debug) then
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az1-0-x------ START' 
         call ga_print(g_Az1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az1-0-x------ END'
      endif ! end-if-debug   
  
          if (ga_nodeid().eq.0) 
     &      write(*,*) 'FA-BEF-uhf_hessv_2e2-cmplx'

            call uhf_hessv_2e2_opt_cmplx(
     &                         acc, 
     &                         g_z(1), 
     &                         g_Az1,
     &                         nsub,
     &                         nvec,
     &                         lifetime)

         if (ga_nodeid().eq.0) 
     &      write(*,*) 'FA-AFT-uhf_hessv_2e2-cmplx'

        endif ! end-if-ncomp
       endif ! end-if-pflg
      end

      subroutine uhf_hessv_2e3_cmplx(
     &                  acc, 
     &                  g_z,
     &                  g_Az1,   ! in: (n1,maxsub) history of Az matrix (large matrix)
     &                  nsub,    ! in: point to (n1,nvec) block to be updated in g_Az1
     &                  nvec,
     &                  limag,
     &                  lifetime)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c Note.- This routine works ONLY for ncomp=2
c        but should work for two cases:
c        lifetime=F  --> g_x_re,g_ax_re
c        lifetime=T  --> g_x_re,g_ax_re g_x_im,g_ax_im

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
c     
      logical limag
      double precision acc      ! [input] required accuracy of products
      integer g_z(2),g_Az(2)
      integer g_Az1
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc1(2),nvir1(2),n,gtype,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istartx,istart2(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,g_Axreim,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,g_tmp1,g_tmp2
      integer g_dens(4),g_fock(4)
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         CalcPerturbedTDPmat1_opt,
     &         get_undosymm_fock,update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim,
     &         getreorim_u,conv2complex4_u
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,nblock,nmul
      logical lifetime,debug
      double precision coef(2,2)
      character*(255) cstemp
      integer ip,im,ipm,ncomp,nsub,
     &        shift,g_xreim(2)
      integer g_h1mat(2),           ! scratch GA array   
     &        g_pmats(2),g_pmata(2),! scratch GA array
     &        cnt,ind,indx(2,2),
     &        nocc_max,nvir_max
      data npol /2/ ! for unrestricted calculations
      data indx /1,2, ! indx(1,1),indx(1,2)
     &           3,4/ ! indx(2,1),indx(2,2)

      debug=.false. ! no printouts

      ncomp=2 ! using two components
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=2 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=4 ! for RE-IM
      endif

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) lifetime,nset,nblock,npol,nmul
 2001   format('(lifetime,nset,nblock,npol,nmul)=(',
     &         L1,',',i3,',',i3,',',i3,',',i3,')')
       endif
      endif ! end-if-debug

c ============= FA-defining scratch GA array g_h1mat() === START
      dims(1)  = nbf
      dims(2)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      do ipm = 1,ncomp
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_h1mat(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo ! end-loop-ipm
c ============= FA-defining scratch GA array g_h1mat() === END
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1

      do ipm = 1,nblock ! =2 or 4 for ncomp=2
c ... allocate g_dens=[g_dens_re g_dens_im]
      dims(1)  = npol*nvec
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
      dims(1)  = nmul*npol*nvec ! if npol=2 nmul=2 to store J+K integrals
                                ! this is done in shell_fock_build2()
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm

      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)
       shift=istartx-1

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'BEF get_dens_reorim-RE'
       endif ! end-if-debug

c ---- Copy g_z --> g_x_reim ------ START
        do ipm=1,ncomp
         if (.not. ga_create(MT_DBL,nocc*nvir,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim(ipm)))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
         call ga_zero(g_xreim(ipm))
         call getreorim_u(g_xreim(ipm),! out : real or im arr
     &                    g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                    shift,       ! in  : shift
     &                    nvir,        ! in  : nr. virtual  MOs
     &                    nocc,        ! in  : nr. occupied MOs
     &                    1)           ! in  : =1 -> re =2 -> im
        enddo ! end-loop-ipm
c ---- Copy g_z --> g_x_reim ------ END

           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    1,      ! in   : =1 1st block RE
     &                    g_xreim,! in   : 
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array  
 
        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug

c ---- Copy g_z --> g_x_reim ------ START

          do ipm=1,ncomp
           call ga_zero(g_xreim(ipm))
           call getreorim_u(g_xreim(ipm),! out : real or im arr
     &                    g_z(ipm),      ! in  : = complx(g_xre,g_xim)
     &                    shift,         ! in  : shift
     &                    nvir,          ! in  : nr. virtual  MOs
     &                    nocc,          ! in  : nr. occupied MOs
     &                    2)             ! in  : =1 -> re =2 -> im
          enddo ! end-loop-ipm

c ---- Copy g_z --> g_x_reim ------ END

           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    2,      ! in   : =2 2nd block IM
     &                    g_xreim,! in   : IM
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array - NOT USED  
     &                    g_h1mat)! in   : scratch GA array   

        endif ! end-if-lifetime
       do ipm = 1,ncomp
        if (.not.ga_destroy(g_xreim(ipm))) call errquit(
     &     'uhf_hessv3_cmplx: ga_destroy failed g_xreim',0,GA_ERR)  
       enddo ! end-loop-ipm
      enddo ! end-loop-iset
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_h1mat(ipm)))  call errquit(
     &      'uhf_hessv_2e3: ga_destroy failed g_h1mat',0,GA_ERR)
      enddo ! end-loop-ipm
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.

      oskel_local = oskel .and. (nvec.eq.1)
      tol2e_local = min(max(acc,itol_floor),itol_ceil)

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .true., ! in : = .true for symm dens
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

      call get_undosymm_fock(
     &                g_fock,  ! in/ou: fock matrix
     &                nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &                nvec,    ! in   : nr. directions (x,y,z)
     &                nbf,     ! in   : nr. basis functions
     &                npol,    ! in   : nr. polarizations
     &                nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &                g_pmats, ! in   : scratch GA array
     &                limag)   ! in   : =.true. imaginary comp. exists

c ------- Remove GA arrays:
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_pmats(ipm))) call errquit(
     &     'uhf_hessv_2e3: ga_destroy failed g_pmats',0,GA_ERR)  
      enddo ! end-loop-ipm

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
          if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az1-0------ START' 
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az1-0----- END'
      endif ! end-if-debug

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 
      call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
c     start loop over components of perturbing field
      nocc1(1)=nalpha
      nocc1(2)=nbeta
      nvir1(1)=nmo-nocc1(1)
      nvir1(2)=nmo-nocc1(2)
      do cnt=1,nset 
       do iset = 1, npol
         if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)

        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,194) nocc_max,nvir_max
 194      format('(noc_max,nvir_max)=(',i4,',',i4,')')
         endif
        endif ! end-if-debug

        if (debug) then
         if (ga_nodeid().eq.0) then
           write(*,192) cnt,iset,nocc,nocc_diff,
     &                  nvir,nvir_diff,istart,istart_diff,
     &                  vbase,vbase_diff
 192       format('(cnt,iset,nocc,nocc_diff,nvir,',
     &             'nvir_diff,istart,istart_diff,vbase,vbase_diff)=(',
     &              i4,',',i4,',',i4,',',i4,',',i4,',',i4,',',
     &              i4,',',i4,',',i4,',',i4,')')
         endif
        endif
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
        do ivec = 1, nvec
         do ipm=1,ncomp ! loop over Fock matrix components +/- here 
          ind=indx(ipm,cnt)

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,ipm,ind
 117        format('XX:(cnt,iset,ivec,ipm,ind)=(',
     &             i3,',',i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug
           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_k        ,blo,bhi) 

           call getreorim1_u1(
     &                      g_Axreim,! out : real or im arr
     &                      g_Az1,   ! in  : = complx(g_xre,g_xim)
     &                      nsub,    ! in  : subblock index
     &                      ipm,     ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      cnt)     ! in  : =1 -> re =2 -> im

c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

            if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-1 update_ax_fock'
            endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim, ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

             if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'AFT-1 update_ax_fock'      
              endif ! end-if-same-spin-contrib
             endif ! end-if-debug

c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

              if (debug) then
               if (ga_nodeid().eq.0) then
                write(*,109) nocc_diff,nvir_diff 
 109            format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
               endif           
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-2 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim,   ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug

           endif ! end-if-diff-spin-contrib
c ---- Move g_Axreim --> g_Az1---- START

               call conv2complex4_u1(
     &                g_Az1,    ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                nsub,     ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                cnt)      ! in : =1 -> re =2 -> im

c ---- Move g_Axreim --> g_Az1---- END
         enddo ! end-loop-ipm
        enddo ! end-loop-ivec

        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)

        enddo ! end-loop-iset (spin A,B)
       enddo ! end-loop-cnt
        if (.not. ga_destroy(g_Axreim)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)

      if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-1------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-1------ END'
      endif ! end-if-debug

      if (.not. ga_destroy(g_j)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_k)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    

      end

      subroutine uhf_hessv_2e2_opt_cmplx(
     &                         acc, 
     &                         g_z,
     &                         g_Az1,
     &                         nsub,
     &                         nvec,
     &                         lifetime)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c Note.- Mimic uhf_hessv_2e2() but including two sets RE and IM.
c        (g_x_re,g_x_im) (g_ax_re,g_ax_im)

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
      integer g_z,    ! in: 
     &        g_Az1,  ! in: history of g_Az 
     &        g_xreim,! Scratch GA array  
     &        g_Axreim! Scratch GA array
      integer nsub,m1,m2,p1,p2
      double precision acc      ! [input] required accuracy of products
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istart2(2),
     &        nocc1(2),nvir1(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,
     &        g_tmp1,g_tmp2,n,gtype
      integer g_dens(2),g_fock(2) ! 2 stands for RE-IM blocks
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim_1,getreorim_u,
     &         getreorim1_u1,conv2complex4_u1
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,ncomp,nblock,nmul,istartx
      logical lifetime,debug
      integer ip,im,ipm,cnt,nocc_max,nvir_max
      data npol /2/ ! for unrestricted calculations

      debug=.false. ! no printouts

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
      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'FA-noskew_uhf=',noskew_uhf 
      endif ! end-if-debug

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         if (.not. ga_create(MT_DBL,nocc*nvir,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
         call ga_zero(g_xreim)
         call getreorim_u(g_xreim,  ! out : real or im arr
     &                    g_z,      ! in  : = complx(g_xre,g_xim)
     &                    istartx-1,! in  : shift
     &                    nvir,     ! in  : nr. virtual  MOs
     &                    nocc,     ! in  : nr. occupied MOs
     &                    1)        ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END

           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    1,         ! in   : =1 1st block RE
     &                    g_xreim,   ! in   : 
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts

        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         call ga_zero(g_xreim)
         call getreorim_u(g_xreim,  ! out : real or im arr
     &                    g_z,      ! in  : = complx(g_xre,g_xim)
     &                    istartx-1,! in  : shift
     &                    nvir,     ! in  : nr. virtual  MOs
     &                    nocc,     ! in  : nr. occupied MOs
     &                    2)        ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END

           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    2,         ! in   : =2 2nd block IM
     &                    g_xreim,   ! in   : IM
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts

        endif ! end-if-lifetime
       if (.not.ga_destroy(g_xreim)) call errquit(
     &     'uhf_hessv2_cmplx: ga_destroy failed g_xreim',0,GA_ERR) 
      enddo ! end-loop-iset

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .false.,! in : =.false. for nonsymm dens
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

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 


      call ga_inquire(g_z,gtype,n,nvec) ! get (n,nvec)
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx2_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
      nocc1(1)=nalpha
      nocc1(2)=nbeta
      nvir1(1)=nmo-nocc1(1)
      nvir1(2)=nmo-nocc1(2)   
c     start loop over components of perturbing field
      do cnt=1,nset 
       do iset = 1, npol
         if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)

        do ivec = 1, nvec

           call getreorim1_u1(
     &                      g_Axreim,! out : real or im arr
     &                      g_Az1,   ! in  : = complx(g_xre,g_xim)
     &                      nsub,    ! in  : subblock index
     &                      1,       ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      cnt)     ! in  : =1 -> re =2 -> im

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,1
 117        format('XX:(cnt,iset,ivec,ipm)=(',
     &             i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug

           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_k        ,blo,bhi) 
c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-1 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim, ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

              if (debug) then 
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-1 update_ax_fock'   
              endif ! end-if-debug
    
            endif ! end-if-same-spin-contrib
c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

             if (debug) then
              if (ga_nodeid().eq.0) then
               write(*,109) nocc_diff,nvir_diff 
 109           format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
              endif           
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-2 update_ax_fock'
             endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim,   ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug

           endif ! end-if-diff-spin-contrib
c ---- Move g_Axreim --> g_Az1---- START
               call conv2complex4_u1(
     &                g_Az1,    ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                nsub,     ! in : subblock index
     &                1,        ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                cnt)      ! in : =1 -> re =2 -> im
c ---- Move g_Axreim --> g_Az1---- END
        enddo ! end-loop-ivec
        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hv2e: tmp1',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hv2e: tmp2',0,
     &       GA_ERR)
        enddo ! end-loop-iset (spin A,B)
       enddo ! end-loop-cnt
        if (.not. ga_destroy(g_Axreim)) call errquit('uhf_hv2e: Axreim',0,
     &       GA_ERR)
      if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-1------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-1------ END'
      endif ! end-if-debug

      if (.not. ga_destroy(g_j)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_k)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'uhf_hv2e: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('uhf_hv2e: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-stop-quest'
        stop
       endif ! end-if-debug

      end
c 000000000000000000000000000000000000000000000000000000000000
c ++++++++++++++++++ udft calc using g_Az1 +++++++++++++ END
c 000000000000000000000000000000000000000000000000000000000000
      subroutine uhf_hessv3_cmplx1(
     &                      acc,     ! in: accuracy of products
     &                      g_z,     ! in : z
     &                      g_Az,
     &                      omega,   ! in:
     &                      limag,   ! in: =.true. includes imaginary part
     &                      lifetime,! in: =.true. includes damping 
     &                      gamwidth,! in: 
     &                      ncomp)   ! in: nr. components
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c Modified from uhf_hessv2 and following rohf_hessv3 (located in rohf_hessv3.F)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"

      logical limag,lifetime
      double precision omega,gamwidth
      integer ipm,ncomp

      double precision acc      ! [input] required accuracy of products

      integer g_x(ncomp),       ! [input] handle to input vectors
     &        g_x_im(ncomp)     ! [input] handle to input vectors
      integer g_ax(ncomp),      ! [input] handle to output products
     &        g_ax_im(ncomp)    ! [input] handle to output products

      integer g_z(ncomp),
     &        g_Az,  
     &        g_xreim,  ! scratch GA arrays
     &        g_Axreim, ! scratch GA arrays
     &        n,n1,
     &        m1,m2,nreim,iset,npol,
     &        shift1(2),nocc1(2),
     &        nvir1(2),pretty
      double complex wls_cmplx
      double precision omg(ncomp),gam(ncomp)
      external conv2complex4,getreorim,getreorim1,
     &         getreorim1_u1,getreorim_u1,
     &         conv2complex4_u1,
     &         uhf_hessv_2e3_cmplx1,
     &         uhf_hessv_2e2_opt_cmplx1

      integer gtype, vlen, nvec, nvecp, ivec
      double precision dnrm,wls,wlsim,coeffw
      integer ilo(2), ihi(2),alo(2),ahi(2)
      logical oprint, olprint
      logical debug
      debug=.false.
      oprint = util_print("hessv",print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
c     Multiply a set of vectors by the level-shifted UHF hessian.
c     
c     Check dimensions    
      if(.not.cuhf_init_flag)
     $     call errquit('uhf_hessv2-dmp: UHF internal block invalid',
     &                  0,INPUT_ERR)

        call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
        call ga_zero(g_Az)
        coeffw=2.0d0 ! u-dft
        omg(1)=-omega
        omg(2)= omega
        gam(1)=-gamwidth
        gam(2)= gamwidth
      if (.not.lifetime) then
        write(*,*) 'uhf_hessv3_cmplx:STOP-not implemented'
        stop
c       no damping: initialize Ax with terms proportional omega
c        wls = lshift - coeffw * omega
c        call ga_dadd( wls, g_x(1), 0.d0, g_ax(1), g_ax(1) )
c        if (ncomp.gt.1) then
c        wls = lshift + coeffw * omega
c        call ga_dadd( wls, g_x(2), 0.d0, g_ax(2), g_ax(2) )
c        endif
      else                    ! lifetime
c        take care of damping here: Re and Im are coupled by gamwidth
         m1=1
         m2=n
         do ipm=1,ncomp
          wls   = lshift + coeffw * omg(ipm)
          wlsim = -coeffw * gam(ipm)
          wls_cmplx=dcmplx(wls,-wlsim)
          call ga_copy_patch('n',g_z(ipm),1 ,n ,1,nvec,
     &                           g_Az    ,m1,m2,1,nvec)
          call ga_scale_patch(g_Az,m1,m2,1,nvec,wls_cmplx)
          m1=m1+n
          m2=m2+n
         enddo ! end-lopp-ipm
      endif                   ! .not.lifetime

      if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az-c------ START' 
            call ga_print(g_Az)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az-c------ END'
      endif ! end-if-debug

c     next: add (e_a - e_i) times A (also called U) matrix to Ax
c 0000000000000000000000000000000000000000000000000000000000
c 0000000000 Adding 1e contrib to g_Az 00000000000000 START
c 0000000000000000000000000000000000000000000000000000000000
       call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
       shift1(1)=0
       shift1(2)=nalpha*(nmo-nalpha)*ncomp
       nocc1(1)=nalpha
       nocc1(2)=nbeta
       nvir1(1)=nmo-nocc1(1)
       nvir1(2)=nmo-nocc1(2)
       npol=2
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('uhf_hessv3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'uhf_hessv3_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('uhf_hessv3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR) 
      do nreim=1,2 ! loop in RE,IM
        do ipm=1,ncomp
         call getreorim_u1(
     &                  g_xreim, ! out : real or im arr
     &                  g_z(ipm),! in  : = complx(g_xre,g_xim)
     &                  npol,    ! in  : nr. polarizations
     &                  nvir1,   ! in  : nr. virtual  MOs
     &                  nocc1,   ! in  : nr. occupied MOs
     &                  nreim)   ! in  : =1 -> re =2 -> im
         call getreorim1_u1(
     &                  g_Axreim,! out : real or im arr
     &                  g_Az,    ! in  : = complx(g_xre,g_xim)
     &                  0,       ! in  : subblock index
     &                  ipm,     ! in  : = 1,2 to access slctd component
     &                  npol,    ! in  : nr. polarizations
     &                  nvir1,   ! in  : nr. virtual  MOs
     &                  nocc1,   ! in  : nr. occupied MOs
     &                  nreim)   ! in  : =1 -> re =2 -> im
         call uhf_hessv_1e(acc, 
     &                     g_xreim, 
     &                     g_Axreim, 
     &                     nvec) 

        call conv2complex4_u1(  ! update g_Az
     &                g_Az,     ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                0,        ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                nreim)    ! in : =1 -> re =2 -> im

        enddo ! end-loop-ipm
      enddo ! end-loop-nreim
        if (.not. ga_destroy(g_xreim))  call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
        if (.not. ga_destroy(g_Axreim)) call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
c 0000000000000000000000000000000000000000000000000000000000
c 0000000000 Adding 1e contrib to g_Az 00000000000000 END
c 0000000000000000000000000000000000000000000000000000000000
c ============== debug g_ax ==================== START
      if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az-d------ START' 
            call ga_print(g_Az)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az-d------ END'
      endif ! end-if-debug
c ============== debug g_ax ==================== END

      if (pflg .gt. 1)then 
        if (ncomp.gt.1) then
            call uhf_hessv_2e3_cmplx1(
     &                  acc, 
     &                  g_z,
     &                  g_Az,    ! in: (n1,maxsub) history of Az matrix (large matrix)
     &                  nvec,
     &                  limag,
     &                  lifetime)
        else                    ! call static 2e code

         if (debug) then
          if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az-0-x------ START' 
          call ga_print(g_Az)
          if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az-0-x------ END'
         endif ! end-if-debug     

            call uhf_hessv_2e2_opt_cmplx1(
     &                         acc, 
     &                         g_z(1), 
     &                         g_Az,
     &                         nvec,
     &                         lifetime)
        endif ! end-if-ncomp
      endif ! end-if-pflg
      end

      subroutine uhf_hessv_2e3_cmplx1(
     &                  acc, 
     &                  g_z,
     &                  g_Az,   ! in: (n1,maxsub) history of Az matrix (large matrix)
     &                  nvec,
     &                  limag,
     &                  lifetime)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c Note.- This routine works ONLY for ncomp=2
c        but should work for two cases:
c        lifetime=F  --> g_x_re,g_ax_re
c        lifetime=T  --> g_x_re,g_ax_re g_x_im,g_ax_im

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
c     
      logical limag
      double precision acc      ! [input] required accuracy of products
      integer g_z(2)
      integer g_Az
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc1(2),nvir1(2),n,gtype,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istartx,istart2(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,g_Axreim,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,g_tmp1,g_tmp2
      integer g_dens(4),g_fock(4)
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         CalcPerturbedTDPmat1_opt,
     &         get_undosymm_fock,update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim,
     &         getreorim_u,conv2complex4_u
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,nblock,nmul
      logical lifetime,debug
      double precision coef(2,2)
      character*(255) cstemp
      integer ip,im,ipm,ncomp,
     &        shift,g_xreim(2)
      integer g_h1mat(2),           ! scratch GA array   
     &        g_pmats(2),g_pmata(2),! scratch GA array
     &        cnt,ind,indx(2,2),
     &        nocc_max,nvir_max
      data npol /2/ ! for unrestricted calculations
      data indx /1,2, ! indx(1,1),indx(1,2)
     &           3,4/ ! indx(2,1),indx(2,2)
c  jbecca START
      logical ldimqm
      integer g_dens_tot(2), g_dim, clo(3), chi(3)
      integer g_dens_tot_i(2)
      integer g_ax_re(2), g_ax_im(2)
c  jbecca END

      debug=.false. ! no printouts

      ncomp=2 ! using two components
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=2 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=4 ! for RE-IM
      endif

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) lifetime,nset,nblock,npol,nmul
 2001   format('(lifetime,nset,nblock,npol,nmul)=(',
     &         L1,',',i3,',',i3,',',i3,',',i3,')')
       endif
      endif ! end-if-debug

c ============= FA-defining scratch GA array g_h1mat() === START
      dims(1)  = nbf
      dims(2)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      do ipm = 1,ncomp
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_h1mat(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo ! end-loop-ipm
c ============= FA-defining scratch GA array g_h1mat() === END
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1

      do ipm = 1,nblock ! =2 or 4 for ncomp=2
c ... allocate g_dens=[g_dens_re g_dens_im]
      dims(1)  = npol*nvec
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
      dims(1)  = nmul*npol*nvec ! if npol=2 nmul=2 to store J+K integrals
                                ! this is done in shell_fock_build2()
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm

      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)
       shift=istartx-1

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'BEF get_dens_reorim-RE'
       endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
        do ipm=1,ncomp
         if (.not. ga_create(MT_DBL,nocc*nvir,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim(ipm)))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
         call ga_zero(g_xreim(ipm))
         call getreorim_u(g_xreim(ipm),! out : real or im arr
     &                    g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                    shift,       ! in  : shift
     &                    nvir,        ! in  : nr. virtual  MOs
     &                    nocc,        ! in  : nr. occupied MOs
     &                    1)           ! in  : =1 -> re =2 -> im

        enddo ! end-loop-ipm
c ---- Copy g_z --> g_x_reim ------ END

           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    1,      ! in   : =1 1st block RE
     &                    g_xreim,! in   : 
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array   

        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
        do ipm=1,ncomp
         call ga_zero(g_xreim(ipm))
         call getreorim_u(g_xreim(ipm),! out : real or im arr
     &                    g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                    shift,       ! in  : shift
     &                    nvir,        ! in  : nr. virtual  MOs
     &                    nocc,     ! in  : nr. occupied MOs
     &                    2)           ! in  : =1 -> re =2 -> im
        enddo ! end-loop-ipm

c ---- Copy g_z --> g_x_reim ------ END
           call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    2,      ! in   : =2 2nd block IM
     &                    g_xreim,! in   : IM
     &                    g_vecs, ! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    iset,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array - NOT USED  
     &                    g_h1mat)! in   : scratch GA array   

        endif ! end-if-lifetime
       do ipm = 1,ncomp
        if (.not.ga_destroy(g_xreim(ipm))) call errquit(
     &     'uhf_hessv3_cmplx: ga_destroy failed g_xreim',0,GA_ERR)  
       enddo ! end-loop-ipm
      enddo ! end-loop-iset
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_h1mat(ipm)))  call errquit(
     &      'uhf_hessv_2e3: ga_destroy failed g_h1mat',0,GA_ERR)
      enddo ! end-loop-ipm
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.

c  DIM/QM jbecca START: Allocate arrays and create real and imaginary
c                    densities for total densities. 
      if (.not.rtdb_get(bgj_get_rtdb_handle(), 'dimqm:lrsp', mt_log,
     $         1, ldimqm)) ldimqm = .false. 
      if (ldimqm) then 
         dims(1)     = nvec
         chunk(1)    = dims(1)
         do ipm = 1, ncomp
         if (.not. nga_create(MT_DBL, 3, dims, 'totdens', chunk,
     $         g_dens_tot(ipm))) call errquit('uhf_hessv3_2e3: DIM dens
     $                  create', 0, GA_ERR)
         if (lifetime) then
            if (.not. nga_create(MT_DBL, 3, dims, 'totdens', chunk,
     $         g_dens_tot_i(ipm))) call errquit('uhf_hessv3_2e3: DIM 
     $             dens create', 0, GA_ERR)
         endif       !lifetime

         alo(1)   =  1
         ahi(1)   =  nvec
c     Copy alpha into total density array

         call nga_copy_patch('n', g_dens(ipm),    alo, ahi,
     $                        g_dens_tot(ipm),    alo, ahi)

         if (lifetime) then
            call nga_copy_patch('n', g_dens(ipm+2),    alo, ahi,
     $                        g_dens_tot_i(ipm),    alo, ahi)
         endif       !lifetime

         alo(1)   =  nvec+1
         ahi(1)   =  2*nvec
         clo(1)   =  1
         chi(1)   =  nvec
         clo(2)   =  1
         chi(2)   =  nbf
         clo(3)   =  1
         chi(3)   =  nbf

c     Add beta density to g_dens_tot
c     divided by 4 to match closed shell static routines. 
c     TODO: make sure that this is also being scaled in a similar fashion
c           in the closed shell FD cases

         call nga_add_patch(0.25d0, g_dens(ipm),     alo,  ahi,
     $                      0.25d0, g_dens_tot(ipm), clo,  chi,
     $                              g_dens_tot(ipm), clo,  chi)
         if (lifetime) then
            call nga_add_patch(0.25d0, g_dens(ipm+2),     alo,  ahi,
     $                         0.25d0, g_dens_tot_i(ipm), clo,  chi,
     $                                 g_dens_tot_i(ipm), clo,  chi)
         endif       !lifetime
         enddo       !ipm

c     reset dimensions back previous values

         dims(1)  =  nmul*npol*nvec
         chunk(1) =  dims(1)
         alo(1)   =  0
         ahi(1)   =  0
      endif       ! ldimqm
c  DIM/QM END

      oskel_local = oskel .and. (nvec.eq.1)
      tol2e_local = min(max(acc,itol_floor),itol_ceil)

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .true., ! in : = .true for symm dens
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

      call get_undosymm_fock(
     &                g_fock,  ! in/ou: fock matrix
     &                nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &                nvec,    ! in   : nr. directions (x,y,z)
     &                nbf,     ! in   : nr. basis functions
     &                npol,    ! in   : nr. polarizations
     &                nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &                g_pmats, ! in   : scratch GA array
     &                limag)   ! in   : =.true. imaginary comp. exists

c ------- Remove GA arrays:
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_pmats(ipm))) call errquit(
     &     'uhf_hessv_2e3: ga_destroy failed g_pmats',0,GA_ERR)       
      enddo ! end-loop-ipm

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
          if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az-0------ START' 
          call ga_print(g_Az)
          if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az-0----- END'
      endif ! end-if-debug

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 
      call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     

c     start loop over components of perturbing field
      nocc1(1)=nalpha
      nocc1(2)=nbeta
      nvir1(1)=nmo-nocc1(1)
      nvir1(2)=nmo-nocc1(2)
      do cnt=1,nset 
       do iset = 1, npol
         if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)

        if (debug) then
         if (ga_nodeid().eq.0) then
         write(*,194) nocc_max,nvir_max
 194     format('(noc_max,nvir_max)=(',i4,',',i4,')')
         endif
        endif ! end-if-debug

        if (debug) then
         if (ga_nodeid().eq.0) then
           write(*,192) cnt,iset,nocc,nocc_diff,
     &                  nvir,nvir_diff,istart,istart_diff,
     &                  vbase,vbase_diff
 192       format('(cnt,iset,nocc,nocc_diff,nvir,',
     &             'nvir_diff,istart,istart_diff,vbase,vbase_diff)=(',
     &              i4,',',i4,',',i4,',',i4,',',i4,',',i4,',',
     &              i4,',',i4,',',i4,',',i4,')')
         endif
        endif
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
        do ivec = 1, nvec
         do ipm=1,ncomp ! loop over Fock matrix components +/- here 
          ind=indx(ipm,cnt)

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,ipm,ind
 117        format('XX:(cnt,iset,ivec,ipm,ind)=(',
     &             i3,',',i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug
           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(ind),alo,ahi,
     &                             g_k        ,blo,bhi) 
           call getreorim1_u1(
     &                      g_Axreim,! out : real or im arr
     &                      g_Az,    ! in  : = complx(g_xre,g_xim)
     &                      0,       ! in  : subblock index
     &                      ipm,     ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      cnt)     ! in  : =1 -> re =2 -> im
c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

            if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-1 update_ax_fock'
            endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim, ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

             if (debug) then
              if (ga_nodeid().eq.0)
     &         write(*,*) 'AFT-1 update_ax_fock'      
              endif ! end-if-same-spin-contrib
             endif ! end-if-debug

c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

              if (debug) then
               if (ga_nodeid().eq.0) then
                write(*,109) nocc_diff,nvir_diff 
 109            format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
               endif           
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-2 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim,   ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug
           endif ! end-if-diff-spin-contrib

c ---- Move g_Axreim --> g_Az---- START
               call conv2complex4_u1(
     &                g_Az,     ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                0,        ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                cnt)      ! in : =1 -> re =2 -> im
c ---- Move g_Axreim --> g_Az---- END
         enddo ! end-loop-ipm
        enddo ! end-loop-ivec
        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        enddo ! end-loop-iset (spin A,B)
       enddo ! end-loop-cnt
c  DIM/QM jbecca START
c  adding in the DIM operator here, the way that g_axreim is handled in 
c  this routine makes things much harder compared to uhf_hessv3_2e3
      if (ldimqm) then
c         write(luout,*)'starting dimqm part of uhf_hessv3_2e3_cmplx'
         do ipm = 1, ncomp
         if (.not. ga_create(MT_DBL, n, nvec, 'g_ax_re', 0,0,
     $         g_ax_re(ipm))) call errquit('uhf_hessv3_2e3: DIM dens
     $                  create', 0, GA_ERR)
         if (.not. ga_create(MT_DBL, n, nvec, 'g_ax_im', 0,0,
     $         g_ax_im(ipm))) call errquit('uhf_hessv3_2e3: DIM dens
     $                  create', 0, GA_ERR)

c         write(luout,*)'done with allocation for',ipm
         call getreorim1_u1(
     &                      g_ax_re(ipm), ! out : real or im arr
     &                      g_Az,    ! in  : = complx(g_xre,g_xim)
     &                      0,       ! in  : subblock index
     &                      ipm,     ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      1)       ! in  : =1 -> re =2 -> im

         call getreorim1_u1(
     &                      g_ax_im(ipm), ! out : real or im arr
     &                      g_Az,    ! in  : = complx(g_xre,g_xim)
     &                      0,       ! in  : subblock index
     &                      ipm,     ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      2)       ! in  : =1 -> re =2 -> im
         enddo      !ipm

c         write(luout,*)'real and imag g_ax arrays made'

         call dimqm_addop_uhf_damp(g_ax_re, g_ax_im, ncomp, limag,
     $         lifetime, g_dens_tot, g_dens_tot_i)

         do ipm   =  1, ncomp
            call conv2complex4_u1(
     &                g_Az,     ! out: = history matrix complex
     &                g_ax_re(ipm), ! in : real      arr
     &                0,        ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                1)      ! in : =1 -> re =2 -> im

            call conv2complex4_u1(
     &                g_Az,     ! out: = history matrix complex
     &                g_ax_im(ipm), ! in : real      arr
     &                0,        ! in : subblock index
     &                ipm,      ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                2)      ! in : =1 -> re =2 -> im

            if (.not. ga_destroy(g_ax_re(ipm))) call errquit
     $               ('uhf_hessv: GA?',0,GA_ERR)
            if (.not. ga_destroy(g_ax_im(ipm))) call errquit
     $               ('uhf_hessv: GA?',0,GA_ERR)
            if (.not. ga_destroy(g_dens_tot(ipm))) call errquit
     $               ('uhf_hessv: GA?',0,GA_ERR)
            if (.not. ga_destroy(g_dens_tot_i(ipm))) call errquit
     $               ('uhf_hessv: GA?',0,GA_ERR)
        enddo     !ipm 
      endif       !ldimqm
c  DIM/QM END
        if (.not. ga_destroy(g_Axreim)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az-1------ START' 
           call ga_print(g_Az)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az-1------ END'
      endif ! end-if-debug

      if (.not. ga_destroy(g_j)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_k)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    
      end

      subroutine uhf_hessv_2e2_opt_cmplx1(
     &                         acc, 
     &                         g_z,
     &                         g_Az,
     &                         nvec,
     &                         lifetime)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c Note.- Mimic uhf_hessv_2e2() but including two sets RE and IM.
c        (g_x_re,g_x_im) (g_ax_re,g_ax_im)

      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
      integer g_z,    ! in: 
     &        g_Az,   ! in: 
     &        g_xreim,! Scratch GA array  
     &        g_Axreim! Scratch GA array
      integer m1,m2
      double precision acc      ! [input] required accuracy of products
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.

      logical oskel_local
      integer ivec
      integer nocc,nvir,g_vecs,
     &        nocc_diff,nvir_diff,g_vecs_diff
      integer npol,nocc_arr(2),istart2(2),
     &        nocc1(2),nvir1(2)
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart,istart_diff,
     &        iend  ,iend_diff,vbase,vbase_diff
      integer g_ax,g_j,g_k,
     &        g_tmp1,g_tmp2,n,gtype
      integer g_dens(2),g_fock(2) ! 2 stands for RE-IM blocks
      integer alo(3),ahi(3), 
     &        blo(2),bhi(2), 
     &        dims(3),chunk(3)
      double precision tol2e_local
      double precision itol_floor,itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         update_ax_fock,
     &         shell_fock_build2,
     &         get_dens_reorim_1,getreorim_u,
     &         getreorim1_u1,conv2complex4_u1,
     &         getreorim1_u1_debug
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
      integer nset,ncomp,nblock,nmul,istartx
      logical lifetime,debug
      integer ip,im,ipm,cnt,nocc_max,nvir_max
      data npol /2/ ! for unrestricted calculations

      debug=.false. ! no printouts

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
      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nocc_arr(1)=nalpha
      nocc_arr(2)=nbeta
      istart2(1)=1
      istart2(2)=nalpha*(nmo-nalpha) + 1

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'FA-noskew_uhf=',noskew_uhf 
      endif ! end-if-debug

      do iset = 1, npol
       g_vecs = g_movecs(iset)
       nocc   = nocc_arr(iset)
       nvir   = nmo-nocc
       istartx= istart2(iset)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         if (.not. ga_create(MT_DBL,nocc*nvir,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
         call ga_zero(g_xreim)
         call getreorim_u(g_xreim,  ! out : real or im arr
     &                    g_z,      ! in  : = complx(g_xre,g_xim)
     &                    istartx-1,! in  : shift
     &                    nvir,     ! in  : nr. virtual  MOs
     &                    nocc,     ! in  : nr. occupied MOs
     &                    1)        ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END
           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    1,         ! in   : =1 1st block RE
     &                    g_xreim,   ! in   : 
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts
        if (lifetime) then

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
         call ga_zero(g_xreim)
         call getreorim_u(g_xreim,  ! out : real or im arr
     &                    g_z,      ! in  : = complx(g_xre,g_xim)
     &                    istartx-1,! in  : shift
     &                    nvir,     ! in  : nr. virtual  MOs
     &                    nocc,     ! in  : nr. occupied MOs
     &                    2)        ! in  : =1 -> re =2 -> im
c ---- Copy g_z --> g_x_reim ------ END
           call get_dens_reorim_1(
     &                    g_dens,    ! in/ou: perturbed density matrix
     &                    2,         ! in   : =2 2nd block IM
     &                    g_xreim,   ! in   : IM
     &                    g_vecs,    ! in   : MO coefficients
     &                    nbf,       ! in   : nr. basis functions
     &                    nmo,       ! in   : nr. MOs
     &                    1,         ! in   : shift nocc-nvirt block
     &                    nocc,      ! in   : nr. occupied MOs
     &                    nvir,      ! in   : nr. virtual  MOs 
     &                    nvec,      ! in   : nr. directions (x,y,z)
     &                    iset,      ! in   : nr. polarizations
     &                    noskew_uhf,! in   : logical var
     &                    debug)     ! in   : = .true. -> debugging printouts
        endif ! end-if-lifetime
       if (.not.ga_destroy(g_xreim)) call errquit(
     &     'uhf_hessv2_cmplx: ga_destroy failed g_xreim',0,GA_ERR) 
      enddo ! end-loop-iset

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock, ! out: Fock    matrices
     &                       g_dens, ! in : density matrices
     &                       geom,   ! in : geom  handle
     &                       basis,  ! in : basis handle
     &                       nbf,    ! in : nr. basis functions
     &                       nvec,   ! in : nr. vecs (x,y,z)
     &                       npol,   ! in : nr. polarizations =1 RDFT =2 UDFT
     &                       ncomp,  ! in : nr. components
     &                       nblock, ! in : nr. of g_dens,g_fock blocks
     &                       .false.,! in : =.false. for nonsymm dens
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
     &    write(*,*) '------- g_Az-0------ START' 
          call ga_print(g_Az)
          if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az-0------ END'
       endif ! end-if-debug

      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens') 

      call ga_inquire(g_z,gtype,n,nvec) ! get (n,nvec)
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx2_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
      nocc1(1)=nalpha
      nocc1(2)=nbeta
      nvir1(1)=nmo-nocc1(1)
      nvir1(2)=nmo-nocc1(2)   
c     start loop over components of perturbing field
      do cnt=1,nset 
       do iset = 1, npol
         if (iset .eq. 1) then
            isetoff     = 0    ! index for J_A terms
            isetoff2    = nvec ! index for K_A terms
            istart      = 1
            nocc        = nalpha
            g_vecs      = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff   = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff     = 2*nvec ! index for J_B terms
            isetoff2    = 3*nvec ! index for K_B terms
            istart      = nalpha*(nmo-nalpha) + 1
            nocc        = nbeta
            g_vecs      = g_movecs(2)
            istart_diff = 1
            nocc_diff   = nalpha
            g_vecs_diff = g_movecs(1)
        endif

        nvir       = nmo - nocc
        nvir_diff  = nmo - nocc_diff
        iend       = istart     +nocc     *nvir      - 1
        iend_diff  = istart_diff+nocc_diff*nvir_diff - 1
        vbase      = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1     ! First virtual
        nnocc      = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
        nocc_max=max(nnocc,nnocc_diff)
        nvir_max=max(nvir ,nvir_diff)
c Fix code: 05-01-12: Dimensioning (g_tmp1,g_tmp2) to
c                     using max(nocc(i)) max(nvir(i))
        if (.not. ga_create(MT_DBL, nbf, nocc_max, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
        if (.not. ga_create(MT_DBL, nvir_max, nocc_max, 
     &                      'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)

        do ivec = 1, nvec
           call getreorim1_u1(
     &                      g_Axreim,! out : real or im arr
     &                      g_Az,    ! in  : = complx(g_xre,g_xim)
     &                      0,       ! in  : subblock index
     &                      1,       ! in  : = 1,2 to access slctd component
     &                      npol,    ! in  : nr. polarizations
     &                      nvir1,   ! in  : nr. virtual  MOs
     &                      nocc1,   ! in  : nr. occupied MOs
     &                      cnt)     ! in  : =1 -> re =2 -> im

          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,117) cnt,iset,ivec,1
 117        format('XX:(cnt,iset,ivec,ipm)=(',
     &             i3,',',i3,',',i3,',',i3,')')
           endif
          endif ! end-if-debug

           alo(1) = isetoff  + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_j        ,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock(cnt),alo,ahi,
     &                             g_k        ,blo,bhi) 
c     Same spin 2-e contributions    
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one,g_j,mone,g_k,g_k) ! K <-- J-K

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'BEF-1 update_ax_fock'
              endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim, ! in/ou: 
     &               g_k,      ! in   : AO Fock matrix Coul or Exch
     &               g_vecs,   ! in   : MO vecs
     &               ivec,     ! in   : ivec-th MO
     &               vbase,    ! in   : virtual base index for g_vecs
     &               nbf,      ! in   : nr. basis functions
     &               nmo,      ! in   : nr. MOs
     &               nocc,     ! in   : nr. occupied MOs
     &               nvir,     ! in   : nr. virtual  MOs
     &               istart,   ! in   : istart in g_ax(ipm)
     &               iend,     ! in   : iend   in g_ax(ipm)
     &               g_tmp1,   ! in   : scratch GA array
     &               g_tmp2,   ! in   : scratch GA array
     &               debug)

              if (debug) then 
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-1 update_ax_fock'   
              endif ! end-if-debug
    
            endif ! end-if-same-spin-contrib
c     Different spin
           if (nocc_diff*nvir_diff .gt. 0) then

             if (debug) then
              if (ga_nodeid().eq.0) then
               write(*,109) nocc_diff,nvir_diff 
 109           format('(nocc_diff,nvir_diff)=(',i3,',',i3,')')
              endif           
              if (ga_nodeid().eq.0)
     &         write(*,*) 'BEF-2 update_ax_fock'
             endif ! end-if-debug

              call update_ax_fock(
     &               g_Axreim,   ! in/ou: 
     &               g_j,        ! in   : AO Fock matrix Coul or Exch
     &               g_vecs_diff,! in   : MO vecs
     &               ivec,       ! in   : ivec-th MO
     &               vbase_diff, ! in   : virtual base index for g_vecs
     &               nbf,        ! in   : nr. basis functions
     &               nmo,        ! in   : nr. MOs
     &               nocc_diff,  ! in   : nr. occupied MOs
     &               nvir_diff,  ! in   : nr. virtual  MOs
     &               istart_diff,! in   : istart in g_ax(ipm)
     &               iend_diff,  ! in   : iend   in g_ax(ipm)
     &               g_tmp1,     ! in   : scratch GA array
     &               g_tmp2,     ! in   : scratch GA array
     &               debug)

              if (debug) then
               if (ga_nodeid().eq.0)
     &          write(*,*) 'AFT-2 update_ax_fock'
              endif ! end-if-debug

           endif ! end-if-diff-spin-contrib
c ---- Move g_Axreim --> g_Az---- START
               call conv2complex4_u1(
     &                g_Az,     ! out: = history matrix complex
     &                g_Axreim, ! in : real      arr
     &                0,        ! in : subblock index
     &                1,        ! in : = 1,2 to access slctd component
     &                npol,     ! in : nr. polarizations
     &                nvir1,    ! in : nr. virtual  MOs
     &                nocc1,    ! in : nr. occupied MOs
     &                cnt)      ! in : =1 -> re =2 -> im
c ---- Move g_Axreim --> g_Az---- END
        enddo ! end-loop-ivec
        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hv2e: tmp1',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hv2e: tmp2',0,
     &       GA_ERR)
        enddo ! end-loop-iset (spin A,B)
       enddo ! end-loop-cnt
        if (.not. ga_destroy(g_Axreim)) call errquit('uhf_hv2e: Axreim',0,
     &       GA_ERR)
      if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az-1------ START' 
           call ga_print(g_Az)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az-1------ END'
      endif ! end-if-debug

      if (.not. ga_destroy(g_j)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_k)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'uhf_hv2e: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('uhf_hv2e: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm   
 
      end
c $Id$
