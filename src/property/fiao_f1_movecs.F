      subroutine fiao_f1_movecs(
     &  rtdb, basis, geom,      ! IN
     &  omega, lstatic,         ! IN
     &  ncomp,  lvelocity,      ! IN
     &  lifetime, gamwidth,     ! IN
     &  g_vecF1, g_vecF1_im,    ! OUT
     &  index4cphf)             ! IN
c $Id$ 
c     =================================================================
c     purpose: Calculate MO vectors perturbed by a static or dynamic
c     external field.  This is a stripped-down and then
c     extended version of hnd_giaox, see the original aoresponse routine
c     for the frequency dependent extension and for the finite lifetime
c     damping.
c     g_vecF1(1:2) should have at least dimensions (nbf, nocc, 3)
c     same for the imaginary components, *if* lifetime.eq..true.
c     Note: the CPHF call terminates the integrals
c     =================================================================
c
c  Written by J. Autschbach, SUNY Buffalo
c  Extension to spin-unrestricted case 
c          by F. Aquino,     Northwestern University 
c          03-15-12

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
#include "util.fh"
#include "bgj.fh"
      
c     ---------------------
c     subroutine arguments:
c     ---------------------
      
      integer rtdb            ! [input] run-time database handle
      integer basis           ! [input] basis handle
      integer geom            ! [input] geometry handle
      integer g_vecF1(2,2)    ! [out] Field perturbed MO coefficients
      integer g_vecF1_im(2,2) ! [out] Field perturbed MOs, Imag
      integer ncomp 
      double precision omega, gamwidth
      logical lstatic, lvelocity, lifetime
      
c     local variables:
      
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, iatom, iocc, ifld, ioff
      
      integer alo(3),ahi(3), 
     &        alo1(3),ahi1(3),shift,
     &        blo(3),bhi(3), clo(3), chi(3)
      integer dlo(3),dhi(3)
      
      integer l_occ, k_occ, l_eval, k_eval
      
      integer g_u(2,2),g_u_im(2,2) ! g_u(ispin,ncomp)
      integer g_rhs(2),g_rhs_im(2) ! needed for 'lifetime'
      
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer vectors1(2)
      
      double precision atn, tol2e, val
      
      double precision a(6)
      character*255 aorespfilename
      character*3 scftyp
      character*16 tag
      character*2 symbol
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      integer nat
      parameter (nat=1)
      logical skip_cphf_ev_fiaof1 
      integer index4cphf

      integer ndir,    ! nr directions (x,y,z)
     &        ntot,    ! sum_{i=1,npol} nocc(i)*nvirt(i)
     &        ispin,
     &        nocc(2), ! store nr occupations 
     &        npol     ! nr of polarizations =1 (  restricted) 
                       !                     =2 (unrestricted) calc
      integer m1,m2
c
      logical     oskel, status, debug, dbgmat, limag
      data tol2e   /1.0d-10/
      integer ipm,read_grhs_fiaof1,nskip
      
c     nwchem file access character strings:
      character*256 cphf_rhs(2), cphf_sol(2), cphf_rhs_im(2),
     &  cphf_sol_im(2), cstemp
      logical lclfld,use_dimqm,ldimqm,lrsp
c     external functions:
      logical  cphf2, cphf3, file_write_ga, file_read_ga, cphf
      external cphf2, cphf3, file_write_ga, file_read_ga, cphf

c  ====================================================================
      ndir=3 ! nr directions (x,y,z)
      debug  = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices
      use_dimqm=util_module_avail("dimqm")

c      debug=.true.

      if (debug) write (luout,*) 'hello from fiao_f1_movecs'

      limag = lvelocity         ! imaginary perturbation ?     
c     -------------------------
c     set-up for CPKS procedure
c     -------------------------
      if (lifetime) then
        if (ga_nodeid().eq.0) write (LuOut,'(a/a,f15.8,a/a,f15.8,a/)')
     &     ' Setting up CPKS',
     &     ' with frequency omega = ',omega,' a.u.',
     &     ' and a damping of       ',gamwidth,' a.u.'
      else
        if (ga_nodeid().eq.0) write (LuOut,'(a,f15.8,a/)')
     &     ' Setting up CPKS'//
     &     ' with frequency omega = ',omega,' a.u.'
      end if
      if (lstatic .and. ga_nodeid().eq.0) write (LuOut,'(a/)')
     &   ' STATIC response' 
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
c     
      oskel = .false.
c
c     Integral initialization 
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
      status = rtdb_parallel(.true.)
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues
      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('fiao_f1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('fiao_f1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('fiao_f1: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)

      call get_nocc(rtdb,   ! in : rtdb handle
     &              nocc,   ! out: nr occupations
     &              npol,   ! out: nr of polarization
     &              nclosed,! in : nr closed shells
     &              nopen,  ! in : nr open shells
     &              nvirt,  ! in : nr virtual MOs
     &              scftyp, ! in : string = UHF or RHF
     &              ntot)   ! out: sum_{i,npol} nocc(i)*nvirt(i)

c      if (ga_nodeid().eq.0) then
c        write(*,10) npol,nocc(1)   ,nocc(2),
c     &              nopen(1)  ,nopen(2),
c     &              nclosed(1),nclosed(2),
c     &              nvirt(1)  ,nvirt(2),scftyp,ntot
c 10    format('npol=',i3,' nocc =(',i3,',',i3,') ',
c     &        'nopen=(',i3,',',i3,') ',
c     &        'nclos=(',i3,',',i3,') ',
c     &        'nvirt=(',i3,',',i3,') ',
c     &        'scftyp=',a,' ntot=',i3)
c      endif
    
c      if (ga_nodeid().eq.0)
c     &  write(*,*) 'STOP-test-urest'
c      stop
     
c ------ define npol ----- END
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('fiao_f1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('fiao_f1: ma_pop_stack failed k_occ',0,MA_ERR)

      if (debug) 
     &   write (luout,*) 'fiao_f1: nclosed,nopen,nvirt',nclosed(1),
     &   nopen(1), nvirt(1)
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c     ipm counts the plus and minus combinations according to
c     Santry's notation. ipm=1: plus, ipm=2: minus
c     ipm = 1 corresponds to F(k,i)(+)
c     ipm = 2 corresponds to F(k,i)(-) = F(i,k)(+)*
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
      do ipm = 1,ncomp
        if(.not.ga_create(MT_DBL,ntot,ndir,
     &     'RHS',-1,-1,g_rhs(ipm)))
     &     call errquit('fiao_f1: ga_create failed g_rhs',0,GA_ERR)
        call ga_zero(g_rhs(ipm))
        if (lifetime) then
          if(.not.ga_create(MT_DBL,ntot,ndir,
     &       'RHS',-1,-1,g_rhs_im(ipm)))
     &       call errquit('fiao_f1: ga_create failed g_rhs',0,GA_ERR)
          call ga_zero(g_rhs_im(ipm))
        end if                  ! lifetime
      end do                    ! ipm
      if (debug) write (luout,*) 'fiao_f1: rhs arrays allocated'
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'FA-BEF-update_rhs_dipole'
c      debug=.true.
c
c     JEM: Updated to include imaginary portion of rhs
      call update_rhs_dipole(
     &                     g_rhs,    ! in/out:
     &                     g_rhs_im, ! in/out:
     &                     vectors,  ! in : MO vectors
     &                     rtdb,     ! in : rtdb  handle
     &                     basis,    ! in : basis handle
     &                     lvelocity,! in : logical var
     &                     nat,      ! in : nr. atoms
     &                     npol,     ! in : nr. polarizations
     &                     nocc,     ! in : nr. occ  shells
     &                     nvirt,    ! in : nr. virt shells
     &                     nbf,      ! in : nr. basis functions
     &                     nmo,      ! in : nr. MOs   
     &                     ncomp,    ! in : nr components of ...
     &                     debug,    ! in : logical for debugging
     &                     omega,    ! in : freq value
     &                     lifetime) ! in : logical for debugging

      call schwarz_tidy()
      call int_terminate()
      if (debug) then
       do ipm = 1,ncomp
        if (ga_nodeid().eq.0) then
         write(*,201) ipm
  201    format('---- g_rhs-AFT-udipole(',i3,')-------- START')
        endif
        call ga_print(g_rhs(ipm))
        if (ga_nodeid().eq.0) then
         write(*,202) ipm
  202    format('---- g_rhs-AFT-udipole(',i3,')-------- END')
        endif
       enddo ! end-loop-ipm
        do ispin=1,npol
          if (ga_nodeid().eq.0) 
     &     write(*,*) '--C--vectors(',ispin,')--- START'
          call ga_print(vectors(ispin))
          if (ga_nodeid().eq.0) 
     &     write(*,*) '--C--vectors(',ispin,')--- END'
        enddo ! end-loop-ispin     
      endif ! end-if-debug
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'FA-AFT-update_rhs_dipole'
c     there is a division by -4 somewhere in the
c     CPKS solver so we scale the rhs by -4 in order to
c     get the correct magnitude of the result back   
      blo(1) = 1
      bhi(1) = ntot
      blo(2) = 1
      bhi(2) = ndir  
      do ipm = 1,ncomp
        call nga_scale_patch(g_rhs(ipm),blo,bhi,-4.0d0)
c       DIM/QM JEM
      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlrsp(lrsp)
         call dimqm_getlclfld(lclfld)
      else
         ldimqm=.false.
         lrsp=.false.
         lclfld=.false.
c
c
c

      endif
        if (ldimqm .and. lrsp .and. lclfld .and. lifetime) then
          call nga_scale_patch(g_rhs_im(ipm),blo,bhi,-4.0d0)
        end if
      end do

c --> Put in rtdb to be read in ga_lkain_2cpl3_damp_cmplx_redmem()
c             or  to be read in ga_lkain_2cpl3_damp_cmplx_redmem1()
      if (.not. rtdb_put(rtdb,'cphf3-aores:guess1', 
     &          mt_int, 1,index4cphf)) call
     $     errquit('fiao_f1: failed to write index4cphf ', 0, RTDB_ERR)

c         if (ga_nodeid().eq.0)         
c     &   write(*,*) 'COMPUTE cphf fiao_f1 data ...'
c ====== writing g_rhs to file ============== START
c     Write ga_rhs to disk 
       call cphf_fname('cphf_rhs',cphf_rhs(1))
       call cphf_fname('cphf_sol',cphf_sol(1))
       if(.not.file_write_ga(cphf_rhs(1),g_rhs(1))) call errquit
     $   ('fiao_f1: could not write cphf_rhs(1)',0, DISK_ERR)
       if (ncomp.gt.1) then
        call cphf_fname('cphf_rhs_2',cphf_rhs(2))
        call cphf_fname('cphf_sol_2',cphf_sol(2))
        if(.not.file_write_ga(cphf_rhs(2),g_rhs(2))) call errquit
     $     ('fiao_f1: could not write cphf_rhs 2',0, DISK_ERR)
       endif
        if (lifetime) then     
c       mainly, there are only zeros written to file here. i.e.
c       we start up the CPKS with zero imaginary parts. Probably
c       not the best guess near resonance but the preconditioner
c       should take some care of that    
        call cphf_fname('cphf_rhs_im',cphf_rhs_im(1))
        call cphf_fname('cphf_sol_im',cphf_sol_im(1))
        if(.not.file_write_ga(cphf_rhs_im(1),g_rhs_im(1)))
     +     call errquit
     $     ('fiao_f1: could not write cphf_rhs_im 1',0, DISK_ERR)      
        if (ncomp.gt.1) then
          call cphf_fname('cphf_rhs_im_2',cphf_rhs_im(2))
          call cphf_fname('cphf_sol_im_2',cphf_sol_im(2))
          if(.not.file_write_ga(cphf_rhs_im(2),g_rhs_im(2)))
     +       call errquit
     $       ('fiao_e2: could not write cphf_rhs_im 2',0, DISK_ERR)
        endif
      endif                     ! lifetime
c ====== writing g_rhs to file ============== END
c
c     Call the appropriate CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.

      if (lvelocity) then       ! we do have skew-sym matrix
        if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $     errquit('fiao_f1: failed to write skew ', 0, RTDB_ERR)
      else                      ! dipole-length not skew-symmetric
        if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.true.)) call
     $     errquit('fiao_f1: failed to write skew ', 0, RTDB_ERR)
      end if
c     we also need to tell the dynamic CPKS solver if the perturbation
c     is imaginary and how many components to treat     
      if (.not. rtdb_put(rtdb, 'cphf:imag', mt_log, 1, limag)) call
     &   errquit('fiao_f1: failed to write limag ', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'cphf:ncomp', mt_int, 1, ncomp)) call
     &   errquit('fiao_f1: failed to write ncomp ', 0, RTDB_ERR)    
c     ----------------------------------------------------
c     Call the CPKS solver. It will solve the sets of
c     equations for each component of the perturbing field
c     simultaneously.
c     ----------------------------------------------------   
      if (debug) write (luout,*) 'fiao_f1: calling CPKS'   

c      if (ga_nodeid().eq.0) then
c       write(*,222) lifetime,lstatic
c 222     format('INPUT-4: (lifetime,lstatic)=(',
c     &         L1,',',L1,')')
c      endif   
  
      if (lstatic .and. .not.lifetime) then
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-BEF cphf2()'
c       ... static CPKS solver:         
        if (.not.cphf2(rtdb)) call errquit
     &       ('fiao_f1: failure in cphf2 ',0, RTDB_ERR)
      else       
c       ... dynamic CPKS solver:
        if (debug) write (luout,*)
     &     'calling cphf3',omega,lifetime,gamwidth

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-BEF cphf3()'

        if (.not.cphf3(rtdb, omega, lifetime, gamwidth))
     +     call errquit
     $     ('fiao_f1: failure in cphf3 ',0, RTDB_ERR)
      endif   
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('fiao_f1: failed to delete skew', 0, RTDB_ERR)
      if (.not. rtdb_delete(rtdb, 'cphf:imag')) call
     &   errquit('fiao_f1: failed to delete limag', 0, RTDB_ERR)
      if (.not. rtdb_delete(rtdb, 'cphf:ncomp')) call
     &   errquit('fiao_f1: failed to delete ncomp', 0, RTDB_ERR)    
      if (debug) write (luout,*) 'fiao_f1: back from CPKS'
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
c ---- Read g_rhs from files --------- START
      blo(1) = 1
      bhi(1) = ntot
      blo(2) = 1
      bhi(2) = ndir
      do ipm=1,ncomp
        call ga_zero(g_rhs(ipm))
        if(.not.file_read_ga(cphf_sol(ipm),g_rhs(ipm)))
     &     call errquit
     $     ('fiao_f1: could not read cphf solution',0, DISK_ERR)  
        if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,203) ipm
  203      format('---- g_rhs-AFT-cphf(',i3,')-------- START')
          endif
          call ga_print(g_rhs(ipm))
          if (ga_nodeid().eq.0) then
           write(*,204) ipm
  204      format('---- g_rhs-AFT-cphf(',i3,')-------- END')
          endif    
        endif ! end-if-debug        
        if (lifetime) then
          if(.not.file_read_ga(cphf_sol_im(ipm),g_rhs_im(ipm)))
     &       call errquit
     $       ('fiao_f1: could not read cphf solution Im',
     &       ipm, DISK_ERR)
          
         if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,205) ipm
  205      format('---- g_rhs_im-AFT-cphf(',i3,')-------- START')
          endif
          call ga_print(g_rhs_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,206) ipm
  206      format('---- g_rhs_im-AFT-cphf(',i3,')-------- END')
          endif            
         endif ! end-if-debug

        endif ! end-if-lifetime (damping)
      enddo ! end-loop-ipm
      
       call get_vecF1(
     &              g_vecF1,    ! out:
     &              g_vecF1_im, ! out:
     &              g_rhs,      ! in :
     &              g_rhs_im,   ! in :
     &              vectors,    ! in : MO vectors
     &              nbf,        ! in : nr. basis functions
     &              nmo,        ! in : nr. MOs
     &              ncomp,      ! in :
     &              npol,       ! in : nr. polarizations
     &              lifetime,   ! in : = (.true.,.false.) with/out damping
     &              nocc,       ! in : nr. occupied MOs
     &              nvirt,      ! in : nr. virtual  MOs
     &              debug)      ! in : = .true. for debugging

      do ispin=1,npol
c ========== visualize (g_vecF1,g_vecF1_im) ====== START
       if (debug) then
        do ipm = 1,ncomp
         if (ga_nodeid().eq.0) then
           write(*,5) ispin,ipm 
 5         format('---- g_vecF1(',i3,',',i3,')-------- START')
         endif
         call ga_print(g_vecF1(ispin,ipm))
         if (ga_nodeid().eq.0) then
           write(*,6) ispin,ipm 
 6         format('---- g_vecF1(',i3,',',i3,')-------- END')
         endif
         if (lifetime) then
          if (ga_nodeid().eq.0) then
           write(*,7) ispin,ipm 
 7         format('---- g_vecF1_im(',i3,',',i3,')-------- START')
          endif
          call ga_print(g_vecF1_im(ispin,ipm))
          if (ga_nodeid().eq.0) then
           write(*,8) ispin,ipm 
 8         format('---- g_vecF1_im(',i3,',',i3,')-------- END')
          endif
         endif ! end-if-lifetime
        enddo ! end-loop-ipm
       endif ! end-if-debug

c ========== visualize (g_vecF1,g_vecF1_im) ====== END
       if (.not.ga_destroy(vectors(ispin))) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
      enddo ! end-loop-ispin
c      RHS arrays are no longer needed
       do ipm=1,ncomp
        if (.not.ga_destroy(g_rhs(ipm))) call 
     &     errquit('fiao_f1: ga_destroy failed g_rhs',0,GA_ERR)
        if (lifetime) then
          if (.not.ga_destroy(g_rhs_im(ipm))) call 
     &       errquit('fiao_f1: ga_destroy failed g_rhs',0,GA_ERR)
        end if
c     delete scratch files:
        if(ga_nodeid().eq.0) then
         call util_file_unlink(cphf_rhs(ipm))
         call util_file_unlink(cphf_sol(ipm))
        if (lifetime) then
         call util_file_unlink(cphf_rhs_im(ipm))
         call util_file_unlink(cphf_sol_im(ipm))
        endif
        endif
       enddo ! end-loop-ipm
      call ga_sync()
      return

 7000 format(/,10x,'Field perturbed MOs cannot be calculated for',
     1      ' UHF or ROHF wave functions at present')
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9999 format(
     1 /,10x,54(1h-),/,
     2 10x,'Calculating external field perturbed MO vectors',/,
     3 10x,54(1h-),/)

      end
