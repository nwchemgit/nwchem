
      subroutine fiao_f1_movecs (rtdb, basis, geom, omega, lstatic,
     &   ncomp,  lvelocity, lifetime, gamwidth,g_vecF1, g_vecF1_im)

c $Id:  $
      
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

c     ---------------------
c     subroutine arguments:
c     ---------------------

      integer rtdb          ! [input] run-time database handle
      integer basis         ! [input] basis handle
      integer geom          ! [input] geometry handle
      integer g_vecF1(2)    ! [out] Field perturbed MO coefficients
      integer g_vecF1_im(2) ! [out] Field perturbed MOs, Imag

      integer ncomp 
      double precision omega, gamwidth
      logical lstatic, lvelocity, lifetime

c     local variables:

      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval

      integer g_dipole, g_rhs(2), g_u(2)
      integer g_u_im(2), g_rhs_im(2) ! needed for 'lifetime'
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val
      double precision jfac(3),kfac(3),a(6),xfac
      character*3 scftyp
      character*16 tag
      character*2 symbol

      double precision origin(3)
      data origin/0d0,0d0,0d0/

      integer nat
      parameter (nat=1)

      logical     oskel, status, debug, dbgmat, limag
      data tol2e   /1.0d-10/

      integer ipm

c     nwchem file access character strings:

      character*256 cphf_rhs(2), cphf_sol(2), cphf_rhs_im(2),
     &   cphf_sol_im(2), cstemp

c     external functions:

      logical  cphf2, cphf3, file_write_ga, file_read_ga, cphf
      external cphf2, cphf3, file_write_ga, file_read_ga, cphf

c  ====================================================================

      debug = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices

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
c     If DFT get part of the exact exchange defined
c
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
c
c     Integral initialization 
c
c      call int_init(rtdb,1,basis)
c      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
      status = rtdb_parallel(.true.)

c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues
c
      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('fiao_f1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('fiao_f1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('fiao_f1: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('fiao_f1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('fiao_f1: ma_pop_stack failed k_occ',0,MA_ERR)
c
c     Get Unperturbed Density Matrix
c
c$$$      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
c$$$     &                      nclosed,nopen,nvirt)
c
c     Error exit if scftyp equals UHF (= ROHF)
c
      if (scftyp.eq.'UHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          call errquit('fiao_f1: incompatible SCF type for Response',
     &       0,INPUT_ERR)
      endif

      if (debug) 
     &   write (luout,*) 'fiao_f1: nclosed,nopen,nvirt',nclosed(1),
     &   nopen(1), nvirt(1)

c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c     ipm counts the plus and minus combinations according to
c     Santry's notation. ipm=1: plus, ipm=2: minus
c     ipm = 1 corresponds to F(k,i)(+)
c     ipm = 2 corresponds to F(k,i)(-) = F(i,k)(+)*

      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nclosed(1)
      ahi(3) = 3

      do ipm = 1,ncomp
        write (cstemp,'(a,i1)') 'U matrix ',ipm
        if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),alo,g_u(ipm)))
     &     call errquit('fiao_f1: nga_create failed g_u',0,GA_ERR)
        call ga_zero(g_u(ipm))

        if (lifetime) then
          write (cstemp,'(a,i1)') 'U matrix Im ',ipm
          if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),alo,
     &       g_u_im(ipm))) 
     &     call errquit('fiao_f1: nga_create failed g_u_im',0,GA_ERR)
          call ga_zero(g_u_im(ipm))
        end if                  ! lifetime
      end do                    ! ipm

      if (debug) write (luout,*) 'fiao_f1: U matrices allocated'

c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)

      do ipm = 1,ncomp
        if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,
     &     'RHS',-1,-1,g_rhs(ipm)))
     &     call errquit('fiao_f1: ga_create failed g_rhs',0,GA_ERR)
        call ga_zero(g_rhs(ipm))

        if (lifetime) then
          if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,
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
c
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
c    
c     Get dipole integrals in GA and transform to MO set (virt,occ)
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'dip matrix',alo,g_dipole)) call 
     &    errquit('fiao_f1: nga_create failed g_dipole',0,GA_ERR)

      if (debug) write (luout,*) 'fiao_f1: dipole matrix allocated'


c     Get H10 in GA, using g_dipole array
c     note: origin has been set to (0,0,0) for multipole integs.

      call ga_zero(g_dipole)
      if (lvelocity) then 
        call int_giao_1ega(basis,basis,g_dipole,'velocity',origin,
     &     nat, oskel)
        call ga_scale (g_dipole, -1d0)
      else
        call int_mpole_1ega(basis,basis,g_dipole,'dipole',origin,
     &     oskel)
      end if

      if (debug) write (luout,*) 'fiao_f1: AO integrals done'

c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs.
c     (the rhs is NOT divided by (e_a - e_i -/+ omega), this
c     will be considered in the CPKS solver, in the precon-
c     ditioner and the 1e part of the "product" routine)

      call giao_aotomo(g_dipole,vectors,nclosed,nvirt,1,3,nbf)
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3

      do ipm = 1,ncomp
        call nga_add_patch(1.0d0,g_rhs(ipm),blo,bhi,1.0d0,g_dipole,alo
     &     ,ahi,g_rhs(ipm),blo,bhi)
      end do

      if (debug) write (luout,*) 'fiao_f1: dipole added to rhs'

c     Cleanup g_dipole as we do not need it right now

      if (.not.ga_destroy(g_dipole)) call 
     &    errquit('fiao_f1: ga_destroy failed g_dipole',0,GA_ERR)


c     there is a division by -4 somewhere in the
c     CPKS solver so we scale the rhs by -4 in order to
c     get the correct magnitude of the result back
      
      do ipm = 1,ncomp
        call nga_scale_patch(g_rhs(ipm),blo,bhi,-4.0d0)
      end do

c     Write ga_rhs to disk 

      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs(1))
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol(1))
      if(.not.file_write_ga(cphf_rhs(1),g_rhs(1))) call errquit
     $  ('fiao_f1: could not write cphf_rhs(1)',0, DISK_ERR)

      if (ncomp.gt.1) then
        call util_file_name('cphf_rhs_2',.true.,.true.,cphf_rhs(2))
        call util_file_name('cphf_sol_2',.true.,.true.,cphf_sol(2))
        if(.not.file_write_ga(cphf_rhs(2),g_rhs(2))) call errquit
     $     ('fiao_f1: could not write cphf_rhs 2',0, DISK_ERR)
      endif

      if (lifetime) then
c       
c       mainly, there are only zeros written to file here. i.e.
c       we start up the CPKS with zero imaginary parts. Probably
c       not the best guess near resonance but the preconditioner
c       should take some care of that
        
        call util_file_name('cphf_rhs_im',.true.,.true.,
     +     cphf_rhs_im(1))
        call util_file_name('cphf_sol_im',.true.,.true.,
     +     cphf_sol_im(1))
        if(.not.file_write_ga(cphf_rhs_im(1),g_rhs_im(1)))
     +     call errquit
     $     ('fiao_f1: could not write cphf_rhs_im 1',0, DISK_ERR)
        
        if (ncomp.gt.1) then
          call util_file_name('cphf_rhs_im_2',.true.,.true.,
     +       cphf_rhs_im(2))
          call util_file_name('cphf_sol_im_2',.true.,.true.,
     +       cphf_sol_im(2))
          if(.not.file_write_ga(cphf_rhs_im(2),g_rhs_im(2)))
     +       call errquit
     $       ('fiao_e2: could not write cphf_rhs_im 2',0, DISK_ERR)
        endif
      endif                     ! lifetime

      if (debug) write (luout,*) 'fiao_f1: cphs rhs written to file'
c
      call schwarz_tidy()
      call int_terminate()
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
c     
      if (debug) write (luout,*) 'fiao_f1: calling CPKS' 
      
      if (lstatic .and. .not.lifetime) then
c       ... static CPKS solver:         
        if (.not.cphf2(rtdb)) call errquit
     &       ('fiao_f1: failure in cphf2 ',0, RTDB_ERR)
      else       
c       ... dynamic CPKS solver:
        if (debug) write (luout,*)
     &     'calling cphf3',omega,lifetime,gamwidth
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

      do ipm = 1,ncomp
        call ga_zero(g_rhs(ipm))
        if(.not.file_read_ga(cphf_sol(ipm),g_rhs(ipm)))
     &     call errquit
     $     ('fiao_f1: could not read cphf solution',0, DISK_ERR)      
        call nga_copy_patch('n',g_rhs(ipm),blo,bhi,g_u(ipm),alo,ahi)
        
        if (lifetime) then
          if(.not.file_read_ga(cphf_sol_im(ipm),g_rhs_im(ipm)))
     &       call errquit
     $       ('fiao_f1: could not read cphf solution Im',
     &       ipm, DISK_ERR)
          call nga_copy_patch('n',g_rhs_im(ipm),blo,bhi,
     &       g_u_im(ipm),alo,ahi)
        endif                   ! lifetime
      end do                    ! ncomp

c     RHS arrays are no longer needed

      do ipm = 1,ncomp
        if (.not.ga_destroy(g_rhs(ipm))) call 
     &     errquit('fiao_f1: ga_destroy failed g_rhs',0,GA_ERR)
        if (lifetime) then
          if (.not.ga_destroy(g_rhs_im(ipm))) call 
     &       errquit('fiao_f1: ga_destroy failed g_rhs',0,GA_ERR)
        end if
      end do

c     delete scratch files:

      call util_file_unlink(cphf_rhs(1))
      call util_file_unlink(cphf_sol(1))
      if (lifetime) then
        call util_file_unlink(cphf_rhs_im(1))
        call util_file_unlink(cphf_sol_im(1))
      endif
      if (ncomp.gt.1) then
        call util_file_unlink(cphf_rhs(2))
        call util_file_unlink(cphf_sol(2))
        if (lifetime) then
          call util_file_unlink(cphf_rhs_im(2))
          call util_file_unlink(cphf_sol_im(2))
        endif
      endif

c     --------------------------------------------
c     allocate main result arrays and compute them
c     --------------------------------------------

      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      do ipm = 1,ncomp
       write (cstemp,'(a,i1)') 'aor vecF1 ',ipm
        if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),
     &    alo,g_vecF1(ipm) ))
     &     call errquit('fiao_f1: nga_create failed vecF1',0,GA_ERR)
        call ga_zero(g_vecF1(ipm))

        if (lifetime) then
          write (cstemp,'(a,i1)') 'aor vecF1_Im ',ipm
          if (.not.nga_create(MT_DBL,3,ahi,trim(cstemp),
     &       alo,g_vecF1_im(ipm) ))
     &       call errquit('fiao_f1: nga_create failed E1_im',0,GA_ERR)
          call ga_zero(g_vecF1_im(ipm))
        end if
      end do                    ! ncomp

c     C1 = C0 * U10

      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
c$$$      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1)) call 
c$$$     &   errquit('fiao_f1: nga_create failed g_d1',0,GA_ERR)
c     
      alo(1) = 1
      alo(2) = 1
      blo(1) = 1
      blo(2) = 1
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(1) = nbf
      dhi(2) = nclosed(1)
      do ifld = 1, 3
        alo(3) = ifld
        ahi(3) = ifld
        blo(3) = ifld
        bhi(3) = ifld
        clo(3) = ifld
        chi(3) = ifld
        dlo(3) = ifld
        dhi(3) = ifld
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nclosed(1)
c       
c       Make C1
c       
        
        do ipm = 1,ncomp
          call nga_matmul_patch(
     &       'n','n',1.0d0,0.0d0,
     &       vectors(1),blo,bhi,
     &       g_u(ipm),alo,ahi,
     &       g_vecF1(ipm),dlo,dhi)  
          
          if (lifetime) then
            call nga_matmul_patch(
     &         'n','n',1.0d0,0.0d0,
     &         vectors(1),blo,bhi,
     &         g_u_im(ipm),alo,ahi,
     &         g_vecF1_im(ipm),dlo,dhi)  
          end if                ! lifetime
          
        end do                  ! ncomp
        
      enddo                     ! ifld

c     that's it for the heavy lifting. Clean up workspace
c     and return. We do NOT deallocate the g_vecF1 arrays, 
c     this must be done in the calling routine after using the
c     data !!!

      do ipm = 1,ncomp
        if (.not.ga_destroy(g_u(ipm))) call 
     &     errquit('fiao_f1: ga_destroy failed g_u',0,GA_ERR)
        if (lifetime) then
          if (.not.ga_destroy(g_u_im(ipm))) call 
     &       errquit('fiao_f1: ga_destroy failed g_u_im',0,GA_ERR)
        end if                  ! lifetime
      end do                    ! ncomp

      if (.not.ga_destroy(vectors(1))) call 
     &   errquit('fiao_f1: ga_destroy failed vectors',0,GA_ERR)

c$$$      if (.not.ga_destroy(g_dens(1))) call 
c$$$     &    errquit('fiao_f1: ga_destroy failed g_dens',0,GA_ERR)

c     All done.
c     At this point, we don't need to terminate the integrals.
c     They were terminated by the cphf calls.

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





