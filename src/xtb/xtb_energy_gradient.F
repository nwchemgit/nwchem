      logical function xtb_energy_gradient(rtdb)
#ifndef USE_TBLITE
      implicit none
#include "errquit.fh"
      integer rtdb
      call errquit('Recompile NWChem with TBLITE support',0,0)
#else
      use mctc_env, only : wp
      use mctc_io, only: structure_type, new
      use tblite_wavefunction_type, only : wavefunction_type, 
     $                                     new_wavefunction
      use tblite_context_type, only: context_type
      use tblite_xtb_calculator, only : xtb_calculator
      use tblite_xtb_gfn2, only : new_gfn2_calculator, export_gfn2_param
      use tblite_xtb_gfn1, only : new_gfn1_calculator, export_gfn1_param
      use tblite_xtb_singlepoint, only : xtb_singlepoint
      use tblite_param, only : param_record
#ifdef USE_OPENMP
      use omp_lib, only : omp_get_max_threads
#endif
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "util_params.fh"
#include "sym.fh"
#include "xtb.fh"

      integer rtdb
      integer geom, i, j, iat
      logical converged, status, oprint, oparam, ovecs
      double precision walltime, cputime

      integer istat
      integer icharge,lcharge,itags,ltags,iact,lact
      integer ncenters,nat,dummybq,atn

      integer iuhf8

      logical oskel, restart
      integer nactive

      character(len=16) aname,elem
      character(len=2) symb
      character(len=8) theory
      character(len=255) temp

      INTGR4 :: iuhf4, verb4, nspin4
      INTGR4,allocatable :: atnum4(:)

      double precision energy, sigma(3,3), dipole(3)
      double precision, allocatable :: coords(:,:),grad(:,:)
      double precision, allocatable :: rcov(:), trans(:,:), cn(:)

      type(structure_type),allocatable    :: mol
      type(xtb_calculator),allocatable    :: calc
      type(wavefunction_type),allocatable :: wfn
      type(context_type),allocatable      :: ctx
      type(param_record),allocatable      :: param

#ifdef USE_OPENMP
      integer maxthreads
      maxthreads = omp_get_max_threads()
      call util_blas_set_num_threads(maxthreads)
#endif

      xtb_energy_gradient = .false.
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      !
      ! Printing options
      !
      call util_print_push()
      call util_print_rtdb_load(rtdb,'xtb')
      oprint = util_print('information',print_low).and.ga_nodeid().eq.0
      oparam = ga_nodeid().eq.0 .and.util_print('params',print_high)
      ovecs = ga_nodeid().eq.0 .and. util_print('vecs',print_default)

      call xtb_rdinput(rtdb)

      if (oprint) then
        write(luout,*)
        call util_print_centered(6,
     $       'NWChem TBlite interface',
     $       40,.true.)
        write(luout,*)   
      endif

      !
      ! Initialize
      !
      if (.not.rtdb_put(rtdb, 'xtb:energy', mt_dbl, 1, 0.0d0))
     $  call errquit('xtb_energy: faild to store energy',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'xtb:converged', mt_log, 1, .false.))
     $  call errquit('xtb_energy: failed to store converged',0,RTDB_ERR)

      !
      ! Load geometry information
      !
      if (.not.geom_create(geom, 'geometry'))
     $  call errquit('xtb_energy: could not create geom',0,GEOM_ERR) 
      if (.not.geom_rtdb_load(rtdb, geom, 'geometry'))
     $  call errquit('xtb_energy: no geometry',0,RTDB_ERR)
      if (.not. geom_ncent(geom, ncenters))
     $  call errquit('xtb_energy: geom_ncent failed',0,GEOM_ERR) 
      if (.not.rtdb_get(rtdb,'charge',mt_dbl,1,rcharge))
     $  rcharge = 0d0


      iuhf4 = uhf
      nspin4 = nspin
      verb4 = verb

C     symmetry or not
      if (.not. rtdb_get(rtdb, 'gradients:use symmetry', mt_log, 1,
     $     oskel)) then
         if (.not. rtdb_get(rtdb, 'scf:skeleton', mt_log, 1,
     $        oskel)) then
            oskel = sym_number_ops(geom) .gt. 0
         end if
      end if
c     
      if (oskel .and. ga_nodeid().eq.0 .and. 
     $     util_print('information', print_default)) then
         write(luout,*) ' Using symmetry'
         call util_flush(luout)
      end if

      !
      ! Allocate space for atomic coordinates
      !
      allocate(coords(3,ncenters),stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to alloc coordinates',0,MA_ERR)
      allocate(grad(3,ncenters),stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to alloc gradient',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,ncenters,'charges',lcharge,icharge))
     $  call errquit('xtb_energy: failed to alloc charges',0,MA_ERR)
      if (.not.ma_push_get(mt_byte,ncenters*16,'tags',ltags,itags))
     $  call errquit('xtb_energy: failed to alloc tags',0,MA_ERR)
      if (.not.ma_push_get(mt_log,ncenters,'active',lact,iact))
     $  call errquit('xtb_energy: failed to alloc tags',0,MA_ERR)

      !
      ! Load geometry
      !
      if (.not. geom_cart_get(geom, ncenters, byte_mb(itags),
     $            coords, dbl_mb(icharge)))
     $  call errquit('xtb_energy: geom_cart_get failed',0,GEOM_ERR)

      dummybq = 0
      nat = 0
      do i=1,ncenters
        call util_set_ma_char_arrays(16, byte_mb(itags+16*(i-1)),
     $                               aname)
        if (aname(1:1).eq.'X' .or. aname(1:1).eq.'x' ) then
        else if(.not.geom_tag_to_element(aname, symb, elem, atn)) then
          if (symb.ne.'bq') then
            call errquit('xtb_energy: center neither atom nor bq',0,
     $                    GEOM_ERR)
          else
            dummybq = dummybq + 1
          endif
        else
          nat = nat + 1
        endif

      enddo 

      if (dummybq.gt.0)
     $  call errquit('xtb_energy: not ready for Bqs',0,GEOM_ERR)

      if (nat.ne.ncenters)
     $  call errquit('xtb_energy: not ready for dummy atoms',0,GEOM_ERR)

      !
      ! Atomic numbers
      !
      allocate(atnum4(ncenters),stat=istat)
      if (istat.ne.0)
     $  call errquit('xtb_energy: could not allocate numbers',0,MA_ERR)
      do i=0,ncenters-1
        atnum4(i+1) = nint(dbl_mb(icharge+i))
      enddo

      energy = 0.0d0
      grad(:,:) = 0.0d0
      sigma(:,:) = 0.0d0

      !
      ! This line is needed for MOLDEN
      !
      if (oprint) then
        write(luout,'(10X,A,I6)') 'No. of atoms     :',ncenters
        write(luout,*)
      endif

      !
      ! Allocate xTB objects
      !
      allocate(mol)
      allocate(param)
      allocate(wfn)
      allocate(calc)
      allocate(ctx)

      if (ga_nodeid().eq.0) then
        status = rtdb_parallel(.false.)

        !
        ! Create mol object
        !
        call new(mol,atnum4,coords,charge=rcharge,uhf=iuhf4)

        !
        ! Select the calculator type and get params
        !
        if (gfn2) then
          theory = 'xTB GFN2'
          call new_gfn2_calculator(calc, mol)
          call export_gfn2_param(param)
        elseif (gfn1) then
          theory = 'xTB GFN1'
          call new_gfn1_calculator(calc, mol)
          call export_gfn1_param(param)
        endif

        !
        ! Initialize wavefunction object
        !
        call new_wavefunction(wfn,mol%nat,calc%bas%nsh,calc%bas%nao,
     $                      nspin4,kt)

        !
        ! Compute starting charges 
        !
        status = rtdb_get(rtdb,'xtb:restart',mt_log,1,restart)
        if (guess.and..not.status) then
          call xtb_guess(mol,wfn%qat,wfn%qsh,wfn%dpat,wfn%qpat,
     $                   param,mol%nat,
     $                   calc%bas%nsh,calc%bas%nsh_at,atnum4)
        elseif(status) then
          call xtb_restart_read(mol, calc, wfn)
        endif

        !
        ! Setup spin-unrestricted calculation
        !
        if (nspin.eq.2) then
          call xtb_spinpol(mol, calc, calc%bas)
        endif

        !
        ! Set other parameters
        !
        calc%mixer_damping = broyden

        !
        ! Print params
        !
        if (oparam) then
          call xtb_print_param(calc, param)
        endif

        !
        ! Obtain energy and gradients
        !
        call xtb_singlepoint(ctx,mol,calc,wfn,acc,energy,grad,sigma,
     $                       verbosity=verb4)


        !
        ! Total energy output
        !
        write(luout,*)
        write(luout,'(A,F15.8)') '          Total xTB energy = ',energy
        write(luout,*)

        !
        ! zero force contributions on inactive atoms
        !
        call grad_active_atoms(rtdb, nat, log_mb(iact), nactive)
        call zero_forces(grad, log_mb(iact), nat)
        
        !
        ! symmetrize
        !
        if (oskel) call sym_grad_symmetrize(geom, grad)

        !
        ! Print Orbital Energies
        !
        if (ovecs) then
          write(luout,*)
          write(luout,*) '      Alpha Orbital Energies '
          do i=1,calc%bas%nao
            write(luout,'(I6,4X,F8.4,4XD18.8)') 
     $        i,wfn%focc(i,1),wfn%emo(i,1)
          enddo
          write(luout,*)
          write(luout,*)

          if (nspin4.gt.1) then
            write(luout,*)
            write(luout,*) '      Beta Orbital Energies '
            do i=1,calc%bas%nao
              write(luout,'(I6,4X,F8.4,4XD18.8)') 
     $          i,wfn%focc(i,2),wfn%emo(i,2)
            enddo
            write(luout,*)
            write(luout,*)
          endif
        endif

        !
        ! Print Forces
        !
        if (ga_nodeid().eq.0 .and.
     $      util_print('forces', print_low)) then
 1000     format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',
     $           15X,'coordinates',24X,'gradient',/,6X,
     $           2(1X,(3(10X,A1))))
 2000     format(1X,I4,1X,A4,2(1X,3(1X,F10.6)))
          write(luout,1000) theory,'x','y','z','x','y','z'
          call util_flush(luout)
          do i=1,nat
            call util_set_ma_char_arrays(16, byte_mb(itags+16*(i-1)),
     $                                   aname)
            write(luout,2000) i,aname,(coords(j,i),j=1,3),
     $                                  (grad(j,i),j=1,3)
          enddo
          call util_flush(luout)
        endif

        !call xtb_restart_write(wfn%density, int(calc%bas%nao))
        call xtb_restart_write(mol, calc, wfn)
        status = rtdb_put(rtdb,'xtb:restart',mt_log,1,.true.)

        dipole(1:3) = sum(wfn%dpat(:,:,1),dim=2)
        status = rtdb_put(rtdb,'task:dipole',mt_dbl,3,dipole(1))
        status = rtdb_parallel(.true.)
      endif

      xtb_energy_gradient = .not.ctx%failed()

      !
      ! Deallocate all xTB objects
      !
      deallocate(ctx)
      deallocate(wfn)
      deallocate(calc)
      deallocate(mol)
      deallocate(param)
     
      if (.not.ma_chop_stack(lcharge))
     $  call errquit('xtb_energy: failed to pop stack',0,MA_ERR)
      if (.not.geom_destroy(geom))
     $  call errquit('xtb_energy: could not destroy geom',0,GEOM_ERR)
      deallocate(coords,atnum4,stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to deallocate',0,MA_ERR)

      if (ga_nodeid().eq.0) then
        status = rtdb_parallel(.false.)
        if (.not.rtdb_put(rtdb,'xtb:energy',mt_dbl,1,energy))
     $    call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'xtb:converged',mt_log,1,.true.))
     $    call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)

        if (.not.rtdb_put(rtdb,'xtb:gradient',mt_dbl,3*nat,grad(1,1)))
     $    call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)
        status = rtdb_parallel(.true.)
      endif

      deallocate(grad,stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to deallocate',0,MA_ERR)


      call util_print_pop()

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif
#endif
      end
