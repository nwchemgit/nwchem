      logical function xtb_energy(rtdb)
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
      use tblite_scf_info, only : scf_info,atom_resolved,shell_resolved
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
#include "xtb.fh"

      integer rtdb
      integer geom, i, iat
      logical converged, status, oprint, ovecs, oparam
      double precision walltime, cputime

      integer istat
      integer icharge,lcharge,itags,ltags
      integer ncenters,nat,dummybq,atn

      integer iuhf8

      character(len=16) aname,elem
      character(len=2) symb
      character(len=8) theory
      character(len=255) temp

      integer*4 :: iuhf4, one4, zero4, nspin4, verb4
      integer*4,allocatable :: atnum(:)

      double precision energy
      double precision, allocatable :: coords(:,:)

      double precision, allocatable :: rcov(:), trans(:,:), cn(:)

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      type(context_type) :: ctx
      type(param_record) :: param
      type(scf_info) :: info

      parameter(one4=1, zero4=0)
#ifdef USE_OPENMP
      integer maxthreads
      maxthreads = omp_get_max_threads()
      call util_blas_set_num_threads(maxthreads)
#endif

      xtb_energy = .false.
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push()
      call util_print_rtdb_load(rtdb,'xtb')
      oprint = util_print('information', print_low) .and. 
     $         ga_nodeid().eq.0
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
      if (.not.rtdb_put(rtdb, 'xtb:energy', mt_dbl, 1, 0d0))
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

      verb4 = verb
      iuhf4 = uhf
      nspin4 = 1
C      if (iuhf4.gt.0) nspin4 = 2

      !
      ! Allocate space for atomic coordinates
      !
      allocate(coords(3,ncenters),stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to alloc coordinates',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,ncenters,'charges',lcharge,icharge))
     $  call errquit('xtb_energy: failed to alloc charges',0,MA_ERR)
      if (.not.ma_push_get(mt_byte,ncenters*16,'tags',ltags,itags))
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

      allocate(atnum(ncenters),stat=istat)
      if (istat.ne.0)
     $  call errquit('xtb_energy: could not allocate numbers',0,MA_ERR)

      do i=0,ncenters-1
        atnum(i+1) = nint(dbl_mb(icharge+i))
      enddo

      energy = 0.0d0

      if (ga_nodeid().eq.0) then
        call new(mol,atnum,coords,charge=rcharge,uhf=iuhf4)
        if (gfn2) then
          call new_gfn2_calculator(calc, mol)
          call export_gfn2_param(param)
        elseif (gfn1) then
          call new_gfn1_calculator(calc, mol)
          call export_gfn1_param(param)
        endif
        call new_wavefunction(wfn,mol%nat,calc%bas%nsh,calc%bas%nao,
     $                      nspin4,kt)


        ! Compute starting charges 
        if (guess) then
          call xtb_guess(mol,wfn%qat,wfn%qsh,param,mol%nat,
     $                   calc%bas%nsh,calc%bas%nsh_at,atnum)
        endif


        if (oparam) then

          info = calc%coulomb%es3%variable_info()
          select case(info%charge)
          case (atom_resolved)
            temp = 'atom-resolved'
          case (shell_resolved)
            temp = 'shell-resolved'
          end select


 5000     format(/,10X,'No. of atomic orbitals: ',I10,
     $         /,10X,'No. of shells         : ',I10,
     $         /,10X,'Electronic temp.      : ',F10.6,
     $         /,10X,'Accuracy              : ',F10.6,
     $         /,10X,'    Integral cutoff   : ',F10.6,
     $         /,10X,'    E. convergence    : ',F10.6,
     $         /,10X,'    W.F. convergence  : ',F10.6,
     $         /,10X,'Broyden damping       : ',F10.6,
     $         /,10X,'Couomb 3rd order      : ',A,
     $         /,10X,'H0-scaling            : ',3(F10.6),
     $         /,10X,'s8                    : ',F10.6,
     $         /,10X,'a1                    : ',F10.6,
     $         /,10X,'a2                    : ',F10.6,
     $         /,10X,'s9                    : ',F10.6,
     $         /,10X,'kExp                  : ',2(F10.6),
     $         /,10X,'alpha                 : ',F10.6,
     $         /,10X,'a3                    : ',F10.6,
     $         /,10X,'a5                    : ',F10.6,
     $         /,10X,'cn-shift              : ',F10.6,
     $         /,10X,'cn-exp                : ',F10.6,
     $         /,10X,'max-rad               : ',F10.6)


          write(luout,5000) calc%bas%nao,calc%bas%nsh,kt/kb_au,acc,
     $      calc%bas%intcut, 1d-6*acc, 2d-5*acc,calc%mixer_damping,
     $      temp(1:len(trim(temp))),
     $      (param%hamiltonian%ksh(i,i),i=0,param%hamiltonian%lmax),
     $      param%dispersion%s8,param%dispersion%a1,param%dispersion%a2,
     $      param%dispersion%s9,param%repulsion%kexp,
     $      param%repulsion%klight,param%charge%gexp,
     $      param%multipole%dmp3,param%multipole%dmp5,
     $      param%multipole%shift,param%multipole%kexp,
     $      param%multipole%rmax
          call util_flush(luout)
        endif


        call xtb_singlepoint(ctx,mol,calc,wfn,acc,energy,
     $                       verbosity=verb4)

        write(luout,*)
        write(luout,'(A,F15.8)') '          Total xTB energy = ',energy
        write(luout,*)

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

        call xtb_deallocate_wavefunction(wfn)
        call xtb_deallocate_structure(mol)
      endif

     
      if (.not.ma_chop_stack(lcharge))
     $  call errquit('xtb_energy: failed to pop stack',0,MA_ERR)
      if (.not.geom_destroy(geom))
     $  call errquit('xtb_energy: could not destroy geom',0,GEOM_ERR)
      deallocate(coords,atnum,stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to deallocate',0,MA_ERR)


      if (.not.rtdb_put(rtdb,'xtb:energy',mt_dbl,1,energy))
     $  call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'xtb:converged',mt_log,1,.true.))
     $  call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)

      xtb_energy = .not.ctx%failed()
      call util_print_pop()

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

#endif
      end
