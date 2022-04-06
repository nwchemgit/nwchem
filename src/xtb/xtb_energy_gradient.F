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
      use tblite_xtb_gfn2, only : new_gfn2_calculator
      use tblite_xtb_gfn1, only : new_gfn1_calculator
      use tblite_xtb_singlepoint, only : xtb_singlepoint
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
#include "sym.fh"
#include "xtb.fh"

      integer rtdb
      integer geom, i, j, iat
      logical converged, status, oprint
      double precision walltime, cputime

      integer istat
      integer icharge,lcharge,itags,ltags,iact,lact
      integer ncenters,nat,dummybq,atn

      integer iuhf8

      logical oskel
      integer nactive

      character(len=16) aname,elem
      character(len=2) symb
      character(len=8) theory

      integer*4 :: iuhf4, one4, zero4, nspin4
      integer*4,allocatable :: atnum(:)

      double precision energy, sigma(3,3)
      double precision, allocatable :: coords(:,:),grad(:,:)

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      type(context_type) :: ctx
      parameter(one4=1, zero4=0)

#ifdef USE_OPENMP
      integer maxthreads
      maxthreads = omp_get_max_threads()
      call util_blas_set_num_threads(maxthreads)
#endif

      xtb_energy_gradient = .false.
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push()
      call util_print_rtdb_load(rtdb,'xtb')

      oprint = util_print('information', print_low) .and. 
     $         ga_nodeid().eq.0

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
      iuhf4 = uhf

      nspin4 = 1
C      if (iuhf4.gt.0) nspin4 = 2

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

      allocate(atnum(ncenters),stat=istat)
      if (istat.ne.0)
     $  call errquit('xtb_energy: could not allocate numbers',0,MA_ERR)

      do i=0,ncenters-1
        atnum(i+1) = nint(dbl_mb(icharge+i))
      enddo

      energy = 0.0d0
      grad(:,:) = 0.0d0
      sigma(:,:) = 0.0d0

      if (ga_nodeid().eq.0) then
        call new(mol,atnum,coords,charge=rcharge,uhf=iuhf4)
        if (gfn2) then
          theory = 'xTB GFN2'
          call new_gfn2_calculator(calc, mol)
        elseif (gfn1) then
          theory = 'xTB GFN1'
          call new_gfn1_calculator(calc, mol)
        endif
        call new_wavefunction(wfn,mol%nat,calc%bas%nsh,calc%bas%nao,
     $                      nspin4,kt)
        call xtb_singlepoint(ctx,mol,calc,wfn,acc,energy,grad,sigma,
     $                       verbosity=one4)


        write(luout,*)
        write(luout,*) '          Total xTB energy = ',energy
        write(luout,*)

c
C     zero force contributions on inactive atoms
        call grad_active_atoms(rtdb, nat, log_mb(iact), nactive)
        call zero_forces(grad, log_mb(iact), nat)
         
C     symmetrize
        if (oskel) call sym_grad_symmetrize(geom, grad)

        if (oprint) then
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

      if (.not.rtdb_put(rtdb,'xtb:gradient',mt_dbl,3*nat,grad))
     $  call errquit('xtb_energy: failed to put in rtdb',0,RTDB_ERR)

      deallocate(grad,stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to deallocate',0,MA_ERR)


      xtb_energy_gradient = .not.ctx%failed()
      call util_print_pop()

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif
#endif
      end
