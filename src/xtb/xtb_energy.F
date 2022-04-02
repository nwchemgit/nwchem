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
      use tblite_xtb_gfn2, only : new_gfn2_calculator
      use tblite_xtb_singlepoint, only : xtb_singlepoint
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"

      integer rtdb
      integer geom, i, iat
      logical converged, status, oprint
      double precision walltime, cputime

      integer istat
      integer icharge,lcharge,itags,ltags
      integer ncenters,nat,dummybq,atn

      integer iuhf8
      double precision rcharge

      character(len=16) aname,elem
      character(len=2) symb

      integer*4 :: iuhf4, one4, zero4
      integer*4,allocatable :: atnum(:)

      double precision energy
      double precision, allocatable :: coords(:,:)

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      type(context_type) :: ctx

      double precision kt, acc
      parameter(kt = 300.0d0*3.166808578545117d-06)
      parameter(acc = 0.01d0)
      parameter(one4=1, zero4=0)


      xtb_energy = .false.
      oprint = ga_nodeid().eq.0

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
      if (.not.rtdb_get(rtdb,'xtb:uhf',mt_int,1,iuhf8))then
        iuhf4 = 0
      else
        iuhf4 = iuhf8
      endif

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

      if (ga_nodeid().eq.0) then
        call new(mol,atnum,coords,charge=rcharge,uhf=iuhf4)
        call new_gfn2_calculator(calc, mol)
        call new_wavefunction(wfn,mol%nat,calc%bas%nsh,calc%bas%nao,
     $                      one4,kt)
        call xtb_singlepoint(ctx,mol,calc,wfn,acc,energy,verbosity=one4)
      endif
     
      if (.not.ma_chop_stack(lcharge))
     $  call errquit('xtb_energy: failed to pop stack',0,MA_ERR)
      if (.not.geom_destroy(geom))
     $  call errquit('xtb_energy: could not destroy geom',0,GEOM_ERR)
      deallocate(coords,atnum,stat=istat)
      if (istat.ne.0) 
     $  call errquit('xtb_energy: failed to deallocate',0,MA_ERR)

      xtb_energy = .not.ctx%failed()

#endif
      end
