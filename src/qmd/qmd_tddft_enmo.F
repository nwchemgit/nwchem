c
c     qmd_tddft_enmo
c
      subroutine qmd_tddft_enmo(rtdb,geom,state,nstates,nroots,ipol,
     $                          nmo,nocc,en,g_v)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb               ! rtdb
      integer geom
      integer state
      integer nstates
      integer nroots
      integer ipol
      integer nmo(2)
      integer nocc(2)
      integer g_v(2)
      double precision en(nstates)
c
      character*30 pname
c
c     preliminaries
      pname = "qmd_tddft_enmo: "
c
      call qmd_tddft_energy(rtdb,state,nstates,nroots,en)
      call qmd_tddft_movecs(rtdb,geom,ipol,nmo,nocc,g_v)
c
      call ga_sync()
c
      return
      end
c
c     qmd_tddft_energy:
c
      subroutine qmd_tddft_energy(rtdb,state,nstates,nroots,en)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb               ! rtdb
      integer state
      integer nstates
      integer nroots
      double precision en(nstates)
c
      integer ien
      integer l_enall, k_enall
c
      character*30 pname
c
c     preliminaries
      pname = "qmd_tddft_energy: "
c
      if (state.eq.1) call qmd_energy(rtdb,'tddft')
c
      if (.not.ma_push_get(mt_dbl,nroots,'enall',l_enall,k_enall))
     $   call errquit(pname//'ma_push_get enall',nroots,MA_ERR)
c
      if(.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,en(1)))
     $   call errquit(pname//'dft energy',1,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'tddft:energy-all',mt_dbl,nroots,
     $                 dbl_mb(k_enall)))
     $   call errquit(pname//'tddft energy',1,RTDB_ERR)
      do ien=2,nstates
        en(ien)=dbl_mb(k_enall+ien-2)+en(1)
      end do
c
      if(.not.ma_pop_stack(l_enall))
     + call errquit(pname//'ma_chop_stack enall',0,MA_ERR)
c
      return
      end
c
c     qmd_tddft_movecs
c
      subroutine qmd_tddft_movecs(rtdb,geom,ipol,nmo,nocc,g_v)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
c
#include "qmd_common.fh"
c
      integer rtdb             ! RTDB handle
      integer geom             ! geometry handle
      integer g_v(2)           ! movecs handle
      integer ipol
      integer nmo(2)
      integer nocc(2)
c
c     local variables
      integer nbf              ! number of basis functions
      integer ipol_verify
      integer nmo_verify(2)
      integer nocc_verify(2)
      integer basis
      integer ispin
      integer l_vals, k_vals
      integer l_occ, k_occ
      character*255 movecs          ! MO vector file name
      character*32 theory
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read_header
      external movecs_read_header
      logical  movecs_read
      external movecs_read
c
      character*30 pname
c
c     preliminaries
      pname='qmd_tddft_movecs: '
c
      if (.not.bas_create(basis,'ao basis'))
     &     call errquit(pname//'failed to create',0,BASIS_ERR)
c
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &     call errquit(pname//'failed to load',0,BASIS_ERR)
c
      if (.not. rtdb_cget(rtdb,'dft:output vectors',1,movecs))
     $     call errquit(pname//'DFT MO vectors not defined',0,
     &    RTDB_ERR)
c
      if (.not.movecs_read_header(movecs,title,basisname,
     $  scftype,nbf,ipol_verify,nmo_verify,2)) call errquit
     $  (pname//'failed to read MO file header',0, DISK_ERR)
c
      if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc_verify))
     &  call errquit(pname//'failed to get noc',0,RTDB_ERR)
c
      if (nmo(1).gt.0) then
        do ispin=1,ipol
         if (nmo(ispin).ne.nmo_verify(ispin))
     &     call errquit(pname//'#MOs changed',0,RTDB_ERR)
         if (nocc(ispin).ne.nocc_verify(ispin))
     &     call errquit(pname//'#electrons changed',0,RTDB_ERR)
        end do
      end if
      nmo=nmo_verify
      nocc=nocc_verify
c
c     Allocate MO energies array
      if(.not.ma_push_get(mt_dbl,nbf,'vals',l_vals,k_vals))
     $   call errquit(pname//'could not allocate evals',nbf,MA_ERR)
c
c     Allocate MO occup array
      if(.not.ma_push_get(mt_dbl,nbf,'occ',l_occ,k_occ))
     $   call errquit(pname//'could not allocate occup',nbf,MA_ERR)
c
      do ispin=1,ipol
        if(.not.movecs_read(movecs,ispin,dbl_mb(k_occ),dbl_mb(k_vals),
     $                      g_v(ispin)))
     $     call errquit(pname//'failed to read MO vectors',0,DISK_ERR)
      end do
c
      if (.not.ma_chop_stack(l_vals)) 
     $  call errquit(pname//'failed chopping MA stack',0, MA_ERR)
      if (.not. ma_verify_allocator_stuff())
     $  call errquit(pname//'ma_verify_allocator_stuff',0,MA_ERR)
      if (.not.bas_destroy(basis))
     $  call errquit(pname//'failed to destroy basis',basis,BASIS_ERR)
c
      return
      end
