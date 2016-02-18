c
c     qmd_tdks_enmo
c
      subroutine qmd_tdks_enmo(rtdb,geom,state,nstates,ipol,nmo,nocc,
     $                         swap_occ,swap_vir,en,g_v,ks_spin)
c
      implicit none
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer rtdb               ! rtdb
      integer geom
      integer state
      integer nstates
      integer ipol
      integer nmo(2)
      integer nocc(2)
      integer swap_occ(nstates)
      integer swap_vir(nstates)
      integer g_v(2)
      integer ks_spin
      double precision en(nstates)
c
      integer swap(2)
c
      character*30 pname
c
c     preliminaries
      pname = "qmd_tdks_enmo: "
c
c     force atomic guess
      if(.not.rtdb_cput(rtdb,'dft:input vectors',1,'atomic'))
     &   call errquit(pname//'failed to set atomic',0,RTDB_ERR)
c
c     turn off noscf
      if(.not.rtdb_put(rtdb,'dft:noscf',mt_log,1,.false.))
     $   call errquit(pname//'noscf failed',0,RTDB_ERR)
c
c     calculate ground state energy
      call qmd_energy(rtdb,'dft')
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,en(1)))
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
      swap_occ(1)=0
      swap_vir(1)=0
c
c     get movecs
      call qmd_tdks_movecs(rtdb,geom,ipol,nmo,nocc,nstates,swap_occ,
     $                     swap_vir,en,g_v,ks_spin)
c
c     if we are not on the ground state, set up movecs for gradient
      if (state.ne.1) then
        call qmd_tdks_swap(rtdb,state,nstates,ks_spin,swap_occ,
     $                     swap_vir)
      end if
c
      return
      end
c
c     qmd_tdks_movecs
c
      subroutine qmd_tdks_movecs(rtdb,geom,ipol,nmo,nocc,nstates,
     $                           swap_occ,swap_vir,en,g_v,ks_spin)
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
      integer nstates
      integer swap_occ(nstates)
      integer swap_vir(nstates)
      integer ks_spin
      double precision en(nstates)
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
      pname='qmd_tdks_movecs: '
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
      if (.not.ma_push_get(mt_dbl,nbf,'vals',l_vals,k_vals))
     $     call errquit(pname//'could not allocate evals',nbf,MA_ERR)
c
c     Allocate MO occup array
      if (.not. ma_push_get(mt_dbl,nbf,'occ',l_occ,k_occ))
     $     call errquit(pname//'could not allocate occup',nbf,MA_ERR)
c
      do ispin=1,ipol
        if(.not.movecs_read(movecs,ispin,dbl_mb(k_occ),dbl_mb(k_vals),
     $                      g_v(ispin)))
     $      call errquit(pname//'failed to read MO vectors',0,DISK_ERR)
c
        if (ks_spin.eq.ispin) then
          call qmd_tdks_energy(nocc(ks_spin),nmo(ks_spin),nstates,
     $                         swap_occ,swap_vir,en,dbl_mb(k_vals))
        end if
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
c
c     qmd_tdks_energy
c
      subroutine qmd_tdks_energy(nocc,nmo,nstates,swap_occ,
     $                           swap_vir,en,movals)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
#include "qmd_common.fh"
c
      integer nmo
      integer nocc
      integer nstates
      integer swap_occ(nstates)
      integer swap_vir(nstates)
      double precision en(nstates)
      double precision movals(nmo)
c
c     local variables
      integer l_exe, k_exe
      integer l_exo, k_exo
      integer l_exv, k_exv
      integer maxsize
      integer iostart
      integer ivend
      integer i, j, n, m
      integer isort
      double precision esort
      double precision exen
c
      character*30 pname
c
c     preliminaries
      pname='qmd_tdks_energy: '
c
      maxsize=0
      iostart=max(nocc-nstates,1)
      ivend=min(nocc+nstates+1,nmo)
      maxsize=(nocc-iostart+1)*(ivend-nocc)
c
      if (.not.ma_push_get(mt_dbl,maxsize,'exe',l_exe,k_exe))
     $     call errquit(pname//'could not allocate exen',maxsize,MA_ERR)
      if (.not.ma_push_get(mt_int,maxsize,'exo',l_exo,k_exo))
     $     call errquit(pname//'could not allocate exen',maxsize,MA_ERR)
      if (.not.ma_push_get(mt_int,maxsize,'exv',l_exv,k_exv))
     $     call errquit(pname//'could not allocate exen',maxsize,MA_ERR)
c
c     calculate excitation energies
      n=0
      do i=iostart,nocc
        do j=nocc+1,ivend
          exen=movals(j)-movals(i)
          dbl_mb(k_exe+n)=exen
          int_mb(k_exo+n)=i
          int_mb(k_exv+n)=j
          n=n+1
        end do
      end do
c
c     sort by excitation energy
      do i=1,nstates-1
        do j=i+1,n
          if (dbl_mb(k_exe+j-1).lt.dbl_mb(k_exe+i-1)) then
c           sort energy
            esort=dbl_mb(k_exe+j-1)
            dbl_mb(k_exe+j-1)=dbl_mb(k_exe+i-1)
            dbl_mb(k_exe+i-1)=esort
c           sort occupied
            isort=int_mb(k_exo+j-1)
            int_mb(k_exo+j-1)=int_mb(k_exo+i-1)
            int_mb(k_exo+i-1)=isort
c           sort virtual
            isort=int_mb(k_exv+j-1)
            int_mb(k_exv+j-1)=int_mb(k_exv+i-1)
            int_mb(k_exv+i-1)=isort
          end if
        end do
        en(i+1)=dbl_mb(k_exe+i-1)+en(1)
        swap_occ(i+1)=int_mb(k_exo+i-1)
        swap_vir(i+1)=int_mb(k_exv+i-1)
      end do
c
      if (.not.ma_chop_stack(l_exe)) 
     $  call errquit(pname//'failed chopping MA stack',0, MA_ERR)
c
      return
      end
c
c     qmd_tdks_swap
c
      subroutine qmd_tdks_swap(rtdb,state,nstates,ks_spin,swap_occ,
     $                         swap_vir)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
      integer state
      integer nstates
      integer ks_spin
      integer swap_occ(nstates)
      integer swap_vir(nstates)
c
      integer swap(2)
c
      character*30 pname
c
c     preliminaries
      pname='qmd_tdks_swap: '
c
c     set up swapped movecs for gradient
      swap(1)=swap_occ(state)
      swap(2)=swap_vir(state)
      if (ks_spin.eq.1) then
        if(.not. rtdb_put(rtdb,'dft:swap alpha',mt_int,2,swap))
     $     call errquit(pname//'swap failed',1,RTDB_ERR)
      else
        if(.not. rtdb_put(rtdb,'dft:swap beta',mt_int,2,swap))
     $     call errquit(pname//'swap failed',2,RTDB_ERR)
      end if
c
c     turn on noscf
      if(.not.rtdb_put(rtdb,'dft:noscf',mt_log,1,.true.))
     $   call errquit(pname//'noscf failed',0,RTDB_ERR)
c
c     unset convergence flag
      if(.not.rtdb_put(rtdb,'dft:converged',mt_log,1,.false.))
     $   call errquit(pname//'converged failed',0,RTDB_ERR)
c
      return
      end
