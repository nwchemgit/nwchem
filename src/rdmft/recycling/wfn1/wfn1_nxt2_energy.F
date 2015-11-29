C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief The main driver for the energy optimization in the next step
C> method
C>
C> This driver routine loads the settings, the Hartree-Fock
C> wavefunction, arranges the memory, drives the energy optimization and
C> subsequently the excited state optimization. The methods utilized are
C> the ones described in the file wfn1_next_step.dox.
C>
      logical function wfn1_nxt2_energy(rtdb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "wfn1_nxt.fh"
c
      integer rtdb !< [Input] The runtime database handle
c
c     Local variables
c
      integer nea !< The number of alpha electrons
      integer neb !< The number of beta electrons
c     integer nbf !< The number of basis functions (see cscf.fh)
      integer nst !< The number of states to calculate (including 
                  !< the ground state)
      integer ist !< The state number
c
      integer l_stat,  k_stat  !< The state function
      integer l_state, k_state !< The excited-state state function
      integer l_grad,  k_grad  !< The state function gradient
c
      character*255 bas_vecs     !< The basis set name on file
      character*255 title_vecs   !< The calculation title on file
      character*20  file_scftype !< The energy expression type on file
      integer k_ehfa, k_ehfb     !< Orbital energy pointers
      integer nbf_vecs           !< The number of basis functions on
                                 !< file
      integer nmo_vecs(2)        !< The number of vectors on file
      integer nsets              !< The number of sets (spin channels?)
      integer ii                 !< Counter
      integer ij                 !< Counter
      integer it                 !< Scratch
      logical status
c
c     Swaps are stored in pairs (occupied1,virtual1), ...,
c     (occupiedn,virtualn)
c
      integer nswapa           !< The number of alpha swaps
      integer nswapb           !< The number of beta swaps
      integer l_swapa, k_swapa !< The alpha orbitals to swap
      integer l_swapb, k_swapb !< The beta orbitals to swap
      integer iswptyp          !< The type of swap lists
      integer iannihilate      !< Orbital to remove electron from
      integer icreate          !< Orbital to put electron into
c
      double precision nuclear_charge !< The total nuclear charge
      double precision etot           !< The total energy
c
      integer lkey        !< The length of the key
      character*32 keyfmt !< The format of key
      character*8  key    !< The key number as a string
      character*32 keya   !< The alpha electron key
      character*32 keyb   !< The beta electron key
c
      character*8  wfn1_init_mode !< The type of starting wavefunction
c
      logical  movecs_read_header
      external movecs_read_header
      logical  movecs_read
      external movecs_read
      logical  int_normalize
      external int_normalize
c
c     Code
c
      wfn1_nxt2_energy = .false.
c
c
c     Figure the geometry out
c     =======================
c     - load the geometry
c     - get the total nuclear charge
c     - get the total molecular charge
c     - compute the number of electrons
c
c
      if (.not. geom_create(geom,'geometry'))
     &  call errquit("wfn1_nxt2_energy: geom_create failed",0,GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     &  call errquit("wfn1_nxt2_energy: geom load failed",0,RTDB_ERR)
      if (.not. geom_nuc_charge(geom,nuclear_charge))
     &  call errquit("wfn1_nxt2_energy: geometry nuclear charge failed",
     &               0,GEOM_ERR)
      if (.not. rtdb_get(rtdb,"charge",MT_DBL,1,charge))
     &  charge = 0.0d0
c
      if (.not. geom_print(geom))
     &  call errquit("wfn1_nxt2_energy: geometry print failed",
     &               0,GEOM_ERR)
c
      nelec = int(nuclear_charge - charge)
      if (nelec.lt.0)
     &   call errquit("wfn1_nxt2_energy: negative number of electrons",
     &                nelec,INPUT_ERR)

c
c
c     Figure the type of wavefunction out
c     ===================================
c     - Is an unrestricted wavefunction or a restricted one
c     - How many electrons are there in each spin channel
c
c
      if (.not. rtdb_cget(rtdb,'scf:scftype',1,scftype)) scftype = 'RHF'
      if (scftype.eq.'RHF'.or.scftype.eq.'ROHF') then
        if (.not. rtdb_get(rtdb,'scf:nopen',MT_INT,1,nopen)) nopen = 0
        if (nopen.ne.0) then
          scftype = 'ROHF'
        else
          scftype = 'RHF'
        endif
        if (nopen.gt.nelec)
     &     call errquit("wfn1_nxt2_energy: more unpaired electrons than"
     &                //" all electrons",nopen,INPUT_ERR)
        if (mod(nelec-nopen,2).ne.0)
     &     call errquit("wfn1_nxt2_energy: inconsistent no. electrons "
     &                //"and spin",nopen,INPUT_ERR)
        nclosed = (nelec-nopen)/2
        nea = nclosed + nopen
        neb = nclosed
      else if (scftype.eq.'UHF') then
        if (.not. rtdb_get(rtdb,'scf:nalpha',MT_INT,1,nea))
     &     call errquit("wfn1_nxt2_energy: no. alpha electrons not "
     &                  //"found",0,INPUT_ERR)
        if (.not. rtdb_get(rtdb,'scf:nbeta',MT_INT,1,neb))
     &     call errquit("wfn1_nxt2_energy: no. beta electrons not "
     &                  //"found",0,INPUT_ERR)
      else
        call errquit("wfn1_nxt2_energy: unexpect scftype: "//scftype,
     &               0,INPUT_ERR)
      endif
c
c
c     Figure the basis set out
c     ========================
c     - Load the basis set
c     - Find the number of basis functions
c     - Normalize the basis functions
c     - Initialize the off-set common block
c
c
      if (.not. bas_create(basis,'ao basis'))
     &  call errquit("wfn1_nxt2_energy: bas_create failed",0,BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &  call errquit("wfn1_nxt2_energy: basis load failed",0,RTDB_ERR)
      if (.not. int_normalize(rtdb,basis))
     &  call errquit("wfn1_nxt2_energy: normalization failed",0,INT_ERR)
      if (.not. bas_numbf(basis,nbf))
     &  call errquit("wfn1_nxt2_energy: no basis functions",0,BASIS_ERR)
      if (.not. bas_name(basis, trans, trans))
     &  call errquit('wfn1_nxt2_energy: bad basis',0, BASIS_ERR)
      if (.not. bas_print(basis))
     &  call errquit("wfn1_nxt2_energy: no basis",0,BASIS_ERR)
      call wfn1_nxt_init(nbf)
c
c
c     Load the Hartree-Fock wavefunctions
c     ===================================
c     - Create Global Arrays to hold the wavefunction
c     - Allocate memory
c     - Load the wavefunctions
c     - Check the data obtained
c
c
      if (.not.ga_create(MT_DBL,nbf,nbf,'alpha mos',-1,-1,g_movecs(1)))
     &   call errquit("wfn1_nxt2_energy: could not allocate alpha mos",
     &                0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,'beta mos',-1,-1,g_movecs(2)))
     &   call errquit("wfn1_nxt2_energy: could not allocate beta mos",
     &                0,GA_ERR)
      call ga_zero(g_movecs(1))
      call ga_zero(g_movecs(2))
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo evals',
     &   l_eval, k_eval)) call errquit
     &   ('wfn1_nxt2_energy: insufficient memory?',
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
      k_ehfa = k_eval
      k_ehfb = k_eval + nbf
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo occ',
     &   l_occ, k_occ)) call errquit
     &   ('wfn1_nxt2_energy: insufficient memory?',
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
c
      if (.not. ma_push_get(MT_INT,2*nbf, 'scf_init: mo irs',
     &   l_irs, k_irs)) call errquit
     &   ('scf_init: insufficient memory?',
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
c
      call ifill(2*nbf, 1, int_mb(k_irs), 1) ! In case not adapting
c
c     Find the MO-vectors
c
      if (.not. rtdb_cget(rtdb,'scf:input vectors',1,movecs_in))
     &   call errquit("wfn1_nxt2_energy: MO-vectors from where?",
     &                0,INPUT_ERR)
      if (.not. rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     &   movecs_out = movecs_in
      call util_file_name_resolve(movecs_in, .false.)
      call util_file_name_resolve(movecs_out, .false.)
c
      status = movecs_read_header(movecs_in, title_vecs, bas_vecs,
     &         file_scftype, nbf_vecs, nsets, nmo_vecs, 2)
      if (status) then
c
c        Check that:
c        1) nbf from file = current nbf
c        2) name of basis from file = current basis name
c        3) nmo from file = current nmo (in case of change of linear
c                           dep)
c
         if ((nbf.eq.nbf_vecs) .and. (bas_vecs.eq.trans) .and.
     &       (nbf.eq.nmo_vecs(1))) then
            if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     &                            dbl_mb(k_ehfa), g_movecs))
     &         call errquit('wfn1_nxt2_energy: movecs_read failed ',
     &                      0, DISK_ERR)
c
            if (file_scftype .eq. 'UHF' .or.
     &          file_scftype .eq. 'MCSCF') then
               if (.not. movecs_read(movecs_in, 2, dbl_mb(k_occ+nbf),
     &                               dbl_mb(k_ehfb),
     &                               g_movecs(2))) then
                  call ga_copy(g_movecs(1), g_movecs(2))
                  call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
                  call dcopy(nbf,dbl_mb(k_ehfa),1,dbl_mb(k_ehfb),1)
               endif
            else
               call ga_copy(g_movecs(1), g_movecs(2))
               do ii = 0, nbf-1
                 dbl_mb(k_occ+ii) = dbl_mb(k_occ+ii)/2.0d0
               enddo
               call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
               call dcopy(nbf,dbl_mb(k_ehfa),1,dbl_mb(k_ehfb),1)
            endif
         else
            call errquit("wfn1_nxt2_energy: no suitable vectors found",
     &                   0, DISK_ERR)
         endif
c        if (.not.rtdb_get(rtdb,'scf:orb:ehfa',MT_DBL,
c    $                     nbf,dbl_mb(k_ehfa))) then
c          write(*,*)'WARNING: could not load scf:orb:ehfa'
c        endif
c        if (.not.rtdb_get(rtdb,'scf:orb:ehfb',MT_DBL,
c    $                     nbf,dbl_mb(k_ehfb))) then
c          write(*,*)'WARNING: could not load scf:orb:ehfb'
c        endif
         if (.not.rtdb_get(rtdb,'mcscf:orb:occa',MT_DBL,
     $                     nbf,dbl_mb(k_occ))) then
c          write(*,*)'WARNING: could not load mcscf:orb:occa'
           call dfill(nbf,0.0d0,dbl_mb(k_occ),1)
           call dfill(nea,1.0d0,dbl_mb(k_occ),1)
         endif
         if (.not.rtdb_get(rtdb,'mcscf:orb:occb',MT_DBL,
     $                     nbf,dbl_mb(k_occ+nbf))) then
c          write(*,*)'WARNING: could not load mcscf:orb:occb'
           call dfill(nbf,0.0d0,dbl_mb(k_occ+nbf),1)
           call dfill(neb,1.0d0,dbl_mb(k_occ+nbf),1)
         endif
      else
         call errquit("wfn1_nxt2_energy: no vectors found",0, DISK_ERR)
      endif
c
c     Load energy expression parameters
c     ---------------------------------
c
c     Read the parameters that define the energy expression from the
c     runtime database. These include:
c     - x_pow: the power of the density matrix in the exchange term
c
      if (.not.rtdb_get(rtdb,'wfn1:power',MT_DBL,1,wfn1_nxt_x_pow))
     &    wfn1_nxt_x_pow = 1.0d0 ! use Hartree-Fock exchange
c
c     - ad: if .true. use automatic differentiation
c
      if (.not.rtdb_get(rtdb,'wfn1:ad',MT_LOG,1,wfn1_nxt_ad))
     &    wfn1_nxt_ad = .false. ! do not use automatic differentiation
c
c     - mode: if "even" distribute electrons evenly across the natural
c             orbitals, if "fermi" use Fermi-smearing to distribute
c             the electrons, if "hf" generate the Hartree-Fock 
c             wavefunction.
c
      if (.not.rtdb_cget(rtdb,'wfn1:mode',1,wfn1_init_mode))
     &    wfn1_init_mode = "even" ! distribute electrons evenly
c
c     - ncore: the number of core occupation functions (i.e. straight
c              unit vectors at the low energy end).
c
      if (.not.rtdb_get(rtdb,'wfn1:ncore',MT_INT,2,wfn1_nxt_ncore)) then
          wfn1_nxt_ncore(1) = 0 ! no core orbitals
          wfn1_nxt_ncore(2) = 0 ! no core orbitals
      endif
c
c     - nvirt: the number of virtual occupation functions (i.e. straight
c              unit vectors at the high energy end).
c
      if (.not.rtdb_get(rtdb,'wfn1:nvirt',MT_INT,2,wfn1_nxt_nvirt)) then
          wfn1_nxt_nvirt(1) = 0 ! no core orbitals
          wfn1_nxt_nvirt(2) = 0 ! no core orbitals
      endif
c
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"state",
     &                     l_stat,k_stat))
     &   call errquit("wfn1_nxt2_energy: could not allocate state",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
      k_stat = k_stat - 1
      call wfn1_nxt2_init_state(wfn1_init_mode,wfn1_nxt_ncore,
     &     wfn1_nxt_nvirt,g_movecs,nbf,nea,neb,
     &     dbl_mb(k_ehfa),dbl_mb(k_ehfb),
     &     dbl_mb(k_stat+wfn1_nxt_aocc()),
     &     dbl_mb(k_stat+wfn1_nxt_anat()),
     &     dbl_mb(k_stat+wfn1_nxt_aoccl()),
     &     dbl_mb(k_stat+wfn1_nxt_anatl()),
     &     dbl_mb(k_stat+wfn1_nxt_aconf()),
     &     dbl_mb(k_stat+wfn1_nxt_bocc()),
     &     dbl_mb(k_stat+wfn1_nxt_bnat()),
     &     dbl_mb(k_stat+wfn1_nxt_boccl()),
     &     dbl_mb(k_stat+wfn1_nxt_bnatl()),
     &     dbl_mb(k_stat+wfn1_nxt_bconf()))
      k_stat = k_stat + 1
cDEBUG
      write(*,*)'*** Initial wavefunction'
      call wfn1_nxt2_print_state(dbl_mb(k_stat),nbf)
cDEBUG
c
      call wfn1_nxt2_min_energy(rtdb,geom,basis,nbf,nea,neb,
     &     dbl_mb(k_stat),etot)
c
c     Store the ground state total energy to keep the task layer happy
c
      if (etot.lt.0.0d0) then
        if (.not.rtdb_put(rtdb,'wfn1-nxt2:energy',MT_DBL,1,etot))
     &    call errquit("wfn1_nxt2_energy: could not store energy",
     &                 0,RTDB_ERR)
      endif
c
c     See if we need to calculate any excited states (i.e. nst > 1)
c
      if (rtdb_get(rtdb,"wfn1:nstate",MT_INT,1,nst)) then
      else if (rtdb_get(rtdb,"wfn1-nxt:state",MT_INT,1,nst)) then
      else if (rtdb_get(rtdb,"wfn1-nxt2:state",MT_INT,1,nst)) then
      else
        nst = 1
      endif
cDEBUG
      write(*,*)"number of states requested: ",nst
cDEBUG
c
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"state-exc",
     &                     l_state,k_state))
     &   call errquit("wfn1_nxt2_energy: could not allocate state-exc",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
      do ist = 2, nst
c
c       Read the orbital swaps from the runtime data base
c
        if (ist.lt.10) then
          keyfmt = "(i1)"
        else if (ist.lt.100) then
          keyfmt = "(i2)"
        else
          call errquit("wfn1_nxt2_energy: too many states",nst,UERR)
        endif
        write(key,keyfmt)ist
        lkey=inp_strlen(key)
        keya="wfn1-nxt2:"//key(1:lkey)//":a"
        keyb="wfn1-nxt2:"//key(1:lkey)//":b"
        lkey=inp_strlen(keya)
cDEBUG
        write(*,*)"state: ",ist," keya: ",keya(1:lkey)
        write(*,*)"state: ",ist," keyb: ",keyb(1:lkey)
cDEBUG
        if (rtdb_ma_get(rtdb,keya(1:lkey),iswptyp,nswapa,l_swapa)) then
          if (.not.ma_get_index(l_swapa,k_swapa)) then
            call errquit("wfn1_nxt2_energy: could not retrieve index",
     &                   l_swapa,MA_ERR)
          endif
        else
          k_swapa = 0
          nswapa  = 0
        endif
        if (iswptyp.ne.MT_INT) then
          call errquit("wfn1_nxt2_energy: wrong type for alpha orbital "
     &               //"swaps",iswptyp,RTDB_ERR)
        endif
cDEBUG
        write(*,*)"state: ",ist," nswapa: ",nswapa
        write(*,*)(int_mb(k_swapa+ii),ii=0,nswapa-1)
cDEBUG
        if (rtdb_ma_get(rtdb,keyb(1:lkey),iswptyp,nswapb,l_swapb)) then
          if (.not.ma_get_index(l_swapb,k_swapb)) then
            call errquit("wfn1_nxt2_energy: could not retrieve index",
     &                   l_swapb,MA_ERR)
          endif
        else
          k_swapb = 0
          nswapb  = 0
        endif
        if (iswptyp.ne.MT_INT) then
          call errquit("wfn1_nxt2_energy: wrong type for beta orbital "
     &               //"swaps",iswptyp,RTDB_ERR)
        endif
cDEBUG
        write(*,*)"state: ",ist," nswapb: ",nswapb
        write(*,*)(int_mb(k_swapb+ii),ii=0,nswapb-1)
cDEBUG
c
c       Do some tests
c
        if (mod(nswapa,2).ne.0) then
          call errquit("wfn1_nxt2_energy: swapa not an integer number "
     &               //"of pairs",nswapa,UERR)
        endif
        if (mod(nswapb,2).ne.0) then
          call errquit("wfn1_nxt2_energy: swapb not an integer number "
     &               //"of pairs",nswapb,UERR)
        endif
        do ii = 1, nswapa/2
          if (int_mb(k_swapa+(ii-1)*2).gt.nea) then
            it = int_mb(k_swapa+(ii-1)*2)
            call errquit("wfn1_nxt2_energy: alpha source orbital not "
     &                 //"occupied",it,UERR)
          endif
          if (int_mb(k_swapa+(ii-1)*2+1).le.nea) then
            it = int_mb(k_swapa+(ii-1)*2+1)
            call errquit("wfn1_nxt2_energy: alpha target orbital not "
     &                 //"unoccupied",it,UERR)
          endif
        enddo
        do ii = 1, nswapb/2
          if (int_mb(k_swapb+(ii-1)*2).gt.neb) then
            it = int_mb(k_swapb+(ii-1)*2)
            call errquit("wfn1_nxt2_energy: beta source orbital not "
     &                 //"occupied",it,UERR)
          endif
          if (int_mb(k_swapb+(ii-1)*2+1).le.neb) then
            it = int_mb(k_swapb+(ii-1)*2+1)
            call errquit("wfn1_nxt2_energy: beta target orbital not "
     &                 //"unoccupied",it,UERR)
          endif
        enddo
        do ii = 1, nswapa/2
          do ij = 1, ii-1
            if(int_mb(k_swapa+(ii-1)*2).eq.int_mb(k_swapa+(ij-1)*2))
     &        call errquit("wfn1_nxt2_energy: alpha source repeated",
     &                     int_mb(k_swapa+(ii-1)*2),UERR)
            if(int_mb(k_swapa+(ii-1)*2+1).eq.int_mb(k_swapa+(ij-1)*2+1))
     &        call errquit("wfn1_nxt2_energy: alpha target repeated",
     &                     int_mb(k_swapa+(ii-1)*2+1),UERR)
          enddo
        enddo
        do ii = 1, nswapb/2
          do ij = 1, ii-1
            if(int_mb(k_swapb+(ii-1)*2).eq.int_mb(k_swapb+(ij-1)*2))
     &        call errquit("wfn1_nxt2_energy: beta source repeated",
     &                     int_mb(k_swapb+(ii-1)*2),UERR)
            if(int_mb(k_swapb+(ii-1)*2+1).eq.int_mb(k_swapb+(ij-1)*2+1))
     &        call errquit("wfn1_nxt2_energy: beta target repeated",
     &                     int_mb(k_swapb+(ii-1)*2+1),UERR)
          enddo
        enddo
c
c       Now create an excited state by applying the permutations
c
        call dcopy(wfn1_nxt_size1(),dbl_mb(k_stat),1,dbl_mb(k_state),1)
        do ii = 1, nswapa/2
c         call wfn1_swap_vectors(nea,nbf,int_mb(k_swapa+(ii-1)*2),
c    &         int_mb(k_swapa+(ii-1)*2+1),
c    &         dbl_mb(k_state-1+wfn1_nxt_aocc()))
          iannihilate = int_mb(k_swapa+(ii-1)*2)
          icreate     = int_mb(k_swapa+(ii-1)*2+1)
c         -2 because of offset implied in wfn1_nxt_aconf() and
c         because iannihilate and icreate run from 1 to nbf
          dbl_mb(k_state-2+wfn1_nxt_aconf()+iannihilate) = 0.0d0
          dbl_mb(k_state-2+wfn1_nxt_aconf()+icreate)     = 1.0d0
        enddo
        do ii = 1, nswapb/2
c         call wfn1_swap_vectors(neb,nbf,int_mb(k_swapb+(ii-1)*2),
c    &         int_mb(k_swapb+(ii-1)*2+1),
c    &         dbl_mb(k_state-1+wfn1_nxt_bocc()))
          iannihilate = int_mb(k_swapb+(ii-1)*2)
          icreate     = int_mb(k_swapb+(ii-1)*2+1)
c         -2 because of offset implied in wfn1_nxt_aconf() and
c         because iannihilate and icreate run from 1 to nbf
          dbl_mb(k_state-2+wfn1_nxt_bconf()+iannihilate) = 0.0d0
          dbl_mb(k_state-2+wfn1_nxt_bconf()+icreate)     = 1.0d0
        enddo
c
c       Clean the memory up
c
        if (nswapa.gt.0) then
          if (.not.ma_free_heap(l_swapa)) then
            call errquit("wfn1_nxt2_energy: could not free l_swapa",
     &                   l_swapa,MA_ERR)
          endif
        endif
        if (nswapb.gt.0) then
          if (.not.ma_free_heap(l_swapb)) then
            call errquit("wfn1_nxt2_energy: could not free l_swapb",
     &                   l_swapb,MA_ERR)
          endif
        endif
cDEBUG
        write(*,*)'just before wfn1_nxt2_min_gradient: state'
        call wfn1_nxt2_print_state(dbl_mb(k_state),nbf)
cDEBUG
c
c       Minimize the gradient of this state
c
        call wfn1_nxt2_min_gradient(rtdb,geom,basis,nbf,nea,neb,
     &       dbl_mb(k_state),etot)
c       call wfn1_nxt2_min_energy(rtdb,geom,basis,nbf,nea,neb,
c    &       dbl_mb(k_state),etot)
c       call wfn1_nxt2_min_excit(rtdb,geom,basis,nbf,nea,neb,
c    &       dbl_mb(k_stat),dbl_mb(k_state),etot)
c       
      enddo
 
      wfn1_nxt2_energy = .true.
c
      return
      end
C>
C> @}

