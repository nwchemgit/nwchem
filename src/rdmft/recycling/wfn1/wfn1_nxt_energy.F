C> \ingroup wfn1_nxt
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
      logical function wfn1_nxt_energy(rtdb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
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
c
      integer l_stat, k_stat !< The state function
      integer l_grad, k_grad !< The state function gradient
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
      logical status
c
      double precision nuclear_charge !< The total nuclear charge
      double precision etot           !< The total energy
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
      wfn1_nxt_energy = .false.
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
     &  call errquit("wfn1_nxt_energy: geom_create failed",0,GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     &  call errquit("wfn1_nxt_energy: geom load failed",0,RTDB_ERR)
      if (.not. geom_nuc_charge(geom,nuclear_charge))
     &  call errquit("wfn1_nxt_energy: geometry nuclear charge failed",
     &               0,GEOM_ERR)
      if (.not. rtdb_get(rtdb,"charge",MT_DBL,1,charge))
     &  charge = 0.0d0
c
      if (.not. geom_print(geom))
     &  call errquit("wfn1_nxt_energy: geometry print failed",
     &               0,GEOM_ERR)
c
      nelec = int(nuclear_charge - charge)
      if (nelec.lt.0)
     &   call errquit("wfn1_nxt_energy: negative number of electrons",
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
     &     call errquit("wfn1_nxt_energy: more unpaired electrons than"
     &                //" all electrons",nopen,INPUT_ERR)
        if (mod(nelec-nopen,2).ne.0)
     &     call errquit("wfn1_nxt_energy: inconsistent no. electrons "
     &                //"and spin",nopen,INPUT_ERR)
        nclosed = (nelec-nopen)/2
        nea = nclosed + nopen
        neb = nclosed
      else if (scftype.eq.'UHF') then
        if (.not. rtdb_get(rtdb,'scf:nalpha',MT_INT,1,nea))
     &     call errquit("wfn1_nxt_energy: no. alpha electrons not "
     &                  //"found",0,INPUT_ERR)
        if (.not. rtdb_get(rtdb,'scf:nbeta',MT_INT,1,neb))
     &     call errquit("wfn1_nxt_energy: no. beta electrons not found",
     &                  0,INPUT_ERR)
      else
        call errquit("wfn1_nxt_energy: unexpect scftype: "//scftype,
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
     &  call errquit("wfn1_nxt_energy: bas_create failed",0,BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &  call errquit("wfn1_nxt_energy: basis load failed",0,RTDB_ERR)
      if (.not. int_normalize(rtdb,basis))
     &  call errquit("wfn1_nxt_energy: normalization failed",0,INT_ERR)
      if (.not. bas_numbf(basis,nbf))
     &  call errquit("wfn1_nxt_energy: no basis functions",0,BASIS_ERR)
      if (.not. bas_name(basis, trans, trans))
     &  call errquit('wfn1_nxt_energy: bad basis',0, BASIS_ERR)
      if (.not. bas_print(basis))
     &  call errquit("wfn1_nxt_energy: no basis",0,BASIS_ERR)
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
     &   call errquit("wfn1_nxt_energy: could not allocate alpha mos",
     &                0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,'beta mos',-1,-1,g_movecs(2)))
     &   call errquit("wfn1_nxt_energy: could not allocate beta mos",
     &                0,GA_ERR)
      call ga_zero(g_movecs(1))
      call ga_zero(g_movecs(2))
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo evals',
     &   l_eval, k_eval)) call errquit
     &   ('wfn1_nxt_energy: insufficient memory?',
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
      k_ehfa = k_eval
      k_ehfb = k_eval + nbf
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo occ',
     &   l_occ, k_occ)) call errquit
     &   ('wfn1_nxt_energy: insufficient memory?',
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
     &   call errquit("wfn1_energy: MO-vectors from where?",0,INPUT_ERR)
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
     &         call errquit('wfn1_nxt_energy: movecs_read failed ',
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
            call errquit("wfn1_nxt_energy: no suitable vectors found",
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
         call errquit("wfn1_nxt_energy: no vectors found",0, DISK_ERR)
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
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"state",
     &                     l_stat,k_stat))
     &   call errquit("wfn1_nxt_energy: could not allocate state",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
      k_stat = k_stat - 1
      call wfn1_nxt_init_state(g_movecs,nbf,
     &     dbl_mb(k_ehfa),dbl_mb(k_ehfb),
     &     dbl_mb(k_stat+wfn1_nxt_aocc()),
     &     dbl_mb(k_stat+wfn1_nxt_anat()),
     &     dbl_mb(k_stat+wfn1_nxt_aoccl()),
     &     dbl_mb(k_stat+wfn1_nxt_anatl()),
     &     dbl_mb(k_stat+wfn1_nxt_bocc()),
     &     dbl_mb(k_stat+wfn1_nxt_bnat()),
     &     dbl_mb(k_stat+wfn1_nxt_boccl()),
     &     dbl_mb(k_stat+wfn1_nxt_bnatl()))
      k_stat = k_stat + 1
cDEBUG
      write(*,*)'*** Initial wavefunction'
      call wfn1_nxt_print_state(dbl_mb(k_stat),nbf)
cDEBUG
c
      call wfn1_nxt_min_energy(rtdb,geom,basis,nbf,nea,neb,
     &     dbl_mb(k_stat),etot)


c
      return
      end
C>
C> \brief Initialize the addressing common block
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. Because both data structures are serialized we need to set
C> the appropriate off-sets to find the various components. We assume that
C> the first element starts at index 1 (just like in regular Fortran).
C>
      subroutine wfn1_nxt_init(nbf)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      iaocc  = 1
      ianat  = iaocc  + nbf*nbf
      iaoccl = ianat  + nbf*nbf
      ianatl = iaoccl + nbf*nbf
      ibocc  = ianatl + nbf*nbf
      ibnat  = ibocc  + nbf*nbf
      iboccl = ibnat  + nbf*nbf
      ibnatl = iboccl + nbf*nbf
      iaconf = ibnatl + nbf*nbf
      ibconf = iaconf + nbf
      isize1 = 8*nbf*nbf+2*nbf
c
      end
C>
C> \brief Fetch alpha occupation function off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the alpha occupation
C> functions.
C>
      integer function wfn1_nxt_aocc()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_aocc = iaocc
c
      return
      end
C>
C> \brief Fetch alpha natural orbital off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the alpha natural
C> orbitals.
C>
      integer function wfn1_nxt_anat()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_anat = ianat
c
      return
      end
C>
C> \brief Fetch alpha occupation function lagrangian off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the alpha occupation
C> function lagrangians.
C>
      integer function wfn1_nxt_aoccl()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_aoccl = iaoccl
c
      return
      end
C>
C> \brief Fetch alpha natural orbital lagrangian off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the alpha natural
C> orbital lagrangians.
C>
      integer function wfn1_nxt_anatl()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_anatl = ianatl
c
      return
      end
C>
C> \brief Fetch beta occupation function off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the beta occupation
C> functions.
C>
      integer function wfn1_nxt_bocc()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_bocc = ibocc
c
      return
      end
C>
C> \brief Fetch beta natural orbital off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the beta natural
C> orbitals.
C>
      integer function wfn1_nxt_bnat()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_bnat = ibnat
c
      return
      end
C>
C> \brief Fetch beta occupation function lagrangian off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the beta occupation
C> function lagrangians.
C>
      integer function wfn1_nxt_boccl()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_boccl = iboccl
c
      return
      end
C>
C> \brief Fetch beta natural orbital lagrangian off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the beta natural
C> orbital lagrangians.
C>
      integer function wfn1_nxt_bnatl()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_bnatl = ibnatl
c
      return
      end
C>
C> \brief Fetch alpha electron configuration off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the alpha 
C> orbital configuration.
C>
      integer function wfn1_nxt_aconf()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_aconf = iaconf
c
      return
      end
C>
C> \brief Fetch beta electron configuration off-set
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the off-set of the beta 
C> orbital configuration.
C>
      integer function wfn1_nxt_bconf()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_bconf = ibconf
c
      return
      end
C>
C> \brief Fetch the size of 1D data structures
C>
C> In wfn1_next_step.dox we describe how the state function and the gradient
C> will be stored. This function returns the size of 1D data structures
C> such as the state function and the gradients
C>
      integer function wfn1_nxt_size1()
      implicit none
c
      integer iaocc  !< The alpha occupation function off-set
      integer ianat  !< The alpha natural orbital off-set
      integer iaoccl !< The alpha occupation function lagrangian off-set
      integer ianatl !< The alpha natural orbital lagrangian off-set
      integer ibocc  !< The beta occupation function off-set
      integer ibnat  !< The beta natural orbital off-set
      integer iboccl !< The beta occupation function lagrangian off-set
      integer ibnatl !< The beta natural orbital lagrangian off-set
      integer iaconf !< The alpha electron configuration
      integer ibconf !< The alpha electron configuration
      integer isize1 !< The total size of a 1D data structure
      common/wfn1_nxt_offset/iaocc,ianat,iaoccl,ianatl,
     +                       ibocc,ibnat,iboccl,ibnatl,
     +                       iaconf,ibconf,
     +                       isize1
c
      wfn1_nxt_size1 = isize1
c
      return
      end
C>
C> @}

