C> \ingroup wfn1
C> @{
C>
C> \brief Initialization and wrapper routine around the actual solver
C>
C> This routine gets the main parameters for the calculation, arranges
C> memory for the core data structures, loads a UHF wavefunction to
C> initialize the wavefunction, and finally calls the solver for the
C> real work.
C>
C> As we now want to solve for multiple states simultaneously the
C> variable nst (indicating the number of states) is particularly 
C> important.
C>
      logical function wfn1_energy(rtdb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
c
      integer rtdb !< [Input] The runtime database handle
c
c     Local variables
c
      integer nst          !< The number of states to solve for
      integer nea          !< The number of alpha electrons
      integer neb          !< The number of beta electrons
      integer nperma       !< The number of alpha electron permutations
      integer npermb       !< The number of beta electron permutations
c     character*20 scftype !< The type of energy expression
      character*20 file_scftype !< The type of energy expression on file
c     character*256 movecs_in !< The input vectors filename
c     character*256 movecs_out !< The output vectors filename
      character*256 title_vecs !< The output vectors filename
c     integer geom         !< The geometry handle
c     integer basis        !< The basis set handle
c     integer nelec        !< The total number of electrons
c     integer nopen        !< The number of open-shell orbitals
c     integer nclosed   !< The number of closed-shell orbitals
      integer nbf_vecs  !< The number of basis functions as stored in
                        !< the movecs file
      character*256 bas_vecs  !< The basis set name
      integer ii        !< Counter
      double precision nuclear_charge !< The sum of all nuclear charges
c     double precision charge         !< The total charge of the molecule
      integer mxst
      double precision tol_act !< The threshold for active natural orbitals
      double precision temperature !< The electron temperature
      double precision power       !< The matrix power
      logical status
      integer nsets
      integer nmo_vecs(2)
      logical uwfn1    ! .true. then unrestricted open shell WFN1
      logical occopt   ! .true. then optimize occupation functions to
                       !        actual occupation numbers
      logical occenopt ! .true. then optimize the energy under the
                       !        constraint that the occupation functions
                       !        generate the actual occupation numbers
      logical anaocc   ! .true. then optimize occupation functions to
                       !        actual occupation numbers analyically
                       !        instead of numerically
      logical oclosed  ! .true. request same permutation for both
                       !        spin channels
c
      integer l_en,  k_en
      integer l_on,  k_on
      integer l_noa, k_noa
      integer l_nob, k_nob
      integer l_ofa, k_ofa
      integer l_ofb, k_ofb
      integer l_noa2, k_noa2
      integer l_nob2, k_nob2
      integer l_ofa2, k_ofa2
      integer l_ofb2, k_ofb2
      integer l_l1a, k_l1a
      integer l_l1b, k_l1b
      integer l_l2a, k_l2a
      integer l_l2b, k_l2b
      integer l_l3a, k_l3a
      integer l_l3b, k_l3b
      integer l_sm1, k_sm1
      integer l_acta, k_acta
      integer l_actb, k_actb
c
      integer l_dnoa, k_dnoa
      integer l_dnob, k_dnob
      integer l_dofa, k_dofa
      integer l_dofb, k_dofb
      integer l_dnoa2, k_dnoa2
      integer l_dnob2, k_dnob2
      integer l_dofa2, k_dofa2
      integer l_dofb2, k_dofb2
      integer l_rnoa, k_rnoa
      integer l_rnob, k_rnob
      integer l_rofa, k_rofa
      integer l_rofb, k_rofb
      integer l_dl1a, k_dl1a
      integer l_dl1b, k_dl1b
      integer l_dl2a, k_dl2a
      integer l_dl2b, k_dl2b
      integer l_dl3a, k_dl3a
      integer l_dl3b, k_dl3b
      integer l_occa, k_occa
      integer l_occb, k_occb
c
      integer l_ohpa, k_ohpa
      integer l_ohpb, k_ohpb
      integer l_ohma, k_ohma
      integer l_ohmb, k_ohmb
      integer l_olpa, k_olpa
      integer l_olpb, k_olpb
      integer l_olma, k_olma
      integer l_olmb, k_olmb
c
      integer l_nhpa, k_nhpa
      integer l_nhpb, k_nhpb
      integer l_nhma, k_nhma
      integer l_nhmb, k_nhmb
      integer l_nlpa, k_nlpa
      integer l_nlpb, k_nlpb
      integer l_nlma, k_nlma
      integer l_nlmb, k_nlmb
c
      integer l_doha, k_doha
      integer l_dohb, k_dohb
      integer l_dola, k_dola
      integer l_dolb, k_dolb
      integer l_dnha, k_dnha
      integer l_dnhb, k_dnhb
      integer l_dnla, k_dnla
      integer l_dnlb, k_dnlb
c
      integer l_ovnoa, k_ovnoa
      integer l_ovnob, k_ovnob
      integer l_ovofa, k_ovofa
      integer l_ovofb, k_ovofb
c
      integer l_ehfa, k_ehfa
      integer l_ehfb, k_ehfb
      integer l_ewfa, k_ewfa
      integer l_ewfb, k_ewfb
      integer l_enta, k_enta
      integer l_entb, k_entb
c
      integer l_h1,  k_h1
      integer l_ov,  k_ov
      integer l_eri, k_eri
      integer l_erix, k_erix
c
      integer l_prma, k_prma !< Memory for permutations
      integer l_prmb, k_prmb !< Memory for permutations
      integer l_sgna, k_sgna !< Memory for sign of permutations
      integer l_sgnb, k_sgnb !< Memory for sign of permutations
      integer l_ovla, k_ovla !< Memory for alpha electron overlaps
      integer l_ovlb, k_ovlb !< Memory for beta electron overlaps
c
      integer ist !< Counter over states
      integer ioa,iob,iva,ivb
c
      logical  movecs_read_header
      external movecs_read_header
      logical  movecs_read
      external movecs_read
      logical  int_normalize
      external int_normalize
      integer  wfn1_factorial
      external wfn1_factorial
c
c     Figure the geometry out
c
      if (.not. geom_create(geom,'geometry'))
     &  call errquit("wfn1_energy: geom_create failed",0,GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     &  call errquit("wfn1_energy: geom load failed",0,RTDB_ERR)
      if (.not. geom_nuc_charge(geom,nuclear_charge))
     &  call errquit("wfn1_energy: geometry nuclear charge failed",
     &               0,GEOM_ERR)
      if (.not. rtdb_get(rtdb,'charge',MT_DBL,1,charge))
     &  charge = 0.0d0
      if (.not. geom_print(geom))
     &  call errquit("wfn1_energy: geometry print failed",0,GEOM_ERR)
c
      nelec = int(nuclear_charge - charge)
      if (nelec.lt.0)
     &   call errquit("wfn1_energy: negative number of electrons",
     &                nelec,INPUT_ERR)
      if (.not. rtdb_get(rtdb,'wfn1:nstate',MT_INT,1,nst))
     &  nst = 1
      if (.not. rtdb_cget(rtdb,'scf:scftype',1,scftype)) scftype = 'RHF'
      if (scftype.eq.'RHF'.or.scftype.eq.'ROHF') then
        if (.not. rtdb_get(rtdb,'scf:nopen',MT_INT,1,nopen)) nopen = 0
        if (nopen.ne.0) then
          scftype = 'ROHF'
        else
          scftype = 'RHF'
        endif
        if (nopen.gt.nelec)
     &     call errquit("wfn1_energy: more unpaired electrons than"
     &                //" all electrons",nopen,INPUT_ERR)
        if (mod(nelec-nopen,2).ne.0) 
     &     call errquit("wfn1_energy: inconsistent no. electrons and"
     &                //" spin",nopen,INPUT_ERR)
        nclosed = (nelec-nopen)/2
        nea = nclosed + nopen
        neb = nclosed 
      else if (scftype.eq.'UHF') then
        if (.not. rtdb_get(rtdb,'scf:nalpha',MT_INT,1,nea))
     &     call errquit("wfn1_energy: no. alpha electrons not found",
     &                  0,INPUT_ERR)
        if (.not. rtdb_get(rtdb,'scf:nbeta',MT_INT,1,neb))
     &     call errquit("wfn1_energy: no. beta electrons not found",
     &                  0,INPUT_ERR)
      else
        call errquit("wfn1_energy: unexpect scftype: "//scftype,
     &               0,INPUT_ERR)
      endif
c
c     Figure the basis set out
c
      if (.not. bas_create(basis,'ao basis'))
     &  call errquit("wfn1_energy: bas_create failed",0,BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &  call errquit("wfn1_energy: basis load failed",0,RTDB_ERR)
      if (.not. int_normalize(rtdb,basis))
     &  call errquit("wfn1_energy: normalization failed",0,INT_ERR)
      if (.not. bas_numbf(basis,nbf))
     &  call errquit("wfn1_energy: no basis functions",0,BASIS_ERR)
      if (.not. bas_name(basis, trans, trans))
     &  call errquit('movecs_write: bad basis',0, BASIS_ERR)
      if (.not. bas_print(basis))
     &  call errquit("wfn1_energy: no basis",0,BASIS_ERR)
c
c     Now go and look for the wavefunction
c
      if (.not.ga_create(MT_DBL,nbf,nbf,'alpha mos',-1,-1,g_movecs(1)))
     &   call errquit("wfn1_energy: could not allocate alpha mos",
     &                0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,'beta mos',-1,-1,g_movecs(2)))
     &   call errquit("wfn1_energy: could not allocate beta mos",
     &                0,GA_ERR)
      call ga_zero(g_movecs(1))
      call ga_zero(g_movecs(2))
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo evals',
     &   l_eval, k_eval)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
      k_ehfa = k_eval
      k_ehfb = k_eval + nbf
      if (.not. ma_push_get(MT_DBL,nbf, 'scf_init: mo ewfa',
     &   l_ewfa, k_ewfa)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,nbf,MT_BYTE), MA_ERR)
      if (.not. ma_push_get(MT_DBL,nbf, 'scf_init: mo ewfb',
     &   l_ewfb, k_ewfb)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,nbf,MT_BYTE), MA_ERR)
      if (.not. ma_push_get(MT_DBL,nbf, 'scf_init: mo enta',
     &   l_enta, k_enta)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,nbf,MT_BYTE), MA_ERR)
      if (.not. ma_push_get(MT_DBL,nbf, 'scf_init: mo entb',
     &   l_entb, k_entb)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,nbf,MT_BYTE), MA_ERR)
c
      if (.not. ma_push_get(MT_DBL,2*nbf, 'scf_init: mo occ',
     &   l_occ, k_occ)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
c
      if (.not. ma_push_get(MT_INT,2*nbf, 'scf_init: mo irs',
     &   l_irs, k_irs)) call errquit
     &   ('scf_init: insufficient memory?', 
     &    ma_sizeof(MT_DBL,2*nbf,MT_BYTE), MA_ERR)
c
      call ifill(2*nbf, 1, int_mb(k_irs), 1) ! In case not adapting
c
c     Fill in the SCF convergence info
c
      if (.not. rtdb_get(rtdb, 'scf:maxiter', MT_INT, 1, maxiter))
     &   maxiter = 30
      if (.not.rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, gnorm_tol))
     &   gnorm_tol = 1.d-5
      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e))
     &   tol2e = min(1.0d-7,gnorm_tol*1d-5)
      if (.not. rtdb_get(rtdb,'wfn1:tol_act',MT_DBL,1,tol_act))
     &   tol_act = 1.0d-5
      if (.not. rtdb_get(rtdb,'wfn1:temperature',MT_DBL,1,temperature))
     &   temperature = 0.000d0
      if (.not. rtdb_get(rtdb,'wfn1:power',MT_DBL,1,power))
     &   power = 1.000d0
      if (.not. rtdb_get(rtdb,'wfn1:uwfn1',MT_LOG,1,uwfn1))
     &   uwfn1 = .true.
      if (.not. rtdb_get(rtdb,'wfn1:anaocc',MT_LOG,1,anaocc))
     &   anaocc = .false.
      if (.not. rtdb_get(rtdb,'wfn1:occopt',MT_LOG,1,occopt))
     &   occopt = .false.
      if (.not. rtdb_get(rtdb,'wfn1:occenopt',MT_LOG,1,occenopt))
     &   occenopt = .false.
      if (.not. rtdb_get(rtdb,'wfn1:oclosed',MT_LOG,1,oclosed))
     &   oclosed = .false.
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
c        3) nmo from file = current nmo (in case of change of linear dep)
c
         if ((nbf.eq.nbf_vecs) .and. (bas_vecs.eq.trans) .and.
     &       (nbf.eq.nmo_vecs(1))) then
            if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     &                            dbl_mb(k_ehfa), g_movecs))
     &         call errquit('scf_movecs_read failed ',0, DISK_ERR)
c
            if (file_scftype .eq. 'UHF' .or.
     &          file_scftype .eq. 'MCSCF') then
               if (.not. movecs_read(movecs_in, 2, dbl_mb(k_occ+nbf),
     &                               dbl_mb(k_ehfb),
     &                               g_movecs(2))) then
                  call ga_copy(g_movecs(1), g_movecs(2))
                  call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
                  call dcopy(nbf,dbl_mb(k_ehfa),1,dbl_mb(k_ehfb),1)
cDEBUG
                  write(*,*)'*** UHF: beta vecs not found'
               else
                  write(*,*)'*** UHF: read beta vecs'
                  do ii = 1, nbf
                    write(*,*)"occ: ",ii,dbl_mb(k_occ+ii-1),
     &                                   dbl_mb(k_occ+nbf+ii-1)
                  enddo
cDEBUG
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
            call errquit("wfn1_energy: no suitable vectors found",
     &                   0, DISK_ERR)
         endif
         if (.not.rtdb_get(rtdb,'scf:orb:ehfa',MT_DBL,
     $                     nbf,dbl_mb(k_ehfa))) then
           write(*,*)'WARNING: could load scf:orb:ehfa'
         endif
         if (.not.rtdb_get(rtdb,'scf:orb:ehfb',MT_DBL,
     $                     nbf,dbl_mb(k_ehfb))) then
           write(*,*)'WARNING: could load scf:orb:ehfb'
         endif
         if (.not.rtdb_get(rtdb,'mcscf:orb:occa',MT_DBL,
     $                     nbf,dbl_mb(k_occ))) then
           write(*,*)'WARNING: could load mcscf:orb:occa'
         endif
         if (.not.rtdb_get(rtdb,'mcscf:orb:occb',MT_DBL,
     $                     nbf,dbl_mb(k_occ+nbf))) then
           write(*,*)'WARNING: could load mcscf:orb:occb'
         endif
      else
         call errquit("wfn1_energy: no vectors found",0, DISK_ERR)
      endif
c
c     Work out permutation stuff
c
      nperma = wfn1_factorial(nea)
      npermb = wfn1_factorial(neb)
      if (.not.ma_push_get(MT_INT,nea*nperma,"perma",l_prma,k_prma))
     &   call errquit("wfn1_energy: could not allocate perma",
     &   ma_sizeof(MT_INT,nea*nperma,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_INT,neb*npermb,"permb",l_prmb,k_prmb))
     &   call errquit("wfn1_energy: could not allocate permb",
     &   ma_sizeof(MT_INT,neb*npermb,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nperma,"signa",l_sgna,k_sgna))
     &   call errquit("wfn1_energy: could not allocate signa",
     &   ma_sizeof(MT_DBL,nperma,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,npermb,"signb",l_sgnb,k_sgnb))
     &   call errquit("wfn1_energy: could not allocate signb",
     &   ma_sizeof(MT_DBL,npermb,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nea*nea,"overlapa",l_ovla,k_ovla))
     &   call errquit("wfn1_energy: could not allocate overlapa",
     &   ma_sizeof(MT_DBL,nea*nea,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,neb*neb,"overlapb",l_ovlb,k_ovlb))
     &   call errquit("wfn1_energy: could not allocate overlapb",
     &   ma_sizeof(MT_DBL,neb*neb,MT_BYTE),MA_ERR)
cDEBUG
      write(*,*)'*** HVD1 A: wfn1_energy'
      call util_flush(6)
      status = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_permutation_table(nea,nperma,int_mb(k_prma))
      call wfn1_permutation_table(neb,npermb,int_mb(k_prmb))
      call wfn1_calc_sign(nea,nperma,int_mb(k_prma),dbl_mb(k_sgna))
      call wfn1_calc_sign(neb,npermb,int_mb(k_prmb),dbl_mb(k_sgnb))
cDEBUG
      write(*,*)'*** HVD1 B: wfn1_energy'
      call util_flush(6)
      status = ma_verify_allocator_stuff()
cDEBUG
c
c     Arrange the required memory
c
      if (.not.ma_push_get(MT_DBL,nst,"En",l_en,k_en))
     &   call errquit("wfn1_energy: could not allocate En",
     &   ma_sizeof(MT_DBL,nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nst,"On",l_on,k_on))
     &   call errquit("wfn1_energy: could not allocate On",
     &   ma_sizeof(MT_DBL,nst,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"noa",l_noa,k_noa))
     &   call errquit("wfn1_energy: could not allocate noa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nob",l_nob,k_nob))
     &   call errquit("wfn1_energy: could not allocate nob",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"ofa",l_ofa,k_ofa))
     &   call errquit("wfn1_energy: could not allocate ofa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"ofb",l_ofb,k_ofb))
     &   call errquit("wfn1_energy: could not allocate ofb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"noa2",l_noa2,k_noa2))
     &   call errquit("wfn1_energy: could not allocate noa2",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nob2",l_nob2,k_nob2))
     &   call errquit("wfn1_energy: could not allocate nob2",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"ofa2",l_ofa2,k_ofa2))
     &   call errquit("wfn1_energy: could not allocate ofa2",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"ofb2",l_ofb2,k_ofb2))
     &   call errquit("wfn1_energy: could not allocate ofb2",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
c
c     if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnoa",l_dnoa,k_dnoa))
c    &   call errquit("wfn1_energy: could not allocate dnoa",
c    &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c     if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnob",l_dnob,k_dnob))
c    &   call errquit("wfn1_energy: could not allocate dnob",
c    &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c     if (.not.ma_push_get(MT_DBL,nbf*nbf,"dofa",l_dofa,k_dofa))
c    &   call errquit("wfn1_energy: could not allocate dofa",
c    &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c     if (.not.ma_push_get(MT_DBL,nbf*nbf,"dofb",l_dofb,k_dofb))
c    &   call errquit("wfn1_energy: could not allocate dofb",
c    &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnoa2",l_dnoa2,k_dnoa2))
     &   call errquit("wfn1_energy: could not allocate dnoa2",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnob2",l_dnob2,k_dnob2))
     &   call errquit("wfn1_energy: could not allocate dnob2",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dofa2",l_dofa2,k_dofa2))
     &   call errquit("wfn1_energy: could not allocate dofa2",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dofb2",l_dofb2,k_dofb2))
     &   call errquit("wfn1_energy: could not allocate dofb2",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l1a",l_l1a,k_l1a))
     &   call errquit("wfn1_energy: could not allocate l1a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l1b",l_l1b,k_l1b))
     &   call errquit("wfn1_energy: could not allocate l1b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l2a",l_l2a,k_l2a))
     &   call errquit("wfn1_energy: could not allocate l2a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l2b",l_l2b,k_l2b))
     &   call errquit("wfn1_energy: could not allocate l2b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l3a",l_l3a,k_l3a))
     &   call errquit("wfn1_energy: could not allocate l3a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"l3b",l_l3b,k_l3b))
     &   call errquit("wfn1_energy: could not allocate l3b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ohpa",l_ohpa,k_ohpa))
     &   call errquit("wfn1_energy: could not allocate ohpa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ohpb",l_ohpb,k_ohpb))
     &   call errquit("wfn1_energy: could not allocate ohpb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ohma",l_ohma,k_ohma))
     &   call errquit("wfn1_energy: could not allocate ohma",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ohmb",l_ohmb,k_ohmb))
     &   call errquit("wfn1_energy: could not allocate ohmb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"olpa",l_olpa,k_olpa))
     &   call errquit("wfn1_energy: could not allocate olpa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"olpb",l_olpb,k_olpb))
     &   call errquit("wfn1_energy: could not allocate olpb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"olma",l_olma,k_olma))
     &   call errquit("wfn1_energy: could not allocate olma",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"olmb",l_olmb,k_olmb))
     &   call errquit("wfn1_energy: could not allocate olmb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"nhpa",l_nhpa,k_nhpa))
     &   call errquit("wfn1_energy: could not allocate nhpa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"nhpb",l_nhpb,k_nhpb))
     &   call errquit("wfn1_energy: could not allocate nhpb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"nhma",l_nhma,k_nhma))
     &   call errquit("wfn1_energy: could not allocate nhma",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"nhmb",l_nhmb,k_nhmb))
     &   call errquit("wfn1_energy: could not allocate nhmb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nlpa",l_nlpa,k_nlpa))
     &   call errquit("wfn1_energy: could not allocate nlpa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nlpb",l_nlpb,k_nlpb))
     &   call errquit("wfn1_energy: could not allocate nlpb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nlma",l_nlma,k_nlma))
     &   call errquit("wfn1_energy: could not allocate nlma",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"nlmb",l_nlmb,k_nlmb))
     &   call errquit("wfn1_energy: could not allocate nlmb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"doha",l_doha,k_doha))
     &   call errquit("wfn1_energy: could not allocate doha",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dohb",l_dohb,k_dohb))
     &   call errquit("wfn1_energy: could not allocate dohb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnha",l_dnha,k_dnha))
     &   call errquit("wfn1_energy: could not allocate dnha",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dnhb",l_dnhb,k_dnhb))
     &   call errquit("wfn1_energy: could not allocate dnhb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dola",l_dola,k_dola))
     &   call errquit("wfn1_energy: could not allocate dola",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dolb",l_dolb,k_dolb))
     &   call errquit("wfn1_energy: could not allocate dolb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dnla",l_dnla,k_dnla))
     &   call errquit("wfn1_energy: could not allocate dnla",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dnlb",l_dnlb,k_dnlb))
     &   call errquit("wfn1_energy: could not allocate dnlb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"sm1",l_sm1,k_sm1))
     &   call errquit("wfn1_energy: could not allocate sm1",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"rnoa",l_rnoa,k_rnoa))
     &   call errquit("wfn1_energy: could not allocate rnoa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"rnob",l_rnob,k_rnob))
     &   call errquit("wfn1_energy: could not allocate rnob",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"rofa",l_rofa,k_rofa))
     &   call errquit("wfn1_energy: could not allocate rofa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"rofb",l_rofb,k_rofb))
     &   call errquit("wfn1_energy: could not allocate rofb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovnoa",l_ovnoa,k_ovnoa))
     &   call errquit("wfn1_energy: could not allocate ovnoa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovnob",l_ovnob,k_ovnob))
     &   call errquit("wfn1_energy: could not allocate ovnob",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovofa",l_ovofa,k_ovofa))
     &   call errquit("wfn1_energy: could not allocate ovofa",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovofb",l_ovofb,k_ovofb))
     &   call errquit("wfn1_energy: could not allocate ovofb",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dnoa",l_dnoa,k_dnoa))
     &   call errquit("wfn1_energy: could not allocate dnoa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dnob",l_dnob,k_dnob))
     &   call errquit("wfn1_energy: could not allocate dnob",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dofa",l_dofa,k_dofa))
     &   call errquit("wfn1_energy: could not allocate dofa",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nst,"dofb",l_dofb,k_dofb))
     &   call errquit("wfn1_energy: could not allocate dofb",
     &   ma_sizeof(MT_DBL,nbf*nbf*nst,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl1a",l_dl1a,k_dl1a))
     &   call errquit("wfn1_energy: could not allocate dl1a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl1b",l_dl1b,k_dl1b))
     &   call errquit("wfn1_energy: could not allocate dl1b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl2a",l_dl2a,k_dl2a))
     &   call errquit("wfn1_energy: could not allocate dl2a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl2b",l_dl2b,k_dl2b))
     &   call errquit("wfn1_energy: could not allocate dl2b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl3a",l_dl3a,k_dl3a))
     &   call errquit("wfn1_energy: could not allocate dl3a",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"dl3b",l_dl3b,k_dl3b))
     &   call errquit("wfn1_energy: could not allocate dl3b",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occa",l_occa,k_occa))
     &   call errquit("wfn1_energy: could not allocate occa",
     &   ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occb",l_occb,k_occb))
     &   call errquit("wfn1_energy: could not allocate occb",
     &   ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_LOG,nbf*1,"acta",l_acta,k_acta))
     &   call errquit("wfn1_energy: could not allocate acta",
     &   ma_sizeof(MT_LOG,nbf*1,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_LOG,nbf*1,"actb",l_actb,k_actb))
     &   call errquit("wfn1_energy: could not allocate actb",
     &   ma_sizeof(MT_LOG,nbf*1,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"h1",l_h1,k_h1))
     &   call errquit("wfn1_energy: could not allocate h1",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ov",l_ov,k_ov))
     &   call errquit("wfn1_energy: could not allocate ov",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"eri",l_eri,k_eri))
     &   call errquit("wfn1_energy: could not allocate eri",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"erix",l_erix,k_erix))
     &   call errquit("wfn1_energy: could not allocate erix",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
c
c     Initialize the arrays: All states start with the same natural
c     orbitals
c
      call ga_get(g_movecs(1),1,nbf,1,nbf,dbl_mb(k_noa),nbf)
      call ga_get(g_movecs(2),1,nbf,1,nbf,dbl_mb(k_nob),nbf)
      do ist = 1, nst-1
        call dcopy(nbf*nbf,dbl_mb(k_noa),1,dbl_mb(k_noa+ist*nbf*nbf),1)
        call dcopy(nbf*nbf,dbl_mb(k_nob),1,dbl_mb(k_nob+ist*nbf*nbf),1)
      enddo
cDEBUG
c     call dfill(nbf*nbf,0.0d0,dbl_mb(k_noa),1)
c     do ii = 0, nbf-1
c        dbl_mb(k_noa+ii*nbf+ii) = 1.0d0
c     enddo
c     call dcopy(nbf*nbf,dbl_mb(k_noa),1,dbl_mb(k_nob),1)
      do ist = 1, nst
        write(*,*)'natural orbitals: state ',ist
        ii = ist - 1
        call hess_hssout(dbl_mb(k_noa+ii*nbf*nbf),nbf,nbf,nbf)
        call hess_hssout(dbl_mb(k_nob+ii*nbf*nbf),nbf,nbf,nbf)
      enddo
c     status = ma_verify_allocator_stuff()
cDEBUG
      call dfill(nbf*nbf*nst,0.0d0,dbl_mb(k_ofa),1)
      do ii = 0, nbf-1
         dbl_mb(k_ofa+ii*nbf+ii) = 1.0d0
      enddo
      do ist = 1, nst-1
         call dcopy(nbf*nbf,dbl_mb(k_ofa),1,dbl_mb(k_ofa+nbf*nbf*ist),1)
      enddo
      call dcopy(nbf*nbf*nst,dbl_mb(k_ofa),1,dbl_mb(k_ofb),1)
      do ist = 1, nst
        call wfn1_swap(nea,neb,nbf,ist,ioa,iva,iob,ivb,oclosed)
        call wfn1_swap_vectors(nea,nbf,ioa,iva,
     +                         dbl_mb(k_ofa+nbf*nbf*(ist-1)))
        call wfn1_swap_vectors(neb,nbf,iob,ivb,
     +                         dbl_mb(k_ofb+nbf*nbf*(ist-1)))
      enddo
      if (.not.uwfn1) then
        call dcopy(nbf*nbf*nst,dbl_mb(k_ofa),1,dbl_mb(k_ofb),1)
        call dcopy(nbf*nbf*nst,dbl_mb(k_noa),1,dbl_mb(k_nob),1)
      endif
cDEBUG
!define DO_HELIUM_HF 1
#ifdef DO_HELIUM
      write(*,*)'*** HELIUM 6-31G'
      write(*,*)'*** correct energy = -2.870162138912'
      nea = 1
      neb = 1
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ofa),1)
      dbl_mb(k_ofa)   =  0.99783183
      dbl_mb(k_ofa+1) =  0.06574935
      dbl_mb(k_ofa+2) = -0.06574935
      dbl_mb(k_ofa+3) =  0.99783183
      call dcopy(nbf*nbf,dbl_mb(k_ofa),1,dbl_mb(k_ofb),1)
      dbl_mb(k_noa)   =  0.589831
      dbl_mb(k_noa+1) =  0.515906
      dbl_mb(k_noa+2) =  1.150974
      dbl_mb(k_noa+3) = -1.185952
      call dcopy(nbf*nbf,dbl_mb(k_noa),1,dbl_mb(k_nob),1)
#endif
#ifdef DO_HELIUM_HF
      write(*,*)'*** HELIUM HF/6-31G'
      write(*,*)'*** correct energy = -2.855160426166'
      nea = 1
      neb = 1
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ofa),1)
      dbl_mb(k_ofa)   =  1.00000000
      dbl_mb(k_ofa+1) =  0.00000000
      dbl_mb(k_ofa+2) =  0.00000000
      dbl_mb(k_ofa+3) =  1.00000000
      call dcopy(nbf*nbf,dbl_mb(k_ofa),1,dbl_mb(k_ofb),1)
      dbl_mb(k_noa)   =  0.592081
      dbl_mb(k_noa+1) =  0.513586
      dbl_mb(k_noa+2) =  1.149818
      dbl_mb(k_noa+3) = -1.186959
      call dcopy(nbf*nbf,dbl_mb(k_noa),1,dbl_mb(k_nob),1)
#endif
#ifdef DO_BERYLIUM
      write(*,*)'*** BERYLIUM 6-31G'
      write(*,*)'*** correct energy = -14.613545269641'
      nea = 2
      neb = 2
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ofa),1)
      dbl_mb(k_ofa)   =  0.99783183
      dbl_mb(k_ofa+1) =  0.06574935
      dbl_mb(k_ofa+2) = -0.06574935
      dbl_mb(k_ofa+3) =  0.99783183
      call dcopy(nbf*nbf,dbl_mb(k_ofa),1,dbl_mb(k_ofb),1)
      dbl_mb(k_noa)   =  0.589831
      dbl_mb(k_noa+1) =  0.515906
      dbl_mb(k_noa+2) =  1.150974
      dbl_mb(k_noa+3) = -1.185952
      call dcopy(nbf*nbf,dbl_mb(k_noa),1,dbl_mb(k_nob),1)
#endif
cDEBUG
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l1a),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l1b),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l2a),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l2b),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l3a),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_l3b),1)
      do ii = 0, nbf-1
c       dbl_mb(k_l1a+ii*nbf+ii) = dbl_mb(k_eval+ii)
c       dbl_mb(k_l1b+ii*nbf+ii) = dbl_mb(k_eval+nbf+ii)
c       dbl_mb(k_l3a+ii*nbf+ii) = dbl_mb(k_eval+ii)
c       dbl_mb(k_l3b+ii*nbf+ii) = dbl_mb(k_eval+nbf+ii)
      enddo
c
c     Go and solve the equations
c
      if (occopt) then
        if (anaocc) then
          write(*,*)
          write(*,*)'wfn1_match_occ: alpha ',ist+1
          call wfn1_match_occ(.true.,nbf,nea,dbl_mb(k_occ),
     &                        dbl_mb(k_ofa),
     &                        dbl_mb(k_ofa2),
     &                        dbl_mb(k_l1a),dbl_mb(k_l1b),
     &                        dbl_mb(k_l2a),dbl_mb(k_l2b),
     &                        dbl_mb(k_l3a),dbl_mb(k_l3b),
     &                        dbl_mb(k_dofa),dbl_mb(k_dofb),.true.)
          write(*,*)
          write(*,*)'wfn1_match_occ: beta  ',ist+1
          call wfn1_match_occ(.true.,nbf,neb,dbl_mb(k_occ+nbf),
     &                        dbl_mb(k_ofb),
     &                        dbl_mb(k_ofb2),
     &                        dbl_mb(k_l1a),dbl_mb(k_l1b),
     &                        dbl_mb(k_l2a),dbl_mb(k_l2b),
     &                        dbl_mb(k_l3a),dbl_mb(k_l3b),
     &                        dbl_mb(k_dofa),dbl_mb(k_dofb),.true.)
          write(*,*)
          do ist = 1, nst-1
            call dcopy(nbf*nbf,dbl_mb(k_ofa),1,
     &                         dbl_mb(k_ofa+ist*nbf*nbf),1)
            call dcopy(nbf*nbf,dbl_mb(k_ofb),1,
     &                         dbl_mb(k_ofb+ist*nbf*nbf),1)
            call wfn1_swap(nea,neb,nbf,ist+1,ioa,iva,iob,ivb,oclosed)
            call wfn1_swap_vectors(nea,nbf,ioa,iva,
     +                             dbl_mb(k_ofa+nbf*nbf*ist))
            call wfn1_swap_vectors(neb,nbf,iob,ivb,
     +                             dbl_mb(k_ofb+nbf*nbf*ist))
          enddo
        endif
        if (occenopt) then
          call wfn1_energyocc_num_doit(rtdb,geom,basis,1.0d-4,nbf,nst,
     &         nea,neb,
     &         dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &         dbl_mb(k_noa),dbl_mb(k_nob),dbl_mb(k_ofa),dbl_mb(k_ofb),
     &         dbl_mb(k_dnoa),dbl_mb(k_dnob),dbl_mb(k_dofa),
     &         dbl_mb(k_dofb),dbl_mb(k_rnoa),dbl_mb(k_rnob),
     &         dbl_mb(k_sm1),dbl_mb(k_en),
     &         dbl_mb(k_dl1a),dbl_mb(k_dl1b),
     &         dbl_mb(k_ehfa),dbl_mb(k_ehfb),
     &         dbl_mb(k_ewfa),dbl_mb(k_ewfb),
     &         dbl_mb(k_enta),dbl_mb(k_entb),power,
     &         dbl_mb(k_occ),dbl_mb(k_occ+nbf),
     &         tol_act)
        endif
      endif
      if (.true.) then
        call wfn1_energy_num_doit(rtdb,geom,basis,gnorm_tol,uwfn1,
     &     nbf,nst,nea,neb,
     &     dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &     dbl_mb(k_noa),dbl_mb(k_nob),dbl_mb(k_ofa),dbl_mb(k_ofb),
     &     dbl_mb(k_noa2),dbl_mb(k_nob2),dbl_mb(k_ofa2),dbl_mb(k_ofb2),
     &     dbl_mb(k_dnoa),dbl_mb(k_dnob),dbl_mb(k_dofa),dbl_mb(k_dofb),
     &     dbl_mb(k_dnoa2),dbl_mb(k_dnob2),
     &     dbl_mb(k_dofa2),dbl_mb(k_dofb2),
     &     dbl_mb(k_rnoa),dbl_mb(k_rnob),dbl_mb(k_rofa),dbl_mb(k_rofb),
     &     dbl_mb(k_sm1),dbl_mb(k_en),dbl_mb(k_on),
     &     dbl_mb(k_dl1a),
     &     dbl_mb(k_ohpa),dbl_mb(k_ohpb),dbl_mb(k_ohma),dbl_mb(k_ohmb),
     &     dbl_mb(k_olpa),dbl_mb(k_olpb),dbl_mb(k_olma),dbl_mb(k_olmb),
     &     dbl_mb(k_nhpa),dbl_mb(k_nhpb),dbl_mb(k_nhma),dbl_mb(k_nhmb),
     &     dbl_mb(k_nlpa),dbl_mb(k_nlpb),dbl_mb(k_nlma),dbl_mb(k_nlmb),
     &     dbl_mb(k_doha),dbl_mb(k_dohb),dbl_mb(k_dola),dbl_mb(k_dolb),
     &     dbl_mb(k_dnha),dbl_mb(k_dnhb),dbl_mb(k_dnla),dbl_mb(k_dnlb),
     &     dbl_mb(k_ovnoa),dbl_mb(k_ovnob),
     &     dbl_mb(k_ovofa),dbl_mb(k_ovofb),
     &     dbl_mb(k_l1a), dbl_mb(k_l1b),dbl_mb(k_l2a),
     &     log_mb(k_acta),log_mb(k_actb),
     &     dbl_mb(k_ehfa),dbl_mb(k_ehfb),
     &     dbl_mb(k_ewfa),dbl_mb(k_ewfb),
     &     dbl_mb(k_enta),dbl_mb(k_entb),temperature,power,
     &     nperma,npermb,int_mb(k_prma),int_mb(k_prmb),
     &     dbl_mb(k_sgna),dbl_mb(k_sgnb),
     &     dbl_mb(k_ovla),dbl_mb(k_ovlb),
     &     tol_act)
      else
        call wfn1_energy_doit(rtdb,geom,basis,gnorm_tol,nbf,nst,nea,neb,
     &       dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &       dbl_mb(k_noa),dbl_mb(k_nob),dbl_mb(k_ofa),
     &       dbl_mb(k_ofb),dbl_mb(k_l1a),dbl_mb(k_l1b),dbl_mb(k_l2a),
     &       dbl_mb(k_l2b),dbl_mb(k_l3a),dbl_mb(k_l3b),dbl_mb(k_en),
     &       dbl_mb(k_dnoa),dbl_mb(k_dnob),
     &       dbl_mb(k_dofa),dbl_mb(k_dofb),
     &       dbl_mb(k_dl1a),dbl_mb(k_dl1b),
     &       dbl_mb(k_dl2a),dbl_mb(k_dl2b),
     &       dbl_mb(k_dl3a),dbl_mb(k_dl3b),
     &       log_mb(k_acta),log_mb(k_actb),
     &       tol_act)
      endif
c
      write(*,*)'Total WFN1 energy = ',dbl_mb(k_en)
      if (.not. rtdb_put(rtdb, "wfn1:energy", MT_DBL, 1, dbl_mb(k_en)))
     &   call errquit("wfn1_energy: could not store energy",
     &                0, RTDB_ERR)
c
      if (.not.ma_pop_stack(l_erix))
     &   call errquit("wfn1_energy: could not deallocate erix",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_eri))
     &   call errquit("wfn1_energy: could not deallocate eri",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ov))
     &   call errquit("wfn1_energy: could not deallocate ov",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_h1))
     &   call errquit("wfn1_energy: could not deallocate h1",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_actb))
     &   call errquit("wfn1_energy: could not deallocate actb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_acta))
     &   call errquit("wfn1_energy: could not deallocate acta",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_occb))
     &   call errquit("wfn1_energy: could not deallocate occb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_occa))
     &   call errquit("wfn1_energy: could not deallocate occa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl3b))
     &   call errquit("wfn1_energy: could not deallocate dl3b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl3a))
     &   call errquit("wfn1_energy: could not deallocate dl3a",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl2b))
     &   call errquit("wfn1_energy: could not deallocate dl2b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl2a))
     &   call errquit("wfn1_energy: could not deallocate dl2a",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl1b))
     &   call errquit("wfn1_energy: could not deallocate dl1b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dl1a))
     &   call errquit("wfn1_energy: could not deallocate dl1a",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dofb))
     &   call errquit("wfn1_energy: could not deallocate dofb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dofa))
     &   call errquit("wfn1_energy: could not deallocate dofa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnob))
     &   call errquit("wfn1_energy: could not deallocate dnob",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnoa))
     &   call errquit("wfn1_energy: could not deallocate dnoa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_ovofb))
     &   call errquit("wfn1_energy: could not deallocate ovofb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ovofa))
     &   call errquit("wfn1_energy: could not deallocate ovofa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ovnob))
     &   call errquit("wfn1_energy: could not deallocate ovnob",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ovnoa))
     &   call errquit("wfn1_energy: could not deallocate ovnoa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_rofb))
     &   call errquit("wfn1_energy: could not deallocate rofb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_rofa))
     &   call errquit("wfn1_energy: could not deallocate rofa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_rnob))
     &   call errquit("wfn1_energy: could not deallocate rnob",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_rnoa))
     &   call errquit("wfn1_energy: could not deallocate rnoa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_sm1))
     &   call errquit("wfn1_energy: could not deallocate sm1",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_dnlb))
     &   call errquit("wfn1_energy: could not deallocate dnlb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnla))
     &   call errquit("wfn1_energy: could not deallocate dnla",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dolb))
     &   call errquit("wfn1_energy: could not deallocate dolb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dola))
     &   call errquit("wfn1_energy: could not deallocate dola",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnhb))
     &   call errquit("wfn1_energy: could not deallocate dnhb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnha))
     &   call errquit("wfn1_energy: could not deallocate dnha",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dohb))
     &   call errquit("wfn1_energy: could not deallocate dohb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_doha))
     &   call errquit("wfn1_energy: could not deallocate doha",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_nlmb))
     &   call errquit("wfn1_energy: could not deallocate nlmb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nlma))
     &   call errquit("wfn1_energy: could not deallocate nlma",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nlpb))
     &   call errquit("wfn1_energy: could not deallocate nlpb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nlpa))
     &   call errquit("wfn1_energy: could not deallocate nlpa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nhmb))
     &   call errquit("wfn1_energy: could not deallocate nhmb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nhma))
     &   call errquit("wfn1_energy: could not deallocate nhma",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nhpb))
     &   call errquit("wfn1_energy: could not deallocate nhpb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nhpa))
     &   call errquit("wfn1_energy: could not deallocate nhpa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_olmb))
     &   call errquit("wfn1_energy: could not deallocate olmb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_olma))
     &   call errquit("wfn1_energy: could not deallocate olma",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_olpb))
     &   call errquit("wfn1_energy: could not deallocate olpb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_olpa))
     &   call errquit("wfn1_energy: could not deallocate olpa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ohmb))
     &   call errquit("wfn1_energy: could not deallocate ohmb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ohma))
     &   call errquit("wfn1_energy: could not deallocate ohma",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ohpb))
     &   call errquit("wfn1_energy: could not deallocate ohpb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ohpa))
     &   call errquit("wfn1_energy: could not deallocate ohpa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_l3b))
     &   call errquit("wfn1_energy: could not deallocate l3b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_l3a))
     &   call errquit("wfn1_energy: could not deallocate l3a",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_l2b))
     &   call errquit("wfn1_energy: could not deallocate l2b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_l2a))
     &   call errquit("wfn1_energy: could not deallocate l2a",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_l1b))
     &   call errquit("wfn1_energy: could not deallocate l1b",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_l1a))
     &   call errquit("wfn1_energy: could not deallocate l1a",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_dofb2))
     &   call errquit("wfn1_energy: could not deallocate dofb2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dofa2))
     &   call errquit("wfn1_energy: could not deallocate dofa2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnob2))
     &   call errquit("wfn1_energy: could not deallocate dnob2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_dnoa2))
     &   call errquit("wfn1_energy: could not deallocate dnoa2",
     &   0,MA_ERR)
c     if (.not.ma_pop_stack(l_dofb))
c    &   call errquit("wfn1_energy: could not deallocate dofb",
c    &   0,MA_ERR)
c     if (.not.ma_pop_stack(l_dofa))
c    &   call errquit("wfn1_energy: could not deallocate dofa",
c    &   0,MA_ERR)
c     if (.not.ma_pop_stack(l_dnob))
c    &   call errquit("wfn1_energy: could not deallocate dnob",
c    &   0,MA_ERR)
c     if (.not.ma_pop_stack(l_dnoa))
c    &   call errquit("wfn1_energy: could not deallocate dnoa",
c    &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_ofb2))
     &   call errquit("wfn1_energy: could not deallocate ofb2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ofa2))
     &   call errquit("wfn1_energy: could not deallocate ofa2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nob2))
     &   call errquit("wfn1_energy: could not deallocate nob2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_noa2))
     &   call errquit("wfn1_energy: could not deallocate noa2",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ofb))
     &   call errquit("wfn1_energy: could not deallocate ofb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ofa))
     &   call errquit("wfn1_energy: could not deallocate ofa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_nob))
     &   call errquit("wfn1_energy: could not deallocate nob",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_noa))
     &   call errquit("wfn1_energy: could not deallocate noa",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_on))
     &   call errquit("wfn1_energy: could not deallocate On",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_en))
     &   call errquit("wfn1_energy: could not deallocate En",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_ovlb))
     &   call errquit("wfn1_energy: could not deallocate ovlb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ovla))
     &   call errquit("wfn1_energy: could not deallocate ovla",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_sgnb))
     &   call errquit("wfn1_energy: could not deallocate sgnb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_sgna))
     &   call errquit("wfn1_energy: could not deallocate sgna",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_prmb))
     &   call errquit("wfn1_energy: could not deallocate prmb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_prma))
     &   call errquit("wfn1_energy: could not deallocate prma",
     &   0,MA_ERR)

c
      if (.not.ma_pop_stack(l_irs))
     &   call errquit("wfn1_energy: could not deallocate irs",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_occ))
     &   call errquit("wfn1_energy: could not deallocate occ",
     &   0,MA_ERR)
c
      if (.not.ma_pop_stack(l_entb))
     &   call errquit("wfn1_energy: could not deallocate entb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_enta))
     &   call errquit("wfn1_energy: could not deallocate enta",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ewfb))
     &   call errquit("wfn1_energy: could not deallocate ewfb",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_ewfa))
     &   call errquit("wfn1_energy: could not deallocate ewfa",
     &   0,MA_ERR)
      if (.not.ma_pop_stack(l_eval))
     &   call errquit("wfn1_energy: could not deallocate eval",
     &   0,MA_ERR)
c
      if (.not.ga_destroy(g_movecs(1)))
     &   call errquit("wfn1_energy: could not destroy g_movecs(1)",
     &   0, GA_ERR)
      if (.not.ga_destroy(g_movecs(2)))
     &   call errquit("wfn1_energy: could not destroy g_movecs(2)",
     &   0, GA_ERR)
c
      wfn1_energy = .true.
      return
      end
C> @}
