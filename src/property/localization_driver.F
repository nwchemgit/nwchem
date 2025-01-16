      subroutine localization_driver (rtdb, basis, geom, loctype)

c     $Id$

c     ==================================================================
c
c     purpose: perform orbital localization, HF and DFT


c     input:   rtdb    - nwchem run time database 
c              basis   - basis set handle
c              geom    - geometry handle
c              loctype - localization type
c                        1: Pipek-Mezey (PM)
c                        2: Boys
c                        3: IAO-IBO
c     
c
c     driver witten by by J. Autschbach, jochena@buffalo.edu
c     localization routines are in src/ddscf/localize.F

c     ==================================================================

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
#include "util.fh"
#include "nwc_const.fh"

c     subroutine arguments:
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer loctype

c     local GA handles:
      integer g_uc(4), g_smat, g_sc, g_t
      integer g_movecs(2), g_cmo(2), g_temp, g_tmp1, g_tmp2
      integer g_sc0, g_iao(2)

c     MA variables:
      integer l_c, k_c, l_sc, k_sc, l_eval, k_eval, l_occ, k_occ
      integer l_dip(3), k_dip(3)
      integer l_pop, k_pop, l_cpop, k_cpop, l_list, k_list
      integer l_iloc, k_iloc

c     other local variables:

      integer loc_opt

      integer nclosed(2), nopen(2), nvirt(2), ncore, nocc, nvir,
     &  nbf, nmo, nloc, i, natoms, nspin

      character*(2) lspin(2)
      data lspin /'_A','_B'/
      character*(10) filename

      integer info, dbg, ispin
      
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)

      double precision atn, tol2e, val, rtmp, sumd2, dx, dy, dz, rtemp
      character*3 scftyp

      logical     oskel, status
      data tol2e   /1.0d-10/

      double precision dummy(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      double precision angstrom
      parameter (angstrom = 0.529177d0)

      character*(256) lmotrans
      logical file_write_ga
      external file_write_ga

      logical debug, master
      logical oprint, have_iao(2), have_mbs
      data have_iao(1), have_iao(2)/.false.,.false./
      integer minbas
      data have_mbs/.false./

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      logical int_normalize
      external int_normalize

      character*(19) pname
      character*(3) ltyp

      character*(16) tag 
      double precision atxyz(3), q_nuc, q_mol

c     ==================================================================

      pname = 'localization_driver'

      dbg = 0
      master =  ga_nodeid().eq.0
      debug = (dbg>0) .and. master ! .true. during development

      if (debug) write(luout,*) 'entering '//pname

      if (loctype.lt.1 .or. loctype.gt.3)
     &   call errquit(pname//': loctype out of range',0, RTDB_ERR)

      oprint = util_print('information', print_low)
      oprint = oprint .and. ga_nodeid().eq.0
      
c     ----------------------
c     start-up: print header
c     ----------------------
      
      if (oprint) then
        write (luout,9000) ! general header
        if (loctype.eq.1) then
          call util_print_centered(luout,
     &       'Pipek-Mezey localization', 25, .true.)
        elseif (loctype.eq.2) then
          call util_print_centered(luout,
     &      'Boys localization', 25, .true.)
        elseif (loctype.eq.3) then
          call util_print_centered(luout,
     &      'IAO-IBO localization', 25, .true.)
        else
          write (luout,*) 'loctype =',loctype
          write (luout,*) 'LocType value not supported. Exiting.'
          return
        end if
      end if

c     retrieve localization option from rtdb if it exists

      if (.not. rtdb_get(rtdb, 'prop:loc_opt', mt_int, 1,
     &                               loc_opt  ))
     &  loc_opt = 0

      if (oprint) then
        if (loc_opt .eq.0) then
          call util_print_centered(luout,
     &      '(occupied orbitals)', 25, .false.)
        else if (loc_opt.eq.1) then
          call util_print_centered(luout,
     &      '(virtual orbitals)', 25, .false.)
        else if (loc_opt.eq.2) then
          call util_print_centered(luout,
     &      '(occ and virt orbitals)', 25, .false.)
        else
          call errquit(pname//': localization option out of range',
     &      loc_opt, INPUT_ERR)
        end if
      end if

c     -------------------------------------------------------
c     assemble some data:
c     MO coeffs, orbital energies, density matrix, occupation
c     numbers, basis set info, matrix elements, etc.
c     -------------------------------------------------------

      if (.not. bas_numbf(basis,nbf)) call
     &   errquit(pname//': could not get nbf',0, BASIS_ERR)
c     allocate dbl_mb(k_occ) = fractional occupation numbers:
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit(pname//': ma_push_get failed k_occ',0,MA_ERR)
c     allocate dbl_mb(leval) = orbital energies:
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit(pname//': ma_push_get failed k_eval',0,MA_ERR)

      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,g_movecs,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)

      if (debug)
     &   write (luout,*) 'loc driver: nclosed,nopen,nvirt',nclosed(1),
     &  nopen(1), nvirt(1),nclosed(2),
     &  nopen(2), nvirt(2)

      natoms = 0

      if (.not. geom_ncent(geom, natoms)) call errquit
     &  (pname//': geom',0, GEOM_ERR)

      if (natoms.lt.1) call errquit
     &  (pname//': could not get natoms from rtdb', 911, UNKNOWN_ERR)

      if (natoms.gt.nw_max_atom) call errquit
     &  (pname//': nw_max_atom too small ', 911, UNKNOWN_ERR)

c     Check if we have a spin-restricted or unrestricted
c     calculation and set nspin accoringly

      nspin = 0

      if (scftyp.eq.'UHF') then

        if (master) write (luout,'(/1x,a)')
     &    'UNRESTRICTED calculation'
        nspin = 2
        if (nclosed(1).ne.0 .or. nclosed(2).ne.0) call
     &    errquit(pname//': UHF but nclosed.ne.0',1, CALC_ERR)

      else if (scftyp.eq.'RHF') then
        nspin = 1
c     for RHF system there should be no open shell orbitals
      if (nopen(1).ne.0) call
     &    errquit(pname//': RHF but nopen.ne.0',nopen(1), CALC_ERR)

      else
c       at present, ROHF, MCSCF not supported
        call errquit(pname//': SCF type not supported: '//
     &    trim(scftyp),1,INPUT_ERR)
      end if

c     Restricted: nclosed and nvirt should add up to nmo
      if (nspin .eq. 1) then
      if ((nclosed(1)+nvirt(1)).ne.nmo) call
     &    errquit(pname//': wrong no. of rhf orbitals',
     &    nmo, CALC_ERR)
c     Unrestrcited: nopen + nvirt = nmo, per spin
      else if (nspin.eq.2) then
        if ((nopen(1)+nvirt(1)+nopen(2)+nvirt(2)).ne.(2*nmo)) call
     &    errquit(pname//': wrong no. of uhf orbitals',
     &    nmo, CALC_ERR)
      end if

      if (nspin.le.0 .or. nspin.gt.2) call
     &  errquit(pname//': nspin wrong',nspin, CALC_ERR)

      if (debug) write (luout,*) 'loc_driver: setup complete'


c     --------------------------------------------
c     setup complete. Now proceed to the localization
c     --------------------------------------------

c
c     Integral initialization
c     
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)

c     =================================================================

      if (loctype.eq.1) then

c       ========================
c       Pipek-Mezey localization
c       ========================

        if (.not. ma_push_get(mt_int, nmo, 'iloc', l_iloc,
     &    k_iloc)) call errquit(pname//': loc:iloc', 0, MA_ERR)

c       AO Overlap Matrix S:

        g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
        call ga_zero(g_smat)
        call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

c       allocate memory

        if (.not. ga_create(MT_DBL, nbf, nmo, 'loc:sc',
     $    nbf, 0, g_sc)) call errquit(pname//': sc',0, GA_ERR)

        if (.not.ga_duplicate(g_sc, g_sc0,'g_sc0')) call
     &     errquit(pname//': ga_create failed g_sc0',1,GA_ERR)

        if (.not. ma_push_get(mt_dbl, 2*nbf, 'sc', l_sc, k_sc))
     $       call errquit(pname//': ma for sc', 0, MA_ERR)
        if (.not. ma_push_get(mt_dbl, 2*nbf, 'c', l_c, k_c))
     $       call errquit(pname//': ma for c', 0, MA_ERR)
        if (.not. ma_push_get(mt_dbl, natoms, 'pop', l_pop, k_pop))
     &    call errquit(pname//': loc:pop', 0, MA_ERR)
         if (.not. ma_push_get(mt_int, natoms,'list', l_list, k_list))
     &     call errquit(pname//': loc:list', 0, MA_ERR)

         do ispin = 1,nspin

           if (nspin.eq.1) then
             nocc = nclosed(ispin)
           else if (nspin.eq.2) then
             nocc = nopen(ispin)
           else
             call errquit(pname//': pmloc nspin out of order',
     &         nspin, UNKNOWN_ERR)
           end if
           nvir = nvirt(ispin)

           if (nspin.gt.1 .and. master) then
             write(luout,"(/1x,6('=')/1x,'Spin ',i1/1x,6('='))") ispin
           end if

c          PM localization needs S*C for this spin:
           call ga_dgemm('n', 'n', nbf, nmo, nbf,
     $       1.0d0, g_smat, g_movecs(ispin), 0.0d0, g_sc)

c          save g_sc in g_sc0 to calculate the loc. transform later

           call ga_copy(g_sc, g_sc0)

           if (loc_opt.eq.0 .or. loc_opt.eq.2) then

             ltyp = 'occ'
             call pm_localization(rtdb, geom, ltyp, basis,
     &         g_movecs(ispin), g_sc, g_sc0,
     &         nocc, nvir, nmo, nbf, natoms,
     &         dbl_mb(k_eval+(ispin-1)*nbf),
     &         dbl_mb(k_occ+(ispin-1)*nbf),
     &         dbl_mb(k_c), dbl_mb(k_sc), int_mb(k_iloc),
     &         dbl_mb(k_pop), int_mb(k_list))

           end if

           if (loc_opt.eq.1 .or. loc_opt.eq.2) then

             ltyp = 'vir'
             call pm_localization(rtdb, geom, ltyp, basis,
     &         g_movecs(ispin), g_sc, g_sc0,
     &         nocc, nvir, nmo, nbf, natoms,
     &         dbl_mb(k_eval+(ispin-1)*nbf),
     &         dbl_mb(k_occ+(ispin-1)*nbf),
     &         dbl_mb(k_c), dbl_mb(k_sc), int_mb(k_iloc),
     &         dbl_mb(k_pop), int_mb(k_list))

           end if

c          write transformation for this spin to scratch file
c          (occ-occ transformation only)

           if (loc_opt.eq.0 .or. loc_opt.eq.2) then

             if (.not. ga_create(MT_DBL, nocc, nocc, 'loc:trans',
     &         nocc, 0, g_t)) call errquit(pname//' pmloc: g_t',
     &         1, GA_ERR)

             call ga_dgemm('t', 'n', nocc, nocc, nbf,
     &         1.0d0, g_sc0, g_movecs(ispin), 0.0d0, g_t)

             if (nspin.eq.1) then
               filename='lmotrans'
             else if (nspin.eq.2) then
               filename='lmotrans'//lspin(ispin)
             else
               call errquit (
     &           pname//' pmloc: nspin value not valid',1,
     &           UNKNOWN_ERR)
             end if

             call util_file_name(filename,.false.,.false.,lmotrans)
             if(.not.file_write_ga(lmotrans,g_t)) call errquit
     &         (pname//' pmloc: could not write lmotrans',1, DISK_ERR)

             if (master) write (luout,'(1x,a/1x,a)')
     &         'PM transformation (occ.) written to file ',
     &         trim(lmotrans)

             if (.not. ga_destroy(g_t)) call errquit(
     &         pname//' pmloc: error destroying g_t',1, GA_ERR)

           end if ! save lmotrans for occ MOs

         end do ! ispin

c        deallocate arrays

c        clean up some temp arrays
         if (.not. ga_destroy(g_sc)) call errquit(
     &      pname//': error destroying g_sc',0, GA_ERR)

         if (.not. ma_chop_stack(l_sc)) call errquit(
     &      pname//': error MA pop sc',0, MA_ERR)

         if (.not. ga_destroy(g_smat)) call errquit(
     &     pname//': error destroying g_smat',1, GA_ERR)
         if (.not. ga_destroy(g_sc0)) call errquit(
     &     pname//': error destroying g_sc0',1, GA_ERR)

c     =================================================================

      else if (loctype.eq.2) then

c       =================
c       Boys localization
c       =================

        if (nspin.gt.1) call
     &    errquit(pname//': Boys loc. not supported for open shell',
     &    nspin, INPUT_ERR)

        ispin = 1

        if (loc_opt.gt.0) call
     &    errquit(pname//': Boys loc. not supported for virtuals',
     &    loc_opt, INPUT_ERR)

        do i = 1,4
          if (.not. ga_create(MT_DBL, nbf, nbf, 'uc',
     $       nbf, 0, g_uc(i))) call errquit(pname//': uc'
     &       ,i,GA_ERR)
          call ga_zero(g_uc(i))
         end do

         if (debug) write(luout,*) 'g_uc created'

c        dipole moment AO matrices -> uc(1) - uc(3)
         call int_dip_ga(basis, basis, g_uc(1), g_uc(2), g_uc(3))

         if (debug) write(luout,*) 'dipole done'

c        calculate dipole-AO times C (mo-coeffs), store in uc(i)
         if (.not. ga_create(MT_DBL, nbf, nmo, 'sc',
     $        nbf, 0, g_sc)) call errquit(pname//': sc',0, GA_ERR)
         do i = 1, 3
            call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $           1.0d0, g_uc(i), g_movecs(ispin), 0.0d0, g_sc)
            call ga_copy_patch('n',g_sc,1,nbf,1,nmo,g_uc(i),1,nbf,1,nmo)
         end do

         if (debug) write(luout,*) 'g_uc complete'

c       AO Overlap Matrix S:

        g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
        call ga_zero(g_smat)
        call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

        if (debug) write(luout,*) 'smat done'

c       store S C in array uc(4)
        call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $     1.0d0, g_smat, g_movecs(ispin), 0.0d0, g_uc(4))

        if (debug) write(luout,*) 'g_uc(4) done'
c     
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'sc', l_sc, k_sc))
     $        call errquit(pname//': ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'c', l_c, k_c))
     $     call errquit(pname//': ma for c', 0, MA_ERR)

         if (.not. ma_push_get(mt_int, nmo, 'iloc', l_iloc,
     &     k_iloc)) call errquit(pname//': loc:iloc', 1, MA_ERR)

         if (debug) write(luout,*) 'MA c, sc complete'

c        localize core and occupied orbitals:

         do i = 1, nclosed(1)
           int_mb(k_iloc - 1 + i) = i
c           iloc(i) = i
         end do
         nloc = nclosed(1)

c        calculate orbital centroids and determine B2 = the sum of 
c        square distances of the centroids from the coordinate origin.
c        The Boys localization maximizes B2. 
c        we can use g_sc as a temp array as it is not further used

         if (.not. ga_destroy(g_sc)) call errquit(
     &      pname//': error destroying g_sc',0, GA_ERR)
         if (.not. ga_create(MT_DBL, nmo, nmo, 'sc',
     $      nbf, 0, g_sc)) call errquit(pname//': sc',0, GA_ERR)

         do i = 1,3
           if (.not. ma_push_get(mt_dbl, nmo, 'sc', l_dip(i), k_dip(i)))
     $        call errquit(pname//': ma for dip', i, MA_ERR)
         end do

         do i = 1,3
           call ga_zero(g_sc)
           call ga_dgemm('t','n',nmo,nmo,nbf, 1.0d0, g_movecs(ispin),
     &       g_uc(i), 0.0d0, g_sc)
           call ga_get_diagonal(g_sc, dbl_mb(k_dip(i)))
         end do

         if (ga_nodeid().eq.0) then
           write (luout,'(/t6,a,t15,a,t27,a,t39,a,t49,a)')
     &        '#','<x>','<y>','<z>','<r>**2'
           write (luout,'(1x,53(''-''))')
         end if
         sumd2 = 0.0d0
         do i = 1,nmo
           dx = dbl_mb(k_dip(1) + i-1) 
           dy = dbl_mb(k_dip(2) + i-1) 
           dz = dbl_mb(k_dip(3) + i-1)
           rtmp = dx**2 + dy**2 + dz**2
           sumd2 = sumd2 + rtmp
           if (ga_nodeid().eq.0) then
             write (luout,'(1x,i5,4f12.4)') i,dx,dy,dz,rtmp
           end if
         end do
c         if (ga_nodeid().eq.0) then
c           write (luout,'(1x,a,2f12.4)') 'B2 before localization, av:',
c     &        sumd2, sumd2/float(nloc)
c         end if

         do i = 3,1,-1
           if (.not. ma_pop_stack(l_dip(i))) call errquit(
     &        pname//': error MA pop dip',i, MA_ERR)
         end do

c        jochen: comment:
c        the Boys localization routine was already available
c        in nwchem
         call localizeFB(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, int_mb(k_iloc), nbf, nmo, g_movecs(ispin), g_uc)

c        calculate orbital centroids again and print information

         do i = 1,3
           if (.not. ma_push_get(mt_dbl, nmo, 'sc', l_dip(i), k_dip(i)))
     $        call errquit(pname//': ma for dip', i, MA_ERR)
         end do

         do i = 1,3
           call ga_zero(g_sc)
           call ga_dgemm('t','n',nmo,nmo,nbf, 1.0d0, g_movecs(ispin),
     &       g_uc(i), 0.0d0, g_sc)
           call ga_get_diagonal(g_sc, dbl_mb(k_dip(i)))
         end do

         if (ga_nodeid().eq.0) then
           write (luout,'(/t6,a,t15,a,t27,a,t39,a,t49,a)')
     &        '#','<x>','<y>','<z>','<r>**2'
           write (luout,'(1x,53(''-''))')
         end if
         sumd2 = 0.0d0
         do i = 1,nmo
           dx = dbl_mb(k_dip(1) + i-1)
           dy = dbl_mb(k_dip(2) + i-1)
           dz = dbl_mb(k_dip(3) + i-1)
           rtmp = dx**2 + dy**2 + dz**2
           sumd2 = sumd2 + rtmp
           if (ga_nodeid().eq.0) then
             write (luout,'(1x,i5,4f12.4)') i,dx,dy,dz,rtmp
           end if
         end do
c         if (ga_nodeid().eq.0) then
c           write (luout,'(1x,a,2f12.4)') 'B2 after localization, av:',
c     &        sumd2, sumd2/float(nloc)
c         end if

         do i = 3,1,-1
           if (.not. ma_pop_stack(l_dip(i))) call errquit(
     &        pname//': error MA pop dip',i, MA_ERR)
         end do

c        clean up  temp arrays:

         if (.not. ga_destroy(g_sc)) call errquit(
     &      pname//': error destroying g_sc',0, GA_ERR)

         if (.not. ma_pop_stack(l_c)) call errquit(
     &      pname//': error MA pop c',0, MA_ERR)

         if (.not. ma_pop_stack(l_sc)) call errquit(
     &      pname//': error MA pop sc',0, MA_ERR)

         do i = 1, 4
            if (.not. ga_destroy(g_uc(i)))
     &        call errquit(pname//': error destroying uc',i, GA_ERR)
         end do

c        smat not needed anymore
         if (.not. ga_destroy(g_smat)) call errquit(
     &     pname//': error destroying g_smat',0, GA_ERR)


c      ===============================================================

       else if (loctype.eq.3) then

c        ===================================================
c        IAO-IBO localization (occupied or virtual orbitals)
c        ===================================================

         g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
         call ga_zero(g_smat)
         call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

c        allocate memory

         if (.not. ma_push_get(mt_dbl, 2*nbf, 'c mo', l_c, k_c))
     &     call errquit(pname//': ma for c', 0, MA_ERR)

         if (.not. ma_push_get(mt_dbl, natoms, 'pop', l_pop, k_pop))
     &     call errquit(pname//': loc:pop', 0, MA_ERR)

         if (.not. ma_push_get(mt_dbl, natoms*nspin,'cpop',
     &     l_cpop,k_cpop))
     &     call errquit(pname//': loc:cpop', 0, MA_ERR)

         if (.not. ma_push_get(mt_int, natoms, 'list', l_list,
     &     k_list)) call errquit(pname//': loc:list', 0, MA_ERR)

         if (.not. ga_create(MT_DBL, nbf, nmo, 'loc:sc0',
     &     nbf, 0, g_sc0)) call errquit(pname//' ibo: sc0',0, GA_ERR)

c        initialize cumulative IAO populations with zeros
         do i=1,natoms*nspin
           dbl_mb(k_cpop+i-1)=0.0d0
         enddo

c        begin loop over spins

         do ispin = 1,nspin

           if (nspin.eq.1) then
             nocc = nclosed(ispin)
           else if (nspin.eq.2) then
             nocc = nopen(ispin)
           else
             call errquit(pname//': nspin out of order',
     &         nspin, UNKNOWN_ERR)
           end if
           nvir = nvirt(ispin)

           if (nspin.gt.1 .and. master) then
             write(luout,"(/1x,6('=')/1x,'Spin ',i1/1x,6('='))") ispin
           end if

           if (debug) write(luout,*) pname//': ispin,nocc,nvir,nmo',
     &       ispin,nocc,nvir,nmo

           if (nocc+nvir.ne.nmo) call errquit(pname//
     &       ': nocc+nvir.ne.nmo. Non-aufbau config?',nocc+nvir,
     &       UNKNOWN_ERR)

c          calculate SC with the original MOs, so we can generate
c          the localization transformation later for this spin

           call ga_dgemm('n', 'n', nbf, nmo, nbf,
     $       1.0d0, g_smat, g_movecs(ispin), 0.0d0, g_sc0)

c          note: upon first run of the IBO loc., we need to create the
c          set of IAOs and the minimal basis (minbas). In a subsequent
c          run, this is not required again. the corresponsing variable
c          have_iao(ispin) is initialized as .false. in the declarations
c          block of this routine. The GA g_iao(ispin) is allocated in
c          routine ibo_localization and deallocated here. Likewise,
c          minbas is created in ibo_localization but needs to be
c          destroyed here.

           have_mbs = .false.
           
c          localize occupied MOs:

           if (loc_opt.eq.0 .or. loc_opt.eq.2) then
             ltyp = 'occ'

             call ibo_localization(rtdb, geom, ltyp, basis,
     &         g_smat,
     &         g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &         dbl_mb(k_eval+(ispin-1)*nbf),
     &         dbl_mb(k_occ+(ispin-1)*nbf),
     &         dbl_mb(k_c),
     &         dbl_mb(k_pop),
     &         dbl_mb(k_cpop+(ispin-1)*natoms),
     &         int_mb(k_list),
     &         have_iao(ispin), g_iao(ispin), have_mbs, minbas)

             if (debug) write (luout,*)
     &         pname//': back from ibo_localization (occ)'
             
           end if

c          localize virtual MOs:

           if(loc_opt.eq.1 .or. loc_opt.eq.2) then             
             ltyp = 'vir'

             call ibo_localization(rtdb, geom, ltyp, basis,
     &         g_smat,
     &         g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &         dbl_mb(k_eval+(ispin-1)*nbf),
     &         dbl_mb(k_occ+(ispin-1)*nbf),
     &         dbl_mb(k_c),
     &         dbl_mb(k_pop),
     &         dbl_mb(k_cpop+(ispin-1)*natoms),
     &         int_mb(k_list),
     &         have_iao(ispin), g_iao(ispin), have_mbs, minbas)

             if (debug) write (luout,*)
     &         pname//': back from ibo_localization (vir)'

           end if

c          destroy minimal basis (iao basis), it was created in
c          routine ibo_localization
           if (.not.bas_destroy(minbas))
     &       call errquit(pname//'iao bas_destroy failed',0,BASIS_ERR)

c          write transformation for this spin to scratch file
c          (occ-occ transformation only)

           if (loc_opt.eq.0 .or. loc_opt.eq.2) then

             if (.not. ga_create(MT_DBL, nocc, nocc, 'loc:trans',
     &         nocc, 0, g_t)) call errquit(pname//' ibo: g_t',1, GA_ERR)

             call ga_dgemm('t', 'n', nocc, nocc, nbf,
     &         1.0d0, g_sc0, g_movecs(ispin), 0.0d0, g_t)

             if (nspin.eq.1) then
               filename='lmotrans'
             else if (nspin.eq.2) then
               filename='lmotrans'//lspin(ispin)
             else
               call errquit (
     &           pname//' ibo: nspin value not valid',1,
     &           UNKNOWN_ERR)
             end if

             call util_file_name(filename,.false.,.false.,lmotrans)
             if(.not.file_write_ga(lmotrans,g_t)) call errquit
     &         (pname//' ibo: could not write lmotrans',1, DISK_ERR)

             if (master) write (luout,'(1x,a/1x,a)')
     &         'IBO transformation (occ.) written to file ',
     &         trim(lmotrans)

             if (.not. ga_destroy(g_t)) call errquit(
     &         pname//' ibo: error destroying g_t',1, GA_ERR)

           end if ! save lmotrans for occ MOs

        end do                  ! ispin
        
c       Print IAO population-based atomic charges:
        
        if (.not. rtdb_get(rtdb, 'charge', mt_dbl, 1, q_mol))
     &    q_mol = 0.0d0
        
        if(master .and. loc_opt.ne.1) then  
          write(luout,'(/10x,a)') "   IAO-based populations    "
          write(luout, '(10x,a)') "----------------------------"
          
          write(luout,'(/10x,a)') "Atom    #        Charge"
          write(luout, '(10x,a)') "---------   ----------------"
          
          do i=1,natoms
c           fetch nuclear charge for this atom
            if (.not.geom_cent_get(geom, i, tag, atxyz, q_nuc)) 
     &        call errquit (
     &        pname//': cannot get atom info iao chg',1,
     &        UNKNOWN_ERR)

c           IAO atomic charges:
            rtmp = -dbl_mb(k_cpop+i-1)
            if (nspin.eq.1) rtmp = rtmp * 2.0d0
            rtmp = rtmp + q_nuc
            if (nspin.eq.2) rtmp = rtmp - dbl_mb(k_cpop+i-1+natoms)
            write(luout, '(t11,a,t14,i6,t23,f10.3)') trim(tag),i, rtmp
            dbl_mb(k_pop+i-1) = rtmp  
          end do                ! i=1,natoms
          write(luout, '(10x,a)') "----------------------------"
          
c         Now the IAO populations are in k_pop MA array.
c         Check if total charge = population; print warning if not.
c         Use q_nuc to store the total charge temporarily
          
          q_nuc = 0.0d0
          do i=1,natoms 
            q_nuc = q_nuc + dbl_mb(k_pop+i-1)
          end do
          write(luout, '(t11, a3, t23, f10.3/)') 'sum', q_nuc
          
          if (dabs(q_mol - q_nuc).gt.1e-8) then 
            write(luout,'(/1x,a,1x,a/)')
     &        "WARNING: sum of IAO-population based atomic",
     &        "charges not equal to total mol. charge in RTDB"
            write(luout,*) "mol charge:", q_mol
            write(luout,*) "sum of atomic charges", q_nuc
          end if
        end if                  ! master

c       save IAOs to an 'movecs' file in case we want to
c       visualize them

        alo(:) = 0
        ahi(:) = 0
        do ispin = 1,nspin
          call ga_inquire (g_iao(ispin), info, alo(ispin), ahi(ispin))
          if (debug) write(luout,*)
     &      'g_iao dims:', alo(ispin), ahi(ispin)
        end do
        
        if (nspin.eq.2 .and. ahi(1).ne.ahi(2)) then
          if (master) then
            write(luout,'(1x,a/1x,a/1x,a)')
     &        'Warning: The IAO sets for the two spins do not have',
     &        'the same column dimension. Unable to save iaos.movecs.',
     &        'You probably want to check what went wrong !'
          end if
        else

          call hnd_vec_write(rtdb,geom,basis,nbf,nclosed,nopen,
     &      nvirt,scftyp,g_iao,dbl_mb(k_occ),
     &      dbl_mb(k_eval),ahi(1), 'iaos.movecs')

          write(tag,'(i0)') ahi(1)
          if (master) write(luout,'(1x,a,a,a)')
     &      'IAOs saved to file iaos.movecs (',trim(tag),' orbitals)'

        end if
        
c       clean up arrays. g_iao(ispin) was allocated in routine
c       ibo_localization

        do ispin = 1,nspin
          if (.not. ga_destroy(g_iao(ispin))) call errquit(
     &      pname//': error destroying g_iao',ispin, GA_ERR)
        end do ! ispin
        if (.not. ga_destroy(g_smat)) call errquit(
     &    pname//' ibo: error destroying g_smat',1, GA_ERR)
        if (.not. ga_destroy(g_sc0)) call errquit(
     &    pname//' ibo: error destroying g_sc0',1, GA_ERR)

       end if ! loctype

c      ================================================================

c     -----------------------------------
c     localization done. Store MO vectors
c     -----------------------------------

      if (debug) write (luout,*) 'calling hnd_vec_write'
      call hnd_vec_write(rtdb,geom,basis,nbf,nclosed,nopen,
     &   nvirt,scftyp,g_movecs,dbl_mb(k_occ),
     &   dbl_mb(k_eval),nmo, 'locorb.movecs')


c     --------
c     clean up
c     --------


      if (.not.ga_destroy(g_movecs(1))) call
     &  errquit(pname//': ga_destroy failed g_movecs',0,GA_ERR)

      if (nspin.eq.2) then
        if (.not.ga_destroy(g_movecs(2))) call
     &    errquit(pname//': ga_destroy failed g_movecs-2',0,GA_ERR)
      end if

      if (.not.ma_chop_stack(l_occ)) call
     &   errquit(pname//': ma_chop_stack failed k_occ',l_occ,MA_ERR)

      call schwarz_tidy()
      call int_terminate()

 7000 continue

c     ---------------------------------------
c     Localization done. return
c     ---------------------------------------
      
      if (oprint) then
        write (LuOut,'(/1x,a)') 'Exiting Localization driver routine'
      endif

c     ==================================================================

      return



 9000 format(//80(1h-)//,
     &  10x,28(1h*),/,10x,'*** ORBITAL LOCALIZATION ***',/,
     &  10x,28(1h*)/)

 9020 format(/,10x,'Spin-unrestricted calculations not yet implemented',
     &      ' for the Localization module. Exiting')
      end
