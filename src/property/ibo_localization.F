
      subroutine ibo_localization(rtdb, geom, ltyp, basis, g_movecs,
     &  nocc, nvir, nmo, nbf, natoms, eval, occ, c, pop, list)

c     =================================================================
c     IAO construction and generation of occupied or virtual IBOs.
c     Written by J. Autschbach, jochena@buffalo.edu
c
c     The formalism is described in G. Knizia, J. Chem. Theory
c     Comput. 2013, 9, 4834 (with additional useful information provided
c     in Senjean et al., J. Chem. Theory Comput. 2021, 17, 1337 and the
c     SI of that article). What is implemented here corresponds to
c     Knizia's original formulation, and the equations for p = 4 in
c     Appendix D as given in the 'bugfixed' version of the article on
c     the author's web page [https://sites.psu.edu/knizia/software/
c     (accessed 02/2023)].

c     g_movecs and arrays occ and eval are for the MOs of a given
c     spin (alpha or beta)

c     The localization routine called from here is in
c     src/ddscf/localize.F

c     =================================================================

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"


c     subroutine arguments:

      integer rtdb, geom, basis
      character*(3) ltyp
      integer g_movecs
      integer nocc, nvir, nmo, nbf, natoms
      double precision eval(nbf), occ(nbf), c(nbf,2)
      double precision pop(natoms)
      integer list(natoms)

c     local GA handles:

      integer g_s2, g_s12, g_p12, g_p21
      integer g_ctilde, g_iao, g_mo
      integer g_temp, g_tmp1, g_tmp2, g_cib, g_u, g_vt
      integer g_smat

c     local variables:

      logical master, debug
      integer dbg, info

      integer ncore, imo, jmo, kmo

      integer n1, n2, nsing, noff

      integer s, nlist, a, i, t, tt, u, bflo, bfhi

      double precision rtemp, qas

      double precision small
      parameter (small=1d-8)

      integer
     &  mnbf, iao_mxprim, iao_mxang, iao_mxcont,
     &  iao_mxnbf_cn, iao_mxnbf_ce, iao_nshells,
     &  minbas

      double precision minval, swap
      integer l_val, k_val

      logical int_normalize
      external int_normalize

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      character*(256) lmotrans
      logical file_write_ga
      external file_write_ga

      character*(16) pname

c     =================================================================

      pname = 'ibo_localization'

      dbg = 0
      master =  ga_nodeid().eq.0
      debug = (dbg>0) .and. master ! .true. during development

      if (ltyp.ne.'occ' .and. ltyp.ne.'vir')  call errquit
     &     (pname//': loc. type unknown', 0, BASIS_ERR)

      if(debug) then
        if (ltyp.eq.'occ') write (luout,*)
     &       'entering occupied IBO localization'
        if (ltyp.eq.'vir') write (luout,*)
     &    'entering virtual IBO localization'
      end if

      if (.not. geom_num_core(rtdb, geom, 'ddscf', ncore)) ncore = 0

      if (debug) write (luout,*) 'ncore = ',ncore



c     ------------------------------
c     construct IAOs in the AO basis
c     ------------------------------

c     Note:
c     Basis 1 is the AO basis used in the SCF calculation
c     Basis 2 is the minimal auxiliary basis

c     AO Overlap Matrix S1 -> g_smat:

      g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
      call ga_zero(g_smat)
      call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

c     auxiliary basis overlap S2 -> g_s2 basis needs to defined in
c     the input as "iao basis". we will now create the basis here and
c     then calculate the overlap and the mixed ao-iao basis overlap
c     S12

      if (.not. bas_create(minbas, 'iao basis'))
     &  call errquit(pname//': cannot create iao bas', 86, BASIS_ERR)

      if (.not. bas_rtdb_load(rtdb, geom, minbas, 'iao basis'))
     &  call errquit(pname//': iao basis not present', 86, BASIS_ERR)

      if (.not. int_normalize(rtdb,minbas))
     &  call errquit(pname//': int_normalize iao', 86, INT_ERR)

      if (.not. bas_numbf(minbas, mnbf))
     &  call errquit(pname//': bas err iao mnbf', 86, BASIS_ERR)

      if (.not. bas_nprim_cn_max(minbas,iao_mxprim))
     &  call errquit(pname//': bas err iao mxprim', 86, BASIS_ERR)

      if (.not. bas_high_angular(minbas,iao_mxang))
     &  call errquit(pname//': bas err iao mxang', 86, BASIS_ERR)

      if (.not. bas_ncontr_cn_max(minbas,iao_mxcont))
     &  call errquit(pname//': bas err iao mxcont', 86, BASIS_ERR)

      if (.not. bas_nbf_cn_max(minbas,iao_mxnbf_cn))
     &  call errquit(pname//': bas err iao mcnbf_cn', 86, BASIS_ERR)

      if (.not. bas_nbf_ce_max(minbas,iao_mxnbf_ce))
     &  call errquit(pname//': bas err iao mxnbf_ce', 86, BASIS_ERR)

      if (.not. bas_numcont(minbas,iao_nshells))
     &  call errquit(pname//': bas err iao nshells', 86, BASIS_ERR)

c     not sure if it makes much sense in case the minimal basis for
c     the IAOs is larger than our regular basis, so let's check that
c     here and exit with an error if nbf < mnbf

      if (nbf.lt.mnbf) call errquit
     &  (pname//': nbf < mnbf. cannot handle!', 66, UNKNOWN_ERR)

c     create overlap for minbas, and the mixed basis-minbas
c     overlap S12.
c     Then we calculate the projectors P12 and P21
c     via solving S1 P12 = S12 for P12 etc.

      g_s2  = ga_create_atom_blocked(geom, minbas, 'loc:s2')
      call ga_zero(g_s2)
      call int_1e_ga(minbas, minbas, g_s2, 'overlap', .false.)

      if (.not. ga_create(MT_DBL, nbf, mnbf, 'loc:s12',
     &  nbf, 0, g_s12)) call errquit(pname//': S12',61, GA_ERR)

      call ga_zero(g_s12)
      call int_1e_ga(basis, minbas, g_s12, 'overlap', .false.)

      if (.not.ga_duplicate(g_s12, g_p12, 'g_p12')) call
     &  errquit(pname//': ga_dup P12',62,GA_ERR)

      call ga_copy(g_s12, g_p12)

#ifdef SCALAPACK
      info = ga_solve(g_smat, g_p12)
      if (debug) write(luout,*)
     &  'info from ga_solve P12:', info
#else
      if (debug) write(luout,*) 'call 1 ma_solve'
      call ma_solve(g_smat, g_p12)
#endif

      if (.not. ga_create(MT_DBL, mnbf, nbf, 'loc:p21',
     &  mnbf, 0, g_p21)) call errquit(pname//': P21',63, GA_ERR)

      call ga_transpose(g_s12, g_p21)
      if (debug) write(luout,*) 'P12 assembled'

c     S12 no longer needed
      if (.not. ga_destroy(g_s12))
     &  call errquit(pname//': ga_destroy failed g_s12',61, GA_ERR)

c : scalapack throws a 'ga_to_SL: rows/cols error' so we're using
c       ma_solve for the time being
c#ifdef SCALAPACK
c        info = ga_solve(g_s2, g_p21)
c        if (debug) write(luout,*)
c     &    'info from ga_solve P21:', info
c#else
      if (debug) write(luout,*) 'call 2 ma_solve'
      call ma_solve(g_s2, g_p21)
c#endif

c     g_s2 no longer needed
      if (.not. ga_destroy(g_s2)) call errquit(
     &  pname//': error destroying g_s2',0, GA_ERR)

      if (dbg>1) then
        call ga_print(g_p12)
        call ga_print(g_p21)
      end if

c     store P12 * P21 in array g_temp; we'll need it soon

      if (.not.ga_duplicate(g_smat, g_temp, 'g_temp')) call
     &  errquit(pname//': ga_dup temp', 71, GA_ERR)

      call ga_dgemm('n', 'n', nbf, nbf, mnbf,
     &  1.0d0, g_p12, g_p21, 0.0d0, g_temp)

      if (debug) write (luout,*) 'p12 p21 -> temp'

c     P21 is no longer needed

      if (.not. ga_destroy(g_p21))
     &  call errquit(pname//': ga_destroy failed g_p21',61, GA_ERR)

c     construct IAOS from occ. MOs:

c     copy the relevant CMOs to g_mo

      if (.not. ga_create(MT_DBL, nbf, nocc, 'loc:g_mo',
     &  nbf, 0, g_mo)) call errquit(pname//': g_mo',0, GA_ERR)

      call ga_copy_patch('n',
     &  g_movecs, 1, nbf, 1, nocc,
     &  g_mo,     1, nbf, 1, nocc)

      if (debug) write (luout,*) 'movecs(occ) -> mo'

c     create C-tilde from Appendix C of Knizia's IBO paper.  g_temp
c     holds P12 * P21; we won't need it after the next matrix
c     multiplication

      if (.not. ga_create(MT_DBL, nbf, nocc, 'loc:ctilde',
     &  nbf, 0, g_ctilde)) call errquit('loc_driver: sc',0, GA_ERR)

      call ga_dgemm('n', 'n', nbf, nocc, nbf,
     &  1.0d0, g_temp, g_mo, 0.0d0, g_ctilde)

      if (.not. ga_destroy(g_temp))
     &  call errquit(pname//': ga_destroy failed g_temp',71, GA_ERR)

      if (debug) write (luout,*) 'pre sym_orth 1'
      call sym_orth(g_ctilde, g_smat, nbf, nocc)

      if (dbg>1) call ga_print(g_ctilde)

c     now assemble the matrix products to form matrix A from Appendix
c     C of the IBO paper. The result will be accumulated in array
c     g_iao which are the IAOs in basis 1. The shape of the array is
c     the same as P12

      if (.not.ga_duplicate(g_smat, g_tmp1, 'g_tmp1')) call
     &  errquit(pname//': ga_dup tmp1', 72, GA_ERR)

      if (.not.ga_duplicate(g_smat, g_tmp2, 'g_tmp2')) call
     &  errquit(pname//': ga_dup tmp2', 72, GA_ERR)

      if (.not.ga_duplicate(g_smat, g_temp, 'g_temp')) call
     &  errquit(pname//': ga_dup temp', 72, GA_ERR)

      if (.not.ga_duplicate(g_p12, g_iao, 'g_iao')) call
     &  errquit(pname//': ga_dup iao', 72, GA_ERR)


c     C * C(T)
      call ga_dgemm('n', 't', nbf, nbf, nocc,
     &  1.0d0, g_mo, g_mo, 0.0d0, g_temp)

c     C * C(T) * S1 -> tmp1
      call ga_dgemm('n', 'n', nbf, nbf, nbf,
     &  1.0d0, g_temp, g_smat, 0.0d0, g_tmp1)

c     tilde-C * tilde-C(T)
      call ga_dgemm('n', 't', nbf, nbf, nocc,
     &  1.0d0, g_ctilde, g_ctilde, 0.0d0, g_temp)

c     tilde-C * tilde-C(T) * S1 -> tmp2
      call ga_dgemm('n', 'n', nbf, nbf, nbf,
     &  1.0d0, g_temp, g_smat, 0.0d0, g_tmp2)

c     C * C(T) * S1  * tilde-C * tilde-C(T) * S1

      call ga_dgemm('n', 'n', nbf, nbf, nbf,
     &  1.0d0, g_tmp1, g_tmp2, 0.0d0, g_temp)

c     multiply prev. result by P12 and store in g_iao

      call ga_dgemm('n', 'n', nbf, mnbf, nbf,
     &  1.0d0, g_temp, g_p12, 0.0d0, g_iao)

c     form C * C(T) -1
      call ga_add_todiag(g_tmp1, -1.0d0)

c     form tilde-C * tilde-C(T) -1
      call ga_add_todiag(g_tmp2, -1.0d0)

c     multiply the two prev. results with each other, then multiply
c     the result with P12 and add to g_iao

      call ga_dgemm('n', 'n', nbf, nbf, nbf,
     &  1.0d0, g_tmp1, g_tmp2, 0.0d0, g_temp)

      call ga_dgemm('n', 'n', nbf, mnbf, nbf,
     &  1.0d0, g_temp, g_p12, 1.0d0, g_iao)

c     g_ctilde no longer needed

      if (.not. ga_destroy(g_ctilde))
     &  call errquit(pname//': ga_destroy failed ctilde',62, GA_ERR)

c     delete temp arrays before calling sym_orth, which allocates its
c     own temp memory

      if (.not. ga_destroy(g_tmp1))
     &  call errquit(pname//': ga_destroy failed g_tmp1',71, GA_ERR)

      if (.not. ga_destroy(g_tmp2))
     &  call errquit(pname//': ga_destroy failed g_tmp2',71, GA_ERR)

      if (.not. ga_destroy(g_temp))
     &  call errquit(pname//': ga_destroy failed g_temp',71, GA_ERR)
c     last step: orthonormalize g_iao

      if (debug) write (luout,*) 'pre sym_orth 2'
      call sym_orth(g_iao, g_smat, nbf, mnbf)

      call ga_sync

      if (dbg>1) call ga_print(g_iao)

c     check orthogonality of IAOs w.r.t. metric S1
c     g_p12 is used for scratch

      if (.not. ga_create(MT_DBL, mnbf, mnbf, 'g_temp',
     &  mnbf, 0, g_temp)) call errquit(pname//': temp',73, GA_ERR)

      call ga_dgemm('n', 'n', nbf, mnbf, nbf,
     &  1.0d0, g_smat, g_iao, 0.0d0, g_p12)

      call ga_dgemm('t', 'n', mnbf, mnbf, nbf,
     &  1.0d0, g_iao, g_p12, 0.0d0, g_temp)

      if (dbg>1) call ga_print(g_temp)

      call ga_add_todiag(g_temp, -1.0d0)

      call ga_maxelt(g_temp, rtemp)
      if (master) then
        write(luout,'(/1x,a,1x,f18.8)')
     &    'IBO loc: largest element in C(iao,T) S C(iao) -1: ',
     &    rtemp
        write(luout,'(1x,a/1x,a/)')
     &    'Significant deviations from zero may indicate',
     &    'elevated numerical noise in the IAO generation'
      end if

      if (.not. ga_destroy(g_temp))
     &  call errquit(pname//': ga_destroy failed g_temp',73, GA_ERR)

c     g_p12 no longer needed

      if (.not. ga_destroy(g_p12))
     &  call errquit(pname//': ga_destroy failed g_p12',61, GA_ERR)

c     ---------------------------------------------------------------
c     IAOs are now in array g_iao. Next, generate occupied or virtual
c     IBOs, depending on the input settings (ltyp)
c     ---------------------------------------------------------------

      if (ltyp.eq.'occ') then

c       -----------------------
c       generate occupied IBOs:
c       -----------------------

c       note: g_mo already holds the occupied MOs

        if (debug) then
          write(luout,*) 'movecs(occ) -> mo'
        end if

c       few more sanity check, just in case
        if (nocc.gt.mnbf) call errquit
     &    (pname//': nocc > mnbf', 66, UNKNOWN_ERR)
        if (nocc.gt.nbf) call errquit
     &    (pname//': nocc > nbf', 66, UNKNOWN_ERR)

c       transform the set of occupied MOs to be localized to the IAO
c       basis via C(iao,T) S C(MO)

        if (.not. ga_create(MT_DBL, mnbf, nocc , 'loc:cib',
     &    mnbf, 0, g_cib))
     &    call errquit(pname//': cib',0, GA_ERR)

        if (.not. ga_create(MT_DBL, nbf, nocc , 'loc:tmp1',
     &    0, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

        call ga_dgemm('n', 'n', nbf, nocc, nbf,
     &    1.0d0, g_smat, g_mo, 0.0d0, g_tmp1)

        call ga_dgemm('t', 'n', mnbf, nocc, nbf,
     &    1.0d0, g_iao, g_tmp1, 0.0d0, g_cib)

        if (.not. ga_destroy(g_tmp1))
     &    call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

        if (dbg>1) call ga_print(g_cib)

        if (debug) write(luout,*) 'cib done'

c       save a copy of the starting MOs for later

        if (.not.ga_duplicate(g_cib, g_tmp2, 'g_tmp2')) call
     &    errquit(pname//': ga_dup cib', 0, GA_ERR)

        call ga_copy(g_cib, g_tmp2)

c       perform localization of the MOs in IAO basis:

        call localizeIBO(minbas, c, g_cib, nocc, nbf, mnbf,
     &    natoms)


      else if (ltyp.eq.'vir') then

c       ----------------------
c       generate virtual IBOs:
c       ----------------------

c       store the virtual MOs in g_mo. The array is allocated, so we
c       destroy it frst

        if (.not. ga_destroy(g_mo)) call errquit(
     &    pname//': error destroying g_mo',1, GA_ERR)

        if (.not. ga_create(MT_DBL, nbf, nvir, 'loc:g_mo',
     &    nbf, 0, g_mo)) call errquit(pname//': g_mo',1, GA_ERR)

        call ga_copy_patch('n',
     &    g_movecs, 1, nbf, nocc+1, nmo,
     &    g_mo,     1, nbf, 1,      nvir)

c       sanity check: for virtual localization, assume no. of
c       virtuals >= (mnbf-nocc)

        if (nvir.lt.(mnbf-nocc)) call errquit
     &    (pname//': nvir < mnbf-nocc. Cannot handle', 66,
     &    UNKNOWN_ERR)

        if (debug) then
          write(luout,*) 'movecs(vir) -> mo'
          write(luout,*) 'nocc, nvir, nmo, nbf, mnbf',
     &      nocc, nvir, nmo, nbf, mnbf
        end if

c       transform the set of virtual MOs to the IAO basis
c       via C(iao,T) S C(MO), store in g_tmp2 (and keep the array)

        if (.not. ga_create(MT_DBL, mnbf, nvir , 'loc:tmp2',
     &    0, 0, g_tmp2)) call errquit(pname//': tmp2',0, GA_ERR)

        if (.not. ga_create(MT_DBL, nbf, nvir , 'loc:tmp1',
     &    0, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

        call ga_dgemm('n', 'n', nbf, nvir, nbf,
     &    1.0d0, g_smat, g_mo, 0.0d0, g_tmp1)

        call ga_dgemm('t', 'n', mnbf, nvir, nbf,
     &    1.0d0, g_iao, g_tmp1, 0.0d0, g_tmp2)

        if (.not. ga_destroy(g_tmp1))
     &    call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

        if (debug .and. master) write(luout,*) 'cib before svd'

c       now perform SVD on the matrix C(vir) that we just calculated and
c       stored in array g_tmp2.  The "eigenvector" matrix U gives the
c       true virtual MOs in the IAO basis that we want to localize

c       #(sing. values) = min(mnbf,nvir) = mnbf , assuming virtuals.
c       The non-zero singular values must number mnbf - nocc, which
c       means we will localize a number of virtuals that's smaller than
c       nvir and also smaller than min(mnbf,nvir).

        nsing=min(mnbf,nvir)

        if (.not. ma_push_get(MT_DBL, nsing, 'svd-vals', l_val, k_val))
     &    call errquit(pname//': svd-vals', nsing, MA_ERR)

        if (.not. ga_create(MT_DBL, mnbf, nsing,'g_u', 0, 0, g_u))
     &    call errquit(pname//': svd-U', mnbf*nsing, GA_ERR)

        if (.not. ga_create(MT_DBL, nsing, nvir, 'g_vt',0, 0, g_vt))
     &    call errquit(pname//': svd-vT', nsing*nvir, GA_ERR)

        if (debug) write(luout,*) 'svd: m,n,nsing=',mnbf,nvir,nsing

        call wrap_dgesvd(g_tmp2, g_u, g_vt, dbl_mb(k_val))

c       wrap_dgesvd leaves the content of the array in arg 1 intact, no
c       matter what the options in the call to dgesvd are inside the
c       routine, so g_tmp2 still contains the initial (redundant)
c       virtual MO set expressed in the IAO basis. We'll use it later.

        if (debug) then
          write(luout,*) pname//' singular values:'
          do imo = 1,mnbf
            write(luout,*) imo, dbl_mb(k_val+imo-1)
          end do
        end if

c       count the non-zero singular values: re-define nsing

        nsing = 0
        do imo = 1,mnbf
          if (dbl_mb(k_val+imo-1) .gt.small) nsing = nsing+1
        end do
        if (master) write(luout,'(1x,a,i10)')
     &    'non-zero singular values: ',nsing

c       we should have as many non-zero singlular values as we have
c       additional IAOs beyond the number of occupied MOs. Check that,
c       and print a warning if the values differ.

        if (nsing .ne. (mnbf-nocc)) then
          if (master)
     &      write(luout,'(/1x,a/1x,a/1x,a/1x,a,i6/1x,a,i6/1x,a)')
     &      'WARNING: number of non-zero singular values differs from',
     &      ' the number of IAOs minus the number of occupied MOs.',
     &      ' Will continue suspiciously!',
     &      ' nsing       : ', nsing,
     &      ' mnbf - nocc : ', mnbf - nocc,
     &      ' The numbers in the two previous lines should be the same.'
        end if

c       array g_cib will hold the left-singular vectors corresponding to
c       non-zero singlular values. Since dgesvd returns them in
c       descending order, and they are by definition positive, we simply
c       take the first nsing vectors

        if (.not. ga_create(MT_DBL, mnbf, nsing,'loc:cib', 0, 0, g_cib))
     &    call errquit(pname//': cib', mnbf*nsing, GA_ERR)

        call ga_copy_patch('n',
     &    g_u,   1, mnbf, 1, nsing,
     &    g_cib, 1, mnbf, 1, nsing )

c       dellocate memory used for SVD

        if (.not. ga_destroy(g_vt))
     &    call errquit(pname//': ga_destroy failed g_vt',0, GA_ERR)

        if (.not. ga_destroy(g_u))
     &    call errquit(pname//': ga_destroy failed g_u',0, GA_ERR)

        if (.not. ma_pop_stack(l_val)) call errquit(
     &    pname//': error MA pop svd-vals',0, MA_ERR)

        if (dbg>1) call ga_print(g_cib)

        if (debug) write(luout,*) 'cib done'

c       perform localization of the MOs in IAO basis:

        call localizeIBO(minbas, c, g_cib, nsing, nbf, mnbf,
     &    natoms)

      end if ! ltyp occ or vir

c     localization interations done

c     calculate the localization transform. The starting MOs were
c     saved in g_tmp2 in the IAO basis
c     CMOs(iao,T) * LMOs(iao) = localization transform -> g_tmp1

      if (ltyp.eq.'occ') then

        if (.not. ga_create(MT_DBL, nocc, nocc, 'loc:tmp1',
     &    0, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

        n1 = nocc ! dimensioning used later
        n2 = nocc
        noff = 0

        call ga_dgemm('t', 'n', nocc, nocc, mnbf,
     &    1.0d0, g_tmp2, g_cib, 0.0d0, g_tmp1)

      if (debug) write(luout,*) 'ibo occ: transformation in g_tmp1'

      elseif (ltyp.eq.'vir') then

        n1 = nvir ! dimensioning used later
        n2 = nsing
        noff = nocc ! offset for virtual orbitals

        if (.not. ga_create(MT_DBL, nvir, nsing, 'loc:tmp1',
     &    0, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

        call ga_dgemm('t', 'n', nvir, nsing, mnbf,
     &    1.0d0, g_tmp2, g_cib, 0.0d0, g_tmp1)

      end if ! ltyp

      if(noff+n1 > nbf) call errquit (pname//
     &  ': noff + n1 > nbf', noff, UNKNOWN_ERR)

      if (.not. ga_destroy(g_tmp2))
     &  call errquit(pname//': ga_destroy failed g_tmp2',0, GA_ERR)

c     transform MO energies. We assume that the Fock matrix was diagonal
c     in the basis of input MOs. We also assume that the loop dimensions
c     do not exceed the capacity of array c(:,:)

      do jmo = 1,n2
        call ga_get(g_tmp1, 1, n1, jmo, jmo, c(1,1), 1)
        rtemp = 0.0d0
        do imo = 1,n1
          rtemp = rtemp + eval(imo+noff) * c(imo,1)**2
        end do ! jmo
        c(jmo,2) = rtemp ! transformed MO energies
      end do ! imo
      do imo = 1,n2
        eval(imo+noff) = c(imo,2)
      end do

      if (debug) write(luout,*) 'transformed MO energies'

c     sort IBOs by energy. We'll have to do this on the master
c     node, or otherwise avoid having multiple nodes working on the
c     same data

      if (master) then
      do imo = 1,n2 -1
        minval = eval(imo+noff)
c       find lowest eval(kmo) below eval(imo) (with offsets)
        kmo = 0
        do jmo = imo+1, n2
          if (eval(jmo+noff).lt.minval) then
            kmo = jmo
            minval = eval(jmo+noff)
          end if
        end do
c       if kmo > 0 we swap imo and jmo
        if (kmo.ne.0) then
          if (debug) write(luout,*) 'swap ',imo,kmo
c         swap eigenvalues
          swap = eval(imo+noff)
          eval(imo+noff) = eval(kmo+noff)
          eval(kmo+noff) = swap
c         swap columns of LMO transformation
          call ga_get(g_tmp1, 1, n1, imo, imo, c(1,1), 1)
          call ga_get(g_tmp1, 1, n1, kmo, kmo, c(1,2), 1)
          call ga_put(g_tmp1, 1, n1, imo, imo, c(1,2), 1)
          call ga_put(g_tmp1, 1, n1, kmo, kmo, c(1,1), 1)
c         swap LMOs in IAO basis, for the analysis below
          call ga_get(g_cib, 1, mnbf, imo, imo, c(1,1), 1)
          call ga_get(g_cib, 1, mnbf, kmo, kmo, c(1,2), 1)
          call ga_put(g_cib, 1, mnbf, imo, imo, c(1,2), 1)
          call ga_put(g_cib, 1, mnbf, kmo, kmo, c(1,1), 1)
        end if
      end do
      end if ! master

      call ga_sync

c     write transformation to scratch file

      call util_file_name('lmotrans',.true.,.true.,lmotrans)
      if(.not.file_write_ga(lmotrans,g_tmp1)) call errquit
     &  (pname//': could not write lmotrans',0, DISK_ERR)

      if (debug) write (luout,*)
     &  pname//': lmos reordered and lmotrans written'

c     transform CMOs to re-ordered LMOs

      if (.not. ga_create(MT_DBL, nbf, n2 , 'loc:tmp2',
     &  nbf, 0, g_tmp2)) call errquit(pname//': tmp2',0, GA_ERR)

      call ga_dgemm('n', 'n', nbf, n2, n1,
     &  1.0d0, g_mo, g_tmp1, 0.0d0, g_tmp2)


      if (debug) write (luout,*) pname//': g_tmp2 now has LMOs'

      if (.not. ga_destroy(g_tmp1))
     &  call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

c     copy the IBOs into the relevant part of movecs

      if (ltyp.eq.'occ') then
        call ga_copy_patch('n',
     &    g_tmp2,   1, nbf, 1, nocc,
     &    g_movecs, 1, nbf, 1, nocc)

      else if (ltyp.eq.'vir') then
        call ga_copy_patch('n',
     &    g_tmp2,   1, nbf, 1, nsing,
     &    g_movecs, 1, nbf, nocc+1, nocc+nsing)
      end if

      if (debug) write (luout,*) pname//': IBOs copied to g_movecs'

      if (.not. ga_destroy(g_tmp2))
     &  call errquit(pname//': ga_destroy failed g_tmp2',0, GA_ERR)

c     Analyze localization of each MO: per LMO, a list of atomic
c     populations is printed in decreasing magnitude, with the
c     polulations in parentheses.  This code is equivalent to the on
c     found in the Pipek-Mezey localization routine

      if (master) then
        write(luout,'(/1x,a/)') 'IAO-IBO localized orbitals'
        do s = 1, n2
          call ga_get(g_cib,  1, mnbf, s, s, c(1,1), 1)
          nlist = 0
          do a = 1, natoms
            if (.not. bas_ce2bfr(minbas, a, bflo, bfhi))
     &        call errquit('localized: basis ', 3,
     &        BASIS_ERR)
            qas  = 0.0d0
            do u = bflo, bfhi
              qas  = qas  + c(u,1)*c(u,1)
            end do
            if (abs(qas) .gt. 0.01d0) then
              nlist = nlist + 1
              list(nlist) = a
              pop(nlist) = qas
            end if
          end do
          do u = 1, nlist
            do t = 1, u-1
              if (abs(pop(t)).lt.abs(pop(u))) then
                rtemp = pop(u)
                pop(u) = pop(t)
                pop(t) = rtemp
                tt = list(u)
                list(u) = list(t)
                list(t) = tt
              end if
            end do
          end do
          write(luout,9002) s, eval(s+noff),
     &      occ(s+noff),(list(a), pop(a), a=1,nlist)
 9002     format(i5, 1x, f14.6,1x, f5.3, 1x,100(2x,i4,'(',f5.2,')'))
        end do
        write(luout,*)
        call util_flush(luout)
      end if

      if (.not. ga_destroy(g_cib))
     &  call errquit(pname//': ga_destroy failed g_cib',0, GA_ERR)

c     deallocate remaining arrays that are no longer needed

      if (.not. ga_destroy(g_iao)) call errquit(
     &  pname//': error destroying g_iao',0, GA_ERR)

      if (.not. ga_destroy(g_mo)) call errquit(
     &  pname//': error destroying g_mo',0, GA_ERR)

c     smat not needed anymore
      if (.not. ga_destroy(g_smat)) call errquit(
     &  pname//': error destroying g_smat',0, GA_ERR)

c     destroy minimal basis (iao basis)
      if (.not.bas_destroy(minbas))
     &  call errquit(pname//'iao bas_destroy failed',0,BASIS_ERR)

      if (ltyp.eq.'occ') then
        n1 = 1
        n2 = nocc
      else if (ltyp.eq.'vir') then
        n1 = nocc + 1
        n2 = nocc + nsing
      end if

      if (master) write(luout,
     &  '(/1x,a,a,a/1x,a/1x,i10,2x,a,1x,i10/)')
     &  'IBO localization (',ltyp,'): IBOs will be stored',
     &  'in file locorb.movecs, number ',n1, 'to', n2

      return

      end



c     =================================================================

      subroutine sym_orth(g_c, g_s, n1, n2)
      
      implicit none
      
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "numerical_constants.fh"
      
c     Given input vectors collected in matrix C, and a metric S,
c     calculate
c     C' = C times [transpose(C) S C]**(-1/2) such that
c     transpose(C') S C' = unit matrix
      
c     matrix S is square symmetric, dimension n1, n1
c     matrix C has dimension n1, n2
      
c     in the calling routine, matrix C may be allocated with more than
c     n2 columns, which is the reason we use ga_copy_patch below
      
      integer g_c !  in/out
      integer g_s !  input
      
      integer n1, n2 ! input
      
      integer g_sc, g_csc
      
      character*(8) pname
      
      logical debug
      
c     -----------------------------------------------------------------
      
      pname = 'sym_orth'
      
      debug = .false. .and. ga_nodeid().eq.0 !during development
      
      if (debug) write(luout,*) 'entering '//pname
      
      if(.not.ga_create(mt_dbl,n1,n2,'g_sc',0,0,g_sc))
     &  call errquit(pname//': error creating g_sc',141, GA_ERR)
      
      if(.not.ga_create(mt_dbl,n2,n2,'g_csc',0,0,g_csc))
     &  call errquit(pname//': error creating g_csc',142, GA_ERR)
      
      call ga_dgemm('n','n',n1, n2, n1, 1.0d0,
     &  g_s, g_c, 0.0d0, g_sc)
      
      call ga_dgemm('t','n',n2, n2, n1, 1.0d0,
     &  g_c, g_sc, 0.0d0, g_csc)
      
      call ga_matpow(g_csc, -FP_Half, 1d-12)
      
      call ga_copy_patch('n',
     &  g_c,  1, n1, 1, n2,
     &  g_sc, 1, n1, 1, n2)
      
      call ga_dgemm('n','n',n1, n2, n2, 1.0d0,
     &  g_sc, g_csc, 0.0d0, g_c)
      
      
      if (.not. ga_destroy(g_sc))
     &  call errquit(pname//': error destroy g_sc',141, GA_ERR)
      
      if (.not. ga_destroy(g_csc))
     &  call errquit(pname//': error destroy g_csc',142, GA_ERR)
      
      if (debug) write(luout,*) 'leaving '//pname
      
      return
      
      end
      
c     =================================================================

      subroutine wrap_dgesvd(g_a, g_u, g_vt, values)
      
      implicit none
      
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      
      integer g_a, g_u, g_vt
      double precision values(*)

c     wrapper routine for call to dgesvd. This is similar to ga_svd_seq
c     in util/ga_it2.F but some of the array dimensions had to be
c     adjusted to match the usage case in the virtual IBO construction
c     where m < n
      
      integer n, m, type, l_a, k_a, l_u, k_u, l_vt, k_vt,
     &  l_work, k_work, lwork, info, nsing

      character*(11) pname

c     -----------------------------------------------------------------

      pname = 'wrap_dgesvd'
      
      call ga_inquire(g_a, type, m, n)

      if (type.ne.MT_DBL) then
        call errquit (pname//': array type not MT_DBL',
     &    type, UNKNOWN_ERR)
      end if
      
      nsing = min(m,n)

      if (ga_nodeid() .eq. 0) then
        
        lwork = 10*max(m,n)
        if (.not. ma_push_get(type, m*n, 'svd', l_a, k_a))
     &    call errquit(pname//': a',m*n, MA_ERR)
        if (.not. ma_push_get(type, m*nsing, 'svd', l_u, k_u))
     &    call errquit(pname//': u',m*nsing, MA_ERR)
        if (.not. ma_push_get(type, nsing*n, 'svd', l_vt, k_vt))
     &    call errquit(pname//': vt',nsing*n, MA_ERR)
        if (.not. ma_push_get(type, lwork, 'svd', l_work, k_work))
     &    call errquit(pname//': work',lwork, MA_ERR)
        
        call ga_get(g_a, 1, m, 1, n, dbl_mb(k_a), m)

        call dgesvd('s','s',m,n,dbl_mb(k_a),m,values,
     &    dbl_mb(k_u),m,dbl_mb(k_vt),nsing,
     &    dbl_mb(k_work),lwork,info)
           
        if (info.ne.0) call errquit(pname//':svd failed',info,MEM_ERR)
        
        call ga_put(g_u,  1, m,     1, nsing, dbl_mb(k_u),  m)
        call ga_put(g_vt, 1, nsing, 1, n,     dbl_mb(k_vt), nsing)

        if (.not. ma_chop_stack(l_a)) call errquit(pname//': svd ma',0,
     &    MA_ERR)
         
      end if
      
      call ga_sync()
      call ga_brdcst(1,values,nsing*8,0)
      call ga_sync()

      end      
