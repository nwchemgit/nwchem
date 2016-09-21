      subroutine  waimfile(rtdb)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
c   This subroutine was designed to creat a WFN/WFX file
c   within a nwchem movecs file.
c   The format used is the same as the requeried for 
c   codes to do topological analyses as Atoms In Molecules 
c   theories and Electron Localization Function.
c
      integer       rtdb
      integer       geom, basis
      character*255 vec_file
      character*255 titlel_vec,  basis_name
      character*32  theory
c
      integer numcont           ! number of mapped contractions
      integer icont             ! contraction index
      integer type              ! type (sp/s/p/d/..)
c      integer nprim            ! no. of primitives per shell
      integer ngeno             ! no. of contractions
      integer sphcart           ! 0/1 for cartesian/shperical
      integer iset
      logical status
c
      integer kprim             ! sum of primitives per shell
      integer primprim          ! sum of the total of primitives used
      integer sprimprim         ! cartesian primitives after transformation from spherical
c
      double precision total_charge
c
      character*(nw_max_path_len) wfnfile
      character*255 title
      character*255 title_vec
      character*20 scftype
      integer nset
      integer nbf
      integer nmo(2)
      integer g_vecs(2)
      integer g_hcore(2)
      integer g_dens(2)
c
      double precision ekin, virial, epot
      double precision etot
c
      double precision forderco(21)
c
      integer nocc(2)
      integer noccs
c
      integer l_icent, k_icent, l_ricent, k_ricent
      integer l_itype, k_itype, l_ritype, k_ritype
      integer ls_ritype, ks_ritype
      integer ls_ricent, ks_ricent
      integer ls_rexp, ks_rexp
      integer ls_rcoef, ks_rcoef
      integer l_exp, k_exp, l_rexp, k_rexp
      integer l_rprim, k_rprim
      integer l_eval(2), k_eval(2)
      integer l_vecs, k_vecs  
      integer l_coef, k_coef, l_rcoef, k_rcoef !basis coefficients
      integer l_co, k_co
      integer ls_co, ks_co
      integer l_tprim, k_tprim
      integer l_nprim, k_nprim
      integer l_occ(2), k_occ(2)
c
      integer idum
      integer i,j,k,l,n,m
      integer naco
c
      integer natoms
c
      double precision coord(3,nw_max_atom)
      double precision xyz(3)
      character*16 tag
      character*80 buf
      double precision qnuc(nw_max_atom)
      character*4 symbol(nw_max_atom)
      integer atn
      character*16 rhf
c
      integer rprims(6)
      integer rprimc(6)
c functions per type or orbital
      data rprims /1, 3, 5, 7, 9, 11/ !wrong 
      data rprimc /1, 3, 6, 10, 15, 21/
c
      logical nowfx
c
      logical movecs_read_header
      external movecs_read_header
c
      logical movecs_read
      external movecs_read
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical int_normalize
      external int_normalize
c
      logical waimfile_printwfn, waimfile_printwfx, waimfile_sph2cart
      external waimfile_printwfn , waimfile_printwfx, waimfile_sph2cart
c
      integer cartbf(-1:5,1:21)
      integer sphebf(-1:5,1:11)
c

      character*4 atmsymbol(109)
      data atmsymbol/'  H ','  HE','  LI','  BE','  B ','  C ',
     *               '  N ','  O ','  F ','  NE','  NA','  MG',
     *               '  AL','  SI','  P ','  S ','  CL','  AR',
     *               '  K ','  CA','  SC','  TI','  V ','  CR',
     *               '  MN','  FE','  CO','  NI','  CU','  ZN',
     *               '  GA','  GE','  AS','  SE','  BR','  KR',
     *               '  RB','  SR','  Y ','  ZR','  NB','  MO',
     *               '  TC','  RU','  RH','  PD','  AG','  CD',
     *               '  IN','  SN','  SB','  TE','  I ','  XE',
     *               '  CS','  BA','  LA','  CE','  PR','  ND',
     *               '  PM','  SM','  EU','  GD','  TB','  DY',
     *               '  HO','  ER','  TM','  YB','  LU','  HF',
     *               '  TA','  W ','  RE','  OS','  IR','  PT',
     *               '  AU','  HG','  TL','  PB','  BI','  PO',
     *               '  AT','  RN','  FR','  RA','  AC','  TH',
     *               '  PA','  U ','  NP','  PU','  AM','  CM',
     *               '  BK','  CF','  ES','  FM','  MD','  NO',
     *               '  LR','  RF','  DB','  SG','  BH','  HS', 
     *               '  MT'/
c L
      cartbf(-1,1) = 1
      cartbf(-1,2) = 2
      cartbf(-1,3) = 3
      cartbf(-1,4) = 4
c 1 S
      cartbf(0, 1) = 1
c 2 PX
      cartbf(1, 1) = 2
c 3 PY
      cartbf(1, 2) = 3
c 4 PZ
      cartbf(1, 3) = 4
c 5 DXX
      cartbf(2, 1) = 5
c 8 DXY
      cartbf(2, 2) = 8
c 9 DXZ
      cartbf(2, 3) = 9
c 6 DYY
      cartbf(2, 4) = 6
c 10 DYZ
      cartbf(2, 5) = 10
c 7 DZZ
      cartbf(2, 6) = 7
c 11 FXXX
      cartbf(3, 1) = 11
c 14 FXXY
      cartbf(3, 2) = 14
c 15 FXXZ
      cartbf(3, 3) = 15
c 17 FXYY
      cartbf(3, 4) = 17
c 20 FXYZ
      cartbf(3, 5) = 20
c 18 FXZZ
      cartbf(3, 6) = 18
c 12 FYYY
      cartbf(3, 7) = 12
c 16 FYYZ
      cartbf(3, 8) = 16
c 19 FYZZ
      cartbf(3, 9) = 19
c 13 FZZZ
      cartbf(3,10) = 13
c 21 GXXXX
      cartbf(4, 1) = 21
c 24 GXXXY
      cartbf(4, 2) = 24
c 25 GXXXZ
      cartbf(4, 3) = 25
c 30 GXXYY
      cartbf(4, 4) = 30
c 33 GXXYZ
      cartbf(4, 5) = 33
c 31 GXXZZ
      cartbf(4, 6) = 31
c 26 GXYYY
      cartbf(4, 7) = 26
c 34 GXYYZ
      cartbf(4, 8) = 34
c 35 GXYZZ
      cartbf(4, 9) = 35
c 28 GXZZZ
      cartbf(4,10) = 28
c 22 GYYYY
      cartbf(4,11) = 22
c 27 GYYYZ
      cartbf(4,12) = 27
c 32 GYYZZ
      cartbf(4,13) = 32
c 29 GYZZZ
      cartbf(4,14) = 29
c 23 GZZZZ
c    cartbf(4,15) = 'zzzz'
      cartbf(4,15) = 23

c 56 HXXXXX (500)
c      cartbf(5, 1) = 'xxxxx'
      cartbf(5, 1) = 56
c 55 HXXXXY (410)
c      cartbf(5, 2) = 'xxxxy'
      cartbf(5, 2) = 55
c 54 HXXXXZ (401)
c      cartbf(5, 3) = 'xxxxz'
      cartbf(5, 3) = 54
c 53 HXXXYY (320)
c      cartbf(5, 4) = 'xxxyy'
      cartbf(5, 4) = 53
c 52 HXXXYZ (311)
c      cartbf(5, 5) = 'xxxyz'
      cartbf(5, 5) = 52
c 51 HXXXZZ (302)
c      cartbf(5, 6) = 'xxxzz'
      cartbf(5, 6) = 51
c 50 HXXYYY (230)
c      cartbf(5, 7) = 'xxyyy'
      cartbf(5, 7) = 50
c 49 HXXYYZ (221)
c      cartbf(5, 8) = 'xxyyz'
      cartbf(5, 8) = 49
c 48 HXXYZZ (212)
c      cartbf(5, 9) = 'xxyzz'
      cartbf(5, 9) = 48
c 47 HXXZZZ (203)
c      cartbf(5,10) = 'xxzzz'
      cartbf(5, 10) = 47
c 46 HXYYYY (140)
c      cartbf(5,11) = 'xyyyy'
      cartbf(5, 11) = 46
c 45 HXYYYZ (131)
c      cartbf(5,12) = 'xyyyz'
      cartbf(5, 12) = 45
c 44 HXYYZZ (122)
c      cartbf(5,13) = 'xyyzz'
      cartbf(5, 13) = 44
c 43 HXYZZZ (113)
c      cartbf(5,14) = 'xyzzz'
      cartbf(5, 14) = 43
c 42 HXZZZZ (104)
c      cartbf(5,15) = 'xzzzz'
      cartbf(5, 15) = 42
c 41 HYYYYY (050)
c      cartbf(5,16) = 'yyyyy'
      cartbf(5, 16) = 41
c 40 HYYYYZ (041)
c      cartbf(5,17) = 'yyyyz'
      cartbf(5, 17) = 40
c 39 HYYYZZ (032)
c      cartbf(5,18) = 'yyyzz'
      cartbf(5, 18) = 39
c 38 HYYZZZ (023)
c      cartbf(5,19) = 'yyzzz'
      cartbf(5, 19) = 38
c 37 HYZZZZ (014)
c      cartbf(5,20) = 'yzzzz'
      cartbf(5, 20) = 37
c 36 HZZZZZ (005)
c      cartbf(5,21) = 'zzzzz'
      cartbf(5, 21) = 36

      sphebf(-1,1) = 1
      sphebf(-1,2) = 2
      sphebf(-1,3) = 3
      sphebf(-1,4) = 4
      sphebf(0, 1) = 1
      sphebf(1, 1) = 2
      sphebf(1, 2) = 3
      sphebf(1, 3) = 4
      sphebf(2, 1) = 5
      sphebf(2, 2) = 6
      sphebf(2, 3) = 7
      sphebf(2, 4) = 8
      sphebf(2, 5) = 9
      sphebf(3, 1) = 10
      sphebf(3, 2) = 11
      sphebf(3, 3) = 12
      sphebf(3, 4) = 13
      sphebf(3, 5) = 14
      sphebf(3, 6) = 15
      sphebf(3, 7) = 16
      sphebf(4, 1) = 17
      sphebf(4, 2) = 18
      sphebf(4, 3) = 19
      sphebf(4, 4) = 20
      sphebf(4, 5) = 21
      sphebf(4, 6) = 22
      sphebf(4, 7) = 23
      sphebf(4, 8) = 24
      sphebf(4, 9) = 25
      sphebf(5, 1) = 26
      sphebf(5, 2) = 37
      sphebf(5, 3) = 38
      sphebf(5, 4) = 39
      sphebf(5, 5) = 40
      sphebf(5, 6) = 41
      sphebf(5, 7) = 42
      sphebf(5, 8) = 43
      sphebf(5, 9) = 44
      sphebf(5, 10) = 45
      sphebf(5, 11) = 46
     
c> open wfn file
c
      buf = ' '
      write(buf,*) '  WFX file creator '
      write(6,*)
      write(6,*)
      call util_print_centered(6,buf,40,.true.)
      write(6,*)
c
c> get title
      if(.not. rtdb_cget(rtdb, 'title', 1, title)) title = 'NWChem Job '

c> wfn / wfx
      if (.not. rtdb_get(rtdb, 'prop:nowfx',mt_log, 1, nowfx))
     $      nowfx=.true.
c theory
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('waimfile: theory not specified',555,
     &       INPUT_ERR)
c

c>  load geometry and symmetry info
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('waimfile : geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('waimfile : no geometry ', 0, RTDB_ERR)
      if (.not.geom_ncent(geom, natoms))
     $   call errquit('waimfile : error code = ',2,GEOM_ERR)
      if (.not. geom_nuc_charge(geom,total_charge))
     $    call errquit('waimfile: failed to read total charge',
     &     0, RTDB_ERR)

c> geometry
      do i=1,natoms
        if (.not.geom_cent_get(geom, i, tag, 
     $                           xyz, qnuc(i)))
     $       call errquit('waimfile : error code = ',3,GEOM_ERR)
        if (tag(1:2) .eq. 'bq' .or. 
     $       ((tag(1:1).eq.'x') .and. (tag(2:2).ne.'e'))) then ! X but not Xe
            symbol(i) = '  Bq '
        else
            symbol(i) = atmsymbol(int(qnuc(i)))
        endif
        coord(1,i)=xyz(1)
        coord(2,i)=xyz(2)
        coord(3,i)=xyz(3)
      enddo


c>  load the basis set and get info about it
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('waimfile : bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) then
         if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $        call errquit('waimfile : no mo or ao basis set', 0,
     &       RTDB_ERR)
      endif
c
      if (.not. int_normalize(rtdb,basis)) call errquit
     &    ('waimfile : int_normalize failed',911, INT_ERR)
c
c> READ vectors from movecs file
c
      if (.not. rtdb_cget(rtdb, 'prop:vectors',  1,
     $     vec_file)) then
         if (.not. rtdb_cget(rtdb, theory(1:inp_strlen(theory))//
     A          ':output vectors', 1, vec_file)) then
           call util_file_name('movecs', .false.,.false.,vec_file )
        endif
      endif
           write(6,*)' using vec ',vec_file,'@@'
      call util_file_name_resolve(vec_file, .false.)
c> movecs 
      if (.not. movecs_read_header(vec_file, title_vec, basis_name,
     &    scftype, nbf, nset, nmo, 2))
     &  call errquit('waimfile : basis set error:', 86, basis_err)
c
c
      if (.not. ma_push_get(mt_dbl,nbf,'occ',l_occ(1),k_occ(1)))
     &    call errquit('waimfile : failed to allocate orb.',nbf,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nbf,'eigenval',l_eval(1),k_eval(1)))
     &    call errquit('waimfile : failed to allocate orb.',nbf,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_read',
     $                                          l_vecs,k_vecs))
     $    call errquit('waimfile : ma failed', nbf, MA_ERR)

c     if (nset.eq.2) then
        if (.not. ma_push_get(mt_dbl,nbf,'occ',l_occ(2),k_occ(2)))
     &    call errquit('waimfile : failed to allocate orb.',nbf,MA_ERR)
        if (.not. ma_push_get(mt_dbl,nbf,'eigenval',
     &                                          l_eval(2),k_eval(2)))
     &    call errquit('waimfile : failed to allocate orb.',nbf,MA_ERR)
c     endif

      g_vecs(1) = ga_create_atom_blocked(geom, basis, 'vecs1')
      call ga_zero(g_vecs(1))
      if (nset.eq.2) then
        g_vecs(2) = ga_create_atom_blocked(geom, basis, 'vecs2')
        call ga_zero(g_vecs(2))
      endif
c
c> determine orb occupation
      nocc(1) = 0
      nocc(2) = 0
      do iset = 1, nset
        if (.not. movecs_read(vec_file, iset, dbl_mb(k_occ(iset)),
     $            dbl_mb(k_eval(iset)), g_vecs(iset))) call errquit
     $              ('waimfile  fragment: failed read fragment MOs',0,
     $               INPUT_ERR)
        do i=1,nbf
c> count electrons
        if ( dbl_mb(k_occ(iset)+i-1).gt.1.0d-1) 
     $                       nocc(iset) = nocc(iset) + 1
        enddo
      enddo

c>  total occupation
      noccs = nocc(1)
      if (nset.gt.1) noccs = noccs + nocc(2)
      

c> total of contractions
      if (.not.bas_numcont(basis,numcont))
     $   call errquit('waimfile : error code = ',2,BASIS_ERR)

c> expand basis representation

      if (.not. ma_push_get(mt_int,numcont,'nprim',l_nprim,k_nprim))
     $    call errquit('waimfile : ma failedi nprim', numcont, MA_ERR)
      if (.not. ma_push_get(mt_int,numcont,'tprim', l_tprim, k_tprim))
     $    call errquit('waimfile : ma failed tprim', numcont, MA_ERR)
      if (.not. ma_push_get(mt_int,numcont,'rprim',l_rprim,k_rprim))
     $    call errquit('waimfile : ma failed rprim', numcont, MA_ERR)
      if (.not. ma_push_get(mt_int,numcont,'icent',l_icent,k_icent))
     &  call errquit('waimfile : failed icent icent.',numcont,MA_ERR)
c
      primprim = 0
      kprim = 0
      do icont = 1, numcont
        if (.not. bas_continfo(basis, icont, type,
     $                    int_mb(k_nprim+icont-1), ngeno,sphcart))
     $    call errquit('waimfile : error code = ',3,BASIS_ERR)
c
        if (sphcart.eq.1) then
c          write (6,*) 'SPHERICAL BASIS WARNING'
           int_mb(k_rprim+icont-1) = rprims(type+1)
        else
           int_mb(k_rprim+icont-1) = rprimc(type+1)
        endif
c
        int_mb(k_tprim+icont-1) = 
     $                 int_mb(k_rprim+icont-1)*int_mb(k_nprim+icont-1)
        primprim = primprim + int_mb(k_tprim+icont-1)
c
        if (.not. bas_cn2ce(basis, icont, i))
     $    call errquit('waimfile : error code = ',5,BASIS_ERR)
c
        int_mb(k_icent+icont-1) = i
        kprim = kprim + int_mb(k_nprim+icont-1)
      enddo
c
c*      do icont = 1, numcont
c*        write(6,*) 'debug icont ', icont
c*        write(6,*) 'debug    nprim ', int_mb(k_nprim+icont-1)
c*        write(6,*) 'debug    rprim ', int_mb(k_rprim+icont-1)
c*        write(6,*) 'debug      tprim ', int_mb(k_tprim+icont-1)
c*      enddo
c
      if (.not. ma_push_get(mt_dbl,kprim,'exponents', l_exp, k_exp))
     $  call errquit('waimfile : ma failed exp', kprim, MA_ERR)
      if (.not. ma_push_get(mt_dbl,kprim,'coeffts', l_coef, k_coef))
     $  call errquit('waimfile : ma failed exp', kprim, MA_ERR)

      kprim = 0
      do icont = 1, numcont
        if (.not.bas_get_exponent(basis, icont,
     $                                     dbl_mb(k_exp + kprim)))
     $      call errquit('waimfile : failed bas_get_exponent',0,
     $       BASIS_ERR)

        if (.not.bas_get_coeff(basis, icont, dbl_mb(k_coef + kprim)))
     $         call errquit('waimfile : failed bas_get_coeff',0,
     &       BASIS_ERR)

        kprim = kprim + int_mb(k_nprim+icont-1)
      enddo
c
      if (.not. ma_push_get( mt_dbl, primprim, 'rexponents', 
     $                                               l_rexp, k_rexp))
     $  call errquit('waimfile : ma failed exp', primprim, MA_ERR)
      if (.not. ma_push_get( mt_dbl, primprim, 'rcoefficients', 
     $                                               l_rcoef, k_rcoef))
     $  call errquit('waimfile : ma failed exp', primprim, MA_ERR)
      if (.not. ma_push_get( mt_int, primprim, 'ricent',
     $                                               l_ricent,k_ricent))
     &  call errquit('waimfile : failed icent.',primprim,MA_ERR)
      if (.not. ma_push_get( mt_int, primprim, 'ritype',
     $                                               l_ritype,k_ritype))
     &  call errquit('waimfile : failed itype.',primprim,MA_ERR)
c
c
c
      k = 0
      l = 1
      m = 0
      do icont = 1, numcont
c> determine number type
        if (.not. bas_continfo(basis, icont, type,
     $                                        idum, ngeno,sphcart))
     $        call errquit('waimfile : error code = ',3,BASIS_ERR)
        n=1
        do i = 1 + k, int_mb(k_rprim+icont-1) + k
          do j = 1  , int_mb(k_nprim+icont-1) 
c> itype   1:numcount to ritype   1:primprim
            if (sphcart.eq.1) then
                int_mb(k_ritype+l-1) = sphebf(type,n)
            else
                int_mb(k_ritype+l-1) = cartbf(type,n)
            endif
c> exp     1:kprim    to rexp     1:primprim
            dbl_mb(k_rexp+l-1) = dbl_mb(k_exp+m+j-1)
c> coeff   1:kprim    to rcoeff   1:primprim
            dbl_mb(k_rcoef+l-1) = dbl_mb(k_coef+m+j-1)
c> icenter 1:nprim    to ricenter 1:primprim
            int_mb(k_ricent+l-1) =  int_mb(k_icent+icont-1)
            l = l + 1
          enddo
            n = n + 1
        enddo
        k = k + int_mb(k_tprim+icont-1) 
        m = m + int_mb(k_nprim+icont-1)
      enddo !icont
c
      if (.not. ma_push_get(mt_dbl,primprim*noccs,'coeffients', 
     $                                                      l_co, k_co))
     $    call errquit('waimfile : ma failed', primprim*noccs, MA_ERR)

      naco = 0
      do iset = 1, nset
        do n = 1, nocc(iset)
          call ga_get(g_vecs(iset),1, nbf,n,n,dbl_mb(k_vecs),1)
          k = 0
          l = 1
          m = 0
          do icont = 1, numcont
            do i = 1 + k, int_mb(k_rprim+icont-1) + k
              m = m + 1
              do j=1,int_mb(k_nprim+icont-1) 
c> expand eigenvector 1:nbf  to 1:primprim
                dbl_mb(k_co + primprim*naco + l-1) = 
     $                                      dbl_mb(k_vecs+m-1)
c> easy human reading
                if ( abs( dbl_mb(k_co + primprim*naco +l-1)).lt.1d-30)
     $            dbl_mb(k_co + primprim*naco +l-1) = 0.0
                l = l + 1
              enddo !j
            enddo  ! i
            k = k + int_mb(k_tprim+icont-1) 
          enddo !icont
          naco = naco + 1 
        enddo
      enddo

c> WFN files is always in catesian
         
      if (sphcart.eq.1) then
         if (.not. waimfile_sph2cart( primprim, 
     $                               noccs,
     $                               nocc,
     $                               nset,
     $                               k_ritype, 
     $                               k_icent, k_ricent,
     $                               k_exp, k_rexp,
     $                               k_coef,  k_rcoef,
     $                               k_co,
     $                               ks_rexp,ks_rcoef,
     $                               ks_ricent,ks_ritype,
     $                               ks_co,
     $                               sprimprim))
     $       call errquit('waimfile: write error wfn',0, GA_ERR)
      endif


cccccccccccccccccccccccccc VIRIAL THEOREM VALUES ccccccccccccccccccccc
c
c Recover potential energy (V) by evaluating kinetic energy (T)
c
c                        E_total = T + V                       
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c
      if (.not. rtdb_get(rtdb,theory(1:inp_strlen(theory))//':energy'
     $                                            , MT_DBL,1,etot))
     $  etot = 0.0d0




      if (etot.ne.0.0d0) then

c>  do kinetic energy
        call int_init(rtdb,1,basis)
c
        g_dens(1) = ga_create_atom_blocked(geom, basis, 'dens1')
c
        call ga_zero(g_dens(1))
        call ga_dgemm('n', 't', nbf, nbf, nocc(1), 1.0d0, g_vecs(1),
     $                                      g_vecs(1), 0.0d0, g_dens(1))

        g_hcore(1)  = ga_create_atom_blocked(geom, basis,'kin')
        call ga_zero(g_hcore(1))
        call int_1e_ga(basis, basis, g_hcore(1), 'kinetic', .false.)
c
        ekin = ga_ddot(g_dens(1),g_hcore(1)) 
c
        if (nset.eq.2) then
          g_dens(2) = ga_create_atom_blocked(geom, basis, 'dens2')
          call ga_zero(g_dens(2))
          call ga_dgemm('n', 't', nbf, nbf, nocc(2), 1.0d0, g_vecs(2),
     $        g_vecs(2), 0.0d0, g_dens(2))
          g_hcore(2)  = ga_create_atom_blocked(geom, basis,'kin')
          call ga_zero(g_hcore(2))
          call int_1e_ga(basis, basis, g_hcore(2), 'kinetic', .false.)
c
          ekin = ekin + ga_ddot(g_dens(2),g_hcore(2))
          if (.not. ga_destroy(g_hcore(2)))
     $       call errquit('waimfile: destroy g_hcore2',0, GA_ERR)
          if (.not. ga_destroy(g_dens(2)))
     $       call errquit('waimfile: destroy g_dens2',0, GA_ERR)
        else
          ekin = 2*ekin
        endif

        if (.not. ga_destroy(g_hcore(1)))
     $       call errquit('waimfile: destroy g_hcore1',0, GA_ERR)
        if (.not. ga_destroy(g_dens(1)))
     $       call errquit('waimfile: destroy g_dens1',0, GA_ERR)
c
        epot = etot - ekin
        virial =  abs(epot/ekin)
c
        if (ga_nodeid() .eq. 0) then
          write (6,*) '             Total Energy E = ', etot
          write (6,*) '           Kinetic Energy V = ', ekin
          write (6,*) '         Potential Energy T = ', epot
          write (6,*) '            Virial ABS(V/T) = ', virial
        endif
c
        call int_terminate()
      else
        if (ga_nodeid() .eq. 0) 
     $       write (6,*) ' Total Energy = Zero (NO FOUND)'
          virial = 0.0d0
      endif

      call ga_sync()

cccccccccccccccccccccccccc VIRIAL THEOREM VALUES end cccccccccccccccccc

c> wavefunction

      rhf='RHF'
      if (nset.gt.1) rhf = 'UHF'


ccccccccccccccccccccc   Print file   cccccccccccccccccccccccc
c>   Print .WFN file :  waimfile_printwfn
c>   Print .WFX file :  waimfile_printwfx

      if (ga_nodeid() .eq. 0) then
      if (sphcart.eq.1) then
      if(nowfx) then 
        call util_file_name('wfn', .false., .false., wfnfile)
        if (.not. waimfile_printwfn(wfnfile,title,natoms,nset,sprimprim,
     $              nocc,int_mb(ks_ricent),int_mb(ks_ritype),
     $              dbl_mb(ks_co), dbl_mb(ks_co+sprimprim*nocc(1)),
     $              dbl_mb(ks_rexp), dbl_mb(ks_rcoef),
     $              coord,qnuc,symbol,
     $              dbl_mb(k_occ(1)),dbl_mb(k_occ(2)),
     $              dbl_mb(k_eval(1)),dbl_mb(k_eval(2)),
     $              etot,virial,rhf))
     $       call errquit('waimfile: write error wfn',0, GA_ERR)
         else
         call util_file_name('wfx', .false., .false., wfnfile)
         if (.not. waimfile_printwfx( wfnfile,title,
     $              natoms,nset,sprimprim,
     $              nocc,int_mb(ks_ricent),int_mb(ks_ritype),
     $              dbl_mb(ks_co), dbl_mb(ks_co+sprimprim*nocc(1)),
     $              dbl_mb(ks_rexp), dbl_mb(ks_rcoef),
     $              coord,qnuc,symbol,
     $              dbl_mb(k_occ(1)),dbl_mb(k_occ(2)),
     $              dbl_mb(k_eval(1)),dbl_mb(k_eval(2)),
     $              etot,virial,rhf,total_charge))
     $       call errquit('waimfile: write error wfn',0, GA_ERR)
      endif !nowfx
      else
      if(nowfx) then 
        call util_file_name('wfn', .false., .false., wfnfile)
        if (.not. waimfile_printwfn( wfnfile,title,natoms,nset,primprim,
     $              nocc,int_mb(k_ricent),int_mb(k_ritype),
     $              dbl_mb(k_co), dbl_mb(k_co+primprim*nocc(1)),
     $              dbl_mb(k_rexp), dbl_mb(k_rcoef),
     $              coord,qnuc,symbol,
     $              dbl_mb(k_occ(1)),dbl_mb(k_occ(2)),
     $              dbl_mb(k_eval(1)),dbl_mb(k_eval(2)),
     $              etot,virial,rhf))
     $       call errquit('waimfile: write error wfn',0, GA_ERR)
         else
         call util_file_name('wfx', .false., .false., wfnfile)
         if (.not. waimfile_printwfx( wfnfile,title,
     $              natoms,nset,primprim,
     $              nocc,int_mb(k_ricent),int_mb(k_ritype),
     $              dbl_mb(k_co), dbl_mb(k_co+primprim*nocc(1)),
     $              dbl_mb(k_rexp), dbl_mb(k_rcoef),
     $              coord,qnuc,symbol,
     $              dbl_mb(k_occ(1)),dbl_mb(k_occ(2)),
     $              dbl_mb(k_eval(1)),dbl_mb(k_eval(2)),
     $              etot,virial,rhf,total_charge))
     $       call errquit('waimfile: write error wfn',0, GA_ERR)
      endif !nowfx
      endif !sphcart
      endif !ga.eq.0

ccccccccccccccccccccc   Print file END   cccccccccccccccccccc
c


      if (nset.gt.1) then
         if (.not. ga_destroy(g_vecs(2)))
     $       call errquit('waimfile: destroy g_movecs',0, GA_ERR)
      endif

      if (.not. ga_destroy(g_vecs(1)))
     $       call errquit('waimfile: destroy g_movecs',0, GA_ERR)

      if(.not.ma_chop_stack(l_occ(1)))
     &   call errquit('waimfile, ma_chop_stack of l_icent failed',911,
     &       ma_err)
c
c*            if (.not. geom_print(geom))
c*     $           call errquit('property: geom_print ?',0, GEOM_ERR)
c*            call util_flush(luout)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('waimfile:error destroying geom and basis handles',911,
     &       GEOM_ERR)
c
      return
      end

      logical function waimfile_printwfn(filename,title,natoms,
     $                                   nset,primprim,
     $                                   nocc,ricent,ritype,
     $                                   coa,cob,
     $                                   rexp,rcoef,coord,qnuc,symbol,
     $                                   occa,occb,evala,evalb,
     $                                   etot, virial,rhf)
      implicit none
      character*80 filename !in
      character*80 title !in
      integer natoms
      integer iset,nset
      integer primprim
      integer nocc(2)
      integer ricent(primprim), ritype(primprim)
      integer i,j
      integer naco
      double precision rexp(primprim),rcoef(primprim)
      double precision coord(3,natoms), qnuc(natoms)
      double precision coa(primprim,nocc(1)),cob(primprim,nocc(2))
      character*4 symbol(natoms)
      double precision occa(*), occb(*)
      double precision evala(*), evalb(*)
      double precision etot, virial
      character*16 rhf
c
      integer unitno
      parameter (unitno = 38)

      waimfile_printwfn=.false.
      write(6,*) '   Name of Final wnf-file  ', filename
cc
      open(unit=unitno, file=filename, status='unknown',
     &       form='formatted')

      write(unitno,1001) title(1:80)
      write(unitno,1002) nocc(1)+nocc(2),primprim,natoms
ccc write geometry
      write(unitno,1003) (symbol(i),i, i,
     $                      coord(1,i),coord(2,i),coord(3,i),
     $                      qnuc(i),i=1,natoms)
      write(unitno,1004) (ricent(i),i=1,primprim)
      write(unitno,1005) (ritype(i),i=1,primprim)
      write(unitno,1006) (rexp(i),  i=1,primprim)
c
c>>occupation and orb. ener.
c>>coeff X eigenvec
      naco=1
      do i = 1, nocc(1)
        write(unitno,1007) i,occa(i),evala(i)
        write(unitno,1008) (coa(j,naco) *rcoef(j),j=1, primprim) 
        naco = naco + 1
      enddo
      if(nset.eq.2) then 
        naco=1
        do i = 1, nocc(2)
           write(unitno,1007) i,occb(i),evalb(i)
           write(unitno,1008) (cob(j,naco) *rcoef(j),j=1, primprim) 
           naco = naco + 1
        enddo
      endif

c>>end message
      write(unitno,1009) 
c>> restricted or unretricted banner
      write(unitno,1010) rhf,etot, virial
c
      close(unit=unitno, status='keep')

      waimfile_printwfn=.true.

ccc>>>>> WFN Format  <<<<<<<<<<<<<<<<<<<
 1001 format (1A80)
 1002 format ('GAUSSIAN',10X,I5,' MOL ORBITALS',1X,I6,' PRIMITIVES',
     &                                                4X,I5,' NUCLEI')
 1003 format (A4,I4,4X,'(CENTRE',I3,')',1X,3F12.8,'  CHARGE =',F5.1)
 1004 format ('CENTRE ASSIGNMENTS',2X,20I3)
 1005 format ('TYPE ASSIGNMENTS',4X,20I3)
 1006 format ('EXPONENTS',1X,1P,5E14.7)
 1007 format ('MO',I3,21X,'OCC NO = ',F12.8,' ORB. ENERGY =',F13.8 )
 1008 format (1P,5E16.8)
 1009 format ('END DATA')
 1010 format (A8,' ENERGY =',F20.10,'   VIRIAL(-V/T)  =',F13.8)
ccc>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      end
      logical function waimfile_printwfx(filename,title,natoms,
     $                                   nset,primprim,
     $                                   nocc,ricent,ritype,
     $                                   coa,cob,
     $                                   rexp,rcoef,coord,qnuc,symbol,
     $                                   occa,occb,evala,evalb,
     $                                   etot, virial,rhf,total_charge)
      implicit none
      character*80 filename !in
      character*80 title !in
      integer natoms
      integer iset,nset
      integer primprim
      integer nocc(2)
      integer ricent(primprim), ritype(primprim)
      integer i,j
      integer naco
      double precision rexp(primprim),rcoef(primprim)
      double precision coord(3,natoms), qnuc(natoms)
      double precision coa(primprim,nocc(1)),cob(primprim,nocc(2))
      character*4 symbol(natoms)
      character*30 string_var
      double precision occa(*), occb(*)
      double precision evala(*), evalb(*)
      double precision etot, virial
      double precision total_charge
      character*16 rhf
c
      integer unitno
      parameter (unitno = 38)

      waimfile_printwfx=.false.


      write(6,*) '   Name of Final wnf-file  ', filename
c
      open(unit=unitno, file=filename, status='unknown',
     &     form='formatted')

      write(unitno,2001) title(1:80)
      write(unitno,2004) ' GTO' 
      write(unitno,2007) natoms
      write(unitno,2010) nocc(1)+nocc(2)
      write(unitno,2013) 0
      if (nset.eq.1) then
          total_charge=total_charge-2*nocc(1)        
          write(unitno,2016) total_charge
      else
          total_charge=total_charge-(nocc(1)+nocc(2))
          write(unitno,2016) total_charge
      endif
      if (nset.eq.1) then
          write(unitno,2019) 2*nocc(1)
          write(unitno,2022) nocc(1)
          write(unitno,2025) nocc(1)
          i = 1
          write(unitno,2028) i  
      else
          write(unitno,2019) nocc(1)+nocc(2)
          write(unitno,2022) nocc(1)
          write(unitno,2025) nocc(2)
          write(unitno,2028) (nocc(1)-nocc(2))+1
      endif 
      write(unitno,2031) 0

c write geometry
      write(unitno,2034)
      do i=1,natoms
        write(string_var, '(I30)') i
        write(unitno,2035) symbol(i),adjustl(string_var)
      enddo
      write(unitno,2036)

      write(unitno,2037)
      write(unitno,2038) (int(qnuc(i)),i=1,natoms)
      write(unitno,2039)

      write(unitno,2040)
      write(unitno,2041) (qnuc(i),i=1,natoms)
      write(unitno,2042)
      
      write(unitno,2043)
      write(unitno,2044) (coord(1,i),coord(2,i),coord(3,i),
     $                                            i=1,natoms)

      write(unitno,2045)
      write(unitno,2046)
      write(unitno,2047) primprim
      write(unitno,2048)
      write(unitno,2049)
      write(unitno,2050) (ricent(i),i=1,primprim)
      write(unitno,2051)
      write(unitno,2052)
      write(unitno,2050) (ritype(i),i=1,primprim)
      write(unitno,2053)
      write(unitno,2054)
      write(unitno,2055) (rexp(i),  i=1,primprim)
      write(unitno,2056)


c>>Occupation Numbers
      write(unitno,2057)    
      write(unitno,2041) (occa(i),i=1,nocc(1))
      if(nset.eq.2) write(unitno,2041) (occb(i),i=1,nocc(2))
      write(unitno,2058)
c Orbital Energies
      write(unitno,2059)
      write(unitno,2041) (evala(i),i=1,nocc(1))
      if(nset.eq.2) write(unitno,2041) (evalb(i),i=1,nocc(2))
      write(unitno,2060)
c Orbital types
      write(unitno,2061)
      if (nset.eq.1) then
         do i = 1,nocc(1)
         write(unitno,2062)
           enddo
      else
        do i = 1,nocc(1)
         write(unitno,2063)
        enddo
        do i = 1,nocc(2)
         write(unitno,2064)
        enddo 
      endif
      write(unitno,2065)   

      write(unitno,2066)
      naco = 1
      do i = 1, nocc(1)
         write(unitno,2067)
         write(unitno,2047) i
         write(unitno,2068)
         write(unitno,2069) (coa(j,naco) *rcoef(j),j=1, primprim) 
         naco = naco + 1
      enddo
      if (nset.eq.2) then 
        do i = 1, nocc(1)
           write(unitno,2067)
           write(unitno,2047) i+nocc(1)
           write(unitno,2068)
           write(unitno,2069) (coa(j,naco) *rcoef(j),j=1, primprim) 
           naco = naco + 1
        enddo
      endif
      write(unitno,2070)
cEnergy print
      write(unitno,2071)
      write(unitno,2041) etot
      write(unitno,2072)
cVirial print
      write(unitno,2073)
      write(unitno,2041) virial
      write(unitno,2074)

      if (nset.eq.1) then 
        write(unitno,2075) 'Restricted SCF'
      else
        write(unitno,2075) 'Unrestricted SCF'
      endif

      close(unit=unitno, status='keep')

      waimfile_printwfx=.true.



c>>>> WFX Format  <<<<<<<<<<<<<<<<

 2001 format ('<Title>',/,1A80,/,'</Title>',/)

 2004 format ('<Keywords>',/,1A5,/,'</Keywords>',/)

 2007 format ('<Number of Nuclei>',/,I5,/,'</Number of Nuclei>',/)

 2010 format ('<Number of Occupied Molecular Orbitals>',/,I5,/,
     $        '</Number of Occupied Molecular Orbitals>',/)

 2013 format ('<Number of Perturbations>',/,I2,/,
     $       ' </Number of Perturbations>',/)

 2016 format ('<Net Charge>',/,F5.2,/,'</Net Charge>',/)

 2019 format ('<Number of Electrons>',/,I3,/,'</Number of Electrons>',/)

 2022 format ('<Number of Alpha Electrons>',/,I5,/,
     $        '</Number of Alpha Electrons>',/)

 2025 format ('<Number of Beta Electrons>',/,I5,/,
     $        '</Number of Beta Electrons>',/)

 2028 format ('<Electronic Spin Multiplicity>',/,I5,/,
     $         '</Electronic Spin Multiplicity>',/)

 2031 format ('<Number of Core Electrons>',/,I5,/,
     $        '</Number of Core Electrons>',/)

 2034 format ('<Nuclear Names>')
 2035 format (A3,A30)
 2036 format ('</Nuclear Names>',/)

 2037 format ('<Atomic Numbers>')
 2038 format (I3)
 2039 format ('</Atomic Numbers>',/)

 2040 format ('<Nuclear Charges>')
 2041 format (E20.12)
 2042 format ('</Nuclear Charges>',/)

 2043 format ('<Nuclear Cartesian Coordinates>')
 2044 format (1X,3E20.12)
 2045 format ('</Nuclear Cartesian Coordinates>',/)

 2046 format ('<Number of Primitives>')
 2047 format (I5)
 2048 format ('</Number of Primitives>',/)

 2049 format ('<Primitive Centers>')
 2050 format (5I20)
 2051 format ('</Primitive Centers>',/)

 2052 format ('<Primitive Types>')
 2053 format ('</Primitive Types>',/)

 2054 format ('<Primitive Exponents>')
 2055 format (5E20.12)
 2056 format ('</Primitive Exponents>',/)

 2057 format ('<Molecular Orbital Occupation Numbers>')
 2058 format ('</Molecular Orbital Occupation Numbers>',/)

 2059 format ('<Molecular Orbital Energies>')
 2060 format ('</Molecular Orbital Energies>',/)

 2061 format ('<Molecular Orbital Spin Types>')
 2062 format (' Alpha and Beta')
 2063 format (' Alpha')
 2064 format (' Beta')
 2065 format ('</Molecular Orbital Spin Types>',/)

 2066 format ('<Molecular Orbital Primitive Coefficients>')
 2067 format ('<MO Number>')
 2068 format ('</MO Number>',/)
 2069 format (4E20.12)
 2070 format ('</Molecular Orbital Primitive Coefficients>',/)

 2071 format ('<Energy = T + Vne + Vee + Vnn>')
 2072 format ('</Energy = T + Vne + Vee + Vnn>',/)

 2073 format ('<Virial Ratio (-V/T)>')
 2074 format ('</Virial Ratio (-V/T)>',/)

 2075 format ('<Model>',/,1A80,/,'</Model>',/)

      end
      logical function waimfile_sph2cart( primprim, 
     $                                   noccs,
     $                                   nocc,
     $                                   nset,
     $                                   k_ritype, 
     $                                   k_icent, k_ricent,
     $                                   k_exp, k_rexp,
     $                                   k_coef,  k_rcoef,
     $                                   k_co,
     $                                   ks_rexp,ks_rcoef,
     $                                   ks_ricent,ks_ritype,
     $                                   ks_co,
     $                                   sprimprim) 
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer primprim
      integer sprimprim
      integer nset,iset
      integer numd,numf,numg,numh
      integer  noccs, nocc(2)
      integer  k_icent, k_ricent
      integer k_ritype
      integer k_exp, k_rexp
      integer k_coef,  k_rcoef 
      integer k_co 
      integer ls_co, ks_co
      integer ls_ritype, ks_ritype
      integer ls_ricent, ks_ricent
      integer ls_rexp, ks_rexp
      integer ls_rcoef, ks_rcoef
      integer i,j,k,l,n,r,h,m,t
      integer np,iorb,forb,nsf
      integer naco
      double precision sqrt3, sqrt5, sqrt7, DUMMY
c:
c: author A. Guevara 7/2013
c:
      waimfile_sph2cart=.false.

      !sqrt3=1.7320508076D+00
      sqrt3=3.0**(1.0/3.0)
        
      numd=0
      numf=0
      numg=0
      numh=0
      do i = 0,primprim-1
         if (int_mb(k_ritype+i).eq.5) numd=numd+1
         if (int_mb(k_ritype+i).eq.10) numf=numf+1
         if (int_mb(k_ritype+i).eq.17) numg=numg+1
         if (int_mb(k_ritype+i).eq.26) numh=numh+1
      enddo
      sprimprim=primprim+numd+3*numf+6*numg+10*numh 
      if (.not. ma_push_get(mt_dbl,sprimprim,'rexponents',
     $                                       ls_rexp, ks_rexp))
     $ call errquit('waimfile_sph2cart : ma failed exp', 
     $                                         sprimprim, MA_ERR)
      if (.not. ma_push_get(mt_dbl,sprimprim,'rcoefficients',
     $                                     ls_rcoef, ks_rcoef))
     $call errquit('waimfile_sph2cart : ma failed exp',
     $                                         sprimprim, MA_ERR)
      if (.not. ma_push_get(mt_int,sprimprim,'ricent',
     $                                     ls_ricent,ks_ricent))
     &call errquit('waimfile_sph2cart : failed icent.',sprimprim,MA_ERR)
      if (.not. ma_push_get(mt_int,sprimprim,'ritype',
     $                                      ls_ritype,ks_ritype))
     & call errquit('waimfile_sph2cart : failed itype.',
     $                                         sprimprim, MA_ERR)


      if (.not. ma_push_get(mt_dbl,sprimprim*noccs,'coeffients',
     $                                                 ls_co, ks_co))
     $    call errquit('waimfile : ma failed ls coeff', 
     $                                        sprimprim*noccs, MA_ERR)

      i=0
      t=0
      do while (i<primprim)
        if ( (int_mb(k_ritype +i).eq.5)  .or.
     $       (int_mb(k_ritype +i).eq.10) .or. 
     $       (int_mb(k_ritype +i).eq.17) .or.
     $       (int_mb(k_ritype +i).eq.26)) then

           if (int_mb(k_ritype +i).eq.5) then
              iorb=5
              forb=10
              nsf=5
           elseif (int_mb(k_ritype +i).eq.10) then
              iorb=11
              forb=20
              nsf=7
           else if (int_mb(k_ritype +i).eq.17) then
              iorb=21
              forb=35
              nsf=9
           else if (int_mb(k_ritype +i).eq.26) then
              iorb=36
              forb=56
              nsf=11
           endif 
           j=i+1
           np=1
           do while (int_mb(k_ritype +i).eq.int_mb(k_ritype +j))
              np=np+1
              j=j+1
           enddo
           do l=iorb,forb
              k=0
              do while (k<np)
                 int_mb(ks_ritype +t+k)=l
                 int_mb(ks_ricent +t+k)=int_mb(k_ricent +i)
                 dbl_mb(ks_rexp+t+k)=dbl_mb(k_rexp +i+k)
                 dbl_mb(ks_rcoef+t+k)=dbl_mb(k_rcoef+i+k)
                 k=k+1
              enddo
              t=np+t
           enddo
           i=i+nsf*np
        else
          int_mb(ks_ritype +t)=int_mb(k_ritype +i)
          int_mb(ks_ricent +t)=int_mb(k_ricent +i)
          dbl_mb(ks_rexp+t)=dbl_mb(k_rexp +i)
          dbl_mb(ks_rcoef+t)=dbl_mb(k_rcoef+i)          
          i=i+1            
          t=t+1
        endif!k_ritype
      enddo !primprim



      naco = 0
      do iset = 1,nset
      do m = 1, nocc(iset)
      i=0
      t=0
      do while (i<primprim)
          if ( (int_mb(k_ritype +i).eq.5)  .or.
     $         (int_mb(k_ritype +i).eq.10) .or.
     $         (int_mb(k_ritype +i).eq.17) .or.
     $         (int_mb(k_ritype +i).eq.26)) then
c d functions
            if (int_mb(k_ritype +i).eq.5) then
               iorb=1
               forb=6
               nsf=5
                j=i+1
                np=1
                do while (int_mb(k_ritype +i).eq.int_mb(k_ritype +j))
                   np=np+1
                   j=j+1
                enddo 
                do l = iorb,forb 
                   k=0
                   do while (k<np)
                      if (l.eq.1) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                  -0.5*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           +0.5d0*sqrt(3d0)*dbl_mb(k_co+primprim*naco+i+k+4*np)
                      elseif (l.eq.2) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                  -0.5*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           -0.5d0*sqrt(3d0)*dbl_mb(k_co+primprim*naco+i+k+4*np)
                      elseif (l.eq.3) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                    dbl_mb(k_co+primprim*naco+i+k+2*np)
                      elseif (l.eq.4) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                    sqrt3*dbl_mb(k_co+primprim*naco+i+k)
                      elseif (l.eq.5) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                   -sqrt3*dbl_mb(k_co+primprim*naco+i+k+3*np)
                      elseif (l.eq.6) then
                        dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                    sqrt3*dbl_mb(k_co+primprim*naco+i+k+np)
                      endif
                      k=k+1
                   enddo
                   t=np+t
                enddo
                i=i+nsf*np
c d functions
c f functions
            else if (int_mb(k_ritype +i).eq.10) then
               iorb=1
               forb=10
               nsf=7
                j=i+1
                np=1
                do while (int_mb(k_ritype +i).eq.int_mb(k_ritype +j))
                  np=np+1
                  j=j+1
                enddo
                do l = iorb,forb
                  k=0
                  do while (k<np)
                    if (l.eq.1) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           6.1237243570D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $           -7.9056941504D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
                    else if (l.eq.2) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -7.9056941504D-01*dbl_mb(k_co+primprim*naco+i+k)
     $           -6.1237243570D-01*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    else if (l.eq.3) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $                dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.4) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k)
     $           -6.1237243570D-01*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    else if (l.eq.5) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -1.5000000000D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
     $           +1.9364916731D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
                    else if (l.eq.6) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -1.5000000000D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
     $           -1.9364916731D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
                    else if (l.eq.7) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           6.1237243570D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $           +2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
                    else if (l.eq.8) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -2.4494897428D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.9) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           2.4494897428D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    else if (l.eq.10) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           3.8729833462D+00*dbl_mb(k_co+primprim*naco+i+k+1*np)   
                    endif
                    k=k+1
                  enddo!k
                  t=np+t
                enddo !l
                i=i+nsf*np
c f functions
c g funtions
            else if (int_mb(k_ritype +i).eq.17) then
               iorb=1
               forb=15
               nsf=9
               j=i+1
               np=1
               do while (int_mb(k_ritype +i).eq.int_mb(k_ritype +j))
                  np=np+1
                  j=j+1
               enddo
               do l = iorb,forb
                 k=0
                 do while (k<np)
                    if (l.eq.1) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         3.7500000000D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $         -5.5901699437D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $         +7.3950997289D-01*dbl_mb(k_co+primprim*naco+i+k+8*np)
                    else if (l.eq.2) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         3.7500000000D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $         +5.5901699437D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $         +7.3950997289D-01*dbl_mb(k_co+primprim*naco+i+k+8*np)
                    else if (l.eq.3) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $              dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.4) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         2.9580398915D+00*dbl_mb(k_co+primprim*naco+i+k)
     $         -1.1180339887D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    else if (l.eq.5) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $         -2.0916500663D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
                    else if (l.eq.6) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         -2.9580398915D+00*dbl_mb(k_co+primprim*naco+i+k)
     $         -1.1180339887D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    else if (l.eq.7) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         -2.0916500663D+00*dbl_mb(k_co+primprim*naco+i+k+np)
     $         -2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.8) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         -3.1622776602D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
                    else if (l.eq.9) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         3.1622776602D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.10) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         7.5000000000D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $         -4.4370598373D+00*dbl_mb(k_co+primprim*naco+i+k+8*np)
                    else if (l.eq.11) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         -3.0000000000D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $         +3.3541019662D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
                    else if (l.eq.12) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         -3.0000000000D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
     $         -3.3541019662D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
                    else if (l.eq.13) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         6.2749501990D+00*dbl_mb(k_co+primprim*naco+i+k+1*np)
     $         -2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.14) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         2.3717082451D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $         +6.2749501990D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
                    else if (l.eq.15) then
                    dbl_mb(ks_co+sprimprim*naco+t+k)=
     $         6.7082039325D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
                    endif
                    k=k+1
                 enddo
                 t=np+t
               enddo
               i=i+nsf*np
c g functions
c h functions
            else if (int_mb(k_ritype +i).eq.26) then
               iorb=1
               forb=21
               nsf=11
               j=i+1
               np=1
               do while (int_mb(k_ritype +i).eq.int_mb(k_ritype +j))
                 np=np+1
                 j=j+1
               enddo
               do l = iorb,forb
                 k=0
                 do while (k<np)
                   if (l.eq.1) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           dbl_mb(k_co+primprim*naco+i+k+5*np)
                    else if (l.eq.2) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           3.8729833462D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.3) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -5.0000000000D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $           -5.1234753830D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
                    else if (l.eq.4) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -4.1833001327D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           -5.8094750193D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.5) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           1.8750000000D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $           +2.5617376915D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
     $           +2.2185299187D+00*dbl_mb(k_co+primprim*naco+i+k+9*np)
                    else if (l.eq.6) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           7.0156076002D-01*dbl_mb(k_co+primprim*naco+i+k)
     $           +5.2291251658D-01*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           +4.8412291828D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.7) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -3.8729833462D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
                    else if (l.eq.8) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           1.0246950766D+01*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.9) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           5.8094750193D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $           +1.2549900398D+01*dbl_mb(k_co+primprim*naco+i+k+8*np)
                    else if (l.eq.10) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -8.8741196746D+00*dbl_mb(k_co+primprim*naco+i+k+1*np)
     $           -5.1234753830D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.11) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -4.8412291828D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $           -1.5687375498D+00*dbl_mb(k_co+primprim*naco+i+k+8*np)
     $           -3.5078038001D+00*dbl_mb(k_co+primprim*naco+i+k+10*np)
                    else if (l.eq.12) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -5.0000000000D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $           +5.1234753830D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
                    else if (l.eq.13) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           1.2549900398D+01*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           -5.8094750193D+00*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.14) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           3.7500000000D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $           -1.3311179512D+01*dbl_mb(k_co+primprim*naco+i+k+9*np)
                    else if (l.eq.15) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -7.0156076002D+00*dbl_mb(k_co+primprim*naco+i+k)
     $           -1.0458250332D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           +9.6824583655D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.16) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           5.8094750193D+00*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $           -4.1833001327D+00*dbl_mb(k_co+primprim*naco+i+k+8*np)
                    else if (l.eq.17) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           8.8741196746D+00*dbl_mb(k_co+primprim*naco+i+k+1*np)
     $           -5.1234753830D+00*dbl_mb(k_co+primprim*naco+i+k+3*np)
                    else if (l.eq.18) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -9.6824583655D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $           -1.0458250332D+00*dbl_mb(k_co+primprim*naco+i+k+8*np)
     $           +7.0156076002D+00*dbl_mb(k_co+primprim*naco+i+k+10*np)
                    else if (l.eq.19) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           1.8750000000D+00*dbl_mb(k_co+primprim*naco+i+k+5*np)
     $           -2.5617376915D+00*dbl_mb(k_co+primprim*naco+i+k+7*np)
     $           +2.2185299187D+00*dbl_mb(k_co+primprim*naco+i+k+9*np)
                    else if (l.eq.20) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           3.5078038001D+00*dbl_mb(k_co+primprim*naco+i+k)
     $           -1.5687375498D+00*dbl_mb(k_co+primprim*naco+i+k+2*np)
     $           +4.8412291828D-01*dbl_mb(k_co+primprim*naco+i+k+4*np)
                    else if (l.eq.21) then
                      dbl_mb(ks_co+sprimprim*naco+t+k)=
     $           -4.8412291828D-01*dbl_mb(k_co+primprim*naco+i+k+6*np)
     $           +5.2291251658D-01*dbl_mb(k_co+primprim*naco+i+k+8*np)
     $           -7.0156076002D-01*dbl_mb(k_co+primprim*naco+i+k+10*np)
                    endif
                    k=k+1
                  enddo
                  t=np+t
                enddo
                i=i+nsf*np
c h functions
            endif
            else
            dbl_mb(ks_co+sprimprim*naco+t)=
     $                   dbl_mb(k_co+primprim*naco+i)
            i=i+1
            t=t+1
            endif!k
        enddo !i
        naco = naco + 1
      enddo !m
      enddo !iset

      waimfile_sph2cart=.true.

      end
cgraveyardc      character*3 atmsymbol(nelements)
cgraveyardc      data atmsymbol/'  H',' He',' Li','  Be','  B','  C',
cgraveyardc     *               '  N','  O','  F','  Ne','  Na','  Mg',
cgraveyardc     *               ' Al',' Si','  P','  S','  Cl',' Ar',
cgraveyardc     *               '  K',' Ca',' Sc',' Ti','  V',' Cr',
cgraveyardc     *               ' Mn',' Fe',' Co',' Ni',' Cu',' Zn',
cgraveyardc     *               ' Ga',' Ge',' As',' Se',' Br',' Kr',
cgraveyardc     *               ' Rb',' Sr','  Y',' Zr',' Nb',' Mo',
cgraveyardc     *               ' Tc',' Ru',' Rh',' Pd',' Ag',' Cd',
cgraveyardc     *               ' In',' Sn',' Sb',' Te','  I',' Xe',
cgraveyardc     *               ' Cs',' Ba',' La',' Ce',' Pr',' Nd',
cgraveyardc     *               ' Pm',' Sm',' Eu',' Gd',' Tb',' Dy',
cgraveyardc     *               ' Ho',' Er',' Tm',' Yb',' Lu',' Hf',
cgraveyardc     *               ' Ta','  W',' Re',' Os',' Ir',' Pt',
cgraveyardc     *               ' Au',' Hg',' Tl',' Pb',' Bi',' Po',
cgraveyardc     *               ' At',' Rn',' Fr',' Ra',' Ac',' Th',
cgraveyardc     *               ' P','  U',' Np',' Pu',' Am','  Cm',
cgraveyardc     *               ' Bk',' Cf',' Es',' Fm',' Md',' No',
cgraveyardc     *               ' Lr',' Rf',' Db',' Sg',' Bh',' Hs', 
cgraveyardc     *               ' MT',' Ds',' Rg'/
cgraveyard
cgraveyardc
cgraveyardc itypes
cgraveyardc     /'l','s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'/
cgraveyardc      -1   0    1    2    3    4    5
c $Id$
