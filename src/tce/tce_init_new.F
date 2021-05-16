      subroutine tce_init_new(rtdb)
c
c $Id$
c
c     Initialize for many-electron calculations.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "dra.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
c
      integer rtdb             ! RTDB handle
      integer nbf_verify       ! nbf read from movecs_read_header
      integer ipol_verify      ! ipol read from movecs_read_header
      character*255 movecs_out ! MO vector file name
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
      character*4 hftype
      character*3 onoff1,onoff2
      integer i,j,n
      integer nclosed, nopen
      double precision a
      character*60 ctype
      character*60 cperturbative
      integer itol2e
      logical nodezero
      logical dummy
      logical rohf
      logical int_normalize
      external int_normalize
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical xc_gotxc
      external xc_gotxc
      character*10 module
      character*8 multchar(8)
      data multchar/'singlet ','doublet ','triplet ','quartet ',
     1              'quintet ','sextet  ','septet  ','octet   '/
       character*26 date
       integer ma_type
c
c
      nodezero=(ga_nodeid().eq.0)
c
c     -----------------------
c     Initialize stop watches
c     -----------------------
c
      do i = 1,max_clocks
        cpusecs(i) = 0.0d0
      enddo
c
c     ---------------------
c     Geometry and symmetry
c     ---------------------
c
      bytes = ma_sizeof(mt_dbl,1,mt_byte)
      if (.not.geom_create(geom,'geometry'))
     1  call errquit('tce_init_new: geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('tce_init_new: no geometry found',0,GEOM_ERR)
      if (.not.rtdb_get(rtdb,'dft:skeleton',mt_log,1,oskel)) then
        oskel=sym_number_ops(geom).gt.0
      endif
      if (.not.rtdb_get(rtdb,'dft:adapt',mt_log,1,oadapt)) then
        oadapt=sym_number_ops(geom).gt.0
      endif
      if (.not.sym_abelian_group(geom))
     1  call errquit('tce_init_new: non-Abelian symmetry not
     2                permitted',0,GEOM_ERR)
c
c     ---------
c     Basis set
c     ---------
c
      if (.not.bas_create(ao_bas_han,'ao basis'))
     1  call errquit('tce_init_new: bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,ao_bas_han,'ao basis'))
     1  call errquit('tce_init_new: no ao basis found',0,BASIS_ERR)
      if (.not.int_normalize(rtdb,ao_bas_han))
     1  call errquit('tce_init_new: int_normalize failed',0,INT_ERR)
      if (.not.bas_numbf(ao_bas_han,nbf))
     1  call errquit('tce_init_new: failed to read nbf',0,BASIS_ERR)
      if (.not.bas_numcont(ao_bas_han,nshells))
     1  call errquit('tce_init_new: failed to read nshells',0,BASIS_ERR)
c
c     ----------
c     MO vectors
c     ----------
c
      if (reference.eq.0) then
        if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_out))
     1    call errquit('tce_init_new: no final MO found',0,RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     1    call errquit('tce_init_new: no final MO found',0,RTDB_ERR)
      endif
c
c     ------------------------------------------------------
c     Number of electrons and orbitals and spin multiplicity
c     ------------------------------------------------------
c
      if (reference.eq.0) then
        rohf = .false.
        if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     1    call errquit('tce_init_new: failed to read ipol',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc))
     1    call errquit('tce_init_new: failed to read nocc',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'dft:mult',mt_int,1,mult))
     1    call errquit('tce_init_new: failed to read mult',0,RTDB_ERR)
        if (ipol.eq.1) then
          restricted = .true.
        else
          restricted = .false.
        endif
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nclosed))
     1    call errquit('tce_init_new: failed to read
     2                  nclosed',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'scf:nopen',mt_int,1,nopen))
     1    call errquit('tce_init_new: failed to read nopen',0,RTDB_ERR)
        nocc(1) = nclosed + nopen
        nocc(2) = nclosed
        mult = nopen + 1
        if(.not.rtdb_cget(rtdb,'scf:scftype',1,hftype))
     1    call errquit('tce_init_new: failed to read
     2                  scftype',0,RTDB_ERR)
        if (hftype.eq.'RHF') then
          restricted = .true.
          rohf = .false.
          ipol = 1
        else if (hftype.eq.'UHF') then
          restricted = .false.
          rohf = .false.
          ipol = 2
        else if(hftype.eq.'ROHF')then
          restricted = .false.
          rohf = .true.
          ipol = 2
        endif
      endif
c     ---------
c     TCE input
c     ---------
c
      if(.not. geom_freeze(rtdb, geom, 'tce', nfc(1))) then
         call errquit('tce_init_new: failed to get number of cores',0,
     1   GEOM_ERR)
      endif
      if (ipol.eq.2) nfc(2) = nfc(1)
      if (.not.rtdb_get(rtdb,'tce:frozen virtual',
     1  mt_int,1,nfv(1))) nfv(1) = 0
      if (ipol.eq.2) nfv(2) = nfv(1)
      if (.not.rtdb_cget(rtdb,'tce:module',1,module)) then
        module='tce'
      endif
      if (.not.rtdb_cget(rtdb,'tce:model2e',1,model2e))
     1  call errquit('tce_init_new: model2e not  specified',0,RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'tce:model',1,model))
     1  call errquit('tce_init_new: no theory specified',0,RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'tce:perturbative',1,perturbative)) then
        perturbative = 'none'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1  call errquit('tce_init_new: failed to write perturbative',0,
     2  RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'tce:ccsdvar',1,ccsd_var)) then
        ccsd_var = 'none'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1  call errquit('tce_init_new: failed to write ccsdvar',0,
     2  RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:thresh',mt_dbl,1,thresh)) then
        thresh = 1.0d-6
        if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_init_new: failed to write thresh',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:zlshift',mt_dbl,1,zlshift)) then
        zlshift = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:zlshift',mt_dbl,1,zlshift))
     1  call errquit('tce_init_new: failed to write zlshift',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:maxiter',mt_int,1,maxiter)) then
        maxiter = 100
        if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_init_new: failed to write maxiter',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ioalg',mt_int,1,ioalg)) then
        ioalg = 2
        if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_init_new: failed to write ioalg',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diis',mt_int,1,diis)) then
        diis = 5
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init_new: failed to write diis',0,RTDB_ERR)
      endif
      if (diis .gt. maxdiis) then
        diis = maxdiis
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init_new: failed to write diis',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:nroots',mt_int,1,nroots)) then
        nroots = 0
        if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1  call errquit('tce_init_new: failed to write nroots',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:hbard',mt_int,1,hbard)) then
        hbard = 500
        if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1  call errquit('tce_init_new: failed to write hbard',0,RTDB_ERR)
      endif
c 4index trans.
      if (.not.rtdb_get(rtdb,'tce:maxs',mt_int,1,imaxsize)) then
        imaxsize = 30
        if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,imaxsize))
     1  call errquit('tce_init_new: failed to write maxs',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ichopx',mt_int,1,ichop)) then
        ichop = 1
        if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichop))
     1  call errquit('tce_init_new: failed to write ichop',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:i4im',mt_int,1,fast2e)) then
        fast2e = 1
        if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,fast2e))
     1  call errquit('tce_init_new: failed to write i4im',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:idiskx',mt_int,1,idisk)) then
        idisk = 0
        if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idisk))
     1  call errquit('tce_init_new: failed to write idiskx',0,RTDB_ERR)
      endif
c CC2 and LR-CCSD
      if (.not.rtdb_cget(rtdb,'tce:ccsdvar',1,ccsd_var)) then
        ccsd_var='none'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1  call errquit('tce_init_new: failed to write ccsdvar',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:target',mt_int,1,target)) then
         target=1
        if (.not.rtdb_put(rtdb,'tce:target',mt_int,1,target))
     1  call errquit('tce_init_new: failed to write target',0,RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'tce:targetsym',1,targetsym)) then
        targetsym='none'
        if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1  call errquit('tce_init_new: failed to write
     2                targetsym',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:symmetry',mt_log,1,symmetry)) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1  call errquit('tce_init_new: failed to write
     2                symmetry',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:left',mt_log,1,left)) then
        left=.false.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_init_new: failed to write left',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:multipole',mt_int,1,multipole)) then
        multipole=0
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1  call errquit('tce_init_new: failed to write
     2                multipole',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diagnt1',mt_log,1,diagnt1)) then
        diagnt1=.false.
        if (.not.rtdb_put(rtdb,'tce:diagnt1',mt_log,1,diagnt1))
     1  call errquit('tce_init_new: failed to write diagnt1',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ampnorms',mt_log,1,ampnorms)) then
        ampnorms=.false.
        if (.not.rtdb_put(rtdb,'tce:ampnorms',mt_log,1,ampnorms))
     1  call errquit('tce_init_new: failed to write
     2                ampnorms',0,RTDB_ERR)
      endif
c
c    Beginning of response property input section
c
      if (.not.rtdb_get(rtdb,'tce:lineresp',mt_log,1,lineresp)) then
        lineresp=.false.
        if (.not.rtdb_put(rtdb,'tce:lineresp',mt_log,1,lineresp))
     1  call errquit('tce_init_new: failed to write
     2                lineresp',0,RTDB_ERR)
      endif
      if (lineresp) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_init_new: failed to write left',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:leftresp',mt_log,1,leftresp)) then
        leftresp=.false.
        if (.not.rtdb_put(rtdb,'tce:leftresp',mt_log,1,leftresp))
     1  call errquit('tce_init_new: failed to write
     2                leftresp',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:quadresp',mt_log,1,quadresp)) then
        quadresp=.false.
        if (.not.rtdb_put(rtdb,'tce:quadresp',mt_log,1,quadresp))
     1  call errquit('tce_init_new: failed to write
     2                quadresp',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get_info(rtdb,'tce:afreq',ma_type,anumfreq,
     1  date)) then
        anumfreq = 1
        if (.not.rtdb_put(rtdb,'tce:anumfreq',mt_int,1,anumfreq))
     1   call errquit('tce_init_new: failed to write
     2                anumfreq',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:afreq',mt_dbl,anumfreq,
     1  afreq)) then
        afreq(1) = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:afreq',mt_dbl,anumfreq,afreq))
     1    call errquit('tce_init_new: failed to write afreq',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get_info(rtdb,'tce:bfreq',ma_type,bnumfreq,
     1  date)) then
        bnumfreq = 1
        if (.not.rtdb_put(rtdb,'tce:bnumfreq',mt_int,1,bnumfreq))
     1   call errquit('tce_init_new: failed to write
     2                bnumfreq',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:bfreq',mt_dbl,bnumfreq,
     1  bfreq)) then
        bfreq(1) = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:bfreq',mt_dbl,bnumfreq,bfreq))
     1    call errquit('tce_init_new: failed to write bfreq',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:respaxis',mt_log,3,respaxis)) then
        respaxis(1) = .true.
        respaxis(2) = .true.
        respaxis(3) = .true.
        if (.not.rtdb_put(rtdb,'tce:respaxis',mt_log,3,respaxis))
     1    call errquit('tce_init_new: failed to write
     2                respaxis',0,RTDB_ERR)
      endif
c
c     End of response property input section
c
      if ((symmetry).and.(targetsym.eq.'none')) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1    call errquit('tce_init_new: failed to write
     2                symmetry',0,RTDB_ERR)
      endif
c
c     ----------------------
c     2-e integral tolerance
c     ----------------------
c
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:itol2e',mt_int,1,itol2e))
     1    call errquit('tce_init_new: failed to read itol2e',0,RTDB_ERR)
        tol2e=10.0d0**(-itol2e)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:tol2e',mt_dbl,1,tol2e))
     1    tol2e = 1.0d-10
      endif
c
c     ----------------------
c     2-e storage
c     ----------------------
c
      intorb=.false.
      if (model2e.eq.'2eorb') intorb=.true.
      if(intorb) then
      if(.not.(intorb.and.(hftype.eq.'RHF'.or.hftype.eq.'ROHF'))) then
        call errquit('tce_init_new: 2eorb needs RHF or ROHF
     2                ref.',0,RTDB_ERR)
      endif
      endif
c
c     --------------------------
c     Determine calculation type
c     --------------------------
c
      activecalc = .false.
c
      if (reference.eq.0) then
        if (xc_gotxc()) write(LuOut,9500)
      endif
      if ((perturbative.ne.'none').and.
     1    (model.ne.'ccsd').and.(model.ne.'ccsdt'))
     2  call errquit('tce_init_new: unknown theory',0,CAPMIS_ERR)
      cperturbative = 'none'
      if (model.eq.'ccd') then
        ctype='Coupled-cluster doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccd'))
     3    call errquit('tce_init_new: ucc/uccd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'lccd') then
        ctype='Linearized coupled-cluster doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'ulccd'))
     3    call errquit('tce_init_new: ucc/ulccd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsd') then
        if (perturbative.eq.'none') then
          ctype='Coupled-cluster singles & doubles'
        else if (perturbative.eq.'(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = '(T)'
        else if (perturbative.eq.'[t]') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = '[T] or +T(CCSD)'
        else if (perturbative.eq.'cr_(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'Completely renormalized (T)'
        else if (perturbative.eq.'lr_(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'locally renormalized (T)'
        else if (perturbative.eq.'lr_(tq1)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'locally renormalized (TQ)-1'
        else if (perturbative.eq.'lr_(tq1p)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'locally renormalized (TQ)'
        else if (perturbative.eq.'cr_[t]') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'Completely renormalized [T]'
        else if (perturbative.eq.'2_t') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = '(2) triples only'
        else if (perturbative.eq.'2_tq') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = '(2) triples & quadruples (nonfactorized)'
        else if (perturbative.eq.'emb1') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'reduced completely renorm. EOMCCSD(T)-I'
        else if (perturbative.eq.'emb2') then
         ctype='Coupled-cluster singles & doubles w/ perturbation'
         cperturbative = 'reduced completely renormalized EOMCCSD(T)-II'
        else if (perturbative.eq.'creom_(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'completely renormalized EOMCCSD(T)'
        endif
        if (model.eq.'ccsd') then
         if(ccsd_var.eq.'cc2') then
          ctype='CC2 calculations'
         end if
         if(ccsd_var.eq.'lr-ccsd') then
          ctype='Locally renormalized CCSD calculations'
         end if
        end if
c
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsd').and.
     3      (module.ne.'uccsd(t)').and.
     4      (module.ne.'uccsd[t]').and.
     6      (module.ne.'cr-ccsd(t)').and.
     6      (module.ne.'lr-ccsd(t)').and.
     6      (module.ne.'lrsd(tq1p)').and.
     6      (module.ne.'lrsd(tq)-1').and.
     6      (module.ne.'creom_(t)').and.
     6      (module.ne.'emb1').and.
     6      (module.ne.'emb2').and.
     7      (module.ne.'cr-ccsd[t]'))
     8    call errquit('tce_init_new: ucc/uccsd module overridden',0,
     9    UNKNOWN_ERR)
      else if (model.eq.'lccsd') then
        ctype='Linearized coupled-cluster singles & doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'ulccsd'))
     3    call errquit('tce_init_new: ucc/ulccsd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsdta') then
        activecalc = .true.
        ctype='Coupled-cluster singles, doubles, & active triples'
        if (perturbative.eq.'2_q') then
          ctype='Coupled-cluster through triples w/ perturbation'
          cperturbative = '(2) quadruples (nonfactorized)'
        endif
         if(intorb) then
          call errquit('tce_init_new: 2eorb with active
     1                  calc.',0,RTDB_ERR)
         endif
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsdt'))
     3    call errquit('tce_init_new: ucc/uccsdt module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsdt') then
        ctype='Coupled-cluster singles, doubles, & triples'
        if (perturbative.eq.'2_q') then
          ctype='Coupled-cluster through triples w/ perturbation'
          cperturbative = '(2) quadruples (nonfactorized)'
        endif
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsdt'))
     3    call errquit('tce_init_new: ucc/uccsdt module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsdtq') then
        ctype='Coupled-cluster singles, doubles, triples, & quadruples'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsdtq'))
     3    call errquit('tce_init_new: ucc/uccsdtq module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'qcisd') then
        ctype='Quadratic configuration interaction singles & doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uqcisd'))
     3    call errquit('tce_init_new: ucc/uqcisd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cisd') then
        ctype='Configuration interaction singles & doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisd'))
     3    call errquit('tce_init_new: uci/ucisd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cisdt') then
        ctype='Configuration interaction singles, doubles, & triples'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisdt'))
     3    call errquit('tce_init_new: uci/ucisdt module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cisdtq') then
        ctype='Configuration interaction singles through quadruples'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisdtq'))
     3    call errquit('tce_init_new: uci/ucisdtq module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt2') then
        ctype='Second-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt2'))
     3    call errquit('tce_init_new: umbpt/umbpt2 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt3') then
        ctype='Third-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt3'))
     3    call errquit('tce_init_new: umbpt/umbpt3 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt4') then
        ctype='Fourth-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt4'))
     3    call errquit('tce_init_new: umbpt/umbpt4 module overridden',0,
     4    UNKNOWN_ERR)
      else
        call errquit('tce_init_new: unknown theory',0,CAPMIS_ERR)
      endif
c
c     -------------------------------------
c     Allocate global arrays for MO vectors
c     -------------------------------------
c
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nbf,'eigenvalues',
     1    l_evl(i),k_evl(i))) call errquit
     2    ('tce_init_new: failed to allocate evl',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,nbf,'occupancy',
     1    l_occ(i),k_occ(i))) call errquit
     2    ('tce_init_new: failed to allocate occ',0,MA_ERR)
        if (.not.ma_push_get(mt_int,nbf,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tce_init_new: failed to allocate irs',0,MA_ERR)
      enddo
      if (.not.movecs_read_header(movecs_out,title,basisname,
     1  scftype,nbf_verify,ipol_verify,nmo,2)) call errquit
     2  ('tce_init_new: failed to read MO file header',0,RTDB_ERR)
      if (nbf.ne.nbf_verify) call errquit
     1  ('tce_init_new: corrupted MO vectors',0,RTDB_ERR)
      if (rohf) nmo(2) = nmo(1)
      do i=1,ipol
        g_movecs(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'MO eigenvectors')
        if (rohf.and.(i.eq.2)) then
          call ga_copy(g_movecs(1),g_movecs(2))
          do j = 1,nbf
            dbl_mb(k_evl(2)+j-1) = dbl_mb(k_evl(1)+j-1)
            dbl_mb(k_occ(2)+j-1) = dbl_mb(k_occ(1)+j-1)
          enddo
        else
          if (.not.movecs_read(movecs_out,i,
     1      dbl_mb(k_occ(i)),dbl_mb(k_evl(i)),g_movecs(i)))
     2      call errquit('tce_init_new: failed to read MO vectors',0,
     3      RTDB_ERR)
        endif
        if (util_print('MO vectors',print_debug)) then
          if (nodezero) call ma_print(dbl_mb(k_evl(i)),
     1      nbf,1,'Eigenvalues')
          if (nodezero) call ma_print(dbl_mb(k_occ(i)),
     1      nbf,1,'Occupancies')
          call ga_print(g_movecs(i))
        endif
      enddo
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_init_new: no alpha irreps found',0,RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_init_new: no alpha irreps found',1,RTDB_ERR)
      endif
      do i=1,nbf
        if (oadapt) then
          int_mb(k_irs(1)+i-1) = int_mb(k_irs(1)+i-1) - 1
        else
          int_mb(k_irs(1)+i-1) = 0
        endif
      enddo
      if (ipol.eq.2) then
        if (reference.eq.0) then
          if (.not.rtdb_get(rtdb,'dft:beta irreps',mt_int,
     1      nbf,int_mb(k_irs(2))))
     2      call errquit('tce_init_new: no beta irreps
     3                    found',0,RTDB_ERR)
        else if (reference.eq.1) then
          if (rohf) then
            do j = 1,nbf
              int_mb(k_irs(2)+j-1) = int_mb(k_irs(1)+j-1) + 1
            enddo
          else
            if (.not.rtdb_get(rtdb,'scf:beta irreps',mt_int,
     1        nbf,int_mb(k_irs(2))))
     2        call errquit('tce_init_new: no beta irreps
     3                      found',1,RTDB_ERR)
          endif
        endif
        do i=1,nbf
          if (oadapt) then
            int_mb(k_irs(2)+i-1) = int_mb(k_irs(2)+i-1) - 1
          else
            int_mb(k_irs(2)+i-1) = 0
          endif
        enddo
      endif
c
c     ---------------
c     Active orbitals
c     ---------------
c
      if(activecalc) THEN
      if (.not.rtdb_get(rtdb,'tce:active_oa',mt_int,1,oactive(1)))
     1  call errquit('tce_init_new: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1  call errquit('tce_init_new: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1  call errquit('tce_init_new: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1  call errquit('tce_init_new: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1  call errquit('tce_init_new: failed to read numact',0,RTDB_ERR)
      if (restricted) then
        oactive(2) = oactive(1)
        vactive(2) = vactive(1)
      endif
      do i = 1, ipol
        if (oactive(i).gt.nocc(i)-nfc(i)) call errquit
     1  ('tce_init_new: reduce the active space',0,RTDB_ERR)
        if (vactive(i).gt.nmo(i)-nocc(i)-nfv(i)) call errquit
     1  ('tce_init_new: reduce the active space',0,RTDB_ERR)
      enddo
        if (numact.le.0.or.numact.gt.3) call errquit
     1  ('tce_init_new: change active excit. level: t3a_lvl',0,RTDB_ERR)
      do i = 1, ipol
        if (oactive(i).le.0) call errquit
     1  ('tce_init_new: incorrect size of active space',0,RTDB_ERR)
        if (vactive(i).le.0) call errquit
     1  ('tce_init_new: incorrect size of active space',0,RTDB_ERR)
      enddo
      end if
      if(.not.activecalc) THEN
       oactive(1) = 0
       oactive(2) = 0
       vactive(1) = 0
       vactive(2) = 0
      end if
c
c     ----------------------
c     Initialize integrators
c     ----------------------
c
      call int_init(rtdb,1,ao_bas_han)
      call schwarz_init(geom,ao_bas_han)
c
c     ----------------------
c     Initialize file system
c     ----------------------
c
      nfiles = 0
      if (ioalg.eq.5) then
        if (dra_init(max_files,1.0d11,1.0d11,
     1    dble(buffer_size)*dble(bytes)).ne.0)
     2    call errquit('tce_init_new: dra initialization failed',0,
     3    DISK_ERR)
      endif
c
c     -------------
c     Printing info
c     -------------
c
      if (nodezero.and.util_print('information',print_low)) then
        call util_print_centered
     1  (LuOut,'NWChem Extensible Many-Electron Theory Module',
     2   40,.true.)
        write(LuOut,*)
        call util_print_centered
     1  (LuOut,'======================================================',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'This portion of the program was automatically',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'generated by a Tensor Contraction Engine (TCE).',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'The development of this portion of the program',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'and TCE was supported by US Department of Energy,',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'Office of Science, Office of Basic Energy Science.',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'TCE is a product of Battelle and PNNL.',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'Please cite: S.Hirata, J.Phys.Chem.A 107, 9887 (2003).',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'======================================================',
     2  40,.false.)
        write(LuOut,*)
        if (title.ne.' ') then
          call util_print_centered(LuOut,title,40,.false.)
          write(LuOut,*)
          write(LuOut,*)
        endif
      endif
      if (nodezero.and.util_print('general information',
     1  print_default)) then
        call util_print_centered
     1    (LuOut,'General Information',20,.true.)
        write(LuOut,9200) ga_nnodes()
        if (ipol.eq.1) then
          if (reference.eq.0) then
            write(LuOut,9010) 'Restricted Kohn-Sham'
          else if (reference.eq.1) then
            write(LuOut,9010) 'Restricted Hartree-Fock'
          endif
          write(LuOut,9030) nocc(1)*2,nocc(1),nocc(1)
        else if (ipol.eq.2) then
          if (reference.eq.0) then
            write(LuOut,9010) 'Unrestricted Kohn-Sham'
          else if (reference.eq.1) then
            if (rohf) then
              write(LuOut,9010) 'Restricted open-shell Hartree-Fock'
            else
              write(LuOut,9010) 'Unrestricted Hartree-Fock'
            endif
          endif
          write(LuOut,9030) nocc(1)+nocc(2),nocc(1),nocc(2)
        endif
        if (ipol.eq.1) then
          write(LuOut,9410) nmo(1)*2,nmo(1),nmo(1)
        else if (ipol.eq.2) then
          write(LuOut,9410) nmo(1)+nmo(2),nmo(1),nmo(2)
        endif
        if (ipol.eq.1) then
          write(LuOut,9430) nfc(1),nfc(1)
        else if (ipol.eq.2) then
          write(LuOut,9430) nfc(1),nfc(2)
        endif
        if (ipol.eq.1) then
          write(LuOut,9440) nfv(1),nfv(1)
        else if (ipol.eq.2) then
          write(LuOut,9440) nfv(1),nfv(2)
        endif
        if(activecalc) then
         write(LuOut,9540) oactive(1),oactive(2)
         write(LuOut,9550) vactive(1),vactive(2)
         write(LuOut,9555) numact
        end if
        write(LuOut,9420) multchar(mult)
        write(LuOut,9050) nbf
        write(LuOut,9060) nshells
        if (oskel) then
          onoff1='on '
        else
          onoff1='off'
        endif
        if (oadapt) then
          onoff2='on '
        else
          onoff2='off'
        endif
        write(LuOut,9040) onoff1,onoff2
        write(LuOut,9100) tol2e
        if (ipol.eq.1) then
          if (nmo(1).lt.nbf) write(LuOut,9400)
        else
          if ((nmo(1).lt.nbf).or.(nmo(2).lt.nbf))
     1      write(LuOut,9400)
        endif
      endif
      if (nodezero.and.util_print('correlation information',
     1  print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'Correlation Information',20,.true.)
        write(LuOut,9000) ctype
        write(LuOut,9230) cperturbative
        write(LuOut,9130) maxiter
        write(LuOut,9110) thresh
        write(LuOut,9115) zlshift
        if (diis .gt. 0) then
          write(LuOut,9150) diis,'-th order DIIS'
        else
          write(LuOut,9160) 'Jacobi'
        endif
        if (nroots .gt. 0) then
          write(LuOut,9170) nroots
          write(LuOut,9180) target
          write(LuOut,9190) targetsym
          if (symmetry) then
            write(LuOut,9210) 'on '
          else
            write(LuOut,9210) 'off'
          endif
          if (left) then
            write(LuOut,9220) 'on '
          else
            write(LuOut,9220) 'off'
          endif
        endif
        if (ioalg.eq.0) then
          write(LuOut,9140) 'Fortran77 Direct Access'
        else if (ioalg.eq.1) then
          write(LuOut,9140) 'Exclusive Access File Library'
        else if (ioalg.eq.2) then
          write(LuOut,9140) 'Global Array Library'
        else if (ioalg.eq.3) then
          write(LuOut,9140) 'Shared File Library'
        else if (ioalg.eq.4) then
          write(LuOut,9140) 'Replicated Exclusive Access Files'
        else if (ioalg.eq.5) then
          write(LuOut,9140) 'Disk Resident Array Algorithm'
        else if (ioalg.eq.6) then
          write(LuOut,9140) 'GA/EAF Replicated Algorithm'
        else
          call errquit('tce_init_new: unknown I/O
     1                  scheme',ioalg,CAPMIS_ERR)
        endif
      endif
      if (nodezero) call util_flush(LuOut)
 9200 format(3x,'   Number of processors :',2x,i4)
 9000 format(3x,'       Calculation type : ',a)
 9230 format(3x,'Perturbative correction : ',a)
 9010 format(3x,'      Wavefunction type : ',a)
 9030 format(3x,'       No. of electrons :',2x,i4,/,
     1       3x,'        Alpha electrons :',2x,i4,/,
     2       3x,'         Beta electrons :',2x,i4)
 9410 format(3x,'        No. of orbitals :',2x,i4,/,
     1       3x,'         Alpha orbitals :',2x,i4,/,
     2       3x,'          Beta orbitals :',2x,i4)
 9430 format(3x,'     Alpha frozen cores :',2x,i4,/,
     1       3x,'      Beta frozen cores :',2x,i4)
 9440 format(3x,'  Alpha frozen virtuals :',2x,i4,/,
     1       3x,'   Beta frozen virtuals :',2x,i4)
 9420 format(3x,'      Spin multiplicity : ',a)
 9040 format(3x,'     Use of symmetry is : ',a3,/,
     1       3x,'   Symmetry adaption is : ',a3)
 9050 format(3x,' Number of AO functions :',2x,i4)
 9060 format(3x,'    Number of AO shells :',2x,i4)
 9100 format(3x,'      Schwarz screening : ',d8.2)
 9110 format(3x,'     Residual threshold : ',d8.2)
 9115 format(3x,'       DIIS level shift : ',d8.2)
 9130 format(3x,'         Max iterations :',2x,i7)
 9140 format(3x,'             I/O scheme : ',a)
 9150 format(3x,'       Amplitude update : ',i2,a)
 9160 format(3x,'       Amplitude update : ',a)
 9170 format(3x,'  No. of excited states :',2x,i4)
 9180 format(3x,'            Target root :',2x,i4)
 9190 format(3x,'        Target symmetry :',1x,a4)
 9210 format(3x,'   Symmetry restriction :',1x,a3)
 9220 format(3x,'Dipole & oscillator str :',1x,a3)
 9400 format(/,2x,
     1  '!! WARNING !! The number of MO is less than the number of AO')
 9500 format(/,2x,
     1  '!! WARNING !! A non-HF reference wave function')
 9540 format(3x,' Alpha active occupieds :',2x,i4,/,
     1       3x,'  Beta active occupieds :',2x,i4)
 9550 format(3x,'  Alpha active virtuals :',2x,i4,/,
     1       3x,'   Beta active virtuals :',2x,i4)
 9555 format(3x,'  T3 active excitation level :',2x,i4)
c
c     ------
c     Return
c     ------
c
      return
      end
