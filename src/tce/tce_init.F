      subroutine tce_init(rtdb)
c
c $Id: tce_init.F,v 1.71 2009/01/28 18:11:42 marat Exp $
c
c     Initialize for many-electron calculations.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
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
#include "tce_restart.fh"
#include "dra.fh"
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
      double precision t_i,dt
      double precision fac,xxx
      double precision pi
      parameter(pi = 3.14159265358979323846264338327950288419D0)
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
     1  call errquit('tce_init: geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('tce_init: no geometry found',0,GEOM_ERR)
      if (.not.rtdb_get(rtdb,'dft:skeleton',mt_log,1,oskel)) then
        oskel=sym_number_ops(geom).gt.0
      endif
      if (.not.rtdb_get(rtdb,'dft:adapt',mt_log,1,oadapt)) then
        oadapt=sym_number_ops(geom).gt.0
      endif
      if (.not.sym_abelian_group(geom))
     1  call errquit('tce_init: non-Abelian symmetry not permitted',0,
     2  GEOM_ERR)
c
c     ---------
c     Basis set
c     ---------
c
      if (.not.bas_create(ao_bas_han,'ao basis'))
     1  call errquit('tce_init: bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,ao_bas_han,'ao basis'))
     1  call errquit('tce_init: no ao basis found',0,BASIS_ERR)
      if (.not.int_normalize(rtdb,ao_bas_han))
     1  call errquit('tce_init: int_normalize failed',0,INT_ERR)
      if (.not.bas_numbf(ao_bas_han,nbf))
     1  call errquit('tce_init: failed to read nbf',0,BASIS_ERR)
      if (.not.bas_numcont(ao_bas_han,nshells))
     1  call errquit('tce_init: failed to read nshells',0,BASIS_ERR)
c
c     ----------
c     MO vectors
c     ----------
c
      if (reference.eq.0) then
        if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_out))
     1    call errquit('tce_init: no final MO found',0,RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     1    call errquit('tce_init: no final MO found',0,RTDB_ERR)
      endif
c
c     ------------------------------------------------------
c     Number of electrons and orbitals and spin multiplicity
c     ------------------------------------------------------
c
      if (reference.eq.0) then
        rohf = .false.
        if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     1    call errquit('tce_init: failed to read ipol',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc))
     1    call errquit('tce_init: failed to read nocc',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'dft:mult',mt_int,1,mult))
     1    call errquit('tce_init: failed to read mult',0,RTDB_ERR)
        if (ipol.eq.1) then
          restricted = .true.
        else
          restricted = .false.
        endif
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nclosed))
     1    call errquit('tce_init: failed to read nclosed',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'scf:nopen',mt_int,1,nopen))
     1    call errquit('tce_init: failed to read nopen',0,RTDB_ERR)
        nocc(1) = nclosed + nopen
        nocc(2) = nclosed
        mult = nopen + 1
        if(.not.rtdb_cget(rtdb,'scf:scftype',1,hftype))
     1    call errquit('tce_init: failed to read scftype',0,RTDB_ERR)
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
         call errquit('tce_init: failed to get number of cores',0,
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
     1  call errquit('tce_init: model2e not  specified',0,RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'tce:model',1,model))
     1  call errquit('tce_init: no theory specified',0,RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'tce:perturbative',1,perturbative)) then
        perturbative = 'none'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1  call errquit('tce_init: failed to write perturbative',0,
     2  RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'tce:ccsdvar',1,ccsd_var)) then
        ccsd_var = 'none'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1  call errquit('tce_init: failed to write ccsdvar',0,
     2  RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:thresh',mt_dbl,1,thresh)) then
        thresh = 1.0d-6
        if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_init: failed to write thresh',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:printtol',mt_dbl,1,printtol)) then
        printtol = 1.0d-1
        if (.not.rtdb_put(rtdb,'tce:printtol',mt_dbl,1,printtol))
     1  call errquit('tce_init: failed to write printtol',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:zlshift',mt_dbl,1,zlshift)) then
        zlshift = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:zlshift',mt_dbl,1,zlshift))
     1  call errquit('tce_init: failed to write zlshift',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:zlshiftl',mt_dbl,1,zlshift)) then
        zlshiftl = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:zlshiftl',mt_dbl,1,zlshiftl))
     1  call errquit('tce_init: failed to write zlshiftl',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:zlshift2',mt_dbl,2,zlshift2)) then
        zlshift2(1) = 0.0d0
        zlshift2(2) = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:zlshift2',mt_dbl,2,zlshift2))
     1  call errquit('tce_init: failed to write zlshift2',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:zlshift3',mt_dbl,2,zlshift3)) then
        zlshift3(1) = 0.0d0 ! denom_power = 2 aka real
        zlshift3(2) = 0.0d0 ! denom_power = 1 aka imag
        if (.not.rtdb_put(rtdb,'tce:zlshift3',mt_dbl,2,zlshift3))
     1  call errquit('tce_init: failed to write zlshift3',0,RTDB_ERR)
      endif
c#ifdef DEBUG_PRINT
c      if(nodezero) write(6,*) 'zlshift  = ',zlshift
c      if(nodezero) write(6,*) 'zlshift2 = ',zlshift2
c      if(nodezero) write(6,*) 'zlshift3 = ',zlshift3
c#endif
      if (.not.rtdb_get(rtdb,'tce:damp',mt_dbl,1,damp)) then
        damp = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:damp',mt_dbl,1,damp))
     1  call errquit('tce_init: failed to write damp',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:damp2',mt_dbl,1,damp2)) then
        damp2 = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:damp2',mt_dbl,1,damp2))
     1  call errquit('tce_init: failed to write damp2',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:damp3',mt_dbl,1,damp3)) then
        damp3 = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:damp3',mt_dbl,1,damp3))
     1  call errquit('tce_init: failed to write damp3',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:maxiter',mt_int,1,maxiter)) then
        maxiter = 100
        if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_init: failed to write maxiter',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ioalg',mt_int,1,ioalg)) then
        ioalg = 2
        if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_init: failed to write ioalg',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:eoms',mt_int,1,eom_solver)) then
        eom_solver = 1
        if (.not.rtdb_put(rtdb,'tce:eoms',mt_int,1,eom_solver))
     1  call errquit('tce_init: failed to write eoms',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diis',mt_int,1,diis)) then
        diis = 5
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init: failed to write diis',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diis2',mt_int,1,diis2)) then
        diis2 = 5
        if (.not.rtdb_put(rtdb,'tce:diis2',mt_int,1,diis2))
     1  call errquit('tce_init: failed to write diis2',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diis3',mt_int,1,diis3)) then
        diis3 = 5
        if (.not.rtdb_put(rtdb,'tce:diis3',mt_int,1,diis3))
     1  call errquit('tce_init: failed to write diis3',0,RTDB_ERR)
      endif
      if (diis .gt. maxdiis) then
        diis = maxdiis
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init: failed to write diis',0,RTDB_ERR)
      endif
      if (diis2 .gt. maxdiis) then
        diis2 = maxdiis
        if (.not.rtdb_put(rtdb,'tce:diis2',mt_int,1,diis2))
     1  call errquit('tce_init: failed to write diis2',0,RTDB_ERR)
      endif
      if (diis3 .gt. maxdiis) then
        diis3 = maxdiis
        if (.not.rtdb_put(rtdb,'tce:diis3',mt_int,1,diis3))
     1  call errquit('tce_init: failed to write diis3',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:cc_ir_alg',mt_int,1,cc_ir_alg)) then
        cc_ir_alg = 1
        if (.not.rtdb_put(rtdb,'tce:cc_ir_alg',mt_int,1,cc_ir_alg))
     2  call errquit('tce_init: failed to write cc_ir_alg',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:guess_lr_real',mt_log,1,
     1                  guess_lr_real)) then
        guess_lr_real = .false.
        if (.not.rtdb_put(rtdb,'tce:guess_lr_real',mt_log,1,
     1                    guess_lr_real))
     2  call errquit('tce_init: failed to write guess_lr_real',0,
     3               RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:guess_qr_real',mt_log,1,
     1                  guess_qr_real)) then
        guess_qr_real = .false.
        if (.not.rtdb_put(rtdb,'tce:guess_qr_real',mt_log,1,
     1                    guess_qr_real))
     2  call errquit('tce_init: failed to write guess_qr_real',0,
     3               RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:guess_ir_real',mt_log,1,
     1                  guess_ir_real)) then
        guess_ir_real = .false.
        if (.not.rtdb_put(rtdb,'tce:guess_ir_real',mt_log,1,
     1                    guess_ir_real))
     2  call errquit('tce_init: failed to write guess_ir_real',0,
     3               RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:guess_ir_imag',mt_log,1,
     1                  guess_ir_imag)) then
        guess_ir_imag = .false.
        if (.not.rtdb_put(rtdb,'tce:guess_ir_imag',mt_log,1,
     1                    guess_ir_imag))
     2  call errquit('tce_init: failed to write guess_ir_imag',0,
     3               RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:nroots',mt_int,1,nroots)) then
        nroots = 0
        if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1  call errquit('tce_init: failed to write nroots',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:maxdiff',mt_dbl,1,maxdiff)) then
        maxdiff = 0.5d0
        if (.not.rtdb_put(rtdb,'tce:maxdiff',mt_dbl,1,maxdiff))
     1  call errquit('tce_init: failed to write maxdiff',0,RTDB_ERR)
c
      endif
      if (.not.rtdb_get(rtdb,'tce:hbard',mt_int,1,hbard)) then
        hbard = 500
        if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1  call errquit('tce_init: failed to write hbard',0,RTDB_ERR)
      endif
c 4index trans.
      if (.not.rtdb_get(rtdb,'tce:maxs',mt_int,1,imaxsize)) then
        imaxsize = 30
        if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,imaxsize))
     1  call errquit('tce_init: failed to write maxs',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ichopx',mt_int,1,ichop)) then
        ichop = 1
        if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichop))
     1  call errquit('tce_init: failed to write ichop',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:i4im',mt_int,1,fast2e)) then
        fast2e = 1
        if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,fast2e))
     1  call errquit('tce_init: failed to write i4im',0,RTDB_ERR)
      endif
c      if ((model2e.ne.'2eorb').and.(fast2e.gt.3)) then
c        fast2e = 1
c        if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,fast2e))
c     1  call errquit('tce_init: failed to write i4im',0,RTDB_ERR)
c      endif
      if (.not.rtdb_get(rtdb,'tce:idiskx',mt_int,1,idisk)) then
c        idisk = 0
        if ( (fast2e.eq.5).or.(fast2e.eq.6) ) then
          idisk = 1
        else
          idisk = 0
        endif
        if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idisk))
     1  call errquit('tce_init: failed to write idiskx',0,RTDB_ERR)
      endif
c CC2 and LR-CCSD
      if (.not.rtdb_cget(rtdb,'tce:ccsdvar',1,ccsd_var)) then
        ccsd_var='none'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1  call errquit('tce_init: failed to write ccsdvar',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:target',mt_int,1,target)) then
         target=1
        if (.not.rtdb_put(rtdb,'tce:target',mt_int,1,target))
     1  call errquit('tce_init: failed to write target',0,RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'tce:targetsym',1,targetsym)) then
        targetsym='none'
        if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1  call errquit('tce_init: failed to write targetsym',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:symmetry',mt_log,1,symmetry)) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1  call errquit('tce_init: failed to write symmetry',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:left',mt_log,1,left)) then
        left=.false.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_init: failed to write left',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:multipole',mt_int,1,multipole)) then
        multipole=0
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1  call errquit('tce_init: failed to write multipole',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:diagnt1',mt_log,1,diagnt1)) then
        diagnt1=.false.
        if (.not.rtdb_put(rtdb,'tce:diagnt1',mt_log,1,diagnt1))
     1  call errquit('tce_init: failed to write diagnt1',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ampnorms',mt_log,1,ampnorms)) then
        ampnorms=.false.
        if (.not.rtdb_put(rtdb,'tce:ampnorms',mt_log,1,ampnorms))
     1  call errquit('tce_init: failed to write ampnorms',0,RTDB_ERR)
      endif
c
c    Beginning of response property input section
c
      if (.not.rtdb_get(rtdb,'tce:lineresp',mt_log,1,lineresp)) then
        lineresp=.false.
        if (.not.rtdb_put(rtdb,'tce:lineresp',mt_log,1,lineresp))
     1  call errquit('tce_init: failed to write lineresp',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:imagresp',mt_log,1,imagresp)) then
        imagresp=.false.
        if (.not.rtdb_put(rtdb,'tce:imagresp',mt_log,1,imagresp))
     1  call errquit('tce_init: failed to write imagresp',0,RTDB_ERR)
      endif
      if (lineresp.or.imagresp) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_init: failed to write left',0,RTDB_ERR)
      endif
      if (imagresp) then
        lineresp=.false.
        if (.not.rtdb_put(rtdb,'tce:lineresp',mt_log,1,lineresp))
     1  call errquit('tce_init: failed to write lineresp',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:leftresp',mt_log,1,leftresp)) then
        leftresp=.false.
        if (.not.rtdb_put(rtdb,'tce:leftresp',mt_log,1,leftresp))
     1  call errquit('tce_init: failed to write leftresp',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:quadresp',mt_log,1,quadresp)) then
        quadresp=.false.
        if (.not.rtdb_put(rtdb,'tce:quadresp',mt_log,1,quadresp))
     1  call errquit('tce_init: failed to write quadresp',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get_info(rtdb,'tce:afreq',ma_type,anumfreq,
     1  date)) then
        anumfreq = 1
        if (.not.rtdb_put(rtdb,'tce:anumfreq',mt_int,1,anumfreq))
     1   call errquit('tce_init: failed to write anumfreq',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:afreq',mt_dbl,anumfreq,
     1  afreq)) then
        afreq(1) = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:afreq',mt_dbl,anumfreq,afreq))
     1    call errquit('tce_init: failed to write afreq',0,RTDB_ERR)
      endif
      if (.not.rtdb_get_info(rtdb,'tce:bfreq',ma_type,bnumfreq,
     1  date)) then
        bnumfreq = 1
        if (.not.rtdb_put(rtdb,'tce:bnumfreq',mt_int,1,bnumfreq))
     1   call errquit('tce_init: failed to write bnumfreq',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:bfreq',mt_dbl,bnumfreq,
     1  bfreq)) then
        bfreq(1) = 0.0d0
        if (.not.rtdb_put(rtdb,'tce:bfreq',mt_dbl,bnumfreq,bfreq))
     1    call errquit('tce_init: failed to write bfreq',0,RTDB_ERR)
      endif
c
      if(imagresp) then !--------------------------
c
      if (.not.rtdb_get(rtdb,'tce:ifreqauto',mt_log,1,ifreqauto)) then
        ifreqauto=.true.
        if (.not.rtdb_put(rtdb,'tce:ifreqauto',mt_log,1,ifreqauto))
     1  call errquit('tce_init: failed to write ifreqauto',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ifreqgrgl',mt_log,1,ifreqgrgl)) then
        ifreqgrgl=.false.
        if (.not.rtdb_put(rtdb,'tce:ifreqgrgl',mt_log,1,ifreqgrgl))
     1  call errquit('tce_init: failed to write ifreqgrgl',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ifreqgrch',mt_log,1,ifreqgrch)) then
        ifreqgrch=.false.
        if (.not.rtdb_put(rtdb,'tce:ifreqgrch',mt_log,1,ifreqgrch))
     1  call errquit('tce_init: failed to write ifreqgrch',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:ifreqgrst',mt_log,1,ifreqgrst)) then
        ifreqgrst=.false.
        if (.not.rtdb_put(rtdb,'tce:ifreqgrch',mt_log,1,ifreqgrst))
     1  call errquit('tce_init: failed to write ifreqgrst',0,RTDB_ERR)
      endif
#ifdef DEBUG_PRINT
       if(nodezero) write(LuOut,*) 'ifreqauto = ',ifreqauto 
       if(nodezero) write(LuOut,*) 'ifreqgrgl = ',ifreqgrgl 
       if(nodezero) write(LuOut,*) 'ifreqgrch = ',ifreqgrch
       if(nodezero) write(LuOut,*) 'ifreqgrst = ',ifreqgrst
#endif
c
c     Check for valid options
c
      if (( ifreqgrgl.and.ifreqgrch ).or.
     1    ( ifreqgrgl.and.ifreqgrst ).or.
     1    ( ifreqgrch.and.ifreqgrst )) then
        if(nodezero) then
           write(LuOut,*) 'These are mutually-exclusive:'
           write(LuOut,*) 'Gauss-Legendre  (ifreqgrgl = T)'
           write(LuOut,*) 'Gauss-Chebyshev (ifreqgrgc = T)'
           write(LuOut,*) 'Stanton         (ifreqgrst = T)'
           call util_flush(LuOut)
        endif
        call errquit('tce_init: user input error',0,INPUT_ERR)
      endif
      if (ifreqauto) then
c
c     Gauss-Chebychev grid identical to Dalton 2.0
c
      if (ifreqgrch) then
#ifdef DEBUG_PRINT
        if(nodezero) write(LuOut,*) 'Gauss-Chebychev grid'
        if(nodezero) call util_flush(LuOut)
#endif
        if (.not.rtdb_get(rtdb,'tce:inumfreq',mt_int,1,inumfreq)) then
          inumfreq = 10
          if (.not.rtdb_put(rtdb,'tce:inumfreq',mt_int,1,inumfreq))
     1     call errquit('tce_init: failed to write inumfreq',0,RTDB_ERR)
        endif
        fac = pi/(4.0d0*inumfreq)
#ifdef DEBUG_PRINT
        if(nodezero) write(LuOut,*) 'fac = ',fac
        if(nodezero) call util_flush(LuOut)
#endif
        ifreq(1) = 0.0d0
        do i = 1, inumfreq-1
          j = (2*inumfreq-2) - 2*(i-1)
          xxx = fac*j
          ifreq(i+1) = (dtan(xxx))**(-2)
#ifdef DEBUG_PRINT
          if(nodezero) write(LuOut,36) i,j,xxx,ifreq(i)
   36 format(1x,'$$$',2x,i10,2x,i10,2x,f14.6,2x,f14.6)
#endif
        enddo 
#ifdef DEBUG_PRINT
        if(nodezero) call util_flush(LuOut)
#endif
c
c     Automatic grid using Stanton's scheme [Phys. Rev. A 49, 1698 (1994)]
c
      elseif (ifreqgrst) then
#ifdef DEBUG_PRINT
       if(nodezero) write(LuOut,*) 'Stanton grid'
#endif
        if (.not.rtdb_get(rtdb,'tce:ifreqw0',mt_dbl,1,ifreqw0)) then
          ifreqw0 = 1.0d0
          if (.not.rtdb_put(rtdb,'tce:ifreqw0',mt_dbl,1,ifreqw0))
     1      call errquit('tce_init: failed to write ifreqw0',0,RTDB_ERR)
        endif
        if (.not.rtdb_get(rtdb,'tce:inumfreq',mt_int,1,inumfreq)) then
          inumfreq = 10
          if (.not.rtdb_put(rtdb,'tce:inumfreq',mt_int,1,inumfreq))
     1     call errquit('tce_init: failed to write inumfreq',0,RTDB_ERR)
        endif
        dt = 1.0d0/(inumfreq-1)
        do i = inumfreq, 1, -1
          t_i = (i-1)*dt
          ifreq(i) = ifreqw0 * (1.0d0 - t_i)/(1.0d0 + t_i)
#ifdef DEBUG_PRINT
          if(nodezero) write(LuOut,37) i,t_i,ifreq(i)
   37 format(1x,'$$$',2x,i10,2x,f14.6,2x,f14.6)
#endif
        enddo
#ifdef DEBUG_PRINT
        if(nodezero) call util_flush(LuOut)
#endif
        if (.not.rtdb_put(rtdb,'tce:ifreq',mt_dbl,inumfreq,ifreq))
     1    call errquit('tce_init: failed to write ifreq',0,RTDB_ERR)
      else ! if (ifreqgrgl) then
c
c     Gauss-Legendre grid with same coefficients as Dalton 2.0
c
#ifdef DEBUG_PRINT
        if(nodezero) write(LuOut,*) 'Gauss-Legendre grid'
        if(nodezero) call util_flush(LuOut)
#endif
        inumfreq = 17
        if (.not.rtdb_put(rtdb,'tce:inumfreq',mt_int,1,inumfreq))
     1     call errquit('tce_init: failed to write inumfreq',0,RTDB_ERR)
        ifreq(1)  = 0.0000000000d0
        ifreq(2)  = 0.0000011354d0
        ifreq(3)  = 0.0000324954d0
        ifreq(4)  = 0.0002074939d0
        ifreq(5)  = 0.0007766098d0
        ifreq(6)  = 0.0022314002d0
        ifreq(7)  = 0.0055272185d0
        ifreq(8)  = 0.0125684274d0
        ifreq(9)  = 0.0273216537d0
        ifreq(10) = 0.0585616091d0
        ifreq(11) = 0.1273031181d0
        ifreq(12) = 0.2894765239d0
        ifreq(13) = 0.7170385627d0
        ifreq(14) = 2.0602366551d0
        ifreq(15) = 7.7110713698d0
        ifreq(16) = 49.2377677794d0
        ifreq(17) = 1409.1898779477d0
        if (.not.rtdb_put(rtdb,'tce:ifreq',mt_dbl,inumfreq,ifreq))
     1    call errquit('tce_init: failed to write ifreq',0,RTDB_ERR)
      endif
c
c     User-specified grid
c
      else
#ifdef DEBUG_PRINT
        if(nodezero) write(LuOut,*) 'User-specified grid'
        if(nodezero) call util_flush(LuOut)
#endif
        if (.not.rtdb_get_info(rtdb,'tce:ifreq',ma_type,inumfreq,
     1    date)) then
          inumfreq = 1
          if (.not.rtdb_put(rtdb,'tce:inumfreq',mt_int,1,inumfreq))
     1     call errquit('tce_init: failed to write inumfreq',0,RTDB_ERR)
        endif
        if (.not.rtdb_get(rtdb,'tce:ifreq',mt_dbl,inumfreq,ifreq)) then
          ifreq(1) = 0.0d0
          if (.not.rtdb_put(rtdb,'tce:ifreq',mt_dbl,inumfreq,ifreq))
     1      call errquit('tce_init: failed to write ifreq',0,RTDB_ERR)
        endif
      endif
      endif !------------------------------
c
      if (anumfreq.gt.maxnumfreq) then
        anumfreq = maxnumfreq
        if (.not.rtdb_put(rtdb,'tce:anumfreq',mt_int,1,anumfreq))
     1    call errquit('tce_init: failed to write anumfreq',0,RTDB_ERR)
      endif
      if (bnumfreq.gt.maxnumfreq) then
        bnumfreq = maxnumfreq
        if (.not.rtdb_put(rtdb,'tce:bnumfreq',mt_int,1,bnumfreq))
     1    call errquit('tce_init: failed to write bnumfreq',0,RTDB_ERR)
      endif
      if (inumfreq.gt.maxnumfreq) then
        inumfreq = maxnumfreq
        if (.not.rtdb_put(rtdb,'tce:inumfreq',mt_int,1,inumfreq))
     1    call errquit('tce_init: failed to write inumfreq',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'tce:respaxis',mt_log,3,respaxis)) then
        respaxis(1) = .true.
        respaxis(2) = .true.
        respaxis(3) = .true.
        if (.not.rtdb_put(rtdb,'tce:respaxis',mt_log,3,respaxis))
     1    call errquit('tce_init: failed to write respaxis',0,RTDB_ERR)
      endif
c
c     these force skipping of alpha, beta and gamma evaluation
c     which is useful during restarts
c
      if (.not.rtdb_get(rtdb,'tce:nopolar',mt_log,1,nopolar)) then
        nopolar=.false.
        if (.not.rtdb_put(rtdb,'tce:nopolar',mt_log,1,nopolar))
     1  call errquit('tce_init: failed to write nopolar',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:nohyper',mt_log,1,nohyper)) then
        nohyper=.false.
        if (.not.rtdb_put(rtdb,'tce:nohyper',mt_log,1,nohyper))
     1  call errquit('tce_init: failed to write nohyper',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:nogamma',mt_log,1,nogamma)) then
        nogamma=.false.
        if (.not.rtdb_put(rtdb,'tce:nogamma',mt_log,1,nogamma))
     1  call errquit('tce_init: failed to write nogamma',0,RTDB_ERR)
      endif
c
c     End of response property input section
c
      if ((symmetry).and.(targetsym.eq.'none')) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1    call errquit('tce_init: failed to write symmetry',0,RTDB_ERR)
      endif
c
c     vectorized CCSD(T)
c
      if (.not.rtdb_get(rtdb,'tce:fast_ccsd_t',mt_log,1,
     1  fast_ccsd_t)) then
        fast_ccsd_t = .false.
        if (.not.rtdb_put(rtdb,'tce:fast_ccsd_t',mt_log,1,fast_ccsd_t))
     1    call errquit('tce_init: failed for fast_ccsd_t',0,RTDB_ERR)
      endif
c
c     Temporary restart input
c
      if (.not.rtdb_get(rtdb,'tce:restart_ccsd_t',mt_int,1,
     1  restart_ccsd_t)) then
        restart_ccsd_t = 0
        if (.not.rtdb_put(rtdb,'tce:restart_ccsd_t',mt_int,1,
     1                    restart_ccsd_t))
     2  call errquit('tce_init: failed to write restart_ccsd_t',0,
     3               RTDB_ERR)
      endif
c
c     Save toggles
c
      if (.not.rtdb_get(rtdb,'tce:save_interval',mt_int,1,
     1                            save_interval)) then
        save_interval = 1000
        if (.not.rtdb_put(rtdb,'tce:save_interval',mt_int,1,
     1      save_interval)) call errquit(
     2     'tce_init: failed to write save_interval',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:restart_ignore_failure',
     &                  mt_log,4,restart_ignore_failure)) then
        restart_ignore_failure(1)=.false. ! all integrals
        restart_ignore_failure(2)=.false. ! all zeroth-order amplitudes (including EOM)
        restart_ignore_failure(3)=.false. ! all response amplitudes
        if (.not.rtdb_put(rtdb,'tce:restart_ignore_failure',
     1                    mt_log,4,restart_ignore_failure))
     1  call errquit('tce_init: failed to write restart_ignore_failure',
     3               0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_integrals',mt_log,5,
     1         save_integrals)) then
        save_integrals(1)=.false. ! 1-e integrals
        save_integrals(2)=.false. ! 2-e integrals
        save_integrals(3)=.false. ! dipole integrals
        save_integrals(4)=.false. ! quadrupole integrals
        save_integrals(5)=.false. ! octopole integrals
        if (.not.rtdb_put(rtdb,'tce:save_integrals',mt_log,5,
     1           save_integrals))
     1  call errquit('tce_init: failed to write save_integrals',
     2               0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_t',mt_log,4,save_t)) then
        do i = 1, 4
          save_t(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:save_t',mt_log,4,save_t))
     1  call errquit('tce_init: failed to write save_t',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_l',mt_log,4,save_l)) then
        do i = 1, 4
          save_l(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:save_l',mt_log,4,save_l))
     1  call errquit('tce_init: failed to write save_l',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_tr',mt_log,12,save_tr)) then
c       Real positive and negative frequency amplitudes, then
c       the real component of imaginary response when doing CC-IR          
c       S(+),D(+),T(+),Q(+),S(-),D(-),T(-),Q(-),S(r),D(r),T(r),Q(r)
        do i = 1, 12
          save_tr(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:save_tr',mt_log,12,save_tr))
     1  call errquit('tce_init: failed to write save_tr',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_lr',mt_log,12,save_lr)) then
        do i = 1, 12
          save_lr(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:save_lr',mt_log,12,save_lr))
     1  call errquit('tce_init: failed to write save_lr',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_x',mt_log,4,save_x)) then
        save_x(1)=.false.
        save_x(2)=.false.
        save_x(3)=.false.
        save_x(4)=.false.
        if (.not.rtdb_put(rtdb,'tce:save_x',mt_log,4,save_x))
     1  call errquit('tce_init: failed to write save_x',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:save_y',mt_log,4,save_y)) then
        save_y(1)=.false.
        save_y(2)=.false.
        save_y(3)=.false.
        save_y(4)=.false.
        if (.not.rtdb_put(rtdb,'tce:save_y',mt_log,4,save_y))
     1  call errquit('tce_init: failed to write save_y',0,RTDB_ERR)
      endif
c
c     Read toggles
c
c moved to tce_energy.F so skips SCF on integral read
c
c      if (.not.rtdb_get(rtdb,'tce:read_integrals',mt_log,5,
c     1         read_integrals)) then
c        read_integrals(1)=.false. ! 1-e integrals
c        read_integrals(2)=.false. ! 2-e integrals
c        read_integrals(3)=.false. ! dipole integrals
c        read_integrals(4)=.false. ! quadrupole integrals
c        read_integrals(5)=.false. ! octopole integrals
c        if (.not.rtdb_put(rtdb,'tce:read_integrals',mt_log,5,
c     1           read_integrals))
c     1  call errquit('tce_init: failed to write read_integrals',
c     2               0,RTDB_ERR)
c      endif
      if (.not.rtdb_get(rtdb,'tce:read_t',mt_log,4,read_t)) then
        do i = 1, 4
          read_t(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:read_t',mt_log,4,read_t))
     1  call errquit('tce_init: failed to write read_t',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:read_l',mt_log,4,read_l)) then
        do i = 1, 4
          read_l(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:read_l',mt_log,4,read_l))
     1  call errquit('tce_init: failed to write read_l',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:read_tr',mt_log,12,read_tr)) then
c       Real positive and negative frequency amplitudes, then
c       the real component of imaginary response when doing CC-IR          
c       S(+),D(+),T(+),Q(+),S(-),D(-),T(-),Q(-),S(r),D(r),T(r),Q(r)
        do i = 1, 12
          read_tr(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:read_tr',mt_log,12,read_tr))
     1  call errquit('tce_init: failed to write read_tr',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:read_lr',mt_log,12,read_lr)) then
c       Real positive and negative frequency amplitudes, then
c       the real component of imaginary response when doing CC-IR          
c       S(+),D(+),T(+),Q(+),S(-),D(-),T(-),Q(-),S(r),D(r),T(r),Q(r)
        do i = 1, 12
          read_lr(i)=.false.
        enddo
        if (.not.rtdb_put(rtdb,'tce:read_lr',mt_log,12,read_lr))
     1  call errquit('tce_init: failed to write read_lr',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:read_x',mt_log,4,read_x)) then
        read_x(1)=.false.
        read_x(2)=.false.
        read_x(3)=.false.
        read_x(4)=.false.
        if (.not.rtdb_put(rtdb,'tce:read_x',mt_log,4,read_x))
     1  call errquit('tce_init: failed to write read_x',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:read_y',mt_log,4,read_y)) then
        read_y(1)=.false.
        read_y(2)=.false.
        read_y(3)=.false.
        read_y(4)=.false.
        if (.not.rtdb_put(rtdb,'tce:read_y',mt_log,4,read_y))
     1  call errquit('tce_init: failed to write read_y',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:nroots_restart',mt_int,1,
     &                  nroots_restart)) then
        nroots_restart = 0
        if (.not.rtdb_put(rtdb,'tce:nroots_restart',mt_int,1,
     1                    nroots_restart))
     2  call errquit('tce_init: failed to write nroots_restart',0,
     3               RTDB_ERR)
      endif
c
c     forces skipping of iterations during restart because for
c     huge basis sets these do not always converge in one iteration
c
      if (.not.rtdb_get(rtdb,'tce:skip_t',mt_log,1,skip_t)) then
        skip_t=.false.
        if (.not.rtdb_put(rtdb,'tce:skip_t',mt_log,1,skip_t))
     1  call errquit('tce_init: failed to write skip_t',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:skip_l',mt_log,1,skip_l)) then
        skip_l=.false.
        if (.not.rtdb_put(rtdb,'tce:skip_l',mt_log,1,skip_l))
     1  call errquit('tce_init: failed to write skip_l',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:skip_tr',mt_log,1,skip_tr)) then
        skip_tr=.false.
        if (.not.rtdb_put(rtdb,'tce:skip_tr',mt_log,1,skip_tr))
     1  call errquit('tce_init: failed to write skip_tr',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tce:skip_lr',mt_log,1,skip_lr)) then
        skip_lr=.false.
        if (.not.rtdb_put(rtdb,'tce:skip_lr',mt_log,1,skip_lr))
     1  call errquit('tce_init: failed to write skip_lr',0,RTDB_ERR)
      endif
c
c     ----------------------
c     2-e integral tolerance
c     ----------------------
c
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:itol2e',mt_int,1,itol2e))
     1    call errquit('tce_init: failed to read itol2e',0,RTDB_ERR)
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
        if(.not.((hftype.eq.'RHF').or.(hftype.eq.'ROHF'))) then
          call errquit('tce_init: 2eorb needs RHF/ROHF ref.',0,RTDB_ERR)
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
     1    (model.ne.'ccsd').and.(model.ne.'ccsdt').and.
     2    (model.ne.'qcisd').and.(model.ne.'lccsd'))
     3  call errquit('tce_init: unknown theory',0,CAPMIS_ERR)
      cperturbative = 'none'
      if (model.eq.'eionly') then
        ctype='Electron integral printing'
        if (module.ne.'tce')
     1    call errquit('tce_init: eionly module overridden',0,
     2    UNKNOWN_ERR)
      else if (model.eq.'ccd') then
        ctype='Coupled-cluster doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccd'))
     3    call errquit('tce_init: ucc/uccd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'lccd') then
        ctype='Linearized coupled-cluster doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'ulccd'))
     3    call errquit('tce_init: ucc/ulccd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'multi') then
        ctype='Multiple CCSD-like methods'
        if (module.ne.'tce')
     1    call errquit('tce_init: ucc/ulccd module overridden',0,
     2    UNKNOWN_ERR)
      else if (model.eq.'ccsd') then
        if (perturbative.eq.'none') then
          ctype='Coupled-cluster singles & doubles'
        else if (perturbative.eq.'(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = '(T)'
        else if (perturbative.eq.'lambda(t)') then
          ctype='Coupled-cluster singles & doubles w/ perturbation'
          cperturbative = 'Lambda (T)'
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
     8    call errquit('tce_init: ucc/uccsd module overridden',0,
     9    UNKNOWN_ERR)
      else if (model.eq.'lccsd') then
        ctype='Linearized coupled-cluster singles & doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'ulccsd'))
     3    call errquit('tce_init: ucc/ulccsd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsdta') then
        activecalc = .true.
        ctype='Coupled-cluster singles, doubles, & active triples'
        if (perturbative.eq.'2_q') then
          ctype='Coupled-cluster through triples w/ perturbation'
          cperturbative = '(2) quadruples (nonfactorized)'
        endif
c         if(intorb) then
c          call errquit('tce_init: 2eorb with active calc.',0,RTDB_ERR)
c         endif
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsdt'))
     3    call errquit('tce_init: ucc/uccsdt module overridden',0,
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
     3    call errquit('tce_init: ucc/uccsdt module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'ccsdtq') then
        ctype='Coupled-cluster singles, doubles, triples, & quadruples'
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uccsdtq'))
     3    call errquit('tce_init: ucc/uccsdtq module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'qcisd') then
        if (perturbative.eq.'none') then
          ctype='Quadratic configuration interaction singles & doubles'
        else if (perturbative.eq.'(t)') then
          ctype='Quadratic configuration interaction w/ perturbation'
          cperturbative = '(T)'
        else if (perturbative.eq.'cr_(t)') then
          ctype='Quadratic configuration interaction w/ perturbation'
          cperturbative = 'Completely renormalized (T)'
        endif
        if ((module.ne.'tce').and.
     1      (module.ne.'ucc').and.
     2      (module.ne.'uqcisd'))
     3    call errquit('tce_init: ucc/uqcisd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cis') then
        ctype='Configuration interaction singles'
        if((module.ne.'tce').and.
     1     (module.ne.'uci').and.
     2     (module.ne.'ucis'))
     3     call errquit('tce_init: uci/ucis module overridden',0,
     4     unknown_err)
      else if (model.eq.'cisd') then
        ctype='Configuration interaction singles & doubles'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisd'))
     3    call errquit('tce_init: uci/ucisd module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cisdt') then
        ctype='Configuration interaction singles, doubles, & triples'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisdt'))
     3    call errquit('tce_init: uci/ucisdt module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'cisdtq') then
        ctype='Configuration interaction singles through quadruples'
        if ((module.ne.'tce').and.
     1      (module.ne.'uci').and.
     2      (module.ne.'ucisdtq'))
     3    call errquit('tce_init: uci/ucisdtq module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt2') then
        ctype='Second-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt2'))
     3    call errquit('tce_init: umbpt/umbpt2 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt3') then
        ctype='Third-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt3'))
     3    call errquit('tce_init: umbpt/umbpt3 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt4') then
        ctype='Fourth-order tensor many-body perturbation theory'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt4'))
     3    call errquit('tce_init: umbpt/umbpt4 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt4sdq') then
        ctype='Fourth-order tensor many-body perturbation theory 
     1 with zero T3 approximation'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt4'))
     3    call errquit('tce_init: umbpt/umbpt4 module overridden',0,
     4    UNKNOWN_ERR)
      else if (model.eq.'mbpt4sdq_t') then
        ctype='Fourth-order tensor many-body perturbation theory
     1 with non-iterative triples'
        if ((module.ne.'tce').and.
     1      (module.ne.'umbpt').and.
     2      (module.ne.'umbpt4'))
     3    call errquit('tce_init: umbpt/umbpt4 module overridden',0,
     4    UNKNOWN_ERR)
      else
        call errquit('tce_init: unknown theory',0,CAPMIS_ERR)
      endif
c
c     -------------------------------------
c     Allocate global arrays for MO vectors
c     -------------------------------------
c
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nbf,'eigenvalues',
     1    l_evl(i),k_evl(i))) call errquit
     2    ('tce_init: failed to allocate evl',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,nbf,'occupancy',
     1    l_occ(i),k_occ(i))) call errquit
     2    ('tce_init: failed to allocate occ',0,MA_ERR)
        if (.not.ma_push_get(mt_int,nbf,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tce_init: failed to allocate irs',0,MA_ERR)
      enddo
      if (.not.movecs_read_header(movecs_out,title,basisname,
     1  scftype,nbf_verify,ipol_verify,nmo,2)) call errquit
     2  ('tce_init: failed to read MO file header',0,RTDB_ERR)
      if (nbf.ne.nbf_verify) call errquit
     1  ('tce_init: corrupted MO vectors',0,RTDB_ERR)
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
     2      call errquit('tce_init: failed to read MO vectors',0,
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
     2    call errquit('tce_init: no alpha irreps found',0,RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_init: no alpha irreps found',1,RTDB_ERR)
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
     2      call errquit('tce_init: no beta irreps found',0,RTDB_ERR)
        else if (reference.eq.1) then
          if (rohf) then
            do j = 1,nbf
              int_mb(k_irs(2)+j-1) = int_mb(k_irs(1)+j-1) + 1
            enddo
          else
            if (.not.rtdb_get(rtdb,'scf:beta irreps',mt_int,
     1        nbf,int_mb(k_irs(2))))
     2        call errquit('tce_init: no beta irreps found',1,RTDB_ERR)
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
     1  call errquit('tce_init: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1  call errquit('tce_init: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1  call errquit('tce_init: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1  call errquit('tce_init: failed to read oactive',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1  call errquit('tce_init: failed to read numact',0,RTDB_ERR)
      if (restricted) then
        oactive(2) = oactive(1)
        vactive(2) = vactive(1)
      endif
      do i = 1, ipol
        if (oactive(i).gt.nocc(i)-nfc(i)) call errquit
     1  ('tce_init: reduce the active space',0,RTDB_ERR)
        if (vactive(i).gt.nmo(i)-nocc(i)-nfv(i)) call errquit
     1  ('tce_init: reduce the active space',0,RTDB_ERR)
      enddo
        if (numact.le.0.or.numact.gt.3) call errquit
     1  ('tce_init: change active excit. level: t3a_lvl',0,RTDB_ERR)
      do i = 1, ipol
        if (oactive(i).le.0) call errquit
     1  ('tce_init: incorrect size of active space',0,RTDB_ERR)
        if (vactive(i).le.0) call errquit
     1  ('tce_init: incorrect size of active space',0,RTDB_ERR)
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
     1    dfloat(buffer_size)*dfloat(bytes)).ne.0)
     2    call errquit('tce_init: dra initialization failed',0,
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
        write(LuOut,9115) 'T(0)',zlshift
        write(LuOut,9115) 'L(0)',zlshiftl
        write(LuOut,9115) 'T(1)',zlshift2(1)
        write(LuOut,9115) 'L(1)',zlshift2(2)
        write(LuOut,9115) 'T(R)',zlshift3(1)
        write(LuOut,9115) 'T(I)',zlshift3(2)
        if (diis .gt. 0) then
          write(LuOut,9150) diis,'-th order DIIS'
        else
          write(LuOut,9160) 'Jacobi'
        endif
        if (lineresp) then
          if (diis2 .gt. 0) then
            write(LuOut,9152) diis2,'-th order DIIS'
          else
            write(LuOut,9162) 'Jacobi'
          endif
        endif
        if (imagresp) then
          if (diis3 .gt. 0) then
            write(LuOut,9153) diis3,'-th order DIIS'
          else
            write(LuOut,9163) 'Jacobi'
          endif
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
            if (multipole.ge.2) then
              write(LuOut,9221) 'on '
            endif
            if (multipole.ge.3) then
              write(LuOut,9222) 'on '
            endif
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
          call errquit('tce_init: unknown I/O scheme',ioalg,CAPMIS_ERR)
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
 9115 format(5x,a4,' DIIS level shift : ',d8.2)
c 9116 format(3x,' CC-LR DIIS level shift : ',d8.2)
c 9117 format(3x,' CC-IR DIIS level shift : ',d8.2,' (imag part)')
c 9118 format(3x,' CC-IR DIIS level shift : ',d8.2,' (real part)')
 9130 format(3x,'         Max iterations :',2x,i7)
 9140 format(3x,'             I/O scheme : ',a)
 9150 format(3x,'CC-T/L Amplitude update : ',i2,a)
 9152 format(3x,' CC-LR Amplitude update : ',i2,a)
 9153 format(3x,' CC-IR Amplitude update : ',i2,a)
 9160 format(3x,'       Amplitude update : ',a)
 9162 format(3x,' CC-LR Amplitude update : ',a)
 9163 format(3x,' CC-IR Amplitude update : ',a)
 9170 format(3x,'  No. of excited states :',2x,i4)
 9180 format(3x,'            Target root :',2x,i4)
 9190 format(3x,'        Target symmetry :',1x,a4)
 9210 format(3x,'   Symmetry restriction :',1x,a3)
 9220 format(3x,'Dipole & oscillator str :',1x,a3)
 9221 format(3x,'...same for quadrupoles :',1x,a3)
 9222 format(3x,'  ...same for octupoles :',1x,a3)
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
#ifndef CCSDTLR
      if ((model.eq.'ccsdt').and.(lineresp)) then
        call errquit('CCSDTLR not compiled (tce_init.F)',0,CAPMIS_ERR)
        return
      endif
#endif
      return
      end
