      subroutine tddft_init(rtdb,geom,oskel,oadapt,ao_bas_han,
     1  nbf_ao,movecs_out,ipol,nocc,nmo,nfc,nfv,numfunc,
     2  xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,nlcfac,xccomb,
     3  lhf,ldft,lhybrid,lhasdft,lhashf,
     4  tda,singlet,triplet,nroots,maxvecs,tol2e,thresh,
     5  maxiter,target,targetsym,symmetry,algorithm,g_movecs,
     6  l_evl,k_evl,l_occ,k_occ,l_irs,k_irs,mult)
c
c $Id: tddft_init.F,v 1.13 2004-10-27 21:18:42 edo Exp $
c
c Set various parameters for TDDFT calculations, which are
c passed as arguments to tddft_init.  Most are read from
c runtime database and are originally set by NWDFT module.
c See also dft_rdinput.F and dft_inpana.F.
c Written by So Hirata, Feb 2002.
c
      implicit none
#include "errquit.fh"
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
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      logical oskel,oadapt     ! Symmetry control parameters
      integer nbf_ao           ! Number of AO basis functions
      integer nbf_ao_verify    ! nbf_ao read from movecs_read_header
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer ipol_verify      ! ipol read from movecs_read_header
      integer numfunc          ! Number of functionals implemented
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nmix(2)          ! Number of symmetry contaminations
      character*255 movecs_out ! MO vector file name
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
      character*40 xname(numfunc)       ! Names of X functionals
      character*40 cname(numfunc)       ! Names of C functionals
      character*40 xcname(numfunc)      ! Names of XC functionals
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      logical lxfac(numfunc)            ! True if local X functional
      logical nlxfac(numfunc)           ! True if gc X functional
      logical lcfac(numfunc)            ! True if local C functional
      logical nlcfac(numfunc)           ! True if gc C functional
      logical xccomb(numfunc)           ! True if composite XC functional
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
      logical tda                       ! True if Tamm-Dancoff approximation
      logical singlet                   ! True if singlet excited states
      logical triplet                   ! True if triplet excited states
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      double precision tol2e   ! 2-e integral tolerance
      double precision thresh  ! Convergence threshold for Davidson iteration
      integer maxiter          ! Maximum number of Davidson iterations
      integer target           ! Target excited state
      character*4 targetsym    ! Symmetry of target excited state
      logical symmetry         ! True if symmetry restriction
      integer algorithm        ! =1 multiple fock build; =2: conventional
      integer g_movecs(2)      ! GA handle for MO vectors
      integer l_evl,k_evl      ! MA handles for eigenvalues (AB)
      integer l_occ,k_occ      ! MA handles for occupancies (AB)
      integer l_irs(2),k_irs(2)! MA handles for irreps (AB)
      integer mult             ! Ground state spin multiplicity
      logical lb94             ! LB94 asymptotic correction
      logical cs00             ! CS00 asymptotic correction
      double precision shift   ! Shift for CS00
C     character*80 vector      ! CI vector filename
c
      character*3 onoff1,onoff2
      integer i,j,n
      integer g_work
      double precision a
      character*18 ctype
      character*4 oname
      integer itol2e
      logical nodezero
      logical dummy
c
      logical int_normalize
      external int_normalize
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      nodezero=(ga_nodeid().eq.0)
c
c ---------------------
c Geometry and symmetry
c ---------------------
      if (.not.geom_create(geom,'geometry'))
     1  call errquit('tddft_init: geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('tddft_init: no geometry found',0, RTDB_ERR)
      if (.not.sym_abelian_group(geom))
     1  call errquit('tddft_init: non-Abelian symmetry not permitted',0,
     2  GEOM_ERR)

      if (.not.rtdb_get(rtdb,'dft:skeleton',mt_log,1,oskel)) then
        oskel=.false.
      endif
      if (oskel) call errquit
     1  ('tddft_init: turn off skeleton symmetry',0, INPUT_ERR)
      if (.not.rtdb_get(rtdb,'dft:adapt',mt_log,1,oadapt)) then
        oadapt=sym_number_ops(geom).gt.0
      endif
c ---------
c Basis set
c ---------
      if (.not.bas_create(ao_bas_han,'ao basis'))
     1  call errquit('tddft_init: bas_create failed',0, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,ao_bas_han,'ao basis'))
     1  call errquit('tddft_init: no ao basis found',0, BASIS_ERR)
      if (.not.int_normalize(rtdb,ao_bas_han))
     1  call errquit('tddft_init: int_normalize failed',0, BASIS_ERR)
      if (.not.bas_numbf(ao_bas_han,nbf_ao))
     1  call errquit('tddft_init: failed to read nbf_ao',0, BASIS_ERR)
c ----------
c MO vectors
c ----------
      if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_out))
     1  call errquit('tddft_init: no final MO found',0, RTDB_ERR)
c ------------------------------------------------------
c Number of electrons and orbitals and spin multiplicity
c ------------------------------------------------------
      if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     1  call errquit('tddft_init: failed to read ipol',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc))
     1  call errquit('tddft_init: failed to read nocc',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:mult',mt_int,1,mult))
     1  call errquit('tddft_init: failed to read mult',0, RTDB_ERR)
c --------------
c XC functionals
c --------------
      if (.not.rtdb_get(rtdb,'dft:xfac',mt_dbl,numfunc,xfac))
     1  call errquit('tddft_init: failed to read xfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:cfac',mt_dbl,numfunc,cfac))
     1  call errquit('tddft_init: failed to read cfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lxfac',mt_log,numfunc,lxfac))
     1  call errquit('tddft_init: failed to read lxfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:nlxfac',mt_log,numfunc,nlxfac))
     1  call errquit('tddft_init: failed to read nlxfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lcfac',mt_log,numfunc,lcfac))
     1  call errquit('tddft_init: failed to read lcfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:nlcfac',mt_log,numfunc,nlcfac))
     1  call errquit('tddft_init: failed to read nlcfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:xccomb',mt_log,numfunc,xccomb))
     1  call errquit('tddft_init: failed to read xccomb',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:cname',numfunc,cname))
     1  call errquit('tddft_init: failed to read xname',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:xcname',numfunc,xcname))
     1  call errquit('tddft_init: failed to read cname',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:xname',numfunc,xname))
     1  call errquit('tddft_init: failed to read xcname',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lb94',mt_log,1,lb94))
     1  call errquit('tddft_init: failed to read lb94',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:cs00',mt_log,1,cs00))
     1  call errquit('tddft_init: failed to read cs00',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:delta_ac',mt_dbl,1,shift))
     1  call errquit('tddft_init: failed to read delta_ac',0, RTDB_ERR)
c -----------
c TDDFT input
c -----------
      if(.not. geom_freeze(rtdb, geom, 'tddft', nfc(1))) then
        call errquit('tddft_init: failed to get number of cores', 0,
     &       GEOM_ERR)
      endif

      if (ipol.eq.2) nfc(2) = nfc(1)
      if (.not.rtdb_get(rtdb,'tddft:frozen virtual',
     1  mt_int,1,nfv(1))) nfv(1) = 0
      if (ipol.eq.2) nfv(2) = nfv(1)
      if (.not.rtdb_get(rtdb,'tddft:nroots',mt_int,1,nroots)) then
        nroots=1
        if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1  call errquit('tddft_init: failed to write nroots',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:maxvecs',mt_int,1,maxvecs)) then
        maxvecs=1000
        if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1  call errquit('tddft_init: failed to write maxvecs',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:tda',mt_log,1,tda)) then
        tda=.false.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1  call errquit('tddft_init: failed to write tda',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:singlet',mt_log,1,singlet)) then
        singlet=.true.
        if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1  call errquit('tddft_init: failed to write singlet',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:triplet',mt_log,1,triplet)) then
        triplet=.true.
        if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1  call errquit('tddft_init: failed to write triplet',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:thresh',mt_dbl,1,thresh)) then
        thresh=1.0d-4
        if (.not.rtdb_put(rtdb,'tddft:thresh',mt_dbl,1,thresh))
     1  call errquit('tddft_init: failed to write thresh',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:maxiter',mt_int,1,maxiter)) then
        maxiter=1000
        if (.not.rtdb_put(rtdb,'tddft:maxiter',mt_int,1,maxiter))
     1  call errquit('tddft_init: failed to write maxiter',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:target',mt_int,1,target)) then
         target=1
        if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,target))
     1  call errquit('tddft_init: failed to write target',0, RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'tddft:targetsym',1,targetsym)) then
        targetsym='none'
        if (.not.rtdb_cput(rtdb,'tddft:targetsym',1,targetsym))
     1  call errquit('tddft_init: failed to write targetsym',0,
     &       RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:symmetry',mt_log,1,symmetry)) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tddft:symmetry',mt_log,1,symmetry))
     1  call errquit('tddft_init: failed to write symmetry',0, RTDB_ERR)
      endif
      if ((symmetry).and.(targetsym.eq.'none')) then
        symmetry=.false.
        if (.not.rtdb_put(rtdb,'tddft:symmetry',mt_log,1,symmetry))
     1    call errquit('tddft_init: failed to write symmetry',0,
     &       RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:algorithm',mt_int,1,algorithm)) then
        algorithm=0
        if (.not.rtdb_put(rtdb,'tddft:algorithm',mt_int,1,algorithm))
     1  call errquit('tddft_init: failed to write algorithm',0,
     &       RTDB_ERR)
      endif
C     if (.not.rtdb_cget(rtdb,'tddft:vector',1,vector)) then
C       call util_file_name('tddft',.false.,.false.,vector)
C       if (.not.rtdb_cput(rtdb,'tddft:vector',1,vector))
C    1  call errquit('tddft_init: failed to write vector',0, RTDB_ERR)
C     endif
c ----------------------
c 2-e integral tolerance
c ----------------------
      if (.not.rtdb_get(rtdb,'dft:itol2e',mt_int,1,itol2e))
     1  call errquit('tddft_init: failed to read itol2e',0, RTDB_ERR)
      tol2e=10.0d0**(-itol2e)
c ----------------------
c Initialize integrators
c ----------------------
      call int_init(rtdb,1,ao_bas_han)
      call schwarz_init(geom,ao_bas_han)
      call scf_get_fock_param(rtdb,tol2e)
c --------------------------
c Determine calculation type
c --------------------------
      a=0.0d0
      do n=2,numfunc
        a=a+dabs(xfac(n))
      enddo
      do n=1,numfunc
        a=a+dabs(cfac(n))
      enddo
      if ((a.gt.1.0d-8).and.(dabs(xfac(1)).le.1.0d-8)) then
c Pure DFT
        ldft=.true.
        lhf=.false.
        lhybrid=.false.
        lhasdft=.true.
        lhashf=.false.
      else if ((a.gt.1.0d-8).and.(dabs(xfac(1)).gt.1.0d-8)) then
c Hybrid DFT & HF
        ldft=.false.
        lhf=.false.
        lhybrid=.true.
        lhasdft=.true.
        lhashf=.true.
      else if ((a.le.1.0d-8).and.(dabs(xfac(1)-1.0d0)).le.1.0d-8) then
c Pure HF
        ldft=.false.
        lhf=.true.
        lhybrid=.false.
        lhasdft=.false.
        lhashf=.true.
      else
c Screened HF
        ldft=.false.
        lhf=.false.
        lhybrid=.false.
        lhasdft=.false.
        lhashf=.true.
      endif
      if (tda) then
        if (ldft.or.lhybrid) then
          ctype='Tamm-Dancoff TDDFT'
        else if (lhf) then
          ctype='CIS               '
        else
          ctype='Screened CIS      '
        endif
      else
        if (ldft.or.lhybrid) then
          ctype='TDDFT             '
        else if (lhf) then
          ctype='TDHF              '
        else
          ctype='Screened TDHF     '
        endif
      endif
c -------------------------------------
c Allocate global arrays for MO vectors
c -------------------------------------
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,'eigenvalues',
     1  l_evl,k_evl)) call errquit
     2  ('tddft_init: failed to allocate evl',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,'occupancy',
     1  l_occ,k_occ)) call errquit
     2  ('tddft_init: failed to allocate occ',0, MA_ERR)
      if (.not.movecs_read_header(movecs_out,title,basisname,
     1  scftype,nbf_ao_verify,ipol_verify,nmo,2)) call errquit
     2  ('tddft_init: failed to read MO file header',0, DISK_ERR)
      if (nbf_ao.ne.nbf_ao_verify) call errquit
     1  ('tddft_init: corrupted MO vectors',0, DISK_ERR)
      if (ipol.ne.ipol_verify) call errquit
     1  ('tddft_init: corrupted MO vectors',0, DISK_ERR)
      do i=1,ipol
        g_movecs(i)=ga_create_atom_blocked
     1    (geom,ao_bas_han,'MO eigenvectors')
        if (.not.movecs_read(movecs_out,i,
     1    dbl_mb(k_occ+(i-1)*nbf_ao),
     2    dbl_mb(k_evl+(i-1)*nbf_ao),g_movecs(i))) call errquit
     3    ('tddft_init: failed to read MO vectors',0, DISK_ERR)
        if (util_print('MO vectors',print_debug)) then
          if (nodezero) call ma_print(dbl_mb(k_evl+(i-1)*nbf_ao),
     1      nbf_ao,1,'Eigenvalues')
          if (nodezero) call ma_print(dbl_mb(k_occ+(i-1)*nbf_ao),
     1      nbf_ao,1,'Occupancies')
          call ga_print(g_movecs(i))
        endif
      enddo
c ---------------------------------------
c Determine symmetry (irreps) of orbitals
c ---------------------------------------
      do i=1,ipol
c
c We re-symmetry-adapt the MO's to get irreps information
c Make sure not to alter g_movecs
c
        if (.not.ma_push_get(mt_int,nbf_ao,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tddft_init: failed to allocate irs',0, MA_ERR)
        if (.not.ga_duplicate(g_movecs(i),g_work,'Work')) call errquit
     1    ('tddft_init: failed to duplicate g_movecs',0, GA_ERR)
        call ga_copy(g_movecs(i),g_work)
        call sym_movecs_adapt(ao_bas_han,1.0d-8,g_work,
     1    int_mb(k_irs(i)),nmix(i))
        if (.not.ga_destroy(g_work)) call errquit
     1    ('tddft_init: failed to destroy g_work',0, GA_ERR)
        if (nodezero.and.util_print('MO vectors',print_debug)) then
          write(LuOut,*) 'Spin ',i
          do j=1,nbf_ao
            call sym_irrepname(geom,int_mb(k_irs(i)+j-1),oname)
            write(LuOut,*) ' orbital',j,oname
          enddo
        endif
      enddo
c --------------
c Adjust maxvecs
c --------------
      maxvecs=min(maxvecs,(nroots+20)*maxiter*2)
      j=0
      do i=1,ipol
        j=j+(nmo(i)-nocc(i))*nocc(i)
      enddo
      maxvecs=min(maxvecs,j)
c -------------
c Printing info
c -------------
      if (nodezero.and.util_print('information',print_low)) then
        call util_print_centered
     1  (LuOut,'NWChem TDDFT Module',40,.true.)
        write(LuOut,*)
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
        write(LuOut,9420) mult
        write(LuOut,9050) nbf_ao
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
          if (nmo(1).lt.nbf_ao) write(LuOut,9400)
        else
          if ((nmo(1).lt.nbf_ao).or.(nmo(2).lt.nbf_ao))
     1      write(LuOut,9400)
        endif
      endif
      if (nodezero.and.util_print('xc information',print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'XC Information',20,.true.)
        do n=1,numfunc
          if (xccomb(n)) write(LuOut,9300) xcname(n)
        enddo
        if (lxfac(1).or.nlxfac(1))
     1    write(LuOut,9310) xname(1),xfac(1),'         '
        do n=2,numfunc
          if (lxfac(n).and.nlxfac(n)) then
            write(LuOut,9310) xname(n),xfac(n),'         '
          else if (lxfac(n).and.(.not.nlxfac(n))) then
            write(LuOut,9310) xname(n),xfac(n),'local    '
          else if ((.not.lxfac(n)).and.nlxfac(n)) then
            write(LuOut,9310) xname(n),xfac(n),'non-local'
          endif
        enddo
        do n=1,numfunc
          if (lcfac(n).and.nlcfac(n)) then
            write(LuOut,9310) cname(n),cfac(n),'         '
          else if (lcfac(n).and.(.not.nlcfac(n))) then
            write(LuOut,9310) cname(n),cfac(n),'local    '
          else if ((.not.lcfac(n)).and.nlcfac(n)) then
            write(LuOut,9310) cname(n),cfac(n),'non-local'
          endif
        enddo
        if (cs00) then
          if (shift.gt.1.0d90) then
          write(LuOut,*)
          write(LuOut,9330)
     1    '  CS with a Zhan-Nichols-Dixon shift        '
          else
          write(LuOut,*)
          write(LuOut,9320)
     1    '  Casida-Salahub correction with a shift    ',
     2    shift,'au'
          endif
        else if (lb94) then
          write(LuOut,*)
          write(LuOut,9330)
     1    '         van Leeuwen-Baerends correction    '
        endif
      endif
      if (nodezero.and.util_print('xc information',print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'TDDFT Information',20,.true.)
        write(LuOut,9000) ctype  
        if (ipol.eq.1) then
          if ((singlet).and.(.not.triplet)) write(LuOut,9080)
          if ((triplet).and.(.not.singlet)) write(LuOut,9090)
          if ((singlet).and.(triplet)) write(LuOut,9170)
          write(LuOut,9030) nocc(1)*2,nocc(1),nocc(1)
        else if (ipol.eq.2) then
          write(LuOut,9020)
          write(LuOut,9030) nocc(1)+nocc(2),nocc(1),nocc(2)
        endif
        write(LuOut,9060) nroots
        write(LuOut,9070) maxvecs
        write(LuOut,9130) maxiter
        write(LuOut,9120) target
        write(LuOut,9150) targetsym
        if (symmetry) then
          write(LuOut,9160) 'on '
        else
          write(LuOut,9160) 'off'
        endif
        if (algorithm.eq.0) then
          write(LuOut,9140) 'Optimal'
        else if (algorithm.eq.1) then
          write(LuOut,9140) 'Incore multiple tensor contraction'
        else if (algorithm.eq.2) then
          write(LuOut,9140) 'Incore single tensor contraction'
        else if (algorithm.eq.3) then
          write(LuOut,9140) 'Disk-based multiple tensor contraction'
        else if (algorithm.eq.4) then
          write(LuOut,9140) 'Disk-based single tensor contraction'
        else
          write(LuOut,9140) 'Unknown'
        endif
        write(LuOut,9110) thresh
      endif
      if (nodezero) call util_flush(LuOut)
 9000 format(3x,'       Calculation type : ',a18)
 9080 format(3x,'      Wavefunction type : Restricted singlets')
 9090 format(3x,'      Wavefunction type : Restricted triplets')
 9170 format(3x,
     1'      Wavefunction type : Restricted singlets & triplets')
 9020 format(3x,'      Wavefunction type : Unrestricted')
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
 9420 format(3x,'      Spin multiplicity :',2x,i4)
 9040 format(3x,'     Use of symmetry is : ',a3,/,
     1       3x,'   Symmetry adaption is : ',a3)
 9050 format(3x,' Number of AO functions :',2x,i4)
 9060 format(3x,'           No. of roots :',2x,i4)
 9070 format(3x,'       Max subspacesize :',2x,i4)
 9100 format(3x,'      Schwarz screening : ',d8.2)
 9110 format(3x,'     Davidson threshold : ',d8.2)
 9120 format(3x,'            Target root :',2x,i4)
 9150 format(3x,'        Target symmetry :',1x,a4)
 9160 format(3x,'   Symmetry restriction :',1x,a3)
 9130 format(3x,'         Max iterations :',2x,i4)
 9140 format(3x,'              Algorithm : ',a)
 9300 format(3x,a40)
 9310 format(3x,a40,1x,f6.2,1x,a9)
 9320 format(3x,a44,1x,f10.6,1x,a2)
 9330 format(3x,a44)
 9400 format(/,2x,
     1  '!! WARNING !! The number of MO is less than the number of AO')
c
c Return
c
      return
      end
