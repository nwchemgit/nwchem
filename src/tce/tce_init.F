      subroutine tce_init(rtdb)
c
c $Id: tce_init.F,v 1.22 2003-01-08 22:59:13 sohirata Exp $
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
#include "tce.fh"
c
      integer rtdb             ! RTDB handle
      integer numfunc          ! Number of functionals implemented
      parameter (numfunc=30)   ! See if it is consistent with nwdft
      integer nbf_verify       ! nbf read from movecs_read_header
      integer ipol_verify      ! ipol read from movecs_read_header
      character*255 movecs_out ! MO vector file name
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      character*4 hftype
      character*3 onoff1,onoff2
      integer i,j,n
      integer nclosed, nopen
      double precision a
      character*60 ctype
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
      character*8 multchar(8)
      data multchar/'singlet ','doublet ','triplet ','quartet ',
     1              'quintet ','sextet  ','septet  ','octet   '/
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
     1  call errquit('tce_init: geom_create failed',0)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('tce_init: no geometry found',0)
      if (.not.rtdb_get(rtdb,'dft:skeleton',mt_log,1,oskel)) then
        oskel=sym_number_ops(geom).gt.0
      endif
      if (.not.rtdb_get(rtdb,'dft:adapt',mt_log,1,oadapt)) then
        oadapt=sym_number_ops(geom).gt.0
      endif
      if (.not.sym_abelian_group(geom))
     1  call errquit('tce_init: non-Abelian symmetry not permitted',0)
c
c     ---------
c     Basis set
c     ---------
c
      if (.not.bas_create(ao_bas_han,'ao basis'))
     1  call errquit('tce_init: bas_create failed',0)
      if (.not.bas_rtdb_load(rtdb,geom,ao_bas_han,'ao basis'))
     1  call errquit('tce_init: no ao basis found',0)
      if (.not.int_normalize(rtdb,ao_bas_han))
     1  call errquit('tce_init: int_normalize failed',0)
      if (.not.bas_numbf(ao_bas_han,nbf))
     1  call errquit('tce_init: failed to read nbf',0)
      if (.not.bas_numcont(ao_bas_han,nshells))
     1  call errquit('tce_init: failed to read nshells',0)
c
c     ----------
c     MO vectors
c     ----------
c
      if (reference.eq.0) then
        if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_out))
     1    call errquit('tce_init: no final MO found',0)
      else if (reference.eq.1) then
        if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     1    call errquit('tce_init: no final MO found',0)
      endif
c
c     ------------------------------------------------------
c     Number of electrons and orbitals and spin multiplicity
c     ------------------------------------------------------
c
      if (reference.eq.0) then
        rohf = .false.
        if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     1    call errquit('tce_init: failed to read ipol',0)
        if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc))
     1    call errquit('tce_init: failed to read nocc',0)
        if (.not.rtdb_get(rtdb,'dft:mult',mt_int,1,mult))
     1    call errquit('tce_init: failed to read mult',0)
        if (ipol.eq.1) then
          restricted = .true.
        else
          restricted = .false.
        endif
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nclosed))
     1    call errquit('tce_init: failed to read nclosed',0)
        if (.not.rtdb_get(rtdb,'scf:nopen',mt_int,1,nopen))
     1    call errquit('tce_init: failed to read nopen',0)
        nocc(1) = nclosed + nopen
        nocc(2) = nclosed
        mult = nopen
        if(.not.rtdb_cget(rtdb,'scf:scftype',1,hftype))
     1    call errquit('tce_init: failed to read scftype',0)
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
c
c     --------------
c     XC functionals
c     --------------
c
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:xfac',mt_dbl,numfunc,xfac))
     1    call errquit('tce_init: failed to read xfac',0)
        if (.not.rtdb_get(rtdb,'dft:cfac',mt_dbl,numfunc,cfac))
     1    call errquit('tce_init: failed to read cfac',0)
      else if (reference.eq.1) then
        do i = 1,numfunc
          xfac(i) = 0.0d0
          cfac(i) = 0.0d0
        enddo
      endif
c
c     ---------
c     TCE input
c     ---------
c
      if (rtdb_get(rtdb,'tce:frozen core:freeze by atoms',
     1  mt_log,1,dummy)) then
        if (.not.geom_num_core(geom,nfc(1)))
     1    call errquit('tce_init: failed to get number of cores',0)
      else if (rtdb_get(rtdb,'cc:frozen core',
     1  mt_int,1,nfc(1))) then
        continue
      else
        nfc(1) = 0
      endif
      if (ipol.eq.2) nfc(2) = nfc(1)
      if (.not.rtdb_get(rtdb,'tce:frozen virtual',
     1  mt_int,1,nfv(1))) nfv(1) = 0
      if (ipol.eq.2) nfv(2) = nfv(1)
      if (.not.rtdb_cget(rtdb,'tce:model',1,model)) then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1  call errquit('tce_init: failed to write model',0)
      endif
      if (.not.rtdb_get(rtdb,'tce:thresh',mt_dbl,1,thresh)) then
        thresh = 1.0d-6
        if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_init: failed to write thresh',0)
      endif
      if (.not.rtdb_get(rtdb,'tce:maxiter',mt_int,1,maxiter)) then
        maxiter = 100
        if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_init: failed to write maxiter',0)
      endif
      if (.not.rtdb_get(rtdb,'tce:ioalg',mt_int,1,ioalg)) then
        ioalg = 2
        if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_init: failed to write ioalg',0)
      endif
      if (.not.rtdb_get(rtdb,'tce:diis',mt_int,1,diis)) then
        diis = 5
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init: failed to write diis',0)
      endif
      if (diis .gt. maxdiis) then
        diis = maxdiis
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_init: failed to write diis',0)
      endif
c
c     ----------------------
c     2-e integral tolerance
c     ----------------------
c
      if (reference.eq.0) then
        if (.not.rtdb_get(rtdb,'dft:itol2e',mt_int,1,itol2e))
     1    call errquit('tce_init: failed to read itol2e',0)
        tol2e=10.0d0**(-itol2e)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:tol2e',mt_dbl,1,tol2e))
     1    tol2e = 1.0d-10
      endif
c
c     --------------------------
c     Determine calculation type
c     --------------------------
c
      if (reference.eq.0) then
        a=0.0d0
        do n=2,numfunc
          a=a+dabs(xfac(n))
        enddo
        do n=1,numfunc
          a=a+dabs(cfac(n))
        enddo
        if ((a.gt.1.0d-8).or.(dabs(xfac(1)-1.0d0)).gt.1.0d-8)
     1    write(LuOut,9500)
      endif
      if (model.eq.'ccd') then
        ctype='Coupled-cluster doubles'
      else if (model.eq.'lccd') then
        ctype='Linearized coupled-cluster doubles'
      else if (model.eq.'ccsd') then
        ctype='Coupled-cluster singles & doubles'
      else if (model.eq.'lccsd') then
        ctype='Linearized coupled-cluster singles & doubles'
      else if (model.eq.'ccsdt') then
        ctype='Coupled-cluster singles, doubles, & triples'
      else if (model.eq.'ccsdtq') then
        ctype='Coupled-cluster singles, doubles, triples, & quadruples'
      else if (model.eq.'qcisd') then
        ctype='Quadratic configuration interaction singles & doubles'
      else if (model.eq.'cisd') then
        ctype='Configuration interaction singles & doubles'
      else if (model.eq.'cisdt') then
        ctype='Configuration interaction singles, doubles, & triples'
      else if (model.eq.'cisdtq') then
        ctype='Configuration interaction singles through quadruples'
      else if (model.eq.'mbpt2') then
        ctype='Second-order tensor many-body perturbation theory'
      else if (model.eq.'mbpt3') then
        ctype='Third-order tensor many-body perturbation theory'
      else if (model.eq.'mbpt4') then
        ctype='Fourth-order tensor many-body perturbation theory'
      else
        call errquit('tce_init: unknown theory',0)
      endif
c
c     -------------------------------------
c     Allocate global arrays for MO vectors
c     -------------------------------------
c
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nbf,'eigenvalues',
     1    l_evl(i),k_evl(i))) call errquit
     2    ('tce_init: failed to allocate evl',0)
        if (.not.ma_push_get(mt_dbl,nbf,'occupancy',
     1    l_occ(i),k_occ(i))) call errquit
     2    ('tce_init: failed to allocate occ',0)
        if (.not.ma_push_get(mt_int,nbf,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tce_init: failed to allocate irs',0)
      enddo
      if (.not.movecs_read_header(movecs_out,title,basisname,
     1  scftype,nbf_verify,ipol_verify,nmo,2)) call errquit
     2  ('tce_init: failed to read MO file header',0)
      if (nbf.ne.nbf_verify) call errquit
     1  ('tce_init: corrupted MO vectors',0)
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
     2      call errquit('tce_init: failed to read MO vectors',0)
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
     2    call errquit('tce_init: no alpha irreps found',0)
      else if (reference.eq.1) then
        if (.not.rtdb_get(rtdb,'scf:alpha irreps',mt_int,
     1    nbf,int_mb(k_irs(1))))
     2    call errquit('tce_init: no alpha irreps found',1)
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
     2      call errquit('tce_init: no beta irreps found',0)
        else if (reference.eq.1) then
          if (rohf) then
            do j = 1,nbf
              int_mb(k_irs(2)+j-1) = int_mb(k_irs(1)+j-1) + 1
            enddo
          else
            if (.not.rtdb_get(rtdb,'scf:beta irreps',mt_int,
     1        nbf,int_mb(k_irs(2))))
     2        call errquit('tce_init: no beta irreps found',1)
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
     1  (LuOut,'generated by the Tensor Contraction Engine (TCE).',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'The development of this portion of the program',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'and TCE was supported by US DOE Basic Energy Sciences.',
     2  40,.false.)
        call util_print_centered
     1  (LuOut,'TCE is a product of Battelle and PNNL.',
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
        write(LuOut,9420) multchar(mult+1)
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
        write(LuOut,9130) maxiter
        write(LuOut,9110) thresh
        if (diis .gt. 0) then
          write(LuOut,9150) diis,'-th order DIIS'
        else
          write(LuOut,9160) 'Jacobi'
        endif
        if (ioalg.eq.0) then
          write(LuOut,9140) 'Global files (Fortran77 direct access)'
        else if (ioalg.eq.1) then
          write(LuOut,9140) 'Global files (C low-level I/O)'
        else if (ioalg.eq.2) then
          write(LuOut,9140) 'Global files (Global array library)'
        else if (ioalg.eq.3) then
          write(LuOut,9140) 'Global files (Shared file library)'
        else if (ioalg.eq.4) then
          write(LuOut,9140) 'Distributed files (C low-level I/O)'
        else
          call errquit('tce_init: unknown I/O scheme',ioalg)
        endif
      endif
      if (nodezero) call util_flush(LuOut)
 9200 format(3x,'   Number of processors :',2x,i4)
 9000 format(3x,'       Calculation type : ',a)
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
 9130 format(3x,'         Max iterations :',2x,i7)
 9140 format(3x,'             I/O scheme : ',a)
 9150 format(3x,'       Amplitude update : ',i2,a)
 9160 format(3x,'       Amplitude update : ',a)
 9400 format(/,2x,
     1  '!! WARNING !! The number of MO is less than the number of AO')
 9500 format(/,2x,
     1  '!! WARNING !! A non-HF reference wave function')
c
c     ------
c     Return
c     ------
c
      return
      end
