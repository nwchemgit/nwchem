       logical function tce_mrcc_energy(rtdb)
c
c Parallel version of multireference CC interface
c Jiri Brabec 2011
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
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_hetio.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"

#ifdef MRCC_METHODS
#include "tce_mrcc.fh"
#endif

      integer rtdb
      logical nodezero

      character*255 modelname
      character*255 filename
      double precision ref     ! Ground state energy
      logical dft_energy
      external dft_energy
      logical scf
      external scf
      logical mcscf
      external mcscf
      integer g_ao1e(2)        ! GA handle for AO Fock matrices
      integer d_ao2e           ! SF handle for AO 2e integrals
ckbn      integer i,j,iref,nref,yref
      integer i,j,iref,yref
c      logical mrccdebug
      integer l_f1_offsetm(maxref),k_f1_offsetm(maxref)
      integer l_v2_offsetm(maxref),k_v2_offsetm(maxref)
      integer size_1em(maxref)
      integer size_2em(maxref)
      logical int2e_file_close
      external int2e_file_close
      character*3 namechunk
      character*3 namechunk2
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      double precision itcpu,itcpu2     ! CPU sec counter
      double precision itwall,itwall2    ! WALL sec counter
      integer d_v2m(maxref)
      logical needt1
      logical needt2
      logical needt3
      logical needt3a
      logical needr3act
      logical needr4act
      integer d_t1m(maxref)             ! SF handle for t1 amplitudes
      integer l_t1_offsetm(maxref)      ! Offset for t1 file
      integer k_t1_offsetm(maxref)      ! Offset for t1 file
      integer size_t1m(maxref)          ! File size in doubles
      integer d_t2m(maxref)             ! SF handle for t2 amplitudes
      integer d_t3m(maxref)             ! SF handle for t3 amplitudes
      integer l_t2_offsetm(maxref)      ! Offset for t2 file
      integer k_t2_offsetm(maxref)      ! Offset for t2 file
      integer size_t2m(maxref)          ! File size in doubles
      integer l_t3_offsetm(maxref)      ! Offset for t3 file
      integer k_t3_offsetm(maxref)      ! Offset for t3 file
      integer size_t3m(maxref)          ! File size in doubles
      integer l_r1_offsetm(maxref),k_r1_offsetm(maxref)
      integer l_r2_offsetm(maxref),k_r2_offsetm(maxref)
      integer l_r3_offsetm(maxref),k_r3_offsetm(maxref)
      integer size_r3m(maxref),d_r3m(maxref)
      integer d_r3am(maxref)
      integer size_r3am(maxref)
      integer l_r3a_offsetm(maxref),k_r3a_offsetm(maxref)
      integer d_r4am(maxref)
      integer size_r4am(maxref)
      integer l_r4a_offsetm(maxref),k_r4a_offsetm(maxref)
      integer d_r2actm(maxref)
      integer size_r2actm(maxref)
      integer l_r2act_offsetm(maxref),k_r2act_offsetm(maxref)
      integer l_e_offsetm(maxref),k_e_offsetm(maxref)
      integer size_em(maxref),size_r1m(maxref),size_r2m(maxref)
      integer d_em(maxref),d_r1m(maxref),d_r2m(maxref)
      integer d_corr,d_corr2
      double precision dbwcorr,dbwcorrf,dbwcorrf2
      double precision dbench1

c t1 localization
c      integer l_t1_local,k_t1_local
c IC logic 
      logical ic_cc,icrunl
c intermediate file
      integer d_c2
c
      double precision r1(maxref)      ! Residual t1
      double precision r2(maxref)      ! Residual t2
      double precision eone,etwo,enrep,energy
      integer dummy            ! Dummy argument for DIIS
      double precision residual,corr
      integer g5b,g6b
      double precision edelta
      double precision norm_r1,norm_r2
      double precision dsummary(5000,2)  !Print summary in the end
      integer k_movecs,l_movecs
      logical lcheckpoints
      integer icheckpoint
ckbn      logical lsave2e
ckbn      logical lread2e,lreadt
ckbn gf -2      integer lbwcorr 
      logical lbwcorr 
      logical lcas
      integer innodes
      integer inodesperref
      integer ifreenodes
      integer inofg, isum
      integer ipg,igg
      logical lstatus
      integer k
      integer isubsize
      integer d_f1mtmp(maxref)
      integer jref
      double precision dfirstterm
      integer istartmk
      double precision ddotr1,ddotr2
ckbn-5
      integer irefpt
      integer imkconv,iiref
      integer irefpth,jrefpth
      double precision pt3_1m(maxref),pt3_2m(maxref),pt3_3m(maxref)
      double precision pt_heff(maxref*maxref)
      integer l_t1_localpt,k_t1_localpt
      integer l_f1_localpt,k_f1_localpt
      integer l_f1_local,k_f1_local
      integer l_t1_local,k_t1_local
      integer pgidtobrdcst 
      integer k_ptheff,l_ptheff
      double precision hcnubycmu(maxref*maxref)
      double precision pt3_off2m(maxref*maxref)
      integer ipt,jpt
      integer itargetx
      double precision itcpu3  ! CPU sec counter
      double precision itwall3 ! WALL sec counter
      double precision itcpu3_1,itcpu3_2
      double precision itwall3_1,itwall3_2
      integer itmpkk,itmpkkk
      double precision dmint1,dmint2
      integer ioff
      double precision dmaxt1,dmaxt2
      double precision dresr1r2,dr1dotr1,dr2dotr2
      integer itype,idim1, idim2(4)
      integer tg_heff, ngheff, dg_heff

c      logical lusesamefock
c (t)
c PNNL: new variables
c === jaguar ===
      integer unitn,recforio
      logical diskint,dlogic
c === eaf ======
      integer eaf_han(maxref)   !handles of eaf files for v2(iref)
      logical eafl,eafa
      character*255 filemr(maxref)
c===============
cjb
C The number of bytes in double
      bytes = 8
#ifdef MRCC_METHODS
      tce_mrcc_energy = .false.
      nodezero=(ga_nodeid().eq.0)
      residual = 0.0d0
      lconverged = .false.
      ddotr1 = 0.0d0
      ddotr2 = 0.0d0
ckbn notes Just to remember IO
c    fortran : Fortran77 direct access,
c    eaf : Exclusive Access File library,
c    ga : Fully incore, Global Array virtual file,
c    sf : Shared File library,
c    replicated : Semi-replicated file on distributed file system with EAF library.
c    dra : Distributed file on distributed file system with DRA library.
c    ga_eaf : Semi-replicated file on distributed file system with EAF library. GA is used to speedup the file reconciliation. 
ckbn notes

ckbn remove below
      write(*,*) "I am here1"
      call ga_sync()
      tce_mrcc_energy = .true.

      if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,epsilon+enrep))
     + call errquit('tce_energy: RTDB problem',0,MA_ERR)

       if (nodezero) call util_flush(LuOut)

ckbn mrcc-r-1 CC      call tce_tidy(rtdb)
ckbn mrcc-r-1 CC      call util_print_pop
ckbn remove below end
      return
#else
      CALL ERRQUIT('MRCC Methods not compiled (tce_mrcc_energy.F) 
     + Compile  with "export MRCC_METHODS=y"',0,CAPMIS_ERR)
#endif

      end
