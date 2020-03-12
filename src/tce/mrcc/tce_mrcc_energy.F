      logical function tce_mrcc_energy(rtdb)
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
#include "tce_mrcc_diis.fh"
#endif

      integer rtdb
      logical se4t
      logical nodezero
      logical no_aposteriori
      character*255 modelname
      character*255 filename
      double precision ref     ! Ground state energy
      logical dft_energy
      external dft_energy
      logical scf
      external scf
ckbn      logical mcscf
ckbn      external mcscf
      integer g_ao1e(2)        ! GA handle for AO Fock matrices
      integer d_ao2e           ! SF handle for AO 2e integrals
ckbn  integer i,j,iref,nref,yref
      integer i,j,iref,yref
c     logical mrccdebug
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
c DEBUG LM2
      double precision lm2wall,lm2cpu
c
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
      integer d_corr
ckbn  integer d_corr2
ckbn      double precision dbwcorr,dbwcorrf,dbwcorrf2
ckbn      double precision dbench1

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
ckbn      double precision dsummary(5000,2)  !Print summary in the end
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
ckbn      double precision pt3_1m(maxref),pt3_2m(maxref),pt3_3m(maxref)
      double precision pt3_3m(maxref)
ckbn      double precision pt_heff(maxref*maxref)
      integer l_t1_localpt,k_t1_localpt
      integer l_f1_localpt,k_f1_localpt
      integer l_f1_local,k_f1_local
      integer l_t1_local,k_t1_local
      integer pgidtobrdcst 
      integer k_ptheff,l_ptheff
ckbn      double precision hcnubycmu(maxref*maxref)
ckbn      double precision pt3_off2m(maxref*maxref)
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
      integer sub_test
c =========== USS variables ===============
      logical lussd,lap_uss
      double precision dbwcorr,dbwcorrf,dbwcorrf2,dbwcorr1,dbwuss
      integer d_r1d(maxref), d_r2xy(maxref),d_corr1
      double precision dsmult
      integer signfact
c =========================================
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

ckbn introduce all checks in this block before going to actual calculation
c
c     ======================
c     Initialize parameters
c     ======================
c
       lsubterm = .false.
       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
     1 mrccdebug = .false.
      if (.not.rtdb_get(rtdb,'mrcc:checkpoints',mt_log,1,lcheckpoints))
     1   lcheckpoints = .false.
      if (.not.rtdb_get(rtdb,'mrcc:icheckpoint',mt_int,1,icheckpoint))
     1   icheckpoint = 0
      if (.not.rtdb_get(rtdb,'mrcc:save2e',mt_log,1,lsave2e))
     1   lsave2e = .false.
      if (.not.rtdb_get(rtdb,'mrcc:read2e',mt_log,1,lread2e))
     1   lread2e = .false.
      if (.not.rtdb_get(rtdb,'mrcc:readt',mt_log,1,lreadt))
     1   lreadt = .false.
      if (.not.rtdb_get(rtdb,'mrcc:forcedegen',mt_log,1,forcedegen))
     1   forcedegen = .false.
      if (.not.rtdb_get(rtdb,'mrcc:bwcorr',mt_log,1,lbwcorr))
     1   lbwcorr = .false.
      if (.not.rtdb_get(rtdb,'mrcc:cas',mt_log,1,lcas))
     2   lcas = .false.
ckbn sub -2
c      if (.not.rtdb_get(rtdb,'mrcc:subgroups',mt_log,1,lusesub))
c     +
         lusesub = .false.
      if (.not.rtdb_get(rtdb,'mrcc:subsize',mt_int,1,isubsize))
     2   isubsize = -1
      if (.not.rtdb_get(rtdb,'mrcc:istartmk',mt_int,1,istartmk))
     2   istartmk = -1
      if (.not.rtdb_get(rtdb,'mrcc:useeaf2el',mt_log,1,luseeaf2e))
     2   luseeaf2e = .false.
      if (.not.rtdb_get(rtdb,'mrcc:fullheff',mt_log,1,lfullheff))
     2   lfullheff = .false.

      if (.not.rtdb_get(rtdb,'mrcc:usesamefock_it',mt_log,1,
     +     lusesamefock_it))
     +   lusesamefock_it = .false.
       if(nodezero) write(*,*) "usesamefock_it",lusesamefock_it

       if (.not.rtdb_get(rtdb,'mrcc:ignorecomplex',mt_log,1,
     +   ignorecomplex)) ignorecomplex= .false.
       if(nodezero) write(*,*) "ignorecomplex",ignorecomplex
     
      if (.not.rtdb_get(rtdb,'mrcc:savet',mt_log,1,lsavet))
     1   lsavet = .false.
       if(nodezero) write(*,*) "savet",lsavet

ckbn-4
      if (.not.rtdb_get(rtdb,'mrcc:se4t',mt_log,1,se4t))
     2   se4t = .false.
      if(nodezero) write(*,*) "se4t", se4t

ckbn-4
      if (.not.rtdb_get(rtdb,'mrcc:no_aposteriori',
     +     mt_log,1,no_aposteriori)) no_aposteriori = .false.
      if(nodezero) write(*,*) "no_aposteriori", no_aposteriori


ckbn-13 placed it here and added to common block
       if (.not.rtdb_get(rtdb,'mrcc:usescffermiv',mt_log,1,lusescffv))
     + lusescffv = .false.
       if(nodezero) write(*,*) "Using fermivacuum by default" 
       lusescffv = .true.

ckbn roottooverlap
       if (.not.rtdb_get(rtdb,'mrcc:rootfromoverlap',mt_log,1,
     +  lrootfromoverlap))
     + lrootfromoverlap = .false.
       if(nodezero)  write(LuOut,*) 
     +  "Using rootfromoverlap ",lrootfromoverlap
       if(lrootfromoverlap) then
       if (.not.rtdb_get(rtdb,'mrcc:iroottooverlap',mt_int,1,
     +  iroottooverlap))
     +  call errquit('No root to overlap ',1,rtdb_err)
       if (.not.rtdb_get(rtdb,'mrcc:iroottooverlapiter',mt_int,1,
     +  iroottooverlapiter))
     +  call errquit('No iteration number to select root to overlap ',1,
     + rtdb_err)
       if(nodezero)
     +   write(LuOut,'(A,I5,A,I5)') 
     +   "Using rootfromoverlap with right eigen vector ",
     +    iroottooverlap," from iteration ", iroottooverlapiter
        lrootfromoverlap = .false. ! Turn it on from iroottooverlapiter
       endif


       if (.not.rtdb_get(rtdb,'mrcc:improvetiling',mt_log,1,limprovet))
     + limprovet = .false.
       if (.not.rtdb_get(rtdb,'mrcc:diistype',mt_int,1,idiis))
     +  idiis=0
       if (.not.rtdb_get(rtdb,'mrcc:zignore',mt_int,1,iignore))
     + iignore = 0
       if (.not.rtdb_get(rtdb,'mrcc:rootmuc',mt_int,1,nrootmuc))
     + nrootmuc = 0
       if (.not.rtdb_get(rtdb,'bwcc:targetroot',mt_int,1,iitarget))
     + iitarget = 0

       if (.not.rtdb_get(rtdb,'tce:tceiop',mt_int,1,recforio)) then
        rec_mem=1
       else
        rec_mem=recforio
       end if
ccccccccccccccccccccc USS correction settings ccccccccccccccccccccccccccccc
      if (.not.rtdb_get(rtdb,'mrcc:ussd',mt_log,1,lussd))
     1   lussd = .false.
      if (.not.rtdb_get(rtdb,'mrcc:uss_ap',mt_log,1,lap_uss))
     1   lap_uss = .false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c IC logic 
       icrunl=.false.
       if (.not.rtdb_get(rtdb,'tce:icrun',mt_log,1,icrunl)) then
        ic_cc=.false.
       else
        ic_cc=icrunl
       end if
c
       if(nodezero) then
        if(ic_cc) then
         write(LuOut,*)'IC version of CCSD is used'
        else
         write(LuOut,*)'standard version of CCSD is used'
        end if  
        call util_flush(LuOut)
       end if
c DISKINT
       dlogic=.false.
       if (.not.rtdb_get(rtdb,'tce:dlogic',mt_log,1,dlogic)) then
        diskint=.false.
       else
        diskint=dlogic
       end if
c EAF
       eafa=.false.
       if (.not.rtdb_get(rtdb,'tce:eafl',mt_log,1,eafa)) then
        eafl=.false.
       else
        eafl=eafa
       end if
c
       if(nodezero) then
        if(diskint.or.luseeaf2e) then
         write(LuOut,*)'DISK-MEMORY option enabled for integrals'
          if(eafl.or.luseeaf2e) then
           write(LuOut,*)'EAF files will be used'
          else 
           write(LuOut,*)'Standard direct access file will be used'
          end if
        else
         write(LuOut,*)'All integrlas stored in memory'
        end if
        call util_flush(LuOut)
       end if 
c *** debug ***
      if(nodezero) then
       write(LuOut,*)'BYTES BYTES',bytes
       call util_flush(6)
      end if

c     =========================
c     Ground-state HF/DFT first
c     =========================
c
      if (.not.rtdb_get(rtdb,'tce:reference',mt_int,1,reference)) then
        reference = 1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_energy: failed to write reference',0,
     2    RTDB_ERR)
      endif
      if (reference.eq.0) then
        if (.not.dft_energy(rtdb)) return
        if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,ref))
     1    call errquit('tce_energy: failed to get dft energy',0,
     2    RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.scf(rtdb)) return
        if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,ref))
     1    call errquit('tce_energy: failed to get scf energy',0,
     2    RTDB_ERR)
ckbn has to be tested      else if (reference.eq.2) then
ckbn has to be tested        if (.not.mcscf(rtdb)) return
ckbn has to be tested        if (.not.rtdb_get(rtdb,'mcscf:energy',mt_dbl,1,ref))
ckbn has to be tested     1   call errquit('tce_energy: failed to get mcscf energy',0,
ckbn has to be tested     2    RTDB_ERR)
      else
        call errquit('tce_energy: illegal reference',reference,
     1    UNKNOWN_ERR)
      endif
c
c
c
c        general model space even if RHF/ROHF orbitals 
c        are employed
c
c
         restricted=.false.
ckbn introduce all checks above this block before going to actual calculation end


ckbn analyze reference for mrcc calculation
c     ===================
c     Print utility start
c     ===================
      call util_print_push
      call util_print_rtdb_load(rtdb,'tce')

c     ==========
c     Initialize
c     ==========

ckbn MR-r      call tce_mrcc_init(rtdb)
ckbn replaced tce_mrcc_init with standard tce_init	!has to be tested
c      call tce_mrcc_init(rtdb)
      
      call tce_init(rtdb)

ckbn  check for maxorb
      if( ((nmo(1)-nfv(1)-nfc(1)) .gt. maxorb) .or. 
     +    ((nmo(2)-nfv(2)-nfc(2)) .gt. maxorb) ) then
       write(LuOut,'(A,I5,A,I5)') 'Maximum number of MOs ',
     +  maxorb,' is lower than current number of MOs',
     +  max((nmo(1)-nfv(1)-nfc(1)),(nmo(2)-nfv(2)-nfc(2)))
       call util_flush(LuOut)
        call errquit
     + ('tce_mrcc_energy: Above maximum number of MOs',0,RTDB_ERR)
      endif

ckbn Exit if diis is greater than maxdiis
      if(diis .gt. maxdiism) then
       write(LuOut,'(A,I5,A,I5)') 'Current diis value', diis,
     + ' is larger than maximum', maxdiism
       call util_flush(LuOut)
        call errquit
     + ('tce_mrcc_energy: Above maximum diis value',0,RTDB_ERR)
      endif

c      write(*,*) "ioalg", ioalg
c      if(ioalg .ne. 2 ) call errquit
c     + ('tce_mrcc_energy: TCE MRCC needs ioalgorithm',0,INPUT_ERR)

ckbn MR-r      if(lcas)call tce_mrcc_create_cas(rtdb)
c      if(lcas)call tce_mrcc_create_cas(rtdb)
      if(lcas)call tce_mrcc_create_cas1(rtdb)
      if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)
ckbn MR-r      call tce_mrcc_readref(rtdb) 
      call tce_mrcc_readref(rtdb) 

      if (.not.rtdb_get(rtdb,'bwcc:targetroot',mt_int,1,itargetx)) then
c       write(*,*) "I am here1","itargetx",itargetx
       if(.not.rtdb_get(rtdb,'mrcc:roottooverlap',mt_int,1,
     +  iroottooverlap)) call errquit
     +  ('tce_mrcc_read: rtdb failed reading root to overlap',
     +  0,RTDB_ERR)
      else
        mkrootold = itargetx
        mkroot    = itargetx
c        write(*,*) "I am here1 mkroot",itargetx,mkrootold,mkroot
      endif
ckbn analyze reference for mrcc calculation end



      if (nodezero) call util_flush(LuOut)

ckbn Test for MA and GA 
ckbn @todo : Here test for memory bottle-neck and exit if it is not met.
ckbn			Use ga_memory_avail()
ckbn 			Add a brdcst test whether everybody else has same value
      call sf_test()

c     ================
c     Arrays for tile
c     ================
      cpu = - util_cpusec()
      wall = - util_wallsec()
      call tce_mrcc_tile(rtdb)                     ! more or less ok
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()

      if(nodezero) then
       write(LuOut,"(/,'MRCC tiling completed in ',2f15.1)") cpu, wall
       call util_flush(LuOut)
      endif

      call tce_mrcc_fillindexarray(rtdb)           ! ok

c     ==============
c     1-e integrals 
c     ==============
      do iref=1,nref
       g_movecs(1) = g_movecsm(iref,1)
       g_movecs(2) = g_movecsm(iref,2)
       k_sym = k_symm(iref)
       k_offset = k_offsetm(iref)
       k_range = k_rangem(iref)        
       k_spin = k_spinm(iref)
       k_movecs_sorted = k_movecs_sortedm(iref)
       noa = nblcks(1,iref)
       nob = nblcks(2,iref)      
       nva = nblcks(3,iref)
       nvb = nblcks(4,iref)
       noab = noa+nob
       nvab = nva+nvb

       cpu  = - util_cpusec()
       wall = - util_wallsec()
       call tce_ao1e_fock2e(rtdb,g_ao1e,lzorafromdft) ! F in AO first, g_ao1e destroyed in the mo1e routine
       cpu  = cpu  + util_cpusec()
       wall = wall + util_wallsec()
       if(nodezero) then
        write(LuOut,'(A,2f15.5)')"tce_ao1e_fock2e", cpu, wall
        call util_flush(LuOut)
       endif

       write(namechunk,"(I3.3)")iref
       call tce_filename('f1m'//namechunk,filename)
       call tce_mo1e_offset(l_f1_offsetm(iref),
     1 k_f1_offsetm(iref),size_1em(iref))


      if((nodezero .and. (.not.lusesub)).or.
     +((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub))
     +write(*,'(A,I5,A,I20)')"F: ",iref," in bytes = ",size_1em(iref)*8

         call createfile(filename,d_f1m(iref),size_1em(iref))
c PNNL: add read option
       if(lread2e) then
          call tce_filenameindexed(iref,'mr_f',filename)
          unitn=76
          call read_tensor(filename,d_f1m(iref),size_1em(iref),unitn)
c         close(unitn)
         call reconcilefile(d_f1m(iref),size_1em(iref))
       else
         cpu  = - util_cpusec()
         wall = - util_wallsec()
         call tce_mo1e(g_ao1e,d_f1m(iref),k_f1_offsetm(iref))
         cpu  = cpu  + util_cpusec()
         wall = wall + util_wallsec()
         if(nodezero) then
           write(LuOut,'(A,2f15.5)')"tce_mo1e", cpu, wall
           call util_flush(LuOut)
         endif
        if(lsave2e) then
           call tce_filenameindexed(iref,'mr_f',filename)
           unitn=76
           call write_tensor(filename,d_f1m(iref),size_1em(iref),unitn)
c          close(unitn)
           call reconcilefile(d_f1m(iref),size_1em(iref))
         endif
       endif
c PNNL: add write option
c PNNL: kill d_f1m(iref)
c         call gatoeaf(d_f1m(iref))
c         call ga_print(d_f1m(iref))

         cpu = - util_cpusec()
         wall = - util_wallsec()
         call mrcc_uhf_energy(g_movecs, eone, etwo, enrep, energy,
     +    lzorafromdft)
         duhfens(iref) = eone+etwo
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()
         if(nodezero) then
           write(LuOut,'(A,2f15.5)')"mrcc_uhf_energy", cpu, wall
           call util_flush(LuOut)
         endif

ckbn zoroa
c        if(nodezero) then
c         call ga_print(d_f1m(iref))
cc         call ga_print(size_1em(iref))
c         write(*,*) "UHF Starts"
c         call ga_print(d_f1m(1))
c         call ga_print(size_1em(1))
c         write(*,*)"E1",energy
c        endif
c        call ga_sync()

c DISKINT
c         if(diskint) then
c          call tce_filenameindexed(iref,'mrfn_r',filename)
c          unitn=76+(iref-1)*2+1
c          call write_tensor(filename,d_f1m(iref),size_1em(iref),unitn)
c          call deletefile(d_f1m(iref))
c         end if
c
      if (.not. int2e_file_close(.false.))
     $     call errquit('tce_mrcc_energy: closing aoints?', 0, INT_ERR)
      enddo

c
c     ==============
c        Subgroups
c     ==============
c
      if (.not.rtdb_get(rtdb,'mrcc:subgroupsize',mt_log,1,lusesub))
     2   lusesub = .false.

          if(nodezero.and.lusesub)
     1    write(LuOut,"(/,'Subgroups will be used',/,
     1                    '======================',/)")
        


c         if(lusesub) then 
cckbn
cckbn          write(LuOut,*) "GA_CLUSTER_NNODES",GA_CLUSTER_NNODES()
cckbn          call ga_sync()
cckbn          call util_flush(LuOut)
cckbn          write(LuOut,*) "GA_CLUSTER_NODEID",GA_CLUSTER_NODEID()
cckbn          call ga_sync()
cckbn          call util_flush(LuOut)
c          sub_test = GA_CLUSTER_NODEID()
cckbn          write(LuOut,*) "GA_CLUSTER_NPROCS",GA_CLUSTER_NPROCS(sub_test)
cckbn          call ga_sync()
cckbn          call util_flush(LuOut)
c          isubsize =  GA_CLUSTER_NNODES()/nref
c          if(nodezero) write(*,*) "isubsize ",isubsize
c          if(isubsize .lt. 1) 
c     +    call errquit("tce_mrcc_energy: isubsize",0,MA_ERR)
cckbn          write(LuOut,*) "Subgroup size 1 ", isubsize
cckbn          call ga_sync()
cckbn          call util_flush(LuOut)
c
c          isubsize =  GA_CLUSTER_NPROCS(sub_test) * isubsize
c          if(nodezero) write(LuOut,*) "Subgroup size", isubsize
c          if(isubsize .lt. 
c     +     ((GA_CLUSTER_NNODES()*GA_CLUSTER_NPROCS(sub_test))/nref))then
c             write(LuOut,*) "Warning: Number of free processors",
c     +        ((GA_CLUSTER_NNODES()*GA_CLUSTER_NPROCS(sub_test))
c     +        - (isubsize*nref))
c           call util_flush(LuOut)
c          endif
c        
c        endif       
c




        innodes = ga_nnodes()
        if(lusesub) then 
c        if(nodezero) write(LuOut,"('NNodes: ',I8)") innodes
         if(nodezero) 
     +    write(LuOut,"('Total Number of cores: ',I8)") innodes
         call util_flush(LuOut)
        endif

      if (.not.ma_push_get(mt_int,innodes*2,"all innodes",
     2  l_innodes,k_innodes))
     3  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,innodes,"sg sizes",
     2  l_sgsizes,k_sgsizes))
     3  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,innodes,"sghandles",
     2  l_sghandles,k_sghandles))
     3  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,nref,"nref affiliation",
     2  l_refafi,k_refafi))
     3  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)

        do i=1,innodes
          int_mb(k_innodes+i-1) = i-1
          int_mb(k_sgsizes+i-1) = 0
          int_mb(k_sghandles+i-1) = 0
        enddo

        if(isubsize.le.0) then

        inodesperref = (innodes/nref)
        if(inodesperref.eq.0) inodesperref = 1

        ifreenodes = mod(innodes,nref) 
        if(innodes.le.nref) ifreenodes = 0

        inofg = min(innodes,nref)

        else

        inodesperref = isubsize
        inofg = innodes/isubsize 
c        inofg = (innodes-
c     +        ((GA_CLUSTER_NNODES()*GA_CLUSTER_NPROCS(sub_test))
c     +        - (isubsize*nref)))/isubsize 
        ifreenodes = innodes - (isubsize*inofg)

        endif


        if(lusesub) then 
        if(nodezero) then
          write(LuOut,"(/,'Nodes per ref: ',I8)")inodesperref
          write(LuOut,"('Free nodes left: ',I8)")ifreenodes
          write(LuOut,"('Total number of non-overlapping groups: '
     1 ,I8,/)")inofg
        endif
        endif

      if (.not.ma_push_get(mt_int,inofg,"g offsets",
     2  l_goffset,k_goffset))
     3  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)

        do i=1,inofg
            int_mb(k_sgsizes+i-1) = inodesperref
        enddo
        do i=1,ifreenodes
            int_mb(k_sgsizes+i-1) = int_mb(k_sgsizes+i-1)+1
        enddo

        isum = 0
        do i=1,inofg
           int_mb(k_goffset+i-1) = isum
           isum = isum + int_mb(k_sgsizes+i-1)
        enddo

        if(lusesub) then
        do i=1,inofg
        int_mb(k_sghandles+i-1) = ga_pgroup_create(int_mb(k_innodes+
     1 int_mb(k_goffset+i-1)),int_mb(k_sgsizes+i-1))
        enddo
        endif

        do i=1,nref
          int_mb(k_refafi+i-1) = int_mb(k_sghandles+mod(i-1,inofg))
        enddo

        k = 0
        do j=1,inofg
        do i=1,int_mb(k_sgsizes+j-1)
           k = k + 1
          int_mb(k_innodes+innodes+k-1) = int_mb(k_sghandles+j-1)
        enddo
        enddo

        if(lusesub) then 
         if(nodezero) then

          write(LuOut,"(/,'Processor distribution',/,
     1                    '======================',/)")

          do i=1,inofg
          write(LuOut,"('Subgroup no. ',I6,' size: ',I6)")
     1 int_mb(k_sghandles+i-1),
     1 ga_pgroup_nnodes(int_mb(k_sghandles+i-1))
          enddo

          write(LuOut,"(/,'Reference static distribution',/,
     1                    '=============================',/)")
          do i=1,nref
          write(LuOut,"('Reference no. ',I6,' will run on: ',I6)")i,
     1 int_mb(k_refafi+i-1)
          enddo
           write(LuOut,*)' '
         endif
        endif

         mypgid = int_mb(k_innodes+ga_nodeid()+ga_nnodes())
c         write(LuOut,'(A,I5,I5,I5,I5)') "mypgid2 ",mypgid,ga_nodeid(),
c     +       ga_pgroup_get_default(),ga_nnodes()
         call util_flush(LuOut) 

c -------
c      call create_mr_evl_sorted(nref,k_f1_offsetm)
      call create_mr_evl_sorted(k_f1_offsetm)

      
#ifdef MRCC_LOCAL_FOCK
c      call ga_print(d_f1m(1))
      if(lusesub) then
       do iref = 1, nref 
c         k_sym = k_symm(iref)
c         k_offset = k_offsetm(iref)
c         k_range = k_rangem(iref)
c         k_spin = k_spinm(iref)
c         k_movecs_sorted = k_movecs_sortedm(iref)
c         k_evl_sorted = k_evl_sortedm(iref)
c         k_active = k_active_tmpm(iref)
c         noa = nblcks(1,iref)
c         nob = nblcks(2,iref)
c         nva = nblcks(3,iref)
c         nvb = nblcks(4,iref)
c         noab = noa+nob
c         nvab = nva+nvb

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1   +innodes+ga_nodeid())).or.(.not.lusesub)) then
         if(.not.MA_PUSH_GET(mt_dbl,size_1em(iref),'f1_local',
     +    l_f1_local,k_f1_local)) call errquit('f1_local',1,MA_ERR)
         call ma_zero(dbl_mb(k_f1_local),size_1em(iref))
         if( ga_pgroup_nodeid(mypgid) .eq. 0) then
          call ga_get(d_f1m(iref),1,size_1em(iref),1,1,
     +                dbl_mb(k_f1_local),1) 
         endif
         call ga_pgroup_brdcst(mypgid,65549,dbl_mb(k_f1_local),
     +    MA_sizeof(MT_DBL,1,MT_BYTE)*size_1em(iref)*1,0)
        endif
       enddo
      endif
#endif
      


ckbn @todo : this is our present memory bottle-neck
ckbn @todo : This part has to be rewritten
c     ==============
c     2-e integrals 
c     ==============
      do iref=1,nref
c       if( ga_pgroup_nodeid(mypgid).eq.0) 
c     +   write(*,'(A,I3,A,I3,A,I3)') "F",int_mb(k_refafi+iref-1),
c     +          "S",int_mb(k_innodes+innodes+ga_nodeid()),"iref",iref
       if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then
        g_movecs(1) = g_movecsm(iref,1)
        g_movecs(2) = g_movecsm(iref,2)
        k_sym = k_symm(iref)
        k_offset = k_offsetm(iref)
        k_range = k_rangem(iref)
        k_spin = k_spinm(iref)
        k_movecs_sorted = k_movecs_sortedm(iref)
        k_active = k_active_tmpm(iref)
        noa = nblcks(1,iref)
        nob = nblcks(2,iref)
        nva = nblcks(3,iref)
        nvb = nblcks(4,iref)
        noab = noa+nob
        nvab = nva+nvb
c  2-e integrals first half
        if(.not.lread2e) then
         cpu = - util_cpusec()
         wall = - util_wallsec()
         if(fast2e.eq.1) then
          call tce_mrcc_ao2e(rtdb,d_ao2e,iref)
         elseif (fast2e.eq.20) then
          if(nodezero) then
           write(LuOut,*)'fast2e 20 4-index transformation is used'
           call util_flush(6)
          end if
         else
          write(LuOut,9123)fast2e
          call errquit('tce_energy: invalid 2emet: ',fast2e,CALC_ERR)
         end if
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()
c        if((nodezero) .or. (ga_pgroup_nodeid(mypgid).eq.0)) then
         if((nodezero .and. (.not.lusesub)).or.
     +       ((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub)) then
         write(LuOut,'(A,I4,A,2f15.2)')"Ref.",iref," Half 2-e",cpu, wall
         call util_flush(LuOut)
         endif
        endif !read2e
c  2-e integrals second half
        write(namechunk,"(I3.3)")iref
        cpu = - util_cpusec()
        wall = - util_wallsec()
        call tce_mo2e_offset(l_v2_offsetm(iref),k_v2_offsetm(iref),
     1                       size_2em(iref))
        call tce_filename('v2'//namechunk,filename)
c       if((nodezero) .or. (ga_pgroup_nodeid(mypgid).eq.0)) 
        if((nodezero .and. (.not.lusesub)).or.
     +       ((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub))
     +   write(LuOut,'(A,A,A,I18)')
     +    "V 2-e ",filename(1:20)," in bytes= ",size_2em(iref)*8
        if((fast2e.eq.1).or.(fast2e.eq.20)) then  ! added 20
         if(lusesub) then
          call pgcreatefile(filename,d_v2m(iref),size_2em(iref),
     1    int_mb(k_innodes+innodes+ga_nodeid()))
         else
          call createfile(filename,d_v2m(iref),size_2em(iref))
         endif
         if(lread2e) then
          call tce_filenameindexed(iref,'mrvn_r',filename)
          unitn=77
          call read_tensor(filename,d_v2m(iref),size_2em(iref),unitn)
c         close(unitn)
c          if(nodezero) write(*,*) "I am here read mrvn_r"
          if(lusesub) then
           call ga_pgroup_sync(mypgid)
          else
           call reconcilefile(d_v2m(iref),size_2em(iref))
          endif
         else
c
c
c
c
c
         if(fast2e.eq.20) then
c
         if(.not.lusesub) 
     1    call errquit('2emet 20 can only be used with subgroups',
     1    0,MA_ERR)
c
        if (.not.rtdb_get(rtdb,'tce:maxs',mt_int,1,imaxsize)) then
         imaxsize = 30
         if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,imaxsize))
     1   call errquit('tce_init: failed to write maxs',0,RTDB_ERR)
        endif
          call atpart_nalength_mrcc
          lm2cpu = - util_cpusec()
          lm2wall = - util_wallsec()
          call new_ga4ind_N5uhf_mrcc_lm2(rtdb,d_v2m(iref),
     1                                k_v2_offsetm(iref),
     1                                size_2em(iref),iref)
          lm2cpu = lm2cpu + util_cpusec()
          lm2wall = lm2wall + util_wallsec()
          if(nodezero) then
           write(LuOut,"('LM2 2-e ints ',2f15.1)") lm2cpu,lm2wall
          end if
         else
          call tce_mrcc_mo2e(rtdb,d_ao2e,d_v2m(iref),k_v2_offsetm(iref),
     1     iref,d_v2m(1),k_v2_offsetm(1))
         end if
c
c
c
c
          if((nodezero .and. (.not.lusesub)).or.
     +       ((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub))
     +    write(LuOut,'(A,I4,A)')"Ref.",iref," 2-e transform. completed"
          if(lusesub) then
           call ga_pgroup_sync(mypgid)
          else
           call reconcilefile(d_v2m(iref),size_2em(iref))
          endif
          if(lsave2e) then
           call tce_filenameindexed(iref,'mrvn_r',filename)
           unitn=77
           call write_tensor(filename,d_v2m(iref),size_2em(iref),unitn)
c          close(unitn)
           if(lusesub) then
            call ga_pgroup_sync(mypgid)
           else
            call reconcilefile(d_v2m(iref),size_2em(iref))
           endif
          endif
         endif
         if(luseeaf2e) then 
          call ga_sync()
          call eaf_w_mrcc(eaf_han(iref),d_v2m(iref),size_2em(iref))
          call ga_sync()
          call deletefile(d_v2m(iref))
         endif
        else
         call errquit('tce_energy: invalid 2emet: ',fast2e,CALC_ERR)
        endif
c        if((nodezero) .or. (ga_pgroup_nodeid(mypgid).eq.0)) 
c     +   write(LuOut,*)'v2-size iref: ',size_2em(iref),iref
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if(((nodezero) .or. (ga_pgroup_nodeid(mypgid).eq.0)) 
     +   .and. mrccdebug) then
          write(LuOut,"('MO 2-e ints ',2f15.1)") cpu, wall
          call util_flush(LuOut)
        endif
       endif ! subgroups
      enddo ! nref
c
c     =========================
c       Allocations of T files
c     =========================
c
      if(nodezero) then
         write(LuOut,"(/,'Integral replication completed.',/,
     1 'Proceeding the allocation of the memory for intermediates.',/)")
      endif

      needt1 = .false.
      needt2 = .false.
      needt3 = .false.
      needt3a= .false.
      needr3act = .false.
      needr4act = .false.

      if (model.eq.'bwccsd') then
         needt1 = .true.
         needt2 = .true.
      elseif (model.eq.'bwccsdt') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
      elseif (model.eq.'mkccsd') then
         needt1 = .true.
         needt2 = .true.
      elseif (model.eq.'mkccsdt') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
      elseif (model.eq.'succsd') then
         needt1 = .true.
         needt2 = .true.
      endif

      do iref=1,nref

c        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
c     1 +innodes+ga_nodeid())).or.(.not.lusesub).or.
c     2 (model.eq.'mkccsd')) then
        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then


         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)

c         k_evl_sorted = k_evl_sortedm(iref)
c         k_active = k_active_tmpm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

      if(needt2) then

        write(namechunk,"(I3.3)")iref
        call tce_filename('t2'//namechunk,filename)
        call tce_t2_offset(l_t2_offsetm(iref),k_t2_offsetm(iref),
     1 size_t2m(iref))

      if((nodezero .and. (.not.lusesub)).or.
     +((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub)) then
       write(LuOut,'(A,I5,A,I20)')
     +  "T2: ",iref," in bytes = ",size_t2m(iref)*8
       call util_flush(LuOut)
      endif
        if(lusesub.and..not.(model.eq.'mkccsd')) then
        call pgcreatefile(filename,d_t2m(iref),size_t2m(iref),
     1 int_mb(k_innodes+innodes+ga_nodeid()))
        else
          call createfile(filename,d_t2m(iref),size_t2m(iref))
          call reconcilefile(d_t2m(iref),size_t2m(iref))
        endif

        if(lreadt) then
           call tce_filenameindexed(iref,'mrt2',filename)
           unitn=79
           call read_tensor(filename,d_t2m(iref),size_t2m(iref),unitn)
c           close(unitn)
           call reconcilefile(d_t2m(iref),size_t2m(iref))
c           if(nodezero) write(*,*) "I am here read mrt2"
        endif

      endif



      if(needt1) then

        write(namechunk,"(I3.3)")iref
        call tce_filename('t1'//namechunk,filename)
        call tce_t1_offset(l_t1_offsetm(iref),k_t1_offsetm(iref),
     1 size_t1m(iref))

      if((nodezero .and. (.not.lusesub)).or.
     +((ga_pgroup_nodeid(mypgid).eq.0) .and. lusesub)) then
       write(LuOut,'(A,I5,A,I20)')
     +  "T1: ",iref," in bytes = ",size_t1m(iref)*8
       call util_flush(LuOut)
      endif
        if(lusesub.and..not.(model.eq.'mkccsd')) then
        call pgcreatefile(filename,d_t1m(iref),size_t1m(iref),
     1 int_mb(k_innodes+innodes+ga_nodeid()))
        else
        call createfile(filename,d_t1m(iref),size_t1m(iref))
        call reconcilefile(d_t1m(iref),size_t1m(iref))
        endif

        if(lreadt) then
           call tce_filenameindexed(iref,'mrt1',filename)
           unitn=78
           call read_tensor(filename,d_t1m(iref),size_t1m(iref),unitn)
c           close(unitn)
           call reconcilefile(d_t1m(iref),size_t1m(iref))
c           if(nodezero) write(*,*) "I am here read mrt1"
        endif

      endif




      if(needt3) then

        write(namechunk,"(I3.3)")iref
        call tce_filename('t3'//namechunk,filename)
        call tce_t3_offset(l_t3_offsetm(iref),k_t3_offsetm(iref),
     1 size_t3m(iref))
       if(nodezero) write(*,'(A,I5,I20)')"T3sz for",iref,size_t3m(iref)
        call createfile(filename,d_t3m(iref),size_t3m(iref))
c        call gatoeaf(d_t3m(iref))
c        call ga_zero(d_t3m(iref))
        call reconcilefile(d_t3m(iref),size_t3m(iref))

      endif
 
      if(lfullheff) then

c        write(namechunk,"(I3.3)")iref

c        call tce_filename('r3act'//namechunk,filename)
      call tce_t3ac_offset(l_r3a_offsetm(iref),k_r3a_offsetm(iref),
     1 size_r3am(iref),iref)
c        call createfile(filename,d_r3am(iref),size_r3am(iref))

c        call tce_filename('r4act'//namechunk,filename)
      call tce_r4ac_offset(l_r4a_offsetm(iref),k_r4a_offsetm(iref),
     1 size_r4am(iref),iref)
c        call createfile(filename,d_r4am(iref),size_r4am(iref))

      endif

      endif !sub

      enddo



c
c     =======================
c     Create residual offsets
c     =======================
c
      do iref=1,nref

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_active = k_active_tmpm(iref)
         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb
  
      call tce_e_offset(l_e_offsetm(iref),k_e_offsetm(iref),
     1 size_em(iref))
      if (needt1) call tce_t1_offset(l_r1_offsetm(iref),
     1 k_r1_offsetm(iref),size_r1m(iref))
      if (needt2) then
         call tce_t2_offset(l_r2_offsetm(iref),
     1 k_r2_offsetm(iref),size_r2m(iref))
      endif
      if (needt3) call tce_t3_offset(l_r3_offsetm(iref),
     1 k_r3_offsetm(iref),size_r3m(iref))
c      if (needr3act) call tce_r3_offset(l_r3act_offsetm(iref),
c     1 k_r3act_offsetm(iref),size_r3actm(iref))
c      if (needr4act) call tce_r4_offset(l_r4act_offsetm(iref),
c     1 k_r4act_offsetm(iref),size_r4actm(iref))

c       if(nodezero)
c     1 write(LuOut,"(/,'Tasks for ref: ',I5,2I8)")iref,
c     2 int_mb(k_r1_offsetm(iref)), int_mb(k_r2_offsetm(iref))
 
      enddo
c
c     ======================
c        k_evl_sorted
c     ======================
c
c      call create_mr_evl_sorted(nref,k_f1_offsetm)
c
c     ======================
c        Heff matrix alloc
c     ======================
c
      if (.not.ma_push_get(mt_dbl,nref*nref,'Heff',l_heff,k_heff))
     1 call errquit('tce_mrcc_energy: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nref*nref,'sqc',l_sqc,k_sqc))
     1 call errquit('tce_sqc: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nref*nref,'sqc',l_sqcl,k_sqcl))
     1 call errquit('tce_sqcl: MA problem',0,MA_ERR)

      if(lusesub) then
        call tce_filename('g_heff',filename)
        call createfile(filename,g_heff,nref*nref)
        if(lbwcorr) then
        call tce_filename('corr',filename)
        call createfile(filename,d_corr,nref)
        endif
      endif

      if(nodezero) then
         write(LuOut,"('Done.',/)")
      endif
c
c     ===================
c       Refs symmetries
c     ===================
c
c      call tce_mrcc_refs_sym(nref)
      call tce_mrcc_refs_sym()
c
c     ================================
c       Set degenerated orb energies
c     ================================
c
c      if(forcedegen)call tce_mrcc_force_orben(nref)
      if(forcedegen)call tce_mrcc_force_orben()



c
c     ==============
c        Methods
c     ==============
c
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
      if(lusesub) then
       do iref = 1, nref 
        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1   +innodes+ga_nodeid())).or.(.not.lusesub)) then
         if(.not.MA_PUSH_GET(mt_dbl,size_t1m(iref),'t1_local',
     +    l_t1_local,k_t1_local)) call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1m(iref))
         if( ga_pgroup_nodeid(mypgid) .eq. 0) then
          call ga_get(d_t1m(iref),1,size_t1m(iref),1,1,
     +                dbl_mb(k_t1_local),1) 
         endif
         call ga_pgroup_brdcst(mypgid,65568,dbl_mb(k_t1_local),
     +    MA_sizeof(MT_DBL,1,MT_BYTE)*size_t1m(iref)*1,0)
        endif
       enddo
      endif
#endif





      if(lusesub) call ga_sync()

c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd begin
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd begin
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd begin


      if (model.eq.'bwccsd') then

       if(nodezero) write (LuOut,"('MR BWCCSD, version 1.0')")

       lsubterm = .true.

       cpu = - util_cpusec()
       wall = - util_wallsec()

       call tce_mrcc_diis_init()

       edelta =0.0d0

       do iter=1,maxiter

        itcpu2 = - util_cpusec()
        itwall2 = - util_wallsec()

c        if(nodezero) write(*,*) "I am here1"

#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
        if(lusesub) then
         do iref = 1, nref 
c         k_sym = k_symm(iref)
c         k_offset = k_offsetm(iref)
c         k_range = k_rangem(iref)
c         k_spin = k_spinm(iref)
c         k_movecs_sorted = k_movecs_sortedm(iref)
c         k_evl_sorted = k_evl_sortedm(iref)
c         k_active = k_active_tmpm(iref)
c         noa = nblcks(1,iref)
c         nob = nblcks(2,iref)
c         nva = nblcks(3,iref)
c         nvb = nblcks(4,iref)
c         noab = noa+nob
c         nvab = nva+nvb

          if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1     +innodes+ga_nodeid())).or.(.not.lusesub)) then
c           if(.not.MA_PUSH_GET(mt_dbl,size_t1m(iref),'t1_local',
c     +     l_t1_local,k_t1_local)) call errquit('t1_local',1,MA_ERR)
           call ma_zero(dbl_mb(k_t1_local),size_t1m(iref))
           if( ga_pgroup_nodeid(mypgid) .eq. 0) then
            call ga_get(d_t1m(iref),1,size_t1m(iref),1,1,
     +                dbl_mb(k_t1_local),1) 
           endif
           call ga_pgroup_brdcst(mypgid,65568,dbl_mb(k_t1_local),
     +         MA_sizeof(MT_DBL,1,MT_BYTE)*size_t1m(iref)*1,0)
          endif
         enddo
        endif
#endif



c        if(.not. lusesub) then
cckbn rootto overlap intialization
cc      if(nodezero) then
cc      write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+0)
cc      write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+1)
cc      write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+0)
cc      write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+1)
cc      endif
c         if(iter.eq.3) then
c          if(iroottooverlap.gt.0) then
c           if(nodezero) write(*,*) "iroottooverlap ", iroottooverlap
c           do i=1,nref
c            bwcoefwanted(i) = dbl_mb(k_sqc+(i-1)*nref+iroottooverlap-1) 
c           enddo
c           lrootfromoverlap = .true.
c           do i=1,nref
c            if(nodezero) write(*,'(A,F17.10)') 
c     +                   "bwcoefwanted",bwcoefwanted(i)
c           enddo
c          endif
c          mkrootold = mkroot
c          endif

ckbn root to overlap intialization
        if(.not. lusesub) then !ckbn can be removed
c         if(nodezero) then
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+0)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+1)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+0)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+1)
c         endif
         if(iter.eq.iroottooverlapiter) then
          if(iroottooverlap.gt.0) then
c           if(nodezero) write(*,*) "iroottooverlap ", iroottooverlap
           if(nodezero) write(LuOut,'(A,I5,A,I5)')
     +      "Overlap selected from iteration ", iroottooverlapiter,
     +      " with respect to root ", iroottooverlap
           do i=1,nref
            bwcoefwanted(i) = dbl_mb(k_sqc+(i-1)*nref+iroottooverlap-1)
           enddo
           lrootfromoverlap = .true.
           do i=1,nref
            if(nodezero) write(*,'(A,I5,F17.10)')
     +       "Overlapping right eigen vector ",i,bwcoefwanted(i)
           enddo
          else
           call errquit('root to overlap not natural number',0,RTDB_ERR)
          endif
          mkrootold = mkroot
         endif


ckbn check if roots swapped
         if(iter.gt.4) then
          if(nodezero) then
c           write(LuOut,'(A,I3,A,I3)')
c     +          "mkroot",mkroot,"mkrootold",mkrootold
           if(mkroot.ne.mkrootold) then
            write(LuOut,'(A,I5,A,I5)') 
     +      "Warning: Eigenvector swap occured from",
     +      mkrootold," to ",mkroot
            call util_flush(LuOut)
           endif
          endif
          mkrootold = mkroot
         endif
        endif !.not. lusesub

c        if(nodezero) write(*,*) "I am here3"

        do iref=1,nref        

         if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1     +innodes+ga_nodeid())).or.(.not.lusesub)) then

          if(lusesub)
     1      ipg = int_mb(k_innodes+ga_nnodes()+ga_nodeid())

          g_movecs(1) = g_movecsm(iref,1)
          g_movecs(2) = g_movecsm(iref,2)

          k_sym = k_symm(iref)
          k_offset = k_offsetm(iref)
          k_range = k_rangem(iref)
          k_spin = k_spinm(iref)
          k_movecs_sorted = k_movecs_sortedm(iref)
          k_evl_sorted = k_evl_sortedm(iref)
          k_active = k_active_tmpm(iref)

          noa = nblcks(1,iref)
          nob = nblcks(2,iref)
          nva = nblcks(3,iref)
          nvb = nblcks(4,iref)

          noab = noa+nob
          nvab = nva+nvb

          write(namechunk,"(I3.3)")iref

c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step5 iter ref',iter,iref
c       call util_flush(6)
c      endif
c kk ---
          if(luseeaf2e) then
           call tce_filename('v2',filename)
           call createfile(filename,d_v2m(iref),size_2em(iref))
           call eaf_r_mrcc(eaf_han(iref),d_v2m(iref),size_2em(iref))
           call ga_sync()
          end if
c
          d_v2orb=d_v2m(iref)
c
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step6 iter ref',iter,iref
c       call util_flush(6)
c      endif
c kk ---
          call tce_filename('e'//namechunk,filename)

          if(lusesub) then
          call pgcreatefile(filename,d_em(iref),size_em(iref),mypgid)
#ifdef MRCC_LOCAL_FOCK
ckbn @todo : localize fock
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
          call cxsd_e(k_f1_local,d_em(iref),k_t1_local,
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_e_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref),ipg)
#else
          call cxsd_e(k_f1_local,d_em(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_e_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref),ipg)
#endif
#else
          call cxsd_e(d_f1m(iref),d_em(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_e_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref),ipg)
#endif

          else
          call createfile(filename,d_em(iref),size_em(iref))
          call mr_ccsd_e(d_f1m(iref),d_em(iref),d_t1m(iref),
     1               d_t2m(iref),d_v2m(iref),
     1               k_f1_offsetm(iref),k_e_offsetm(iref),
     2               k_t1_offsetm(iref),k_t2_offsetm(iref),
     2               k_v2_offsetm(iref))
        endif
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step6 iter ref',iter,iref
c       write(6,*)'ga mem used6',ga_inquire_memory()
c       write(6,*)'ga mem avail6',ga_memory_avail()
c       call util_flush(6)
c      endif
c kk ---

          call tce_filename('r1'//namechunk,filename)
       if(lusesub) then
          call pgcreatefile(filename,d_r1m(iref),size_r1m(iref),mypgid)
#ifdef MRCC_LOCAL_FOCK
ckbn @did : localize fock
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
          call cxsd_t1(k_f1_local,d_r1m(iref),k_t1_local,
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),mypgid)
#else
          call cxsd_t1(k_f1_local,d_r1m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),mypgid)
#endif
#else
          call cxsd_t1(d_f1m(iref),d_r1m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),mypgid)
#endif
       else
          call createfile(filename,d_r1m(iref),size_r1m(iref))
          call mr_ccsd_t1(d_f1m(iref),d_r1m(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref))
       endif
ccx          call reconcilefile(d_r1m(iref),size_r1m(iref))
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step7 iter ref',iter,iref
c       call util_flush(6)
c      endif
c kk ---



          call tce_filename('r2'//namechunk,filename)
       if(lusesub) then
          call pgcreatefile(filename,d_r2m(iref),size_r2m(iref),mypgid)
#ifdef MRCC_LOCAL_FOCK
ckbn @did : localize fock
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
          call cxsd_t2(k_f1_local,d_r2m(iref),k_t1_local,
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),size_r2m(iref),mypgid)
#else
          call cxsd_t2(k_f1_local,d_r2m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),size_r2m(iref),mypgid)
#endif
#else
          call cxsd_t2(d_f1m(iref),d_r2m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),size_r2m(iref),mypgid)
#endif
       else
          call createfile(filename,d_r2m(iref),size_r2m(iref))
          call mr_ccsd_t2(d_f1m(iref),d_r2m(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref),size_r2m(iref))
       endif
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step9 iter ref',iter,iref
c       call util_flush(6)
c      endif
c kk ---
c
c     - T1/X1 LOCALIZATION --
#if 0       
          if(ic_cc) then
           if(.not.MA_POP_STACK(l_t1_local))
     &       call errquit('l_t1_local',2,MA_ERR)
        end if
#endif
c -----------------------
c
c kk ---

c          if(nodezero) write(*,*) "I am here4"

c_blocked_for_now      if(lfullheff) then
c_blocked_for_now
c_blocked_for_now          itcpu = - util_cpusec()
c_blocked_for_now          itwall = - util_wallsec()
c_blocked_for_now
c_blocked_for_now          call tce_filename('r3act'//namechunk,filename)
c_blocked_for_now          call createfile(filename,d_r3am(iref),size_r3am(iref))
c_blocked_for_nowc kk--
c_blocked_for_now          call ccsd_activet3(d_f1m(iref),d_r3am(iref),d_t1m(iref),
c_blocked_for_now     1                 d_t2m(iref),d_v2m(iref),
c_blocked_for_now     1                 k_f1_offsetm(iref),k_r3a_offsetm(iref),
c_blocked_for_now     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
c_blocked_for_now     2                 k_v2_offsetm(iref))
c_blocked_for_now          call reconcilefile(d_r3am(iref),size_r3am(iref))
c_blocked_for_now
c_blocked_for_now          itcpu = itcpu + util_cpusec()
c_blocked_for_now          itwall = itwall + util_wallsec()
c_blocked_for_now
c_blocked_for_now          if(nodezero) then
c_blocked_for_now            write(LuOut,"('R3 active: ',2f15.1)") itcpu, itwall
c_blocked_for_now            call util_flush(LuOut)
c_blocked_for_now          endif
c_blocked_for_now
c_blocked_for_now          itcpu = - util_cpusec()
c_blocked_for_now          itwall = - util_wallsec()
c_blocked_for_now
c_blocked_for_now          call tce_filename('r4act'//namechunk,filename)
c_blocked_for_now          call createfile(filename,d_r4am(iref),size_r4am(iref))
c_blocked_for_now
c_blocked_for_nowc kk--
c_blocked_for_now          call ccsd_activet4(d_r4am(iref),d_t1m(iref),d_t2m(iref),
c_blocked_for_now     1                 d_v2m(iref),k_r4a_offsetm(iref),
c_blocked_for_now     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
c_blocked_for_now     3                 k_v2_offsetm(iref))
c_blocked_for_now          call reconcilefile(d_r4am(iref),size_r4am(iref))
c_blocked_for_now
c_blocked_for_now          itcpu = itcpu + util_cpusec()
c_blocked_for_now          itwall = itwall + util_wallsec()
c_blocked_for_now
c_blocked_for_now          if(nodezero) then
c_blocked_for_now            write(LuOut,"('R4 active: ',2f15.1)") itcpu, itwall
c_blocked_for_now            call util_flush(LuOut)
c_blocked_for_now          endif
c_blocked_for_nowc kk--
c_blocked_for_now        endif

c          if (nodezero) call util_flush(LuOut)

c          call tce_residual_t1(d_r1m(iref),k_r1_offsetm(iref),r1(iref))
c          call tce_residual_t2(d_r2m(iref),k_r2_offsetm(iref),r2(iref))
c
c PNNL: kill global arrays for d_f1m(iref) & d_v2m(iref)
c DISKINT
          if(diskint.or.luseeaf2e) then
c           call deletefile(d_f1m(iref))
            call deletefile(d_v2m(iref))
          end if

         endif ! subgroups
c
        enddo !iref

        if(lusesub) call ga_sync()

c        do iref=1,nref
c        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
c     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then
c        if((lusesub.and.(ga_pgroup_nodeid(int_mb(k_innodes
c     1 +innodes+ga_nodeid())).eq.0)).or.
c     2 ((.not.lusesub).and.nodezero)) then
c            call tce_mrcc_debug_pfile(d_r2m(iref),size_r2m(iref),
c     1 'R2er',iter,iref)
c            call tce_mrcc_debug_pfile(d_r1m(iref),size_r1m(iref),
c     1 'R1er',iter,iref)
c        endif
c        endif
c          enddo
c
c     ==========================
c        Construction of Heff
c     ==========================
c
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step11 iter',iter
c       write(6,*)'ga mem used11',ga_inquire_memory()
c       write(6,*)'ga mem avail11',ga_memory_avail()
c       call util_flush(6)
c      endif
c kk ---
        do i=1,nref*nref
          dbl_mb(k_heff+i-1) = 0.0d0
        enddo

c          itcpu = - util_cpusec()
c          itwall = - util_wallsec()

c        if(nodezero) write(*,*) "I am here5"

        call tce_heff(d_em,k_e_offsetm,k_r1_offsetm,
     1       k_r2_offsetm,k_r3a_offsetm,k_r4a_offsetm,d_r1m,d_r2m,
     2       d_r3am,d_r4am,.true.,.true.,lfullheff,lfullheff,rtdb)

c        if(nodezero) write(*,*) "I am here6"


c          itcpu = itcpu + util_cpusec()
c          itwall = itwall + util_wallsec()

c          if(nodezero.and.lusesub) then
c            write(LuOut,"('HEFF: ',2f15.1)") itcpu, itwall
c            call util_flush(LuOut)
c          endif


c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step12 iter',iter
c       write(6,*)'ga mem used12',ga_inquire_memory()
c       write(6,*)'ga mem avail12',ga_memory_avail()
c       call util_flush(6)
c      endif
c kk ---
c          itcpu = - util_cpusec()
c          itwall = - util_wallsec()

c        if(nodezero) write(*,*) "I am here7"

c       call ga_sync()
        call tce_diagonalize_heff(rtdb,iter)
c       call ga_sync()
ckbn subg -4
c        write(LuOut,'(A,F17.10,I5,I5)') "ckbn epsilon",
c     + epsilon,iter,ga_nodeid()
c        call util_flush(LuOut)

c          itcpu = itcpu + util_cpusec()
c          itwall = itwall + util_wallsec()

c          if(nodezero) then
c            write(LuOut,"('DIAG: ',2f15.1)") itcpu, itwall
c            call util_flush(LuOut)
c          endif

c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step13 iter',iter
c       write(6,*)'ga mem used13',ga_inquire_memory()
c       write(6,*)'ga mem avail13',ga_memory_avail()
c       call util_flush(6)
c      endif
c kk ---

c        if(nodezero) write(*,*) "I am here8"

        if (nodezero) write(LuOut,"(/,'BWCC iter. #',I4,2f24.13,2f24.13,
     1     2f18.13)")iter,epsilon,enrep+epsilon,epsilon-edelta
ckbn        dsummary(iter,1) = enrep+epsilon
ckbn        dsummary(iter,2) = epsilon-duhfens(mkroot)
c
c     ==========================
c             Converged
c     ==========================
c
ckbn -4
c        if(nodezero) write(LuOut,"('ddot Rckbn:',3F16.12)")
c     1    dsqrt(ddotr1),dsqrt(ddotr2),(sqrt(ddotr2)+sqrt(ddotr1))

ckbn sub-2
ckbn sub-2        
c        call ga_sync()

        if(((abs(epsilon-edelta).lt.thresh).and.
     1    ((sqrt(ddotr2)+sqrt(ddotr1)).lt.thresh)).or.lconverged) then

         if(.not.lconverged) then
          lconverged = .true.
          if(nodezero) write (LuOut,
     +      "(/,'The MR-BWCCSD iteration has converged.',/)")

ckbnttttt
c         if(.not. lusesub) then
          if(lsavet) then
           if(nodezero) then
            write(LuOut,*) "Stored amplitudes are before ap correction"
            call util_flush(LuOut)
           endif
           do iref=1,nref
            call tce_filenameindexed(iref,'mrt1',filename)
            unitn=78
            call write_tensor(filename,d_t1m(iref),size_t1m(iref),unitn)
c           close(unitn)
            call tce_filenameindexed(iref,'mrt2',filename)
            unitn=79
            call write_tensor(filename,d_t2m(iref),size_t2m(iref),unitn)
c           close(unitn)
           enddo
          endif
c         endif
ckbnttttt

c                    Universal SS Correction
c =====================================================================
c
          if(lap_uss) then
             dbwcorrf = 0.0d0
             dbwuss = 0.0d0

             do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)

            noa = nblcks(1,iref)
            nob = nblcks(2,iref)
            nva = nblcks(3,iref)
            nvb = nblcks(4,iref)
c
            noab = noa+nob
            nvab = nva+nvb

          itcpu3 = - util_cpusec()
          itwall3 = - util_wallsec()

               if(.not.lusesub) then

               call tce_filename('corr',filename)
               call createfile(filename,d_corr,size_em(iref))
C
               call tce_filename('corr1',filename)
               call createfile(filename,d_corr1,size_em(iref))

               call  c1_d1(d_r1m(iref),d_t1m(iref),d_corr,
     1 k_r1_offsetm(iref),k_t1_offsetm(iref),k_e_offsetm(iref))
        call c2_d2(d_r2m(iref),d_t2m(iref),d_corr,k_r2_offsetm(iref),
     1 k_t2_offsetm(iref),k_e_offsetm(iref))

               call get_block(d_corr,dbwcorr,1,0)

          dbwcorrf = dbwcorrf+dbwcorr*dbl_mb(k_sqc+(iref-1)*nref+
     1 mkroot-1)*dbl_mb(k_sqcl+(iref-1)*nref+mkroot-1)

       call deletefile(d_corr)

      if(.not.lussd) then

        do jref=1,nref
c
      call tce_filename('r1d',filename)
      call createfile(filename,d_r1d(jref),size_r1m(jref))
      call tce_filename('r2xy',filename)
      call createfile(filename,d_r2xy(jref),size_r2m(jref))
      call tce_filename('corr1',filename)
      call createfile(filename,d_corr1,size_em(iref))
! 
       if(iref.ne.jref) then

      call tce_uss_offdiagonal_1(d_r1m,d_r2m,k_r1_offsetm,
     1 k_r2_offsetm,iref,jref,d_r1d(jref),
     1 d_r2xy(jref))

       call tce_uss_offdiagonal_2(d_r1m,d_r2m,k_r1_offsetm,
     1 k_r2_offsetm,iref,jref,d_r2xy(jref),
     1 d_r1d(jref))
c
       call tce_uss_offdiagonal_3(d_t2m(iref),d_v2m(iref),
     2                 k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),d_r1d(jref),d_r2xy(jref),
     2                 iref,jref,k_r1_offsetm,k_r2_offsetm)


       call tce_uss_offdiagonal_4(d_t2m(iref),d_v2m(iref),
     2                 k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),d_r1d(jref),d_r2xy(jref),
     2                 iref,jref,k_r1_offsetm,k_r2_offsetm)


      call signf(iref,jref,signfact)
      dsmult=1.0d0
      if(mod(signfact,2).eq.1)then
      dsmult=-1.0d0
      endif

      call tce_uss_T1_R1(d_r1d(jref),d_t1m(jref),d_corr1,
     1k_r1_offsetm(jref),k_t1_offsetm(jref),k_e_offsetm(jref),iref,jref)

      call tce_uss_T2_R2(d_r2xy(jref),d_t2m(jref),d_corr1,
     1k_r2_offsetm(jref),k_t2_offsetm(jref),k_e_offsetm(jref),iref,jref)

       call get_block(d_corr1,dbwcorr1,1,0)
      dbwuss = dbwuss +dsmult*dbwcorr1 * dbl_mb(k_sqc+(iref-1)*nref+
     1 mkroot-1)*dbl_mb(k_sqcl+(jref-1)*nref+mkroot-1)
c
      endif
!
       call deletefile(d_r1d)
       call deletefile(d_r2xy)
       call deletefile(d_corr1)
        end do ! jref
        endif
        endif!ussd
        endif
       enddo !iref

      if(nodezero)write(6,"(/,'USSD Corr         :',F18.12)")dbwcorrf
      if(nodezero)write(6,"(/,'Total USSD ENERGY : ',2F18.12)")
     1 dbwcorrf+epsilon+enrep
c
      if(.not.lussd) then

      if(nodezero)write(6,"(/,'USS(pt) Corr      : ',F18.12)")
     1 dbwuss+dbwcorrf
      if(nodezero)write(6,"(/,'BW-MRCCSD USS(pt) : ',F18.12,/)")
     1 epsilon+enrep+dbwcorrf+dbwuss

       endif

        endif

c
c =====================================================================

cxxx
ckbn-2 avoid a posteriori correction
c          if(nodezero)write(*,'(A,3F17.12)') 
c     +        "Tckbn ", epsilon,edelta,abs(epsilon-edelta)
          if(nodezero) 
     +      write (LuOut,"('No a posteriori correction to Energy: Iter',
     +      I4,3f24.13)") iter, epsilon, enrep+epsilon, epsilon-edelta
ckbn subg -2
          call ga_sync()
          if(no_aposteriori) goto 1333
cxxx

          goto 1222  ! a posteriori

         endif ! not lconverged
ckbnjjjj
ckbn*******************************************************************
ckbn**************************BWCCSD(T)********************************
ckbn*******************************************************************
 1333    continue ! After avoiding a posteriori correction

cx       if(lconverged)write(*,*) "Before", lconverged,ga_nodeid()
       if(se4t .and. lconverged )  then
ckbn subg -4
c       call dfill(maxref,0.d0,pt3_1m,1)
c       call dfill(maxref,0.d0,pt3_2m,1)
       call dfill(maxref,0.d0,pt3_3m,1)
       if(ga_pgroup_nodeid(mypgid).eq.0) call util_flush(LuOut)
       call ga_sync()
       if(lusesub)call ga_pgroup_sync(mypgid)
          itcpu3 = - util_cpusec()
          itwall3 = - util_wallsec()
          itcpu3_1 =  util_cpusec()
          itwall3_1 =  util_wallsec()

        if(lusesub) then
        if( ga_pgroup_nodeid(mypgid).eq.0) then
         if (.not.MA_PUSH_GET(mt_dbl,nref*nref,'(T) ptheff',l_ptheff,
     1        k_ptheff)) call errquit('ptheff: MA error',3,MA_ERR)
         call dfill(nref*nref,0.d0,dbl_mb(k_ptheff),1)
        endif
        endif

c        if(lusesub) call ga_print(g_heff)

        do irefpt=1,nref

        if((int_mb(k_refafi+irefpt-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then

         k_sym = k_symm(irefpt)
         k_offset = k_offsetm(irefpt)
         k_range = k_rangem(irefpt)
         k_spin = k_spinm(irefpt)
         k_movecs_sorted = k_movecs_sortedm(irefpt)
         k_evl_sorted = k_evl_sortedm(irefpt)
         k_active = k_active_tmpm(irefpt)

         noa = nblcks(1,irefpt)
         nob = nblcks(2,irefpt)
         nva = nblcks(3,irefpt)
         nvb = nblcks(4,irefpt)
         noab = noa+nob
         nvab = nva+nvb
       if(lusesub) then

        if (.not.MA_PUSH_GET(mt_dbl,size_r1m(irefpt),'t1_localpt',
     1      l_t1_localpt,k_t1_localpt))
     1      call errquit('t1_localpt',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_localpt),size_r1m(irefpt))
c    copy d_t1 ==> k_t1_local
         if( ga_pgroup_nodeid(mypgid) .eq. 0) then
          call ga_get(d_t1m(irefpt),1,size_r1m(irefpt),1,1,
     &                dbl_mb(k_t1_localpt),1) 
        endif
          call ga_pgroup_brdcst(mypgid,65540,dbl_mb(k_t1_localpt),
     +     MA_sizeof(MT_DBL,1,MT_BYTE)*size_r1m(irefpt)*1,0)

        if (.not.MA_PUSH_GET(mt_dbl,size_1em(irefpt),'f1_localpt',
     1      l_f1_localpt,k_f1_localpt))
     1      call errquit('f1_localpt',1,MA_ERR)
        call ma_zero(dbl_mb(k_f1_localpt),size_1em(irefpt))
c    copy d_f1 ==> k_f1_local
         if( ga_pgroup_nodeid(mypgid) .eq. 0) then
          call ga_get(d_f1m(irefpt),1,size_1em(irefpt),1,1,
     &                dbl_mb(k_f1_localpt),1) 
        endif
          call ga_pgroup_brdcst(mypgid,65541,dbl_mb(k_f1_localpt),
     +     MA_sizeof(MT_DBL,1,MT_BYTE)*size_1em(irefpt)*1,0)

       endif
       
       if(lusesub) then
        call tce_mrcc_ccsdpt_subg(k_t1_localpt,k_t1_offsetm(irefpt),
     +                           d_t2m(irefpt),k_t2_offsetm(irefpt),
     +                           d_v2m(irefpt),k_v2_offsetm(irefpt),
     +                           k_f1_localpt,
     +                           k_f1_offsetm(irefpt),
c     +                           pt3_1m(irefpt),pt3_2m(irefpt),
     +                           pt3_3m(irefpt),irefpt)

         else
          call tce_mrcc_ccsdpt(d_t1m(irefpt),k_t1_offsetm(irefpt),
     +                           d_t2m(irefpt),k_t2_offsetm(irefpt),
     +                           d_v2m(irefpt),k_v2_offsetm(irefpt),
     +                           d_f1m(irefpt),
     +                           k_f1_offsetm(irefpt),
c     +                           pt3_1m(irefpt),pt3_2m(irefpt),
     +                           pt3_3m(irefpt),irefpt)
       endif

        if(lusesub) then
c        if( ga_pgroup_nodeid(mypgid) .eq. 0) write(LuOut,'(A,I10)') 
c     +     "I am here my pgroup_nodeid ",ga_pgroup_nodeid(mypgid)
        call util_flush(LuOut)
        if( ga_pgroup_nodeid(mypgid) .eq. 0) then
         dbl_mb(k_ptheff+((irefpt-1)*nref+(irefpt-1)))=pt3_3m(irefpt)
        endif
        endif

        if(lusesub) then
         if(.not.MA_POP_STACK(l_f1_localpt))
     &      call errquit('l_t1_localpt',2,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_localpt))
     &      call errquit('l_t1_localpt',2,MA_ERR)
        endif

c        if(lusesub) then
c        if( ga_pgroup_nodeid(mypgid).eq.0) then
c         write(LuOut,'(A,I5,F15.10)')"Correction to Heff from iref " ,
c     +     irefpt,pt3_3m(irefpt)
c        call util_flush(LuOut)
c        endif
c        else 
c        if(nodezero) then
c         write(LuOut,'(A,I5,F15.10)')"Correction to Heff from iref " ,
c     +     irefpt,pt3_3m(irefpt)
c        call util_flush(LuOut)
c        endif 
c        endif ! lusesub

       endif ! .not. lusesub
       end do

       if(lusesub) then
        if( ga_pgroup_nodeid(mypgid).eq.0)
     +    call ga_acc(g_heff,1,nref*nref,1,1,dbl_mb(k_ptheff),1,1.0d0)
       endif
       call ga_sync()

      if(.not.lusesub) then 
      do i=1,nref
       dbl_mb(k_heff+((i-1)*nref+(i-1)))=
     +       dbl_mb(k_heff+((i-1)*nref+(i-1)))+pt3_3m(i)
      enddo
      endif

        if(lusesub) then
        if( ga_pgroup_nodeid(mypgid).eq.0) then
         if (.not.MA_POP_STACK(l_ptheff)) 
     1    call errquit('ccsd_t',3,MA_ERR)
        endif
        endif

c        if(lusesub) call ga_print(g_heff)

       if(nodezero) then
        write(LuOut,"(/,'Final1 MR-BWCCSD energy before (T) correction:
     +     ',2F18.12)")epsilon+enrep
       endif

        call tce_diagonalize_heff(rtdb,iter)

       if(nodezero) then
        write(LuOut,"(/,'Final2 MR-BWCCSD energy after (T) correction:
     +     ',2F18.12)")epsilon+enrep
       endif

         itcpu3 = itcpu3 + util_cpusec()
         itwall3 = itwall3 + util_wallsec()
         itcpu3_2 =  util_cpusec()
         itwall3_2 = util_wallsec()
         call ga_dgop(mt_dbl,itcpu3 ,1,'max')
         call ga_dgop(mt_dbl,itwall3,1,'max')
         if(nodezero) then
            write(LuOut,"('Iter cpu (t) ',2f15.5)") itcpu3, itwall3
           write(LuOut,"('Iter cpu (t)n ',2f15.5)") (itcpu3_2-itcpu3_1),
     +                     (itwall3_2-itwall3_1)
            call util_flush(LuOut)
         endif

       endif

ckbn*******************************************************************

cxxx
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()

         do iref=1,nref
          if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1      +innodes+ga_nodeid())).or.(.not.lusesub)) then
           call deletefile(d_em(iref))
           call deletefile(d_r1m(iref))
           call deletefile(d_r2m(iref))
           if(lfullheff) then
            call deletefile(d_r3am(iref))
            call deletefile(d_r4am(iref))
           endif
          endif !sub
         enddo ! nref
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step14 iter',iter
c       call util_flush(6)
c      endif
c kk ---

         call tce_mrcc_diis_tidy()

         if(nodezero) then
c     +     write(LuOut,"(/,'Final MR-BWCCSD energy and ecorr: ',2F18.12,
c     1     2F18.12)")epsilon+enrep,epsilon-duhfens(1)
c          write(LuOut,"('A posteriori correction: ',2F18.12,
c     1     2F18.12,/)")epsilon-edelta
ckbn          call tce_mrcc_print_summary(iter,dsummary)

          write(LuOut,*)""
          write(LuOut,'(A)')
     +"=================================================="//
     +"==================================================="
          write(LuOut,*)""
          if( .not. se4t ) then
           if(.not. no_aposteriori) write(LuOut,9080) 
     +      "BW-MRCCSD with a posteriori correction ",
     +      epsilon+enrep
           if(no_aposteriori) write(LuOut,9080)
     +      "BW-MRCCSD without a posteriori correction ",
     +      epsilon+enrep
          else
           if(.not. no_aposteriori) write(LuOut,9080)
     +      "BW-MRCCSD(T) with a posteriori correction ",
     +      epsilon+enrep
           if(no_aposteriori) write(LuOut,9080)
     +      "BW-MRCCSD(T) without a posteriori correction ",
     +      epsilon+enrep
          endif
          if(iter.le.6) write(LuOut,*) "Converged below 6 iterations!"//
     +     " Please rerun with higher convergence threshold."
          write(LuOut,*)""
          write(LuOut,'(A)')
     +"=================================================="//
     +"==================================================="
          write(LuOut,*)""
          write(LuOut,9220) cpu, wall
          call util_flush(LuOut)
         endif ! nodezero 

         tce_mrcc_energy = .true.
         if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,epsilon+enrep))
     1        call errquit('tce_energy: RTDB problem',0,MA_ERR)
         if (nodezero) call util_flush(LuOut)

c      call ga_print(d_f1m(1))
c      call ga_print(d_f1m(2))
c      call ga_print(d_f1m(3))
c      call ga_print(d_f1m(4))
         goto 1000

        endif ! energy convergence check
c
c     ==========================
c             Update T's  
c     ==========================
c
 1222   continue
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step15 iter',iter
c       call util_flush(6)
c      endif
c kk ---

        do iref=1,nref
c        write(6,*)'TEST1',ga_nodeid(),iref
         if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1    +ga_nnodes()+ga_nodeid())).or.(.not.lusesub)) then

          g_movecs(1) = g_movecsm(iref,1)
          g_movecs(2) = g_movecsm(iref,2)

          k_sym = k_symm(iref)
          k_offset = k_offsetm(iref)
          k_range = k_rangem(iref)
          k_spin = k_spinm(iref)
          k_movecs_sorted = k_movecs_sortedm(iref)
          k_evl_sorted = k_evl_sortedm(iref)

          noa = nblcks(1,iref)
          nob = nblcks(2,iref)
          nva = nblcks(3,iref)
          nvb = nblcks(4,iref)

          noab = noa+nob
          nvab = nva+nvb

          mrccshift = epsilon - dbl_mb(k_heff+iref-1+(iref-1)*nref) ! Global denominator shift
c         write(6,*)'BWSHIFT: ',mrccshift,ga_nodeid(),int_mb(k_innodes+
c     1 ga_nnodes()+ga_nodeid())
          if(lconverged) mrccshift = 0.0d0 ! a posteriori corr

c
c     ==========================
c       Add DCl,L and UL terms
c     ==========================
c
          if(.not.lconverged) then      ! a posteriori corr
c           write(*,*) "I am here 1 dcl ul"

           itcpu = - util_cpusec()
           itwall = - util_wallsec()

         if(lusesub) then
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
          call xalf_t1t1_1(d_r1m(iref),k_r1_offsetm(iref),k_t1_local,
     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
#else
          call xalf_t1t1_1(d_r1m(iref),k_r1_offsetm(iref),d_t1m(iref),
     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
#endif
         else
           call mr_half_t1t1_1(d_r1m(iref),k_r1_offsetm(iref),
     +           d_t1m(iref),k_t1_offsetm(iref),d_r2m(iref),
     +           k_r2_offsetm(iref))
           call reconcilefile(d_r2m(iref),size_r2m(iref))
         endif

         if(lusesub) then
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
          call xcl_t1t1_1(k_t1_local,k_t1_offsetm(iref),k_t1_local,
     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
#else
          call xcl_t1t1_1(d_t1m(iref),k_t1_offsetm(iref),d_t1m(iref),
     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
          call ga_pgroup_sync(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
#endif
         else
           call dcl_t1t1_1(d_t1m(iref),k_t1_offsetm(iref),d_t1m(iref),
     &           k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
           call reconcilefile(d_r2m(iref),size_r2m(iref))
         endif

           itcpu = itcpu + util_cpusec()   
           itwall = itwall + util_wallsec()


          endif ! not lconverged
c kk ---
c      if(nodezero) then
c       write(6,*)'mrcc step16 iter',iter
c       call util_flush(6)
c      endif
c kk ---
c        write(6,*)'TEST2',ga_nodeid(),iref
          itcpu = - util_cpusec()
          itwall = - util_wallsec()

c          residual = residual + sqrt(ga_ddot(d_r2m(iref),d_r2m(iref)))/
c     1 size_t2m(iref)
cnewdii      call tce_mrcc_diis_new(.false.,iter,.true.,.true.,.true.,.false.,
cnewdii     1  d_r1m,d_t1m,k_t1_offsetm,size_t1m,
cnewdii     2  d_r2m,d_t2m,k_t2_offsetm,size_t2m,
cnewdii     3  d_r3m,d_t3m,k_t3_offsetm,size_t3m,
cnewdii     4  dummy,dummy,dummy,dummy,
cnewdii     5  rtdb,iref)

          call tce_mrcc_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1m(iref),d_t1m(iref),k_t1_offsetm(iref),
     1                  size_t1m(iref),d_r2m(iref),d_t2m(iref),
     2                  k_t2_offsetm(iref),size_t2m(iref),
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy,rtdb,iref)


          itcpu = itcpu + util_cpusec()
          itwall = itwall + util_wallsec()



         endif !sub
        enddo !iref


ckbn #ifdef MRCC_LOCAL_T1 as update was done on d_t1(iref) in tce_jacobi_sub_t1 zero global part
        call tce_internal_t_zero(d_t1m,d_t2m,k_t1_offsetm,
     1        k_t2_offsetm,.false.,dummy,dummy,rtdb)

c
c     ==========================
c             DDOT conv.
c     ==========================
c

c        if(.not.lusesub) then
c
c        call tce_internal_t_zero(d_r1m,d_r2m,k_t1_offsetm,
c     1 k_t2_offsetm,nref,.false.,dummy,dummy,rtdb)

c          ddotr1 = 0.0d0
c          ddotr2 = 0.0d0

c        do iref=1,nref
c
c          ddotr1 = ga_ddot(d_r1m(iref),d_r1m(iref))
c          ddotr2 = ga_ddot(d_r2m(iref),d_r2m(iref))

c        enddo
c
c         if(nodezero)write(LuOut,"('ddot R:',2F16.12)")sqrt(ddotr1),
c     1 sqrt(ddotr2)

c        endif


        do iref=1,nref
         if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1    +ga_nnodes()+ga_nodeid())).or.(.not.lusesub)) then
          call deletefile(d_em(iref))
          call deletefile(d_r1m(iref))
          call deletefile(d_r2m(iref))
          if(lfullheff) then
           call deletefile(d_r3am(iref))
           call deletefile(d_r4am(iref))
          endif
          if (nodezero) call util_flush(LuOut)
         endif
        enddo


        edelta = epsilon

        itcpu2 = itcpu2 + util_cpusec()
        itwall2 = itwall2 + util_wallsec()

        if(nodezero) then
          write(LuOut,9222) itcpu2, itwall2
          call util_flush(LuOut)
        endif

        if(lcheckpoints) then
         if(mod(iter,icheckpoint).eq.0) then
          do iref=1,nref
           call tce_filenameindexed(iref,'mrt1',filename)
           unitn=76
           call write_tensor(filename,d_t1m(iref),size_t1m(iref),unitn)
c          close(unitn)
           call tce_filenameindexed(iref,'mrt2',filename)
           unitn=77
           call write_tensor(filename,d_t2m(iref),size_t2m(iref),unitn)
c          close(unitn)
          enddo
         endif
        endif

       enddo !maxiter

       call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)

 1000  continue
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd end
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd end
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd end
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd end
c bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd bwccsd end

c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd begin 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd begin 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd begin 

c
c     ===============
c     Mk CCSD method
c     ===============
c
      else if (model.eq.'mkccsd') then

       if(nodezero) write (LuOut,"('MR MkCCSD, version 1.0')")

ckbn
       if(lusesub) then
        if(nodezero) write (LuOut,
     +    "('In this version MR MkCCSD in subgroup doesnt work')")
        call errquit('Subgroups not present in MR MKCCSD',1,MA_ERR)
       endif

       cpu = - util_cpusec()
       wall = - util_wallsec()

       call tce_mrcc_diis_init()

       edelta =0.0d0

c      call ga_sync()

c -----------------
c  Loop over iters
c -----------------
c
       do iter=1,maxiter

        itcpu = - util_cpusec()
        itwall = - util_wallsec()


cckbn root to overlap intialization
c        if(.not. lusesub) then !ckbn can be removed
cc         if(nodezero) then
cc         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+0)
cc         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+1)
cc         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+0)
cc         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+1)
cc         endif
c         if(iter.eq.3) then
c          if(iroottooverlap.gt.0) then
c           if(nodezero) write(*,*) "iroottooverlap ", iroottooverlap
c           do i=1,nref
c            bwcoefwanted(i) = dbl_mb(k_sqc+(i-1)*nref+iroottooverlap-1) 
c           enddo
c           lrootfromoverlap = .true.
c           do i=1,nref
c            if(nodezero)
c     +        write(*,'(A,F17.10)')"bwcoefwanted",bwcoefwanted(i)
c           enddo
c          endif
c          mkrootold = mkroot
c         endif

ckbn root to overlap intialization
        if(.not. lusesub) then !ckbn can be removed
c         if(nodezero) then
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+0)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(0)*nref+1)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+0)
c         write(*,*)"ckbn iter",iter,dbl_mb(k_sqc+(1)*nref+1)
c         endif
         if(iter.eq.iroottooverlapiter) then
          if(iroottooverlap.gt.0) then
c           if(nodezero) write(*,*) "iroottooverlap ", iroottooverlap
           if(nodezero) write(LuOut,'(A,I5,A,I5)')
     +      "Overlap selected from iteration ", iroottooverlapiter,
     +      " with respect to root ", iroottooverlap
           do i=1,nref
            bwcoefwanted(i) = dbl_mb(k_sqc+(i-1)*nref+iroottooverlap-1)
           enddo
           lrootfromoverlap = .true.
           do i=1,nref
            if(nodezero) write(*,'(A,I5,F17.10)')
     +       "Overlapping right eigen vector ",i,bwcoefwanted(i)
           enddo
          else
           call errquit('root to overlap not natural number',0,RTDB_ERR)
          endif
          mkrootold = mkroot
         endif



ckbn check if roots swapped
         if(iter.gt.4) then
          if(nodezero) then
c           write(LuOut,'(A,I3,A,I3)')"mkroot",mkroot,
c     +      "mkrootold",mkrootold
           call util_flush(LuOut)
           if(mkroot.ne.mkrootold) then
             write(LuOut,'(A,I5,A,I5)') 
     +      "Warning: Eigenvector swap occured from",
     +      mkrootold," to ",mkroot
           endif
          endif
          mkrootold = mkroot
         endif
         call util_flush(LuOut)
        endif ! ckbn can be removed

c
c -----------------
c  Loop over refs
c -----------------
c
        do iref=1,nref
c

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then

        if(lusesub)
     1   ipg = int_mb(k_innodes+ga_nnodes()+ga_nodeid())

         g_movecs(1) = g_movecsm(iref,1)
         g_movecs(2) = g_movecsm(iref,2)

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

          write(namechunk,"(I3.3)")iref
c
c - T1/X1 LOCALIZATION ----------
        if(ic_cc) then
        if (.not.MA_PUSH_GET(mt_dbl,size_r1m(iref),'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_r1m(iref))
c    copy d_t1 ==> l_t1_local
        call ga_get(d_t1m(iref),1,size_r1m(iref),1,1,
     &       dbl_mb(k_t1_local),1)
        end if

c kk ---
c PNNL: create global arrays for d_f1m(iref) & d_v2m(iref)
c PNNL: read on d_f1m(iref) d_v2m(iref) from a disk
c PNNL: set d_v2orb to the current handle of d_v2m(iref)
c
          d_v2orb=d_v2m(iref)

          call tce_filename('e'//namechunk,filename)
          if(lusesub) then
          call pgcreatefile(filename,d_em(iref),size_em(iref),mypgid)
          call cxsd_e(d_f1m(iref),d_em(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_e_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref),ipg)
          else
          call createfile(filename,d_em(iref),size_em(iref))
          call mr_ccsd_e(d_f1m(iref),d_em(iref),d_t1m(iref),
     1                d_t2m(iref),d_v2m(iref),
     1                k_f1_offsetm(iref),k_e_offsetm(iref),
     2                k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                k_v2_offsetm(iref))
          endif
          call tce_filename('r1'//namechunk,filename)
       if(lusesub) then
          call pgcreatefile(filename,d_r1m(iref),size_r1m(iref),mypgid)
          call cxsd_t1(d_f1m(iref),d_r1m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),mypgid)
       else
          call createfile(filename,d_r1m(iref),size_r1m(iref))
          call mr_ccsd_t1(d_f1m(iref),d_r1m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r1_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref))
       endif


          call tce_filename('r2'//namechunk,filename)
       if(lusesub) then
          call pgcreatefile(filename,d_r2m(iref),size_r2m(iref),mypgid)
          call cxsd_t2(d_f1m(iref),d_r2m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),size_r2m(iref),mypgid)
       else
          call createfile(filename,d_r2m(iref),size_r2m(iref))
          call mr_ccsd_t2(d_f1m(iref),d_r2m(iref),d_t1m(iref),
     1                 d_t2m(iref),d_v2m(iref),
     1                 k_f1_offsetm(iref),k_r2_offsetm(iref),
     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),size_r2m(iref))
       endif

c - T1/X1 LOCALIZATION --
        if(ic_cc) then
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',2,MA_ERR)
        end if
c -----------------------

c_blocked_for_now      if(lfullheff) then
c_blocked_for_now
c_blocked_for_now          call tce_filename('r3act'//namechunk,filename)
c_blocked_for_now          call createfile(filename,d_r3am(iref),size_r3am(iref))
c_blocked_for_nowc kk--
c_blocked_for_now          call ccsd_activet3(d_f1m(iref),d_r3am(iref),d_t1m(iref),
c_blocked_for_now     1                 d_t2m(iref),d_v2m(iref),
c_blocked_for_now     1                 k_f1_offsetm(iref),k_r3a_offsetm(iref),
c_blocked_for_now     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
c_blocked_for_now     2                 k_v2_offsetm(iref))
c_blocked_for_now          call reconcilefile(d_r3am(iref),size_r3am(iref))
c_blocked_for_now
c_blocked_for_now          call tce_filename('r4act'//namechunk,filename)
c_blocked_for_now          call createfile(filename,d_r4am(iref),size_r4am(iref))
c_blocked_for_now
c_blocked_for_nowc kk--
c_blocked_for_now          call ccsd_activet4(d_r4am(iref),d_t1m(iref),d_t2m(iref),
c_blocked_for_now     1                 d_v2m(iref),k_r4a_offsetm(iref),
c_blocked_for_now     2                 k_t1_offsetm(iref),k_t2_offsetm(iref),
c_blocked_for_now     3                 k_v2_offsetm(iref))
c_blocked_for_now          call reconcilefile(d_r4am(iref),size_r4am(iref))
c_blocked_for_nowc kk--
c_blocked_for_now        endif

          if (nodezero) call util_flush(LuOut)

         if(diskint) then
c          call deletefile(d_f1m(iref))
          call deletefile(d_v2m(iref))
         end if

        endif ! sub filter

        enddo !iref

        if(lusesub) call ga_sync()

c
c     ==========================
c        Construction of Heff
c     ==========================
c
        do i=1,nref*nref
          dbl_mb(k_heff+i-1) = 0.0d0
        enddo

        do i=1,nref*nref
          dbl_mb(k_sqc+i-1) = 0.0d0
        enddo


c        call tce_heff(nref,d_em,k_e_offsetm,k_r1_offsetm,
        call tce_heff(d_em,k_e_offsetm,k_r1_offsetm,
     1 k_r2_offsetm,k_r3a_offsetm,k_r4a_offsetm,d_r1m,d_r2m,
     2 d_r3am,d_r4am,.true.,.true.,lfullheff,lfullheff,rtdb)


c        call tce_diagonalize_heff(nref,rtdb)
        call tce_diagonalize_heff(rtdb,iter)

        if (nodezero) write(LuOut,"(/,'MkCC iter. #',I4,2f24.13,2f24.13,
     1 2f18.13)")iter,epsilon,enrep+epsilon,epsilon-edelta
        if (nodezero) call util_flush(LuOut)

ckbn        dsummary(iter,1) = enrep+epsilon
ckbn        dsummary(iter,2) = epsilon-duhfens(1)
c
c     ==========================
c             Converged
c     ==========================
c
ckbn diis ***
c         if(nodezero)write(LuOut,"('ddot Rckbn:',i5,4F16.12)")iref,
c     1 dsqrt(ddotr1),dsqrt(ddotr2),
c     1 sqrt(ddotr1),sqrt(ddotr2)
        if(((abs(epsilon-edelta).lt.thresh).and.
     1 ((sqrt(ddotr2)+sqrt(ddotr1)).lt.thresh)).or.lconverged) then

          if(.not.lconverged) then
             lconverged = .true.
             if(nodezero) then
             write(LuOut,"(/,'The MR-MkCCSD iteration has converged.',
     1 /)")
             call util_flush(LuOut)
             end if
c             goto 1223
          endif
ckbnttttt
        if(lsavet) then
          if(nodezero) 
     +     write(LuOut,*) "Stored amplitudes"
          call util_flush(LuOut)
          do iref=1,nref
           call tce_filenameindexed(iref,'mrt1',filename)
           unitn=78
           call write_tensor(filename,d_t1m(iref),size_t1m(iref),unitn)
c           close(unitn)
           call tce_filenameindexed(iref,'mrt2',filename)
           unitn=79
           call write_tensor(filename,d_t2m(iref),size_t2m(iref),unitn)
c           close(unitn)
          enddo
        endif
ckbnttttt

c  ==================================
c Universal State-Selective Correction
c  ==================================
c
c
          if(lap_uss) then
             dbwcorrf = 0.0d0
             dbwuss = 0.0d0

             do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)

            noa = nblcks(1,iref)
            nob = nblcks(2,iref)
            nva = nblcks(3,iref)
            nvb = nblcks(4,iref)
c
            noab = noa+nob
            nvab = nva+nvb
               if(.not.lusesub) then

               call tce_filename('corr',filename)
               call createfile(filename,d_corr,size_em(iref))
C
               call tce_filename('corr1',filename)
               call createfile(filename,d_corr1,size_em(iref))

               call  c1_d1(d_r1m(iref),d_t1m(iref),d_corr,
     1 k_r1_offsetm(iref),k_t1_offsetm(iref),k_e_offsetm(iref))
        call c2_d2(d_r2m(iref),d_t2m(iref),d_corr,k_r2_offsetm(iref),
     1 k_t2_offsetm(iref),k_e_offsetm(iref))

               call get_block(d_corr,dbwcorr,1,0)

          dbwcorrf = dbwcorrf+dbwcorr*dbl_mb(k_sqc+(iref-1)*nref+
     1 mkroot-1)*dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

       call deletefile(d_corr)

      if(.not.lussd) then

        do jref=1,nref
c
      call tce_filename('r1d',filename)
      call createfile(filename,d_r1d(jref),size_r1m(jref))
      call tce_filename('r2xy',filename)
      call createfile(filename,d_r2xy(jref),size_r2m(jref))
      call tce_filename('corr1',filename)
      call createfile(filename,d_corr1,size_em(iref))
! 
       if(iref.ne.jref) then

      call tce_uss_offdiagonal_1(d_r1m,d_r2m,k_r1_offsetm,
     1 k_r2_offsetm,iref,jref,d_r1d(jref),
     1 d_r2xy(jref))

       call tce_uss_offdiagonal_2(d_r1m,d_r2m,k_r1_offsetm,
     1 k_r2_offsetm,iref,jref,d_r2xy(jref),
     1 d_r1d(jref))
c
       call tce_uss_offdiagonal_3(d_t2m(iref),d_v2m(iref),
     2                 k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),d_r1d(jref),d_r2xy(jref),
     2                 iref,jref,k_r1_offsetm,k_r2_offsetm)


       call tce_uss_offdiagonal_4(d_t2m(iref),d_v2m(iref),
     2                 k_t2_offsetm(iref),
     2                 k_v2_offsetm(iref),d_r1d(jref),d_r2xy(jref),
     2                 iref,jref,k_r1_offsetm,k_r2_offsetm)


      call signf(iref,jref,signfact)
      dsmult=1.0d0
      if(mod(signfact,2).eq.1)then
      dsmult=-1.0d0
      endif

      call tce_uss_T1_R1(d_r1d(jref),d_t1m(jref),d_corr1,
     1k_r1_offsetm(jref),k_t1_offsetm(jref),k_e_offsetm(jref),iref,jref)

      call tce_uss_T2_R2(d_r2xy(jref),d_t2m(jref),d_corr1,
     1k_r2_offsetm(jref),k_t2_offsetm(jref),k_e_offsetm(jref),iref,jref)

       call get_block(d_corr1,dbwcorr1,1,0)
      dbwuss = dbwuss +dsmult*dbwcorr1 * dbl_mb(k_sqc+(iref-1)*nref+
     1 mkroot-1)*dbl_mb(k_sqcl+(jref-1)*nref+mkroot-1)


      endif
!
       call deletefile(d_r1d)
       call deletefile(d_r2xy)
       call deletefile(d_corr1)
        end do ! jref
        endif
        endif!ussd
        endif
       enddo !iref
c

c      if(nodezero)write(luOut,*)"USS_OFFDIAGONAL",dbwuss
      if(nodezero)write(6,"(/,'USSD Corr        :',F18.12)")
     1 dbwcorrf
      if(nodezero)write(6,"(/,'Mk-MRCCSD USSD   :',F18.12)")
     1 dbwcorrf+epsilon+enrep

      if(.not.lussd) then
 
      if(nodezero)write(6,"(/,'USS(pt) Corr     :',F18.12)")
     1 dbwuss+dbwcorrf
      if(nodezero)write(6,"(/,'Mk-MRCCSD USS(pt):',F18.12,/)")
     1 epsilon+enrep+dbwcorrf+dbwuss

       endif

          endif

c
c
c ============END of main USS code =====================

ckbn*******************************************************************
ckbn**************************MkCCSD(T)********************************
ckbn*******************************************************************
ckbn   SE4T for MkCC
c       if((se4t .and. lconverged) .or. (uncoup_pt3 .and. lconverged)) 
       if((se4t .and. lconverged) ) then
        call ga_sync()
        if(lusesub) call errquit(
     +   'tce_mrcc_energy: Subgroups not supported in Mk-MRCCSD(T)',30,
     +    INPUT_ERR)

          itcpu3 = - util_cpusec()
          itwall3 = - util_wallsec()
c       call dfill(maxref,0.d0,pt3_1m,1)
c       call dfill(maxref,0.d0,pt3_2m,1)
       call dfill(maxref,0.d0,pt3_3m,1)

        do irefpt=1,nref

         k_sym = k_symm(irefpt)
         k_offset = k_offsetm(irefpt)
         k_range = k_rangem(irefpt)
         k_spin = k_spinm(irefpt)
         k_movecs_sorted = k_movecs_sortedm(irefpt)
         k_evl_sorted = k_evl_sortedm(irefpt)
         k_active = k_active_tmpm(irefpt)
         noa = nblcks(1,irefpt)
         nob = nblcks(2,irefpt)
         nva = nblcks(3,irefpt)
         nvb = nblcks(4,irefpt)
         noab = noa+nob
         nvab = nva+nvb

c          if(nodezero) write(*,*) "Call tce_mrcc_ccsdpt"
          call tce_mrcc_ccsdpt(d_t1m(irefpt),k_t1_offsetm(irefpt),
     +                           d_t2m(irefpt),k_t2_offsetm(irefpt),
     +                           d_v2m(irefpt),k_v2_offsetm(irefpt),
     +                           d_f1m(irefpt),
     +                           k_f1_offsetm(irefpt),
c     +                           pt3_1m(irefpt),pt3_2m(irefpt),
     +                           pt3_3m(irefpt),irefpt)
c       if(nodezero) then
cckbn       write(*,*)"E4T  ","iref",irefpt," ",pt3_1m(irefpt) 
c        write(*,'(A,I5,F15.10)')"Correction to Heff from iref " ,
c     +    irefpt,pt3_3m(irefpt)
c       endif 


        end do

      do i=1,nref
       dbl_mb(k_heff+((i-1)*nref+(i-1)))=
     +       dbl_mb(k_heff+((i-1)*nref+(i-1)))+pt3_3m(i)
      enddo

       if(nodezero) then
        write(LuOut,"(/,'Final1 MR-MkCCSD energy before (T) correction:
     +     ',2F18.12)")epsilon+enrep
       endif

        call ga_sync()
        call tce_diagonalize_heff(rtdb,iter)

       if(nodezero) then
        write(LuOut,"(/,'Final2 MR-MkCCSD energy after (T) correction:
     +     ',2F18.12)")epsilon+enrep
       endif
         itcpu3 = itcpu3 + util_cpusec()
         itwall3 = itwall3 + util_wallsec()
         if(nodezero) then
            write(LuOut,"('Iter cpu (t) ',2f15.5)") itcpu3, itwall3
            call util_flush(LuOut)
         endif

       endif ! se4t
ckbn*******************************************************************



          do iref=1,nref
        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then
           call deletefile(d_em(iref))
           call deletefile(d_r1m(iref))
           call deletefile(d_r2m(iref))
           if(lfullheff)call deletefile(d_r3am(iref))
           if(lfullheff)call deletefile(d_r4am(iref))
        endif
          enddo

          call tce_mrcc_diis_tidy()

          if(nodezero) then
ckbn           call tce_mrcc_print_summary(iter,dsummary)

           write(LuOut,*) ""
           write(LuOut,'(A)')
     +"=================================================="//
     +"==================================================="
           write(LuOut,*)""
           if(.not. se4t) write(LuOut,9080)
     +     "Mk-MRCCSD ",epsilon+enrep

           if(se4t) write(LuOut,9080)
     +     "Mk-MRCCSD(T) ",epsilon+enrep

          if(iter.le.6) write(LuOut,*) "Converged below 6 iterations!"//
     +     " Please rerun with higher convergence threshold."

           write(LuOut,*) ""
           write(LuOut,'(A)')
     +"=================================================="//
     +"==================================================="
           write(LuOut,*)""
           call util_flush(LuOut)
          endif

           tce_mrcc_energy = .true.
           if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,epsilon+enrep))
     1        call errquit('tce_energy: RTDB problem',0,MA_ERR)

          if (nodezero) call util_flush(LuOut)

           goto 1001

        endif
c
c     ==============================
c             Update R's and T's  
c     ==============================
c
 1223   continue

        call tce_mrcc_diis_save(.false.,iter,
     1 .true.,.true.,.false.,.false.,
     1  d_r1m,d_t1m,k_t1_offsetm,size_t1m,
     2  d_r2m,d_t2m,k_t2_offsetm,size_t2m,
     3  dummy,dummy,dummy,dummy,
     4  dummy,dummy,dummy,dummy,
     5  rtdb)

cjb update r's
        do iref=1,nref

        if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then

         g_movecs(1) = g_movecsm(iref,1)
         g_movecs(2) = g_movecsm(iref,2)

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_evl_sorted = k_evl_sortedm(iref)
         k_active = k_active_tmpm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb


c
c     ==========================
c       Add DCl,L and UL terms
c     ==========================
c
c_blocked         if(iter.le.istartmk) then
c_blocked
c_blocked         lsubterm = .true.
c_blocked         diis = 0

c_blocked         mrccshift = epsilon - dbl_mb(k_heff+iref-1+(iref-1)*nref) ! Global denominator shift

c_blocked_for_now         if(lusesub) then
c_blocked_for_now          call xalf_t1t1_1(d_r1m(iref),k_r1_offsetm(iref),d_t1m(iref),
c_blocked_for_now     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
c_blocked_for_now          call ga_pgroup_sync(mypgid)
c_blocked_for_now         else
c_blocked          call half_t1t1_1(d_r1m(iref),k_r1_offsetm(iref),d_t1m(iref),
c_blocked     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
c_blocked          call reconcilefile(d_r2m(iref),size_r2m(iref))
c_blocked_for_now         endif
c_blocked_for_now
c_blocked_for_now         if(lusesub) then
c_blocked_for_now          call xcl_t1t1_1(d_t1m(iref),k_t1_offsetm(iref),d_t1m(iref),
c_blocked_for_now     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
c_blocked_for_now          call ga_pgroup_sync(mypgid)
c_blocked_for_now         else
c_blocked          call dcl_t1t1_1(d_t1m(iref),k_t1_offsetm(iref),d_t1m(iref),
c_blocked     & k_t1_offsetm(iref),d_r2m(iref),k_r2_offsetm(iref))
c_blocked          call reconcilefile(d_r2m(iref),size_r2m(iref))
c_blocked_for_now         endif
c_blocked_for_now
c_blocked         else

        lsubterm = .false.

         if(iter.eq.(istartmk+1)) then
      if (.not.rtdb_get(rtdb,'tce:diis',mt_int,1,diis)) then
        diis = 5
         endif
         endif

         cpu = - util_cpusec()
         wall = - util_wallsec()

c         call tce_mrcc_c1(size_t1m,k_t1_offsetm,d_t1m,nref,iref,
c     1 d_r1m,k_r1_offsetm)
         call tce_mrcc_c1(size_t1m,k_t1_offsetm,d_t1m,iref,
     1 d_r1m,k_r1_offsetm)

         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()

c         if(iref.eq.1) then
c            write(LuOut,"('C1 ',2f15.1,I4)") cpu, wall,iref
c            call util_flush(LuOut)
c         endif

         cpu = - util_cpusec()
         wall = - util_wallsec()

c         call tce_mrcc_c2(size_t2m,k_t1_offsetm,d_t1m,nref,iref,
c     1 k_t2_offsetm,d_t2m,d_r2m,k_r2_offsetm,size_t1m)
         call tce_mrcc_c2(size_t2m,k_t1_offsetm,d_t1m,iref,
     1 k_t2_offsetm,d_t2m,d_r2m,k_r2_offsetm,size_t1m)

         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()

c         if(nodezero) then
c         if(iref.eq.1) then
c            write(LuOut,"('C2 ',2f15.1,I4)") cpu, wall,iref
c            call util_flush(LuOut)
c         endif
c         endif

c_blocked        endif
        endif
        enddo

cjb continue in updating of T's 

        if(lusesub)call ga_sync()

        call tce_internal_t_zero(d_r1m,d_r2m,k_t1_offsetm,
     1 k_t2_offsetm,.false.,dummy,dummy,rtdb)

c             DDOT conv.
        if(.not.lusesub) then
          ddotr1 = 0.0d0
          ddotr2 = 0.0d0
         do iref=1,nref
          ddotr1 = ddotr1 + ga_ddot(d_r1m(iref),d_r1m(iref)) *
     2             dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
     3             dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

          ddotr2 = ddotr2 + ga_ddot(d_r2m(iref),d_r2m(iref)) *
     2             dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)*
     3             dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

c          if(nodezero) write(LuOut,"('ddot R:',i5,2F16.12)")iref,
c     1                      dsqrt(ddotr1),dsqrt(ddotr2)
         enddo
        endif

          if(nodezero) write(LuOut,"('ddot R:',2F16.12)")
     1                      dsqrt(ddotr1),dsqrt(ddotr2)


        call tce_mrcc_diis_new(.false.,iter,
     1       .true.,.true.,.false.,.false.,
     1       d_r1m,d_t1m,k_t1_offsetm,size_t1m,
     2       d_r2m,d_t2m,k_t2_offsetm,size_t2m,
     3       dummy,dummy,dummy,dummy,
     4       dummy,dummy,dummy,dummy,
     5       rtdb,iref)
ckbn
c        call ga_sync()


        do iref=1,nref
         if((int_mb(k_refafi+iref-1).eq.mypgid).or.(.not.lusesub)) then
          call deletefile(d_em(iref))
          call deletefile(d_r1m(iref))
          call deletefile(d_r2m(iref))
          if (lfullheff) call deletefile(d_r3am(iref))
          if (lfullheff) call deletefile(d_r4am(iref))
         endif
        enddo !iref

        edelta = epsilon

        if(lusesub)call ga_sync()

        itcpu = itcpu + util_cpusec()
        itwall = itwall + util_wallsec()
        if(nodezero) then
          write(LuOut,"('Iter cpu ',2f15.1,I4)") itcpu, itwall,iter
          call util_flush(LuOut)
        endif

        if(lcheckpoints) then
         if(mod(iter,icheckpoint).eq.0) then
          do iref=1,nref
           call tce_filenameindexed(iref,'mrt1',filename)
           unitn=76
           call write_tensor(filename,d_t1m(iref),size_t1m(iref),unitn)
c           close(unitn)
           call tce_filenameindexed(iref,'mrt2',filename)
           unitn=77
           call write_tensor(filename,d_t2m(iref),size_t2m(iref),unitn)
c           close(unitn)
          enddo
         endif
        endif

       enddo !maxiter

       call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)

 1001  continue

c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd end 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd end 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd end 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd end 
c mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd mkccsd end 

      else
       call errquit('tce_mrcc_energy: unkown method',INPUT_ERR)
      endif !methods 

c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup begin
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup begin
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup begin

ckbn------------------------------- clean up--------------------------------
ckbn If there is some final cleanup issue check it carefully it may be a bug
ckbn------------------------------- clean up--------------------------------
c
c     ===============
c     Destroy a mutex
c     ===============
c
c      if (.not.ga_destroy_mutexes())
c     1  call errquit('tce_energy: GA problem',1,GA_ERR)
c
c     ==================
c     Deallocate arrays
c     ==================
c
c      if(nodezero) write(*,*) "I am here1"
c      call ga_print(d_f1m(1))
c      call ga_print(d_f1m(2))
c      call ga_print(d_f1m(3))
c      call ga_print(d_f1m(4))

      do iref=1,nref

        if (.not.ga_destroy(g_movecsm(iref,1)))
     1    call errquit('tce_mrcc_energy: GA problem',1,GA_ERR)
        if (.not.ga_destroy(g_movecsm(iref,2)))
     1    call errquit('tce_mrcc_energy: GA problem',1,GA_ERR)

      enddo

#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
       if(lusesub) then
        do iref=1,nref
         if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     +      +innodes+ga_nodeid())).or.(.not.lusesub)) then
          if(.not.MA_POP_STACK(l_t1_local)) 
     +     call errquit('t1_local',1,MA_ERR)
         endif
        enddo
       endif
#endif


      if (.not.ma_pop_stack(l_sqcl))
     1 call errquit("tce_mrcc_energy: MA problem",331,MA_ERR)

      if (.not.ma_pop_stack(l_sqc))
     1 call errquit("tce_mrcc_energy: MA problem",332,MA_ERR)

      if (.not.ma_pop_stack(l_heff))
     1 call errquit("tce_mrcc_energy: MA problem",58,MA_ERR)

      if(lusesub) then
          lstatus = ga_destroy(g_heff)
      endif





      if((.not.diskint).and.(.not.luseeaf2e)) then
       do iref=1,nref
        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then
          call deletefile(d_v2m(iref))
        endif
       enddo 
      end if
c
      do iref=1,nref
        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub).or.
     2 (model.eq.'mkccsd')) then
         if(needt2) call deletefile(d_t2m(iref))
         if(needt1) call deletefile(d_t1m(iref))
         if(needt3) call deletefile(d_t3m(iref))
        endif
      enddo

      do i=1,nref

c       if(needr4act) then
c         if (.not.ma_pop_stack(l_r4act_offsetm(nref-i+1)))
cc     1    call errquit("tce_energy: MA problem",61,MA_ERR)
c       endif
       if(needt3) then
         if (.not.ma_pop_stack(l_r3_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",60,MA_ERR)
       endif
       if(needt2) then
         if (.not.ma_pop_stack(l_r2_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",58,MA_ERR)
        endif
        if(needt1) then
         if (.not.ma_pop_stack(l_r1_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",57,MA_ERR)
        endif

       if (.not.ma_pop_stack(l_e_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",59,MA_ERR)

      enddo

      do i=1,nref

        if((int_mb(k_refafi+nref-i).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub).or.
     2 (model.eq.'mkccsd')) then

       if(lfullheff) then
         if (.not.ma_pop_stack(l_r4a_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",56,MA_ERR)
         if (.not.ma_pop_stack(l_r3a_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",56,MA_ERR)
       endif
        if(needt1) then
         if (.not.ma_pop_stack(l_t1_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",55,MA_ERR)
        endif
       if(needt2) then
         if (.not.ma_pop_stack(l_t2_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",56,MA_ERR)
        endif
       if(needt3) then
         if (.not.ma_pop_stack(l_t3_offsetm(nref-i+1)))
     1    call errquit("tce_energy: MA problem",88,MA_ERR)
        endif

        endif

      enddo

c      if(nodezero) write(*,*) "I am here2"
c      call ga_print(d_f1m(1))
c      call ga_print(d_f1m(2))
c      call ga_print(d_f1m(3))
c      call ga_print(d_f1m(4))

c      do i=1,nref

c       if(lfullheff) then
c         if (.not.ma_pop_stack(l_r4a_offsetm(nref-i+1)))
c     1    call errquit("tce_energy: MA problem",88,MA_ERR)
c        endif
c       if(lfullheff) then
c         if (.not.ma_pop_stack(l_r3a_offsetm(nref-i+1)))
c     1    call errquit("tce_energy: MA problem",56,MA_ERR)
c        endif

c      enddo

      do i=1,nref

        if((int_mb(k_refafi+nref-i).eq.int_mb(k_innodes
     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then

      if (.not.ma_pop_stack(l_v2_offsetm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",54,MA_ERR)

        endif

      enddo

c      if(nodezero) write(*,*) "I am here3"
c      call ga_print(d_f1m(1))
c      call ga_print(d_f1m(2))
c      call ga_print(d_f1m(3))
c      call ga_print(d_f1m(4))
c PNNL
c#ifdef MRCC_LOCAL_FOCK
c      if(lusesub) then
c       do iref=1,nref
c        call deletefile(d_f1m(iref))
c       enddo
c      else
c       do iref=1,nref
c        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
c     +     +innodes+ga_nodeid())).or.(.not.lusesub)) then
c         call deletefile(d_f1m(iref))
c        endif
c       enddo
c      endif
c#else
      do iref=1,nref
c        if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
c     1 +innodes+ga_nodeid())).or.(.not.lusesub)) then
         call deletefile(d_f1m(iref))
c        endif
      enddo
c#endif

#ifdef MRCC_LOCAL_FOCK
ckbn @todo : localize fock
       if(lusesub) then
        do iref=1,nref
         if((int_mb(k_refafi+iref-1).eq.int_mb(k_innodes
     +      +innodes+ga_nodeid())).or.(.not.lusesub)) then
             if(.not.MA_POP_STACK(l_f1_local)) 
     +        call errquit('f1_local',1,MA_ERR)
         endif
        enddo
       endif
#endif

      if (.not.ma_pop_stack(l_goffset))
     1    call errquit("tce_mrcc_energy: MA problem",611,MA_ERR)

      if (.not.ma_pop_stack(l_refafi))
     1    call errquit("tce_mrcc_energy: MA problem",612,MA_ERR)

      if (.not.ma_pop_stack(l_sghandles))
     1    call errquit("tce_mrcc_energy: MA problem",613,MA_ERR)

      if (.not.ma_pop_stack(l_sgsizes))
     1    call errquit("tce_mrcc_energy: MA problem",614,MA_ERR)

      if (.not.ma_pop_stack(l_innodes))
     1    call errquit("tce_mrcc_energy: MA problem",615,MA_ERR)
 
      do i=1,nref

      if (.not.ma_pop_stack(l_f1_offsetm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",44,MA_ERR)

      enddo

      do i=1,nref

      if (.not.ma_pop_stack(l_alpham(nref-i+1)))
     1  call errquit("tce_energy: MA problem",4,MA_ERR)
      if (.not.ma_pop_stack(l_offsetm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",5,MA_ERR)
      if (.not.ma_pop_stack(l_rangem(nref-i+1)))
     1  call errquit("tce_energy: MA problem",6,MA_ERR)
      if (.not.ma_pop_stack(l_symm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",7,MA_ERR)
      if (.not.ma_pop_stack(l_spinm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",8,MA_ERR)

      enddo

      do i=1,ipol

      if (.not.ma_pop_stack(l_irs_tile(ipol-i+1)))
     1  call errquit("tce_energy: MA problem",9,MA_ERR)
      if (.not.ma_pop_stack(l_evl_tile(ipol-i+1)))
     1  call errquit("tce_energy: MA problem",9,MA_ERR)

      enddo

      do i=1,nref

      if (.not.ma_pop_stack(l_active_tmpm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",9,MA_ERR)
      if (.not.ma_pop_stack(l_range_tmpm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",9,MA_ERR)
      if (.not.ma_pop_stack(l_sym_tmpm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",10,MA_ERR)
      if (.not.ma_pop_stack(l_spin_tmpm(nref-i+1)))
     1  call errquit("tce_energy: MA problem",11,MA_ERR)

      enddo

      do i=1,nref
       j = nref - i + 1
      if(.not.ma_pop_stack(l_isactive(j)))
     1  call errquit("tce_mrcc_tile: ma problem",95,MA_ERR)
      if(.not.ma_pop_stack(l_mo_indexm(j)))
     1  call errquit("tce_mrcc_tile: ma problem",85,MA_ERR)
      if (.not.ma_pop_stack(l_evl_sortedm(j)))
     1  call errquit("tce_mrcc_tile: MA problem",9,MA_ERR)
      if (.not.ma_pop_stack(l_irs_sortedm(j)))
     1  call errquit("tce_mrcc_tile: MA problem",10,MA_ERR)
      if (.not.ma_pop_stack(l_spin_sortedm(j)))
     1  call errquit("tce_mrcc_tile: MA problem",11,MA_ERR)
      if (.not.ma_pop_stack(l_movecs_sortedm(j)))
     1  call errquit("tce_mrcc_tile: MA problem",12,MA_ERR)

      end do

c PNNL
c      do iref=1,nref
c         call deletefile(d_f1m(iref))
c      enddo
c      if(.not.diskint) then
c       do iref=1,nref
c          call deletefile(d_v2m(iref))
c       enddo
c      end if
c
c      do iref=1,nref
c         if(needt1) call deletefile(d_t1m(iref))
c         if(needt2) call deletefile(d_t2m(iref))
c         if(needt3) call deletefile(d_t3m(iref))
c      enddo

c
c     =========
c     Terminate
c     =========
c
ckbn
      call ga_sync()
      call tce_tidy(rtdb)
      call util_print_pop

c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup end
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup end
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup end
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup end
c cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup cleanup end

      return

 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9100 format(1x,i4,2f18.13,2f8.1,2f18.13)
 9123 format(1x,'4-index algorithm nr.',i4,1x,'is used')
 9220 format(1x,'MRCC cpu & wall time / sec',2f15.1)
 9222 format(1x,'Residues and T update / sec',2f15.1)

 9080 format(1x,A,'(Hartrees) = ',f25.15)

#else
      CALL ERRQUIT('MRCC Methods not compiled (tce_mrcc_energy.F) 
     + Compile  with "export MRCC_METHODS=y"',0,CAPMIS_ERR)
#endif

      end

c $Id$
