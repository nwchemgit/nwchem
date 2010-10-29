      subroutine intb_chk_sum(rtdb,basis,print_ints)
* $Id$
* does checksum for 2 electron eri blocking integral interface only
*
* geom and basis data structures initiated and destroyed outside 
* this routine
      implicit none
c::includes
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "errquit.fh"
*
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basis      ! [input] basis set handle 
      logical print_ints ! [input] also print integral and label values to luout
*
c::local
      integer h_eri, i_eri, h_scr, i_scr         ! handle and index for eri and scratch arrays
      integer i_ilab, h_ilab, i_jlab, h_jlab     ! handles and index for labels buffers
      integer i_klab, h_klab, i_llab, h_llab     ! handles and index for labels buffers
      integer i_iqtet, h_iqtet, i_jqtet, h_jqtet ! handles and index for quartet buffers
      integer i_kqtet, h_kqtet, i_lqtet, h_lqtet ! handles and index for quartet buffers
      integer i_iat, h_iat, i_jat, h_jat         ! handles and index for atom cache
      integer i_kat, h_kat, i_lat, h_lat         ! handles and index for atom cache
      integer h_stats, i_stats                   ! handle and index for stats of checksum
      integer h_q4, k_q4                         ! handle and index for q4 symmetry factors
      integer atom_blk_sz                        ! number of atom quartets in atom cache
c
      character*255 name69
      logical print69  ! print integrals to unit 69
      logical dumm69   ! used to read print69 from rtdb
      logical all_at_once ! logical for single block pass
      logical status      ! local logical value
c
      integer geom   ! geometry handle used to load basis set
      integer nat    ! number of atoms
      integer ncont  ! total number of contractions
      integer nbf    ! total number of basis functions
      integer maxg   ! maximum memory for eri quartet
      integer memscr ! maximum scratch for eri quartet (nwchem)
      integer max_shell_q ! max number of shell quartets
      integer max_atom_q  ! max number of atom quartets
      integer max_ints    ! max number of integrals 
      integer memb        ! memory used for buffers of integrals and labels  
      integer int_tmp           ! local integer temp variable
      integer neri
      integer num_quart
      integer num_ints_total
      double precision dbl_tmp  ! local dp temp variable
      double precision zerotol
c
      logical intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      logical debug    ! internal debug print enabled if .true.
      data debug /.true./
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c set zero tolerance
      zerotol = 1.0d-08
c
c
c set print69 flag
c
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:b2e4cprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.b2e4c',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c get information from geom/basis objects
c      
*.. get geometry handle
      if (.not.bas_geom(basis,geom)) call errquit
     &    ('intb_chk_sum:bas_geom failed',911, BASIS_ERR)
*.. get number of atoms
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('intb_chk_sum:geom_ncent failed',911, GEOM_ERR)
*.. get number of contractions
      if (.not.bas_numcont(basis,ncont)) call errquit
     &    ('intb_chk_sum:bas_numcont failed ',911, BASIS_ERR)
*.. get number of basis functions
      if (.not.bas_numbf(basis,nbf)) call errquit
     &    ('intb_chk_sum:bas_numbf failed ',911, BASIS_ERR)
*.. maximum number of integrals N**4/8
      max_ints = nbf*nbf*nbf*nbf/8
      int_tmp = nat*(nat+1)/2
      max_atom_q = int_tmp*(int_tmp+1)/2
      int_tmp = ncont*(ncont+1)/2
      max_shell_q = int_tmp*(int_tmp+1)/2
      write(luout,*)' maximum number of atom  quartets: ',max_atom_q
      write(luout,*)' maximum number of shell quartets: ',max_shell_q
      write(luout,*)' maximum number of integrals     : ',max_ints
*.. check memory to see if we can do it all in one pass
      memb = ma_inquire_stack(MT_DBL)
      call int_mem_2e4c(maxg,memscr)
*
      int_tmp = 2*max_ints                    ! eri buffer
      int_tmp = int_tmp + ((4*max_ints)/2)    ! labels 
      int_tmp = int_tmp + ((max_shell_q*4)/2) ! quartet buffers
      int_tmp = int_tmp + memscr + 15         ! scratch and stats
      if (debug) then
        write(luout,*)' memory on stack : ',memb
        write(luout,*)' memory needed   : ',int_tmp
      endif
      all_at_once = memb.ge.int_tmp
      if (rtdb_get(rtdb,'intcsum:b2e4c:all_at_once',
     &    MT_LOG,1,status)) then
        write(luout,*)' all_at_once forced to be ',status
        all_at_once = status
      endif
      if (.not.all_at_once) then
        atom_blk_sz = 20
        if(rtdb_get(rtdb,'intcsum:b2e4c:atom_size',
     &      MT_INT,1,int_tmp)) then
          write(luout,*)' atom cahce block size set to ',int_tmp
          atom_blk_sz = int_tmp
        endif
        max_shell_q = atom_blk_sz*50
        if(rtdb_get(rtdb,'intcsum:b2e4c:quartet_size',
     &      MT_INT,1,int_tmp)) then
          write(luout,*)' quartet cache buffer size set to ',int_tmp
          max_shell_q = int_tmp
        endif
        memb = max_shell_q*81
        if(rtdb_get(rtdb,'intcsum:b2e4c:buffer_size',
     &      MT_INT,1,int_tmp)) then
          write(luout,*)' eri/label buffer size set to ',int_tmp
          memb = int_tmp
        endif
        memscr = memscr + memscr/10
        if(rtdb_get(rtdb,'intcsum:b2e4c:scratch_size',
     &      MT_INT,1,int_tmp)) then
          write(luout,*)' scratch array size set to ',int_tmp
          memscr = int_tmp
        endif
        if (debug) then
          write(luout,*)' atom cahce block size is     :',atom_blk_sz
          write(luout,*)' eri/label buffer size is     :',memb
          write(luout,*)' quartet cache block size is  :',max_shell_q
          write(luout,*)' scratch array size is        :',memscr
        endif
      endif
*
      if (all_at_once) then
c
c do it all at once
c
        write(luout,*)'intb_chk_sum: all at once computation done '
*
        memb = max_ints + max_ints/10  ! add 10 percent
        status = ma_push_get(MT_DBL,memb,'eri buffer',h_eri,i_eri)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'i label buffer',h_ilab,i_ilab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'j label buffer',h_jlab,i_jlab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'k label buffer',h_klab,i_klab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'l label buffer',h_llab,i_llab)
        status = status.and.
     &      ma_push_get(MT_DBL,memscr,'scratch buffer',h_scr,i_scr)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'i cont buffer',
     &      h_iqtet,i_iqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'j cont buffer',
     &      h_jqtet,i_jqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'k cont buffer',
     &      h_kqtet,i_kqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'l cont buffer',
     &      h_lqtet,i_lqtet)
        status = status.and.
     &      ma_push_get(MT_DBL,max_shell_q,'q4 cont buffer',
     &      h_q4,k_q4)
        status = status.and.
     &      ma_push_get(MT_DBL,(5*3),'summary stats',h_stats,i_stats)
        if(.not.status) call errquit
     &      ('intb_chk_sum: ma_push_get failed somewhere',911,
     &       MA_ERR)
* zero memory
        call ifill(memb,0,int_mb(i_ilab),1)
        call ifill(memb,0,int_mb(i_jlab),1)
        call ifill(memb,0,int_mb(i_klab),1)
        call ifill(memb,0,int_mb(i_llab),1)
        call ifill(max_shell_q,0,int_mb(i_iqtet),1)
        call ifill(max_shell_q,0,int_mb(i_jqtet),1)
        call ifill(max_shell_q,0,int_mb(i_kqtet),1)
        call ifill(max_shell_q,0,int_mb(i_lqtet),1)
        call dfill(max_shell_q,0.0d00,dbl_mb(k_q4),1)
        call dfill((5*3),0.0d00,dbl_mb(i_stats),1)
        call dfill(memb,0.0d00,dbl_mb(i_eri),1)
        call dfill(memscr,0.0d00,dbl_mb(i_scr),1)
*
        call intb_all_at_once_gen(
     &      geom,basis,nat,max_shell_q,
     &      int_mb(i_iqtet),int_mb(i_jqtet),
     &      int_mb(i_kqtet),int_mb(i_lqtet),num_quart,
     &      dbl_mb(k_q4),debug)

        num_ints_total = 0
        status = intb_init4c(basis,int_mb(i_iqtet),int_mb(i_jqtet),
     &      basis,int_mb(i_kqtet),int_mb(i_lqtet),num_quart,
     $      dbl_mb(k_q4),.true.,
     &      memscr,dbl_mb(i_scr),memb,dbl_tmp)
        if (.not.status) call errquit
     &      ('intb_chk_sum: init blocking error',911, INT_ERR)
00001   continue
        status =
     &      intb_2e4c(
     &      basis,int_mb(i_iqtet),int_mb(i_jqtet),
     &      basis,int_mb(i_kqtet),int_mb(i_lqtet),num_quart,
     &      dbl_mb(k_q4),.true.,
     &      zerotol,.false.,
     &      int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab),
     &      dbl_mb(i_eri),memb,neri,memscr,dbl_mb(i_scr))
        call intb_simple_canon(int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab),
     &      dbl_mb(i_eri), neri, zerotol)
*        write(6,*)' stsum_inve called 1'
        call stsum_inve(neri,
     &      dbl_mb(i_eri),
     &      int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab))
        call stsum2(neri,num_ints_total,
     &      dbl_mb(i_eri),
     &      int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab),
     &      dbl_mb(i_stats),
     &      print_ints,print69,'intb: 2e4c')
        if (status) then
          if (debug) write(luout,*)
     &        ' intb_chk_sum: more integrals goto 00001'
          goto 00001
        endif
c
c print stats
c
        call prt_stat_sum('all at once blocked two electron',
     &      dbl_mb(i_stats),num_ints_total)
        call stsum_inve(-1,
     &      dbl_mb(i_eri),
     &      int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab))
      else
c
c multiple passes required
c
        write(luout,*)'intb_chk_sum: multi-pass computation done '
        status = ma_push_get(MT_DBL,memb,'eri buffer',h_eri,i_eri)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'i label buffer',h_ilab,i_ilab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'j label buffer',h_jlab,i_jlab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'k label buffer',h_klab,i_klab)
        status = status.and.
     &      ma_push_get(MT_INT,memb,'l label buffer',h_llab,i_llab)
        status = status.and.
     &      ma_push_get(MT_DBL,memscr,'scratch buffer',h_scr,i_scr)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'i cont buffer',
     &      h_iqtet,i_iqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'j cont buffer',
     &      h_jqtet,i_jqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'k cont buffer',
     &      h_kqtet,i_kqtet)
        status = status.and.
     &      ma_push_get(MT_INT,max_shell_q,'l cont buffer',
     &      h_lqtet,i_lqtet)
        status = status.and.
     &      ma_push_get(MT_DBL,max_shell_q,'q4 cont buffer',
     &      h_q4,k_q4)
        status = status.and.
     &      ma_push_get(MT_DBL,(5*3),'summary stats',h_stats,i_stats)
        status = status.and.
     &      ma_push_get(MT_INT,atom_blk_sz,'i atom buffer',
     &      h_iat,i_iat)
        status = status.and.
     &      ma_push_get(MT_INT,atom_blk_sz,'j atom buffer',
     &      h_jat,i_jat)
        status = status.and.
     &      ma_push_get(MT_INT,atom_blk_sz,'k atom buffer',
     &      h_kat,i_kat)
        status = status.and.
     &      ma_push_get(MT_INT,atom_blk_sz,'l atom buffer',
     &      h_lat,i_lat)
        if(.not.status) call errquit
     &      ('intb_chk_sum: ma_push_get failed somewhere',911,
     &       MA_ERR)
* zero memory
        call ifill(memb,0,int_mb(i_ilab),1)
        call ifill(memb,0,int_mb(i_jlab),1)
        call ifill(memb,0,int_mb(i_klab),1)
        call ifill(memb,0,int_mb(i_llab),1)
        call ifill(max_shell_q,0,int_mb(i_iqtet),1)
        call ifill(max_shell_q,0,int_mb(i_jqtet),1)
        call ifill(max_shell_q,0,int_mb(i_kqtet),1)
        call ifill(max_shell_q,0,int_mb(i_lqtet),1)
        call ifill(atom_blk_sz,0,int_mb(i_iat),1)
        call ifill(atom_blk_sz,0,int_mb(i_jat),1)
        call ifill(atom_blk_sz,0,int_mb(i_kat),1)
        call ifill(atom_blk_sz,0,int_mb(i_lat),1)
        call dfill(max_shell_q,0.0d00,dbl_mb(k_q4),1)
        call dfill((5*3),0.0d00,dbl_mb(i_stats),1)
        call dfill(memb,0.0d00,dbl_mb(i_eri),1)
        call dfill(memscr,0.0d00,dbl_mb(i_scr),1)
*
        call intb_atom_blocks_chksum(
     &      memb,memscr,max_shell_q,atom_blk_sz,nat,
     &      geom,basis,
     &      dbl_mb(i_eri),int_mb(i_ilab),int_mb(i_jlab),
     &      int_mb(i_klab),int_mb(i_llab),
     &      int_mb(i_iqtet),int_mb(i_jqtet),
     &      int_mb(i_kqtet),int_mb(i_lqtet),
     &      dbl_mb(k_q4),
     &      int_mb(i_iat),int_mb(i_jat),
     &      int_mb(i_kat),int_mb(i_lat),
     &      dbl_mb(i_scr),dbl_mb(i_stats),zerotol,
     &      print_ints,print69,debug)
      endif
*free memory in order
      status = .true.
      if (.not.all_at_once) then
        status = status .and. ma_pop_stack(h_lat)
        status = status .and. ma_pop_stack(h_kat)
        status = status .and. ma_pop_stack(h_jat)
        status = status .and. ma_pop_stack(h_iat)
      endif
      status = status .and. ma_pop_stack(h_stats)
      status = status .and. ma_pop_stack(h_q4)
      status = status .and. ma_pop_stack(h_lqtet)
      status = status .and. ma_pop_stack(h_kqtet)
      status = status .and. ma_pop_stack(h_jqtet)
      status = status .and. ma_pop_stack(h_iqtet)
      status = status .and. ma_pop_stack(h_scr)
      status = status .and. ma_pop_stack(h_llab)
      status = status .and. ma_pop_stack(h_klab)
      status = status .and. ma_pop_stack(h_jlab)
      status = status .and. ma_pop_stack(h_ilab)
      status = status .and. ma_pop_stack(h_eri)
      if (.not.status) call errquit
     &    ('intb_chk_sum: ma_pop_stack failed',911,
     &       MA_ERR)
* close file      
      close(unit=69,status='keep')
      call ga_sync()
      end
      subroutine intb_atom_blocks_chksum(
     &    leri,memscr,max_q,a_blk,nat,
     &    geom,basis,eri,ilb,jlb,klb,llb,
     &    iq,jq,kq,lq,q4,
     &    iiat,jjat,kkat,llat,scr,stats,zerotol,
     &    printem,print69,debug)
      implicit none
c::includes
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical intb_init4c, intb_2e4c
      external intb_init4c, intb_2e4c
c::passed
      integer leri,memscr,max_q,a_blk,nat,basis,geom
      double precision eri(leri)
      integer ilb(leri),jlb(leri),klb(leri),llb(leri)
      integer iq(max_q),jq(max_q),kq(max_q),lq(max_q)
      integer iiat(a_blk),jjat(a_blk),kkat(a_blk),llat(a_blk)
      double precision q4(max_q)
      double precision scr(memscr),stats(5,3), zerotol
      logical printem,print69,debug
c::local      
      integer num_ints_total, num_q, iat, kat, jat, lat, lathi
      integer neri, nat_cache, num_q_max
      double precision dbl_tmp
      logical status, more_integrals, last_quartet
#define ATOM_BLOCKED_SCF
#if defined(ATOM_BLOCKED_SCF)
      character*30 my_date
      integer kbk, ibk, jbk 
      integer kathi, katlo, iathi, iatlo, jathi, jatlo
#endif
c
      num_ints_total = 0
      num_q_max = -1
#if defined(ATOM_BLOCKED_SCF)
c---------------------------------------------|
*-- define block chunk sizes for atom loops --*
c---------------------------------------------|
      kbk = 2          ! katom block size     |
      ibk = 2          ! iatom block size     |
      jbk = 2          ! jatom block size     |
c---------------------------------------------|

      nat_cache = 0
      do kathi = nat, 1, -kbk
        katlo = max(1,kathi-kbk+1)
        call util_date(my_date)
        write(6,'(a,i4,a,i4,a,a)')
     &      ' doing atom block [',kathi,'->',katlo,'] ',
     &      my_date(1:26)
        call util_flush(6)
        do iathi = nat, 1, -ibk
          iatlo = max(1,iathi-ibk+1)
          do jathi = nat, 1, -jbk
            jatlo = max(1,jathi-jbk+1)
            do kat = kathi,katlo,-1
              do iat = iathi, iatlo, -1
                do jat = jathi,jatlo,-1
                  if (iat.ge.jat) then
                    if (iat.ge.kat) then
                      lathi = kat
                      if (iat.eq.kat) lathi = jat
                      do lat = 1,lathi
                        nat_cache = nat_cache + 1
                        iiat(nat_cache) = iat
                        jjat(nat_cache) = jat
                        kkat(nat_cache) = kat
                        llat(nat_cache) = lat
                        last_quartet = (kat.eq.1)
                        last_quartet = last_quartet.and.(iat.eq.1)
                        last_quartet = last_quartet.and.(jat.eq.1)
                        last_quartet = last_quartet.and.(lat.eq.lathi)
#else
      nat_cache = 0
      do kat = nat, 1, -1
        do iat = nat, kat, -1
          do jat = iat, 1, -1
            lathi = kat
            if (iat.eq.kat) lathi = jat
            do lat = 1,lathi
              nat_cache = nat_cache + 1
              iiat(nat_cache) = iat
              jjat(nat_cache) = jat
              kkat(nat_cache) = kat
              llat(nat_cache) = lat
              last_quartet = (kat.eq.1)
              last_quartet = last_quartet.and.(iat.eq.kat)
              last_quartet = last_quartet.and.(jat.eq.1)
              last_quartet = last_quartet.and.(lat.eq.lathi)
#endif
                        if ((.not.last_quartet).and.
     &                      nat_cache.lt.a_blk) goto 00001 ! skip until atom block filled
c---------------------------------------------------------------------
c same part for both loop structures - start
c---------------------------------------------------------------------
*.. do group of atoms work              
              num_q = 0
              call intb_quart_from_atoms(
     &            geom, basis,iiat,jjat,kkat,llat,
     &            nat_cache,max_q,iq,jq,kq,lq,num_q,q4,debug)
*rak:              write(6,*)'iiat  =',iiat
*rak:              write(6,*)'jjat  =',jjat
*rak:              write(6,*)'kkat  =',kkat
*rak:              write(6,*)'llat  =',llat
*rak:              num_q_max = max(num_q,num_q_max)
*rak:              write(6,*)'num_q =',num_q, num_q_max
*rak:              write(6,*)
*rak:     &          ' max_q, leri, memscr, a_blk, nat, neri 1',
*rak:     &            max_q, leri, memscr, a_blk, nat, neri
*rak:              status = ma_verify_allocator_stuff()
              status = intb_init4c(basis,iq,jq,basis,kq,lq,num_q,
     &            q4, .true.,memscr,scr,leri,dbl_tmp)
*rak:              write(6,*)
*rak:     &          ' max_q, leri, memscr, a_blk, nat, neri 2',
*rak:     &            max_q, leri, memscr, a_blk, nat, neri
*rak:              status = ma_verify_allocator_stuff()
              if (.not.status) call errquit
     &            ('intb_atom_blocks_chksum: init blocking error',911,
     &       INT_ERR)
00002         continue
              more_integrals = intb_2e4c(
     &            basis, iq, jq, basis, kq, lq, num_q, q4, .true.,
     &            zerotol,.false.,
     &            ilb,jlb,klb,llb,eri,leri,neri,memscr,scr)
*rak:              write(6,*)
*rak:     &          ' max_q, leri, memscr, a_blk, nat, neri 3',
*rak:     &            max_q, leri, memscr, a_blk, nat, neri
*rak:              status = ma_verify_allocator_stuff()
              call intb_simple_canon(ilb,jlb,klb,llb,eri,neri,zerotol)
*rak:              write(6,*)
*rak:     &          ' max_q, leri, memscr, a_blk, nat, neri 4',
*rak:     &            max_q, leri, memscr, a_blk, nat, neri
*rak:              status = ma_verify_allocator_stuff()
              call stsum2(neri,num_ints_total,eri,ilb,jlb,klb,llb,
     &            stats,printem,print69,'intb: 2e4c')
*              write(6,*)' stsum_inve called 2'
              call stsum_inve(neri,eri,ilb,jlb,klb,llb)
*rak:              write(6,*)
*rak:     &          ' max_q, leri, memscr, a_blk, nat, neri 5',
*rak:     &            max_q, leri, memscr, a_blk, nat, neri
*rak:              status = ma_verify_allocator_stuff()
              if (more_integrals) then
                if(debug) write(luout,*)
     &              'intb_atom_blocks_chksum: more ints goto 00002'
                goto 00002
              endif
              nat_cache = 0
00001         continue
c---------------------------------------------------------------------
c same part for both loop structures - end
c---------------------------------------------------------------------
#if defined(ATOM_BLOCKED_SCF)
                      enddo
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
#else
            enddo
          enddo
        enddo
      enddo
#endif
*
      call prt_stat_sum('multi-passed blocked two electron',
     &    stats,num_ints_total)
      call stsum_inve(-1,eri,ilb,jlb,klb,llb)
*
      end
      subroutine intb_quart_from_atoms(
     &    geom,basis,iiat,jjat,kkat,llat,num_at,
     &    max_q,iq,jq,kq,lq,num_q,q4,debug)
      implicit none
c
c generates all shell quartets for given atom quartets
c
c::includes
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#if defined(SUM_WITH_SYM)
#include "sym.fh"
#endif
c::passed
      integer geom
      integer basis
      integer num_at
      integer iiat(num_at),jjat(num_at),kkat(num_at),llat(num_at)
      integer max_q
      integer iq(max_q),jq(max_q),kq(max_q),lq(max_q)
      double precision q4(max_q)
      integer num_q
      logical debug
c::local
      integer count
      integer iat,jat,kat,lat,iatom
      integer islo,jslo,kslo,lslo
      integer ishi,jshi,kshi,lshi
      integer ish,jsh,ksh,lsh
      integer kstop,jstop,lstop
      double precision q4_atom
      logical oaij,oakl,oaikjl,status
c
*.. reset debug to false
      if (debug) debug = .false.
c
      count = 0
      do iatom = 1,num_at
        iat = iiat(iatom)
        jat = jjat(iatom)
        kat = kkat(iatom)
        lat = llat(iatom)
#if defined(SUM_WITH_SYM)
        if (.not.sym_atom_quartet(geom,iat,jat,kat,lat,q4_atom))
     &      goto 00001
#else
        q4_atom  = 1.0d00  
#endif
        status = .true.
        status = status .and. bas_ce2cnr(basis, iat, islo, ishi)
        status = status .and. bas_ce2cnr(basis, jat, jslo, jshi)
        status = status .and. bas_ce2cnr(basis, kat, kslo, kshi)
        status = status .and. bas_ce2cnr(basis, lat, lslo, lshi)
        if (.not.status) call errquit
     &      ('intb_quart_from_atoms: bas_ce2cnr failed',911, BASIS_ERR)
        oaij = iat.eq.jat
        oakl = kat.eq.lat
        oaikjl = iat.eq.kat .and. jat.eq.lat
        do ish = islo,ishi
          jstop = jshi
          if (oaij) jstop = ish
          do jsh = jslo,jstop
            kstop = kshi
            if (oaikjl) kstop = ish
            do ksh = kslo, kstop
              lstop = lshi
              if (oakl) lstop = ksh
              if (oaikjl.and. (ish.eq.ksh)) lstop = jsh
              do lsh = lslo, lstop
                if ((count+1).le.max_q) then
                  count = count + 1
                  iq(count) = ish
                  jq(count) = jsh
                  kq(count) = ksh
                  lq(count) = lsh
                  q4(count) = q4_atom
                  if (q4_atom.gt.1.0d00) then
                    write(6,*)' 1 ',ish,jsh,ksh,lsh,q4_atom
                  endif
                  if (debug) 
     &                write(luout,10000)count,
     &                iat,jat,kat,lat,ish,jsh,ksh,lsh
                else
                  write(luout,*)' count = ',(count+1)
                  write(luout,*)' max_q = ',max_q
                  call errquit
     &                ('intb_quart_from_atoms: failed',911, INT_ERR)
                endif
              enddo
            enddo
          enddo
        enddo
#if defined(SUM_WITH_SYM)
00001   continue
#endif
      enddo
*
      num_q = count
*
10000 format(1x,'quartet ',i6,' (atoms|<',i3,i3,'|',i3,i3,
     &    '>) (shells|<',i3,i3,'|',i3,i3,'>)')
      end
      subroutine intb_all_at_once_gen(geom,basis,nat,
     &    max_q,iq,jq,kq,lq,num_q,q4,debug)
c
c generate all quartets ddscf style
c
      implicit none
c::includes
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#if defined(SUM_WITH_SYM)
#include "sym.fh"
#endif
c::passed
      integer geom
      integer basis
      integer nat
      integer max_q
      integer iq(max_q)
      integer jq(max_q)
      integer kq(max_q)
      integer lq(max_q)
      integer num_q
      double precision q4(max_q)
      logical debug
c::local      
      double precision q4_atom
      integer count
      integer iat,jat,kat,lat
      integer islo,jslo,kslo,lslo
      integer ishi,jshi,kshi,lshi
      integer ish,jsh,ksh,lsh
      integer kstop,jstop,lstop
      integer lathi 
      logical oaij,oakl,oaikjl,status
c
*   reset debug to off
      if (debug) debug = .false.
c
      count = 0
      do kat = nat, 1, -1
        do iat = nat, kat, -1
          do jat = iat, 1, -1
            lathi = kat
            if (iat.eq.kat) lathi = jat
            do lat = 1,lathi
              status = .true.
              status = status .and. bas_ce2cnr(basis, iat, islo, ishi)
              status = status .and. bas_ce2cnr(basis, jat, jslo, jshi)
              status = status .and. bas_ce2cnr(basis, kat, kslo, kshi)
              status = status .and. bas_ce2cnr(basis, lat, lslo, lshi)
              if (.not.status) call errquit
     &            ('intb_all_at_once_gen: bas_ce2cnr failed',911,
     &       BASIS_ERR)
              oaij = iat.eq.jat
              oakl = kat.eq.lat
              oaikjl = iat.eq.kat .and. jat.eq.lat
#if defined(SUM_WITH_SYM)
              if (.not.sym_atom_quartet(geom,iat,jat,kat,lat,q4_atom))
     &            goto 00001
#else
              q4_atom = 1.0d00
#endif
              do ish = islo,ishi
                jstop = jshi
                if (oaij) jstop = ish
                do jsh = jslo,jstop
                  kstop = kshi
                  if (oaikjl) kstop = ish
                  do ksh = kslo, kstop
                    lstop = lshi
                    if (oakl) lstop = ksh
                    if (oaikjl.and. (ish.eq.ksh)) lstop = jsh
                    do lsh = lslo, lstop
                      if ((count+1).le.max_q) then
                        count = count + 1
                        iq(count) = ish
                        jq(count) = jsh
                        kq(count) = ksh
                        lq(count) = lsh
                        q4(count) = q4_atom
                        if (q4_atom.gt.1.0d00) then
                          write(6,*)' 1 ',ish,jsh,ksh,lsh,q4_atom
                        endif
                        if (debug) 
     &                      write(luout,10000)count,
     &                      iat,jat,kat,lat,ish,jsh,ksh,lsh
                      else
                        write(luout,*)' count = ',(count+1)
                        write(luout,*)' max_q = ',max_q
                        call errquit
     &                      ('intb_all_at_once_gen: failed',911,
     &       INT_ERR)
                      endif
                      
                    enddo
                  enddo
                enddo
              enddo
00001         continue
c... end of atom loops
            enddo
          enddo
        enddo
      enddo
c
      num_q = count
c
10000 format(1x,'quartet ',i6,' (atoms|<',i3,i3,'|',i3,i3,
     &    '>) (shells|<',i3,i3,'|',i3,i3,'>)')
c
      end
      subroutine intb_simple_canon(ilb,jlb,klb,llb,eri,neri,zerotol)
      implicit none
c
c canonicalization routine and zero removal routine
c
      integer neri               ! number of integrals
      integer ilb(neri)          ! i labels
      integer jlb(neri)          ! j labels
      integer klb(neri)          ! k labels
      integer llb(neri)          ! l labels
      double precision eri(neri) ! two electron integral values
      double precision zerotol   ! throw away integrals with this magnitude or less
c
      logical oij,okl,oijkl
      integer count
      integer ii,jj,kk,ll
      integer i,j,isym2
      isym2(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)          
*
      if (neri.gt.1) return
      count = 0
      do 00100 i=1,neri
        if (abs(eri(i)).gt.zerotol) then
          ii = ilb(i)
          jj = jlb(i)
          kk = klb(i)
          ll = llb(i)
          oij   = ii.ge.jj
          okl   = kk.ge.ll
          oijkl = (isym2(ii,jj)).ge.(isym2(kk,ll))
          if (oij.and.okl.and.oijkl) then
            count = count + 1
            if (count.lt.i) then
              ilb(count) = ii
              jlb(count) = jj
              klb(count) = kk
              llb(count) = ll
              eri(count) = eri(i)
            endif
          endif
        endif
00100 continue
*
      neri = count
*
      end

*      integer isym2,isym4                                 
*      ISYM2(I,J)=
*      ISYM4(I,J,K,L)=MAX(ISYM2(I,J),ISYM2(K,L))*          
*     &               (MAX(ISYM2(I,J),ISYM2(K,L))-1)/2+    
*     &               MIN(ISYM2(I,J),ISYM2(K,L))           
