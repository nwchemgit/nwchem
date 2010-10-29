      logical function raktest(rtdb)
      implicit none
#include "errquit.fh"
c $Id$
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "stdio.fh"
c::functions
      logical task_hondo_deriv_check
      external task_hondo_deriv_check
      logical task_ecp_deriv_check, task_ecp_print_integrals
      external task_ecp_deriv_check, task_ecp_print_integrals
      logical task_computeSld, task_printsoints, task_pderiv
      external task_computeSld, task_printsoints, task_pderiv
      logical task_dddd, task_accy, raktask_intdd, raktask_ecppe
      external task_dddd, task_accy, raktask_intdd, raktask_ecppe
      logical raktask_intdd_3c
      external raktask_intdd_3c
      logical raktask_geomcalc
      external raktask_geomcalc
      logical raktask_fullsc
      external raktask_fullsc
      logical rak_justrunvib
      external rak_justrunvib
c::passed
      integer rtdb          ! rtdb handle
c::local
      integer me
      integer raktask, rak_tmp
c
      call ga_sync()
      call ga_sync()
      raktask = 0
      if (rtdb_get(rtdb,'raktask',MT_INT,1,rak_tmp))
     &    raktask = rak_tmp
c
      call ga_sync()
      call ga_sync()
      me = ga_nodeid()
c
      if (raktask.eq.0) then    !...................................   0
        if (me.eq.0) then
          write(luout,*)' default raktest task '
          write(luout,*)' test semi empirical interface '
        endif
        call raktest_semi(rtdb)
        raktest = .true.
      else if (raktask.eq.1) then !.................................   1
        if (me.eq.0)write(luout,*)' raktest task 1 stepper test'
        call raktest_stpr(rtdb)
        raktest = .true.
      else if (raktask.eq.2) then !.................................   2
        if (me.eq.0)write(luout,*)' raktest task 2 check int_init'
        call raktest_init(rtdb)
        raktest = .true.
      else if (raktask.eq.3) then !.................................   3
        if (me.eq.0)write(luout,*)' raktest task 3 check intd_init'
        call raktest_initd(rtdb)
        raktest = .true.
      else if (raktask.eq.4) then !.................................   4
        if (me.eq.0)write(luout,*)' raktest check 3ctr nai'
        call raktest_3ctr(rtdb)
        raktest = .true.
      else if (raktask.eq.5) then !.................................   5
        if (me.eq.0)write(luout,*)' test of general contraction code '
        call raktest_gc(rtdb)
        raktest = .true.
      else if (raktask.eq.6) then !.................................   6
        if (me.eq.0)write(luout,*)' test of orbital printing code '
        call raktest_printorb(rtdb)
        raktest = .true.
      else if (raktask.eq.7) then !.................................   7
        if (me.eq.0)write(luout,*)' test of writing geom objects out '
        call raktest_geomwrt(rtdb)
        raktest = .true.
      else if (raktask.eq.8) then !.................................   8
        if (me.eq.0)write(luout,*)' test of spcart stuff '
        call raktest_spcart(rtdb)
        raktest = .true.
      else if (raktask.eq.9) then !.................................   9
        if (me.eq.0)write(luout,*)' test of spcart stuff all in one'
        call raktest_test9(rtdb)
        raktest = .true.
      else if (raktask.eq.10) then !................................. 10
        if (me.eq.0)write(luout,*)' test of ecp stuff '
        call raktest_ecp(rtdb)
        raktest = .true.
      else if (raktask.eq.11) then !................................. 11
        if (me.eq.0)write(luout,*)' bug in integrals test '
        call raktest_bug(rtdb)
        raktest = .true.
      else if (raktask.eq.12) then !................................. 12
        if (me.eq.0)write(luout,*)' geometry printing routines '
        call raktest_geomprt(rtdb)
        raktest = .true.
      else if (raktask.eq.13) then !................................. 13
        if (me.eq.0)write(luout,*)' test 3 center derivatives '
        call raktest_3cd(rtdb)
        raktest = .true.
      else if (raktask.eq.14) then !................................. 14
        if (me.eq.0)write(luout,*)' test derivative overlap '
        call raktest_ovd(rtdb)
        raktest = .true.
      else if (raktask.eq.15) then !................................. 15
        if (me.eq.0)write(luout,*)
     &      ' test blocking 2e derivative integral interface'
        call raktest_bd2e(rtdb)
        raktest = .true.
      else if (raktask.eq.16) then !................................. 16
        if (me.eq.0)write(luout,*)' raktest: disk test code '
        call raktest_diskspeed(rtdb)
        raktest = .true.
      else if (raktask.eq.17) then !................................. 17
        if (me.eq.0)write(luout,*)
     &      ' raktest: compare 2e integrals from nwchem & texas'
        call raktest_2ecompare(rtdb)
        raktest = .true.
      else if (raktask.eq.18) then !................................. 18
        if (me.eq.0)write(luout,*)
     &      ' raktest: compute overlap and linear dependence '
        raktest = task_computeSld(rtdb)
      else if (raktask.eq.19) then !................................. 19
        if (me.eq.0)write(luout,*)
     &      ' raktest: task print SO integrals '
        raktest = task_printSOints(rtdb)
      else if (raktask.eq.20) then !................................. 20
        if (me.eq.0)write(luout,*)
     &      ' raktest: task test periodic deriv '
        raktest = task_pderiv(rtdb)
      else if (raktask.eq.21) then !................................. 21
        if (me.eq.0)write(luout,*)
     &      ' raktest: task test dddd bug'
        raktest = task_dddd(rtdb)
      else if (raktask.eq.22) then !................................. 22
        if (me.eq.0)write(luout,*)
     &      ' raktest: accuracy test for ints/grads '
        raktest = task_accy(rtdb)
      else if (raktask.eq.23) then !................................. 23
        if (me.eq.0)write(luout,*)
     &      ' raktest: second derivative code '        
        raktest = raktask_intdd(rtdb)
      else if (raktask.eq.24) then !................................. 24
        if (me.eq.0)write(luout,*)
     &      ' raktest: ecp PE code test'        
        raktest = raktask_ecppe(rtdb)
      else if (raktask.eq.25) then !................................. 25
        if (me.eq.0) write(luout,*)
     &      ' raktest: ecp deriv code test'        
        raktest = task_ecp_deriv_check(rtdb)
      else if (raktask.eq.26) then !................................. 26
        if (me.eq.0) write(luout,*)
     &      ' raktest: ecp print integrals'        
        raktest = task_ecp_print_integrals(rtdb)
      else if (raktask.eq.27) then !................................. 27
        if (me.eq.0) write(luout,*)
     &      ' raktest: hondo deriv code test'        
        raktest = task_hondo_deriv_check(rtdb)
      else if (raktask.eq.28) then !................................. 28
        raktest = raktask_geomcalc(rtdb)
      else if (raktask.eq.29) then !................................. 29
        raktest = raktask_fullsc(rtdb)
      else if (raktask.eq.30) then !................................. 30
        if (me.eq.0) write(luout,*)
     &      ' raktest: check 2e3c second derivatives'        
        raktest = raktask_intdd_3c(rtdb)
      else if (raktask.eq.31) then !................................. 31
        if (me.eq.0) write(luout,*)
     &      ' raktest: rerun vib'        
        raktest = rak_justrunvib(rtdb)
      else
        if (me.eq.0) then
          write(luout,*)' unknown raktask number :',raktask
          call errquit('raktest: fatal error',911, INPUT_ERR)
        endif
      endif
c
      end
      subroutine raktest_geomwrt(rtdb)
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "geom.fh"
c
      integer rtdb
c
      character*40 new_geom_name
      integer geom
      integer igeom
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911,
     &       RTDB_ERR)
c
      do igeom = 1,110
	new_geom_name = ' '
        if (abs(igeom) .lt. 10) then
           write(new_geom_name,'(''g-'',i1,''-step'')') 
     $          abs(igeom)
        else if (abs(igeom) .lt. 100) then
           write(new_geom_name,'(''g-'',i2,''-step'')') 
     $          abs(igeom)
        else if (abs(igeom) .lt. 1000) then
           write(new_geom_name,'(''g-'',i3,''-step'')') 
     $          abs(igeom)
        else
           write(new_geom_name,'(''g-'',i4,''-step'')') 
     $          abs(igeom)
        endif

        call sym_geom_project(geom, 1d-4)

        if (.not.geom_rtdb_store(rtdb,geom,new_geom_name))
     &      call errquit
     &      ('raktest_geomwrt: geom_rtdb_store (of copy) failed',911,
     &       RTDB_ERR)
        write(luout,*)' stored geometry ',new_geom_name
      enddo
      end
      subroutine raktest_printorb(rtdb)
      implicit none
c
#include "rtdb.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
c
      integer rtdb ! [input]
c
      
c
      end
      subroutine gen_bf_tag(basis,i_bf,bf_tag)
      implicit none
#include "bas.fh"
c
c generate a string that tells all about the basis function      
c structure
c  bf_tag(1:16)  = geom_tag() or bas_tag ! user symbol for basis function
c  bf_tag(17:17) = ' '                   ! blank  
c  bf_tag(18:18) = type                  ! l, s, p, d, ....
c  bf_tag(19:19) = ' '                   ! blank  
c  bf_tag(20:27) = xyz_tag()             ! xyz's for bf
c  
      integer basis
      integer i_bf
      character*27 bf_tag
      integer cont
      integer center
      character*1 ch_type(-1:5)
      data ch_type /'l','s','p','d','f','g','h'/
c
c map bf -> cn
c map bf -> ce -> tag
c map cn -> type
c map cn -> bfr -> ic 
      if (.not. bas_bf2cn(basis,i_bf,cont)) stop 'ceq'
      if (.not. bas_bf2ce(basis,i_bf,center)) stop 'ceq'
c      
      end
      subroutine int_xyz_tag(lval,ic,xyz_tag,l_tag)
      implicit none
#include "errquit.fh"
      integer lval    ! [input] l value
      integer ic      ! [input] cartesean component
      integer l_tag ! [input] length of xyz_tag character array
      character*(*) xyz_tag ! [output] left justified
c
      integer nxyz(3)
      character*1 pxyz(3)
      integer ixyz, i, j
      data pxyz /'x','y','z'/
      save pxyz
c
      if (lval.eq.0) then
        xyz_tag(1:3) = ' s '
        ixyz = 4
c      
      elseif (lval.eq.-1) then
        if (ic.eq.1) then
          xyz_tag(1:3) = ' s '
        elseif (ic.eq.2) then
          xyz_tag(1:3) = ' x '
        elseif (ic.eq.2) then
          xyz_tag(1:3) = ' y '
        elseif (ic.eq.2) then
          xyz_tag(1:3) = ' z '
        else
          call errquit('int_xyz_tag: error on lval=-1,ic=',ic,
     &       INPUT_ERR)
        endif
        ixyz = 4
      elseif (lval.gt.0) then
        call defNxyz(lval)
        call getNxyz(lval,ic,nxyz)
c
        ixyz = 1 
c     
        do i=1,3
          do j=1,nxyz(i)
            xyz_tag(ixyz:ixyz) = pxyz(i)
            ixyz = ixyz + 1
          enddo
        enddo
      else
        call errquit('int_xyz_tag: error on lval=',lval, INPUT_ERR)
      endif
      do i = ixyz, l_tag
        xyz_tag(i:i) = ' '
      enddo
      end
c...............................................................................
      subroutine raktest_spcart(rtdb)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer rtdb ! [input] rtdb handle
c
      integer geom, basis
      integer sp_basis, nshell_sp
      integer nbf, nbfsq, nbf_sp, nbf_chk, nshell
      integer max1e, max2e, mscr1, mscr2, m_scr, m_buf
      integer h_cart_s, h_sph_s, h_scr, h_buf, h_2bfr
      integer k_cart_s, k_sph_s, k_scr, k_buf, k_2bfr
      integer h_cart_s2, h_sph_s2
      integer k_cart_s2, k_sph_s2
      integer h_eri_1, h_eri_2, h_eri_sp_1, h_eri_sp_2
      integer k_eri_1, k_eri_2, k_eri_sp_1, k_eri_sp_2
      double precision norm_cart, norm_sph
      double precision ddot
      external ddot
      logical status
c
      logical int_normalize, int_norm_2c
      external int_normalize, int_norm_2c
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &      ('geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &      ('geom_rtdb_load failed',911, RTDB_ERR)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
c
      write(6,*)' geom/basis loaded'
c
      write(6,*)' raw basis '
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0, BASIS_ERR)
      write(6,*)' first normalize'
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0, BASIS_ERR)
      write(6,*)' second normalize'
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 2'
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0, BASIS_ERR)
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &      ('numbf failed',911, BASIS_ERR)
c
      nbfsq = nbf*nbf
      if (.not.ma_push_get(mt_dbl,nbfsq,'square cart overlap',
     &      h_cart_s, k_cart_s)) call errquit
     &      (' cart overlap ma failed ',911, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square spher overlap',
     &      h_sph_s, k_sph_s)) call errquit
     &      (' sph overlap ma failed ',911, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square cart overlap 2',
     &      h_cart_s2, k_cart_s2)) call errquit
     &      (' cart2 overlap ma failed ',911, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square spher overlap 2',
     &      h_sph_s2, k_sph_s2)) call errquit
     &      (' sph2 overlap ma failed ',911, MA_ERR)
c
      if (.not.bas_numcont(basis,nshell)) call errquit
     &      ('numcont error',911, BASIS_ERR)
c
      call int_init(rtdb,1,basis)
      call int_mem(max1e,max2e,mscr1,mscr2)
      m_buf = max(max1e*2,max2e*2)
      m_scr = max(mscr1*2,mscr2)
      m_buf = m_buf + (m_buf*110)/100
      m_scr = m_scr + (m_scr*110)/100
      
c
      if (.not.ma_push_get(mt_dbl,m_scr,'scr for 1e',h_scr,k_scr))
     &      call errquit('ma scr failed',911, MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,m_buf,'buf for 1e',h_buf,k_buf))
     &      call errquit('ma buf failed',911, MA_ERR)
c
      if (.not.ma_push_get(mt_int,2*nshell,'buf for sp cn2bfr',
     &      h_2bfr,k_2bfr))
     &      call errquit('ma buf failed',911, MA_ERR)
      call rak_tospbfr(basis,nshell,nbf_chk,nbf_sp,int_mb(k_2bfr))
c
      if (nbf_chk.ne.nbf) then
        write(6,*)' nbf not right ',nbf_chk, nbf
      endif
c
      write(6,*)' nbf    ',nbf
      write(6,*)' nbf_sp ',nbf_sp
c
      if (.not.ma_push_get(mt_dbl,(nbf*nbf*nbf*nbf),'eri cart 1',
     &      h_eri_1,k_eri_1)) call errquit('ma failed',911, MA_ERR)
      if (.not.ma_push_get(mt_dbl,(nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &    'eri_sp 1',
     &      h_eri_sp_1,k_eri_sp_1)) call errquit('ma failed',911,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,(nbf*nbf*nbf*nbf),
     &      'eri cart 2',
     &      h_eri_2,k_eri_2)) call errquit('ma failed',922, MA_ERR)
      if (.not.ma_push_get(mt_dbl,(nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &      'eri_sp  2',
     &      h_eri_sp_2,k_eri_sp_2)) call errquit('ma failed',922,
     &       MA_ERR)
c
c
      call rak_ovlap_test_sp(basis,nbf,nbf_sp,nshell,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_cart_s),dbl_mb(k_sph_s),
     &      int_mb(k_2bfr))
c
      call rak_2el_test_sp(basis,nbf,nbf_sp,nshell,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      int_mb(k_2bfr),
     &      dbl_mb(k_eri_sp_1),dbl_mb(k_eri_1))
c
      call rak_ovlap(basis,nbf,nshell,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_cart_s2),.false.,'cartcart')

      call rak_2el(basis,nbf,nshell,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_eri_2),
     &      .false., ' cartcart ')
      if (.not.bas_create(sp_basis,'ao sp_basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,sp_basis,'ao sp_basis'))
     &      call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
c
      write(6,*)' geom/sp_basis loaded'
c
      write(6,*)' raw sp_basis '
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0, BASIS_ERR)
      write(6,*)' first normalize'
      if (.not.int_normalize(rtdb,sp_basis)) stop ' norm error 1'
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0, BASIS_ERR)
      if (.not.bas_numbf(sp_basis,nbf_sp)) call errquit
     &      ('numbf failed',911, BASIS_ERR)
      if (.not.bas_numcont(sp_basis,nshell_sp)) call errquit
     &      ('numcont error',911, BASIS_ERR)
      write(6,*)' sp_basis b4 rak_ovlap'
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0, BASIS_ERR)
      call rak_core(sp_basis,nbf_sp,nshell_sp,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_sph_s2),.true.,'spsp')
      call rak_ovlap(sp_basis,nbf_sp,nshell_sp,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_sph_s2),.false.,'spsp')
      call rak_2el(sp_basis,nbf_sp,nshell_sp,
     &      dbl_mb(k_scr),m_scr,
     &      dbl_mb(k_buf),m_buf,
     &      dbl_mb(k_eri_sp_2),
     &      .true.,' spsp ')
c
      call print_diff_vec((nbf*nbf),
     &    dbl_mb(k_cart_s),
     &    dbl_mb(k_cart_s2),
     &    1.0d-05,' cart s ')
      call print_diff_vec((nbf_sp*nbf_sp),
     &    dbl_mb(k_sph_s),
     &    dbl_mb(k_sph_s2),
     &    1.0d-05,' spher s ')
      call daxpy((nbf*nbf),-1.0d00,
     &      dbl_mb(k_cart_s2),1,
     &      dbl_mb(k_cart_s),1)
      norm_cart = ddot((nbf*nbf),dbl_mb(k_cart_s),1,dbl_mb(k_cart_s),1)
      call daxpy((nbf*nbf),-1.0d00,
     &      dbl_mb(k_sph_s2),1,
     &      dbl_mb(k_sph_s),1)
      norm_sph = ddot((nbf*nbf),dbl_mb(k_sph_s),1,dbl_mb(k_sph_s),1)
c
      write(6,*)'1e diff norm_cart:',norm_cart
      write(6,*)'1e diff norm_sph :',norm_sph
c
      call print_diff_vec((nbf*nbf*nbf*nbf),
     &    dbl_mb(k_eri_1),
     &    dbl_mb(k_eri_2),
     &    1.0d-05,' eri cart ')
      call print_diff_vec((nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &    dbl_mb(k_eri_sp_1),
     &    dbl_mb(k_eri_sp_2),
     &    1.0d-05,' eri spherical ')
      call daxpy((nbf*nbf*nbf*nbf),-1.0d00,
     &      dbl_mb(k_eri_2),1,
     &      dbl_mb(k_eri_1),1)
      norm_cart = ddot((nbf*nbf*nbf*nbf),
     &      dbl_mb(k_eri_1),1,dbl_mb(k_eri_1),1)
      call daxpy((nbf_sp*nbf_sp*nbf_sp*nbf_sp),-1.0d00,
     &      dbl_mb(k_eri_sp_2),1,
     &      dbl_mb(k_eri_sp_1),1)
      norm_sph = ddot((nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &      dbl_mb(k_eri_sp_1),1,dbl_mb(k_eri_sp_1),1)
c
      write(6,*)'2e diff norm_cart:',norm_cart
      write(6,*)'2e diff norm_sph :',norm_sph
c
      call int_terminate()
c
      status = .true.
      status = status.and.ma_pop_stack(h_eri_sp_2)
      status = status.and.ma_pop_stack(h_eri_2)
      status = status.and.ma_pop_stack(h_eri_sp_1)
      status = status.and.ma_pop_stack(h_eri_1)
      status = status.and.ma_pop_stack(h_2bfr)
      status = status.and.ma_pop_stack(h_buf)
      status = status.and.ma_pop_stack(h_scr)
      status = status.and.ma_pop_stack(h_sph_s2)
      status = status.and.ma_pop_stack(h_cart_s2)
      status = status.and.ma_pop_stack(h_sph_s)
      status = status.and.ma_pop_stack(h_cart_s)
c
      if (.not.status) call errquit('ma pop fail',911, MA_ERR)
c
      if(.not.bas_destroy(basis)) call errquit
     &      ('basis bas_destroy failed',911, BASIS_ERR)
      if(.not.bas_destroy(sp_basis)) call errquit
     &      ('sp_basis bas_destroy failed',911, BASIS_ERR)
      if(.not.geom_destroy(geom)) call errquit
     &      ('geom_destroy failed',911, GEOM_ERR)
c
      end
*.......................................................................
      subroutine rak_2el(basis,nbf,nshell,
     &    scr,mscr,buf,mbuf,eri,print_int,msg)
      implicit none
#include "bas.fh"
#include "stdio.fh"
#include "util.fh"
      integer basis,nbf,nshell,mscr,mbuf
      double precision eri(*), scr(mscr), buf(mbuf)
      logical print_int
      character*(*) msg
c
      integer ish, jsh, ksh, lsh
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer count, indx
      logical stat_indx
      integer ii,jj,kk,ll
c
      integer i,j,k,l,isym2,isym4
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
     &               (max(isym2(i,j),isym2(k,l))-1)/2+
     &               min(isym2(i,j),isym2(k,l))
c
      write(6,*)' 2el ',msg
c
      call dfill((nbf*nbf*nbf*nbf),0.0d00,eri,1)
c
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &      stop 'cn2bfr error i'
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &        stop 'cn2bfr error j'
          do ksh = 1,ish
            if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &          stop 'cn2bfr error k'
            do lsh = 1,ksh
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi))
     &            stop 'cn2bfr error l'
              call int_2e4c
     &            (basis,ish,jsh,basis,ksh,lsh,mscr,scr,mbuf,buf)
              count = 0
              do i=ilo,ihi
                do j=jlo,jhi
                  do k=klo,khi
                    do l=llo,lhi
                      count = count + 1
                      indx = isym4(i,j,k,l)
                      stat_indx = .false.
                      if (stat_indx) then
                        write(6,*)'indx:elel:shells',msg,
     &                      indx,ish,jsh,ksh,lsh
                        write(6,*)'indx:elel:labels',msg,
     &                      indx,i,j,k,l
                      endif
                      if (print_int.and.(abs(buf(count)).gt.0.0d00))
     &                    then
                        call int_canon(i,j,k,l,ii,jj,kk,ll)
                        write(69,*)ii,jj,kk,ll,buf(count),' 2el ',msg
                      endif
                      eri(indx) = buf(count)
                    enddo
                  enddo
                enddo
              enddo
c.... end of shell loops
99999         continue
            enddo
          enddo
        enddo
      enddo

      end
*.......................................................................
      subroutine rak_2el_test_sp(basis,nbf, nbf_sp, nshell,
     &    scr,mscr,buf,mbuf,cn2bfr_sp,eri_sp,eri)
      implicit none
#include "bas.fh"
#include "stdio.fh"
#include "util.fh"
      integer nbf, nbf_sp, nshell, mscr, mbuf, basis
      integer cn2bfr_sp(2,nshell)
      double precision scr(mscr), buf(mbuf)
      double precision eri_sp(*), eri(*)
c
      integer ish, jsh, ksh, lsh
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer ilosp, jlosp, klosp, llosp
      integer ihisp, jhisp, khisp, lhisp
      integer inbf, jnbf, knbf, lnbf
      integer inbf_sp, jnbf_sp, knbf_sp, lnbf_sp
      integer itype, jtype, ktype, ltype
      integer igen, jgen, kgen, lgen
      integer iatom, jatom, katom, latom
      integer count, indx, lshtop
      logical stat_indx
*      integer junk2
      integer junk1,junk3
      double precision ttrans_w, ttrans_c
      double precision tcomp_w, tcomp_c
      double precision tadd_w, tadd_c
c
*      logical spcart_init, spcart_terminate
*      external spcart_init, spcart_terminate
c
      integer i,j,k,l,isym2,isym4
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
     &               (max(isym2(i,j),isym2(k,l))-1)/2+
     &               min(isym2(i,j),isym2(k,l))
c
*      if (.not.spcart_init(5,.true.,.false.))
*     &    stop ' spcart_init failed'
c
      call dfill((nbf_sp*nbf_sp*nbf_sp*nbf_sp),0.0d00,eri_sp,1)
      call dfill((nbf*nbf*nbf*nbf),0.0d00,eri,1)
      ttrans_w = 0.0d00
      ttrans_c = 0.0d00
      tcomp_w  = 0.0d00
      tcomp_c  = 0.0d00
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &      stop 'cn2bfr error i'
        if (.not.bas_continfo
     &      (basis,ish,itype,junk1,igen,junk3))
     &      stop 'bas_continfo error i'
        if (.not.bas_cn2ce(basis,ish,iatom))
     &      stop 'bas_cn2ce error i'
        ilosp = cn2bfr_sp(1,ish)
        ihisp = cn2bfr_sp(2,ish)
        inbf    = ihi-ilo + 1
        inbf_sp = ihisp-ilosp + 1
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &        stop 'cn2bfr error j'
          if (.not.bas_continfo
     &        (basis,jsh,jtype,junk1,jgen,junk3))
     &        stop 'bas_continfo error j'
          if (.not.bas_cn2ce(basis,jsh,jatom))
     &        stop 'bas_cn2ce error j'
          jlosp = cn2bfr_sp(1,jsh)
          jhisp = cn2bfr_sp(2,jsh)
          jnbf    = jhi-jlo + 1
          jnbf_sp = jhisp-jlosp + 1
          do ksh = 1,ish
            if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &          stop 'cn2bfr error k'
            if (.not.bas_continfo
     &          (basis,ksh,ktype,junk1,kgen,junk3))
     &          stop 'bas_continfo error k'
            if (.not.bas_cn2ce(basis,ksh,katom))
     &          stop 'bas_cn2ce error k'
            klosp = cn2bfr_sp(1,ksh)
            khisp = cn2bfr_sp(2,ksh)
            knbf    = khi-klo + 1
            knbf_sp = khisp-klosp + 1
            lshtop = ksh
            if (ksh.eq.ish) lshtop = jsh
            do lsh = 1,lshtop
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi))
     &            stop 'cn2bfr error l'
              if (.not.bas_continfo
     &            (basis,lsh,ltype,junk1,lgen,junk3))
     &            stop 'bas_continfo error l'
              if (.not.bas_cn2ce(basis,lsh,latom))
     &            stop 'bas_cn2ce error l'
              llosp = cn2bfr_sp(1,lsh)
              lhisp = cn2bfr_sp(2,lsh)
              lnbf    = lhi-llo + 1
              lnbf_sp = lhisp-llosp + 1
              tadd_c = util_cpusec()
              tadd_w = util_wallsec()
              call int_2e4c
     &            (basis,ish,jsh,basis,ksh,lsh,mscr,scr,mbuf,buf)
c
              tadd_c = util_cpusec() - tadd_c
              tadd_w = util_wallsec() - tadd_w
              if (tadd_c.gt.0.0d00) tcomp_c = tcomp_c + tadd_c
              if (tadd_w.gt.0.0d00) tcomp_w = tcomp_w + tadd_w
*rak:              write(luout,10000)ish,jsh,ksh,lsh,
*rak:     &            itype,jtype,ktype,ltype,
*rak:     &            iatom,jatom,katom,latom
              count = 0
              do i=ilo,ihi
                do j=jlo,jhi
                  do k=klo,khi
                    do l=llo,lhi
                      count = count + 1
                      indx = isym4(i,j,k,l)
                      stat_indx = indx.eq.5566
                      stat_indx = stat_indx.or.
     &                    (indx.ge.5772.and.indx.le.5784)
                      stat_indx = stat_indx.or.
     &                    (indx.ge.5887.and.indx.le.5901)
                      stat_indx = .false.
                      if (stat_indx) then
                        write(6,*)'indx:cart:shells',
     &                      indx,ish,jsh,ksh,lsh
                        write(6,*)'indx:cart:labels',
     &                      indx,i,j,k,l
                      endif
                      eri(indx) = buf(count)
                    enddo
                  enddo
                enddo
              enddo
*
              tadd_c = util_cpusec()
              tadd_w = util_wallsec()
              call spcart_bra2etran(buf,scr,
     &            jnbf,inbf,jnbf_sp,inbf_sp,
     &            jtype,itype,jgen,igen,
     &            (knbf*lnbf),.false.)
              call spcart_ket2etran(buf,scr,
     &            lnbf,knbf,lnbf_sp,knbf_sp,
     &            ltype,ktype,lgen,kgen,
     &            (inbf_sp*jnbf_sp),.false.)
c
              tadd_c = util_cpusec() - tadd_c
              tadd_w = util_wallsec() - tadd_w
              if (tadd_c.gt.0.0d00) ttrans_c = ttrans_c + tadd_c
              if (tadd_w.gt.0.0d00) ttrans_w = ttrans_w + tadd_w
*
*rak:              write(luout,10000)ish,jsh,ksh,lsh,
*rak:     &            itype,jtype,ktype,ltype,
*rak:     &            iatom,jatom,katom,latom
*rak:              count = 0
*rak:              do i=ilosp,ihisp
*rak:                do j=jlosp,jhisp
*rak:                  do k=klosp,khisp
*rak:                    do l=llosp,lhisp
*rak:                      count = count + 1
*rak:                      if (abs(buf(count)).gt.1.0d-07) then
*rak:                        write(luout,10002)i,j,k,l,buf(count),count
*rak:                      endif
*rak:                    enddo
*rak:                  enddo
*rak:                enddo
*rak:              enddo
              count = 0
              do i=ilosp,ihisp
                do j=jlosp,jhisp
                  do k=klosp,khisp
                    do l=llosp,lhisp
                      count = count + 1
                      indx = isym4(i,j,k,l)
                      stat_indx = indx.eq.5566
                      stat_indx = stat_indx.or.
     &                    (indx.ge.5772.and.indx.le.5784)
                      stat_indx = stat_indx.or.
     &                    (indx.ge.5887.and.indx.le.5901)
                      stat_indx = .false.
                      if (stat_indx) then
                        write(6,*)'indx:sp:shells',
     &                      indx,ish,jsh,ksh,lsh
                        write(6,*)'indx:sp:labels',
     &                      indx,i,j,k,l
                      endif
                      eri_sp(indx) = buf(count)
                    enddo
                  enddo
                enddo
              enddo
c.... end of shell loops
99999         continue
            enddo
          enddo
        enddo
      enddo
c
*      if (.not.spcart_terminate()) stop 'term error'
      write(luout,*)' total compute time ( cpu): ',tcomp_c
      write(luout,*)' total compute time (wall): ',tcomp_w
      write(luout,*)' total tranfrm time ( cpu): ',ttrans_c
      write(luout,*)' total tranfrm time (wall): ',ttrans_w
      if (tcomp_c.gt.1.0d-30)
     &    write(luout,'(1x,a,f10.2)')
     &    '      %    overhead ( cpu): ',
     &    (ttrans_c/tcomp_c*100.0d00)
      if (tcomp_w.gt.1.0d-30)
     &    write(luout,'(1x,a,f10.2)')
     &    '      %    overhead (wall): ',
     &    (ttrans_w/tcomp_w*100.0d00)
c
10000 format(
     &    'Shells <',i5,i5,'|',i5,i5,'>',5x,
     &    'Types  {',i5,i5,'|',i5,i5,'}',5x,
     &    'Atoms  (',i5,i5,'|',i5,i5,')')
10001 format('(',i5,i5,'|',i5,i5,') =',1pd20.10,' cart',1x,i10)
10002 format('[',i5,i5,'|',i5,i5,'] =',1pd20.10,' sphr',1x,i10)
c
      end
*.......................................................................
      subroutine rak_tospbfr(basis,nshell,nbf_chk,nbf_sp,cn2bfr)
      implicit none
#include "errquit.fh"
#include "bas.fh"
      integer nshell, nbf_chk, nbf_sp, basis
      integer cn2bfr(2,nshell)
      integer type, nprim, ngen, spsp, ish
c
      nbf_chk = 0
      nbf_sp  = 0
      do ish = 1,nshell
        if(.not.bas_continfo(basis,ish,type,nprim,ngen,spsp))
     &        call errquit(' continfo failed ',911, BASIS_ERR)
        cn2bfr(1,ish) = nbf_sp + 1
        cn2bfr(2,ish) = nbf_sp + 2*type+1
        nbf_chk = nbf_chk + (type+1)*(type+2)/2
        nbf_sp  = nbf_sp + 2*type+1
*        write(6,10000)ish,type,cn2bfr(1,ish),cn2bfr(2,ish)
      enddo
*10000 format('<ish:type><',i3,':',i3,'>  range ',i3,' to ',i3)
      end
*-------------------------------------------------------------------------------
      subroutine rak_ovlap(basis,nbf,nshell,
     &      scr,mscr,buf,mbuf,s,print_int,msg)
      implicit none
c
#include "mafdecls.fh"
#include "bas.fh"
c
      integer nbf, nshell, mscr, mbuf, basis
      double precision scr(mscr), buf(mbuf)
      double precision s(nbf,nbf)
      logical print_int
      character*(*) msg
c
      integer count
      integer ibf, ish, ilow, ihi, nbfi
      integer jbf, jsh, jlow, jhi, nbfj
      double precision value
      logical do_print
c
      call dfill((nbf*nbf),0.0d00,s,1)
      call dfill(mscr,0.0d00,scr,1)
      call dfill(mbuf,0.0d00,buf,1)
c
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'dead i'
        nbfi    = ihi    - ilow + 1
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'dead j'
          nbfj    = jhi    - jlow + 1
c
          call int_1eov(basis,ish,basis,jsh,mscr,scr,mbuf,buf)
c
          count = 0
          do ibf = ilow,ihi
            do jbf = jlow, jhi
              count = count + 1
              value = buf(count)
              s(ibf,jbf) = value
              s(jbf,ibf) = value
              do_print = print_int
              do_print = do_print .and. (ibf.ge.jbf)
              do_print = do_print .and. (abs(value).gt.0.0d00)
              if (do_print) then
                write(69,*)ibf, jbf, value, ' ovlap ',msg
              endif
            enddo
          enddo

c. close jsh/ish loops
        enddo
      enddo
      write(6,*)' generic overlap matrix (nbf=',nbf,') <',msg,'>'
      call output(s,1,nbf,1,nbf,nbf,nbf,1)
      end
*-------------------------------------------------------------------------------
      subroutine rak_core(basis,nbf,nshell,
     &      scr,mscr,buf,mbuf,s,print_int,msg)
      implicit none
c
#include "mafdecls.fh"
#include "bas.fh"
c
      integer nbf, nshell, mscr, mbuf, basis
      double precision scr(mscr), buf(mbuf)
      double precision s(nbf,nbf)
      logical print_int
      character*(*) msg
c
      integer count
      integer ibf, ish, ilow, ihi, nbfi
      integer jbf, jsh, jlow, jhi, nbfj
      double precision value
      logical do_print
c
      call dfill((nbf*nbf),0.0d00,s,1)
      call dfill(mscr,0.0d00,scr,1)
      call dfill(mbuf,0.0d00,buf,1)
c
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'dead i'
        nbfi    = ihi    - ilow + 1
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'dead j'
          nbfj    = jhi    - jlow + 1
c
          call int_1eh1(basis,ish,basis,jsh,mscr,scr,mbuf,buf)
c
          count = 0
          do ibf = ilow,ihi
            do jbf = jlow, jhi
              count = count + 1
              value = buf(count)
              s(ibf,jbf) = value
              s(jbf,ibf) = value
              do_print = print_int
              do_print = do_print .and. (ibf.ge.jbf)
              do_print = do_print .and. (abs(value).gt.0.0d00)
              if (do_print) then
                write(69,*)ibf, jbf, value, ' h1 ',msg
              endif
            enddo
          enddo

c. close jsh/ish loops
        enddo
      enddo
      write(6,*)' generic h1 matrix (nbf=',nbf,') <',msg,'>'
      call output(s,1,nbf,1,nbf,nbf,nbf,1)
      end
*-------------------------------------------------------------------------------
      subroutine rak_ovlap_test_sp(basis,nbf,nbf_sp,nshell,
     &      scr,mscr,buf,mbuf,s,s_sp,cn2bfr_sp)
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "bas.fh"
c
      integer nbf, nbf_sp, nshell, mscr, mbuf, basis
      integer cn2bfr_sp(2,nshell)
      double precision scr(mscr), buf(mbuf)
      double precision s(nbf,nbf), s_sp(nbf_sp,nbf_sp)
c
      integer ibf, ish, ilow, ihi, ilow_sp, ihi_sp, nbfi, nbfi_sp
      integer jbf, jsh, jlow, jhi, jlow_sp, jhi_sp, nbfj, nbfj_sp
      integer typei, typej, nprim, igen, jgen, spsp
      integer nint, nint_sp, count, hi_ang, st_ang
      integer ii, jj
      double precision value
*rak:      double precision pi, fact
c
*      logical spcart_init, spcart_terminate
*      external spcart_init, spcart_terminate
c
*      write(6,*)'inside rak_ovlap'
      if (.not.bas_high_angular(basis,hi_ang)) stop ' dead ang'
      st_ang = hi_ang/2
*      if (.not.spcart_init(st_ang,.true.,.false.)) stop ' dead sp'
*      if (.not.spcart_init(hi_ang,.true.,.false.)) stop ' dead sp'
c
*      write(6,*)' mscr = ',mscr
*      write(6,*)' mbuf = ',mbuf
c
      call dfill((nbf_sp*nbf_sp),0.0d00,s_sp,1)
      call dfill((nbf*nbf),0.0d00,s,1)
      call dfill(mscr,0.0d00,scr,1)
      call dfill(mbuf,0.0d00,buf,1)

      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'dead i'
        ilow_sp = cn2bfr_sp(1,ish)
        ihi_sp  = cn2bfr_sp(2,ish)
        nbfi    = ihi    - ilow + 1
        nbfi_sp = ihi_sp - ilow_sp + 1
        if(.not.bas_continfo(basis,ish,typei,nprim,igen,spsp))
     &        call errquit(' continfo failed ',911, BASIS_ERR)
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'dead j'
          jlow_sp = cn2bfr_sp(1,jsh)
          jhi_sp  = cn2bfr_sp(2,jsh)
          nbfj    = jhi    - jlow + 1
          nbfj_sp = jhi_sp - jlow_sp + 1
          if(.not.bas_continfo(basis,jsh,typej,nprim,jgen,spsp))
     &          call errquit(' continfo failed ',911, BASIS_ERR)
c
          nint    = nbfi*nbfj
          nint_sp = nbfi_sp*nbfj_sp
c
*rak:          write(6,*)' '
*rak:          write(6,*)'<ish,ilow,ihi,nbfi,typei>',ish,ilow,ihi,nbfi,typei
*rak:          write(6,*)'<jsh,jlow,jhi,nbfj,typej>',jsh,jlow,jhi,nbfj,typej
*rak:          write(6,*)' nint     = ',nint,' nint(sp) = ',nint_sp
*rak:          write(6,*)'<ish,ilowsp,ihisp,nbfisp,typei>',
*rak:     &          ish,ilow_sp,ihi_sp,nbfi_sp,typei
*rak:          write(6,*)'<jsh,jlowsp,jhisp,nbfjsp,typej>',
*rak:     &          jsh,jlow_sp,jhi_sp,nbfj_sp,typej
*rak:          write(6,*) ' ma broke 1'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma broke 1'
*rak:          write(6,*)' '
c
          call int_1eov(basis,ish,basis,jsh,mscr,scr,mbuf,buf)
c
          count = 0
          do ibf = ilow,ihi
            do jbf = jlow, jhi
              count = count + 1
              value = buf(count)
              s(ibf,jbf) = value
              s(jbf,ibf) = value
            enddo
          enddo
#define NEW_WAY
#if defined(NEW_WAY)
          call spcart_tran1e(buf,scr,
     &        nbfj,nbfi,typej,jgen,
     &        nbfj_sp,nbfi_sp,typei,igen,.false.)
#else
*rak:          write(6,*) ' ma broke 2'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma broke 2'
c.... buf is now -- buf(nbfj,nbfi)
          write(6,*)' integral buffer cart,cart '
          call output(buf,1,nbfj,1,nbfi,nbfj,nbfi,1)
          call spcart_a_s(buf,scr,nbfj,typei,.false.,.false.)
c.... scr is now -- scr(nbfj,nbfi_sp)
          write(6,*)' integral buffer  cart,sph'
          call output(scr,1,nbfj,1,nbfi_sp,nbfj,nbfi_sp,1)
*rak:          write(6,*) ' ma broke 3'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma broke 3'
          call spcart_s_a(scr,buf,nbfi_sp,typej,.false.,.false.)
c.... buf is now -- buf(nbfj_sp,nbfi_sp)
          write(6,*)' integral buffer  sph,sph'
          call output(buf,1,nbfj_sp,1,nbfi_sp,nbfj_sp,nbfi_sp,1)
*rak:          write(6,*) ' ma broke 4'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma broke 4'
#endif
          count = 0
          do ibf = ilow_sp,ihi_sp
            do jbf = jlow_sp, jhi_sp
              count = count + 1
              value = buf(count)
              s_sp(ibf,jbf) = value
              s_sp(jbf,ibf) = value
            enddo
          enddo
*rak:          write(6,*) ' ma broke 5'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma broke 5'
        enddo
      enddo
*      write(6,*)' loops done '
      write(6,*)' cartesian overlap matrix '
      call output(s,1,nbf,1,nbf,nbf,nbf,1)
      write(6,*)' nbf    ',nbf
      write(6,*)' nbf_sp ',nbf_sp
      write(6,*)' spherical overlap matrix '
      call output(s_sp,1,nbf_sp,1,nbf_sp,nbf_sp,nbf_sp,1)
      write(6,*)' nbf    ',nbf
      write(6,*)' nbf_sp ',nbf_sp
c
      count = 0 
      do ii = 1,nbf_sp
        do jj = 1,(ii-1)
          if (abs(s_sp(ii,jj)).gt.1.0d-10) then
            count = count + 1
            write(6,'(a,i3,a,i3,a,1pd30.20,i6)')
     &            'non diag element.gt.1.0d-10 s_sp(',
     &            ii,',',jj,') = ',
     &            s_sp(ii,jj),count 
          endif
        enddo
      enddo
      do ii = 1,nbf_sp
        if (abs(s_sp(ii,ii)-1.0d00).gt.1.0d-05)
     &        write(6,'(a,i3,a,i3,a,1pd30.20)')
     &        ' diagonal element.ne.1 s_sp(',ii,',',ii,') = ',
     &        s_sp(ii,ii) 
      enddo
*rak:      do ii = 2, nbf_sp
*rak:        do jj = 1,(ii-1)
*rak:          write(6,'(a,2i5,2f15.8)')
*rak:     &          ' ratios ',ii,jj,(s_sp(ii,ii)/s_sp(jj,jj)),
*rak:     &          (s_sp(jj,jj)/s_sp(ii,ii))
*rak:        enddo
*rak:      enddo
*rak:      PI=2.0d00*acos(0.0d00)
*rak:      write(6,*)pi,(pi-3.1415926535898D0)
*rak:      do ii = 1, nbf_sp
*rak:*rak:        write(6,*)' indx, val, val**2 ',ii,s_sp(ii,ii),
*rak:*rak:     &        (s_sp(ii,ii)*s_sp(ii,ii))
*rak:*rak:        write(6,*)' indx, val, 1/val ',ii,s_sp(ii,ii),
*rak:*rak:     &        (1.0d00/s_sp(ii,ii))
*rak:*rak:        write(6,*)' indx, val, sqrt ',ii,s_sp(ii,ii),
*rak:*rak:     &        (sqrt(s_sp(ii,ii)))
*rak:*rak:        write(6,*)' indx, val, /sqrt(2) ',ii,s_sp(ii,ii),
*rak:*rak:     &        (s_sp(ii,ii)/sqrt(2.0d00))
*rak:*rak:        write(6,*)' indx, val, /sqrt(3) ',ii,s_sp(ii,ii),
*rak:*rak:     &        (s_sp(ii,ii)/sqrt(3.0d00))
*rak:        fact = 1.0d00
*rak:        if (ii.ge.6.and.ii.le.10)  fact = pi*8.0d00/5.0d00
*rak:        if (ii.ge.11.and.ii.le.17) fact = pi*8.0d00/7.0d00
*rak:        if (ii.ge.18.and.ii.le.26) fact = pi*8.0d00/9.0d00
*rak:        write(6,*)' indx, val, *fact',ii,s_sp(ii,ii),
*rak:     &        (s_sp(ii,ii)*fact)
*rak:        write(6,*)' indx, val, *pi',ii,s_sp(ii,ii),
*rak:     &        (s_sp(ii,ii)*pi)
*rak:*rak:        write(6,*)' indx, val, /pi',ii,s_sp(ii,ii),
*rak:*rak:     &        (s_sp(ii,ii)/pi)
*rak:        write(6,*)' '
*rak:      enddo
*      if (.not.spcart_terminate()) stop ' dead sp term'
      end
      subroutine raktest_3ctr(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      logical int_normalize
      external int_normalize
c
c test hf3 nai type routines
      integer rtdb
      integer geom,basis, basis_id
      integer nshell, memscr, membuf
      integer h_scr, k_scr, h_buf, k_buf
      integer ish, jsh, ucont
      integer li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer nint_out
      logical status
      character*255 mo_basis, geom_name
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
     &    mo_basis = 'ao basis'
      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
     &    geom_name = 'geometry'
c
      if(.not.geom_create(geom,geom_name))call errquit
     &    ('raktest_3ctr: geom create error',911, GEOM_ERR)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('raktest_3ctr: basis create error',911, BASIS_ERR)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('raktest_3ctr: geom load ',911, RTDB_ERR)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('raktest_3ctr: basis load ',911, RTDB_ERR)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(rtdb,basis)) call errquit
     &    ('raktest_3ctr: error normalizing ',911, INT_ERR)
c
      call int_init(rtdb,1,basis)
      memscr = 100 000
      membuf = 1000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911, MA_ERR)
      if (.not.ma_push_get(mt_dbl,membuf,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911, MA_ERR)
c
      do ish = 1,nshell
        do jsh = 1,ish
          write(6,*)' ============= shells <',ish,'|',jsh,'>',
     &        '==================== start =========='
          write(6,*)' '
          
          ucont = (sf_ibs_cn2ucn(ish,basis_id))
          Li      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          i_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          i_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          i_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          i_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          i_cent  = (sf_ibs_cn2ce(ish,basis_id))
          i_geom  = ibs_geom(basis_id)
c
          ucont = (sf_ibs_cn2ucn(jsh,basis_id))
          Lj      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          j_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          j_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          j_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          j_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          j_cent  = (sf_ibs_cn2ce(jsh,basis_id))
          j_geom  = ibs_geom(basis_id)
          
          call hf1tmp(
     &          coords(1,i_cent,i_geom),
     &          dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(i_icfp,basis_id)), i_prim, i_gen, Li,
     &          coords(1,j_cent,j_geom),
     &          dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(j_icfp,basis_id)), j_prim, j_gen, Lj,
     &          coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
     &          dbl_mb(k_scr),dbl_mb(k_scr),dbl_mb(k_buf),membuf,
     &          .false., .false., .true., .false., .false.,
     &          dbl_mb(k_scr), memscr)
          write(6,*)' i = c center '
          call hf3pot(
     &          coords(1,i_cent,i_geom),
     &          dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(i_icfp,basis_id)), i_prim, i_gen, Li,
     &          coords(1,j_cent,j_geom),
     &          dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(j_icfp,basis_id)), j_prim, j_gen, Lj,
     &          coords(1,i_cent,i_geom),0.0d00, 1.0d00, 1, 1, 0,
     &          dbl_mb(k_buf), membuf, nint_out, .false.,
     &          dbl_mb(k_scr), memscr)
          write(6,*)' j = c center '
          call hf3pot(
     &          coords(1,i_cent,i_geom),
     &          dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(i_icfp,basis_id)), i_prim, i_gen, Li,
     &          coords(1,j_cent,j_geom),
     &          dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(j_icfp,basis_id)), j_prim, j_gen, Lj,
     &          coords(1,j_cent,j_geom),0.0d00, 1.0d00, 1, 1, 0,
     &          dbl_mb(k_buf), membuf, nint_out, .false.,
     &          dbl_mb(k_scr), memscr)
          write(6,*)' i = c center swap'
          call hf3pot(
     &          coords(1,i_cent,i_geom),0.0d00, 1.0d00, 1, 1, 0,
     &          coords(1,j_cent,j_geom),
     &          dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(j_icfp,basis_id)), j_prim, j_gen, Lj,
     &          coords(1,i_cent,i_geom),
     &          dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(i_icfp,basis_id)), i_prim, i_gen, Li,
     &          dbl_mb(k_buf), membuf, nint_out, .false.,
     &          dbl_mb(k_scr), memscr)
          write(6,*)' j = c center swap'
          call hf3pot(
     &          coords(1,i_cent,i_geom),
     &          dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(i_icfp,basis_id)), i_prim, i_gen, Li,
     &          coords(1,j_cent,j_geom),0.0d00, 1.0d00, 1, 1, 0,
     &          coords(1,j_cent,j_geom),
     &          dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &          dbl_mb(mb_exndcf(j_icfp,basis_id)), j_prim, j_gen, Lj,
     &          dbl_mb(k_buf), membuf, nint_out, .false.,
     &          dbl_mb(k_scr), memscr)
          write(6,*)' ============= shells <',ish,'|',jsh,'>',
     &          '====================  end  =========='
          write(6,*)' '
        enddo
      enddo
c      
      call int_terminate()
      status = ma_pop_stack(h_buf)
      status = status.and.ma_pop_stack(h_scr)
      if (.not.status) call errquit('pop failed',911, MA_ERR)
      status = bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      if (.not.status) call errquit('b/g destroy failed',911, MA_ERR)
      end
      subroutine raktest_stpr(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
      integer stpr_walk
      external stpr_walk
c
      integer geom
      integer nat
      integer k_grad, h_grad
      logical flag
      double precision energy
c
      if (ga_nodeid().eq.0)  then
        if (.not. geom_create(geom, 'geometry'))
     &         call errquit('raktest_stpr: geom_create?', 911, GEOM_ERR)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &         call errquit('raktest_stpr: no geometry ', 911, RTDB_ERR)
c get number of atoms = nat
        if (.not. geom_ncent(geom,nat))
     &      call errquit('raktest_stpr: geom_ncent?',911, GEOM_ERR)
        if (.not. geom_destroy(geom))
     &      call errquit('raktest_stpr: geom_destroy?',911, GEOM_ERR)
        if (.not.
     &         MA_Push_Get(MT_DBL,(3*nat),'stpr fake gradient',
     &         h_grad,k_grad))
     &         call errquit
     &         ('raktest_stpr: allocation for gradient failed?',911,
     &       MA_ERR)
        dbl_mb((k_grad+   0)) =  -0.008697D00
        dbl_mb((k_grad+   1)) =  -0.004076D00
        dbl_mb((k_grad+   2)) =   0.004591D00
        dbl_mb((k_grad+   3)) =   0.000926D00
        dbl_mb((k_grad+   4)) =   0.000922D00
        dbl_mb((k_grad+   5)) =  -0.001206D00
        dbl_mb((k_grad+   6)) =   0.000639D00
        dbl_mb((k_grad+   7)) =   0.000263D00
        dbl_mb((k_grad+   8)) =  -0.003188D00
        dbl_mb((k_grad+   9)) =   0.007578D00
        dbl_mb((k_grad+  10)) =  -0.002686D00
        dbl_mb((k_grad+  11)) =   0.005190D00
        dbl_mb((k_grad+  12)) =  -0.000872D00
        dbl_mb((k_grad+  13)) =   0.000692D00
        dbl_mb((k_grad+  14)) =  -0.001103D00
        dbl_mb((k_grad+  15)) =  -0.000782D00
        dbl_mb((k_grad+  16)) =   0.000717D00
        dbl_mb((k_grad+  17)) =  -0.001393D00
        dbl_mb((k_grad+  18)) =  -0.008078D00
        dbl_mb((k_grad+  19)) =   0.006782D00
        dbl_mb((k_grad+  20)) =   0.000074D00
        dbl_mb((k_grad+  21)) =   0.000549D00
        dbl_mb((k_grad+  22)) =  -0.000865D00
        dbl_mb((k_grad+  23)) =  -0.002881D00
        dbl_mb((k_grad+  24)) =   0.006377D00
        dbl_mb((k_grad+  25)) =   0.001920D00
        dbl_mb((k_grad+  26)) =   0.000909D00
        dbl_mb((k_grad+  27)) =  -0.000695D00
        dbl_mb((k_grad+  28)) =  -0.000362D00
        dbl_mb((k_grad+  29)) =  -0.000411D00
        dbl_mb((k_grad+  30)) =  -0.001585D00
        dbl_mb((k_grad+  31)) =  -0.000732D00
        dbl_mb((k_grad+  32)) =  -0.000523D00
        dbl_mb((k_grad+  33)) =  -0.001062D00
        dbl_mb((k_grad+  34)) =  -0.000824D00
        dbl_mb((k_grad+  35)) =  -0.000303D00
        dbl_mb((k_grad+  36)) =  -0.022866D00
        dbl_mb((k_grad+  37)) =   0.011419D00
        dbl_mb((k_grad+  38)) =  -0.010230D00
        dbl_mb((k_grad+  39)) =   0.013343D00
        dbl_mb((k_grad+  40)) =  -0.017703D00
        dbl_mb((k_grad+  41)) =  -0.018748D00
        dbl_mb((k_grad+  42)) =   0.002823D00
        dbl_mb((k_grad+  43)) =   0.000704D00
        dbl_mb((k_grad+  44)) =   0.000638D00
        dbl_mb((k_grad+  45)) =   0.001783D00
        dbl_mb((k_grad+  46)) =   0.000141D00
        dbl_mb((k_grad+  47)) =   0.000499D00
        dbl_mb((k_grad+  48)) =   0.003863D00
        dbl_mb((k_grad+  49)) =   0.007955D00
        dbl_mb((k_grad+  50)) =  -0.003668D00
        dbl_mb((k_grad+  51)) =  -0.000253D00
        dbl_mb((k_grad+  52)) =  -0.000335D00
        dbl_mb((k_grad+  53)) =   0.003427D00
        dbl_mb((k_grad+  54)) =   0.003510D00
        dbl_mb((k_grad+  55)) =  -0.004396D00
        dbl_mb((k_grad+  56)) =  -0.000280D00
        dbl_mb((k_grad+  57)) =  -0.001050D00
        dbl_mb((k_grad+  58)) =   0.000598D00
        dbl_mb((k_grad+  59)) =   0.000407D00
        dbl_mb((k_grad+  60)) =  -0.001190D00
        dbl_mb((k_grad+  61)) =   0.000896D00
        dbl_mb((k_grad+  62)) =   0.000453D00
        dbl_mb((k_grad+  63)) =  -0.000536D00
        dbl_mb((k_grad+  64)) =   0.000599D00
        dbl_mb((k_grad+  65)) =   0.000487D00
        dbl_mb((k_grad+  66)) =   0.008229D00
        dbl_mb((k_grad+  67)) =   0.006460D00
        dbl_mb((k_grad+  68)) =  -0.000069D00
        dbl_mb((k_grad+  69)) =  -0.000490D00
        dbl_mb((k_grad+  70)) =  -0.000613D00
        dbl_mb((k_grad+  71)) =  -0.003029D00
        dbl_mb((k_grad+  72)) =  -0.006348D00
        dbl_mb((k_grad+  73)) =   0.002016D00
        dbl_mb((k_grad+  74)) =   0.000844D00
        dbl_mb((k_grad+  75)) =   0.001549D00
        dbl_mb((k_grad+  76)) =  -0.000898D00
        dbl_mb((k_grad+  77)) =  -0.000570D00
        dbl_mb((k_grad+  78)) =   0.000697D00
        dbl_mb((k_grad+  79)) =  -0.000361D00
        dbl_mb((k_grad+  80)) =  -0.000393D00
        dbl_mb((k_grad+  81)) =   0.001078D00
        dbl_mb((k_grad+  82)) =  -0.000863D00
        dbl_mb((k_grad+  83)) =  -0.000271D00
        dbl_mb((k_grad+  84)) =  -0.007499D00
        dbl_mb((k_grad+  85)) =  -0.002478D00
        dbl_mb((k_grad+  86)) =   0.005769D00
        dbl_mb((k_grad+  87)) =   0.000808D00
        dbl_mb((k_grad+  88)) =   0.000701D00
        dbl_mb((k_grad+  89)) =  -0.001448D00
        dbl_mb((k_grad+  90)) =   0.000827D00
        dbl_mb((k_grad+  91)) =   0.000750D00
        dbl_mb((k_grad+  92)) =  -0.001107D00
        dbl_mb((k_grad+  93)) =   0.007929D00
        dbl_mb((k_grad+  94)) =  -0.004455D00
        dbl_mb((k_grad+  95)) =   0.003384D00
        dbl_mb((k_grad+  96)) =  -0.000671D00
        dbl_mb((k_grad+  97)) =   0.000441D00
        dbl_mb((k_grad+  98)) =  -0.003097D00
        dbl_mb((k_grad+  99)) =  -0.000863D00
        dbl_mb((k_grad+ 100)) =   0.000987D00
        dbl_mb((k_grad+ 101)) =  -0.001353D00
        dbl_mb((k_grad+ 102)) =   0.007929D00
        dbl_mb((k_grad+ 103)) =   0.004455D00
        dbl_mb((k_grad+ 104)) =  -0.003384D00
        dbl_mb((k_grad+ 105)) =  -0.000671D00
        dbl_mb((k_grad+ 106)) =  -0.000441D00
        dbl_mb((k_grad+ 107)) =   0.003097D00
        dbl_mb((k_grad+ 108)) =  -0.000863D00
        dbl_mb((k_grad+ 109)) =  -0.000987D00
        dbl_mb((k_grad+ 110)) =   0.001353D00
        dbl_mb((k_grad+ 111)) =  -0.007499D00
        dbl_mb((k_grad+ 112)) =   0.002477D00
        dbl_mb((k_grad+ 113)) =  -0.005768D00
        dbl_mb((k_grad+ 114)) =   0.000808D00
        dbl_mb((k_grad+ 115)) =  -0.000701D00
        dbl_mb((k_grad+ 116)) =   0.001448D00
        dbl_mb((k_grad+ 117)) =   0.000827D00
        dbl_mb((k_grad+ 118)) =  -0.000751D00
        dbl_mb((k_grad+ 119)) =   0.001107D00
        dbl_mb((k_grad+ 120)) =   0.008230D00
        dbl_mb((k_grad+ 121)) =  -0.006460D00
        dbl_mb((k_grad+ 122)) =   0.000068D00
        dbl_mb((k_grad+ 123)) =  -0.000491D00
        dbl_mb((k_grad+ 124)) =   0.000613D00
        dbl_mb((k_grad+ 125)) =   0.003029D00
        dbl_mb((k_grad+ 126)) =  -0.006349D00
        dbl_mb((k_grad+ 127)) =  -0.002015D00
        dbl_mb((k_grad+ 128)) =  -0.000844D00
        dbl_mb((k_grad+ 129)) =   0.001078D00
        dbl_mb((k_grad+ 130)) =   0.000863D00
        dbl_mb((k_grad+ 131)) =   0.000270D00
        dbl_mb((k_grad+ 132)) =   0.001549D00
        dbl_mb((k_grad+ 133)) =   0.000898D00
        dbl_mb((k_grad+ 134)) =   0.000570D00
        dbl_mb((k_grad+ 135)) =   0.000697D00
        dbl_mb((k_grad+ 136)) =   0.000361D00
        dbl_mb((k_grad+ 137)) =   0.000393D00
        dbl_mb((k_grad+ 138)) =   0.003863D00
        dbl_mb((k_grad+ 139)) =  -0.007955D00
        dbl_mb((k_grad+ 140)) =   0.003667D00
        dbl_mb((k_grad+ 141)) =  -0.000253D00
        dbl_mb((k_grad+ 142)) =   0.000335D00
        dbl_mb((k_grad+ 143)) =  -0.003427D00
        dbl_mb((k_grad+ 144)) =   0.003510D00
        dbl_mb((k_grad+ 145)) =   0.004397D00
        dbl_mb((k_grad+ 146)) =   0.000279D00
        dbl_mb((k_grad+ 147)) =  -0.001190D00
        dbl_mb((k_grad+ 148)) =  -0.000897D00
        dbl_mb((k_grad+ 149)) =  -0.000453D00
        dbl_mb((k_grad+ 150)) =  -0.000536D00
        dbl_mb((k_grad+ 151)) =  -0.000599D00
        dbl_mb((k_grad+ 152)) =  -0.000487D00
        dbl_mb((k_grad+ 153)) =  -0.001050D00
        dbl_mb((k_grad+ 154)) =  -0.000598D00
        dbl_mb((k_grad+ 155)) =  -0.000407D00
        dbl_mb((k_grad+ 156)) =   0.013343D00
        dbl_mb((k_grad+ 157)) =   0.017702D00
        dbl_mb((k_grad+ 158)) =   0.018749D00
        dbl_mb((k_grad+ 159)) =   0.002823D00
        dbl_mb((k_grad+ 160)) =  -0.000703D00
        dbl_mb((k_grad+ 161)) =  -0.000639D00
        dbl_mb((k_grad+ 162)) =   0.001783D00
        dbl_mb((k_grad+ 163)) =  -0.000141D00
        dbl_mb((k_grad+ 164)) =  -0.000499D00
        dbl_mb((k_grad+ 165)) =  -0.022865D00
        dbl_mb((k_grad+ 166)) =  -0.011418D00
        dbl_mb((k_grad+ 167)) =   0.010231D00
        dbl_mb((k_grad+ 168)) =   0.000385D00
        dbl_mb((k_grad+ 169)) =   0.023722D00
        dbl_mb((k_grad+ 170)) =   0.016418D00
        dbl_mb((k_grad+ 171)) =   0.022263D00
        dbl_mb((k_grad+ 172)) =  -0.013938D00
        dbl_mb((k_grad+ 173)) =  -0.000893D00
        dbl_mb((k_grad+ 174)) =   0.000883D00
        dbl_mb((k_grad+ 175)) =  -0.000287D00
        dbl_mb((k_grad+ 176)) =  -0.002376D00
        dbl_mb((k_grad+ 177)) =  -0.022844D00
        dbl_mb((k_grad+ 178)) =  -0.010601D00
        dbl_mb((k_grad+ 179)) =  -0.000400D00
        dbl_mb((k_grad+ 180)) =  -0.000446D00
        dbl_mb((k_grad+ 181)) =  -0.001743D00
        dbl_mb((k_grad+ 182)) =  -0.001937D00
        dbl_mb((k_grad+ 183)) =   0.020689D00
        dbl_mb((k_grad+ 184)) =  -0.010232D00
        dbl_mb((k_grad+ 185)) =   0.017869D00
        dbl_mb((k_grad+ 186)) =  -0.022845D00
        dbl_mb((k_grad+ 187)) =   0.010601D00
        dbl_mb((k_grad+ 188)) =   0.000400D00
        dbl_mb((k_grad+ 189)) =  -0.000446D00
        dbl_mb((k_grad+ 190)) =   0.001743D00
        dbl_mb((k_grad+ 191)) =   0.001937D00
        dbl_mb((k_grad+ 192)) =   0.022264D00
        dbl_mb((k_grad+ 193)) =   0.013938D00
        dbl_mb((k_grad+ 194)) =   0.000892D00
        dbl_mb((k_grad+ 195)) =   0.000883D00
        dbl_mb((k_grad+ 196)) =   0.000287D00
        dbl_mb((k_grad+ 197)) =   0.002376D00
        dbl_mb((k_grad+ 198)) =   0.020690D00
        dbl_mb((k_grad+ 199)) =   0.010232D00
        dbl_mb((k_grad+ 200)) =  -0.017869D00
        dbl_mb((k_grad+ 201)) =   0.002430D00
        dbl_mb((k_grad+ 202)) =  -0.022402D00
        dbl_mb((k_grad+ 203)) =   0.003835D00
        dbl_mb((k_grad+ 204)) =  -0.001212D00
        dbl_mb((k_grad+ 205)) =  -0.000669D00
        dbl_mb((k_grad+ 206)) =   0.002178D00
        dbl_mb((k_grad+ 207)) =  -0.002379D00
        dbl_mb((k_grad+ 208)) =   0.025844D00
        dbl_mb((k_grad+ 209)) =  -0.005825D00
        dbl_mb((k_grad+ 210)) =   0.001510D00
        dbl_mb((k_grad+ 211)) =   0.000491D00
        dbl_mb((k_grad+ 212)) =  -0.001989D00
        dbl_mb((k_grad+ 213)) =   0.000386D00
        dbl_mb((k_grad+ 214)) =  -0.023721D00
        dbl_mb((k_grad+ 215)) =  -0.016418D00
        dbl_mb((k_grad+ 216)) =   0.022354D00
        dbl_mb((k_grad+ 217)) =  -0.009889D00
        dbl_mb((k_grad+ 218)) =  -0.001034D00
        dbl_mb((k_grad+ 219)) =   0.000445D00
        dbl_mb((k_grad+ 220)) =  -0.001737D00
        dbl_mb((k_grad+ 221)) =  -0.001918D00
        dbl_mb((k_grad+ 222)) =  -0.023384D00
        dbl_mb((k_grad+ 223)) =  -0.006755D00
        dbl_mb((k_grad+ 224)) =   0.020970D00
        dbl_mb((k_grad+ 225)) =  -0.000504D00
        dbl_mb((k_grad+ 226)) =   0.001171D00
        dbl_mb((k_grad+ 227)) =   0.000301D00
        dbl_mb((k_grad+ 228)) =  -0.019164D00
        dbl_mb((k_grad+ 229)) =   0.008952D00
        dbl_mb((k_grad+ 230)) =  -0.018541D00
        dbl_mb((k_grad+ 231)) =  -0.023385D00
        dbl_mb((k_grad+ 232)) =   0.006755D00
        dbl_mb((k_grad+ 233)) =  -0.020970D00
        dbl_mb((k_grad+ 234)) =  -0.000503D00
        dbl_mb((k_grad+ 235)) =  -0.001171D00
        dbl_mb((k_grad+ 236)) =  -0.000301D00
        dbl_mb((k_grad+ 237)) =   0.022354D00
        dbl_mb((k_grad+ 238)) =   0.009888D00
        dbl_mb((k_grad+ 239)) =   0.001034D00
        dbl_mb((k_grad+ 240)) =   0.000445D00
        dbl_mb((k_grad+ 241)) =   0.001737D00
        dbl_mb((k_grad+ 242)) =   0.001918D00
        dbl_mb((k_grad+ 243)) =  -0.019164D00
        dbl_mb((k_grad+ 244)) =  -0.008951D00
        dbl_mb((k_grad+ 245)) =   0.018541D00
        dbl_mb((k_grad+ 246)) =  -0.002378D00
        dbl_mb((k_grad+ 247)) =  -0.025843D00
        dbl_mb((k_grad+ 248)) =   0.005825D00
        dbl_mb((k_grad+ 249)) =   0.001510D00
        dbl_mb((k_grad+ 250)) =  -0.000491D00
        dbl_mb((k_grad+ 251)) =   0.001990D00
        dbl_mb((k_grad+ 252)) =   0.002430D00
        dbl_mb((k_grad+ 253)) =   0.022402D00
        dbl_mb((k_grad+ 254)) =  -0.003835D00
        dbl_mb((k_grad+ 255)) =  -0.001213D00
        dbl_mb((k_grad+ 256)) =   0.000669D00
        dbl_mb((k_grad+ 257)) =  -0.002178D00
        dbl_mb((k_grad+ 258)) =  -0.008078D00
        dbl_mb((k_grad+ 259)) =  -0.006781D00
        dbl_mb((k_grad+ 260)) =  -0.000074D00
        dbl_mb((k_grad+ 261)) =   0.000549D00
        dbl_mb((k_grad+ 262)) =   0.000865D00
        dbl_mb((k_grad+ 263)) =   0.002881D00
        dbl_mb((k_grad+ 264)) =   0.006378D00
        dbl_mb((k_grad+ 265)) =  -0.001920D00
        dbl_mb((k_grad+ 266)) =  -0.000909D00
        dbl_mb((k_grad+ 267)) =  -0.000695D00
        dbl_mb((k_grad+ 268)) =   0.000362D00
        dbl_mb((k_grad+ 269)) =   0.000411D00
        dbl_mb((k_grad+ 270)) =  -0.001585D00
        dbl_mb((k_grad+ 271)) =   0.000732D00
        dbl_mb((k_grad+ 272)) =   0.000524D00
        dbl_mb((k_grad+ 273)) =  -0.001062D00
        dbl_mb((k_grad+ 274)) =   0.000824D00
        dbl_mb((k_grad+ 275)) =   0.000303D00
        dbl_mb((k_grad+ 276)) =   0.007577D00
        dbl_mb((k_grad+ 277)) =   0.002686D00
        dbl_mb((k_grad+ 278)) =  -0.005190D00
        dbl_mb((k_grad+ 279)) =  -0.000782D00
        dbl_mb((k_grad+ 280)) =  -0.000716D00
        dbl_mb((k_grad+ 281)) =   0.001394D00
        dbl_mb((k_grad+ 282)) =  -0.000872D00
        dbl_mb((k_grad+ 283)) =  -0.000692D00
        dbl_mb((k_grad+ 284)) =   0.001103D00
        dbl_mb((k_grad+ 285)) =  -0.008696D00
        dbl_mb((k_grad+ 286)) =   0.004077D00
        dbl_mb((k_grad+ 287)) =  -0.004592D00
        dbl_mb((k_grad+ 288)) =   0.000639D00
        dbl_mb((k_grad+ 289)) =  -0.000263D00
        dbl_mb((k_grad+ 290)) =   0.003187D00
        dbl_mb((k_grad+ 291)) =   0.000926D00
        dbl_mb((k_grad+ 292)) =  -0.000923D00
        dbl_mb((k_grad+ 293)) =   0.001207D00
        energy = -317.5656589D00
        flag  = .true.
c put scf:converged = logical true
        if (.not. rtdb_put(rtdb,'scf:converged', MT_LOG, 1, flag))
     &      call errquit
     &      ('raktest_stpr: failed to read converged in rtdb', 911,
     &       RTDB_ERR)
c put scf:energy   =  real value
        if (.not. rtdb_put(rtdb,'scf:energy', MT_DBL, 1, energy))
     &      call errquit
     &      ('raktest_stpr: failed to read energy in rtdb', 911,
     &       RTDB_ERR)
c put scf:gradient = 3*nat reals
        if (.not. rtdb_put(rtdb, 'scf:gradient', MT_DBL, 
     &       (3*nat),dbl_mb(k_grad)))
     &      call errquit
     &        ('raktest_stpr: reading gradients failed',911, RTDB_ERR)
c free memory
        if (.not. ma_pop_stack(h_grad))
     &      call errquit('raktest_stpr: pop failed',911, MA_ERR)
      endif
      call ga_sync()
      call ga_sync()
      call ga_sync()
      if (stpr_walk(rtdb).eq.1) then
        write(6,*)' walk converged'
      else
        write(6,*)' walk NOT converged'
      endif
      end
      subroutine raktest_initd(rtdb)
c raktest = 4
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer geom
      integer mx1e, mxg, mxs1, mxs2
c
      logical status
c
      integer nbas, bases(6), mynbas
c
      if (.not.bas_rtdb_in(rtdb))
     &    call errquit('raktest4: error loading known basis sets',911,
     &       BASIS_ERR)
c
      write(6,*)' number of basis sets in rtdb ',nbasis_rtdb
c
      do 00100 nbas = 1,nbasis_rtdb
        write(6,*)' basis ',nbas,' is ',bs_names_rtdb(nbas)
00100 continue
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911, RTDB_ERR)
c
      mynbas = 0
      do 00200 nbas = 1,nbasis_rtdb
        if (bs_names_rtdb(nbas).ne.'ecp basis') then
          mynbas = mynbas + 1
          if(.not.bas_create(bases(mynbas),bs_names_rtdb(mynbas)))
     &        call errquit('raktest4: bas_create choked',911, BASIS_ERR)
          if(.not.
     &        bas_rtdb_load
     &        (rtdb,geom,bases(mynbas),bs_names_rtdb(mynbas)))
     &        call errquit('raktest4: bas_rtdb_load choked',911,
     &       RTDB_ERR)
          status = bas_print(bases(mynbas))
          status = gbs_map_print(bases(mynbas))
        endif
00200 continue
      call intd_init(rtdb,mynbas,bases)
c
      call int_mem(mx1e, mxg, mxs1, mxs2)
      write(6,*)' one electron buffer size        :',mx1e
      write(6,*)' two electron buffer size        :',mxg
      write(6,*)' one electron scratch buffer size:',mxs1
      write(6,*)' two electron scratch buffer size:',mxs2
c
      call int_mem_print()
c
      do nbas = 1,mynbas
        if (.not.bas_destroy(bases(nbas))) call errquit
     &        ('raktest_initd: bas_destroy failed',911, BASIS_ERR)
      enddo
      if (.not.geom_destroy(geom)) call errquit
     &      ('raktest_initd: _destroy failed',911, GEOM_ERR)
c
      call intd_terminate()
c
      end
      subroutine raktest_init(rtdb)
c raktest = 5
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer geom
      integer mx1e, mxg, mxs1, mxs2
c
      logical status
c
      integer nbas, bases(6)
c
      if (.not.bas_rtdb_in(rtdb))
     &    call errquit('raktest4: error loading known basis sets',911,
     &       BASIS_ERR)
c
      write(6,*)' number of basis sets in rtdb ',nbasis_rtdb
c
      do 00100 nbas = 1,nbasis_rtdb
        write(6,*)' basis ',nbas,' is ',bs_names_rtdb(nbas)
00100 continue
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911, RTDB_ERR)
c
      do 00200 nbas = 1,nbasis_rtdb
        if(.not.bas_create(bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_create choked',911, BASIS_ERR)
        if(.not.
     &      bas_rtdb_load
     &      (rtdb,geom,bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_rtdb_load choked',911, RTDB_ERR)
        status = bas_print(bases(nbas))
        status = gbs_map_print(bases(nbas))
00200 continue
      call int_init(rtdb,nbasis_rtdb,bases)
c
      call int_mem(mx1e, mxg, mxs1, mxs2)
      write(6,*)' one electron buffer size        :',mx1e
      write(6,*)' two electron buffer size        :',mxg
      write(6,*)' one electron scratch buffer size:',mxs1
      write(6,*)' two electron scratch buffer size:',mxs2
c
      call int_mem_print()
c
      do nbas = 1,nbasis_rtdb
        if (.not.bas_destroy(bases(nbas))) call errquit
     &        ('raktest_init: bas_destroy failed',911, BASIS_ERR)
      enddo
      if (.not.geom_destroy(geom)) call errquit
     &      ('raktest_init: _destroy failed',911, GEOM_ERR)
c
      call int_terminate()
c
      end

      Subroutine hf1mkrtmp(Axyz,Bxyz,Cxyz,zan,ncenters,
     &                  alpha,Pxyz,RS,PC,ff,R,
     &                  R0,R0C,IJK,NPP,Lp,Lp3,CENTER)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical CENTER

      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)

c--> Cartesian Coordinates

      Dimension Axyz(3),Bxyz(3)

c--> Nuclear Cartesian Coordinates & Charges

      Dimension Cxyz(3,ncenters),zan(ncenters)

c--> Exponents

      Dimension alpha(2,NPP)

c--> Auxiliary Function Integrals & Index

      Dimension R0(NPP,Lp3),R0C(ncenters,NPP,Lp3),IJK(0:Lp,0:Lp,0:Lp)

c--> Scratch Space

      Dimension Pxyz(3,NPP),RS(NPP),PC(NPP,3),ff(2,NPP),R(NPP,0:Lp,Lp3)

c
c Define the auxiliary function integrals necessary to compute the nuclear
c attraction integrals (NAIs). These integrals are scaled by an appropriate
c factor, RS, defined as
c
c         / a + b \ 1/2
c   RS = | ------- |
c         \  PI/4 /
c
c The scale factor for the Hermite expansion coefficients is assumed to be
c
c         /   PI  \ 3/2     /   a b   __2 \
c   ES = | ------- |    EXP| - -----  AB   |
c         \ a + b /         \  a + b      /
c
c Therefore,
c
c            2 PI        /   a b   __2 \
c   ES RS = -------  EXP| - -----  AB   |
c            a + b       \  a + b      /
c
c******************************************************************************

      do 100 mp = 1,NPP
        do 100 j = 1,Lp3
          R0(mp,j) = 0.D0
  100 continue

      do 110 mp = 1,NPP

c Define the center "P".

       a = alpha(1,mp)
       b = alpha(2,mp)

       f1 = a/(a+b)
       f2 = b/(a+b)

       Pxyz(1,mp) = f1*Axyz(1) + f2*Bxyz(1)
       Pxyz(2,mp) = f1*Axyz(2) + f2*Bxyz(2)
       Pxyz(3,mp) = f1*Axyz(3) + f2*Bxyz(3)

c Define the scaling factor.

       RS(mp) = sqrt((a+b)*PI4)

  110 continue

c Sum over all centers.

      do 150 ic = 1,ncenters

c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.

       do 120 m = 1,NPP

        alpha_t = alpha(1,m) + alpha(2,m)

        ff(1,m) = RS(m)
        ff(2,m) = -2.D0*alpha_t

        PCx = Pxyz(1,m) - Cxyz(1,ic)
        PCy = Pxyz(2,m) - Cxyz(2,ic)
        PCz = Pxyz(3,m) - Cxyz(3,ic)

        R(m,0,1) = alpha_t*(PCx**2 + PCy**2 + PCz**2)

        PC(m,1) = PCx
        PC(m,2) = PCy
        PC(m,3) = PCz

  120  continue

c Evaluate the incomplete gamma function.

       call igamma(R,NPP,Lp)

c Define the initial auxiliary functions (i.e., R000j, j=1,Lr).

       do 135 j = 0,Lp
        do 130 m = 1,NPP
         R(m,j,1) = ff(1,m)*R(m,j,1)
         ff(1,m) = ff(1,m)*ff(2,m)
  130   continue
  135  continue

c Recursively build the remaining auxiliary functions (i.e., RIJKj, j=0).

       call hfmkr(R,IJK,PC,NPP,Lp,Lp3)

c Transfer to R0 array.

       if( CENTER )then
        do 141 n = 1,Lp3
         do 140 m = 1,NPP
          R0C(ic,m,n) = -zan(ic)*R(m,0,n)
          R0(m,n) = R0(m,n) + R0C(ic,m,n)
  140    continue
  141   continue
       else
        do 146 n = 1,Lp3
         do 145 m = 1,NPP
c*          R0(m,n) = R0(m,n) - zan(ic)*R(m,0,n)
          R0C(ic,m,n) = R(m,0,n)
  145    continue
  146   continue
       end if

  150 continue

      end
      Subroutine hf1tmp(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,zan,ncenters,
     &               bO2I,bKEI,bNAI,Nint,O2I,KEI,NAI,canAB,
     &               DryRun,W0,maxW0)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
#include "errquit.fh"

      Logical O2I,KEI,NAI,canAB

      Logical GenCon,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates & Charges

      Dimension Cxyz(3,ncenters),zan(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      Dimension bO2I(Nint),bKEI(Nint),bNAI(Nint)

c--> Scratch Space.

      Dimension W0(maxW0)
c
c Compute the overlap, kinetic energy, and nuclear attraction integrals for 
c two shells of contracted Gaussians functions. This driver is NOT capable of 
c evaluating integral derivatives.
c
c******************************************************************************
#if defined(INTDEBUG)
      write(6,*)' inside hf1 '
      write(6,*)' npa,nca,la = ',npa,nca,la
      write(6,*)' npb,ncb,lb = ',npb,ncb,lb
      write(6,*)' ncenters   = ',ncenters
      write(6,*)' NINT       = ',nint
      write(6,*)' maxW0      = ',maxw0
      write(6,*)' <canAB:DryRun>-<',canab,':',dryrun,'>'
      write(6,*)' <o2i:kei:nai>-<',o2i,':',kei,':',nai,'>'
      write(6,'(a,3(2x,1pd20.10))')' Axyz =',Axyz
      write(6,'(a,3(2x,1pd20.10))')' Bxyz =',Bxyz
      write(6,'(a,100(3(2x,1pd20.10/)))')' Cxyz =',Cxyz
      do iiii = 1,npa
        write(6,'(a,i3,a,2(2x,1pd20.10))')
     &         'Aprims:Acoeffs:(',iiii,') =',Aprims(iiii),
     &         Acoefs(iiii,1)
      enddo
      do iiii = 1,npb
        write(6,'(a,i3,a,2(2x,1pd20.10))')
     &         'Bprims:Bcoeffs:(',iiii,') =',Bprims(iiii),
     &         Bcoefs(iiii,1)
      enddo
#endif
*      if (.not.dryrun) then
*        call hf_print('hf1: a shell',axyz,aprims,acoefs,npa,nca,la)
*        call hf_print('hf1: b shell',bxyz,bprims,bcoefs,npb,ncb,lb)
*      endif
      MXD = 0
      if (KEI) call errquit('hf1tmp: only for pot ints',911, INT_ERR)
      if (O2I) call errquit('hf1tmp: only for pot ints',911, INT_ERR)

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB

      GenCon = NCP.ne.1

      if( GenCon )then
       write(*,*) 'HF1: Not yet ready for general contraction.'
       stop
      end if

c To determine all the Hermite expansion coefficients required to evaluate
c the kinetic energy integrals, increment the angular momenta by one.

      if( KEI )then
       Li = 1
      else
       Li = 0
      end if

c Define the angular momentum of the overlap distribution.

      Lp = La + Lb

c Increment "Lp" to account for the order of differentiation.

      Lp = Lp + MXD

c Define the accumulated number of angular momentum functions <= Lp.

      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6

c Define the prefactor of the overlap distribution "P".

c Assign pointers to scratch space.
 
      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
 
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp
 
      if( i_left.ge.i_right )then
 
       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0
 
       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp
 
       stop
 
      end if
 
      if( DryRun )then

       MaxMem = i_left + (maxW0 - (i_right-1))
       NPP = NPA*NPB

      else

       call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &            Bxyz,Bprims,Bcoefs,NPB,NCB,
     &            GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)

      end if

      if (npp.eq.0) then
        if (O2I) call dfill(Nint,0.0d00,bO2I,1)
        if (KEI) call dfill(Nint,0.0d00,bKEI,1)
        if (NAI) call dfill(Nint,0.0d00,bNAI,1)
        return
      endif
c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_pf   = i_Ep     + 3*NPP*(MXD+1)*lprod
      i_left = i_pf     + 2*NPP - 1

      if( i_left.ge.i_right )then

       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0

       write(*,*) 'From the right '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'pf    :  ',i_pf
       write(*,*) 'From the left '
       write(*,*) 'ESp   :  ',i_ESp

       stop

      end if

      if( DryRun )then

       MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )

      else

       do 100 nd = 0,MXD
        call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &             nd,NPP,MXD,La+Li,Lb+Li)
  100  continue

      end if
       
c Compute the 2-center overlap integrals, <a|S|b>.

      if( O2I )then
        if( .not. DryRun )then
          call hf2oi(W0(i_Ep),bO2I,Nint,NPP,La,Lb,Li,canAB)
        end if
      end if

c Compute kinetic energy integrals, <a|T|b>.

      if( KEI )then

c Assign pointers to scratch space.

       i_Ti  = i_Ep + 3*NPP*(MXD+1)*lprod
       i_top = i_Ti + NPP - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1:  Insufficient scratch space.'
        write(*,*) '       needed    ',i_top
        write(*,*) '       allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'Ti    :  ',i_Ti

        stop

       end if

       if( DryRun )then

        MaxMem = max( MaxMem, i_top )

       else

        call hfkei(W0(i_ALPHAp),W0(i_Ep),bKEI,W0(i_Ti),
     &             Nint,NPP,La,Lb,Li,canAB)
       end if

      end if
       
c Compute nuclear attraction integrals, <a|V|b>.

      if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

       i_R0  = i_Ep  + 3*NPP*(MXD+1)*lprod
       i_IJK = i_R0  + NPP*Lp3*ncenters
       i_P   = i_IJK + (Lp+1)**3
       i_RS  = i_P   + NPP*3
       i_PC  = i_RS  + NPP
       i_ff  = i_PC  + NPP*3
       i_Rj  = i_ff  + NPP*2
       i_top = i_Rj  + NPP*(Lp+1)*Lp3 - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1:  Insufficient scratch space.'
        write(*,*) '       needed    ',i_top
        write(*,*) '       allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'R0    :  ',i_R0
        write(*,*) 'IJK   :  ',i_IJK
        write(*,*) 'P     :  ',i_P
        write(*,*) 'RS    :  ',i_RS
        write(*,*) 'PC    :  ',i_PC
        write(*,*) 'ff    :  ',i_ff
        write(*,*) 'Rj    :  ',i_Rj

        stop

       end if

       if( DryRun )then

        MaxMem = max( MaxMem, i_top )

       else

        call hf1mkrtmp(Axyz,Bxyz,Cxyz,zan,ncenters,
     &              W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &              W0(i_Rj),W0(i_R0),W0(i_R0),W0(i_IJK),
     &              NPP,Lp,Lp3,.FALSE.)

        call hfnaitmp(W0(i_Ep),W0(i_R0),W0(i_IJK),bNAI,
     &             Nint,NPP,La,Lb,Li,Lp,Lp3,canAB,ncenters)

       end if

      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
      Subroutine hfnaitmp(E,R0,IJK,Vab,Nint,NPP,
     &    La,Lb,Li,Lp,Lp3,canAB,ncenters)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Auxiliary Function Integrals & Index

      Dimension R0(ncenters,NPP,Lp3),IJK(0:Lp,0:Lp,0:Lp)

c--> Nuclear Attraction Integrals

      Dimension Vab(Nint)

c--> Scratch Space

      Dimension Nxyz(3)
c
c Compute the nuclear attraction integrals.
c
c     formula:
c           __
c           \    Ia,Ib    Ja,Jb    Ka,Kb
c     Vab = /  Ex     * Ey     * Ez     * R
c           --   Ip       Jp       Kp      Ip,Jp,Kp
c        Ip,Jp,Kp
c
c******************************************************************************

c Initialize the block of NAIs.

      do 10 nn = 1,Nint
       Vab(nn) = 0.D0
   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c      loop over centers
      do 23456 icic = 1,ncenters
c Loop over shell components.

      nn = 0

      do 50 ma = 1,La2

c Define the angular momentum indices for shell "A".

       call getNxyz(La,ma,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       if( canAB )then
        mb_limit = ma
       else
        mb_limit = Lb2
       end if

       do 40 mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

        call getNxyz(Lb,mb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        nn = nn + 1

        do 30 Ip = 0,Ia+Ib
        do 30 Jp = 0,Ja+Jb
        do 30 Kp = 0,Ka+Kb

         np = IJK(Ip,Jp,Kp)

         do 20 mp = 1,NPP
          Vab(nn) = Vab(nn) + (E(1,mp,Ip,Ia,Ib)*
     &                         E(2,mp,Jp,Ja,Jb)*
     &                         E(3,mp,Kp,Ka,Kb))*R0(icic,mp,np)
   20    continue

   30   continue

   40  continue

   50 continue

      write(6,*)'==================================================',
     &    '=============================='
      write(6,*)' for center ',icic,' modified V is '
      do inn = 1,nn
        if (abs(Vab(inn)).gt.1.0d-07) then
          write(6,*)' Vab(',inn,' ) =',vab(inn),' of ',nn
        endif
      enddo
      write(6,*)'==================================================',
     &    '=============================='
23456 continue
      end
      subroutine raktest_gc(rtdb)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer rtdb
c
      integer geom, basisseg, basisns, bases(2)
      integer mx1e, mxg, mxs1, mxs2, membuf, memscr
      integer h_buf, k_buf, h_scr, k_scr
      integer h_s_seg, k_s_seg, h_s_ns, k_s_ns
      integer h_diff, k_diff
      integer nbf_seg, ncont_seg, nbf_ns, ncont_ns
      integer icount, niter
      double precision norm, tov_seg, tov_ns
      double precision thresh_norm
      logical status
      logical FF, FT
      parameter (FF=.false.,FT=.true.)
c
      logical int_normalize
      external int_normalize
c
      thresh_norm = 1.0d-06
c
      write(luout,*)' ============================================ '
      write(luout,*)' ================            ================ '
      write(luout,*)' ================ raktest_gc ================ '
      write(luout,*)' ================            ================ '
      write(luout,*)' ============================================ '
c
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('raktest_gc: geom_create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &      call errquit('raktest_gc: geom_rtdb_load failed',911,
     &       RTDB_ERR)
c
      if(.not.bas_create(basisseg,'aobasisseg'))
     &      call errquit('raktest_gc: bas_create choked',911, BASIS_ERR)
      if(.not.
     &      bas_rtdb_load(rtdb,geom,basisseg,'aobasisseg'))
     &      call errquit('raktest_gc: bas_rtdb_load choked',911,
     &       RTDB_ERR)
c
      if(.not.bas_create(basisns,'aobasisns'))
     &      call errquit('raktest_gc: bas_create choked',911, BASIS_ERR)
      if(.not.
     &      bas_rtdb_load(rtdb,geom,basisns,'aobasisns'))
     &      call errquit('raktest_gc: bas_rtdb_load choked',911,
     &       RTDB_ERR)
c
      bases(1) = basisseg
      bases(2) = basisns
c
      status = geom_print(geom)
      status = bas_print(basisseg)
      status = gbs_map_print(basisseg)
      status = bas_print(basisns)
      status = gbs_map_print(basisns)
c
c... get memory requirements for segmented basis
      call int_init(rtdb,1,basisseg)
      call int_mem(mx1e, mxg, mxs1, mxs2)
      write(luout,*)' seg: one electron buffer size        :',mx1e
      write(luout,*)' seg: two electron buffer size        :',mxg
      write(luout,*)' seg: one electron scratch buffer size:',mxs1
      write(luout,*)' seg: two electron scratch buffer size:',mxs2
      call int_terminate()
c
c... get memory requirements for non-segmented basis
      call int_init(rtdb,1,basisns)
      call int_mem(mx1e, mxg, mxs1, mxs2)
      write(luout,*)'  ns: one electron buffer size        :',mx1e
      write(luout,*)'  ns: two electron buffer size        :',mxg
      write(luout,*)'  ns: one electron scratch buffer size:',mxs1
      write(luout,*)'  ns: two electron scratch buffer size:',mxs2
      call int_terminate()
      call int_init(rtdb,2,bases)
      if(.not.int_normalize(rtdb,basisseg)) call errquit
     &      ('raktest_gc: int_normalize seg failed',911, INT_ERR)
      if(.not.int_normalize(rtdb,basisns)) call errquit
     &      ('raktest_gc: int_normalize ns failed',911, INT_ERR)
      write(6,*)' after normalization'
      status = bas_print(basisseg)
      status = gbs_map_print(basisseg)
      status = bas_print(basisns)
      status = gbs_map_print(basisns)
c
      call int_mem(mx1e, mxg, mxs1, mxs2)
      write(luout,*)' both: one electron buffer size        :',mx1e
      write(luout,*)' both: two electron buffer size        :',mxg
      write(luout,*)' both: one electron scratch buffer size:',mxs1
      write(luout,*)' both: two electron scratch buffer size:',mxs2
c
      membuf = max(mx1e,mxg,10000)
      memscr = max(mxs1,mxs2,100000)
c
      if(.not.ma_push_get(mt_dbl,membuf,'int buff',h_buf,k_buf))
     &      call errquit('raktest_gc: ma failed 1',911, MA_ERR)
      if(.not.ma_push_get(mt_dbl,memscr,'int scr',h_scr,k_scr))
     &      call errquit('raktest_gc: ma failed 2',911, MA_ERR)
c
      if(.not.bas_numcont(basisseg,ncont_seg))
     &      call errquit('raktest_gc: bas_numcont failed',911,
     &       BASIS_ERR)
      write(luout,*)' ao basis seg number of contractions    :',
     &      ncont_seg
      if(.not.bas_numbf(basisseg,nbf_seg))
     &      call errquit('raktest_gc: bas_numbf failed',911,
     &       BASIS_ERR)
      write(luout,*)' ao basis seg number of basis functions :',
     &      nbf_seg
      if(.not.bas_numcont(basisns,ncont_ns))
     &      call errquit('raktest_gc: bas_numcont failed',911,
     &       BASIS_ERR)
      write(luout,*)' ao basis noseg number of contractions    :',
     &      ncont_ns
      if(.not.bas_numbf(basisns,nbf_ns))
     &      call errquit('raktest_gc: bas_numbf failed',911,
     &       BASIS_ERR)
      write(luout,*)' ao basis noseg number of basis functions :',
     &      nbf_ns
c
      if (nbf_ns.ne.nbf_seg) call errquit
     &      ('raktest_gc: nbf_seg.ne.nbf_ns',(nbf_ns-nbf_seg),
     &       UNKNOWN_ERR)
c
      if(.not.ma_push_get(mt_dbl,(nbf_seg*nbf_seg),'overlap seg',
     &      h_s_seg,k_s_seg))
     &      call errquit('raktest_gc: ma failed 3',911, MA_ERR)
      if(.not.ma_push_get(mt_dbl,(nbf_ns*nbf_ns),'overlap noseg',
     &      h_s_ns,k_s_ns))
     &      call errquit('raktest_gc: ma failed 4',911, MA_ERR)
      if(.not.ma_push_get(mt_dbl,(nbf_ns*nbf_ns),'difference matrix',
     &      h_diff,k_diff))
     &      call errquit('raktest_gc: ma failed 4',911, MA_ERR)

c
c set niter number of iterations
c
      niter = 20
      write(6,'(/,a,i5)')
     &    'number of instances computed for each seg/noseg set:',niter
c
c zero ma segments

      write(6,*)' dfill k_diff',k_diff
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_diff),1)
      write(6,*)' dfill k_ns',k_s_ns
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns),1)
      write(6,*)' dfill k_seg',k_s_seg
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg),1)
      write(6,*)' dfill k_scr',k_scr
      call dfill(memscr,           0.0d00,dbl_mb(k_scr),1)
      write(6,*)' dfill k_buf',k_buf
      call dfill(membuf,           0.0d00,dbl_mb(k_buf),1)
c
c overlap check
      tov_seg = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisseg,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_seg,ncont_seg,dbl_mb(k_s_seg),
     &        FT,FF,FF,'  segment overlap')
      enddo
      tov_seg = util_cpusec() - tov_seg
      tov_ns  = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisns,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_ns,ncont_ns,dbl_mb(k_s_ns),
     &        FT,FF,FF,'nosegment overlap')

      enddo
      tov_ns = util_cpusec() - tov_ns
c
      call dcopy((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_diff),1)
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_diff),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_diff),1,dbl_mb(k_diff),1)
      if (norm.gt.thresh_norm)
     &    call print_diff_gc(nbf_ns,
     &    dbl_mb(k_diff),
     &    dbl_mb(k_s_seg),
     &    dbl_mb(k_s_ns),
     &    thresh_norm,
     &    'overlap')
      write(luout,*)' '
      write(luout,*)' time for segmented overlap         :',tov_seg
      write(luout,*)' time for non-segmented overlap     :',tov_ns
      write(luout,'(a,f10.2)')
     &    ' % speedup                          :',
     &    (tov_seg-tov_ns)/tov_seg*100.0d00
      write(luout,*)'raktest_gc: overlap difference norm :',norm
      write(luout,*)' '
c
c zero ma segments
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_diff),1)
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns),1)
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg),1)
      call dfill(memscr,           0.0d00,dbl_mb(k_scr),1)
      call dfill(membuf,           0.0d00,dbl_mb(k_buf),1)
c kinetic energy check
      tov_seg = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisseg,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_seg,ncont_seg,dbl_mb(k_s_seg),
     &        FF,FT,FF,'  segment kinetic')
      enddo
      tov_seg = util_cpusec() - tov_seg
      tov_ns  = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisns,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_ns,ncont_ns,dbl_mb(k_s_ns),
     &        FF,FT,FF,'nosegment kinetic')

      enddo
      tov_ns = util_cpusec() - tov_ns
c
      call dcopy((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_diff),1)
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_diff),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_diff),1,dbl_mb(k_diff),1)
      if (norm.gt.thresh_norm)
     &    call print_diff_gc(nbf_ns,
     &    dbl_mb(k_diff),
     &    dbl_mb(k_s_seg),
     &    dbl_mb(k_s_ns),
     &    thresh_norm,
     &    'kinetic')
      write(luout,*)' '
      write(luout,*)' time for segmented kinetic         :',tov_seg
      write(luout,*)' time for non-segmented kinetic     :',tov_ns
      write(luout,'(a,f10.2)')
     &    ' % speedup                          :',
     &    (tov_seg-tov_ns)/tov_seg*100.0d00
      write(luout,*)'raktest_gc: kinetic difference norm :',norm
      write(luout,*)' '
c
c zero ma segments
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_diff),1)
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns),1)
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg),1)
      call dfill(memscr,           0.0d00,dbl_mb(k_scr),1)
      call dfill(membuf,           0.0d00,dbl_mb(k_buf),1)
c potential check
      tov_seg = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisseg,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_seg,ncont_seg,dbl_mb(k_s_seg),
     &        FF,FF,FT,'  segment potential')
      enddo
      tov_seg = util_cpusec() - tov_seg
      tov_ns  = util_cpusec()
      do icount = 1,niter
        call raktest_bs_gc(basisns,
     &        memscr,dbl_mb(k_scr),
     &        membuf,dbl_mb(k_buf),
     &        nbf_ns,ncont_ns,dbl_mb(k_s_ns),
     &        FF,FF,FT,'nosegment potential')

      enddo
      tov_ns = util_cpusec() - tov_ns
c
      call dcopy((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_diff),1)
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_diff),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_diff),1,dbl_mb(k_diff),1)
      if (norm.gt.thresh_norm)
     &    call print_diff_gc(nbf_ns,
     &    dbl_mb(k_diff),
     &    dbl_mb(k_s_seg),
     &    dbl_mb(k_s_ns),
     &    thresh_norm,
     &    'potential')
      write(luout,*)' '
      write(luout,*)' time for segmented potential         :',tov_seg
      write(luout,*)' time for non-segmented potential     :',tov_ns
      write(luout,'(a,f10.2)')
     &    ' % speedup                          :',
     &    (tov_seg-tov_ns)/tov_seg*100.0d00
      write(luout,*)'raktest_gc: potential difference norm :',norm
      write(luout,*)' '
c
      status =              ma_pop_stack(h_diff)
      status = status .and. ma_pop_stack(h_s_ns)
      status = status .and. ma_pop_stack(h_s_seg)
      status = status .and. ma_pop_stack(h_scr)
      status = status .and. ma_pop_stack(h_buf)
      if (.not.status) call errquit('raktest_gc: pop error',911, MA_ERR)
c
      if (.not.bas_destroy(basisseg)) call errquit
     &    ('raktest_gc: bas_destroy failed ?',911, BASIS_ERR)
      if (.not.bas_destroy(basisns)) call errquit
     &    ('raktest_gc: bas_destroy failed ?',911, BASIS_ERR)
      if (.not.geom_destroy(geom)) call errquit
     &    ('raktest_gc: geom_destroy failed ?',911, GEOM_ERR)
c
      call int_terminate()
c
      end
      subroutine raktest_bs_gc(basis,nscr,scr,nbuf,buf,nbf,ncont,S,
     &      OV,KE,PE,msg)
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c
#include "bas.fh"
c
      integer basis,nscr,nbuf,nbf,ncont
      double precision scr(nscr), buf(nbuf)
      double precision S(nbf,nbf)
      logical OV,KE,PE
      character*(*) msg
c
      double precision val
      integer nint
      integer ish, ilo, ihi, ibf
      integer jsh, jlo, jhi, jbf
      integer icount
c
c      check validity
      icount = 0   
      if (OV) icount = icount + 1
      if (KE) icount = icount + 1
      if (PE) icount = icount + 1
      if (icount.eq.0) then
        write(luout,*)' no integral set defined '
        call dfill ((nbf*nbf),0.0d00,S,1)
        return
      elseif(icount.ne.1) then
        write(luout,*)' OV =',OV
        write(luout,*)' KE =',KE
        write(luout,*)' PE =',PE
        stop 'error'
      endif
c
      do ish = 1,ncont
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) call errquit
     &          ('raktest_gc: cn2bfr failed',911, BASIS_ERR)
*          write(luout,*)'ish = ',ish,' bfr =',ilo,ihi
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     &          ('raktest_gc: cn2bfr failed',911, BASIS_ERR)
*          write(luout,*)'jsh = ',jsh,' bfr =',jlo,jhi
          nint = (ihi-ilo+1)*(jhi-jlo+1)
*          write(luout,*)'number of integrals = ',nint
          if (OV) then
            call int_1eov(basis,ish,basis,jsh,
     &          nscr,scr,nbuf,buf)
          elseif (KE) then
            call int_1eke(basis,ish,basis,jsh,
     &          nscr,scr,nbuf,buf)
          elseif (PE) then
            call int_1epe(basis,ish,basis,jsh,
     &          nscr,scr,nbuf,buf)
          else
            write(luout,*)' no integral set defined '
            call dfill ((nbf*nbf),0.0d00,S,1)
            return
          endif
          icount = 0
          nint = 0
          do ibf = ilo, ihi
            do jbf = jlo, jhi
              nint = nint+1
              val = buf(nint)
              s(ibf,jbf) = val
              s(jbf,ibf) = val
              if (abs(val).gt.1.0d-07) then
                icount = icount + 1
*                write(69,10000)ibf,jbf,val,msg,' ints ',icount
*                write(luout,10000)ibf,jbf,val,msg,' ints ',icount
              endif
            enddo
          enddo
        enddo
      enddo
c
*      write(luout,*)' raktest_gc matrix ',msg
*      call output(S,1,nbf,1,nbf,nbf,nbf,1)
c
10000 format(i5,i5,1pd20.10,1x, a,a,i12)
c
      end
      subroutine print_diff_gc(nbf,diff,seg,ns,ths,msg)
      implicit none
#include "stdio.fh"
      integer nbf
      double precision diff(nbf,nbf)
      double precision seg(nbf,nbf)
      double precision ns(nbf,nbf)
      double precision ths
      character*(*) msg
c
      integer i,j
c
      write(luout,*)' differeces from ',msg
c
      do i=1,nbf
        do j=1,nbf
          if (abs(diff(i,j)).gt.ths) then
            write(6,10000)i,j,diff(i,j),seg(i,j),ns(i,j),msg
          endif
        enddo
      enddo
10000 format('<',i3,'|',i3,'> <diff',1pd20.10,'>  <seg',1pd20.10,
     &    '>  <noseg',1pd20.10,'>',1x,a)
      end
      subroutine print_diff_vec(n,a,b,ths,msg)
      implicit none
      integer n
      double precision a(n), b(n)
      double precision ths
      character*(*) msg
c
      integer i, count, nza, nzb
      double precision diff
c
      write(6,*)' print_diff_vec <<',msg,'>>'
      write(6,*)' print_diff_vec threshold:',ths
c
      nza = 0
      nzb = 0
      count = 0
      do i = 1,n
        if (a(i).ne.0.0d00) nza = nza + 1
        if (b(i).ne.0.0d00) nzb = nzb + 1
        diff = a(i) - b(i)
        if (abs(diff).gt.ths) then
          count = count + 1
          write(6,'(1x,i8,a,d12.6,a,d12.6,a,d12.6)')
     &        i,' th element delta:',diff,
     &        ' a:',a(i),' b:',b(i)
        endif
      enddo
      write(6,*)' print_diff_vec: number of different elements    :',
     &    count
      write(6,*)' print_diff_vec: number of nonzero elements in a :',
     &    nza
      write(6,*)' print_diff_vec: number of nonzero elements in b :',
     &    nzb
      end
      subroutine raktest_test9(rtdb)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
c      
      integer rtdb
c
      integer basis, geom, count
      integer i, j, k, l
      integer ish, jsh, ksh, lsh
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer k_buf, h_buf, k_scr, h_scr
      integer max1e, max2e, mscr1, mscr2, m_buf, m_scr
      integer inshell(4)
      logical status
      logical int_normalize
      external int_normalize
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911, RTDB_ERR)
c
      if (.not.bas_create(basis,'ao sp_basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao sp_basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911, BASIS_ERR)
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
c
      if (.not.geom_print(geom)) stop ' print error'
      if (.not.bas_print(basis)) stop ' print error'
c
      call int_init(rtdb,1,basis)
      call int_mem(max1e,max2e,mscr1,mscr2)
      m_buf = max(max1e*2,max2e*2)
      m_scr = max(mscr1*2,mscr2)
      m_buf = m_buf + (m_buf*110)/100
      m_scr = m_scr + (m_scr*110)/100
c
      if (.not.ma_push_get(mt_dbl,m_scr,'scr for 1e',h_scr,k_scr))
     &      call errquit('ma scr failed',911, MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,m_buf,'buf for 1e',h_buf,k_buf))
     &      call errquit('ma buf failed',911, MA_ERR)
c
      ish = 1
      jsh = 1
      ksh = 1
      lsh = 1
      if (rtdb_get(rtdb,'rak9:shells',mt_int,4,inshell)) then
        ish = inshell(1)
        jsh = inshell(2)
        ksh = inshell(3)
        lsh = inshell(4)
      else
        write(6,*)'rak9:shells not set on rtdb'
      endif
c
      write(6,*)' rak9 for shells ',ish,jsh,ksh,lsh
c
      call int_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &    m_scr,dbl_mb(k_scr),m_buf,dbl_mb(k_buf))
c
      if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &    stop 'cn2bfr error i'
      if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &    stop 'cn2bfr error j'
      if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &    stop 'cn2bfr error k'
      if (.not.bas_cn2bfr(basis,lsh,llo,lhi))
     &    stop 'cn2bfr error l'
      count = 0
      do i=ilo,ihi
        do j=jlo,jhi
          do k=klo,khi
            do l=llo,lhi
              write(6,*)i,j,k,l,dbl_mb(k_buf+count)
              write(69,*)i,j,k,l,dbl_mb(k_buf+count)
              count = count + 1
            enddo
          enddo
        enddo
      enddo
c
      status = .true.
      status = status.and.ma_pop_stack(h_buf)
      status = status.and.ma_pop_stack(h_scr)
c
      if (.not.status) call errquit('ma pop fail',911, MA_ERR)
      call int_terminate()
c
      if (.not.bas_destroy(basis)) call errquit
     &    ('9: bas_destroy failed ?',911, BASIS_ERR)
      if (.not.geom_destroy(geom)) call errquit
     &    ('9: geom_destroy failed ?',911, GEOM_ERR)
c
      end
      subroutine raktest_ecp(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
      integer rtdb
c
      integer geom, basis, ecpid
c
      logical int_normalize, int_ecp_init
      external int_normalize, int_ecp_init
c
*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 1?'
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911,
     &       RTDB_ERR)
      if (.not.geom_print(geom)) stop ' print error'
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
      if (.not.bas_print(basis)) stop ' print error'
*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 2?'
c
      if (.not.bas_get_ecp_handle(basis,ecpid)) stop 'get/ecp/handle'
*ecp:      if (.not.bas_create(ecpid,'ecp basis')) call errquit
*ecp:     &      ('bas_create failed',911, BASIS_ERR)
*ecp:      if (.not.bas_rtdb_load(rtdb,geom,ecpid,'ecp basis'))
*ecp:     &    call errquit
*ecp:     &      ('bas_rtdb_load failed',911, RTDB_ERR)
      if (.not.bas_print(ecpid)) stop ' print error'
c
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      if (.not.bas_print(basis)) stop ' print error'
      if (.not.gbs_map_print(basis)) stop ' gbs map print error'
      if (.not.bas_print(ecpid)) stop ' print error'
      if (.not.gbs_map_print(basis)) stop ' gbs map print error'
c
      if (.not.int_ecp_init(ecpid,0,0)) stop ' error in int_ecp_init'
      call int_ecp_terminate()
c
      if (.not.bas_destroy(basis)) stop ' bas_dest error 1'
      if (.not.bas_destroy(ecpid)) stop ' bas_dest error 2'
      if (.not.geom_destroy(geom)) stop ' geom_dest error 1'
      if (.not.bas_version()) stop ' bas_version error'
c
      end
      subroutine raktest_bug(rtdb)
      implicit none
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"      
#include "geomP.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      logical int_normalize
      external int_normalize
c
      integer rtdb
c
      integer basao, basecp
      integer indao, indecp
      integer geom
      integer maxg, maxs, maxg1, maxs1, maxg2, maxs2
      integer hao, kao, hscr, kscr, h3p, k3p
      logical ignore
c
      integer icont, ucont
      integer atype, aprim, ie_a, ic_a, a_cent
      integer btype, bprim, ie_b, ic_b, b_cent
      integer ctype, cprim, ie_c, ic_c, c_cent
      integer ac_prim, bc_prim
      integer ir_c
      integer i,ir,icnt,j
      integer l2a, l2b, nint
      logical FF, FT
      double precision xyza(3),xyzb(3),xyzc(3),xyzall(3,3)
      double precision zan(3)
c
      integer maxprim
      parameter (maxprim=100)
      integer rc(maxprim)
      double precision ea(maxprim),ca(maxprim)
      double precision eb(maxprim),cb(maxprim)
      double precision ec(maxprim),cc(maxprim)
      double precision eac(maxprim),cac(maxprim)
      double precision ebc(maxprim),cbc(maxprim)
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      ignore = geom_create(geom,'geometry')
      ignore = ignore.and.geom_rtdb_load(rtdb,geom,'geometry')
      ignore = ignore.and.bas_create(basao,'ao basis')
      ignore = ignore.and.bas_rtdb_load(rtdb,geom,basao,'ao basis')
      ignore = ignore.and.int_normalize(rtdb,basao)
      ignore = ignore.and.geom_print(geom)
      ignore = ignore.and.bas_print(basao)
      ignore = ignore.and.gbs_map_print(basao)
      ignore = ignore.and.bas_get_ecp_handle(basao,basecp)
      ignore = ignore.and.ecp_print(basecp)
      if (.not.ignore) stop ' something failed'
c
      call int_init(rtdb,1,basao)
      call int_mem_1e(maxg1,maxs1)
      call int_mem_2e4c(maxg2,maxs2)
      maxg = max(maxg1,maxg2)
      maxs = max(maxs1,maxs2)
c
      ignore =
     &    ma_push_get(mt_dbl,maxg,'aobuf',hao,kao)
      ignore =
     &    ma_push_get(mt_dbl,maxg,'3pbuf',h3p,k3p)
      ignore=ignore.and.
     &    ma_push_get(mt_dbl,maxs,'scr',hscr,kscr)

      indao  = basao  + basis_handle_offset
      indecp = basecp + basis_handle_offset
c
      call dfill(maxprim,0.0d00,ea,1)
      call dfill(maxprim,0.0d00,ca,1)
      call dfill(maxprim,0.0d00,eb,1)
      call dfill(maxprim,0.0d00,cb,1)
      call dfill(maxprim,0.0d00,ec,1)
      call dfill(maxprim,0.0d00,cc,1)
      call dfill(maxprim,0.0d00,eac,1)
      call dfill(maxprim,0.0d00,cac,1)
      call dfill(maxprim,0.0d00,ebc,1)
      call dfill(maxprim,0.0d00,cbc,1)
      call ifill(maxprim,0,rc,1)
      call dfill(3,0.0d00,xyza,1)
      call dfill(3,0.0d00,xyzb,1)
      call dfill(3,0.0d00,xyzc,1)
      call dfill(3*3,0.0d00,xyzall,1)
      call dfill(3,-1.0d00,zan,1)
c
c   for cont1|ecp|cont1 it fails
c      
c fill a 
      icont = 1
      ucont = (sf_ibs_cn2ucn(icont,indao))
      atype = infbs_cont(Cont_Type,ucont,indao)
      aprim = infbs_cont(Cont_Nprim,ucont,indao)
      if (infbs_cont(Cont_Ngen,ucont,indao).ne.1)
     &    stop ' a ngen wrong '
      ie_a  = infbs_cont(Cont_Iexp,ucont,indao)
      ic_a  = infbs_cont(Cont_Icfp,ucont,indao)
      a_cent = sf_ibs_cn2ce(icont,indao)
      call dcopy(3,coords(1,a_cent,geom),1,xyza,1)
      if (aprim.le.maxprim) then
        call dcopy(aprim,dbl_mb(mb_exndcf(ie_a,indao)),1,ea,1)
        call dcopy(aprim,dbl_mb(mb_exndcf(ic_a,indao)),1,ca,1)
        write(6,*)' ea '
        call output(ea,1,aprim,1,1,aprim,1,1)
        write(6,*)' ca '
        call output(ca,1,aprim,1,1,aprim,1,1)
      else
        stop ' aprim error '
      endif
c fill b
      icont = 1
      ucont = (sf_ibs_cn2ucn(icont,indao))
      btype = infbs_cont(Cont_Type,ucont,indao)
      bprim = infbs_cont(Cont_Nprim,ucont,indao)
      if (infbs_cont(Cont_Ngen,ucont,indao).ne.1)
     &    stop ' b ngen wrong '
      ie_b  = infbs_cont(Cont_Iexp,ucont,indao)
      ic_b  = infbs_cont(Cont_Icfp,ucont,indao)
      b_cent = sf_ibs_cn2ce(icont,indao)
      call dcopy(3,coords(1,b_cent,geom),1,xyzb,1)
      if (bprim.le.maxprim) then
        call dcopy(bprim,dbl_mb(mb_exndcf(ie_b,indao)),1,eb,1)
        call dcopy(bprim,dbl_mb(mb_exndcf(ic_b,indao)),1,cb,1)
        write(6,*)' eb '
        call output(eb,1,bprim,1,1,bprim,1,1)
        write(6,*)' cb '
        call output(cb,1,bprim,1,1,bprim,1,1)
      else
        stop ' bprim error '
      endif
c fill c
      icont = 1
      ucont = sf_ibs_cn2ucn(icont,indecp)
      ctype = infbs_cont(Cont_Type,ucont,indecp)
      cprim = infbs_cont(Cont_Nprim,ucont,indecp)
      if (infbs_cont(Cont_Ngen,ucont,indecp).ne.1)
     &    stop 'c ngen wrong'
      ie_c =  infbs_cont(Cont_iexp, ucont,indecp)
      ic_c =  infbs_cont(Cont_icfp, ucont,indecp)
      ir_c =  infbs_cont(Cont_irexp,ucont,indecp)
      c_cent = sf_ibs_cn2ce(icont,indecp)
      call dcopy(3,coords(1,c_cent,geom),1,xyzc,1)
c
      icnt = 0
      do i = 1,cprim
        ir = int(sf_exndcf((ir_c-1+i),indecp) + 0.00001d0)
        if (ir.eq.1) then
          if ((icnt+1).le.maxprim) then
            icnt = icnt + 1
            ec(icnt) = sf_exndcf((ie_c-1+i),indecp)
            cc(icnt) = sf_exndcf((ic_c-1+i),indecp)
            rc(icnt) = ir
          else
            stop 'cprim > maxprim'
          endif
        endif
      enddo
      cprim = icnt
      write(6,*)' ec '
      call output(ec,1,cprim,1,1,cprim,1,1)
      write(6,*)' cc '
      call output(cc,1,cprim,1,1,cprim,1,1)
      write(6,*)' rc '
      do i = 1,cprim
        write(6,*)i,rc(i)
      enddo
c
c copy all coords
      call dcopy(3,xyza,1,xyzall(1,1),1)
      call dcopy(3,xyzb,1,xyzall(1,2),1)
      call dcopy(3,xyzc,1,xyzall(1,3),1)
c print coords
      write(6,*)' xyza ',xyza
      write(6,*)' xyzb ',xyzb
      write(6,*)' xyzc ',xyzc
      write(6,*)' xyz all'
      call output(xyzall,1,3,1,3,3,3,1)
c
c
      if (a_cent.eq.b_cent.and.a_cent.eq.c_cent) then
        write(6,*)' one center case'
      else
        stop ' not one center case'
      endif
c
c form ac exponents
      icnt = 0
      do i = 1,aprim
        do j = 1,cprim
          icnt = icnt + 1
          eac(icnt) = ea(i)+ec(j)
          cac(icnt) = ca(i)*cc(j)
        enddo
      enddo
      ac_prim = icnt
      write(6,*)' eac'
      call output(eac,1,ac_prim,1,1,ac_prim,1,1)
      write(6,*)' cac'
      call output(cac,1,ac_prim,1,1,ac_prim,1,1)
c form bc exponents
      icnt = 0
      do i = 1,bprim
        do j = 1,cprim
          icnt = icnt + 1
          ebc(icnt) = eb(i)+ec(j)
          cbc(icnt) = cb(i)*cc(j)
        enddo
      enddo
      bc_prim = icnt
      write(6,*)' ebc'
      call output(ebc,1,bc_prim,1,1,bc_prim,1,1)
      write(6,*)' cbc'
      call output(cbc,1,bc_prim,1,1,bc_prim,1,1)
c
      FF = .false.
      FT = .true.
      call dfill(maxg,0.0d00,dbl_mb(kao),1)
      call dfill(maxs,0.0d00,dbl_mb(kscr),1)
      call hf1(
     &    xyza,eac,cac,ac_prim,1,atype,
     &    xyzb,eb,cb,bprim,1,btype,
     &    xyza,zan,1,
     &    dbl_mb(kao),dbl_mb(kao),dbl_mb(kao),maxg,FF,FF,FT,FF,
     &    FF,dbl_mb(kscr),maxs)
c
      call dfill(maxg,0.0d00,dbl_mb(k3p),1)
      call dfill(maxs,0.0d00,dbl_mb(kscr),1)
      call hf3pot(
     &    xyza,ea,ca,aprim,1,atype,
     &    xyzb,eb,cb,bprim,1,btype,
     &    xyzc,ec,cc,cprim,1,0,
     &    dbl_mb(k3p),maxg,nint,
     &    FF,dbl_mb(kscr),maxs)
c
      l2a = (atype+1)*(atype+2)/2
      l2b = (btype+1)*(btype+2)/2
      write(6,*)'from hf1'
      call intintp(dbl_mb(kao),l2a,l2b,'from hf1')
      write(6,*)'from hf3pot'
      call intintp(dbl_mb(k3p),l2a,l2b,'from hf3pot')
c      
      call int_terminate()
      ignore = bas_destroy(basao)
c
      ignore = ma_pop_stack(hscr)  
      ignore = ma_pop_stack(h3p)  
      ignore = ma_pop_stack(hao)  
c
      end
      subroutine intintp(z,r,c,msg)
      implicit none
      integer r,c
      double precision z(r,c)
      character*(*) msg
c
      integer ir, ic
      do ir = 1,r
        do ic = 1,c
          if (abs(z(ir,ic)).gt.1.0d-10) then
            write(6,10000)ir,ic,z(ir,ic),msg
          endif
        enddo
      enddo
10000 format(1x,'(',i5,',',i5,')',1x,1pd20.10,a)
      end
      subroutine raktest_geomprt(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
      integer rtdb
c
      integer geom
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomprt:create error',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &    call errquit('raktest_geomprt:load error',911, RTDB_ERR)
      if (.not.geom_print_distances(geom))
     &    call errquit('raktest_geomprt:print_distance error',911,
     &       GEOM_ERR)
      if (.not.geom_print_angles(geom))
     &    call errquit('raktest_geomprt:print_angles error',911,
     &       GEOM_ERR)
      if (.not.geom_print_dihedrals(geom))
     &    call errquit('raktest_geomprt:print_dihedrals error',911,
     &       GEOM_ERR)
      if (.not.geom_destroy(geom))
     &    call errquit('raktest_geomprt:destory',911,
     &       GEOM_ERR)
      end
      subroutine raktest_3cd(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb
c
      integer geom, basis
      integer nat, nbf, lmax, bufsz, dbufsz, nbftrip, szscr
      integer h_bovp,h_bovm,h_dbuf,h_fd3ov,h_3ov,h_scr,h_c,h_cp,h_cm
      integer k_bovp,k_bovm,k_dbuf,k_fd3ov,k_3ov,k_scr,k_c,k_cp,k_cm
      integer h_fdp3ov, h_fdm3ov
      integer k_fdp3ov, k_fdm3ov
      integer stackleft
c
      logical int_normalize
      external int_normalize
c
*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 1?'
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911,
     &       RTDB_ERR)
      if (.not.geom_print(geom)) stop ' print error'
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
      if (.not.bas_print(basis)) stop ' print error'
      if (.not.rtdb_print(rtdb,.false.)) stop ' rtdb_p 2?'
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      call int_init(rtdb,1,basis)
c
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      if (.not.bas_numbf(basis,nbf)) stop 'bas_numbf fe'
      if (.not.bas_high_angular(basis,Lmax)) stop 'bas_ha fe'
c
      bufsz = (Lmax+1)*(Lmax+2)/2
      bufsz = bufsz*bufsz*bufsz
      dbufsz = 9*bufsz
      nbftrip = nbf*nbf*nbf
      szscr  = 60000
      if (MA_verify_allocator_stuff()) write(6,*)' maok (0)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <0> stack left ',stackleft,
     &    ' next:',bufsz
      if (.not.ma_push_get(mt_dbl,bufsz,' buf ovlap +',
     &    h_bovp,k_bovp)) stop ' ma_pg 1 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (1)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <1> stack left ',stackleft,
     &    ' next:',bufsz
      if (.not.ma_push_get(mt_dbl,bufsz,' buf ovlap -',
     &    h_bovm,k_bovm)) stop ' ma_pg 2 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (2)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <2> stack left ',stackleft,
     &    ' next:',dbufsz
      if (.not.ma_push_get(mt_dbl,dbufsz,'deriv buf ',
     &    h_dbuf,k_dbuf)) stop ' ma_pg 3 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (3)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <3> stack left ',stackleft,
     &    ' next:',(nbftrip*nat*3)
      if (.not.ma_push_get(mt_dbl,(nbftrip*nat*3),' fd 3ov matrix ',
     &    h_fd3ov, k_fd3ov)) stop ' ma_pg 4 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (4)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <4> stack left ',stackleft,
     &    ' next:',(nbftrip*nat*3)
      if (.not.ma_push_get(mt_dbl,(nbftrip*nat*3),' fd + 3ov matrix ',
     &    h_fdp3ov, k_fdp3ov)) stop ' ma_pg 5 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (5)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <5> stack left ',stackleft,
     &    ' next:',(nbftrip*nat*3)
      if (.not.ma_push_get(mt_dbl,(nbftrip*nat*3),' fd - 3ov matrix ',
     &    h_fdm3ov, k_fdm3ov)) stop ' ma_pg 6 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (6)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <6> stack left ',stackleft,
     &    ' next:',(nbftrip*nat*3)
      if (.not.ma_push_get(mt_dbl,(nbftrip*nat*3),' 3ov matrix ',
     &    h_3ov, k_3ov)) stop ' ma_pg 7 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (7)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <7> stack left ',stackleft,
     &    ' next:',szscr
      if (.not.ma_push_get(mt_dbl,szscr,' scratch buffer',
     &    h_scr, k_scr)) stop ' ma_pg 8 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (8)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <8> stack left ',stackleft,
     &    ' next:',(3*nat)
      if (.not.ma_push_get(mt_dbl,(3*nat),' coords ',
     &    h_c, k_c)) stop ' ma_pg 9 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (9)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <9> stack left ',stackleft,
     &    ' next:',(3*nat)
      if (.not.ma_push_get(mt_dbl,(3*nat),' coords +',
     &    h_cp, k_cp)) stop ' ma_pg 10 fail '
      if (MA_verify_allocator_stuff()) write(6,*)' maok (10)'
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <10> stack left ',stackleft,
     &    ' next:',(3*nat)
      if (MA_verify_allocator_stuff()) write(6,*)' maok (11)'
      if (.not.ma_push_get(mt_dbl,(3*nat),' coords -',
     &    h_cm, k_cm)) stop ' ma_pg 11 fail '
      stackleft = ma_inquire_stack(mt_dbl)
      write(6,*)' <11> stack left ',stackleft
c
      call raktest_3cd_1(geom,basis,nbf,nat,szscr,bufsz,dbufsz,
     &    dbl_mb(k_bovp),  dbl_mb(k_bovm),   dbl_mb(k_dbuf),
     &    dbl_mb(k_fd3ov), dbl_mb(k_fdp3ov), dbl_mb(k_fdm3ov),
     &    dbl_mb(k_3ov),  dbl_mb(k_scr),
     &    dbl_mb(k_c), dbl_mb(k_cp), dbl_mb(k_cm))
c
      if (.not.ma_pop_stack(h_cm)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_cp)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_c)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_scr)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_3ov)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_fdm3ov)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_fdp3ov)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_fd3ov)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_dbuf)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_bovm)) stop ' ma pop fail'
      if (.not.ma_pop_stack(h_bovp)) stop ' ma pop fail'
c
      call int_terminate
      if (.not.bas_destroy(basis)) stop ' bas_dest error 1'
      if (.not.geom_destroy(geom)) stop ' geom_dest error 1'
      if (.not.bas_version()) stop ' bas_version error'
      end
      subroutine raktest_3cd_1(geom,basisin,nbf,nat,szscr,bufsz,dbufsz,
     &    bufp, bufm, dbuf, fdo, fdp, fdm, o, scr, c, cp, cm)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
      integer geom
      integer basisin
      integer nbf
      integer nat
      integer szscr
      integer bufsz
      integer dbufsz
      double precision bufp(bufsz), bufm(bufsz), dbuf(dbufsz)
      double precision fdo(nbf,nbf,nbf,3,nat)
      double precision fdp(nbf,nbf,nbf,3,nat)
      double precision fdm(nbf,nbf,nbf,3,nat)
      double precision o(nbf,nbf,nbf,3,nat)
      double precision scr(szscr)
      double precision c(3,nat), cp(3,nat), cm(3,nat)
c
      double precision delta, delta2i, norm
      double precision ddot
      external ddot
      integer ucont
      integer ish, ilow, ihi, Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer jsh, jlow, jhi, Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      integer ksh, klow, khi, Lk, k_prim, k_gen, k_iexp, k_icfp, k_cent
      integer ibf, jbf, kbf
      integer bs
      integer nshell, nint, count
      integer zcent(3), ixyz, zc
      integer pass
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      call dfill((3*nat*nbf*nbf*nbf),0.0d00,fdm,1)
      call dfill((3*nat*nbf*nbf*nbf),0.0d00,fdp,1)
      call dfill((3*nat*nbf*nbf*nbf),0.0d00,fdo,1)
      call dfill((3*nat*nbf*nbf*nbf),0.0d00,  o,1)
      delta = 0.00001d00
      delta2i = 1.0d00/(2.0d00*delta)
      write(6,*)' nbf = ',nbf
      write(6,*)' nat = ',nat
      write(6,*)' delta,2inverse ',delta,delta2i
      bs = basisin + Basis_Handle_Offset
      nshell = ncont_tot_gb(bs)
      write(6,*)' nshell',nshell
c
      call dcopy (3*nat,coords(1,1,geom),1,c,1)
c
      pass = 0
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basisin,ish,ilow,ihi)) stop 'cn2bfr i'
        ucont   = (sf_ibs_cn2ucn(ish,bs))
        Li      = infbs_cont(CONT_TYPE ,ucont,bs)
        i_prim  = infbs_cont(CONT_NPRIM,ucont,bs)
        i_gen   = infbs_cont(CONT_NGEN ,ucont,bs)
        i_iexp  = infbs_cont(CONT_IEXP ,ucont,bs)
        i_icfp  = infbs_cont(CONT_ICFP ,ucont,bs)
        i_cent  = (sf_ibs_cn2ce(ish,bs))
        do jsh = 1,nshell
          if (.not.bas_cn2bfr(basisin,jsh,jlow,jhi)) stop 'cn2bfr j'
          ucont   = (sf_ibs_cn2ucn(jsh,bs))
          Lj      = infbs_cont(CONT_TYPE ,ucont,bs)
          j_prim  = infbs_cont(CONT_NPRIM,ucont,bs)
          j_gen   = infbs_cont(CONT_NGEN ,ucont,bs)
          j_iexp  = infbs_cont(CONT_IEXP ,ucont,bs)
          j_icfp  = infbs_cont(CONT_ICFP ,ucont,bs)
          j_cent  = (sf_ibs_cn2ce(jsh,bs))
          do ksh = 1,nshell
            if (.not.bas_cn2bfr(basisin,ksh,klow,khi)) stop 'cn2bfr k'
            ucont   = (sf_ibs_cn2ucn(ksh,bs))
            Lk      = infbs_cont(CONT_TYPE ,ucont,bs)
            k_prim  = infbs_cont(CONT_NPRIM,ucont,bs)
            k_gen   = infbs_cont(CONT_NGEN ,ucont,bs)
            k_iexp  = infbs_cont(CONT_IEXP ,ucont,bs)
            k_icfp  = infbs_cont(CONT_ICFP ,ucont,bs)
            k_cent  = (sf_ibs_cn2ce(ksh,bs))
            pass = pass + 1
            write(6,*)' pass ',pass
            if (i_cent.eq.j_cent.and.j_cent.eq.k_cent) goto 00100
c
            nint = int_nbf_x(Li)*int_nbf_x(Lj)*int_nbf_x(Lk)
            write(6,*)' nint = ',nint
c
c* icenter x +
            call dcopy(nat*3,c,1,cp,1)
            cp(1,i_cent) = cp(1,i_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* icenter x -
            call dcopy(nat*3,c,1,cm,1)
            cm(1,i_cent) = cm(1,i_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,1,i_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,1,i_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c* icenter y +
            call dcopy(nat*3,c,1,cp,1)
            cp(2,i_cent) = cp(2,i_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* icenter y -
            call dcopy(nat*3,c,1,cm,1)
            cm(2,i_cent) = cm(2,i_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,2,i_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,2,i_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c            
c* icenter z +
            call dcopy(nat*3,c,1,cp,1)
            cp(3,i_cent) = cp(3,i_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* icenter z -
            call dcopy(nat*3,c,1,cm,1)
            cm(3,i_cent) = cm(3,i_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,3,i_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,3,i_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c
c* jcenter x +
            call dcopy(nat*3,c,1,cp,1)
            cp(1,j_cent) = cp(1,j_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* jcenter x -
            call dcopy(nat*3,c,1,cm,1)
            cm(1,j_cent) = cm(1,j_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,1,j_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,1,j_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c* jcenter y +
            call dcopy(nat*3,c,1,cp,1)
            cp(2,j_cent) = cp(2,j_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* jcenter y -
            call dcopy(nat*3,c,1,cm,1)
            cm(2,j_cent) = cm(2,j_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,2,j_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,2,j_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c            
c* jcenter z +
            call dcopy(nat*3,c,1,cp,1)
            cp(3,j_cent) = cp(3,j_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* jcenter z -
            call dcopy(nat*3,c,1,cm,1)
            cm(3,j_cent) = cm(3,j_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,3,j_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,3,j_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c
c* kcenter x +
            call dcopy(nat*3,c,1,cp,1)
            cp(1,k_cent) = cp(1,k_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* kcenter x -
            call dcopy(nat*3,c,1,cm,1)
            cm(1,k_cent) = cm(1,k_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,1,k_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,1,k_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c* kcenter y +
            call dcopy(nat*3,c,1,cp,1)
            cp(2,k_cent) = cp(2,k_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* kcenter y -
            call dcopy(nat*3,c,1,cm,1)
            cm(2,k_cent) = cm(2,k_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,2,k_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,2,k_cent) = bufm(count)
                enddo
              enddo
            enddo
c            
c            
c* kcenter z +
            call dcopy(nat*3,c,1,cp,1)
            cp(3,k_cent) = cp(3,k_cent)+delta
            call hf3ois(
     &          cp(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cp(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cp(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufp,nint,.false.,.false.,scr,szscr)
c* kcenter z -
            call dcopy(nat*3,c,1,cm,1)
            cm(3,k_cent) = cm(3,k_cent)-delta
            call hf3ois(
     &          cm(1,i_cent), dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, Li,
     &          cm(1,j_cent), dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, Lj,
     &          cm(1,k_cent), dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, Lk,
     &          bufm,nint,.false.,.false.,scr,szscr)
            count = 0
            do ibf = ilow,ihi
              do jbf = jlow,jhi
                do kbf = klow,khi
                  count = count + 1
                  fdp(ibf,jbf,kbf,3,k_cent) = bufp(count)
                  fdm(ibf,jbf,kbf,3,k_cent) = bufm(count)
                enddo
              enddo
            enddo
c
            call hfd3ois(c(1,i_cent),dbl_mb(mb_exndcf(i_iexp,bs)),
     &          dbl_mb(mb_exndcf(i_icfp,bs)),i_prim, 1, Li,
     &          c(1,j_cent),dbl_mb(mb_exndcf(j_iexp,bs)),
     &          dbl_mb(mb_exndcf(j_icfp,bs)),j_prim, 1, Lj,
     &          c(1,k_cent),dbl_mb(mb_exndcf(k_iexp,bs)),
     &          dbl_mb(mb_exndcf(k_icfp,bs)),k_prim, 1, Lk,
     &          dbuf,nint,.false.,scr,szscr)
            call raktest_3cd_2(i_cent,j_cent,k_cent,nint,dbuf,zcent)
c
            count = 0
            do zc = 1,3
              if (zcent(zc).ne.0) then
                do ixyz = 1,3
                  do ibf = ilow,ihi
                    do jbf = jlow,jhi
                      do kbf = klow,khi
                        count = count + 1
                        o(ibf,jbf,kbf,ixyz,zcent(zc)) = dbuf(count)
                      enddo
                    enddo
                  enddo
                enddo
              endif
            enddo
c            
00100       continue
c
          enddo
        enddo
      enddo
c
      call dcopy((nbf*nbf*nbf*3*nat),fdp,1,fdo,1)
      call daxpy((nbf*nbf*nbf*3*nat),-1.0d00,fdm,1,fdo,1)
      call dscal((nbf*nbf*nbf*3*nat),delta2i,fdo,1)
      call printboth_3cd(nbf,nat,fdo,o)
      call dcopy((nbf*nbf*nbf*3*nat),fdo,1,fdm,1)
      call daxpy((nbf*nbf*nbf*3*nat),-1.0d00,o,1,fdm,1)
      norm = ddot((nbf*nbf*nbf*3*nat),fdm,1,fdm,1)
c
      write(6,*)' difference norm: ',norm
c
      end
      subroutine raktest_3cd_2(ic,jc,kc,nint,dbuf,idcent)
      implicit none
      integer ic, jc, kc, nint
      integer idcent(3)
      double precision dbuf(nint,9)
c
      if ((ic.ne.jc).and.(ic.ne.kc).and.(jc.ne.kc)) then
        idcent(1) = ic
        idcent(2) = jc
        idcent(3) = kc
        goto 90000
      endif
      if (ic.eq.jc) then
        call daxpy(3*nint,1.0d00,dbuf(1,4),1,dbuf(1,1),1)
        call dcopy(3*nint,dbuf(1,7),1,dbuf(1,4),1)
        call dfill(3*nint,0.0d00,dbuf(1,7),1)
        write(6,*)' moving'
        idcent(1) = ic
        idcent(2) = kc
        idcent(3) = 0
      else if (ic.eq.kc) then
        call daxpy(3*nint,1.0d00,dbuf(1,7),1,dbuf(1,1),1)
        call dfill(3*nint,0.0d00,dbuf(1,7),1)
        idcent(1) = ic
        idcent(2) = jc
        idcent(3) = 0
      else if (jc.eq.kc) then
        call daxpy(3*nint,1.0d00,dbuf(1,7),1,dbuf(1,4),1)
        call dfill(3*nint,0.0d00,dbuf(1,7),1)
        idcent(1) = ic
        idcent(2) = jc
        idcent(3) = 0
      else
        write(6,*)ic,jc,kc
        stop ' how did I get here'
      endif
90000 continue
      write(6,*)'idcent',idcent
      end
      subroutine printboth_3cd(nbf,nat,fdo,o)
      implicit none
      integer nbf
      integer nat
      double precision fdo(nbf,nbf,nbf,3,nat)
      double precision   o(nbf,nbf,nbf,3,nat)
c
      double precision diff, thresh
      integer i,j,k,ixyz,ic, count
c
      thresh =  1.0d-06
*      thresh =  -1.0d00
c
      count = 0
      do ic = 1,nat
        do ixyz=1,3
          do k=1,nbf
            do j=1,nbf
              do i=1,nbf
                count = count + 1
                diff = fdo(i,j,k,ixyz,ic) - o(i,j,k,ixyz,ic)
                diff = abs(diff)
                if (diff.gt.thresh) then
                  write(6,10000)count,i,j,k,ixyz,ic,
     &                fdo(i,j,k,ixyz,ic),
     &                o(i,j,k,ixyz,ic), diff
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'<',i6,'>','(i=',i3,',j=',
     &    i3,',k=',i3,',x=',i3,',at=',i3,') fd=',
     &    1pd12.5,' o=',1pd12.5,' diff=',1pd12.5)
      end
      subroutine raktest_ovd(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb
c
      integer maxg1, maxs1
      integer h_buf, h_scr, h_cm, h_o, h_fdo, h_fdp, h_fdm
      integer k_buf, k_scr, k_cm, k_o, k_fdo, k_fdp, k_fdm
      integer geom, basis
      integer nat, nbf,  bufsz, szscr, osz
c
      logical int_normalize
      external int_normalize
c
*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 1?'
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_ovd: geom_create failed?',911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_ovd: geom_rtdb_load -ref failed',911, RTDB_ERR)
      if (.not.geom_print(geom)) stop ' print error'
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('raktest_ovd:bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &    call errquit
     &      ('raktest_ovd:bas_rtdb_load failed',911, RTDB_ERR)
      if (.not.bas_print(basis)) stop ' print error'
      if (.not.gbs_map_print(basis)) stop ' print error'
      if (.not.rtdb_print(rtdb,.false.)) stop ' rtdb_p 2?'
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      call intd_init(rtdb,1,basis)
      call int_mem_print()
c
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      if (.not.bas_numbf(basis,nbf)) stop 'bas_numbf fe'
      write(6,*)' total number of basis functions:nbf:',nbf
      call int_mem_1e(maxg1,maxs1)      
      bufsz = 2*maxg1
      szscr = 2*maxs1
      osz   = nbf*nbf*3*nat
      if (.not.ma_push_get(mt_dbl,bufsz,'buf',h_buf,k_buf))
     &    stop ' ma buf failed'
      if (.not.ma_push_get(mt_dbl,szscr,'scr',h_scr,k_scr))
     &    stop ' ma scr failed'
      if (.not.ma_push_get(mt_dbl,(3*nat),'cm',h_cm,k_cm))
     &    stop ' ma cm failed'
      if (.not.ma_push_get(mt_dbl,osz,'o',h_o,k_o))
     &    stop ' ma o failed'
      if (.not.ma_push_get(mt_dbl,osz,'fdo',h_fdo,k_fdo))
     &    stop ' ma fdo failed'
      if (.not.ma_push_get(mt_dbl,osz,'fdp',h_fdp,k_fdp))
     &    stop ' ma fdp failed'
      if (.not.ma_push_get(mt_dbl,osz,'fdm',h_fdm,k_fdm))
     &    stop ' ma fdm failed'
      call raktest_ovd1(dbl_mb(k_buf),bufsz,dbl_mb(k_scr),szscr,
     &    dbl_mb(k_fdo),dbl_mb(k_o),dbl_mb(k_fdp),dbl_mb(k_fdm),
     &    nbf,nat,dbl_mb(k_cm),basis,geom)
      call intd_terminate()
      if (.not.ma_pop_stack(h_fdm)) stop 'h_fdm pop error'
      if (.not.ma_pop_stack(h_fdp)) stop 'h_fdp pop error'
      if (.not.ma_pop_stack(h_fdo)) stop 'h_fdo pop error'
      if (.not.ma_pop_stack(h_o))   stop 'h_o pop error'
      if (.not.ma_pop_stack(h_cm))  stop 'h_cm pop error'
      if (.not.ma_pop_stack(h_scr)) stop 'h_scr pop error'
      if (.not.ma_pop_stack(h_buf)) stop 'h_buf pop error'
      if (.not.bas_destroy(basis)) stop ' bas_destroy failed'
      if (.not.geom_destroy(geom)) stop ' geom_destroy failed'
      end
      subroutine raktest_ovd1(buf,lbuf,scr,lscr,fdo,o,fdp,fdm,
     &    nbf,nat,cmaster,basis,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "inp.fh"
#include "ecp_nwc.fh"
      double precision ddot
      external ddot
c
      integer lbuf, lscr
      double precision buf(lbuf),scr(lscr)
      integer nbf, nat
      double precision fdo(nbf,nbf,3,nat)
      double precision fdp(nbf,nbf,3,nat),fdm(nbf,nbf,3,nat)
      double precision o(nbf,nbf,3,nat) 
      double precision cmaster(3,nat)
      integer basis, geom
*      integer idatom(2)
c
      double precision delta, delta2i, norm
      integer ncont, count, zatom, ratom, ixyz
      integer ish,ilo,ihi,ibf
      integer jsh,jlo,jhi,jbf
      character*256 date_string
      integer lds
c
      if (.not.bas_numcont(basis,ncont)) stop ' bas_numcont'
      write(6,*)' buf size ', lbuf
      write(6,*)' lscr ',lscr
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lscr,0.0d00,scr,1)
      call dfill((nbf*nbf*3*nat),0.0d00,fdo,1)
      call dfill((nbf*nbf*3*nat),0.0d00,fdp,1)
      call dfill((nbf*nbf*3*nat),0.0d00,fdm,1)
      call dfill((nbf*nbf*3*nat),0.0d00,o,1)
      do ish = 1,ncont
        call util_date(date_string)
        lds = inp_strlen(date_string)
        write(6,*)' direct: ish: ',ish,' of ',ncont,
     &      '  ',date_string(1:lds)
        call util_flush(6)
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'
        do jsh = 1,ncont
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr i'
*          call intd_1eov(basis,ish,basis,jsh,lscr,scr,lbuf,buf,idatom)
          call intd_1eh1(basis,ish,basis,jsh,lscr,scr,lbuf,buf)
          count = 0
*          do zatom = 1,2
*            if (idatom(zatom).gt.0) then
*              ratom = idatom(zatom)
          do zatom = 1,nat
            ratom = zatom
            do ixyz = 1,3
              do ibf = ilo,ihi
                do jbf = jlo,jhi
                  count = count + 1
                  o(ibf,jbf,ixyz,ratom) = buf(count)
                enddo
              enddo
            enddo
*            endif    
          enddo
        enddo
      enddo
c
      call dcopy((3*nat),coords(1,1,geom),1,cmaster,1)
      delta = 0.00001d00
      delta2i = 1.0d00/(2.0d00*delta)
      do zatom = 1,nat
        do ixyz = 1,3
        call util_date(date_string)
        lds = inp_strlen(date_string)
        write(6,*)' finite diff: atom: ',zatom,' of ',nat,
     &      ' coord(1:x,2:y,3:z)  ',ixyz,
     &      '  ',date_string(1:lds)
        call util_flush(6)
*--
* +
          call dcopy((3*nat),cmaster,1,coords(1,1,geom),1)
          coords(ixyz,zatom,geom) = coords(ixyz,zatom,geom) + delta
c.. get coordinates for ecp centers.
          if (.not.geom_coords_ecp(geom,dbl_mb(k_xyzecp),n_ecp))
     &        call errquit('ecpdebug: geom_coords_ecp failed',911,
     &       GEOM_ERR)
          do ish = 1,ncont
            call util_date(date_string)
            lds = inp_strlen(date_string)
            write(6,*)' finite diff: ish: ',ish,' of ',ncont,
     &          '  ',date_string(1:lds)
            call util_flush(6)
            if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'
            do jsh = 1,ncont
              if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr i'
*              call int_1eov(basis,ish,basis,jsh,lscr,scr,lbuf,buf)
              call int_1eh1(basis,ish,basis,jsh,lscr,scr,lbuf,buf)
              count = 0
              do ibf = ilo,ihi
                do jbf = jlo,jhi
                  count = count + 1
                  fdp(ibf,jbf,ixyz,zatom) = buf(count)
                enddo
              enddo
            enddo
          enddo
* -
          call dcopy((3*nat),cmaster,1,coords(1,1,geom),1)
          coords(ixyz,zatom,geom) = coords(ixyz,zatom,geom) - delta
c.. get coordinates for ecp centers.
          if (.not.geom_coords_ecp(geom,dbl_mb(k_xyzecp),n_ecp))
     &        call errquit('ecpdebug: geom_coords_ecp failed',911,
     &       GEOM_ERR)
          do ish = 1,ncont
            if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'
            do jsh = 1,ncont
              if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr i'
*              call int_1eov(basis,ish,basis,jsh,lscr,scr,lbuf,buf)
              call int_1eh1(basis,ish,basis,jsh,lscr,scr,lbuf,buf)
              count = 0
              do ibf = ilo,ihi
                do jbf = jlo,jhi
                  count = count + 1
                  fdm(ibf,jbf,ixyz,zatom) = buf(count)
                enddo
              enddo
            enddo
          enddo
*--
        enddo
      enddo
      call dcopy((nbf*nbf*3*nat),fdp,1,fdo,1)
      call daxpy((nbf*nbf*3*nat),-1.0d00,fdm,1,fdo,1)
      call dscal((nbf*nbf*3*nat),delta2i,fdo,1)
      call printboth_ovd(nbf,nat,fdo,o)
      call dcopy((nbf*nbf*3*nat),fdo,1,fdm,1)
      call daxpy((nbf*nbf*3*nat),-1.0d00,o,1,fdm,1)
      norm = ddot((nbf*nbf*3*nat),fdm,1,fdm,1)
c
      write(6,*)' difference norm: ',norm
      end
      subroutine printboth_ovd(nbf,nat,fdo,o)
      implicit none
      integer nbf
      integer nat
      double precision fdo(nbf,nbf,3,nat)
      double precision   o(nbf,nbf,3,nat)
c
      double precision diff, thresh
      integer i,j,ixyz,ic, count, nval_good, n_zero
c
      integer o_res, fd_res
      integer  is_this_val_okay
      external is_this_val_okay
c
      thresh =  1.0d-06
c
      n_zero = 0
      nval_good = 0
      count = 0
      do ic = 1,nat
        do ixyz=1,3
          do j=1,nbf
            do i=1,nbf
              count = count + 1
              o_res  = is_this_val_okay(o(i,j,ixyz,ic)) 
*              write(6,*)' o_res',o_res
              if (o_res .eq. 1) then
                write(6,10001)i,j,ixyz,ic
              elseif (o_res .eq. 2) then
                write(6,10003)i,j,ixyz,ic
              elseif (o_res.eq.0) then
                n_zero = n_zero + 1
              elseif (o_res.eq.3) then
                nval_good = nval_good + 1
              endif
              fd_res = is_this_val_okay(fdo(i,j,ixyz,ic)) 
*              write(6,*)' fd_res',fd_res
              if (fd_res .eq. 1) then
                write(6,10002)i,j,ixyz,ic
              elseif (fd_res .eq. 2) then
                write(6,10004)i,j,ixyz,ic
              endif
              diff = fdo(i,j,ixyz,ic) - o(i,j,ixyz,ic)
              diff = abs(diff)
              if (diff.gt.thresh) then
                write(6,10000)count,i,j,ixyz,ic,
     &              fdo(i,j,ixyz,ic),
     &              o(i,j,ixyz,ic), diff
*              else
*                if (abs(o(i,j,ixyz,ic)).gt.thresh)
*     &              nval_good = nval_good + 1
              endif
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'<',i6,'>','(i=',i3,',j=',
     &    i3,',x=',i3,',at=',i3,') fd=',
     &    1pd12.5,' o=',1pd12.5,' diff=',1pd12.5)
      write(6,*)' zero values',n_zero
      write(6,*)' good values',nval_good
      write(6,*)'         sum',(n_zero+nval_good)
      write(6,*)' count     =',count
10001 format(' calculated  value (',i3,',',i3,',',i3,',',
     &    i3,') is a NAN')
10002 format(' finite diff value (',i3,',',i3,',',i3,',',
     &    i3,') is a NAN')
10003 format(' calculated  value (',i3,',',i3,',',i3,',',
     &    i3,') is an INFINITY')
10004 format(' finite diff value (',i3,',',i3,',',i3,',',
     &    i3,') is an INFINITY')
      end
      subroutine raktest_bd2e(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb
c
      integer basis
      integer geom
      integer nat
      integer nbf, nq_tot
      integer maxg, maxs, bufsz
      integer h_txs, h_nw, h_g, h_l, h_scr
      integer k_txs, k_nw, k_g, k_l, k_scr
      double precision norm
      logical int_normalize
      external int_normalize
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911, RTDB_ERR)
      if (.not.geom_print(geom)) stop ' print error'
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
      if (.not.bas_print(basis)) stop ' print error'
      if (.not.rtdb_print(rtdb,.false.)) stop ' rtdb_p 2?'
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
c
      call intd_init(rtdb,1,basis)
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      if (.not.bas_numbf(basis,nbf)) stop 'bas_numbf fe'
      if (.not.bas_numcont(basis,nq_tot)) stop 'bas_numcont fe'
      call int_mem_2e4c(maxg,maxs)
      maxg = maxg*2
      maxs = maxs*2
      bufsz = nbf*nbf*nbf*nbf*nat*3
c
      if (.not.ma_push_get(mt_dbl,bufsz,'texas',h_txs,k_txs))
     &    stop ' ma texas'
      call dfill(bufsz,0.0d00,dbl_mb(k_txs),1)
      if (.not.ma_push_get(mt_dbl,bufsz,'nwchem',h_nw,k_nw))
     &    stop ' ma nwchem'
      call dfill(bufsz,0.0d00,dbl_mb(k_nw),1)
      if (.not.ma_push_get(mt_dbl,12*maxg,'g',h_g,k_g))
     &    stop ' ma g'
      call dfill(12*maxg,0.0d00,dbl_mb(k_g),1)
      if (.not.ma_push_get(mt_int,4*maxg,'labels',h_l,k_l))
     &    stop ' ma labs'
      call ifill(4*maxg,0,int_mb(k_l),1)
      if (.not.ma_push_get(mt_dbl,maxs,'scratch',h_scr,k_scr))
     &    stop ' ma scratch'
      call dfill(maxs,0.0d00,dbl_mb(k_scr),1)
      write(6,*)' computing with texas '
      call raktest_bd2e_getg(nq_tot,
     &    nbf,nat,rtdb,basis,geom,
     &    dbl_mb(k_txs),
     &    maxg,dbl_mb(k_g),int_mb(k_l),
     &    maxs,dbl_mb(k_scr))
c
      call intd_terminate()
      call int_app_set_no_texas(rtdb)
      call intd_init(rtdb,1,basis)
      write(6,*)' computing with texas off'
      call raktest_bd2e_getg(nq_tot,
     &    nbf,nat,rtdb,basis,geom,
     &    dbl_mb(k_nw),
     &    maxg,dbl_mb(k_g),int_mb(k_l),
     &    maxs,dbl_mb(k_scr))
      call intd_terminate()
      call int_app_unset_no_texas(rtdb)
      call raktest_bd2e_print(basis,nbf,nat,
     &    dbl_mb(k_txs),dbl_mb(k_nw))
      write(6,*)' computing norm'
      call daxpy(((nbf**4)*nat*3),-1.0d00,
     &    dbl_mb(k_txs),1,dbl_mb(k_nw),1)
      norm = ddot(((nbf**4)*nat*3),dbl_mb(k_nw),1,dbl_mb(k_nw),1)
      write(6,*)' difference norm ',norm
c
      if (.not.ma_pop_stack(h_scr)) stop 'ma pop failed for h_scr'
      if (.not.ma_pop_stack(h_l)) stop 'ma pop failed for h_l'
      if (.not.ma_pop_stack(h_g)) stop 'ma pop failed for h_g'
      if (.not.ma_pop_stack(h_nw)) stop 'ma pop failed for h_nw'
      if (.not.ma_pop_stack(h_txs)) stop 'ma pop failed for h_txs'
      end
      subroutine raktest_bd2e_print(basis,nbf,nat,txs,nw)
      implicit none
#include "bas.fh"
      integer nbf,nat
      double precision txs(3,nat,nbf,nbf,nbf,nbf)
      double precision nw(3,nat,nbf,nbf,nbf,nbf)
      integer basis
c
      integer i,j,k,l,xyz,atom
      integer ia,ja,ka,la
      double precision tval 
      double precision nval
      double precision dval
      double precision thresh
      logical printit
      thresh = 1.0d-08
      do i = 1,nbf
        do j = 1,nbf
          do k = 1,nbf
            do l = 1,nbf
              do atom = 1,nat
                do xyz = 1,3
                  tval = txs(xyz,atom,i,j,k,l)
                  nval = nw(xyz,atom,i,j,k,l)
                  dval = abs((tval-nval))
                  printit = dval.gt.thresh
                  if (printit) then
                    write(6,10000)xyz,atom,i,j,k,l,
     &                  tval,nval,dval
                    if (.not.bas_bf2ce(basis,i,ia)) stop 'bf2ce i'
                    if (.not.bas_bf2ce(basis,j,ja)) stop 'bf2ce j'
                    if (.not.bas_bf2ce(basis,k,ka)) stop 'bf2ce k'
                    if (.not.bas_bf2ce(basis,l,la)) stop 'bf2ce l'
                    write(6,10001)ia,ja,ka,la
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'deri [',6i5,'] ',1x,1pd20.10,1x,1pd20.10,1x,1pd20.10)
10001 format(1x,6x,10x,4i5)
      end
      subroutine raktest_bd2e_getg(nq_total,nbf,nat,rtdb,basis,geom,
     &    deri,lbuf,buf,labels,lscr,scr)
      implicit none
#include "errquit.fh"
#include "bas.fh"
      integer nq_total, nbf, nat, rtdb, basis,geom, lbuf, lscr
      double precision deri(3,nat,nbf,nbf,nbf,nbf)
      double precision buf(lbuf), scr(lscr)
      integer labels(lbuf,4)
c
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
c
      integer nqs
      parameter (nqs = 20)
      double precision q4(nqs)
      integer iq(nqs),jq(nqs),kq(nqs),lq(nqs)
      integer qi,qj,qk,ql, q
      double precision dbl_dum
      integer ii,jj,kk,ll,iatom,jatom,katom,latom,xyz,int
      integer count, nintout
      logical more, lastq
c

      q = 0
      do qi = 1,nq_total
        do qj = 1,qi
          do qk = 1,qj
            do ql = 1,qk
              lastq =           qi.eq.nq_total
              lastq = lastq.and.qj.eq.qi
              lastq = lastq.and.qk.eq.qj
              lastq = lastq.and.ql.eq.qk
              if ((q+1).gt.nqs.or.lastq) then
                if (.not.intbd_init4c(
     &              basis,iq,jq,basis,kq,lq,q,q4,.false.,
     &              lscr,scr,lbuf,dbl_dum)) call errquit
     &              ('intbd_init4c failed',911, INT_ERR)
00010           continue
                more = intbd_2e4c(basis,iq,jq,basis,kq,lq,q,q4,.false.,
     &              1.0d-8,.false.,
     &              labels(1,1),labels(1,2),labels(1,3),labels(1,4),
     &              buf,lbuf,nintout,lscr,scr)
                count = 0
                do int = 1,nintout
                  ii = labels(int,1)
                  jj = labels(int,2)
                  kk = labels(int,3)
                  ll = labels(int,4)
                  if (.not.bas_bf2ce(basis,ii,iatom)) stop 'bf2ce i'
                  if (.not.bas_bf2ce(basis,jj,jatom)) stop 'bf2ce j'
                  if (.not.bas_bf2ce(basis,kk,katom)) stop 'bf2ce k'
                  if (.not.bas_bf2ce(basis,ll,latom)) stop 'bf2ce l'
                  do xyz = 1,3
                    count = count + 1
                    deri(xyz,iatom,ii,jj,kk,ll) = buf(count)
                  enddo
                  do xyz = 1,3
                    count = count + 1
                    deri(xyz,jatom,ii,jj,kk,ll) = buf(count)
                  enddo
                  do xyz = 1,3
                    count = count + 1
                    deri(xyz,katom,ii,jj,kk,ll) = buf(count)
                  enddo
                  do xyz = 1,3
                    count = count + 1
                    deri(xyz,latom,ii,jj,kk,ll) = buf(count)
                  enddo
                  if (iatom.eq.jatom) then
                    do xyz = 1,3
                      deri(xyz,iatom,ii,jj,kk,ll) =
     &                    deri(xyz,iatom,ii,jj,kk,ll) +
     &                    deri(xyz,jatom,ii,jj,kk,ll) 
                      deri(xyz,jatom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif
                  if (iatom.eq.katom) then
                    do xyz = 1,3
                      deri(xyz,iatom,ii,jj,kk,ll) =
     &                    deri(xyz,iatom,ii,jj,kk,ll) +
     &                    deri(xyz,katom,ii,jj,kk,ll) 
                      deri(xyz,katom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif
                  if (iatom.eq.latom) then
                    do xyz = 1,3
                      deri(xyz,iatom,ii,jj,kk,ll) =
     &                    deri(xyz,iatom,ii,jj,kk,ll) +
     &                    deri(xyz,latom,ii,jj,kk,ll) 
                      deri(xyz,latom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif
                  if (jatom.eq.katom) then
                    do xyz = 1,3
                      deri(xyz,jatom,ii,jj,kk,ll) =
     &                    deri(xyz,jatom,ii,jj,kk,ll) +
     &                    deri(xyz,katom,ii,jj,kk,ll) 
                      deri(xyz,katom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif
                  if (jatom.eq.latom) then
                    do xyz = 1,3
                      deri(xyz,jatom,ii,jj,kk,ll) =
     &                    deri(xyz,jatom,ii,jj,kk,ll) +
     &                    deri(xyz,latom,ii,jj,kk,ll) 
                      deri(xyz,latom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif
                  if (katom.eq.latom) then
                    do xyz = 1,3
                      deri(xyz,katom,ii,jj,kk,ll) =
     &                    deri(xyz,katom,ii,jj,kk,ll) +
     &                    deri(xyz,latom,ii,jj,kk,ll) 
                      deri(xyz,latom,ii,jj,kk,ll) = 0.0d00
                    enddo
                  endif

                enddo
                if (more) goto 00010
                q = 1
              else
                q = q + 1
              endif
              iq(q) = qi
              jq(q) = qj
              kq(q) = qk
              lq(q) = ql
            enddo
          enddo
        enddo
      enddo
c
      end
      subroutine raktest_diskspeed(rtdb)
      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      integer default_size
      parameter (default_size = 512)
      integer default_count
      parameter (default_count = 10)
      integer default_iterations
      parameter (default_iterations = 10)
      integer size, count, iters
      integer h_io, k_io
      integer h_data, k_data
      integer data_size
c
      if (.not.rtdb_get(rtdb,'rak:disksize',mt_int,1,size))
     &    size = default_size
      if (.not.rtdb_get(rtdb,'rak:diskcount',mt_int,1,count))
     &    count = default_count
      if (.not.rtdb_get(rtdb,'rak:diskiters',mt_int,1,iters))
     &    iters = default_iterations
c
      if (ga_nodeid().eq.0) then
        write(luout,*)
     &      ' size of buffer is            : ',size, ' doubles'
        write(luout,*)
     &      '                              : ',(size*8), ' bytes'
        write(luout,*)
     &      ' number of buffers put out is :',count
        write(luout,*)
     &      ' number of iterations is      :',iters
        write(luout,*)
     &      ' output and reading of        :',
     &      8*size*count*iters,' bytes'
      endif
c
      if (.not.ma_push_get(mt_dbl,size,'io buffer',h_io,k_io))
     &    stop ' ma_get of io buffer failed'
      data_size = 5*ga_nnodes()
      if (.not.ma_push_get(mt_dbl,data_size,'data',h_data,k_data))
     &    stop ' ma get of data buffer failed'
c
      call raktest_diskspeed_fill(size,dbl_mb(k_io))
      call raktest_diskspeed_write(dbl_mb(k_io),size,
     &    count,iters,dbl_mb(k_data),'rakdisk',.true.)
      call raktest_diskspeed_read(dbl_mb(k_io),size,
     &    count,iters,dbl_mb(k_data),'rakdisk',.false.)
c
      if (.not.ma_pop_stack(h_data)) stop ' ma pop error '
      if (.not.ma_pop_stack(h_io)) stop ' ma pop error '
c
      end
      subroutine raktest_diskspeed_read(buf,size,
     &    count,iters,adata,filename_stub,save_file)
      implicit none
#include "stdio.fh"
#include "eaf.fh"
#include "util.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer size
      integer count
      integer iters
      double precision buf(size)
      double precision adata(5,*)
      character*(*)filename_stub
      logical save_file
c
      character*256 filename
      integer i, it
      integer fd
      integer eaf_rv
      double precision offset, bytes_tot
      double precision timec, timew, ratec, ratew
      integer inode
      integer adata_size
      integer mtype
      double precision timew_min, timew_max, timew_ave
      double precision ratew_min, ratew_max, ratew_ave
c
      call ga_sync()
      call util_file_name(filename_stub,.false.,.true.,filename)
c
      eaf_rv = eaf_open(filename,EAF_R,fd)
      if (eaf_rv.ne.0) then
        call eaf_errmsg(eaf_rv)
        stop 'eaf open error'
      endif
c
      timec = util_cpusec()
      timew = util_wallsec()
      do it = 1,iters
        offset = 0.0d00
        do i = 1,count
          eaf_rv = eaf_read(fd,offset,buf,8*size)
          if (eaf_rv.ne.0) then
            call eaf_errmsg(eaf_rv)
            stop ' eaf read error'
          endif
          offset = offset + dble(8*size)
        enddo
      enddo
      timec = util_cpusec() - timec
      timew = util_wallsec() - timew
      if (ga_nodeid().eq.0) then
        write(luout,*)' EAF info for node 0 only '
        call eaf_print_stats(fd)
      endif
      call ga_sync()
c
      eaf_rv = eaf_close(fd)
      if (eaf_rv.ne.0) then
        call eaf_errmsg(eaf_rv)
        stop 'eaf close error'
      endif
c
      if (.not.save_file) call util_file_unlink(filename)
c      
      bytes_tot = dble(count)*dble(8)*dble(size)
      bytes_tot = bytes_tot*dble(iters)
      bytes_tot = bytes_tot*1.0d-6
      ratec = bytes_tot/timec
      ratew = bytes_tot/timew
      adata_size = 5*ga_nnodes()
      call dfill(adata_size,0.0d00,adata,1)
      inode = ga_nodeid() + 1
      adata(1,inode) = bytes_tot
      adata(2,inode) = timec
      adata(3,inode) = ratec
      adata(4,inode) = timew
      adata(5,inode) = ratew
      mtype = 1234 + MSGDBL
      call ga_dgop(mtype,adata,adata_size,'+')
      if (ga_nodeid().eq.0) then
        do inode = 1,ga_nnodes()
          write(luout,*)
     &        ' read  statistics for node :',(inode-1)
          write(luout,10000)
     &        ' read  total Mbytes        :',adata(1,inode)
          write(luout,10000)
     &        ' read  cpu  time (s)       :',adata(2,inode)
          write(luout,10000)
     &        ' read  cpu  rate (mb/s)    :',adata(3,inode)
          write(luout,10000)
     &        ' read  wall time (s)       :',adata(4,inode)
          write(luout,10000)
     &        ' read  wall rate (mb/s)    :',adata(5,inode)
          write(luout,*)' '
          call util_flush(luout)
        enddo
      endif
      call ga_sync()
      call ga_sync()
      timew_min = timew
      call ga_dgop(5651,timew_min,1,'min')
      timew_max = timew
      call ga_dgop(5652,timew_max,1,'max')
      timew_ave = timew
      call ga_dgop(5653,timew_ave,1,'+')
      timew_ave = timew_ave/dble(ga_nnodes())
      ratew_min = ratew
      call ga_dgop(5654,ratew_min,1,'min')
      ratew_max = ratew
      call ga_dgop(5655,ratew_max,1,'max')
      ratew_ave = ratew
      call ga_dgop(5656,ratew_ave,1,'+')
      ratew_ave = ratew_ave/dble(ga_nnodes())
      call ga_sync()
      call ga_sync()
      call ga_sync()
      if (ga_nodeid().eq.0) then
        write(luout,10000)' read  minimum wall time   :',timew_min
        write(luout,10000)' read  maximum wall time   :',timew_max
        write(luout,10000)' read  average wall time   :',timew_ave
        write(luout,10000)' read  minimum wall rate   :',ratew_min
        write(luout,10000)' read  maximum wall rate   :',ratew_max
        write(luout,10000)' read  average wall rate   :',ratew_ave
      endif
c
10000 format(1x,a,f10.5)
      end
      subroutine raktest_diskspeed_write(buf,size,
     &    count,iters,adata,filename_stub,save_file)
      implicit none
#include "stdio.fh"
#include "eaf.fh"
#include "util.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer size
      integer count
      integer iters
      double precision buf(size)
      double precision adata(5,*)
      character*(*)filename_stub
      logical save_file
c
      character*256 filename
      integer i, it
      integer fd
      integer eaf_rv
      double precision offset, bytes_tot
      double precision timec, timew, ratec, ratew
      integer inode
      integer adata_size
      integer mtype
      double precision timew_min, timew_max, timew_ave
      double precision ratew_min, ratew_max, ratew_ave
c
      call ga_sync()
      call util_file_name(filename_stub,.false.,.true.,filename)
c
      eaf_rv = eaf_open(filename,EAF_RW,fd)
      if (eaf_rv.ne.0) then
        call eaf_errmsg(eaf_rv)
        stop 'eaf open error'
      endif
c
      timec = util_cpusec()
      timew = util_wallsec()
      do it = 1,iters
        offset = 0.0d00
        do i = 1,count
          eaf_rv = eaf_write(fd,offset,buf,8*size)
          if (eaf_rv.ne.0) then
            call eaf_errmsg(eaf_rv)
            stop ' eaf write error'
          endif
          offset = offset + dble(8*size)
        enddo
      enddo
      timec = util_cpusec() - timec
      timew = util_wallsec() - timew
      if (ga_nodeid().eq.0) then
        write(luout,*)' EAF info for node 0 only '
        call eaf_print_stats(fd)
      endif
      call ga_sync()
c
      eaf_rv = eaf_close(fd)
      if (eaf_rv.ne.0) then
        call eaf_errmsg(eaf_rv)
        stop 'eaf close error'
      endif
c
      if (.not.save_file) call util_file_unlink(filename)
c      
      bytes_tot = dble(count)*dble(8)*dble(size)
      bytes_tot = bytes_tot*dble(iters)
      bytes_tot = bytes_tot*1.0d-6
      ratec = bytes_tot/timec
      ratew = bytes_tot/timew
      adata_size = 5*ga_nnodes()
      call dfill(adata_size,0.0d00,adata,1)
      inode = ga_nodeid() + 1
      adata(1,inode) = bytes_tot
      adata(2,inode) = timec
      adata(3,inode) = ratec
      adata(4,inode) = timew
      adata(5,inode) = ratew
      mtype = 2134 + MSGDBL
      call ga_dgop(mtype,adata,adata_size,'+')
      if (ga_nodeid().eq.0) then
        do inode = 1,ga_nnodes()
          write(luout,*)
     &        ' write statistics for node :',(inode-1)
          write(luout,10000)
     &        ' write total Mbytes        :',adata(1,inode)
          write(luout,10000)
     &        ' write cpu  time (s)       :',adata(2,inode)
          write(luout,10000)
     &        ' write cpu  rate (mb/s)    :',adata(3,inode)
          write(luout,10000)
     &        ' write wall time (s)       :',adata(4,inode)
          write(luout,10000)
     &        ' write wall rate (mb/s)    :',adata(5,inode)
          write(luout,*)' '
          call util_flush(luout)
        enddo
      endif
      call ga_sync()
      call ga_sync()
      timew_min = timew
      call ga_dgop(5657,timew_min,1,'min')
      timew_max = timew
      call ga_dgop(5658,timew_max,1,'max')
      timew_ave = timew
      call ga_dgop(5659,timew_ave,1,'+')
      timew_ave = timew_ave/dble(ga_nnodes())
      ratew_min = ratew
      call ga_dgop(5660,ratew_min,1,'min')
      ratew_max = ratew
      call ga_dgop(5661,ratew_max,1,'max')
      ratew_ave = ratew
      call ga_dgop(5662,ratew_ave,1,'+')
      ratew_ave = ratew_ave/dble(ga_nnodes())
      call ga_sync()
      call ga_sync()
      call ga_sync()
      if (ga_nodeid().eq.0) then
        write(luout,10000)' write minimum wall time   :',timew_min
        write(luout,10000)' write maximum wall time   :',timew_max
        write(luout,10000)' write average wall time   :',timew_ave
        write(luout,10000)' write minimum wall rate   :',ratew_min
        write(luout,10000)' write maximum wall rate   :',ratew_max
        write(luout,10000)' write average wall rate   :',ratew_ave
        call util_flush(luout)
      endif
c
10000 format(1x,a,f10.5)
      end
      subroutine raktest_diskspeed_fill(size,buf)
      implicit none
      integer size
      double precision buf(size)
c
      integer i
c
      do i = 1,size
        buf(size) = sqrt(dble(i))
      enddo
      end
      subroutine raktest_2ecompare(rtdb)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
c:functions
      logical int_normalize
      external int_normalize
c:passed
      integer rtdb
c:local
      integer geom, basis
      integer szbuf, Lmax, nshell
      integer maxg, maxscr
      integer size_dbls, h_dbls, k_dbls
      integer size_ints, h_ints, k_ints
      integer k_bufnw, k_buftxs, k_scr, k_lab, k_eri
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_geomwrt: geom_create failed?',911,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911, RTDB_ERR)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
      if (ga_nodeid().eq.0) then
        if (.not.geom_print(geom)) stop ' print error'
        if (.not.bas_print(basis)) stop ' print error'
        if (.not.gbs_map_print(basis)) stop ' gbs_map_print 2?'
      endif
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont fe'
      if (.not.bas_high_angular(basis,Lmax)) stop 'bas_ha fe'
c
      szbuf = (Lmax+1)*(Lmax+2)/2
      szbuf = szbuf**4
c
      call rak_init(rtdb,1,basis)
      call intb_mem_2e4c(maxg,maxscr)
      maxscr = maxscr + szbuf
      maxscr = maxscr + maxscr/5 + 1
      maxscr = max(51000,maxscr)
      call int_terminate()
c
      if (szbuf.lt.maxg) then
        write(luout,*)' szbuf =',szbuf
        write(luout,*)' maxg  =',maxg
        call errquit('raktest_2ecompare:fatal error',911, UNKNOWN_ERR)
      endif
c
      size_dbls = 3*szbuf+maxscr
      size_ints = 4*szbuf
      if (ga_nodeid().eq.0) then
        write(luout,*)' raktest: maxscr   :',maxscr
        write(luout,*)' raktest: szbuf    :',szbuf
        write(luout,*)' raktest: size_dbls:',size_dbls
        write(luout,*)' raktest: size_ints:',size_ints
        call util_flush(luout)
      endif
      if (.not.ma_push_get(mt_dbl,size_dbls,'dbls scr',h_dbls,k_dbls))
     &    call errquit('raktest_2ecompare:fatal error:ma dbls',911,
     &       MA_ERR)
      if (.not.ma_push_get(mt_int,size_ints,'ints scr',h_ints,k_ints))
     &    call errquit('raktest_2ecompare:fatal error:ma ints',911,
     &       MA_ERR)
      k_bufnw  = k_dbls
      k_buftxs = k_bufnw + szbuf
      k_eri    = k_buftxs + szbuf
      k_scr    = k_eri + szbuf
      k_lab    = k_ints
c
      call ga_sync()
      
      call raktest_2ecompare_a(
     &    dbl_mb(k_bufnw),
     &    dbl_mb(k_buftxs),
     &    dbl_mb(k_scr),
     &    dbl_mb(k_eri),
     &    int_mb(k_lab),
     &    szbuf,maxscr,nshell,basis,geom,rtdb)
c
      if (.not.ma_pop_stack(h_ints)) stop ' ma pop error'
      if (.not.ma_pop_stack(h_dbls)) stop ' ma pop error'
      if (.not.bas_destroy(basis)) stop 'bas_destroy fe'
      if (.not.geom_destroy(geom)) stop 'geom destroy fe'
c
      end
      subroutine raktest_2ecompare_a(bufnw, buftxs, scr, eri,
     &    labs, szbuf, maxscr, nsh, basis, geom, rtdb)
      implicit none
#include "bas.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
c:functions
      logical intb_init4c, intb_2e4c
      external intb_init4c, intb_2e4c
c:passed
      integer rtdb
      integer szbuf
      integer maxscr
      integer nsh
      integer basis
      integer geom
      double precision eri(szbuf)
      double precision bufnw(szbuf)
      double precision buftxs(szbuf)
      integer labs(szbuf,4)
      double precision scr(maxscr)
c:local
      integer me, nproc
      integer ish, jsh, ksh, lsh
      integer tish, tjsh, tksh, tlsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision q4, dummy, norm
      integer nint_t, nint_n
      logical more
      double precision zero,thresh
      double precision nwmax, txsmax
      double precision prev_time, now_time, delta_time
      integer count, wrong
      integer task_count
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      zero = 1.0d-09
      thresh = zero*0.01d00
c
*      call setdbg(1)
      call int_app_set_no_spint(rtdb)
      call rak_init(rtdb,1,basis)
c
      wrong = 0
      count = 0
      task_count = (me-1)
      do ish = nsh,1,-1
        if (ish.eq.nsh) then
          delta_time = 0.0d0
          prev_time = util_wallsec()
        else
          now_time = util_wallsec()
          delta_time = now_time - prev_time
          prev_time = now_time
        endif
        call ga_sync()
        if (me.eq.0) then
          write(luout,10001)ish,delta_time
          call util_flush(luout)
        endif
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'        
        do jsh = ish,1,-1
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr j'        
          do ksh = jsh,1,-1
            if (.not.bas_cn2bfr(basis,ksh,klo,khi)) stop 'cn2bfr k'        
            do lsh = ksh,1,-1
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi)) stop 'cn2bfr l'        
              count = count + 1
              task_count = task_count + 1
              if (mod(task_count,nproc).eq.0) then
*                write(luout,*)ish,jsh,ksh,lsh,':',
*     &              task_count,count,me
*                call util_flush(luout)
c... do default
                tish = ish
                tjsh = jsh
                tksh = ksh
                tlsh = lsh
                call ifill(4*szbuf,0,labs,1)
                call dfill(szbuf,0.0d00,buftxs,1)
*                call rak_init(rtdb,1,basis)
                if (.not.intb_init4c(
     &              basis,tish,tjsh,basis,tksh,tlsh,1,
     &              q4,.false.,maxscr,scr,szbuf,dummy))
     &              stop 'intb_init failed 1'
00100           continue
                call dfill (szbuf,0.0d00,eri,1)
*                write(6,*)ish,jsh,ksh,lsh,':',' intb_2e4c 1',me
*                call util_flush(luout)
                more = intb_2e4c(
     &              basis,tish,tjsh,basis,tksh,tlsh,1,
     &              q4,.false.,zero,.false.,
     &              labs(1,1),labs(1,2),labs(1,3),labs(1,4),
     &              eri,szbuf,nint_t,maxscr,scr)
                call raktest_2ecompare_cp(buftxs,eri,nint_t,
     &              labs(1,1),labs(1,2),labs(1,3),labs(1,4),
     &              ilo,ihi,jlo,jhi,klo,khi,llo,lhi,nwmax)
                if (more) then
*                  write(luout,*)ish,jsh,ksh,lsh,':',' more 2 ?',me
*                  call util_flush(luout)
                  goto 00100
                endif
*                call int_terminate()
c... do no texas .eg., force nwchem integrals
                tish = ish
                tjsh = jsh
                tksh = ksh
                tlsh = lsh
                call ifill(4*szbuf,0,labs,1)
                call dfill(szbuf,0.0d00,bufnw,1)
*                call int_app_set_no_texas(rtdb)
*                call rak_init(rtdb,1,basis)
                call rak_notxs()
                if (.not.intb_init4c(
     &              basis,tish,tjsh,basis,tksh,tlsh,1,
     &              q4,.false.,maxscr,scr,szbuf,dummy))
     &              stop 'intb_init failed 2'
00200           continue
                call dfill (szbuf,0.0d00,eri,1)
*                write(6,*)ish,jsh,ksh,lsh,':',' intb_2e4c 2',me
*                call util_flush(luout)
                more = intb_2e4c(
     &              basis,tish,tjsh,basis,tksh,tlsh,1,
     &              q4,.false.,zero,.false.,
     &              labs(1,1),labs(1,2),labs(1,3),labs(1,4),
     &              eri,szbuf,nint_n,maxscr,scr)
                call raktest_2ecompare_cp(bufnw,eri,nint_n,
     &              labs(1,1),labs(1,2),labs(1,3),labs(1,4),
     &              ilo,ihi,jlo,jhi,klo,khi,llo,lhi,txsmax)
                if (more) then
*                  write(luout,*)ish,jsh,ksh,lsh,':',' more 2 ?',me
*                  call util_flush(luout)
                  goto 00200
                endif
                call rak_nonotxs()
*                call int_terminate()
*                call int_app_unset_no_texas(rtdb)
c
*                write(6,*)ish,jsh,ksh,lsh,':',' dcopy',me
*                call util_flush(luout)
                call dcopy(szbuf,bufnw,1,eri,1)
                call daxpy(szbuf,-1.0d00,buftxs,1,eri,1)
                norm = ddot(szbuf,eri,1,eri,1)
                if (norm.gt.thresh)then
                  wrong = wrong + 1
                  call raktest_2ecompare_mdiff(szbuf,bufnw,buftxs,nwmax)
                  write(luout,10000)
     &                me,wrong,count,ish,jsh,ksh,lsh,
     &                nint_t,nint_n,norm,nwmax
                  call util_flush(luout)
                endif
                if (nint_t.ne.nint_n.and.norm.gt.thresh)then
                  wrong = wrong + 1
                  call raktest_2ecompare_mdiff(szbuf,bufnw,buftxs,nwmax)
                  if (nwmax.gt.zero) then
                    write(luout,10000)
     &                  me,wrong,count,ish,jsh,ksh,lsh,
     &                  nint_t,nint_n,norm,
     &                  nwmax
                    call util_flush(luout)
                  endif
                endif
c
              endif
            enddo
          enddo
        enddo
      enddo
c
*      call ga_sync()
c
      delta_time = util_wallsec() - prev_time
      if (me.eq.0) then
        write(luout,10001) 0, delta_time
        call util_flush(luout)
      endif
      call ga_dgop(5662,wrong,1,'+')
      if (me.eq.0) write(luout,10002) wrong,count
c
      call int_terminate()
      call int_app_unset_no_spint(rtdb)
c
10000 format(1x,i3,':',i5,' of',i8,'(',i4,i4,'|',i4,i4,') [t:',i5,
     &    '|n:',i5,']  norm=',2(1pd15.6))
10001 format(1x,60('-'),'>',i5,1x,f8.2 )
10002 format(1x,i5,' of',i8,' quartets are bad')
c
      end
      subroutine rak_notxs()
      implicit none
* 
* user_* variables determine .true. a user set some value for the 
*        specific integral code.  
* def_* variables is the value that the user set.
*
* this means that if the user does not want to run the sp integral code
* he/she would set "int:cando_sp" false and the values of would be
* user_cando_sp = .true. and def_cando_sp = .false.
*
* to test then use: 
*
* if(user_cando_sp.and.(.not.def_cando_sp) then
*    do not do anything with sp code
* endif
* 
* or 
* 
* if (.not.((user_cando_sp.and.(.not.def_cando_sp)))) call sp_code
*
*
* Ricky A. Kendall, HPCCG, EMSL, PNNL 
*
      logical user_cando_sp  ! did user set a value for sp 
      logical user_cando_nw  ! did user set a value for nw 
      logical user_cando_txs ! did user set a value for txs
      logical def_cando_sp   ! default user setable value for cando_sp
      logical def_cando_nw   ! default user setable value for cando_nw
      logical def_cando_txs  ! default user setable value for cando_txs
c
      logical app_stored_txs   ! value stored in int_app_set_no_texas
      logical app_stored_spint ! value stored in int_app_set_no_spint
      integer rtdbIused
c
      common /clcando/ user_cando_sp, user_cando_nw, user_cando_txs,
     &    def_cando_sp, def_cando_nw, def_cando_txs,
     &    app_stored_txs, app_stored_spint,
     &    rtdbIused
c
      user_cando_txs = .true.
      def_cando_txs = .false.
c
      end
      subroutine rak_nonotxs()
      implicit none
* 
* user_* variables determine .true. a user set some value for the 
*        specific integral code.  
* def_* variables is the value that the user set.
*
* this means that if the user does not want to run the sp integral code
* he/she would set "int:cando_sp" false and the values of would be
* user_cando_sp = .true. and def_cando_sp = .false.
*
* to test then use: 
*
* if(user_cando_sp.and.(.not.def_cando_sp) then
*    do not do anything with sp code
* endif
* 
* or 
* 
* if (.not.((user_cando_sp.and.(.not.def_cando_sp)))) call sp_code
*
*
* Ricky A. Kendall, HPCCG, EMSL, PNNL 
*
      logical user_cando_sp  ! did user set a value for sp 
      logical user_cando_nw  ! did user set a value for nw 
      logical user_cando_txs ! did user set a value for txs
      logical def_cando_sp   ! default user setable value for cando_sp
      logical def_cando_nw   ! default user setable value for cando_nw
      logical def_cando_txs  ! default user setable value for cando_txs
c
      logical app_stored_txs   ! value stored in int_app_set_no_texas
      logical app_stored_spint ! value stored in int_app_set_no_spint
      integer rtdbIused
c
      common /clcando/ user_cando_sp, user_cando_nw, user_cando_txs,
     &    def_cando_sp, def_cando_nw, def_cando_txs,
     &    app_stored_txs, app_stored_spint,
     &    rtdbIused
c
      user_cando_txs = .false.
      def_cando_txs = .false.
c
      end
      subroutine raktest_2ecompare_mdiff(szbuf,bufnw,buftxs,valmax)
      implicit none
      integer szbuf
      double precision bufnw(szbuf),buftxs(szbuf),valmax
c
      double precision diff
      integer i
      valmax = 0.0d00
      do i = 1,szbuf
        diff = bufnw(i)-buftxs(i)
        valmax = max(abs(diff),valmax)
      enddo
      end
      subroutine raktest_2ecompare_cp(buf,eri,nint,
     &            ilb,jlb,klb,llb,
     &            ilo,ihi,jlo,jhi,klo,khi,llo,lhi,valmax)
      implicit none
      integer nint
      integer ilb(nint),jlb(nint),klb(nint),llb(nint)
      integer ilo,ihi,jlo,jhi,klo,khi,llo,lhi
      double precision eri(nint)
      double precision buf(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      double precision valmax
c
      integer n,i,j,k,l
      valmax = 0.0d00
      do n = 1,nint
        i = ilb(n)
        j = jlb(n)
        k = klb(n)
        l = llb(n)
        buf(l,k,j,i) = eri(n)
        valmax = max(valmax,abs(eri(n)))
      enddo
      end
      subroutine rak_init(rtdb, nbas, bases)
c
c initializes integral code to data structers for a integral computation
c: no print
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "apiP.fh"
* 
* user_* variables determine .true. a user set some value for the 
*        specific integral code.  
* def_* variables is the value that the user set.
*
* this means that if the user does not want to run the sp integral code
* he/she would set "int:cando_sp" false and the values of would be
* user_cando_sp = .true. and def_cando_sp = .false.
*
* to test then use: 
*
* if(user_cando_sp.and.(.not.def_cando_sp) then
*    do not do anything with sp code
* endif
* 
* or 
* 
* if (.not.((user_cando_sp.and.(.not.def_cando_sp)))) call sp_code
*
*
* Ricky A. Kendall, HPCCG, EMSL, PNNL 
*
      logical user_cando_sp  ! did user set a value for sp 
      logical user_cando_nw  ! did user set a value for nw 
      logical user_cando_txs ! did user set a value for txs
      logical def_cando_sp   ! default user setable value for cando_sp
      logical def_cando_nw   ! default user setable value for cando_nw
      logical def_cando_txs  ! default user setable value for cando_txs
c
      logical app_stored_txs   ! value stored in int_app_set_no_texas
      logical app_stored_spint ! value stored in int_app_set_no_spint
      integer rtdbIused
c
      common /clcando/ user_cando_sp, user_cando_nw, user_cando_txs,
     &    def_cando_sp, def_cando_nw, def_cando_txs,
     &    app_stored_txs, app_stored_spint,
     &    rtdbIused
      
c

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
c::functions
      logical spcart_init
      external spcart_init
      logical int_ecp_init
      external int_ecp_init
c::passed
c:tex-\begin{verbatim}
      integer rtdb        ! [input] run time data base handle
      integer nbas        ! [input] number of basis sets to be used
      integer bases(nbas) ! [input] basis set handles
c:tex-\end{verbatim}
c::local
      integer txs_mem_min ! memory from texas
      integer ibas, ang2use, angm, type
      logical status
      integer nqmax_texas  ! maximum number of quartets in texas blocking interface
      parameter (nqmax_texas = 10000)
c
c      block data api_data
c
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
c     
      call int_mem_zero()
c
      DCexp     = 0.0D00
      DCcoeff   = 1.0D00
      val_int_acc = 0.0d00
c
      intd_memthresh = 0
      numd_tot       = 0
      numd_okay      = 0
      numd_red       = 0
c
      if(init_int.eq.1) then
        write(6,*)' warning nested int_inits'
        write(6,*)' int_init already called '
        call util_flush(6)
      endif
c
c initialize type-> nbf maps
c
      int_nbf_x(-1) = 4
      int_nbf_s(-1) = 4
      do type = 0,int_nbf_max_ang
        int_nbf_x(type) = (type+1)*(type+2)/2
        int_nbf_s(type) = 2*type+1
      enddo
     
c
c initialize cando information from rtdb
c
      user_cando_sp   = .false.
      user_cando_nw   = .false.
      user_cando_txs  = .false.
      def_cando_sp    = .false.
      def_cando_nw    = .false.
      def_cando_txs   = .false.
c
      if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
        user_cando_sp = .true.
        def_cando_sp  = status
*        if (ga_nodeid().eq.0) then
*          write(6,*)
*     &        ' int_init: cando_sp set to always be ',def_cando_sp
*          call util_flush(6)
*        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
        user_cando_nw = .true.
        def_cando_nw  = status
*rak:        if (ga_nodeid().eq.0) then
*rak:          write(6,*)
*rak:     &        ' int_init: cando_nw set to always be ',def_cando_nw
*rak:          call util_flush(6)
*rak:        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
        user_cando_txs = .true.
        def_cando_txs  = status
*rak:        if (ga_nodeid().eq.0) then
*rak:          write(6,*)
*rak:     &        ' int_init: cando_txs set to always be ',def_cando_txs
*rak:          call util_flush(6)
*rak:        endif
      endif
* sanity checking: e.g., you only want to turn off a particular integral
* code never always turn it on.
*
      if (def_cando_sp.or.def_cando_nw.or.def_cando_txs) then
        if (ga_nodeid().eq.0) then
          write(6,*)' you are trying to turn an integral code on ? '
          write(6,*)' sp  ', def_cando_sp
          write(6,*)' nw  ', def_cando_nw
          write(6,*)' txs ', def_cando_txs
          call util_flush(6)
        endif
        call errquit
     &      ('int_init: logic error with user cando settings',911,
     &       INT_ERR)
      endif
c
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'int_init')
00100 continue

      if (.not.status) then
        write(6,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(6,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('int_init: basis handles hosed ',nbas, INT_ERR)
      endif
*      write(6,*)' int_init: basis set handles valid '
c
c check for both sp and gc shells
c
      call int_bothsp_gc_check(bases,nbas,'int_init')
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('int_init: angm error',angm, INPUT_ERR)
        ang2use = max(ang2use,angm)
00300 continue
*
* test for higher than h functions  0123456
      if (ang2use.ge.6) call errquit
     &    ('only basis sets with s through h functions are allowed',
     &    911, BASIS_ERR)
*
      call defNxyz(ang2use)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init(ang2use,.true.,.false.))) then
        call errquit('int_init: spcart_init failed',911, INT_ERR)
      endif
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(rtdb,bases,nbas)
      call sp_init(nbas,bases)
      call init70               ! To generate tables etc.
      call int_acc_std()
* def u=f d=f -> f.and.!f -> f -> e = t
* no txs u=t d=f -> t.and.!f -> t -> e = f
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_init(rtdb,nbas,bases,nqmax_texas,txs_mem_min,
     *                  'scfd_int')
      endif
c
c See if any basis has an attached ECP
c
      any_ecp = .false.
      ecp_bsh = 0
      do ibas = 1,nbas
        if (bas_get_ecp_handle(bases(ibas),ecp_bsh)) then
          any_ecp = .true.
          goto 00001
        endif
      enddo
00001 continue
      if (any_ecp) then
        if (.not.ecp_check_handle(ecp_bsh,'int_init')) call errquit
     &        ('int_init: ecp handle is invalid fatal error',911,
     &       INT_ERR)
        if (.not.int_ecp_init(ecp_bsh,0,0)) call errquit
     &        ('int_init: int_ecp_init failed ',911, INT_ERR)
      endif
      init_int = 1
      end

