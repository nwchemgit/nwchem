      subroutine raktest(rtdb)
      implicit none
c $Id: raktest.F,v 1.56 1996-07-26 16:25:31 d3g681 Exp $
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb          ! rtdb handle
c::local
      integer raktask, rak_tmp
c
      raktask = 0
      if (rtdb_get(rtdb,'raktask',MT_INT,1,rak_tmp))
     &    raktask = rak_tmp
c
      goto (00010,00020,00030,00040,00050,00060,00070,00080,
     &    00090, 00100) raktask
c
      write(6,*)' default raktest task '
      write(6,*)' test semi empirical interface '
      call raktest_semi(rtdb)
      goto 00001
c
00010 write(6,*)' raktest task 1 stepper test'
      call raktest_stpr(rtdb)
      goto 00001
00020 write(6,*)' raktest task 2 check int_init'
      call raktest_init(rtdb)
      goto 00001
00030 write(6,*)' raktest task 3 check intd_init'
      call raktest_initd(rtdb)
      goto 00001
00040 write(6,*)' raktest check 3ctr nai'
      call raktest_3ctr(rtdb)
      goto 00001
00050 write(6,*)' test of general contraction code '
      call raktest_gc(rtdb)
      goto 00001
00060 write(6,*)' test of orbital printing code '
      call raktest_printorb(rtdb)
      goto 00001
00070 write(6,*)' test of writing geom objects out '
      call raktest_geomwrt(rtdb)
      goto 00001
00080 write(6,*)' test of spcart stuff '
      call raktest_spcart(rtdb)
      goto 00001
00090 write(6,*)' test of spcart stuff all in one'
      call raktest_test9(rtdb)
      goto 00001
00100 write(6,*)' test of ecp stuff '
      call raktest_ecp(rtdb)
      goto 00001
c
00001 continue
      end
      subroutine raktest_geomwrt(rtdb)
      implicit none
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
     &    call errquit('raktest_geomwrt: geom_create failed?',911)
c
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911)
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
     &      ('stpr_walk: geom_rtdb_store (of copy) failed',911)
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
          call errquit('int_xyz_tag: error on lval=-1,ic=',ic)          
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
        call errquit('int_xyz_tag: error on lval=',lval)
      endif
      do i = ixyz, l_tag
        xyz_tag(i:i) = ' '
      enddo
      end
c...............................................................................
      subroutine raktest_spcart(rtdb)
      implicit none 
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
     &      ('geom create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &      ('geom_rtdb_load failed',911)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &      ('bas_rtdb_load failed',911)
c
      write(6,*)' geom/basis loaded'
c
      write(6,*)' raw basis '
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0)
      write(6,*)' first normalize'
      if (.not.int_normalize(basis)) stop ' norm error 1'
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0)
      write(6,*)' second normalize'
      if (.not.int_normalize(basis)) stop ' norm error 2'
      if (.not. bas_print(basis))
     $      call errquit(' basis print failed', 0)
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &      ('numbf failed',911)
c
      nbfsq = nbf*nbf
      if (.not.ma_push_get(mt_dbl,nbfsq,'square cart overlap',
     &      h_cart_s, k_cart_s)) call errquit
     &      (' cart overlap ma failed ',911)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square spher overlap',
     &      h_sph_s, k_sph_s)) call errquit
     &      (' sph overlap ma failed ',911)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square cart overlap 2',
     &      h_cart_s2, k_cart_s2)) call errquit
     &      (' cart2 overlap ma failed ',911)
      if (.not.ma_push_get(mt_dbl,nbfsq,'square spher overlap 2',
     &      h_sph_s2, k_sph_s2)) call errquit
     &      (' sph2 overlap ma failed ',911)
c
      if (.not.bas_numcont(basis,nshell)) call errquit
     &      ('numcont error',911)
c
      call int_init(rtdb,1,basis)
      call int_mem(max1e,max2e,mscr1,mscr2)
      m_buf = max(max1e*2,max2e*2)
      m_scr = max(mscr1*2,mscr2)
      m_buf = m_buf + (m_buf*110)/100
      m_scr = m_scr + (m_scr*110)/100
      
c
      if (.not.ma_push_get(mt_dbl,m_scr,'scr for 1e',h_scr,k_scr))
     &      call errquit('ma scr failed',911)
c
      if (.not.ma_push_get(mt_dbl,m_buf,'buf for 1e',h_buf,k_buf))
     &      call errquit('ma buf failed',911)
c
      if (.not.ma_push_get(mt_int,2*nshell,'buf for sp cn2bfr',
     &      h_2bfr,k_2bfr))
     &      call errquit('ma buf failed',911)
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
     &      h_eri_1,k_eri_1)) call errquit('ma failed',911)
      if (.not.ma_push_get(mt_dbl,(nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &    'eri_sp 1',
     &      h_eri_sp_1,k_eri_sp_1)) call errquit('ma failed',911)
      if (.not.ma_push_get(mt_dbl,(nbf*nbf*nbf*nbf),
     &      'eri cart 2',
     &      h_eri_2,k_eri_2)) call errquit('ma failed',922)
      if (.not.ma_push_get(mt_dbl,(nbf_sp*nbf_sp*nbf_sp*nbf_sp),
     &      'eri_sp  2',
     &      h_eri_sp_2,k_eri_sp_2)) call errquit('ma failed',922)
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
     &      ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,sp_basis,'ao sp_basis'))
     &      call errquit
     &      ('bas_rtdb_load failed',911)
c
      write(6,*)' geom/sp_basis loaded'
c
      write(6,*)' raw sp_basis '
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0)
      write(6,*)' first normalize'
      if (.not.int_normalize(sp_basis)) stop ' norm error 1'
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0)
      if (.not.bas_numbf(sp_basis,nbf_sp)) call errquit
     &      ('numbf failed',911)
      if (.not.bas_numcont(sp_basis,nshell_sp)) call errquit
     &      ('numcont error',911)
      write(6,*)' sp_basis b4 rak_ovlap'
      if (.not. bas_print(sp_basis))
     $      call errquit(' sp_basis print failed', 0)
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
      if (.not.status) call errquit('ma pop fail',911)
c
      if(.not.bas_destroy(basis)) call errquit
     &      ('basis bas_destroy failed',911)
      if(.not.bas_destroy(sp_basis)) call errquit
     &      ('sp_basis bas_destroy failed',911)
      if(.not.geom_destroy(geom)) call errquit
     &      ('geom_destroy failed',911)
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
#include "bas.fh"
      integer nshell, nbf_chk, nbf_sp, basis
      integer cn2bfr(2,nshell)
      integer type, nprim, ngen, spsp, ish
c
      nbf_chk = 0
      nbf_sp  = 0
      do ish = 1,nshell
        if(.not.bas_continfo(basis,ish,type,nprim,ngen,spsp))
     &        call errquit(' continfo failed ',911)
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
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'fuck i'
        nbfi    = ihi    - ilow + 1
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'fuck j'
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
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'fuck i'
        nbfi    = ihi    - ilow + 1
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'fuck j'
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
      if (.not.bas_high_angular(basis,hi_ang)) stop ' fuck ang'
      st_ang = hi_ang/2
*      if (.not.spcart_init(st_ang,.true.,.false.)) stop ' fuck sp'
*      if (.not.spcart_init(hi_ang,.true.,.false.)) stop ' fuck sp'
c
*      write(6,*)' mscr = ',mscr
*      write(6,*)' mbuf = ',mbuf
c
      call dfill((nbf_sp*nbf_sp),0.0d00,s_sp,1)
      call dfill((nbf*nbf),0.0d00,s,1)
      call dfill(mscr,0.0d00,scr,1)
      call dfill(mbuf,0.0d00,buf,1)

      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'fuck i'
        ilow_sp = cn2bfr_sp(1,ish)
        ihi_sp  = cn2bfr_sp(2,ish)
        nbfi    = ihi    - ilow + 1
        nbfi_sp = ihi_sp - ilow_sp + 1
        if(.not.bas_continfo(basis,ish,typei,nprim,igen,spsp))
     &        call errquit(' continfo failed ',911)
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'fuck j'
          jlow_sp = cn2bfr_sp(1,jsh)
          jhi_sp  = cn2bfr_sp(2,jsh)
          nbfj    = jhi    - jlow + 1
          nbfj_sp = jhi_sp - jlow_sp + 1
          if(.not.bas_continfo(basis,jsh,typej,nprim,jgen,spsp))
     &          call errquit(' continfo failed ',911)
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
*rak:          write(6,*) ' ma fucked 1'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 1'
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
*rak:          write(6,*) ' ma fucked 2'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 2'
c.... buf is now -- buf(nbfj,nbfi)
          write(6,*)' integral buffer cart,cart '
          call output(buf,1,nbfj,1,nbfi,nbfj,nbfi,1)
          call spcart_a_s(buf,scr,nbfj,typei,.false.,.false.)
c.... scr is now -- scr(nbfj,nbfi_sp)
          write(6,*)' integral buffer  cart,sph'
          call output(scr,1,nbfj,1,nbfi_sp,nbfj,nbfi_sp,1)
*rak:          write(6,*) ' ma fucked 3'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 3'
          call spcart_s_a(scr,buf,nbfi_sp,typej,.false.,.false.)
c.... buf is now -- buf(nbfj_sp,nbfi_sp)
          write(6,*)' integral buffer  sph,sph'
          call output(buf,1,nbfj_sp,1,nbfi_sp,nbfj_sp,nbfi_sp,1)
*rak:          write(6,*) ' ma fucked 4'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 4'
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
*rak:          write(6,*) ' ma fucked 5'
*rak:          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 5'
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
*      if (.not.spcart_terminate()) stop ' fuck sp term'
      end
      subroutine raktest_3ctr(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "bas.fh"
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
     &    ('raktest_3ctr: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('raktest_3ctr: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('raktest_3ctr: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('raktest_3ctr: basis load ',911)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(basis)) call errquit
     &    ('raktest_3ctr: error normalizing ',911)
c
      call int_init(rtdb,1,basis)
      memscr = 100 000
      membuf = 1000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911)
      if (.not.ma_push_get(mt_dbl,membuf,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911)
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
      if (.not.status) call errquit('pop failed',911)
      status = bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      if (.not.status) call errquit('b/g destroy failed',911)
      end
      subroutine raktest_stpr(rtdb)
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
      logical  stpr_walk
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
     &         call errquit('raktest_stpr: geom_create?', 911)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &         call errquit('raktest_stpr: no geometry ', 911)
c get number of atoms = nat
        if (.not. geom_ncent(geom,nat))
     &      call errquit('raktest_stpr: geom_ncent?',911)
        if (.not. geom_destroy(geom))
     &      call errquit('raktest_stpr: geom_destroy?',911)
        if (.not.
     &         MA_Push_Get(MT_DBL,(3*nat),'stpr fake gradient',
     &         h_grad,k_grad))
     &         call errquit
     &         ('raktest_stpr: allocation for gradient failed?',911)
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
     &      ('raktest_stpr: failed to read converged in rtdb', 911)
c put scf:energy   =  real value
        if (.not. rtdb_put(rtdb,'scf:energy', MT_DBL, 1, energy))
     &      call errquit
     &      ('raktest_stpr: failed to read energy in rtdb', 911)
c put scf:gradient = 3*nat reals
        if (.not. rtdb_put(rtdb, 'scf:gradient', MT_DBL, 
     &       (3*nat),dbl_mb(k_grad)))
     &      call errquit
     &        ('raktest_stpr: reading gradients failed',911)
c free memory
        if (.not. ma_pop_stack(h_grad))
     &      call errquit('raktest_stpr: pop failed',911)
      endif
      call ga_sync()
      call ga_sync()
      call ga_sync()
      if (stpr_walk('scf',rtdb)) then
        write(6,*)' walk converged'
      else
        write(6,*)' walk NOT converged'
      endif
      end
      subroutine raktest_initd(rtdb)
c raktest = 4
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "basP.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer geom
      integer mx1e, mxg, mxs1, mxs2
c
c
      integer nbas, bases(6)
c
      if (.not.bas_rtdb_in(rtdb))
     &    call errquit('raktest4: error loading known basis sets',911)
c
      write(6,*)' number of basis sets in rtdb ',nbasis_rtdb
c
      do 00100 nbas = 1,nbasis_rtdb
        write(6,*)' basis ',nbas,' is ',bs_names_rtdb(nbas)
00100 continue
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911)
c
      do 00200 nbas = 1,nbasis_rtdb
        if(.not.bas_create(bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_create choked',911)
        if(.not.
     &      bas_rtdb_load
     &      (rtdb,geom,bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_rtdb_load choked',911)
00200 continue
      call intd_init(rtdb,nbasis_rtdb,bases)
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
     &        ('raktest_initd: bas_destroy failed',911)
      enddo
      if (.not.geom_destroy(geom)) call errquit
     &      ('raktest_initd: _destroy failed',911)
c
      call int_terminate()
c
      end
      subroutine raktest_init(rtdb)
c raktest = 5
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "basP.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer geom
      integer mx1e, mxg, mxs1, mxs2
c
c
      integer nbas, bases(6)
c
      if (.not.bas_rtdb_in(rtdb))
     &    call errquit('raktest4: error loading known basis sets',911)
c
      write(6,*)' number of basis sets in rtdb ',nbasis_rtdb
c
      do 00100 nbas = 1,nbasis_rtdb
        write(6,*)' basis ',nbas,' is ',bs_names_rtdb(nbas)
00100 continue
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &    call errquit('raktest4: geom_create failed',911)
c
      do 00200 nbas = 1,nbasis_rtdb
        if(.not.bas_create(bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_create choked',911)
        if(.not.
     &      bas_rtdb_load
     &      (rtdb,geom,bases(nbas),bs_names_rtdb(nbas)))
     &      call errquit('raktest4: bas_rtdb_load choked',911)
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
     &        ('raktest_init: bas_destroy failed',911)
      enddo
      if (.not.geom_destroy(geom)) call errquit
     &      ('raktest_init: _destroy failed',911)
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
      if (KEI) KEI = .false.
      if (O2I) O2I = .false.

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
     &      call errquit('raktest_gc: geom_create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &      call errquit('raktest_gc: geom_rtdb_load failed',911)
c
      if(.not.bas_create(basisseg,'aobasisseg'))
     &      call errquit('raktest_gc: bas_create choked',911)
      if(.not.
     &      bas_rtdb_load(rtdb,geom,basisseg,'aobasisseg'))
     &      call errquit('raktest_gc: bas_rtdb_load choked',911)
c
      if(.not.bas_create(basisns,'aobasisns'))
     &      call errquit('raktest_gc: bas_create choked',911)
      if(.not.
     &      bas_rtdb_load(rtdb,geom,basisns,'aobasisns'))
     &      call errquit('raktest_gc: bas_rtdb_load choked',911)
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
      if(.not.int_normalize(basisseg)) call errquit
     &      ('raktest_gc: int_normalize seg failed',911)
      if(.not.int_normalize(basisns)) call errquit
     &      ('raktest_gc: int_normalize ns failed',911)
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
     &      call errquit('raktest_gc: ma failed 1',911)
      if(.not.ma_push_get(mt_dbl,memscr,'int scr',h_scr,k_scr))
     &      call errquit('raktest_gc: ma failed 2',911)
c
      if(.not.bas_numcont(basisseg,ncont_seg))
     &      call errquit('raktest_gc: bas_numcont failed',911)
      write(luout,*)' ao basis seg number of contractions    :',
     &      ncont_seg
      if(.not.bas_numbf(basisseg,nbf_seg))
     &      call errquit('raktest_gc: bas_numbf failed',911)
      write(luout,*)' ao basis seg number of basis functions :',
     &      nbf_seg
      if(.not.bas_numcont(basisns,ncont_ns))
     &      call errquit('raktest_gc: bas_numcont failed',911)
      write(luout,*)' ao basis noseg number of contractions    :',
     &      ncont_ns
      if(.not.bas_numbf(basisns,nbf_ns))
     &      call errquit('raktest_gc: bas_numbf failed',911)
      write(luout,*)' ao basis noseg number of basis functions :',
     &      nbf_ns
c
      if (nbf_ns.ne.nbf_seg) call errquit
     &      ('raktest_gc: nbf_seg.ne.nbf_ns',(nbf_ns-nbf_seg))
c
      if(.not.ma_push_get(mt_dbl,(nbf_seg*nbf_seg),'overlap seg',
     &      h_s_seg,k_s_seg))
     &      call errquit('raktest_gc: ma failed 3',911)
      if(.not.ma_push_get(mt_dbl,(nbf_ns*nbf_ns),'overlap noseg',
     &      h_s_ns,k_s_ns))
     &      call errquit('raktest_gc: ma failed 4',911)
      if(.not.ma_push_get(mt_dbl,(nbf_ns*nbf_ns),'difference matrix',
     &      h_diff,k_diff))
     &      call errquit('raktest_gc: ma failed 4',911)

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
      if (.not.status) call errquit('raktest_gc: pop error',911)
c
      if (.not.bas_destroy(basisseg)) call errquit
     &    ('raktest_gc: bas_destroy failed ?',911)
      if (.not.bas_destroy(basisns)) call errquit
     &    ('raktest_gc: bas_destroy failed ?',911)
      if (.not.geom_destroy(geom)) call errquit
     &    ('raktest_gc: geom_destroy failed ?',911)
c
      call int_terminate()
c
      end
      subroutine raktest_bs_gc(basis,nscr,scr,nbuf,buf,nbf,ncont,S,
     &      OV,KE,PE,msg)
      implicit none
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
     &          ('raktest_gc: cn2bfr failed',911)
*          write(luout,*)'ish = ',ish,' bfr =',ilo,ihi
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     &          ('raktest_gc: cn2bfr failed',911)
*          write(luout,*)'jsh = ',jsh,' bfr =',jlo,jhi
          nint = (ihi-ilo+1)*(jhi-jlo+1)
*          write(luout,*)'number of integrals = ',nint
          if (OV) then
            call int_1eov(basis,ish,basis,jsh,
     &          nscr,scr,nint,buf)
          elseif (KE) then
            call int_1eke(basis,ish,basis,jsh,
     &          nscr,scr,nint,buf)
          elseif (PE) then
            call int_1epe(basis,ish,basis,jsh,
     &          nscr,scr,nint,buf)
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
     &    call errquit('raktest_geomwrt: geom_create failed?',911)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911)
c
      if (.not.bas_create(basis,'ao sp_basis')) call errquit
     &      ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao sp_basis'))
     &    call errquit
     &      ('bas_rtdb_load failed',911)
      if (.not.int_normalize(basis)) stop ' norm error 1'
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
     &      call errquit('ma scr failed',911)
c
      if (.not.ma_push_get(mt_dbl,m_buf,'buf for 1e',h_buf,k_buf))
     &      call errquit('ma buf failed',911)
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
      if (.not.status) call errquit('ma pop fail',911)
      call int_terminate()
c
      if (.not.bas_destroy(basis)) call errquit
     &    ('9: bas_destroy failed ?',911)
      if (.not.geom_destroy(geom)) call errquit
     &    ('9: geom_destroy failed ?',911)
c
      end
      subroutine raktest_ecp(rtdb)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
      integer rtdb
*ecp:c
*ecp:      integer geom, basis, ecpid
*ecp:c
*ecp:      logical int_normalize, int_ecp_init
*ecp:      external int_normalize, int_ecp_init
*ecp:c
*ecp:*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 1?'
*ecp:c
*ecp:      if (.not.geom_create(geom,'geometry'))
*ecp:     &    call errquit('raktest_geomwrt: geom_create failed?',911)
*ecp:      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
*ecp:     &    ('raktest_geomwrt: geom_rtdb_load -ref failed',911)
*ecp:      if (.not.geom_print(geom)) stop ' print error'
*ecp:c
*ecp:      if (.not.bas_create(basis,'ao basis')) call errquit
*ecp:     &      ('bas_create failed',911)
*ecp:      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
*ecp:     &    call errquit
*ecp:     &      ('bas_rtdb_load failed',911)
*ecp:      if (.not.bas_print(basis)) stop ' print error'
*ecp:*      if (.not.rtdb_print(rtdb,.true.)) stop ' rtdb_p 2?'
*ecp:c
*ecp:      if (.not.bas_create(ecpid,'ecp basis')) call errquit
*ecp:     &      ('bas_create failed',911)
*ecp:      if (.not.bas_rtdb_load(rtdb,geom,ecpid,'ecp basis'))
*ecp:     &    call errquit
*ecp:     &      ('bas_rtdb_load failed',911)
*ecp:      if (.not.bas_print(ecpid)) stop ' print error'
*ecp:c
*ecp:      if (.not.int_normalize(basis)) stop ' norm error 1'
*ecp:      if (.not.bas_print(basis)) stop ' print error'
*ecp:      if (.not.gbs_map_print(basis)) stop ' gbs map print error'
*ecp:      if (.not.bas_print(ecpid)) stop ' print error'
*ecp:      if (.not.gbs_map_print(basis)) stop ' gbs map print error'
*ecp:c
*ecp:      if (.not.int_ecp_init(ecpid)) stop ' error in int_ecp_init'
*ecp:      call int_ecp_terminate()
*ecp:c
*ecp:      if (.not.bas_destroy(basis)) stop ' bas_dest error 1'
*ecp:      if (.not.bas_destroy(ecpid)) stop ' bas_dest error 2'
*ecp:      if (.not.geom_destroy(geom)) stop ' geom_dest error 1'
*ecp:      if (.not.bas_version()) stop ' bas_version error'
*ecp:c
      end
