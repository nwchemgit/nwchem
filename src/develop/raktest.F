      subroutine raktest(rtdb)
      implicit none
c $Id: raktest.F,v 1.44 1996-01-22 21:51:06 d3e129 Exp $
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
      goto (00100,00200,00300,00400,00500) raktask
c
      write(6,*)' default raktest task '
      write(6,*)' test of general contraction code '
      call raktest_gc(rtdb)
      goto 00001
c
00100 write(6,*)' raktest task 1 stepper test'
      call raktest_stpr(rtdb)
      goto 00001
00200 write(6,*)' raktest task 2 check intd_init'
      call raktest_initd(rtdb)
      goto 00001
00300 write(6,*)' raktest task 3 check int_init'
      call raktest_init(rtdb)
      goto 00001
00400 write(6,*)' raktest check 3ctr nai'
      call raktest_3ctr(rtdb)
      goto 00001
00500 write(6,*)' test of spcart stuff '
      call raktest_spcart(rtdb)
c
00001 continue
      end
      subroutine raktest_spcart(rtdb)
      implicit none 
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer rtdb ! [input] rtdb handle
c
      integer geom, basis
      integer nbf, nbfsq, nbf_sp, nbf_chk, nshell
      integer max1e, max2e, mscr1, mscr2
      integer h_cart_s, h_sph_s, h_scr, h_buf, h_2bfr
      integer k_cart_s, k_sph_s, k_scr, k_buf, k_2bfr
      logical status
c
      logical int_normalize
      external int_normalize
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
      if (.not.int_normalize(basis)) stop ' norm error '
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
c
      if (.not.bas_numcont(basis,nshell)) call errquit
     &      ('numcont error',911)
c
      call int_init(rtdb,1,basis)
      call int_mem(max1e,max2e,mscr1,mscr2)
      max1e = max1e*2
      mscr1 = mscr1*2
c
      if (.not.ma_push_get(mt_dbl,mscr1,'scr for 1e',h_scr,k_scr))
     &      call errquit('ma scr failed',911)
c
      if (.not.ma_push_get(mt_dbl,max1e,'buf for 1e',h_buf,k_buf))
     &      call errquit('ma buf failed',911)
c
      if (.not.ma_push_get(mt_int,2*nshell,'buf for sp cn2bfr',
     &      h_2bfr,k_2bfr))
     &      call errquit('ma buf failed',911)
c
      call rak_tospbfr(basis,nshell,nbf_chk,nbf_sp,int_mb(k_2bfr))
c
      if (nbf_chk.ne.nbf) then
        write(6,*)' nbf not right ',nbf_chk, nbf
      endif
c
      write(6,*)' nbf    ',nbf
      write(6,*)' nbf_sp ',nbf_sp
c
      call rak_ovlap_test_sp(basis,nbf,nbf_sp,nshell,
     &      dbl_mb(k_scr),mscr1,
     &      dbl_mb(k_buf),max1e,
     &      dbl_mb(k_cart_s),dbl_mb(k_sph_s),
     &      int_mb(k_2bfr))
c
      call int_terminate()
      status = .true.
      status = status.and.ma_pop_stack(h_2bfr)
      status = status.and.ma_pop_stack(h_buf)
      status = status.and.ma_pop_stack(h_scr)
      status = status.and.ma_pop_stack(h_sph_s)
      status = status.and.ma_pop_stack(h_cart_s)
c
      if (.not.status) call errquit('ma pop fail',911)
c
      end
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
        write(6,10000)ish,type,cn2bfr(1,ish),cn2bfr(2,ish)
      enddo
10000 format('<ish:type><',i3,':',i3,'>  range ',i3,' to ',i3)
      end
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
      integer typei, typej, nprim, ngen, spsp
      integer nint, nint_sp, count, hi_ang, st_ang
      integer ii, jj
      double precision value
c
      logical spcart_init, spcart_terminate
      external spcart_init, spcart_terminate
c
      write(6,*)'inside rak_ovlap'
      if (.not.bas_high_angular(basis,hi_ang)) stop ' fuck ang'
      st_ang = hi_ang/2
      if (.not.spcart_init(st_ang,.true.,.false.)) stop ' fuck sp'
      if (.not.spcart_init(hi_ang,.true.,.false.)) stop ' fuck sp'
c
      write(6,*)' mscr = ',mscr
      write(6,*)' mbuf = ',mbuf
c
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilow,ihi)) stop 'fuck i'
        ilow_sp = cn2bfr_sp(1,ish)
        ihi_sp  = cn2bfr_sp(2,ish)
        nbfi    = ihi    - ilow + 1
        nbfi_sp = ihi_sp - ilow_sp + 1
        if(.not.bas_continfo(basis,ish,typei,nprim,ngen,spsp))
     &        call errquit(' continfo failed ',911)
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlow,jhi)) stop 'fuck j'
          jlow_sp = cn2bfr_sp(1,jsh)
          jhi_sp  = cn2bfr_sp(2,jsh)
          nbfj    = jhi    - jlow + 1
          nbfj_sp = jhi_sp - jlow_sp + 1
          if(.not.bas_continfo(basis,jsh,typej,nprim,ngen,spsp))
     &          call errquit(' continfo failed ',911)
c
          nint    = nbfi*nbfj
          nint_sp = nbfi_sp*nbfj_sp
c
          write(6,*)' '
          write(6,*)'<ish,ilow,ihi,nbfi,typei>',ish,ilow,ihi,nbfi,typei
          write(6,*)'<jsh,jlow,jhi,nbfj,typej>',jsh,jlow,jhi,nbfj,typej
          write(6,*)' nint     = ',nint,' nint(sp) = ',nint_sp
          write(6,*)'<ish,ilowsp,ihisp,nbfisp,typei>',
     &          ish,ilow_sp,ihi_sp,nbfi_sp,typei
          write(6,*)'<jsh,jlowsp,jhisp,nbfjsp,typej>',
     &          jsh,jlow_sp,jhi_sp,nbfj_sp,typej
          write(6,*) ' ma fucked 1'
          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 1'
          write(6,*)' '
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
          write(6,*) ' ma fucked 2'
          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 2'
c.... buf is now -- buf(nbfj,nbfi)
          write(6,*)' integral buffer cart,cart '
          call output(buf,1,nbfj,1,nbfi,nbfj,nbfi,1)
          call spcart_a_s(buf,scr,nbfj,typei,.false.)
c.... scr is now -- scr(nbfj,nbfi_sp)
          write(6,*)' integral buffer  cart,sph'
          call output(scr,1,nbfj,1,nbfi_sp,nbfj,nbfi_sp,1)
          write(6,*) ' ma fucked 3'
          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 3'
          call spcart_s_a(scr,buf,nbfi_sp,typej,.false.)
c.... buf is now -- buf(nbfj_sp,nbfi_sp)
          write(6,*)' integral buffer  sph,sph'
          call output(buf,1,nbfj_sp,1,nbfi_sp,nbfj_sp,nbfi_sp,1)
          write(6,*) ' ma fucked 4'
          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 4'
          count = 0
          do ibf = ilow_sp,ihi_sp
            do jbf = jlow_sp, jhi_sp
              count = count + 1
              value = buf(count)
              s_sp(ibf,jbf) = value
              s_sp(jbf,ibf) = value
            enddo
          enddo
          write(6,*) ' ma fucked 5'
          if (.not.ma_verify_allocator_stuff()) stop ' ma fucked 5'
        enddo
      enddo
      write(6,*)' loops done '
      write(6,*)' cartesian overlap matrix '
      call output(s,1,nbf,1,nbf,nbf,nbf,1)
      write(6,*)' spherical overlap matrix '
      call output(s_sp,1,nbf_sp,1,nbf_sp,nbf_sp,nbf_sp,1)
      do ii = 2, nbf_sp
        do jj = 1,(ii-1)
          write(6,'(a,2i5,2f15.8)')
     &          ' ratios ',ii,jj,(s_sp(ii,ii)/s_sp(jj,jj)),
     &          (s_sp(jj,jj)/s_sp(ii,ii))
        enddo
      enddo
      if (.not.spcart_terminate()) stop ' fuck sp term'
      end
*rak:      subroutine raktest0(rtdb)
*rak:      implicit none
*rak:c test0 "default raktest"
*rak:#include "global.fh"
*rak:      integer rtdb
*rak:c::functions
*rak:      logical gen_apphess_at
*rak:      external gen_apphess_at
*rak:c
*rak:      double precision delta
*rak:c
*rak:#if defined(SUN)
*rak:      real etime,tarray(2)
*rak:      external etime
*rak:#endif
*rak:      real t1
*rak:      t1 = 0.0
*rak:c
*rak:c      logical status
*rak:c
*rak:      if(ga_nodeid().eq.0) then
*rak:        write(6,*) ' rak test'
*rak:      endif
*rak:      delta = 0.01d00
*rak:#if defined(SUN)
*rak:      t1 = etime(tarray)
*rak:#endif
*rak:      if(.not.gen_apphess_at(rtdb,delta)) call errquit
*rak:     &    ('raktest:gen_apphess_at failed',911)
*rak:#if defined(SUN)
*rak:      t1 = etime(tarray) - t1 
*rak:#endif
*rak:      if(ga_nodeid().eq.0) then
*rak:        write(6,'(1x,a,f8.2)')
*rak:     &      'time for rak test',t1
*rak:      endif
*rak:      end
*rak:      subroutine raktest2(rtdb)
*rak:      implicit none
*rak:c test0 "default raktest"
*rak:#include "global.fh"
*rak:      integer rtdb
*rak:c::functions
*rak:      logical gen_apphess_opt
*rak:      external gen_apphess_opt
*rak:c
*rak:      double precision delta
*rak:c
*rak:#if defined(SUN)
*rak:      real etime,tarray(2)
*rak:      external etime
*rak:#endif
*rak:      real t1
*rak:      t1 = 0.0
*rak:c
*rak:c      logical status
*rak:c
*rak:      if(ga_nodeid().eq.0) then
*rak:        write(6,*) ' rak test'
*rak:      endif
*rak:      delta = 0.01d00
*rak:#if defined(SUN)
*rak:      t1 = etime(tarray)
*rak:#endif
*rak:      if(.not.gen_apphess_opt(rtdb,delta)) call errquit
*rak:     &    ('raktest:2:gen_apphess_opt failed',911)
*rak:#if defined(SUN)
*rak:      t1 = etime(tarray) - t1 
*rak:#endif
*rak:      if(ga_nodeid().eq.0) then
*rak:        write(6,'(1x,a,f8.2)')
*rak:     &      'time for rak test 2',t1
*rak:      endif
*rak:      end
*rak:      subroutine raktest1(rtdb)
*rak:      implicit none
*rak:      integer rtdb
*rak:      call deriv3c(rtdb)
*rak:      call deriv2c(rtdb)
*rak:      end
*rak:      subroutine deriv3c(rtdb)
*rak:      implicit none
*rak:#include "mafdecls.fh"
*rak:#include "rtdb.fh"
*rak:#include "context.fh"
*rak:#include "geom.fh"
*rak:#include "bas.fh"
*rak:c
*rak:      double precision ddot
*rak:      external ddot
*rak:c
*rak:      integer rtdb
*rak:c
*rak:      character*100 mo_basis, geom_name
*rak:      integer basis, geom
*rak:      integer nbf, nshells, nat
*rak:      integer size2e, size2e_d
*rak:      integer h2e3d,i2e3d, h2e3d_fd,i2e3d_fd
*rak:      integer hxyz,ixyz,nbuf,nscr,hscr,iscr,hbuf,ibuf
*rak:      logical status
*rak:c
*rak:      if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
*rak:     &    mo_basis = 'ao basis'
*rak:      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
*rak:     &    geom_name = 'geometry'
*rak:c
*rak:      if(.not.geom_create(geom,geom_name))call errquit
*rak:     &    ('deriv3c: geom create error',911)
*rak:      if(.not.bas_create(basis,mo_basis))call errquit
*rak:     &    ('deriv3c: basis create error',911)
*rak:c
*rak:      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
*rak:     &    ('deriv3c: geom load ',911)
*rak:      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
*rak:     &    ('deriv3c: basis load ',911)
*rak:c
*rak:      if(.not.bas_numbf(basis,nbf)) call errquit
*rak:     &    ('deriv3c: bas_numbf failed',911)
*rak:      if(.not.bas_numcont(basis,nshells)) call errquit
*rak:     &       ('deriv3c:bas_numcont failed',911)
*rak:      if(.not.geom_ncent(geom,nat)) call errquit
*rak:     &       ('deriv3c:geom_ncent failed',911)
*rak:c
*rak:      size2e   = nbf*nbf*nbf
*rak:      size2e_d = size2e*3*nat
*rak:      write(6,*)' size2e   = ',size2e
*rak:      write(6,*)' size2e_d = ',size2e_d
*rak:c
*rak:      call intd_init(rtdb,1,basis)
*rak:*      call int_normalize(basis)
*rak:      call int_norm_2c(basis)
*rak:c
*rak:      if(.not.ma_set_error_print(.true.)) stop ' fuck 1'
*rak:      if(.not.ma_verify_allocator_stuff()) stop ' fuck 2'
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:      if (.not.ma_push_get(mt_dbl,size2e_d,'2e3c d real',
*rak:     &    h2e3d,i2e3d))
*rak:     &    call errquit('deriv3c: ma fail 1',911)
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:      if (.not.ma_push_get(mt_dbl,size2e_d,'2e3c d fd',
*rak:     &    h2e3d_fd,i2e3d_fd))
*rak:     &    call errquit('deriv3c: ma fail 2',911)
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:      if (.not.ma_push_get(mt_dbl,(3*nat),'xyzbuf',
*rak:     &    hxyz,ixyz))
*rak:     &    call errquit('deriv3c: ma fail 4',911)
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:      call int_mem_2e4c(nbuf,nscr)
*rak:      if (.not.ma_push_get(mt_dbl,nscr,'scratch',
*rak:     &    hscr,iscr))
*rak:     &    call errquit('deriv3c: ma fail 5',911)
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:      if (.not.ma_push_get(mt_dbl,(2*nbuf),'buffer scratch',
*rak:     &    hbuf,ibuf))
*rak:     &    call errquit('deriv3c: ma fail 6',911)
*rak:      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
*rak:      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
*rak:c
*rak:      call dfill(size2e_d, 0.0d00,dbl_mb(i2e3d)     ,1)
*rak:      call dfill(size2e_d, 0.0d00,dbl_mb(i2e3d_fd)  ,1)
*rak:      call dfill(nscr,     0.0d00,dbl_mb(iscr)      ,1)
*rak:      call dfill(2*nbuf,   0.0d00,dbl_mb(ibuf)      ,1)
*rak:      call dfill((3*nat),  0.0d00,dbl_mb(ixyz)      ,1)
*rak:c      
*rak:c      
*rak:      call build_3ed(dbl_mb(i2e3d),nat,size2e,nbf,nshells,
*rak:     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,basis,geom)
*rak:      call dfill(nscr,      0.0d00,dbl_mb(iscr)      ,1)
*rak:      call dfill(2*nbuf,    0.0d00,dbl_mb(ibuf)      ,1)
*rak:      call build_3edfd(dbl_mb(i2e3d_fd),nat,size2e,nbf,
*rak:     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,dbl_mb(ixyz),
*rak:     &    1.0d-05,basis,geom,nshells)
*rak:      call daxpy(size2e_d,-1.0d00,dbl_mb(i2e3d),1,dbl_mb(i2e3d_fd),1)
*rak:c
*rak:      call print_2d3c(dbl_mb(i2e3d),dbl_mb(i2e3d_fd),nbf,nat)
*rak:c
*rak:      write(6,*)
*rak:      write(6,*)' 2e3c deriv diff norm ',
*rak:     &    ddot(size2e_d,dbl_mb(i2e3d_fd),1,dbl_mb(i2e3d_fd),1)
*rak:c
*rak:      call intd_terminate()
*rak:      status = ma_pop_stack(hbuf)
*rak:      status = status.and.ma_pop_stack(hscr)
*rak:      status = status.and.ma_pop_stack(hxyz)
*rak:      status = status.and.ma_pop_stack(h2e3d_fd)
*rak:      status = status.and.ma_pop_stack(h2e3d)
*rak:      if (.not.status) call errquit('pop failed',911)
*rak:      status = bas_destroy(basis)
*rak:      status = status.and.geom_destroy(geom)
*rak:      if (.not.status) call errquit('b/g destroy failed',911)
*rak:      end
*rak:      subroutine build_3edfd(x2e3dfd,nat,size2e,nbf,scr,nscr,
*rak:     &    buf,nbuf,xyz,delta,basis,geom,nshells)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:      integer nat, size2e, nbf, nscr, nbuf, basis, geom, nshells
*rak:      double precision x2e3dfd(nbf,nbf,nbf,3,nat)
*rak:      double precision scr(nscr)
*rak:      double precision buf(nbuf,2)
*rak:      double precision xyz(3,nat)
*rak:      double precision delta
*rak:c
*rak:      double precision scale
*rak:      integer ish, ibflo, ibfhi, ibf
*rak:      integer jsh, jbflo, jbfhi, jbf
*rak:      integer ksh, kbflo, kbfhi, kbf
*rak:      integer datom, dxyz
*rak:      integer nint, icount
*rak:      integer int_nint
*rak:      external int_nint
*rak:c
*rak:c.. store original coords in xyz
*rak:      call dcopy((3*nat),coords(1,1,geom),1,xyz,1)
*rak:      scale = 1.0d00/(2.0d00*delta)
*rak:c
*rak:      do 00100 ish = 1,nshells
*rak:        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &      call errquit('cn2bfr error',ish)
*rak:        do 00200 jsh = 1,ish
*rak:          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &        call errquit('cn2bfr error',jsh)
*rak:          do 00300 ksh = 1,jsh
*rak:            if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
*rak:     &          call errquit('cn2bfr error',ksh)
*rak:            nint = int_nint(basis,ish,basis,jsh,basis,ksh,basis,0)
*rak:            do 00400 datom = 1,nat
*rak:              do 00500 dxyz = 1,3
*rak:c... copy original coords in to coords array
*rak:                call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:c... add delta to coordinate
*rak:                coords(dxyz,datom,geom)=coords(dxyz,datom,geom)+delta
*rak:                call int_2e3c(basis,ish,basis,jsh,ksh,
*rak:     &              nscr,scr,nint,buf(1,1))
*rak:c... copy original coords in to coords array
*rak:                call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:c... subtract delta to coordinate
*rak:                coords(dxyz,datom,geom)=coords(dxyz,datom,geom)-delta
*rak:                call int_2e3c(basis,ish,basis,jsh,ksh,
*rak:     &              nscr,scr,nint,buf(1,2))
*rak:                call daxpy(nint,-1.0d00,buf(1,2),1,buf(1,1),1)
*rak:                call dscal(nint,scale,buf(1,1),1)
*rak:                icount = 0
*rak:                do 00600 ibf = ibflo,ibfhi
*rak:                  do 00700 jbf = jbflo,jbfhi
*rak:                    do 00800 kbf = kbflo, kbfhi
*rak:                      icount = icount + 1
*rak:                      x2e3dfd(ibf,jbf,kbf,dxyz,datom) =
*rak:     &                    buf(icount,1)
*rak:00800               continue
*rak:00700             continue
*rak:00600           continue
*rak:00500         continue
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:c      
*rak:c... copy original coords in to coords array
*rak:      call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:c
*rak:      end
*rak:      subroutine build_3ed(x2e3d,nat,size2e,nbf,nshells,
*rak:     &    scr,nscr,buf,nbuf,basis,geom)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:      integer nat, size2e, nbf, nshells, nscr, nbuf, basis, geom
*rak:      double precision x2e3d(nbf,nbf,nbf,3,nat)
*rak:      double precision scr(nscr)
*rak:      double precision buf(nbuf)
*rak:c
*rak:      integer idatom(4)
*rak:      integer ish, ibflo, ibfhi, num_ibf, ibf
*rak:      integer jsh, jbflo, jbfhi, num_jbf, jbf
*rak:      integer ksh, kbflo, kbfhi, num_kbf, kbf
*rak:      integer nint, icount, id, dxyz
*rak:c
*rak:      integer int_nint
*rak:      external int_nint
*rak:c
*rak:      do 00100 ish = 1,nshells
*rak:        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &      call errquit('cn2bfr error',ish)
*rak:        num_ibf = ibfhi-ibflo+1
*rak:        do 00200 jsh = 1,ish
*rak:          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &        call errquit('cn2bfr error',jsh)
*rak:          num_jbf = jbfhi-jbflo+1
*rak:          do 00300 ksh = 1,jsh
*rak:            if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
*rak:     &          call errquit('cn2bfr error',ksh)
*rak:            num_kbf = kbfhi-kbflo+1
*rak:            nint = int_nint(basis,ish,basis,jsh,basis,ksh,basis,0)
*rak:            call intd_2e3c(basis,ish,basis,jsh,ksh,nscr,scr,(nint*3*4),
*rak:     &          buf,idatom)
*rak:            write(6,*)' idatom 2e3c',idatom
*rak:            icount = 0
*rak:            do 00400 id = 1,4
*rak:              if (idatom(id).gt.0) then
*rak:                do 00500 dxyz = 1,3
*rak:                  do 00600 ibf = ibflo,ibfhi
*rak:                    do 00700 jbf = jbflo,jbfhi
*rak:                      do 00800 kbf = kbflo,kbfhi
*rak:                        icount = icount + 1
*rak:                        x2e3d(ibf,jbf,kbf,dxyz,idatom(id))=buf(icount)
*rak:00800                 continue
*rak:00700               continue
*rak:00600             continue
*rak:00500           continue
*rak:              else
*rak:                icount = icount + 3*num_ibf*num_jbf*num_kbf
*rak:              endif
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:      end
*rak:      subroutine deriv2c(rtdb)
*rak:      implicit none
*rak:#include "mafdecls.fh"
*rak:#include "rtdb.fh"
*rak:#include "context.fh"
*rak:#include "geom.fh"
*rak:#include "bas.fh"
*rak:c
*rak:      double precision ddot
*rak:      external ddot
*rak:c
*rak:      integer rtdb
*rak:c
*rak:      character*100 mo_basis, geom_name
*rak:      integer basis, geom
*rak:      integer nbf, nshells, nat
*rak:      integer size2e, size2e_d
*rak:      integer h2e2d,i2e2d, h2e2d_fd,i2e2d_fd, h2e2d_diff,i2e2d_diff
*rak:      integer hxyz,ixyz,nbuf,nscr,hscr,iscr,hbuf,ibuf
*rak:      logical status
*rak:c
*rak:      if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
*rak:     &    mo_basis = 'ao basis'
*rak:      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
*rak:     &    geom_name = 'geometry'
*rak:c
*rak:      if(.not.geom_create(geom,geom_name))call errquit
*rak:     &    ('deriv2c: geom create error',911)
*rak:      if(.not.bas_create(basis,mo_basis))call errquit
*rak:     &    ('deriv2c: basis create error',911)
*rak:c
*rak:      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
*rak:     &    ('deriv2c: geom load ',911)
*rak:      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
*rak:     &    ('deriv2c: basis load ',911)
*rak:c
*rak:      if(.not.bas_numbf(basis,nbf)) call errquit
*rak:     &    ('deriv2c: bas_numbf failed',911)
*rak:      if(.not.bas_numcont(basis,nshells)) call errquit
*rak:     &       ('deriv2c:bas_numcont failed',911)
*rak:      if(.not.geom_ncent(geom,nat)) call errquit
*rak:     &       ('deriv2c:geom_ncent failed',911)
*rak:c
*rak:      size2e   = nbf*(nbf+1)/2
*rak:      size2e_d = size2e*3*nat
*rak:c
*rak:      call intd_init(rtdb,1,basis)
*rak:*      call int_normalize(basis)
*rak:      call int_norm_2c(basis)
*rak:c
*rak:      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d real',
*rak:     &    h2e2d,i2e2d))
*rak:     &    call errquit('deriv2c: ma fail 1',911)
*rak:      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d fd',
*rak:     &    h2e2d_fd,i2e2d_fd))
*rak:     &    call errquit('deriv2c: ma fail 2',911)
*rak:      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d diff',
*rak:     &    h2e2d_diff,i2e2d_diff))
*rak:     &    call errquit('deriv2c: ma fail 3',911)
*rak:      if (.not.ma_push_get(mt_dbl,(3*nat),'xyzbuf',
*rak:     &    hxyz,ixyz))
*rak:     &    call errquit('deriv2c: ma fail 4',911)
*rak:      call int_mem_2e4c(nbuf,nscr)
*rak:      if (.not.ma_push_get(mt_dbl,nscr,'scratch',
*rak:     &    hscr,iscr))
*rak:     &    call errquit('deriv2c: ma fail 5',911)
*rak:      if (.not.ma_push_get(mt_dbl,(2*nbuf),'buffer scratch',
*rak:     &    hbuf,ibuf))
*rak:     &    call errquit('deriv2c: ma fail 6',911)
*rak:c
*rak:      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d)     ,1)
*rak:      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d_fd)  ,1)
*rak:      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d_diff),1)
*rak:      call dfill(nscr,     0.0d00,dbl_mb(iscr)      ,1)
*rak:      call dfill(2*nbuf,   0.0d00,dbl_mb(ibuf)      ,1)
*rak:      call dfill((3*nat),  0.0d00,dbl_mb(ixyz)      ,1)
*rak:c      
*rak:      call build_2ed(dbl_mb(i2e2d),nat,size2e,nbf,nshells,
*rak:     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,basis,geom)
*rak:      call dfill(nscr,      0.0d00,dbl_mb(iscr)      ,1)
*rak:      call dfill(2*nbuf,    0.0d00,dbl_mb(ibuf)      ,1)
*rak:      call build_2edfd(dbl_mb(i2e2d_fd),nat,size2e,nbf,
*rak:     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,dbl_mb(ixyz),
*rak:     &    1.0d-05,basis,geom,nshells)
*rak:      call dcopy(size2e_d,dbl_mb(i2e2d_fd),1,dbl_mb(i2e2d_diff),1)
*rak:      call daxpy(size2e_d,-1.0d00,dbl_mb(i2e2d),1,dbl_mb(i2e2d_diff),1)
*rak:c
*rak:      call print_2d2c(dbl_mb(i2e2d),dbl_mb(i2e2d_fd),
*rak:     &    dbl_mb(i2e2d_diff),size2e,nat)
*rak:c
*rak:      write(6,*)
*rak:      write(6,*)' 2e2c deriv diff norm ',
*rak:     &    ddot(size2e_d,dbl_mb(i2e2d_diff),1,dbl_mb(i2e2d_diff),1)
*rak:c
*rak:      call intd_terminate()
*rak:      status = ma_pop_stack(hbuf)
*rak:      status = status.and.ma_pop_stack(hscr)
*rak:      status = status.and.ma_pop_stack(hxyz)
*rak:      status = status.and.ma_pop_stack(h2e2d_diff)
*rak:      status = status.and.ma_pop_stack(h2e2d_fd)
*rak:      status = status.and.ma_pop_stack(h2e2d)
*rak:      if (.not.status) call errquit('pop failed',911)
*rak:      status = bas_destroy(basis)
*rak:      status = status.and.geom_destroy(geom)
*rak:      if (.not.status) call errquit('b/g destroy failed',911)
*rak:      end
*rak:      subroutine build_2ed(x2e2d,nat,size2e,nbf,nshells,
*rak:     &    scr,nscr,buf,nbuf, basis, geom)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:      integer basis, geom
*rak:      integer nat
*rak:      integer size2e
*rak:      integer nbf
*rak:      integer nshells
*rak:      double precision x2e2d(size2e,3,nat)
*rak:      integer nscr, nbuf, idatom(2), id, dxyz, index
*rak:      double precision scr(nscr),buf(nbuf)      
*rak:c
*rak:      integer ish, ibflo, ibfhi, num_ibf, ibf
*rak:      integer jsh, jbflo, jbfhi, num_jbf, jbf
*rak:      integer nint, icount
*rak:      integer int_nint
*rak:      external int_nint
*rak:      integer min, max
*rak:      intrinsic min,max
*rak:      integer i,j,isym2
*rak:      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
*rak:c
*rak:      do 00100 ish=1,nshells
*rak:        do 00200 jsh=1,ish
*rak:          nint = int_nint(basis,ish,basis,jsh,basis,0,basis,0)
*rak:          call intd_2e2c(basis,ish,basis,jsh,nscr,scr,(nint*3*4),
*rak:     &        buf,idatom)
*rak:          write(6,*)' idatom 2e2c', idatom
*rak:          if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &        call errquit('cn2bfr error',ish)
*rak:          num_ibf = ibfhi-ibflo+1
*rak:          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &        call errquit('cn2bfr error',ish)
*rak:          num_jbf = jbfhi-jbflo+1
*rak:          icount = 0
*rak:          do 00300 id = 1,2
*rak:            if (idatom(id).gt.0) then
*rak:              do 00400 dxyz = 1,3
*rak:                do 00500 ibf = ibflo,ibfhi
*rak:                  do 00600 jbf = jbflo, jbfhi
*rak:                    icount = icount + 1
*rak:                    index = isym2(ibf,jbf)
*rak:                    x2e2d(index,dxyz,idatom(id)) = buf(icount)
*rak:00600             continue
*rak:00500           continue
*rak:00400         continue
*rak:            else
*rak:              icount = icount + 3*num_ibf*num_jbf
*rak:            endif
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:      end
*rak:      subroutine build_2edfd(x2e2d_fd,nat,size2e,nbf,
*rak:     &    scr,nscr,buf,nbuf,xyz,delta,basis,geom,nshells)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:c
*rak:      integer nat,size2e,nbf, nscr, nbuf, basis, geom, nshells
*rak:      double precision x2e2d_fd(size2e,3,nat)
*rak:      double precision scr(nscr)
*rak:      double precision buf(nbuf,2)
*rak:      double precision xyz(3,nat)
*rak:      double precision delta, scale
*rak:      integer ish, ibf, ibflo, ibfhi
*rak:      integer jsh, jbf, jbflo, jbfhi
*rak:      integer nint, icount, index, datom, dxyz
*rak:      integer int_nint
*rak:      external int_nint
*rak:      integer min, max
*rak:      intrinsic min,max
*rak:      integer i,j,isym2
*rak:      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
*rak:c
*rak:c
*rak:c.. store original coords in xyz
*rak:      call dcopy((3*nat),coords(1,1,geom),1,xyz,1)
*rak:      scale = 1.0d00/(2.0d00*delta)
*rak:c
*rak:      do 00100 ish = 1,nshells
*rak:        do 00200 jsh = 1,ish
*rak:          nint = int_nint(basis,ish,basis,jsh,basis,0,basis,0)
*rak:          do 00300 datom = 1,nat
*rak:            do 00400 dxyz = 1,3
*rak:c... copy original coords in to coords array
*rak:              call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:c... add delta to coordinate
*rak:              coords(dxyz,datom,geom) = coords(dxyz,datom,geom) + delta
*rak:              call int_2e2c(basis,ish,basis,jsh,nscr,scr,nint,buf(1,1))
*rak:c... copy original coords in to coords array
*rak:              call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:c... subtract delta to coordinate
*rak:              coords(dxyz,datom,geom) = coords(dxyz,datom,geom) - delta
*rak:              call int_2e2c(basis,ish,basis,jsh,nscr,scr,nint,buf(1,2))
*rak:              call daxpy(nint,-1.0d00,buf(1,2),1,buf(1,1),1)
*rak:              call dscal(nint,scale,buf(1,1),1)
*rak:              icount = 0
*rak:              if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &            call errquit('cn2bfr error',ish)
*rak:              if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &            call errquit('cn2bfr error',ish)
*rak:              do 00500 ibf = ibflo, ibfhi
*rak:                do 00600 jbf = jbflo, jbfhi
*rak:                  icount = icount+1
*rak:                  index  = isym2(ibf,jbf)
*rak:                  x2e2d_fd(index,dxyz,datom) = buf(icount,1)
*rak:00600           continue
*rak:00500         continue
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:c      
*rak:c... copy original coords in to coords array
*rak:      call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
*rak:      end
*rak:      subroutine print_2d3c(x2e3d,diff,nbf,nat)
*rak:      implicit none
*rak:      integer nbf,nat
*rak:      double precision x2e3d(nbf,nbf,nbf,3,nat)
*rak:      double precision diff(nbf,nbf,nbf,3,nat)
*rak:c
*rak:      double precision thresh
*rak:      integer datom,dxyz, ibf,jbf,kbf
*rak:      logical printit
*rak:c
*rak:      thresh = 1.0d-08
*rak:c
*rak:      do 00100 datom = 1, nat
*rak:        do 00200 dxyz = 1,3
*rak:          write(6,*)
*rak:          write(6,*)' atom = ',datom,' xyz = ',dxyz
*rak:          write(6,*)'--------------------------------------------------'
*rak:          do 00300 ibf = 1,nbf
*rak:            do 00400 jbf = 1,nbf
*rak:              do 00500 kbf = 1,nbf
*rak:                printit = abs(x2e3d(ibf,jbf,kbf,dxyz,datom)).gt.thresh
*rak:                if (printit) then
*rak:                  write(6,10000)ibf,jbf,kbf,datom,dxyz,
*rak:     &                x2e3d(ibf,jbf,kbf,dxyz,datom),
*rak:     &                diff(ibf,jbf,kbf,dxyz,datom)
*rak:                endif
*rak:00500         continue
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:10000 format('d2e3c <',i3,',',i3,',',i3,'|',i3,'|',i3,' > val =',
*rak:     &    d20.10,' diff =',d20.10)
*rak:c
*rak:      end
*rak:      subroutine print_2d2c(x2e2d,x2e2d_fd,diff,size2e,nat)
*rak:      implicit none
*rak:      integer size2e,nat
*rak:      double precision x2e2d(size2e,3,nat)
*rak:      double precision x2e2d_fd(size2e,3,nat)
*rak:      double precision diff(size2e,3,nat)
*rak:c
*rak:      double precision thresh
*rak:      integer datom,dxyz, index
*rak:      logical printit
*rak:c
*rak:      thresh = 1.0d-08
*rak:c
*rak:      do 00100 datom = 1, nat
*rak:        do 00200 dxyz = 1,3
*rak:          write(6,*)
*rak:          write(6,*)' atom = ',datom,' xyz = ',dxyz
*rak:          write(6,*)'--------------------------------------------------'
*rak:          do 00300 index = 1,size2e
*rak:            printit = abs(x2e2d(index,dxyz,datom)).gt.thresh
*rak:            printit = printit .or.
*rak:     &          abs(x2e2d_fd(index,dxyz,datom)).gt.thresh
*rak:            if (printit) then
*rak:              write(6,10000)index,datom,dxyz,
*rak:     &            x2e2d(index,dxyz,datom),
*rak:     &            x2e2d_fd(index,dxyz,datom),
*rak:     &            diff(index,dxyz,datom)
*rak:            endif
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:10000 format('d2e2c <',i5,'|',i3,'|',i3,' > val =',d20.10,
*rak:     &    ' fd =',d20.10,' diff =',d20.10)
*rak:c
*rak:      end
*rak:      logical function gen_apphess_at(rtdb,delta)
*rak:      implicit none
*rak:#include "mafdecls.fh"
*rak:#include "global.fh"
*rak:#include "geom.fh"
*rak:#include "bas.fh"
*rak:#include "rtdb.fh"
*rak:#include "context.fh"
*rak:#include "tcgmsg.fh"
*rak:#include "msgtypesf.h"
*rak:#include "msgids.fh"
*rak:c::functions
*rak:c::passed
*rak:      integer rtdb
*rak:      double precision delta
*rak:c::local
*rak:      integer whoiam, master
*rak:      integer geom
*rak:      integer nat, rank_hess
*rak:      integer basis
*rak:      integer h_grad0, k_grad0 ! handle/index for central gradient
*rak:      integer h_gradp, k_gradp ! handle/index for delta gradient
*rak:      integer h_hess, k_hess ! handle/index for hessian matrix
*rak:      integer iatom, ixyz
*rak:      logical ignore_status
*rak:      double precision xyz(3), chg
*rak:      character*5 tmpname
*rak:      character*16 tag_name
*rak:      character*40 new_geom_name
*rak:      character*255 mo_basis, fd_basis, mo_basis_ref
*rak:c
*rak:      whoiam = ga_nodeid()
*rak:      master = 0
*rak:      ignore_status = rtdb_parallel(.false.)
*rak:      if (whoiam.eq.master) then
*rak:c
*rak:c.. create/load reference geometry
*rak:        if (.not.geom_create(geom,'geometry')) call errquit
*rak:     &      ('gen_apphess_at:geom_create failed?',1)
*rak:        if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
*rak:     &      ('gen_apphess_at:geom_rtdb_load failed?',2)
*rak:c..  get the number of atoms
*rak:        if (.not. geom_ncent(geom,nat)) call errquit
*rak:     &      ('gen_apphess_at:geom_ncent failed?',3)
*rak:c
*rak:c.. copy reference geometry and store on rtdb
*rak:        if (.not.
*rak:     &      geom_rtdb_store(rtdb,geom,'reference'))
*rak:     &      call errquit('gen_apphess_at: geom_rtdb_store failed',5)
*rak:c
*rak:c.. determine what fd basis name is
*rak:c
*rak:        if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
*rak:     &      mo_basis = 'ao basis'
*rak:        mo_basis_ref = mo_basis
*rak:        if (.not.context_rtdb_match(rtdb,'fd basis',fd_basis))
*rak:     &      fd_basis = 'fd basis'
*rak:c
*rak:c.. try to load "fd basis"
*rak:c.. create basis handle
*rak:        if(.not.bas_create(basis,fd_basis)) call errquit
*rak:     &      ('gen_apphess_at: bas_create failed?',7)
*rak:        if(bas_rtdb_load(rtdb,geom,basis,fd_basis)) then
*rak:          continue
*rak:        else
*rak:c...  if fd basis not defined use current mo_basis
*rak:          fd_basis = mo_basis
*rak:        endif
*rak:c
*rak:c set "ao basis" to point to "fd basis"
*rak:        if (.not.rtdb_cput(rtdb,'ao basis',1,fd_basis)) call errquit
*rak:     &      ('gen_apphess_at: rtdb_cpu failed?',8)
*rak:c
*rak:        if (.not.(bas_destroy(basis))) call errquit
*rak:     &      ('gen_apphess_at: bas_destroy failed?',9)
*rak:        if (.not.(geom_destroy(geom))) call errquit
*rak:     &      ('gen_apphess_at: geom_destroy failed?',9)
*rak:c set "geometry" to "reference"
*rak:        new_geom_name='reference'
*rak:        if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
*rak:     &      call errquit
*rak:     &        ('gen_apphess_at: rtdb_cput failed ?',911)
*rak:      endif
*rak:c all nodes work here 
*rak:c
*rak:c..      broadcast number of atoms
*rak:c
*rak:      call ga_brdcst(Msg_RAK+MSGINT, nat, mitob(1), 0)
*rak:c
*rak:      call ga_sync()
*rak:c
*rak:      ignore_status = rtdb_parallel(.true.)
*rak:c      
*rak:      call scf_gradient(rtdb)
*rak:c
*rak:      ignore_status = rtdb_parallel(.false.)
*rak:c      
*rak:      call ga_sync()
*rak:      if (whoiam.eq.master) then
*rak:        rank_hess = 3*nat
*rak:        if (.not. MA_Push_Get(MT_DBL,rank_hess,
*rak:     &      'cent nuclear gradient vector',
*rak:     &      h_grad0,k_grad0)) call errquit
*rak:     &      ('gen_apphess_at: ma_push_get failed?',10)
*rak:        if (.not. MA_Push_Get(MT_DBL,rank_hess,
*rak:     &      'delta nuclear gradient vector',
*rak:     &      h_gradp,k_gradp)) call errquit
*rak:     &      ('gen_apphess_at: ma_push_get failed?',11)
*rak:        if (.not. MA_Push_Get(MT_DBL,(rank_hess*rank_hess),
*rak:     &      'nuclear hessian matrix',
*rak:     &      h_hess,k_hess)) call errquit
*rak:     &      ('gen_apphess_at: ma_push_get failed?',12)
*rak:c
*rak:        call dfill((rank_hess*rank_hess),0.0d00,
*rak:     &      dbl_mb(k_hess),1)
*rak:c
*rak:c.. get central gradient
*rak:        if(.not. rtdb_get(rtdb,'scf:gradients',MT_DBL,
*rak:     &      rank_hess,dbl_mb(k_grad0)))
*rak:     &      call errquit('gen_apphess_at: rtdb_get failed',13)
*rak:c
*rak:      endif
*rak:c
*rak:c... all nodes doing work here
*rak:      call ga_sync()
*rak:      call util_flush(6)
*rak:      do 00100 iatom = 1,nat
*rak:        do 00200 ixyz = 1,3
*rak:          if(whoiam.eq.master) then
*rak:c... master node only
*rak:            if (.not.geom_create(geom,'reference')) call errquit
*rak:     &          ('gen_apphess_at:geom_create failed?',14)
*rak:            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at:geom_rtdb_load failed?',15)
*rak:            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at:geom_cent_get failed?',16)
*rak:            xyz(ixyz) = xyz(ixyz) + delta
*rak:            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at:geom_cent_get failed?',17)
*rak:c.. copy modified geometry and store on rtdb
*rak:            write(tmpname,'(i3,1x,i1)')iatom,ixyz
*rak:            new_geom_name = 'reference '//tmpname
*rak:            if (.not.
*rak:     &          geom_rtdb_store(rtdb,geom,new_geom_name))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at: geom_rtdb_store failed',19)
*rak:            if (.not.geom_destroy(geom))
*rak:     &          call errquit('gen_apphess_at: geom_destroy failed?',20)
*rak:c set "geometry" to "reference atom/xyz"
*rak:            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at: rtdb_cput failed ?',911)
*rak:          endif
*rak:          call ga_sync()
*rak:          ignore_status = rtdb_parallel(.true.)
*rak:          call scf_gradient(rtdb)
*rak:          ignore_status = rtdb_parallel(.false.)
*rak:          call ga_sync()
*rak:          if (whoiam.eq.master) then
*rak:            if(.not. rtdb_get(rtdb,'scf:gradients',MT_DBL,
*rak:     &          rank_hess,dbl_mb(k_gradp)))
*rak:     &          call errquit('gen_apphess_at: rtdb_get failed',21)
*rak:            call fd_upd_hess(dbl_mb(k_hess),dbl_mb(k_grad0),
*rak:     &          dbl_mb(k_gradp),1.0d00,delta,nat,iatom,ixyz)
*rak:          endif
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:      if (whoiam.eq.master) then
*rak:        write(6,*)' finite difference hessian delta = ',delta
*rak:        call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
*rak:     &      rank_hess,rank_hess,1)
*rak:        call wrt_fd_from_sq(dbl_mb(k_hess),rank_hess,'hess_ascii')
*rak:        write(6,*)' triangle hessian written to hess_ascii'
*rak:        
*rak:c.. redefine ao basis to it's original setting
*rak:        if (.not.rtdb_cput(rtdb,'ao basis',1,mo_basis_ref))
*rak:     &      call errquit
*rak:     &      ('gen_apphess_at: rtdb_cpu failed?',22)
*rak:c set "geometry" to "geometry"
*rak:            if(.not.rtdb_cput(rtdb,'geometry',1,'geometry'))
*rak:     &          call errquit
*rak:     &          ('gen_apphess_at: rtdb_cput failed ?',911)
*rak:c
*rak:        if(.not.ma_pop_stack(h_hess)) call errquit
*rak:     &      ('gen_apphess_at: ma_pop_stack(h_hess) failed?',27)
*rak:        if(.not.ma_pop_stack(h_gradp)) call errquit
*rak:     &      ('gen_apphess_at: ma_pop_stack(h_gradp) failed?',27)
*rak:        if(.not.ma_pop_stack(h_grad0)) call errquit
*rak:     &      ('gen_apphess_at: ma_pop_stack(h_grad0) failed?',27)
*rak:
*rak:      endif
*rak:      call ga_sync()
*rak:      ignore_status = rtdb_parallel(.true.)
*rak:      gen_apphess_at = .true.
*rak:      end
*rak:      subroutine wrt_fd_from_sq(hess,rank_hess,filename)
*rak:      implicit none
*rak:c
*rak:      integer rank_hess
*rak:      double precision hess(rank_hess,rank_hess)
*rak:      character*(*) filename
*rak:c
*rak:      integer i, j, lu
*rak:c
*rak:      lu = 69
*rak:      open(unit=lu,file=filename,
*rak:     &    form='formatted',
*rak:     &    access='sequential',
*rak:     &    status='unknown')
*rak:c
*rak:      do 00100 i = 1,rank_hess
*rak:        do 00200 j = 1,i
*rak:          write(lu,10000)hess(i,j)
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:10000 format(1x,1pd20.10)
*rak:c
*rak:      close(unit=lu,status='keep')
*rak:c
*rak:      end
*rak:      subroutine fd_upd_hess(hess,gradm,gradp,s_delta,delta,nat,
*rak:     &    iatom_t,ixyz_t)
*rak:      implicit none
*rak:c::passed
*rak:      integer nat
*rak:      integer iatom_t, ixyz_t
*rak:      double precision hess(3,nat,3,nat)
*rak:      double precision gradm(3,nat)
*rak:      double precision gradp(3,nat)
*rak:      double precision delta, s_delta
*rak:c::local
*rak:      integer iatom,ixyz
*rak:      double precision rdelta, value
*rak:c
*rak:c 
*rak:c finite difference  [g(x+delta) - g(x-delta)]/(s_delta*delta) (s_delta = 2.0)
*rak:c central difference [g(x+delta) - g(x)]/(s_delta*delta) (s_delta = 1.0)
*rak:c 
*rak:c
*rak:      rdelta = 1.0d00/(s_delta*delta)
*rak:      do 00100 iatom = 1,nat
*rak:        do 00200 ixyz = 1,3
*rak:          value = rdelta*(gradp(ixyz,iatom)-gradm(ixyz,iatom))
*rak:          hess(ixyz_t,iatom_t,ixyz,iatom) = value
*rak:          hess(ixyz,iatom,ixyz_t,iatom_t) = value
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:      end
*rak:      logical function gen_apphess_opt(rtdb,delta)
*rak:      implicit none
*rak:c
*rak:c function to generate approximate hessian
*rak:c
*rak:#include "bas.fh"
*rak:#include "mafdecls.fh"
*rak:#include "global.fh"
*rak:#include "geom.fh"
*rak:#include "rtdb.fh"
*rak:#include "context.fh"
*rak:      integer rtdb
*rak:      double precision delta
*rak:c
*rak:      integer whoiam, master
*rak:      logical ignore_status
*rak:      integer geom,nat,basis
*rak:      character*40 mo_basis,mo_basis_ref,fd_basis,new_geom_name
*rak:      character*40 movecs_in_ref
*rak:c
*rak:      gen_apphess_opt=.false.
*rak:      whoiam = ga_nodeid()
*rak:      master = 0
*rak:      ignore_status = rtdb_parallel(.false.)
*rak:      if (whoiam.eq.master) then
*rak:c
*rak:c.. create/load reference geometry
*rak:        if (.not.geom_create(geom,'geometry')) call errquit
*rak:     &      ('gen_apphess_at:geom_create failed?',1)
*rak:        if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
*rak:     &      ('gen_apphess_at:geom_rtdb_load failed?',2)
*rak:c..  get the number of atoms
*rak:        if (.not. geom_ncent(geom,nat)) call errquit
*rak:     &      ('gen_apphess_at:geom_ncent failed?',3)
*rak:c
*rak:c.. copy reference geometry and store on rtdb
*rak:        if (.not.
*rak:     &      geom_rtdb_store(rtdb,geom,'reference'))
*rak:     &      call errquit('gen_apphess_at: geom_rtdb_store failed',5)
*rak:c
*rak:c.. determine what fd basis name is
*rak:c
*rak:        if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
*rak:     &      mo_basis = 'ao basis'
*rak:        mo_basis_ref = mo_basis
*rak:        if (.not.context_rtdb_match(rtdb,'fd basis',fd_basis))
*rak:     &      fd_basis = 'fd basis'
*rak:c
*rak:c.. try to load "fd basis"
*rak:c.. create basis handle
*rak:        if(.not.bas_create(basis,fd_basis)) call errquit
*rak:     &      ('gen_apphess_at: bas_create failed?',7)
*rak:        if(bas_rtdb_load(rtdb,geom,basis,fd_basis)) then
*rak:          continue
*rak:        else
*rak:c...  if fd basis not defined use current mo_basis
*rak:          fd_basis = mo_basis
*rak:          write(6,*) ' fd basis is the ao basis This will be expensive'
*rak:          call util_flush(6)
*rak:        endif
*rak:c
*rak:c.. set "ao basis" to point to "fd basis"
*rak:        if (.not.rtdb_cput(rtdb,'ao basis',1,fd_basis)) call errquit
*rak:     &      ('gen_apphess_at: rtdb_cpu failed?',8)
*rak:c... destroty basis and geom handles
*rak:        if (.not.(bas_destroy(basis))) call errquit
*rak:     &      ('gen_apphess_at: bas_destroy failed?',9)
*rak:        if (.not.(geom_destroy(geom))) call errquit
*rak:     &      ('gen_apphess_at: geom_destroy failed?',9)
*rak:c set "geometry" to "reference"
*rak:        new_geom_name='reference'
*rak:        if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
*rak:     &      call errquit
*rak:     &        ('gen_apphess_at: rtdb_cput failed ?',911)
*rak:      endif
*rak:c... store movecs initial state
*rak:      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in_ref))
*rak:     $     movecs_in_ref = 'atomic'
*rak:c... set initial guess for fd basis to atomic
*rak:      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, 'atomic'))
*rak:     &     call errquit ('gen_apphess_opt: rtdb_cput failed?',911)
*rak:c
*rak:      ignore_status = rtdb_parallel(.true.)
*rak:      call gen_apphess_at(rtdb,delta)
*rak:      call scf_opt(rtdb)
*rak:c
*rak:c.. reset stepper so it knows to restart the calculation
*rak:      call stpr_walk_reset
*rak:      ignore_status = rtdb_parallel(.false.)
*rak:c.. set "ao basis" to point to %val<mo_basis_ref
*rak:      if (.not.rtdb_cput(rtdb,'ao basis',1,mo_basis_ref))
*rak:     &    call errquit ('gen_apphess_opt: rtdb_cput failed?',8)
*rak:c
*rak:c... reset input vectors option here to original guess if appropriate
*rak:      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, movecs_in_ref))
*rak:     &     call errquit ('gen_apphess_opt: rtdb_cput failed?',911)
*rak:c
*rak:      call scf_opt(rtdb)
*rak:c
*rak:      gen_apphess_opt = .true.
*rak:      end
*rak:      subroutine size_integ(rtdb,geom_name, basis_name)
*rak:      implicit none
*rak:#include "geom.fh"
*rak:#include "bas.fh"
*rak:#include "mafdecls.fh"
*rak:c
*rak:      integer rtdb
*rak:      integer basis, geom
*rak:      character*(*) geom_name, basis_name
*rak:c
*rak:      integer max1e,max2e,scr1,scr2,maxuse,scruse
*rak:      double precision xscruse
*rak:c
*rak:      logical int_normalize
*rak:      double precision ddot
*rak:      external int_normalize
*rak:      external ddot
*rak:c
*rak:      call hf_print_set(0)
*rak:c
*rak:      write(6,*)' rtdb      :',rtdb
*rak:      write(6,*)' geom_name : <',geom_name,'>'
*rak:      write(6,*)' basis_name: <',basis_name,'>'
*rak:c
*rak:      if (.not. geom_create(geom, geom_name))
*rak:     $     call errquit('size_integ: geom_create failed?', 0)
*rak:      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
*rak:     $     call errquit('size_integ: geom_load failed', 0)
*rak:      if (.not. geom_print(geom))
*rak:     $     call errquit('size_integ: geom_print failed', 0)
*rak:      if (.not. bas_create(basis, basis_name))
*rak:     $     call errquit('size_integ: basis create failed', 0)
*rak:      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
*rak:     $      call errquit('size_integ: basis load failed', 0)
*rak:      if (.not. bas_print(basis))
*rak:     $      call errquit('size_integ: basis print failed', 0)
*rak:      if (.not. gbs_map_print(basis))
*rak:     $     call errquit('size_integ: gbs_map_print failed', 0)
*rak:c
*rak:c normalize basis set
*rak:c
*rak:      if (.not.int_normalize(basis))
*rak:     &       call errquit('size_integ: basis norm. failed', 0)
*rak:      write(6,*)' after normalization'
*rak:      if (.not. bas_print(basis))
*rak:     $      call errquit('size_integ: basis print failed', 0)
*rak:c
*rak:c      write(6,*)' with integrals  printed '
*rak:c      call int_chk_sum(basis,.true.)
*rak:      call int_init(rtdb,1,basis)
*rak:      call int_terminate()
*rak:      call intd_init(rtdb,1,basis)
*rak:      call int_mem(max1e,max2e,scr1,scr2)
*rak:      maxuse=max(max1e,max2e)
*rak:      scruse=max(scr1,scr2)
*rak:      scruse=scruse+maxuse
*rak:      scruse=scruse*8
*rak:      xscruse=dble(scruse)*1.0d-06
*rak:      write(6,'(1x,a,f10.3,a)')' int deriv memory ',xscruse,' Mbytes'
*rak:      call int_terminate()
*rak:      if (.not.(
*rak:     &       (bas_destroy(basis)).and.(geom_destroy(geom))
*rak:     &       ))
*rak:     &       call errquit('rak:error destroying gbs pair',0)
*rak:c
*rak:      end
*rak:      subroutine derivtest(rtdb,geom_name, basis_name)
*rak:      implicit none
*rak:#include "geom.fh"
*rak:#include "bas.fh"
*rak:#include "mafdecls.fh"
*rak:c
*rak:      integer rtdb
*rak:      integer basis, geom
*rak:      integer nscr, hscr, iscr
*rak:      integer nbf, nshells, nat, maxg, size2e,size2ed
*rak:      integer h2e, i2e, h2efd, i2efd, h2ep, i2ep, h2em, i2em 
*rak:      integer h2ediff, i2ediff
*rak:      integer hbuf2e, ibuf2e
*rak:      integer buf1esz, hbuf1e, ibuf1e
*rak:      integer sa_size, hsa, isa, hsafd, isafd
*rak:      integer hta, ita, htafd, itafd
*rak:      integer threen
*rak:      integer hxyz, ixyz
*rak:      character*(*) geom_name, basis_name
*rak:      logical int_normalize
*rak:      double precision ddot
*rak:      external int_normalize
*rak:      external ddot
*rak:c
*rak:      call hf_print_set(0)
*rak:c
*rak:      write(6,*)' rtdb      :',rtdb
*rak:      write(6,*)' geom_name : <',geom_name,'>'
*rak:      write(6,*)' basis_name: <',basis_name,'>'
*rak:c
*rak:      if (.not. geom_create(geom, geom_name))
*rak:     $     call errquit('derivtest: geom_create failed?', 0)
*rak:      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
*rak:     $     call errquit('derivtest: geom_load failed', 0)
*rak:      if (.not. geom_print(geom))
*rak:     $     call errquit('derivtest: geom_print failed', 0)
*rak:      if (.not. bas_create(basis, basis_name))
*rak:     $     call errquit('derivtest: basis create failed', 0)
*rak:      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
*rak:     $      call errquit('derivtest: basis load failed', 0)
*rak:      if (.not. bas_print(basis))
*rak:     $      call errquit('derivtest: basis print failed', 0)
*rak:      if (.not. gbs_map_print(basis))
*rak:     $     call errquit('derivtest: gbs_map_print failed', 0)
*rak:c
*rak:c normalize basis set
*rak:c
*rak:      if (.not.int_normalize(basis))
*rak:     &       call errquit('derivtest: basis norm. failed', 0)
*rak:      write(6,*)' after normalization'
*rak:      call intd_init(rtdb,1,basis)
*rak:      if (.not. bas_print(basis))
*rak:     $      call errquit('derivtest: basis print failed', 0)
*rak:c
*rak:c      write(6,*)' with integrals  printed '
*rak:c      call int_chk_sum(basis,.true.)
*rak:      write(6,*)' without integrals  printed '
*rak:      call int_chk_sum(basis,.false.)
*rak:
*rak:c
*rak:      write(6,*)' derivative test '
*rak:c
*rak:*      call hf_print_set(1)
*rak:*c
*rak:*      nscr = 1 000 000
*rak:*      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
*rak:*     &       call errquit('raktest: allocation scr failed',nscr)
*rak:*      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:*      call onedint_fd(basis,geom,dbl_mb(iscr),nscr)
*rak:*      if (.not.MA_Pop_Stack(hscr))
*rak:*     &       write(6,*)'raktest: pop on hscr failed'
*rak:c
*rak:      
*rak:      if(.not.bas_numbf(basis,nbf))
*rak:     &       stop 'raktest:bas_numbf failed'
*rak:      if(.not.bas_numcont(basis,nshells))
*rak:     &       stop 'raktest:bas_numcont failed'
*rak:      if(.not.geom_ncent(geom,nat))
*rak:     &       stop 'raktest:geom_ncent failed'
*rak:c
*rak:      sa_size = (nbf*nbf)*(3*nat)
*rak:c..... space for sa
*rak:      if (.not.MA_Push_Get(MT_dbl,sa_size,'sa buffer',hsa,isa))
*rak:     &       call errquit('raktest: allocation for sa failed')
*rak:c..... space for safd
*rak:      if (.not.MA_Push_Get(MT_dbl,sa_size,'safd buffer',hsafd,isafd))
*rak:     &       call errquit('raktest: allocation for safd failed')
*rak:      threen = 3*nat
*rak:      if (.not.MA_Push_Get(MT_dbl,2*threen,'xyzxyz buffer',hxyz,ixyz))
*rak:     &       call errquit('raktest: allocation for safd failed')
*rak:      buf1esz = nbf*nbf*3
*rak:      if (.not.MA_Push_Get(MT_dbl,buf1esz,'safd buffer',hbuf1e,ibuf1e))
*rak:     &       call errquit('raktest: allocation for buf1e failed')
*rak:c..... scratch space
*rak:      nscr = 100 000
*rak:      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
*rak:     &       call errquit('raktest: allocation scr failed',nscr)
*rak:c..... space for Ta&Tafd so to sum it to Va&Fafd
*rak:      if (.not.MA_Push_Get(MT_Dbl,sa_size,'Ta',hta,ita))
*rak:     &       call errquit('raktest: allocation scr failed',sa_size)
*rak:      if (.not.MA_Push_Get(MT_Dbl,sa_size,'Ta fd',htafd,itafd))
*rak:     &       call errquit('raktest: allocation scr failed',sa_size)
*rak:      call dfill(sa_size,0.0,dbl_mb(ita),1)
*rak:      call dfill(sa_size,0.0,dbl_mb(itafd),1)
*rak:      call dfill(sa_size,0.0,dbl_mb(isa),1)
*rak:      call dfill(sa_size,0.0,dbl_mb(isafd),1)
*rak:      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
*rak:      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call buildSa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,.true.)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:      call buildSanew(dbl_mb(isafd),nat,nshells,nbf,basis,
*rak:     &       dbl_mb(iscr),nscr,.true.)
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
*rak:      write(6,*)' Sanew difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:      call buildSafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
*rak:     &       dbl_mb(ibuf1e+2*nbf*nbf),
*rak:     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
*rak:      write(6,*)' Sa difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:c
*rak:c.... do T
*rak:      call dfill(sa_size,0.0,dbl_mb(isa),1)
*rak:      call dfill(sa_size,0.0,dbl_mb(isafd),1)
*rak:      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
*rak:      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call buildTa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,.true.)
*rak:      call dcopy(sa_size,dbl_mb(isa),1,dbl_mb(ita),1) 
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:      call buildTafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
*rak:     &       dbl_mb(ibuf1e+2*nbf*nbf),
*rak:     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
*rak:      call dcopy(sa_size,dbl_mb(isafd),1,dbl_mb(itafd),1) 
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
*rak:      write(6,*)' Ta difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:c.. do V
*rak:      call dfill(sa_size,0.0,dbl_mb(isa),1)
*rak:      call dfill(sa_size,0.0,dbl_mb(isafd),1)
*rak:      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
*rak:      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call buildVa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,.true.)
*rak:      call daxpy(sa_size,1.0d00,dbl_mb(isa),1,dbl_mb(ita),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:      call buildVafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
*rak:     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
*rak:     &       dbl_mb(ibuf1e+2*nbf*nbf),
*rak:     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
*rak:      call daxpy(sa_size,1.0d00,dbl_mb(isafd),1,dbl_mb(itafd),1)
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
*rak:      call printdmat('Va difference matrix',dbl_mb(isafd),nat,nbf)
*rak:      write(6,*)' Va difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:c... print H derivs
*rak:      call printdmat(' Ha ',dbl_mb(ita),nat,nbf)
*rak:      call printdmat(' Ha <<<FD>>> ',dbl_mb(itafd),nat,nbf)
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(ita),1,dbl_mb(itafd),1)
*rak:      write(6,*)' Ha difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(itafd),1,dbl_mb(itafd),1)
*rak:c
*rak:      call buildHa(dbl_mb(itafd),nat,nshells,nbf,basis,
*rak:     &       dbl_mb(iscr),nscr,.true.)
*rak:      call printdmat(' Ha New',dbl_mb(itafd),nat,nbf)
*rak:      call daxpy(sa_size,-1.0d00,dbl_mb(ita),1,dbl_mb(itafd),1)
*rak:      write(6,*)' Ha new difference norm =',
*rak:     &       ddot(sa_size,dbl_mb(itafd),1,dbl_mb(itafd),1)
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:c
*rak:      if (.not.MA_Pop_Stack(htafd))
*rak:     &       call errquit('raktest: ma_pop_stack htafd failed',911)
*rak:      if (.not.MA_Pop_Stack(hta))
*rak:     &       call errquit('raktest: ma_pop_stack hta failed',911)
*rak:      if (.not.MA_Pop_Stack(hscr))
*rak:     &       call errquit('raktest: ma_pop_stack hscr failed',911)
*rak:      if (.not.MA_Pop_Stack(hbuf1e))
*rak:     &       call errquit('raktest: ma_pop_stack hbu1e failed',911)
*rak:      if (.not.MA_Pop_Stack(hxyz))
*rak:     &       call errquit('raktest: ma_pop_stack hxyz failed',911)
*rak:      if (.not.MA_Pop_Stack(hsafd))
*rak:     &       call errquit('raktest: ma_pop_stack hsafd failed',911)
*rak:      if (.not.MA_Pop_Stack(hsa))
*rak:     &       call errquit('raktest: ma_pop_stack hsa failed',911)
*rak:c
*rak:c
*rak:c
*rak:c... do two electron integral stuff
*rak:c
*rak:      call int_mem_2e4c(maxg,nscr)
*rak:      write(6,*)' nscr    => ',nscr
*rak:      write(6,*)' maxg    => ',maxg
*rak:      nscr = nscr*12    ! 8x scratch for derivatives
*rak:      maxg = 12*maxg    ! quadruple buffer size for derivatives
*rak:      size2e = nbf*(nbf+1)/2
*rak:      size2e = size2e*(size2e+1)/2
*rak:      size2ed= size2e*3*nat
*rak:      write(6,*)' nscr    => ',nscr
*rak:      write(6,*)' maxg    => ',maxg
*rak:      write(6,*)' nbf     => ',nbf
*rak:      write(6,*)' nat     => ',nat
*rak:      write(6,*)' size2e  => ',size2e
*rak:      write(6,*)' size2ed => ',size2ed
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e d',h2ediff,i2ediff))
*rak:     &       call errquit('raktest: allocation h2ediff failed',size2ed)
*rak:      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e d',h2e,i2e))
*rak:     &       call errquit('raktest: allocation h2e failed',size2ed)
*rak:      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e fd',h2efd,i2efd))
*rak:     &       call errquit('raktest: allocation h2efd failed',size2ed)
*rak:      if (.not.MA_Push_Get(MT_Dbl,size2e,'2e +buf',h2ep,i2ep))
*rak:     &       call errquit('raktest: allocation h2ep failed',size2e)
*rak:      if (.not.MA_Push_Get(MT_Dbl,size2e,'2e -buf',h2em,i2em))
*rak:     &       call errquit('raktest: allocation h2em failed',size2e)
*rak:      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr 2ed',hscr,iscr))
*rak:     &       call errquit('raktest: allocation hscr failed',nscr)
*rak:      if (.not.MA_Push_Get(MT_Dbl,maxg,'buf2e 2ed',hbuf2e,ibuf2e))
*rak:     &       call errquit('raktest: allocation hbuf2e failed',maxg)
*rak:      threen = 3*nat
*rak:      if (.not.MA_Push_Get(MT_dbl,2*threen,'xyzxyz buffer',hxyz,ixyz))
*rak:     &       call errquit('raktest: allocation for safd failed')
*rak:      call dfill(size2ed,0.0,dbl_mb(i2e),1)
*rak:      call dfill(size2ed,0.0,dbl_mb(i2efd),1)
*rak:      call dfill(size2e,0.0,dbl_mb(i2ep),1)
*rak:      call dfill(size2e,0.0,dbl_mb(i2em),1)
*rak:      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
*rak:      
*rak:      call build_d2e(geom,basis,nbf,nshells,nat,size2e,dbl_mb(i2e),
*rak:     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
*rak:*      call print2de(nat,size2e,dbl_mb(i2e),1.0d-05,nbf,'anal:d2e')
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
*rak:      call build_d2efd(geom,basis,nbf,nshells,nat,size2e,
*rak:     &       dbl_mb(i2efd),
*rak:     &       dbl_mb(i2ep),dbl_mb(i2em),dbl_mb(ixyz),
*rak:     &       dbl_mb(ixyz+threen),
*rak:     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
*rak:*      call print2de(nat,size2e,dbl_mb(i2efd),1.0d-05,nbf,'fd:d2e')
*rak:      call dcopy(size2ed,dbl_mb(i2efd),1,dbl_mb(i2ediff),1)
*rak:      call daxpy(size2ed,-1.0d00,dbl_mb(i2e),1,dbl_mb(i2ediff),1)
*rak:      write(6,*)' 2ed difference norm = ',
*rak:     &       ddot(size2ed,dbl_mb(i2ediff),1,dbl_mb(i2ediff),1)
*rak:*      call print2de_diff(nat,size2e,dbl_mb(i2ediff),
*rak:*     &       dbl_mb(i2e),dbl_mb(i2efd),1.0d-05,nbf)
*rak:      call dfill(size2ed,0.0,dbl_mb(i2efd),1)
*rak:      call dfill(nscr,0.0,dbl_mb(iscr),1)
*rak:      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
*rak:      call build_d2e_new(basis,nbf,nshells,nat,size2e,dbl_mb(i2efd),
*rak:     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
*rak:      call dcopy(size2ed,dbl_mb(i2efd),1,dbl_mb(i2ediff),1)
*rak:      call daxpy(size2ed,-1.0d00,dbl_mb(i2e),1,dbl_mb(i2ediff),1)
*rak:      write(6,*)' new 2ed difference norm = ',
*rak:     &       ddot(size2ed,dbl_mb(i2ediff),1,dbl_mb(i2ediff),1)
*rak:
*rak:      call print2de_diff(nat,size2e,dbl_mb(i2ediff),
*rak:     &       dbl_mb(i2e),dbl_mb(i2efd),1.0d-05,nbf)
*rak:
*rak:      if (.not.MA_Pop_Stack(hxyz))
*rak:     &       call errquit('raktest: ma_pop_stack hxyz failed',911)
*rak:      if (.not.MA_Pop_Stack(hbuf2e))
*rak:     &       call errquit('raktest: ma_pop_stack hbuf2e failed',911)
*rak:      if (.not.MA_Pop_Stack(hscr))
*rak:     &       call errquit('raktest: ma_pop_stack hscr failed',911)
*rak:      if (.not.MA_Pop_Stack(h2em))
*rak:     &       call errquit('raktest: ma_pop_stack h2em failed',911)
*rak:      if (.not.MA_Pop_Stack(h2ep))
*rak:     &       call errquit('raktest: ma_pop_stack h2ep failed',911)
*rak:      if (.not.MA_Pop_Stack(h2efd))
*rak:     &       call errquit('raktest: ma_pop_stack h2efd failed',911)
*rak:      if (.not.MA_Pop_Stack(h2e))
*rak:     &       call errquit('raktest: ma_pop_stack h2e failed',911)
*rak:      if (.not.MA_Pop_Stack(h2ediff))
*rak:     &       call errquit('raktest: ma_pop_stack h2ediff failed',911)
*rak:00911 call int_terminate()
*rak:      if (.not.(
*rak:     &       (bas_destroy(basis)).and.(geom_destroy(geom))
*rak:     &       ))
*rak:     &       call errquit('rak:error destroying gbs pair',0)
*rak:      call MA_summarize_allocated_blocks()
*rak:      end
*rak:c.......................................................................
*rak:      Subroutine onedint_fd(basisin,geom,W0,maxW0)
*rak:
*rak:      implicit none 
*rak:
*rak:      Logical canAB,spherical
*rak:      integer basisin, basis, geom
*rak:
*rak:c      double precision EPS
*rak:c      Parameter (EPS=1.D-15)
*rak:
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c----#include "dimensions.fh"
*rak:c----#include "atoms.fh"
*rak:c----#include "AO_basis_set.fh"
*rak:c----#include "timer.fh"
*rak:
*rak:c--> Scratch Space
*rak:
*rak:      integer ltime
*rak:      parameter (ltime= 69)
*rak:      integer maxW0, itime, i, nshells, i_top
*rak:      integer ucont, nshbfa, nshbfb, nint, ncenters
*rak:      integer iscr, nscr, i2, i1
*rak:      integer i_s0, i_s2, i_s1, i_sa
*rak:      integer i_t0, i_t2, i_t1, i_ta
*rak:      integer i_v0, i_v2, i_v1, i_va
*rak:      integer isha, ictra, igeoa, ipa, npa, ica, nca, la
*rak:      integer ishb, ictrb, igeob, ipb, npb, icb, ncb, lb
*rak:      double precision time(10)
*rak:      Double Precision W0(maxW0),Txyz(3,2)
*rak:      double precision delta
*rak:      double precision adifv, adifs, adift
*rak:c
*rak:c Compute analytic 1-electron integrals derivatives and compare to finite
*rak:c difference evaluation.
*rak:c
*rak:c******************************************************************************
*rak:
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:c Set timer.
*rak:
*rak:      itime = 1
*rak:
*rak:      do 10 i = 1,10
*rak:       time(i) = 0.D0
*rak:   10 continue
*rak:
*rak:      spherical = .FALSE.
*rak:
*rak:c Compute the 1-electron integrals (S, T, V).
*rak:
*rak:c      open(1,file='1e')
*rak:
*rak:      if (.not.bas_numcont(basisin,nshells))
*rak:     &       stop ' bas_numcont:error'
*rak:      
*rak:      ncenters = ncenter(geom)
*rak:      do 190 isha = 1,nshells
*rak:
*rak:       do 180 ishb = 1,isha
*rak:  
*rak:c         write(6,*)' coords at <isha:ishb> <',isha,':',ishb,'>'
*rak:c         if (.not. geom_print(geom)) stop 'onedint_fd:geom_print'
*rak:c Set the canonicalization switches.
*rak:
*rak:        if( isha.eq.ishb )then
*rak:         canAB = .false.
*rak:        else
*rak:         canAB = .FALSE.
*rak:        end if
*rak:
*rak:c Define the basis sets on each of the centers.
*rak:
*rak:        ictra = ibs_cn2ce(isha,basis)
*rak:        igeoa = ibs_geom(basis)
*rak:        ucont = ibs_cn2ucn(isha,basis)
*rak:        ipa   = infbs_cont(CONT_IEXP,ucont,basis)
*rak:        npa   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:        ica   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:        nca   = infbs_cont(CONT_NGEN,ucont,basis)
*rak:        la    = infbs_cont(CONT_TYPE,ucont,basis)
*rak:
*rak:        ictrb = ibs_cn2ce(ishb,basis)
*rak:        igeob = igeoa                       ! always same basis
*rak:        ucont = ibs_cn2ucn(ishb,basis)
*rak:        ipb   = infbs_cont(CONT_IEXP,ucont,basis)
*rak:        npb   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:        icb   = infbs_cont(CONT_ICFP,ucont,basis) 
*rak:        ncb   = infbs_cont(CONT_NGEN,ucont,basis) 
*rak:        lb    = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:
*rak:c Define the number of integrals in this block of shells.
*rak: 
*rak:        if( spherical )then
*rak:         nshbfa = 2*la+1
*rak:         nshbfb = 2*lb+1
*rak:        else
*rak:         nshbfa = ((la+1)*(la+2))/2
*rak:         nshbfb = ((lb+1)*(lb+2))/2
*rak:       end if
*rak: 
*rak:        if( canAB )then
*rak:         nint = (nshbfa*(nshbfa+1))/2
*rak:        else
*rak:         nint = nshbfa*nshbfb
*rak:        end if
*rak: 
*rak:c Assign pointers to scratch space.
*rak:
*rak:        i_S0  = 1
*rak:        i_S1  = i_S0 + nint*6
*rak:        i_S2  = i_S1 + nint
*rak:        i_top = i_S2 + nint - 1
*rak:
*rak:        i_T0  = i_top + 1
*rak:        i_T1  = i_T0 + nint*6
*rak:        i_T2  = i_T1 + nint
*rak:        i_top = i_T2 + nint - 1
*rak:
*rak:        i_V0   = i_top  + 1
*rak:        i_V1   = i_V0   + nint*(3*ncenters)
*rak:        i_V2   = i_V1   + nint
*rak:        i_top  = i_V2   + nint
*rak:
*rak:        if( i_top.gt.maxW0 )then
*rak:
*rak:         write(*,*) 'ONEDINT:  Insufficient scratch space.'
*rak:         write(*,*) '          needed    ',i_top
*rak:         write(*,*) '          allocated ',maxW0
*rak:
*rak:         write(*,*) 'S0:  ',i_S0
*rak:         write(*,*) 'S1:  ',i_S1
*rak:         write(*,*) 'S2:  ',i_S2
*rak:
*rak:         write(*,*) 'T0:  ',i_T0
*rak:         write(*,*) 'T1:  ',i_T1
*rak:         write(*,*) 'T2:  ',i_T2
*rak:
*rak:         write(*,*) 'V0:  ',i_V0
*rak:         write(*,*) 'V1:  ',i_V1
*rak:         write(*,*) 'V2:  ',i_V2
*rak:
*rak:         stop
*rak:
*rak:        end if
*rak:
*rak:        iscr = i_top + 1
*rak:        nscr = maxW0 - i_top
*rak:
*rak:c Analytic derivatives.
*rak:
*rak:        i2  = 1
*rak:        i1  = 1
*rak:
*rak:        i_Sa = i_S0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:        i_Ta = i_T0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:        i_Va = i_V0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:
*rak:                     itime = 1
*rak:
*rak:c        call hf1d(xyz(1,ictra),zprims(ipa),ccoef(ica),npa,nca,la,
*rak:c     &            xyz(1,ictrb),zprims(ipb),ccoef(icb),npb,ncb,lb,
*rak:c     &            xyz,zan,ncenters,ictra,ictrb,
*rak:c     &            W0(i_Sa),W0(i_Ta),W0(i_Va),nint,
*rak:c     &            .TRUE.,.TRUE.,.TRUE.,canAB,W0(iscr),nscr)
*rak:c
*rak:c        write(6,*)' hfd1: <1> nscr:',nscr
*rak:        call hf1d(
*rak:     &   coords(1,ictra,geom),exndcf(ipa,basis),
*rak:     &                      exndcf(ica,basis),npa,nca,la,ictra,
*rak:     &   coords(1,ictrb,geom),exndcf(ipb,basis),
*rak:     &                      exndcf(icb,basis),npb,ncb,lb,ictrb,
*rak:     &   coords(1,1,geom),charge(1,geom),ncenters,
*rak:     &          W0(i_Sa),W0(i_Ta),W0(i_Va),nint,
*rak:     &          .TRUE.,.TRUE.,.TRUE.,canAB,.false.,W0(iscr),nscr)
*rak:        
*rak:                     itime = 0
*rak:
*rak:c Numeric derivatives of 2-ctr OIs and KEIs.
*rak:
*rak:c Set finite difference perturbation.
*rak:
*rak:        DELTA = 1.D-06
*rak:
*rak:c Set coordinates for centers.
*rak:
*rak:        Txyz(1,1) = coords(1,ictra,geom) ! xyz(1,ictra)
*rak:        Txyz(2,1) = coords(2,ictra,geom) ! xyz(2,ictra)
*rak:        Txyz(3,1) = coords(3,ictra,geom) ! xyz(3,ictra)
*rak:
*rak:        Txyz(1,2) = coords(1,ictrb,geom) ! xyz(1,ictrb)
*rak:        Txyz(2,2) = coords(2,ictrb,geom) ! xyz(2,ictrb)
*rak:        Txyz(3,2) = coords(3,ictrb,geom) ! xyz(3,ictrb)
*rak:
*rak:c Loop over the centers.
*rak:
*rak:        do 120 i2 = 1,2
*rak:
*rak:c Cartesian components of the 1st derivative.
*rak:
*rak:         do 110 i1 = 1,3
*rak:
*rak:c Plus delta.
*rak:
*rak:          Txyz(i1,i2) = Txyz(i1,i2) + DELTA
*rak:
*rak:c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
*rak:c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
*rak:c       &             xyz,zan,ncenters,
*rak:c       &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
*rak:c       &             .TRUE.,.TRUE.,.FALSE.,canAB,W0(iscr),nscr)
*rak:c          write(6,*)' hf1: <1> nscr:',nscr
*rak:          call hf1(
*rak:     &   txyz(1,1),exndcf(ipa,basis),
*rak:     &                      exndcf(ica,basis),npa,nca,la,
*rak:     &   txyz(1,2),exndcf(ipb,basis),
*rak:     &                      exndcf(icb,basis),npb,ncb,lb,
*rak:     &   coords(1,1,geom),charge(1,geom),ncenters,
*rak:     &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
*rak:     &             .TRUE.,.TRUE.,.FALSE.,canAB,
*rak:     &           .false.,W0(iscr),nscr)
*rak:
*rak:          Txyz(i1,i2) = Txyz(i1,i2) - DELTA
*rak:
*rak:c Minus delta.
*rak:
*rak:          Txyz(i1,i2) = Txyz(i1,i2) - DELTA
*rak:
*rak:c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
*rak:c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
*rak:c       &             xyz,zan,ncenters,
*rak:c       &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
*rak:c       &             .TRUE.,.TRUE.,.FALSE.,canAB,W0(iscr),nscr)
*rak:c          write(6,*)' hf1: <2> nscr:',nscr
*rak:          call hf1(
*rak:     &   txyz(1,1),exndcf(ipa,basis),
*rak:     &                      exndcf(ica,basis),npa,nca,la,
*rak:     &   txyz(1,2),exndcf(ipb,basis),
*rak:     &                      exndcf(icb,basis),npb,ncb,lb,
*rak:     &   coords(1,1,geom),charge(1,geom),ncenters,
*rak:     &           W0(i_S2),W0(i_T2),W0(i_V2),nint,
*rak:     &           .TRUE.,.TRUE.,.FALSE.,canAB,
*rak:     &           .false.,W0(iscr),nscr)
*rak:          Txyz(i1,i2) = Txyz(i1,i2) + DELTA
*rak:
*rak:c Compute the 1st derivative of the integrals via a central difference formula.
*rak:
*rak:          call cdif(W0(i_S1),W0(i_S2),DELTA,nint)
*rak:          call cdif(W0(i_T1),W0(i_T2),DELTA,nint)
*rak:
*rak:c Check accuracy of derivatives.
*rak:
*rak:          i_Sa = i_S0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:          i_Ta = i_T0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:
*rak:          adifS = 0.D0
*rak:          adifT = 0.D0
*rak:          do 100 i = 1,nint
*rak:           adifS = max( adifS, abs( W0(i_Sa+(i-1)) - W0(i_S1+(i-1)) ) )
*rak:           adifT = max( adifT, abs( W0(i_Ta+(i-1)) - W0(i_T1+(i-1)) ) )
*rak:  100     continue
*rak:
*rak:          if( adifS.gt.1.D-08)then
*rak:           write(21,'(a,e16.8)') 'max. abs. dif. ',adifS
*rak:           write(21,'(3i5)') isha,ictra,la
*rak:           write(21,'(3i5)') ishb,ictrb,lb
*rak:           write(21,'(3i5)') i1,i2
*rak:           write(21,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Sa+(i-1)),  W0(i_S1+(i-1)),
*rak:     &                    (W0(i_Sa+(i-1)) - W0(i_S1+(i-1))), i=1,nint)
*rak:           write(6,'(a,e16.8)') 'S:max. abs. dif. ',adifS
*rak:           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
*rak:           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
*rak:           write(6,'(a,3i5)') 'i1,i2',i1,i2
*rak:           write(6,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Sa+(i-1)),  W0(i_S1+(i-1)),
*rak:     &                    (W0(i_Sa+(i-1)) - W0(i_S1+(i-1))), i=1,nint)
*rak:c          stop
*rak:          end if
*rak:
*rak:          if( adifT.gt.1.D-08)then
*rak:           write(31,'(a,e16.8)') 'max. abs. dif. ',adifT
*rak:           write(31,'(3i5)') isha,ictra,la
*rak:           write(31,'(3i5)') ishb,ictrb,lb
*rak:           write(31,'(3i5)') i1,i2
*rak:           write(31,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Ta+(i-1)),  W0(i_T1+(i-1)),
*rak:     &                    (W0(i_Ta+(i-1)) - W0(i_T1+(i-1))), i=1,nint)
*rak:           write(6,'(a,e16.8)') 'T max. abs. dif. ',adifT
*rak:           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
*rak:           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
*rak:           write(6,'(a,3i5)') 'i1,i2',i1,i2
*rak:           write(6,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Ta+(i-1)),  W0(i_T1+(i-1)),
*rak:     &                    (W0(i_Ta+(i-1)) - W0(i_T1+(i-1))), i=1,nint)
*rak:c          stop
*rak:          end if
*rak:
*rak:  110    continue
*rak:
*rak:  120   continue
*rak:
*rak:
*rak:c Numeric derivatives of NAIs.
*rak:
*rak:c Loop over the centers.
*rak:
*rak:        do 150 i2 = 1,ncenters
*rak:
*rak:c Cartesian components of the 1st derivative.
*rak:
*rak:         do 140 i1 = 1,3
*rak:
*rak:c         if( (i2.eq.ictra) .or. (i2.eq.ictrb) ) go to 150
*rak:
*rak:c Plus delta.
*rak:
*rak:          if( i2.eq.ictra )then
*rak:           Txyz(i1,1) = Txyz(i1,1) + DELTA
*rak:          else if( i2.eq.ictrb )then
*rak:           Txyz(i1,2) = Txyz(i1,2) + DELTA
*rak:          end if
*rak:
*rak:c          xyz(i1,i2) = xyz(i1,i2) + DELTA
*rak:          coords(i1,i2,geom) = coords(i1,i2,geom) + DELTA
*rak:
*rak:c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
*rak:c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
*rak:c       &             xyz,zan,ncenters,
*rak:c       &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
*rak:c       &             .FALSE.,.FALSE.,.TRUE.,canAB,W0(iscr),nscr)
*rak:c          write(6,*)' hf1: <3> nscr:',nscr
*rak:          call hf1(
*rak:     &   txyz(1,1),exndcf(ipa,basis),
*rak:     &                      exndcf(ica,basis),npa,nca,la,
*rak:     &   txyz(1,2),exndcf(ipb,basis),
*rak:     &                      exndcf(icb,basis),npb,ncb,lb,
*rak:     &   coords(1,1,geom),charge(1,geom),ncenters,
*rak:     &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
*rak:     &             .FALSE.,.FALSE.,.TRUE.,canAB,
*rak:     &           .false.,W0(iscr),nscr)
*rak:          if( i2.eq.ictra )then
*rak:           Txyz(i1,1) = Txyz(i1,1) - DELTA
*rak:          else if( i2.eq.ictrb )then
*rak:           Txyz(i1,2) = Txyz(i1,2) - DELTA
*rak:          end if
*rak:
*rak:c          xyz(i1,i2) = xyz(i1,i2) - DELTA
*rak:          coords(i1,i2,geom) = coords(i1,i2,geom) - DELTA
*rak:
*rak:c Minus delta.
*rak:
*rak:          if( i2.eq.ictra )then
*rak:           Txyz(i1,1) = Txyz(i1,1) - DELTA
*rak:          else if( i2.eq.ictrb )then
*rak:           Txyz(i1,2) = Txyz(i1,2) - DELTA
*rak:          end if
*rak:
*rak:c          xyz(i1,i2) = xyz(i1,i2) - DELTA
*rak:          coords(i1,i2,geom) = coords(i1,i2,geom) - DELTA
*rak:
*rak:c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
*rak:c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
*rak:c       &             xyz,zan,ncenters,
*rak:c       &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
*rak:c       &             .FALSE.,.FALSE.,.TRUE.,canAB,W0(iscr),nscr)
*rak:c          write(6,*)' hf1: <4> nscr:',nscr
*rak:          call hf1(
*rak:     &   txyz(1,1),exndcf(ipa,basis),
*rak:     &                      exndcf(ica,basis),npa,nca,la,
*rak:     &   txyz(1,2),exndcf(ipb,basis),
*rak:     &                      exndcf(icb,basis),npb,ncb,lb,
*rak:     &   coords(1,1,geom),charge(1,geom),ncenters,
*rak:     &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
*rak:     &             .FALSE.,.FALSE.,.TRUE.,canAB,
*rak:     &           .false.,W0(iscr),nscr)
*rak:
*rak:          if( i2.eq.ictra )then
*rak:           Txyz(i1,1) = Txyz(i1,1) + DELTA
*rak:          else if( i2.eq.ictrb )then
*rak:           Txyz(i1,2) = Txyz(i1,2) + DELTA
*rak:          end if
*rak:
*rak:c          xyz(i1,i2) = xyz(i1,i2) + DELTA
*rak:          coords(i1,i2,geom) = coords(i1,i2,geom) + DELTA
*rak:
*rak:c Compute the 1st derivative of the integrals via a central difference formula.
*rak:
*rak:          call cdif(W0(i_V1),W0(i_V2),DELTA,nint)
*rak:
*rak:c Check accuracy of derivatives.
*rak:
*rak:          i_Va = i_V0 + ( ((i2-1)*3 + i1) - 1 )*nint
*rak:
*rak:          adifV = 0.D0
*rak:          do 130 i = 1,nint
*rak:           adifV = max( adifV, abs( W0(i_Va+(i-1)) - W0(i_V1+(i-1)) ) )
*rak:  130     continue
*rak:
*rak:          if( adifV.gt.1.D-08)then
*rak:           write(41,'(a,e16.8)') 'max. abs. dif. ',adifV
*rak:           write(41,'(3i5)') isha,ictra,la
*rak:           write(41,'(3i5)') ishb,ictrb,lb
*rak:           write(41,'(3i5)') i1,i2
*rak:           write(41,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Va+(i-1)),  W0(i_V1+(i-1)),
*rak:     &                    (W0(i_Va+(i-1)) - W0(i_V1+(i-1))), i=1,nint)
*rak:           write(6,'(a,e16.8)') 'V max. abs. dif. ',adifV
*rak:           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
*rak:           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
*rak:           write(6,'(a,3i5)') 'i1,i2',i1,i2
*rak:           write(6,'(i10,3e22.14)')
*rak:     &                 (i, W0(i_Va+(i-1)),  W0(i_V1+(i-1)),
*rak:     &                    (W0(i_Va+(i-1)) - W0(i_V1+(i-1))), i=1,nint)
*rak:c          stop
*rak:          end if
*rak:
*rak:  140    continue
*rak:
*rak:  150   continue
*rak:
*rak:c  c Determine the indices for this block of integrals.
*rak:c  
*rak:c          ioffseta = kloc(isha) - 1
*rak:c          ioffsetb = kloc(ishb) - 1
*rak:c  
*rak:c          nn = -1
*rak:c  
*rak:c          do 170 ma = 1,nshbfa
*rak:c  
*rak:c           if( canAB )then
*rak:c            mb_limit = ma
*rak:c           else
*rak:c            mb_limit = nshbfb
*rak:c           end if
*rak:c  
*rak:c           do 160 mb = 1,mb_limit
*rak:c  
*rak:c            nn = nn + 1
*rak:c  
*rak:c            i = ioffseta + ma
*rak:c            j = ioffsetb + mb
*rak:c  
*rak:c            ij = (i*(i-1))/2 + j
*rak:c  
*rak:c    160    continue
*rak:c  
*rak:c    170   continue
*rak:
*rak:  180  continue
*rak:
*rak:  190 continue
*rak:
*rak:      close(1)
*rak:
*rak:      itime = 1
*rak:c      if( itime.eq.1 )then
*rak:c       write(ltime,9999) 'hfset    ',time(1)
*rak:c       write(ltime,9999) 'hfmke    ',time(2)
*rak:c       write(ltime,9999) 'O2Is     ',time(3)
*rak:c       write(ltime,9999) 'KEIs     ',time(4)
*rak:c       write(ltime,9999) 'NAIs     ',time(5)
*rak:c      end if
*rak:
*rak: 9999 format('  time spent in ',a10,20x,f20.6)
*rak:
*rak:      end
*rak:c.......................................................................
*rak:      Subroutine cdif(f1,f2,delta,nfcns)
*rak:
*rak:      Implicit real*8 (a-h,o-z)
*rak:      Implicit integer (i-n)
*rak:
*rak:      Dimension f1(nfcns),f2(nfcns)
*rak:c
*rak:c Compute 1st derivatives via a central difference formula.
*rak:c
*rak:c******************************************************************************
*rak:
*rak:      do 10 n = 1,nfcns
*rak:       f1(n) = (f1(n) - f2(n))/(2.D0*delta)
*rak:   10 continue
*rak:
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildSa(Sa,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision sa(nbf,nbf,3*nat),scr(nscr)
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, ids, ide, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF, compute
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildSa: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      call dfill(nbf*nbf*3*nat,0.0,sa,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh = 1,ish
*rak:          iatom = ibs_cn2ce(ish,basis)
*rak:          jatom = ibs_cn2ce(jsh,basis)
*rak:          compute = iatom.ne.jatom
*rak:          if (compute) then
*rak:            ucont = ibs_cn2ucn(ish,basis)
*rak:            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:            
*rak:            ucont = ibs_cn2ucn(jsh,basis)
*rak:            jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:                
*rak:            nint = (itype+1)*(itype+2)/2
*rak:            nint = nint * (jtype+1)*(jtype+2)/2
*rak:            used = nint*6
*rak:            iscr = used + 1
*rak:            call hf1d(
*rak:     &             coords(1,iatom,geom),
*rak:     &             exndcf(iexp,basis),
*rak:     &             exndcf(icf,basis),
*rak:     &             inp,igen,itype,iatom,
*rak:c
*rak:     &             coords(1,jatom,geom),
*rak:     &             exndcf(jexp,basis),
*rak:     &             exndcf(jcf,basis),
*rak:     &             jnp,jgen,jtype,jatom,
*rak:c
*rak:     &             coords(1,1,geom),charge(1,geom),nat,
*rak:     &             scr,scr,scr,nint,
*rak:     &             FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
*rak:            icount = 0
*rak:            ids = (iatom-1)*3+1
*rak:            ide = ids + 2
*rak:            do 00300 id = ids,ide
*rak:              do 00400 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00500 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  sa(ibf,jbf,id) = scr(icount)
*rak:                  sa(jbf,ibf,id) = scr(icount)
*rak:00500           continue
*rak:00400         continue
*rak:00300       continue
*rak:            ids = (jatom-1)*3+1
*rak:            ide = ids + 2
*rak:            do 00600 id = ids,ide
*rak:              do 00700 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &               ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  sa(ibf,jbf,id) = scr(icount)
*rak:                  sa(jbf,ibf,id) = scr(icount)
*rak:00800           continue
*rak:00700         continue
*rak:00600       continue
*rak:          endif
*rak:00200   continue
*rak:00100 continue
*rak:
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Sa <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Sa ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildSafd(Sa,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision sa(nbf,nbf,3*nat),scr(nscr)
*rak:      double precision bplus(nbf,nbf),bminus(nbf,nbf)
*rak:      double precision b(nbf,nbf)
*rak:      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
*rak:      double precision delta
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildSafd: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      delta = 1.0d-06
*rak:      call dfill(nbf*nbf*3*nat,0.0,sa,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      id = 0
*rak:      do 00050 datom = 1,nat
*rak:        do 00060 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,*)' Sa:fd <atom:',datom,'> <xyz:',dxyz,'>'
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
*rak:          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
*rak:          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*rak:*          call print_coords('   coords',coords(1,1,geom),nat)
*rak:*          call print_coords(' + coords',xyzxyzp,nat)
*rak:*          call print_coords(' - coords',xyzxyzm,nat)
*rak:          do 00100 ish = 1,ncont
*rak:            do 00200 jsh = 1,ncont
*rak:              ucont = ibs_cn2ucn(ish,basis)
*rak:              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              iatom = ibs_cn2ce(ish,basis)
*rak:              
*rak:              ucont = ibs_cn2ucn(jsh,basis)
*rak:              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              jatom = ibs_cn2ce(jsh,basis)
*rak:              
*rak:              nint = (itype+1)*(itype+2)/2
*rak:              nint = nint * (jtype+1)*(jtype+2)/2
*rak:              used = nint
*rak:              iscr = used + 1
*rak:c
*rak:c... do center
*rak:              call hf1(
*rak:     &               coords(1,iatom,geom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               coords(1,jatom,geom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               coords(1,1,geom),charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  b(jbf,ibf) = scr(icount)
*rak:00401           continue
*rak:00301         continue
*rak:c
*rak:c... do plus
*rak:              call dfill(nscr,0.0,scr,1)
*rak:              call hf1(
*rak:     &               xyzxyzp(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzp(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzp,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bplus(jbf,ibf) = scr(icount)
*rak:00400           continue
*rak:00300         continue
*rak:c
*rak:c........ do minus
*rak:              call dfill(nscr,0.0,scr,1)
*rak:c
*rak:              call hf1(
*rak:     &               xyzxyzm(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzm(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzm,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bminus(jbf,ibf) = scr(icount)
*rak:00600           continue
*rak:00500         continue
*rak:00200       continue
*rak:00100     continue
*rak:*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
*rak:          call cdif(bplus,bminus,delta,nbf*nbf)
*rak:          call dcopy(nbf*nbf,bplus,1,sa(1,1,id),1)
*rak:00060   continue
*rak:00050 continue
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Sa <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Sa FD ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c.......................................................................
*rak:      subroutine print_coords(msg,xyzxyz,nat)
*rak:      implicit none
*rak:      character*(*) msg
*rak:      integer nat, i, j
*rak:      double precision xyzxyz(3,nat)
*rak:c
*rak:      write(6,*)' '
*rak:      write(6,*)' '
*rak:      write(6,10000)
*rak:      write(6,'(1x,a)')msg
*rak:      write(6,10000)
*rak:      do 00100 i=1,nat
*rak:        write(6,10001)i,(xyzxyz(j,i),j=1,3)
*rak:00100 continue
*rak:      write(6,10000)
*rak:      return
*rak:10000 format(100('-'))
*rak:10001 format(1x,'<atom:',i3,'>  <x:',1pd20.10,'>  <y:',1pd20.10,
*rak:     &       '>  <z:',1pd20.10,'>')
*rak:      end
*rak:c.......................................................................
*rak:      subroutine myoutput (z,rowlow,rowhi,collow,colhi,rowdim,coldim,
*rak:     $     nctl,msg,zerotol)
*rak:c.......................................................................
*rak:c myoutput prints a real*8 matrix in formatted form with numbered rows
*rak:c and columns.  the input is as follows;
*rak:c        matrix(*,*).........matrix to be output
*rak:c        rowlow..............row number at which output is to begin
*rak:c        rowhi...............row number at which output is to end
*rak:c        collow..............column number at which output is to begin
*rak:c        colhi...............column number at which output is to end
*rak:c        rowdim..............row dimension of matrix(*,*)
*rak:c        coldim..............column dimension of matrix(*,*)
*rak:c        nctl................carriage control flag; 1 for single space
*rak:c                                                   2 for double space
*rak:c                                                   3 for triple space
*rak:c        msg.................informative message
*rak:c        zerotol.............information
*rak:c the parameters that follow matrix are all of type integer.  the
*rak:c program is set up to handle 5 columns/page with a 1p5d24.15 format for
*rak:c the columns.  if a different number of columns is required, change
*rak:c formats 1000 and 2000, and initialize kcol with the new number of
*rak:c columns.
*rak:c author;  nelson h.f. beebe, quantum theory project, university of
*rak:c          florida, gainesville
*rak:c modified: RJ Harrison - 8 columns 
*rak:c modified: Rick A. Kendall, HPCC group, Pacific Northwest Laboratory
*rak:c.......................................................................
*rak:      implicit none
*rak:      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
*rak:      integer nctl, i, j, last, k
*rak:      double precision z(rowdim,coldim)
*rak:      double precision zerotol
*rak:      character*(*) msg
*rak:      character*8 asa(3), column, ctl, blank
*rak:      data column/'column'   /,asa/' ','00000000'  ,
*rak:     1     '--------'  /,blank/' '/
*rak:      data kcol/8/
*rak:c
*rak:      write(6,4000)msg
*rak: 4000 format(1x,80('-'),/,2x,a,/,1x,80('-'),/)
*rak:c
*rak:      do 11 i=rowlow,rowhi
*rak:         do 10 j=collow,colhi
*rak:            if (abs(z(i,j)).ge.zerotol) go to 15
*rak: 10      continue
*rak: 11   continue
*rak:      write (6,3000)
*rak: 3000 format (/' zero matrix'/)
*rak:      go to 3
*rak: 15   continue
*rak:      ctl = blank
*rak:      if ((nctl.le.3).and.(nctl.gt.0)) ctl = asa(nctl)
*rak:      if (rowhi.lt.rowlow) go to 3
*rak:      if (colhi.lt.collow) go to 3
*rak:      last = min(colhi,collow+kcol-1)
*rak:      do 2 begin = collow,colhi,kcol
*rak:         write (6,1000) (column,i,i = begin,last)
*rak:         do 1 k = rowlow,rowhi
*rak:            do 4 i=begin,last
*rak:               if (abs(z(k,i)).ge.zerotol) go to 5
*rak: 4          continue
*rak:            go to 1
*rak: 5          write (6,2000) ctl,k,(z(k,i), i = begin,last)
*rak: 1       continue
*rak:         last = min(last+kcol,colhi)
*rak: 2    continue
*rak: 3    return
*rak:* kcol = 4
*rak:* 1000 format (/1h ,16x,3(a6,i3,2x),(a6,i3))
*rak:* 2000 format (a1,3hrow,i4,2x,4f17.11)
*rak:* kcol = 8
*rak: 1000 format (/1h ,11x,7(a3,i3,3x),(a3,i3))
*rak: 2000 format (a1,'row',i4,1x,8f9.4)
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildTa(Ta,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision Ta(nbf,nbf,3*nat),scr(nscr)
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, ids, ide, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF, compute
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildTa: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      call dfill(nbf*nbf*3*nat,0.0,Ta,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh = 1,ish
*rak:          iatom = ibs_cn2ce(ish,basis)
*rak:          jatom = ibs_cn2ce(jsh,basis)
*rak:          compute = iatom.ne.jatom
*rak:          if (compute) then
*rak:            ucont = ibs_cn2ucn(ish,basis)
*rak:            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:            
*rak:            ucont = ibs_cn2ucn(jsh,basis)
*rak:            jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:                
*rak:            nint = (itype+1)*(itype+2)/2
*rak:            nint = nint * (jtype+1)*(jtype+2)/2
*rak:            used = nint*6
*rak:            iscr = used + 1
*rak:            call hf1d(
*rak:     &             coords(1,iatom,geom),
*rak:     &             exndcf(iexp,basis),
*rak:     &             exndcf(icf,basis),
*rak:     &             inp,igen,itype,iatom,
*rak:c
*rak:     &             coords(1,jatom,geom),
*rak:     &             exndcf(jexp,basis),
*rak:     &             exndcf(jcf,basis),
*rak:     &             jnp,jgen,jtype,jatom,
*rak:c
*rak:     &             coords(1,1,geom),charge(1,geom),nat,
*rak:     &             scr,scr,scr,nint,
*rak:     &             FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
*rak:            icount = 0
*rak:            ids = (iatom-1)*3+1
*rak:            ide = ids + 2
*rak:            do 00300 id = ids,ide
*rak:              do 00400 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00500 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  Ta(ibf,jbf,id) = scr(icount)
*rak:                  Ta(jbf,ibf,id) = scr(icount)
*rak:00500           continue
*rak:00400         continue
*rak:00300       continue
*rak:            ids = (jatom-1)*3+1
*rak:            ide = ids + 2
*rak:            do 00600 id = ids,ide
*rak:              do 00700 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  Ta(ibf,jbf,id) = scr(icount)
*rak:                  Ta(jbf,ibf,id) = scr(icount)
*rak:00800           continue
*rak:00700         continue
*rak:00600       continue
*rak:          endif
*rak:00200   continue
*rak:00100 continue
*rak:
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Ta <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(Ta(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Ta ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildTafd(Ta,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision Ta(nbf,nbf,3*nat),scr(nscr)
*rak:      double precision bplus(nbf,nbf),bminus(nbf,nbf)
*rak:      double precision b(nbf,nbf)
*rak:      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
*rak:      double precision delta
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildTafd: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      delta = 1.0d-06
*rak:      call dfill(nbf*nbf*3*nat,0.0,Ta,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      id = 0
*rak:      do 00050 datom = 1,nat
*rak:        do 00060 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,*)' Ta:fd <atom:',datom,'> <xyz:',dxyz,'>'
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
*rak:          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
*rak:          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*rak:*          call print_coords('   coords',coords(1,1,geom),nat)
*rak:*          call print_coords(' + coords',xyzxyzp,nat)
*rak:*          call print_coords(' - coords',xyzxyzm,nat)
*rak:          do 00100 ish = 1,ncont
*rak:            do 00200 jsh = 1,ncont
*rak:              ucont = ibs_cn2ucn(ish,basis)
*rak:              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              iatom = ibs_cn2ce(ish,basis)
*rak:              
*rak:              ucont = ibs_cn2ucn(jsh,basis)
*rak:              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              jatom = ibs_cn2ce(jsh,basis)
*rak:              
*rak:              nint = (itype+1)*(itype+2)/2
*rak:              nint = nint * (jtype+1)*(jtype+2)/2
*rak:              used = nint
*rak:              iscr = used + 1
*rak:c
*rak:c... do center
*rak:              call hf1(
*rak:     &               coords(1,iatom,geom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               coords(1,jatom,geom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               coords(1,1,geom),charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  b(jbf,ibf) = scr(icount)
*rak:00401           continue
*rak:00301         continue
*rak:c
*rak:c... do plus
*rak:              call dfill(nscr,0.0,scr,1)
*rak:              call hf1(
*rak:     &               xyzxyzp(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzp(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzp,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bplus(jbf,ibf) = scr(icount)
*rak:00400           continue
*rak:00300         continue
*rak:c
*rak:c........ do minus
*rak:              call dfill(nscr,0.0,scr,1)
*rak:c
*rak:              call hf1(
*rak:     &               xyzxyzm(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzm(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzm,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bminus(jbf,ibf) = scr(icount)
*rak:00600           continue
*rak:00500         continue
*rak:00200       continue
*rak:00100     continue
*rak:*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
*rak:          call cdif(bplus,bminus,delta,nbf*nbf)
*rak:          call dcopy(nbf*nbf,bplus,1,Ta(1,1,id),1)
*rak:00060   continue
*rak:00050 continue
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Ta <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(Ta(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Ta FD ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildVa(Va,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision Va(nbf,nbf,3*nat),scr(nscr)
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildVa: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      call dfill(nbf*nbf*3*nat,0.0,Va,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh = 1,ish
*rak:          iatom = ibs_cn2ce(ish,basis)
*rak:          jatom = ibs_cn2ce(jsh,basis)
*rak:
*rak:          ucont = ibs_cn2ucn(ish,basis)
*rak:          inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:          igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:          iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:          icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:          itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:          
*rak:          ucont = ibs_cn2ucn(jsh,basis)
*rak:          jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:          jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:          jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:          jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:          jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:          
*rak:          nint = (itype+1)*(itype+2)/2
*rak:          nint = nint * (jtype+1)*(jtype+2)/2
*rak:          used = nint*3*nat
*rak:          iscr = used + 1
*rak:          call hf1d(
*rak:     &           coords(1,iatom,geom),
*rak:     &           exndcf(iexp,basis),
*rak:     &           exndcf(icf,basis),
*rak:     &           inp,igen,itype,iatom,
*rak:c
*rak:     &           coords(1,jatom,geom),
*rak:     &           exndcf(jexp,basis),
*rak:     &           exndcf(jcf,basis),
*rak:     &           jnp,jgen,jtype,jatom,
*rak:c
*rak:     &           coords(1,1,geom),charge(1,geom),nat,
*rak:     &           scr,scr,scr,nint,
*rak:     &           FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
*rak:          icount = 0
*rak:          id = 0
*rak:          do 00300 datom = 1,nat
*rak:            do 00400 dxyz = 1,3
*rak:              id = id + 1
*rak:              do 00500 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  Va(ibf,jbf,id) = scr(icount)
*rak:                  Va(jbf,ibf,id) = scr(icount)
*rak:00600           continue
*rak:00500         continue
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Va <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(Va(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Va ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildVafd(Va,nat,ncont,nbf,basisin,geom,
*rak:     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basisin,nscr,geom
*rak:      double precision Va(nbf,nbf,3*nat),scr(nscr)
*rak:      double precision bplus(nbf,nbf),bminus(nbf,nbf)
*rak:      double precision b(nbf,nbf)
*rak:      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
*rak:      double precision delta
*rak:      logical printit
*rak:c
*rak:      integer ish,jsh, basis, id, icount
*rak:      integer testgeom, ucont, nint, iscr, used, jbf, ibf
*rak:      integer datom,dxyz
*rak:      integer iatom,inp,igen,iexp,icf,itype
*rak:      integer jatom,jnp,jgen,jexp,jcf,jtype
*rak:      logical FT, FF
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:      testgeom = ibs_geom(basis)
*rak:      if (testgeom.ne.geom) stop 'buildVafd: geom does not match'
*rak:c
*rak:      FT = .true.
*rak:      FF = .false.
*rak:      delta = 1.0d-06
*rak:      call dfill(nbf*nbf*3*nat,0.0,Va,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      id = 0
*rak:      do 00050 datom = 1,nat
*rak:        do 00060 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,*)' Va:fd <atom:',datom,'> <xyz:',dxyz,'>'
*rak:#if defined(SUN)
*rak:      call flush(6)
*rak:#endif
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
*rak:          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
*rak:          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*rak:*          call print_coords('   coords',coords(1,1,geom),nat)
*rak:*          call print_coords(' + coords',xyzxyzp,nat)
*rak:*          call print_coords(' - coords',xyzxyzm,nat)
*rak:          do 00100 ish = 1,ncont
*rak:            do 00200 jsh = 1,ncont
*rak:              ucont = ibs_cn2ucn(ish,basis)
*rak:              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              iatom = ibs_cn2ce(ish,basis)
*rak:              
*rak:              ucont = ibs_cn2ucn(jsh,basis)
*rak:              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:              jatom = ibs_cn2ce(jsh,basis)
*rak:              
*rak:              nint = (itype+1)*(itype+2)/2
*rak:              nint = nint * (jtype+1)*(jtype+2)/2
*rak:              used = nint
*rak:              iscr = used + 1
*rak:c
*rak:c... do center
*rak:              call hf1(
*rak:     &               coords(1,iatom,geom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               coords(1,jatom,geom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               coords(1,1,geom),charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  b(jbf,ibf) = scr(icount)
*rak:00401           continue
*rak:00301         continue
*rak:c
*rak:c... do plus
*rak:              call dfill(nscr,0.0,scr,1)
*rak:              call hf1(
*rak:     &               xyzxyzp(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzp(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzp,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bplus(jbf,ibf) = scr(icount)
*rak:00400           continue
*rak:00300         continue
*rak:c
*rak:c........ do minus
*rak:              call dfill(nscr,0.0d00,scr,1)
*rak:c
*rak:              call hf1(
*rak:     &               xyzxyzm(1,iatom),
*rak:     &               exndcf(iexp,basis),
*rak:     &               exndcf(icf,basis),
*rak:     &               inp,igen,itype,
*rak:c
*rak:     &               xyzxyzm(1,jatom),
*rak:     &               exndcf(jexp,basis),
*rak:     &               exndcf(jcf,basis),
*rak:     &               jnp,jgen,jtype,
*rak:c
*rak:     &               xyzxyzm,charge(1,geom),nat,
*rak:     &               scr,scr,scr,nint,
*rak:     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
*rak:              icount = 0
*rak:              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
*rak:     &               ibs_cn2bfr(2,ish,basis)
*rak:                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                 ibs_cn2bfr(2,jsh,basis)
*rak:                  icount = icount + 1
*rak:                  bminus(jbf,ibf) = scr(icount)
*rak:00600           continue
*rak:00500         continue
*rak:00200       continue
*rak:00100     continue
*rak:*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*rak:*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*rak:*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
*rak:          call cdif(bplus,bminus,delta,nbf*nbf)
*rak:          call dcopy(nbf*nbf,bplus,1,Va(1,1,id),1)
*rak:00060   continue
*rak:00050 continue
*rak:      if (.not.printit) return
*rak:      id = 0
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' Va <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(Va(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' Va FD ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:      subroutine printdmat(msg,Xa,nat,nbf)
*rak:      implicit none
*rak:c
*rak:      character*(*) msg
*rak:      integer nat,nbf
*rak:      double precision Xa(nbf,nbf,3*nat)
*rak:c
*rak:      integer datom, dxyz, id
*rak:c
*rak:      id = 0
*rak:      do 00100 datom = 1,nat
*rak:        do 00200 dxyz = 1,3
*rak:          id = id + 1
*rak:          write(6,'(//,a,a,i2,a,i2,a)')
*rak:     &           msg,' <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(Xa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           msg,1.0d-05)
*rak:00200   continue
*rak:00100 continue
*rak:      return
*rak:      end
*rak:      subroutine build_d2e(geom,basisin,nbf,ncont,nat,intsize,d2e,
*rak:     &       maxg,gbuf,nscr,scr)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c::passed
*rak:      integer geom,basisin,nbf,nat,intsize,maxg,nscr,ncont
*rak:      double precision d2e(intsize,3,nat)
*rak:      double precision gbuf(maxg)
*rak:      double precision scr(nscr)
*rak:c::local
*rak:      integer ucont,basis,nint,iduse,ipfrom,ipto,icount,id
*rak:      integer idatom(4)
*rak:      integer index,ixyz
*rak:c
*rak:C     SET-UP FUNCTIONS
*rak:C
*rak:      integer ibf, ish, iatom, inp, igen, iexp, icf, itype
*rak:      integer jbf, jsh, jatom, jnp, jgen, jexp, jcf, jtype
*rak:      integer kbf, ksh, katom, knp, kgen, kexp, kcf, ktype
*rak:      integer lbf, lsh, latom, lnp, lgen, lexp, lcf, ltype
*rak:c
*rak:c...  statement function
*rak:C
*rak:C     SET-UP FUNCTIONS
*rak:      integer min,max
*rak:      intrinsic min,max
*rak:      integer i,j,k,l,isym2,isym4
*rak:      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
*rak:      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
*rak:     &               (max(isym2(i,j),isym2(k,l))-1)/2+
*rak:     &               min(isym2(i,j),isym2(k,l))
*rak:c
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:c      
*rak:      do 00100 ish = 1,ncont
*rak:c.
*rak:        iatom = ibs_cn2ce(ish,basis)
*rak:        ucont = ibs_cn2ucn(ish,basis)
*rak:        inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:        igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:        iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:        icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:        itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:        do 00200 jsh= 1,ish
*rak:c.
*rak:          jatom = ibs_cn2ce(jsh,basis)
*rak:          ucont = ibs_cn2ucn(jsh,basis)
*rak:          jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:          jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:          jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:          jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:          jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:          do 00300 ksh = 1,jsh
*rak:c.
*rak:            katom = ibs_cn2ce(ksh,basis)
*rak:            ucont = ibs_cn2ucn(ksh,basis)
*rak:            knp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            kgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            kexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            kcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            ktype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:            do 00400 lsh = 1,ksh
*rak:c.
*rak:              latom = ibs_cn2ce(lsh,basis)
*rak:              ucont = ibs_cn2ucn(lsh,basis)
*rak:              lnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              lgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              lexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              lcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              ltype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:              nint = (itype+1)*(itype+2)/2
*rak:              nint = nint*(jtype+1)*(jtype+2)/2
*rak:              nint = nint*(ktype+1)*(ktype+2)/2
*rak:              nint = nint*(ltype+1)*(ltype+2)/2
*rak:c
*rak:              if (
*rak:     &               iatom.eq.jatom.and.
*rak:     &               jatom.eq.katom.and.
*rak:     &               katom.eq.latom
*rak:     &               ) then
*rak:                continue
*rak:              else
*rak:                call hf2d(
*rak:     &                 coords(1,iatom,geom),exndcf(iexp,basis),
*rak:     &                 exndcf(icf,basis),inp,igen,itype,iatom,
*rak:c
*rak:     &                 coords(1,jatom,geom),exndcf(jexp,basis),
*rak:     &                 exndcf(jcf,basis),jnp,jgen,jtype,jatom,
*rak:c
*rak:     &                 coords(1,katom,geom),exndcf(kexp,basis),
*rak:     &                 exndcf(kcf,basis),knp,kgen,ktype,katom,
*rak:c
*rak:     &                 coords(1,latom,geom),exndcf(lexp,basis),
*rak:     &                 exndcf(lcf,basis),lnp,lgen,ltype,latom,
*rak:c
*rak:     &                 gbuf,nint,.false.,.false.,.false.,.false.,
*rak:     &                 scr,nscr)
*rak:c 
*rak:                call logic_atom(idatom,iatom,jatom,katom,latom)
*rak:c do sums
*rak:                do 00500 iduse=2,4
*rak:                  if(idatom(iduse).gt.0) then 
*rak:                    continue
*rak:                  else
*rak:                    ipto   = (abs(idatom(iduse)) -1)*nint*3 + 1
*rak:                    ipfrom = (iduse-1)*nint*3 + 1
*rak:                    call daxpy(nint*3,1.0d00,
*rak:     &                     gbuf(ipfrom),1,gbuf(ipto),1)
*rak:                  endif
*rak:00500           continue
*rak:c
*rak:c ... fill integral buffers              
*rak:                icount = 0
*rak:                do 00600 iduse = 1,4
*rak:                  if (idatom(iduse).lt.1) then
*rak:                    icount = icount + nint*3
*rak:                  else
*rak:                    id = idatom(iduse)
*rak:                    do 00700 ixyz = 1,3
*rak:                      do 00800 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &                       ibs_cn2bfr(2,ish,basis)
*rak:                        do 00900 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                         ibs_cn2bfr(2,jsh,basis)
*rak:                          do 01000 kbf = ibs_cn2bfr(1,ksh,basis),
*rak:     &                           ibs_cn2bfr(2,ksh,basis)
*rak:                            do 01100 lbf = ibs_cn2bfr(1,lsh,basis),
*rak:     &                             ibs_cn2bfr(2,lsh,basis)
*rak:                              index = isym4(ibf,jbf,kbf,lbf)
*rak:                              icount = icount + 1
*rak:                              d2e(index,ixyz,id) = gbuf(icount)
*rak:01100                       continue
*rak:01000                     continue
*rak:00900                   continue
*rak:00800                 continue
*rak:00700               continue
*rak:                  endif
*rak:00600           continue
*rak:c
*rak:              endif
*rak:c
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:      end
*rak:      subroutine logic_atom(idat,iat,jat,kat,lat)
*rak:      implicit none
*rak:      integer iat,jat,kat,lat
*rak:      integer idat(4)
*rak:c
*rak:      idat(1) = iat
*rak:      idat(2) = jat
*rak:      idat(3) = kat
*rak:      idat(4) = lat
*rak:      if (iat.eq.jat) idat(2) = -1
*rak:      if (iat.eq.kat) idat(3) = -1
*rak:      if (iat.eq.lat) idat(4) = -1
*rak:      if (jat.eq.kat) then
*rak:        if(idat(2).gt.0) then
*rak:          idat(3) = -2
*rak:        else
*rak:          idat(3) = idat(2)
*rak:        endif
*rak:      endif
*rak:      if (jat.eq.lat) then
*rak:        if(idat(2).gt.0) then
*rak:          idat(4) = -2
*rak:        else
*rak:          idat(4) = idat(2)
*rak:        endif
*rak:      endif
*rak:      if (kat.eq.lat) then
*rak:        if(idat(3).gt.0) then
*rak:          idat(4) = -3
*rak:        else
*rak:          idat(4) = idat(3)
*rak:        endif
*rak:      endif
*rak:      end
*rak:      subroutine build_d2efd(geom,basisin,nbf,ncont,nat,intsize,d2e,
*rak:     &       bp,bm,xyzp,xyzm,
*rak:     &       maxg,gbuf,nscr,scr)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geom.fh"
*rak:#include "geomP.fh"
*rak:#include "geobasmapP.fh"
*rak:c::passed
*rak:      integer geom,basisin,nbf,nat,intsize,maxg,nscr,ncont
*rak:      double precision d2e(intsize,3,nat)
*rak:      double precision gbuf(maxg)
*rak:      double precision scr(nscr)
*rak:      double precision bp(intsize),bm(intsize)
*rak:      double precision xyzp(3,nat),xyzm(3,nat)
*rak:c:: local
*rak:      integer basis
*rak:      double precision delta
*rak:      integer dxyz,datom,nint,ucont,index,icount
*rak:      integer ibf, ish, iatom, inp, igen, iexp, icf, itype
*rak:      integer jbf, jsh, jatom, jnp, jgen, jexp, jcf, jtype
*rak:      integer kbf, ksh, katom, knp, kgen, kexp, kcf, ktype
*rak:      integer lbf, lsh, latom, lnp, lgen, lexp, lcf, ltype
*rak:c
*rak:c
*rak:c...  statement function
*rak:C
*rak:C     SET-UP FUNCTIONS
*rak:      integer min,max
*rak:      intrinsic min,max
*rak:      integer i,j,k,l,isym2,isym4
*rak:      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
*rak:      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
*rak:     &               (max(isym2(i,j),isym2(k,l))-1)/2+
*rak:     &               min(isym2(i,j),isym2(k,l))
*rak:c
*rak:      delta = 1.0d-06
*rak:      basis = basisin + BASIS_HANDLE_OFFSET
*rak:c      
*rak:      do 00100 datom = 1,nat
*rak:        do 00200 dxyz = 1,3
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzp,1)
*rak:          call dcopy(3*nat,coords(1,1,geom),1,xyzm,1)
*rak:          xyzp(dxyz,datom) =  xyzp(dxyz,datom) + delta
*rak:          xyzm(dxyz,datom) =  xyzm(dxyz,datom) - delta
*rak:          call print_coords('2e +',xyzp,nat)
*rak:          call print_coords('2e -',xyzm,nat)
*rak:c
*rak:          do 00300 ish = 1,ncont
*rak:c.
*rak:            iatom = ibs_cn2ce(ish,basis)
*rak:            ucont = ibs_cn2ucn(ish,basis)
*rak:            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:            igen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:            icf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:            itype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:            do 00400 jsh = 1,ish
*rak:c.
*rak:              jatom = ibs_cn2ce(jsh,basis)
*rak:              ucont = ibs_cn2ucn(jsh,basis)
*rak:              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:              do 00500 ksh = 1,jsh
*rak:c.
*rak:                katom = ibs_cn2ce(ksh,basis)
*rak:                ucont = ibs_cn2ucn(ksh,basis)
*rak:                knp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:                kgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:                kexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:                kcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:                ktype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:                do 00600 lsh = 1,ksh
*rak:c.
*rak:                  latom = ibs_cn2ce(lsh,basis)
*rak:                  ucont = ibs_cn2ucn(lsh,basis)
*rak:                  lnp   = infbs_cont(CONT_NPRIM,ucont,basis)
*rak:                  lgen  = infbs_cont(CONT_NGEN,ucont,basis)
*rak:                  lexp  = infbs_cont(CONT_IEXP,ucont,basis)
*rak:                  lcf   = infbs_cont(CONT_ICFP,ucont,basis)
*rak:                  ltype = infbs_cont(CONT_TYPE,ucont,basis) 
*rak:c.
*rak:                  nint = (itype+1)*(itype+2)/2
*rak:                  nint = nint*(jtype+1)*(jtype+2)/2
*rak:                  nint = nint*(ktype+1)*(ktype+2)/2
*rak:                  nint = nint*(ltype+1)*(ltype+2)/2
*rak:c
*rak:c... do plus
*rak:                  call dfill(nscr,0.0d00,scr,1)
*rak:                  call dfill(maxg,0.0d00,gbuf,1)
*rak:                  call hf2(
*rak:     &                   xyzp(1,iatom),exndcf(iexp,basis),
*rak:     &                   exndcf(icf,basis),inp,igen,itype,
*rak:c
*rak:     &                   xyzp(1,jatom),exndcf(jexp,basis),
*rak:     &                   exndcf(jcf,basis),jnp,jgen,jtype,
*rak:c
*rak:     &                   xyzp(1,katom),exndcf(kexp,basis),
*rak:     &                   exndcf(kcf,basis),knp,kgen,ktype,
*rak:c
*rak:     &                   xyzp(1,latom),exndcf(lexp,basis),
*rak:     &                   exndcf(lcf,basis),lnp,lgen,ltype,
*rak:c
*rak:     &                   gbuf,nint,.false.,.false.,.false.,.false.,
*rak:     &                   scr,nscr)
*rak:c
*rak:                  icount = 0
*rak:                  do 00700 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &                   ibs_cn2bfr(2,ish,basis)
*rak:                    do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                     ibs_cn2bfr(2,jsh,basis)
*rak:                      do 00900 kbf = ibs_cn2bfr(1,ksh,basis),
*rak:     &                       ibs_cn2bfr(2,ksh,basis)
*rak:                        do 01000 lbf = ibs_cn2bfr(1,lsh,basis),
*rak:     &                         ibs_cn2bfr(2,lsh,basis)
*rak:                          index = isym4(ibf,jbf,kbf,lbf)
*rak:                          icount = icount + 1
*rak:                          bp(index) = gbuf(icount)
*rak:01000                   continue
*rak:00900                 continue
*rak:00800               continue
*rak:00700             continue
*rak:c
*rak:c... do minus
*rak:                  call dfill(nscr,0.0d00,scr,1)
*rak:                  call dfill(maxg,0.0d00,gbuf,1)
*rak:                  call hf2(
*rak:     &                   xyzm(1,iatom),exndcf(iexp,basis),
*rak:     &                   exndcf(icf,basis),inp,igen,itype,
*rak:c
*rak:     &                   xyzm(1,jatom),exndcf(jexp,basis),
*rak:     &                   exndcf(jcf,basis),jnp,jgen,jtype,
*rak:c
*rak:     &                   xyzm(1,katom),exndcf(kexp,basis),
*rak:     &                   exndcf(kcf,basis),knp,kgen,ktype,
*rak:c
*rak:     &                   xyzm(1,latom),exndcf(lexp,basis),
*rak:     &                   exndcf(lcf,basis),lnp,lgen,ltype,
*rak:c
*rak:     &                   gbuf,nint,.false.,.false.,.false.,.false.,
*rak:     &                   scr,nscr)
*rak:c
*rak:                  icount = 0
*rak:                  do 00701 ibf = ibs_cn2bfr(1,ish,basis),
*rak:     &                   ibs_cn2bfr(2,ish,basis)
*rak:                    do 00801 jbf = ibs_cn2bfr(1,jsh,basis),
*rak:     &                     ibs_cn2bfr(2,jsh,basis)
*rak:                      do 00901 kbf = ibs_cn2bfr(1,ksh,basis),
*rak:     &                       ibs_cn2bfr(2,ksh,basis)
*rak:                        do 01001 lbf = ibs_cn2bfr(1,lsh,basis),
*rak:     &                         ibs_cn2bfr(2,lsh,basis)
*rak:                          index = isym4(ibf,jbf,kbf,lbf)
*rak:                          icount = icount + 1
*rak:                          bm(index) = gbuf(icount)
*rak:01001                   continue
*rak:00901                 continue
*rak:00801               continue
*rak:00701             continue
*rak:c
*rak:00600           continue
*rak:00500         continue
*rak:00400       continue
*rak:00300     continue
*rak:          call cdif(bp,bm,delta,intsize)
*rak:          call dcopy(intsize,bp,1,d2e(1,dxyz,datom),1)
*rak:00200   continue
*rak:00100 continue
*rak:c
*rak:      end
*rak:      subroutine print2de_diff(nat,size2e,diff,d2e,d2efd,thresh,nbf)
*rak:      implicit none 
*rak:      integer nat,size2e,index
*rak:      double precision diff (size2e,3,nat)
*rak:      double precision d2e  (size2e,3,nat)
*rak:      double precision d2efd(size2e,3,nat)
*rak:      double precision thresh, ra2fd,  rfd2a
*rak:      integer ibf,jbf,kbf,lbf
*rak:      integer nbf
*rak:c
*rak:      integer dxyz,datom
*rak:
*rak:      do 00100 datom = 1,nat
*rak:        do 00200 dxyz = 1,3
*rak:          do 00300 index = 1,size2e
*rak:            if (abs(diff(index,dxyz,datom)).gt.thresh) then
*rak:              call getfromi(index,nbf,ibf,jbf,kbf,lbf)
*rak:              if (abs(d2efd(index,dxyz,datom)).gt.1.0d-05) then
*rak:                ra2fd = d2e(index,dxyz,datom)/d2efd(index,dxyz,datom)
*rak:              else
*rak:                ra2fd = -565.6589d00
*rak:              endif
*rak:              if (abs(d2e(index,dxyz,datom)).gt.1.0d-05) then
*rak:                rfd2a = d2efd(index,dxyz,datom)/d2e(index,dxyz,datom)
*rak:              else
*rak:                rfd2a = -6.0201023d00
*rak:              endif
*rak:              write(6,10000)
*rak:     &               ibf,jbf,kbf,lbf,
*rak:     &               index,dxyz,datom,diff(index,dxyz,datom),
*rak:     &               d2e(index,dxyz,datom),d2efd(index,dxyz,datom),
*rak:     &               ra2fd,rfd2a
*rak:#if defined(SUN)
*rak:              call flush(6)
*rak:#endif              
*rak:            endif
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:10000 format(' <',i2,i2,'|',i2,i2,'||',i5,'|xyz:',i2,'|a:',i2,'>  ',
*rak:     &       '<diff:',1pd14.6,'> ','<anal:',1pd14.6,'> ',
*rak:     &       '<fd:',1pd14.6,'> <ra2fd:',1pd14.6,'> <rfd2a:',
*rak:     &       1pd14.6,'>')
*rak:      end
*rak:      subroutine getfromi(index,nbf,i,j,k,l)
*rak:      implicit none
*rak:      integer index,nbf,i,j,k,l,lhi
*rak:c
*rak:      integer icount
*rak:c
*rak:      icount = 1
*rak:      do 00100 i=1,nbf
*rak:        do 00200 j=1,i
*rak:          do 00300 k=1,i
*rak:            lhi = k
*rak:            if (k.eq.i) lhi = j
*rak:            do 00400 l=1,lhi
*rak:              if (icount.eq.index) return
*rak:              icount = icount + 1
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:      write(6,*)' index/nbf =',index,nbf
*rak:      stop ' getfromi fucked '
*rak:      end
*rak:      subroutine print2de(nat,size2e,d2e,thresh,nbf,msg)
*rak:      implicit none 
*rak:      integer nat,size2e,index
*rak:      double precision d2e  (size2e,3,nat)
*rak:      double precision thresh
*rak:      integer ibf,jbf,kbf,lbf
*rak:      character*(*)msg
*rak:      integer nbf
*rak:c
*rak:      integer dxyz,datom
*rak:
*rak:      do 00100 datom = 1,nat
*rak:        do 00200 dxyz = 1,3
*rak:          do 00300 index = 1,size2e
*rak:            if (abs(d2e(index,dxyz,datom)).gt.thresh) then
*rak:              call getfromi(index,nbf,ibf,jbf,kbf,lbf)
*rak:              write(6,10000) ibf,jbf,kbf,lbf,
*rak:     &               index,dxyz,datom,msg,d2e(index,dxyz,datom)
*rak:#if defined(SUN)
*rak:              call flush(6)
*rak:#endif              
*rak:            endif
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:10000 format(' <',i2,i2,'|',i2,i2,'||',i5,'|xyz:',i2,'|a:',i2,'>  ',
*rak:     &       '<',a,':',1pd14.6,'> ')
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildSanew(Sa,nat,ncont,nbf,basis,
*rak:     &       scr,nscr,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basis,nscr,idatom(2)
*rak:      double precision sa(nbf,nbf,3,nat),scr(nscr)
*rak:      logical printit
*rak:c
*rak:      integer used, iscr, newnscr, nbfmax
*rak:      integer ish, ibf, ibflo, ibfhi
*rak:      integer jsh, jbf, jbflo, jbfhi
*rak:      integer id, dxyz, datom, icount
*rak:c
*rak:      call dfill(nbf*nbf*3*nat,0.0,sa,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      if(.not.bas_nbf_cn_max(basis,nbfmax))
*rak:     &       call errquit('raktest:',922)
*rak:c
*rak:      call hf_print_set(1)
*rak:      used = nbfmax*nbfmax*3*2
*rak:      iscr = used + 1
*rak:      newnscr = nscr - used
*rak:
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh = 1,ish
*rak:
*rak:          call intd_1eov(basis,ish,basis,jsh,newnscr,scr(iscr),
*rak:     &           used,scr,idatom)
*rak:
*rak:          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &           call errquit('error',911)
*rak:          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &           call errquit('error',911)
*rak:          
*rak:          icount = 0
*rak:          do 00300 id = 1,2
*rak:            if (idatom(id).gt.0) then
*rak:              do 00400 dxyz = 1,3
*rak:                do 00500 ibf = ibflo,ibfhi
*rak:                  do 00600 jbf = jbflo,jbfhi
*rak:                    icount = icount + 1
*rak:                    sa(ibf,jbf,dxyz,idatom(id)) = scr(icount)
*rak:                    sa(jbf,ibf,dxyz,idatom(id)) = scr(icount)
*rak:00600             continue
*rak:00500           continue
*rak:00400         continue
*rak:            else
*rak:              icount = icount + 3*(ibfhi-ibflo+1)*(jbfhi-jbflo+1)
*rak:            endif
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:
*rak:      if (.not.printit) return
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' SaNEW <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(sa(1,1,dxyz,datom),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' SaNEW ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      call hf_print_set(0)
*rak:      end
*rak:c.......................................................................
*rak:      subroutine buildHa(Ha,nat,ncont,nbf,basis,
*rak:     &       scr,nscr,printit)
*rak:      implicit none
*rak:c
*rak:#include "bas.fh"
*rak:c
*rak:      integer nat,ncont,nbf,basis,nscr
*rak:      double precision ha(nbf,nbf,3,nat),scr(nscr)
*rak:      logical printit
*rak:c
*rak:      integer used, iscr, newnscr, nbfmax
*rak:      integer ish, ibf, ibflo, ibfhi
*rak:      integer jsh, jbf, jbflo, jbfhi
*rak:      integer dxyz, datom, icount
*rak:c
*rak:      call dfill(nbf*nbf*3*nat,0.0,ha,1)
*rak:      call dfill(nscr,0.0,scr,1)
*rak:c
*rak:      if(.not.bas_nbf_cn_max(basis,nbfmax))
*rak:     &       call errquit('raktest:',922)
*rak:c
*rak:      used = nbfmax*nbfmax*3*nat
*rak:      iscr = used + 1
*rak:      newnscr = nscr - used
*rak:c
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh = 1,ish
*rak:c
*rak:          call intd_1eh1(basis,ish,basis,jsh,newnscr,scr(iscr),
*rak:     &           used,scr)
*rak:c
*rak:          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &           call errquit('error',911)
*rak:          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &           call errquit('error',911)
*rak:c          
*rak:          icount = 0
*rak:          do 00300 datom = 1,nat
*rak:              do 00400 dxyz = 1,3
*rak:                do 00500 ibf = ibflo,ibfhi
*rak:                  do 00600 jbf = jbflo,jbfhi
*rak:                    icount = icount + 1
*rak:                    ha(ibf,jbf,dxyz,datom) = scr(icount)
*rak:                    ha(jbf,ibf,dxyz,datom) = scr(icount)
*rak:00600             continue
*rak:00500           continue
*rak:00400         continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:
*rak:      if (.not.printit) return
*rak:      do 00900 datom = 1,nat
*rak:        do 01000 dxyz = 1,3
*rak:          write(6,'(//,a,i2,a,i2,a)')
*rak:     &           ' HaNEW <atom:',datom,'>  <xyz:',dxyz,'>'
*rak:          call myoutput(ha(1,1,dxyz,datom),1,nbf,1,nbf,nbf,nbf,1,
*rak:     &           ' HaNEW ',1.0d-05)
*rak:01000   continue
*rak:00900 continue
*rak:      end
*rak:c......................................................................
*rak:      subroutine build_d2e_new(basis,nbf,ncont,nat,intsize,d2e,
*rak:     &       maxg,gbuf,nscr,scr)
*rak:      implicit none
*rak:#include "bas.fh"
*rak:c::passed
*rak:      integer basis,nbf,nat,intsize,maxg,nscr,ncont
*rak:      double precision d2e(intsize,3,nat)
*rak:      double precision gbuf(maxg)
*rak:      double precision scr(nscr)
*rak:c::local
*rak:      integer idatom(4)
*rak:      integer dxyz,id,index,icount
*rak:      integer ish,ibf,ibflo,ibfhi
*rak:      integer jsh,jbf,jbflo,jbfhi
*rak:      integer ksh,kbf,kbflo,kbfhi
*rak:      integer lsh,lbf,lbflo,lbfhi
*rak:c
*rak:C     SET-UP FUNCTIONS
*rak:c
*rak:c...  statement function
*rak:C
*rak:C     SET-UP FUNCTIONS
*rak:      integer min,max
*rak:      intrinsic min,max
*rak:      integer i,j,k,l,isym2,isym4
*rak:      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
*rak:      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
*rak:     &               (max(isym2(i,j),isym2(k,l))-1)/2+
*rak:     &               min(isym2(i,j),isym2(k,l))
*rak:c
*rak:c      
*rak:      do 00100 ish = 1,ncont
*rak:        do 00200 jsh= 1,ish
*rak:          do 00300 ksh = 1,jsh
*rak:            do 00400 lsh = 1,ksh
*rak:c
*rak:              call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
*rak:     &               nscr,scr,maxg,gbuf,idatom)
*rak:c
*rak:          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
*rak:     &           call errquit('error',911)
*rak:          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
*rak:     &           call errquit('error',911)
*rak:          if(.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
*rak:     &           call errquit('error',911)
*rak:          if(.not.bas_cn2bfr(basis,lsh,lbflo,lbfhi))
*rak:     &           call errquit('error',911)
*rak:c
*rak:c ... fill integral buffers              
*rak:                icount = 0
*rak:                do 00600 id = 1,4
*rak:                  if (idatom(id).gt.0) then
*rak:                    do 00700 dxyz = 1,3
*rak:                      do 00800 ibf = ibflo,ibfhi
*rak:                        do 00900 jbf = jbflo,jbfhi
*rak:                          do 01000 kbf = kbflo,kbfhi
*rak:                            do 01100 lbf = lbflo,lbfhi
*rak:                              index = isym4(ibf,jbf,kbf,lbf)
*rak:                              icount = icount + 1
*rak:                              d2e(index,dxyz,idatom(id)) = gbuf(icount)
*rak:01100                       continue
*rak:01000                     continue
*rak:00900                   continue
*rak:00800                 continue
*rak:00700               continue
*rak:                  else
*rak:                    icount = icount +
*rak:     &                     3*(ibfhi-ibflo+1)*(jbfhi-jbflo+1)*
*rak:     &                     (kbfhi-kbflo+1)*(lbfhi-lbflo+1)
*rak:                  endif
*rak:00600           continue
*rak:c
*rak:00400       continue
*rak:00300     continue
*rak:00200   continue
*rak:00100 continue
*rak:      end
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
c put scf:gradients = 3*nat reals
        if (.not. rtdb_put(rtdb, 'scf:gradients', MT_DBL, 
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
#if defined(INT_TIMING)
#include "int_tim.fh"
#include "pstat.fh"
#endif
c
c Compute the overlap, kinetic energy, and nuclear attraction integrals for 
c two shells of contracted Gaussians functions. This driver is NOT capable of 
c evaluating integral derivatives.
c
c******************************************************************************
#if defined(INT_TIMING)
      if((.not.dryrun).and.o_int_time.and.(int_time_level.ge.1))
     &   call pstat_on(intt_arr(intt_hf1))
#endif
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
#if defined(INT_TIMING)
      if((.not.dryrun).and.o_int_time.and.(int_time_level.ge.1))
     &   call pstat_off(intt_arr(intt_hf1))
#endif
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
      integer nbf_seg, ncont_seg, nbf_ns, ncont_ns
      integer icount, niter
      double precision norm, tov_seg, tov_ns
      logical status
      logical FF, FT
      parameter (FF=.false.,FT=.true.)
c
      logical int_normalize
      external int_normalize
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
c
c set niter number of iterations
c
      niter = 10
      write(6,'(/,a,i5)')
     &    'number of instances computed for each seg/noseg set:',niter
c
c zero ma segments
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns))
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg))
      call dfill(memscr,           0.0d00,dbl_mb(k_scr))
      call dfill(membuf,           0.0d00,dbl_mb(k_buf))
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
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_s_ns),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_s_ns),1)
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
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns))
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg))
      call dfill(memscr,           0.0d00,dbl_mb(k_scr))
      call dfill(membuf,           0.0d00,dbl_mb(k_buf))
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
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_s_ns),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_s_ns),1)
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
      call dfill((nbf_ns*nbf_ns),  0.0d00,dbl_mb(k_s_ns))
      call dfill((nbf_seg*nbf_seg),0.0d00,dbl_mb(k_s_seg))
      call dfill(memscr,           0.0d00,dbl_mb(k_scr))
      call dfill(membuf,           0.0d00,dbl_mb(k_buf))
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
      call daxpy((nbf_seg*nbf_seg),(-1.0d00),
     &      dbl_mb(k_s_seg),1,dbl_mb(k_s_ns),1)
      norm = ddot((nbf_ns*nbf_ns),dbl_mb(k_s_ns),1,dbl_mb(k_s_ns),1)
      write(luout,*)' '
      write(luout,*)' time for segmented potential         :',tov_seg
      write(luout,*)' time for non-segmented potential     :',tov_ns
      write(luout,'(a,f10.2)')
     &    ' % speedup                          :',
     &    (tov_seg-tov_ns)/tov_seg*100.0d00
      write(luout,*)'raktest_gc: potential difference norm :',norm
      write(luout,*)' '
c
      status =              ma_pop_stack(h_s_ns)
      status = status .and. ma_pop_stack(h_s_seg)
      status = status .and. ma_pop_stack(h_scr)
      status = status .and. ma_pop_stack(h_buf)
      if (.not.status) call errquit('raktest_gc: pop error',911)
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
*              if (abs(val).gt.1.0d-07) then
*                icount = icount + 1
*                write(69,10000)ibf,jbf,val,msg,' ints ',icount
*                write(luout,10000)ibf,jbf,val,msg,' ints ',icount
*              endif
            enddo
          enddo
        enddo
      enddo
c
*      write(luout,*)' raktest_gc overlap matrix ',msg
*      call output(S,1,nbf,1,nbf,nbf,nbf,1)
c
10000 format(i5,i5,1pd20.10,1x, a,a,i12)
c
      end

