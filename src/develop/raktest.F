      subroutine raktest(rtdb)
      implicit none
c $Id: raktest.F,v 1.33 1995-02-07 00:48:29 d3e129 Exp $
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
      goto (00100,00200,00300,00400) raktask
c
      write(6,*)' default raktest task '
      call raktest0(rtdb)
      goto 00001
c
00100 write(6,*)' raktest task 1 gen_hess at'
      call raktest1(rtdb)
      goto 00001
c
00200 write(6,*)' raktest task 2 gen_hess opt/at'
      call raktest2(rtdb)
      goto 00001
c
00300 write(6,*)' raktest task 3 stepper test'
      call raktest3(rtdb)
      goto 00001
00400 write(6,*)' raktest task 4 check intd_init'
      call raktest4(rtdb)
      goto 00001
c
00001 continue
      end
      subroutine raktest0(rtdb)
      implicit none
c test0 "default raktest"
#include "global.fh"
      integer rtdb
c::functions
      logical gen_apphess_at
      external gen_apphess_at
c
      double precision delta
c
#if defined(SUN)
      real etime,tarray(2)
      external etime
#endif
      real t1
      t1 = 0.0
c
c      logical status
c
      if(ga_nodeid().eq.0) then
        write(6,*) ' rak test'
      endif
      delta = 0.01d00
#if defined(SUN)
      t1 = etime(tarray)
#endif
      if(.not.gen_apphess_at(rtdb,delta)) call errquit
     &    ('raktest:gen_apphess_at failed',911)
#if defined(SUN)
      t1 = etime(tarray) - t1 
#endif
      if(ga_nodeid().eq.0) then
        write(6,'(1x,a,f8.2)')
     &      'time for rak test',t1
      endif
      end
      subroutine raktest2(rtdb)
      implicit none
c test0 "default raktest"
#include "global.fh"
      integer rtdb
c::functions
      logical gen_apphess_opt
      external gen_apphess_opt
c
      double precision delta
c
#if defined(SUN)
      real etime,tarray(2)
      external etime
#endif
      real t1
      t1 = 0.0
c
c      logical status
c
      if(ga_nodeid().eq.0) then
        write(6,*) ' rak test'
      endif
      delta = 0.01d00
#if defined(SUN)
      t1 = etime(tarray)
#endif
      if(.not.gen_apphess_opt(rtdb,delta)) call errquit
     &    ('raktest:2:gen_apphess_opt failed',911)
#if defined(SUN)
      t1 = etime(tarray) - t1 
#endif
      if(ga_nodeid().eq.0) then
        write(6,'(1x,a,f8.2)')
     &      'time for rak test 2',t1
      endif
      end
      subroutine raktest1(rtdb)
      implicit none
      integer rtdb
      call deriv3c(rtdb)
      call deriv2c(rtdb)
      end
      subroutine deriv3c(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "bas.fh"
c
      double precision ddot
      external ddot
c
      integer rtdb
c
      character*100 mo_basis, geom_name
      integer basis, geom
      integer nbf, nshells, nat
      integer size2e, size2e_d
      integer h2e3d,i2e3d, h2e3d_fd,i2e3d_fd
      integer hxyz,ixyz,nbuf,nscr,hscr,iscr,hbuf,ibuf
      logical status
c
      if (.not.context_rtdb_match(rtdb,'mo basis',mo_basis))
     &    mo_basis = 'mo basis'
      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
     &    geom_name = 'geometry'
c
      if(.not.geom_create(geom,geom_name))call errquit
     &    ('deriv3c: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('deriv3c: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('deriv3c: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('deriv3c: basis load ',911)
c
      if(.not.bas_numbf(basis,nbf)) call errquit
     &    ('deriv3c: bas_numbf failed',911)
      if(.not.bas_numcont(basis,nshells)) call errquit
     &       ('deriv3c:bas_numcont failed',911)
      if(.not.geom_ncent(geom,nat)) call errquit
     &       ('deriv3c:geom_ncent failed',911)
c
      size2e   = nbf*nbf*nbf
      size2e_d = size2e*3*nat
      write(6,*)' size2e   = ',size2e
      write(6,*)' size2e_d = ',size2e_d
c
      call intd_init(rtdb,1,basis)
*      call int_normalize(basis)
      call int_norm_2c(basis)
c
      if(.not.ma_set_error_print(.true.)) stop ' fuck 1'
      if(.not.ma_verify_allocator_stuff()) stop ' fuck 2'
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
      if (.not.ma_push_get(mt_dbl,size2e_d,'2e3c d real',
     &    h2e3d,i2e3d))
     &    call errquit('deriv3c: ma fail 1',911)
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
      if (.not.ma_push_get(mt_dbl,size2e_d,'2e3c d fd',
     &    h2e3d_fd,i2e3d_fd))
     &    call errquit('deriv3c: ma fail 2',911)
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
      if (.not.ma_push_get(mt_dbl,(3*nat),'xyzbuf',
     &    hxyz,ixyz))
     &    call errquit('deriv3c: ma fail 4',911)
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
      call int_mem_2e4c(nbuf,nscr)
      if (.not.ma_push_get(mt_dbl,nscr,'scratch',
     &    hscr,iscr))
     &    call errquit('deriv3c: ma fail 5',911)
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
      if (.not.ma_push_get(mt_dbl,(2*nbuf),'buffer scratch',
     &    hbuf,ibuf))
     &    call errquit('deriv3c: ma fail 6',911)
      write(6,*)' stack available ',ma_inquire_stack(mt_dbl)
      write(6,*)' heap  available ',ma_inquire_heap(mt_dbl)
c
      call dfill(size2e_d, 0.0d00,dbl_mb(i2e3d)     ,1)
      call dfill(size2e_d, 0.0d00,dbl_mb(i2e3d_fd)  ,1)
      call dfill(nscr,     0.0d00,dbl_mb(iscr)      ,1)
      call dfill(2*nbuf,   0.0d00,dbl_mb(ibuf)      ,1)
      call dfill((3*nat),  0.0d00,dbl_mb(ixyz)      ,1)
c      
c      
      call build_3ed(dbl_mb(i2e3d),nat,size2e,nbf,nshells,
     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,basis,geom)
      call dfill(nscr,      0.0d00,dbl_mb(iscr)      ,1)
      call dfill(2*nbuf,    0.0d00,dbl_mb(ibuf)      ,1)
      call build_3edfd(dbl_mb(i2e3d_fd),nat,size2e,nbf,
     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,dbl_mb(ixyz),
     &    1.0d-05,basis,geom,nshells)
      call daxpy(size2e_d,-1.0d00,dbl_mb(i2e3d),1,dbl_mb(i2e3d_fd),1)
c
      call print_2d3c(dbl_mb(i2e3d),dbl_mb(i2e3d_fd),nbf,nat)
c
      write(6,*)
      write(6,*)' 2e3c deriv diff norm ',
     &    ddot(size2e_d,dbl_mb(i2e3d_fd),1,dbl_mb(i2e3d_fd),1)
c
      call intd_terminate()
      status = ma_pop_stack(hbuf)
      status = status.and.ma_pop_stack(hscr)
      status = status.and.ma_pop_stack(hxyz)
      status = status.and.ma_pop_stack(h2e3d_fd)
      status = status.and.ma_pop_stack(h2e3d)
      if (.not.status) call errquit('pop failed',911)
      status = bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      if (.not.status) call errquit('b/g destroy failed',911)
      end
      subroutine build_3edfd(x2e3dfd,nat,size2e,nbf,scr,nscr,
     &    buf,nbuf,xyz,delta,basis,geom,nshells)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "geomP.fh"
      integer nat, size2e, nbf, nscr, nbuf, basis, geom, nshells
      double precision x2e3dfd(nbf,nbf,nbf,3,nat)
      double precision scr(nscr)
      double precision buf(nbuf,2)
      double precision xyz(3,nat)
      double precision delta
c
      double precision scale
      integer ish, ibflo, ibfhi, ibf
      integer jsh, jbflo, jbfhi, jbf
      integer ksh, kbflo, kbfhi, kbf
      integer datom, dxyz
      integer nint, icount
      integer int_nint
      external int_nint
c
c.. store original coords in xyz
      call dcopy((3*nat),coords(1,1,geom),1,xyz,1)
      scale = 1.0d00/(2.0d00*delta)
c
      do 00100 ish = 1,nshells
        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &      call errquit('cn2bfr error',ish)
        do 00200 jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &        call errquit('cn2bfr error',jsh)
          do 00300 ksh = 1,jsh
            if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
     &          call errquit('cn2bfr error',ksh)
            nint = int_nint(basis,ish,basis,jsh,basis,ksh,basis,0)
            do 00400 datom = 1,nat
              do 00500 dxyz = 1,3
c... copy original coords in to coords array
                call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
c... add delta to coordinate
                coords(dxyz,datom,geom)=coords(dxyz,datom,geom)+delta
                call int_2e3c(basis,ish,basis,jsh,ksh,
     &              nscr,scr,nint,buf(1,1))
c... copy original coords in to coords array
                call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
c... subtract delta to coordinate
                coords(dxyz,datom,geom)=coords(dxyz,datom,geom)-delta
                call int_2e3c(basis,ish,basis,jsh,ksh,
     &              nscr,scr,nint,buf(1,2))
                call daxpy(nint,-1.0d00,buf(1,2),1,buf(1,1),1)
                call dscal(nint,scale,buf(1,1),1)
                icount = 0
                do 00600 ibf = ibflo,ibfhi
                  do 00700 jbf = jbflo,jbfhi
                    do 00800 kbf = kbflo, kbfhi
                      icount = icount + 1
                      x2e3dfd(ibf,jbf,kbf,dxyz,datom) =
     &                    buf(icount,1)
00800               continue
00700             continue
00600           continue
00500         continue
00400       continue
00300     continue
00200   continue
00100 continue
c      
c... copy original coords in to coords array
      call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
c
      end
      subroutine build_3ed(x2e3d,nat,size2e,nbf,nshells,
     &    scr,nscr,buf,nbuf,basis,geom)
      implicit none
#include "bas.fh"
      integer nat, size2e, nbf, nshells, nscr, nbuf, basis, geom
      double precision x2e3d(nbf,nbf,nbf,3,nat)
      double precision scr(nscr)
      double precision buf(nbuf)
c
      integer idatom(4)
      integer ish, ibflo, ibfhi, num_ibf, ibf
      integer jsh, jbflo, jbfhi, num_jbf, jbf
      integer ksh, kbflo, kbfhi, num_kbf, kbf
      integer nint, icount, id, dxyz
c
      integer int_nint
      external int_nint
c
      do 00100 ish = 1,nshells
        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &      call errquit('cn2bfr error',ish)
        num_ibf = ibfhi-ibflo+1
        do 00200 jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &        call errquit('cn2bfr error',jsh)
          num_jbf = jbfhi-jbflo+1
          do 00300 ksh = 1,jsh
            if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
     &          call errquit('cn2bfr error',ksh)
            num_kbf = kbfhi-kbflo+1
            nint = int_nint(basis,ish,basis,jsh,basis,ksh,basis,0)
            call intd_2e3c(basis,ish,basis,jsh,ksh,nscr,scr,(nint*3*4),
     &          buf,idatom)
            write(6,*)' idatom 2e3c',idatom
            icount = 0
            do 00400 id = 1,4
              if (idatom(id).gt.0) then
                do 00500 dxyz = 1,3
                  do 00600 ibf = ibflo,ibfhi
                    do 00700 jbf = jbflo,jbfhi
                      do 00800 kbf = kbflo,kbfhi
                        icount = icount + 1
                        x2e3d(ibf,jbf,kbf,dxyz,idatom(id))=buf(icount)
00800                 continue
00700               continue
00600             continue
00500           continue
              else
                icount = icount + 3*num_ibf*num_jbf*num_kbf
              endif
00400       continue
00300     continue
00200   continue
00100 continue
c
      end
      subroutine deriv2c(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "bas.fh"
c
      double precision ddot
      external ddot
c
      integer rtdb
c
      character*100 mo_basis, geom_name
      integer basis, geom
      integer nbf, nshells, nat
      integer size2e, size2e_d
      integer h2e2d,i2e2d, h2e2d_fd,i2e2d_fd, h2e2d_diff,i2e2d_diff
      integer hxyz,ixyz,nbuf,nscr,hscr,iscr,hbuf,ibuf
      logical status
c
      if (.not.context_rtdb_match(rtdb,'mo basis',mo_basis))
     &    mo_basis = 'mo basis'
      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
     &    geom_name = 'geometry'
c
      if(.not.geom_create(geom,geom_name))call errquit
     &    ('deriv2c: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('deriv2c: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('deriv2c: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('deriv2c: basis load ',911)
c
      if(.not.bas_numbf(basis,nbf)) call errquit
     &    ('deriv2c: bas_numbf failed',911)
      if(.not.bas_numcont(basis,nshells)) call errquit
     &       ('deriv2c:bas_numcont failed',911)
      if(.not.geom_ncent(geom,nat)) call errquit
     &       ('deriv2c:geom_ncent failed',911)
c
      size2e   = nbf*(nbf+1)/2
      size2e_d = size2e*3*nat
c
      call intd_init(rtdb,1,basis)
*      call int_normalize(basis)
      call int_norm_2c(basis)
c
      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d real',
     &    h2e2d,i2e2d))
     &    call errquit('deriv2c: ma fail 1',911)
      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d fd',
     &    h2e2d_fd,i2e2d_fd))
     &    call errquit('deriv2c: ma fail 2',911)
      if (.not.ma_push_get(mt_dbl,size2e_d,'2e2c d diff',
     &    h2e2d_diff,i2e2d_diff))
     &    call errquit('deriv2c: ma fail 3',911)
      if (.not.ma_push_get(mt_dbl,(3*nat),'xyzbuf',
     &    hxyz,ixyz))
     &    call errquit('deriv2c: ma fail 4',911)
      call int_mem_2e4c(nbuf,nscr)
      if (.not.ma_push_get(mt_dbl,nscr,'scratch',
     &    hscr,iscr))
     &    call errquit('deriv2c: ma fail 5',911)
      if (.not.ma_push_get(mt_dbl,(2*nbuf),'buffer scratch',
     &    hbuf,ibuf))
     &    call errquit('deriv2c: ma fail 6',911)
c
      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d)     ,1)
      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d_fd)  ,1)
      call dfill(size2e_d, 0.0d00,dbl_mb(i2e2d_diff),1)
      call dfill(nscr,     0.0d00,dbl_mb(iscr)      ,1)
      call dfill(2*nbuf,   0.0d00,dbl_mb(ibuf)      ,1)
      call dfill((3*nat),  0.0d00,dbl_mb(ixyz)      ,1)
c      
      call build_2ed(dbl_mb(i2e2d),nat,size2e,nbf,nshells,
     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,basis,geom)
      call dfill(nscr,      0.0d00,dbl_mb(iscr)      ,1)
      call dfill(2*nbuf,    0.0d00,dbl_mb(ibuf)      ,1)
      call build_2edfd(dbl_mb(i2e2d_fd),nat,size2e,nbf,
     &    dbl_mb(iscr),nscr,dbl_mb(ibuf),nbuf,dbl_mb(ixyz),
     &    1.0d-05,basis,geom,nshells)
      call dcopy(size2e_d,dbl_mb(i2e2d_fd),1,dbl_mb(i2e2d_diff),1)
      call daxpy(size2e_d,-1.0d00,dbl_mb(i2e2d),1,dbl_mb(i2e2d_diff),1)
c
      call print_2d2c(dbl_mb(i2e2d),dbl_mb(i2e2d_fd),
     &    dbl_mb(i2e2d_diff),size2e,nat)
c
      write(6,*)
      write(6,*)' 2e2c deriv diff norm ',
     &    ddot(size2e_d,dbl_mb(i2e2d_diff),1,dbl_mb(i2e2d_diff),1)
c
      call intd_terminate()
      status = ma_pop_stack(hbuf)
      status = status.and.ma_pop_stack(hscr)
      status = status.and.ma_pop_stack(hxyz)
      status = status.and.ma_pop_stack(h2e2d_diff)
      status = status.and.ma_pop_stack(h2e2d_fd)
      status = status.and.ma_pop_stack(h2e2d)
      if (.not.status) call errquit('pop failed',911)
      status = bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      if (.not.status) call errquit('b/g destroy failed',911)
      end
      subroutine build_2ed(x2e2d,nat,size2e,nbf,nshells,
     &    scr,nscr,buf,nbuf, basis, geom)
      implicit none
#include "bas.fh"
      integer basis, geom
      integer nat
      integer size2e
      integer nbf
      integer nshells
      double precision x2e2d(size2e,3,nat)
      integer nscr, nbuf, idatom(2), id, dxyz, index
      double precision scr(nscr),buf(nbuf)      
c
      integer ish, ibflo, ibfhi, num_ibf, ibf
      integer jsh, jbflo, jbfhi, num_jbf, jbf
      integer nint, icount
      integer int_nint
      external int_nint
      integer min, max
      intrinsic min,max
      integer i,j,isym2
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      do 00100 ish=1,nshells
        do 00200 jsh=1,ish
          nint = int_nint(basis,ish,basis,jsh,basis,0,basis,0)
          call intd_2e2c(basis,ish,basis,jsh,nscr,scr,(nint*3*4),
     &        buf,idatom)
          write(6,*)' idatom 2e2c', idatom
          if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &        call errquit('cn2bfr error',ish)
          num_ibf = ibfhi-ibflo+1
          if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &        call errquit('cn2bfr error',ish)
          num_jbf = jbfhi-jbflo+1
          icount = 0
          do 00300 id = 1,2
            if (idatom(id).gt.0) then
              do 00400 dxyz = 1,3
                do 00500 ibf = ibflo,ibfhi
                  do 00600 jbf = jbflo, jbfhi
                    icount = icount + 1
                    index = isym2(ibf,jbf)
                    x2e2d(index,dxyz,idatom(id)) = buf(icount)
00600             continue
00500           continue
00400         continue
            else
              icount = icount + 3*num_ibf*num_jbf
            endif
00300     continue
00200   continue
00100 continue
      end
      subroutine build_2edfd(x2e2d_fd,nat,size2e,nbf,
     &    scr,nscr,buf,nbuf,xyz,delta,basis,geom,nshells)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "geomP.fh"
c
      integer nat,size2e,nbf, nscr, nbuf, basis, geom, nshells
      double precision x2e2d_fd(size2e,3,nat)
      double precision scr(nscr)
      double precision buf(nbuf,2)
      double precision xyz(3,nat)
      double precision delta, scale
      integer ish, ibf, ibflo, ibfhi
      integer jsh, jbf, jbflo, jbfhi
      integer nint, icount, index, datom, dxyz
      integer int_nint
      external int_nint
      integer min, max
      intrinsic min,max
      integer i,j,isym2
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c
c.. store original coords in xyz
      call dcopy((3*nat),coords(1,1,geom),1,xyz,1)
      scale = 1.0d00/(2.0d00*delta)
c
      do 00100 ish = 1,nshells
        do 00200 jsh = 1,ish
          nint = int_nint(basis,ish,basis,jsh,basis,0,basis,0)
          do 00300 datom = 1,nat
            do 00400 dxyz = 1,3
c... copy original coords in to coords array
              call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
c... add delta to coordinate
              coords(dxyz,datom,geom) = coords(dxyz,datom,geom) + delta
              call int_2e2c(basis,ish,basis,jsh,nscr,scr,nint,buf(1,1))
c... copy original coords in to coords array
              call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
c... subtract delta to coordinate
              coords(dxyz,datom,geom) = coords(dxyz,datom,geom) - delta
              call int_2e2c(basis,ish,basis,jsh,nscr,scr,nint,buf(1,2))
              call daxpy(nint,-1.0d00,buf(1,2),1,buf(1,1),1)
              call dscal(nint,scale,buf(1,1),1)
              icount = 0
              if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &            call errquit('cn2bfr error',ish)
              if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &            call errquit('cn2bfr error',ish)
              do 00500 ibf = ibflo, ibfhi
                do 00600 jbf = jbflo, jbfhi
                  icount = icount+1
                  index  = isym2(ibf,jbf)
                  x2e2d_fd(index,dxyz,datom) = buf(icount,1)
00600           continue
00500         continue
00400       continue
00300     continue
00200   continue
00100 continue
c      
c... copy original coords in to coords array
      call dcopy((3*nat),xyz,1,coords(1,1,geom),1)
      end
      subroutine print_2d3c(x2e3d,diff,nbf,nat)
      implicit none
      integer nbf,nat
      double precision x2e3d(nbf,nbf,nbf,3,nat)
      double precision diff(nbf,nbf,nbf,3,nat)
c
      double precision thresh
      integer datom,dxyz, ibf,jbf,kbf
      logical printit
c
      thresh = 1.0d-08
c
      do 00100 datom = 1, nat
        do 00200 dxyz = 1,3
          write(6,*)
          write(6,*)' atom = ',datom,' xyz = ',dxyz
          write(6,*)'--------------------------------------------------'
          do 00300 ibf = 1,nbf
            do 00400 jbf = 1,nbf
              do 00500 kbf = 1,nbf
                printit = abs(x2e3d(ibf,jbf,kbf,dxyz,datom)).gt.thresh
                if (printit) then
                  write(6,10000)ibf,jbf,kbf,datom,dxyz,
     &                x2e3d(ibf,jbf,kbf,dxyz,datom),
     &                diff(ibf,jbf,kbf,dxyz,datom)
                endif
00500         continue
00400       continue
00300     continue
00200   continue
00100 continue
c
10000 format('d2e3c <',i3,',',i3,',',i3,'|',i3,'|',i3,' > val =',
     &    d20.10,' diff =',d20.10)
c
      end
      subroutine print_2d2c(x2e2d,x2e2d_fd,diff,size2e,nat)
      implicit none
      integer size2e,nat
      double precision x2e2d(size2e,3,nat)
      double precision x2e2d_fd(size2e,3,nat)
      double precision diff(size2e,3,nat)
c
      double precision thresh
      integer datom,dxyz, index
      logical printit
c
      thresh = 1.0d-08
c
      do 00100 datom = 1, nat
        do 00200 dxyz = 1,3
          write(6,*)
          write(6,*)' atom = ',datom,' xyz = ',dxyz
          write(6,*)'--------------------------------------------------'
          do 00300 index = 1,size2e
            printit = abs(x2e2d(index,dxyz,datom)).gt.thresh
            printit = printit .or.
     &          abs(x2e2d_fd(index,dxyz,datom)).gt.thresh
            if (printit) then
              write(6,10000)index,datom,dxyz,
     &            x2e2d(index,dxyz,datom),
     &            x2e2d_fd(index,dxyz,datom),
     &            diff(index,dxyz,datom)
            endif
00300     continue
00200   continue
00100 continue
c
10000 format('d2e2c <',i5,'|',i3,'|',i3,' > val =',d20.10,
     &    ' fd =',d20.10,' diff =',d20.10)
c
      end
      logical function gen_apphess_at(rtdb,delta)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "context.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
c::functions
c::passed
      integer rtdb
      double precision delta
c::local
      integer whoiam, master
      integer geom
      integer nat, rank_hess
      integer basis
      integer h_grad0, k_grad0 ! handle/index for central gradient
      integer h_gradp, k_gradp ! handle/index for delta gradient
      integer h_hess, k_hess ! handle/index for hessian matrix
      integer iatom, ixyz
      logical ignore_status
      double precision xyz(3), chg
      character*5 tmpname
      character*16 tag_name
      character*40 new_geom_name
      character*255 mo_basis, fd_basis, mo_basis_ref
c
      whoiam = ga_nodeid()
      master = 0
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
c
c.. create/load reference geometry
        if (.not.geom_create(geom,'geometry')) call errquit
     &      ('gen_apphess_at:geom_create failed?',1)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &      ('gen_apphess_at:geom_rtdb_load failed?',2)
c..  get the number of atoms
        if (.not. geom_ncent(geom,nat)) call errquit
     &      ('gen_apphess_at:geom_ncent failed?',3)
c
c.. copy reference geometry and store on rtdb
        if (.not.
     &      geom_rtdb_store(rtdb,geom,'reference'))
     &      call errquit('gen_apphess_at: geom_rtdb_store failed',5)
c
c.. determine what fd basis name is
c
        if (.not.context_rtdb_match(rtdb,'mo basis',mo_basis))
     &      mo_basis = 'mo basis'
        mo_basis_ref = mo_basis
        if (.not.context_rtdb_match(rtdb,'fd basis',fd_basis))
     &      fd_basis = 'fd basis'
c
c.. try to load "fd basis"
c.. create basis handle
        if(.not.bas_create(basis,fd_basis)) call errquit
     &      ('gen_apphess_at: bas_create failed?',7)
        if(bas_rtdb_load(rtdb,geom,basis,fd_basis)) then
          continue
        else
c...  if fd basis not defined use current mo_basis
          fd_basis = mo_basis
        endif
c
c set "mo basis" to point to "fd basis"
        if (.not.rtdb_cput(rtdb,'mo basis',1,fd_basis)) call errquit
     &      ('gen_apphess_at: rtdb_cpu failed?',8)
c
        if (.not.(bas_destroy(basis))) call errquit
     &      ('gen_apphess_at: bas_destroy failed?',9)
        if (.not.(geom_destroy(geom))) call errquit
     &      ('gen_apphess_at: geom_destroy failed?',9)
c set "geometry" to "reference"
        new_geom_name='reference'
        if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &      call errquit
     &        ('gen_apphess_at: rtdb_cput failed ?',911)
      endif
c all nodes work here 
c
c..      broadcast number of atoms
c
      call ga_brdcst(Msg_RAK+MSGINT, nat, mitob(1), 0)
c
      call ga_sync()
c
      ignore_status = rtdb_parallel(.true.)
c      
      call scf_gradient(rtdb)
c
      ignore_status = rtdb_parallel(.false.)
c      
      call ga_sync()
      if (whoiam.eq.master) then
        rank_hess = 3*nat
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'cent nuclear gradient vector',
     &      h_grad0,k_grad0)) call errquit
     &      ('gen_apphess_at: ma_push_get failed?',10)
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'delta nuclear gradient vector',
     &      h_gradp,k_gradp)) call errquit
     &      ('gen_apphess_at: ma_push_get failed?',11)
        if (.not. MA_Push_Get(MT_DBL,(rank_hess*rank_hess),
     &      'nuclear hessian matrix',
     &      h_hess,k_hess)) call errquit
     &      ('gen_apphess_at: ma_push_get failed?',12)
c
        call dfill((rank_hess*rank_hess),0.0d00,
     &      dbl_mb(k_hess),1)
c
c.. get central gradient
        if(.not. rtdb_get(rtdb,'scf:gradients',MT_DBL,
     &      rank_hess,dbl_mb(k_grad0)))
     &      call errquit('gen_apphess_at: rtdb_get failed',13)
c
      endif
c
c... all nodes doing work here
      call ga_sync()
      call util_flush(6)
      do 00100 iatom = 1,nat
        do 00200 ixyz = 1,3
          if(whoiam.eq.master) then
c... master node only
            if (.not.geom_create(geom,'reference')) call errquit
     &          ('gen_apphess_at:geom_create failed?',14)
            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &          call errquit
     &          ('gen_apphess_at:geom_rtdb_load failed?',15)
            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('gen_apphess_at:geom_cent_get failed?',16)
            xyz(ixyz) = xyz(ixyz) + delta
            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('gen_apphess_at:geom_cent_get failed?',17)
c.. copy modified geometry and store on rtdb
            write(tmpname,'(i3,1x,i1)')iatom,ixyz
            new_geom_name = 'reference '//tmpname
            if (.not.
     &          geom_rtdb_store(rtdb,geom,new_geom_name))
     &          call errquit
     &          ('gen_apphess_at: geom_rtdb_store failed',19)
            if (.not.geom_destroy(geom))
     &          call errquit('gen_apphess_at: geom_destroy failed?',20)
c set "geometry" to "reference atom/xyz"
            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &          call errquit
     &          ('gen_apphess_at: rtdb_cput failed ?',911)
          endif
          call ga_sync()
          ignore_status = rtdb_parallel(.true.)
          call scf_gradient(rtdb)
          ignore_status = rtdb_parallel(.false.)
          call ga_sync()
          if (whoiam.eq.master) then
            if(.not. rtdb_get(rtdb,'scf:gradients',MT_DBL,
     &          rank_hess,dbl_mb(k_gradp)))
     &          call errquit('gen_apphess_at: rtdb_get failed',21)
            call fd_upd_hess(dbl_mb(k_hess),dbl_mb(k_grad0),
     &          dbl_mb(k_gradp),1.0d00,delta,nat,iatom,ixyz)
          endif
00200   continue
00100 continue
c
      if (whoiam.eq.master) then
        write(6,*)' finite difference hessian delta = ',delta
        call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
        call wrt_fd_from_sq(dbl_mb(k_hess),rank_hess,'hess_ascii')
        write(6,*)' triangle hessian written to hess_ascii'
        
c.. redefine mo basis to it's original setting
        if (.not.rtdb_cput(rtdb,'mo basis',1,mo_basis_ref))
     &      call errquit
     &      ('gen_apphess_at: rtdb_cpu failed?',22)
c set "geometry" to "geometry"
            if(.not.rtdb_cput(rtdb,'geometry',1,'geometry'))
     &          call errquit
     &          ('gen_apphess_at: rtdb_cput failed ?',911)
c
        if(.not.ma_pop_stack(h_hess)) call errquit
     &      ('gen_apphess_at: ma_pop_stack(h_hess) failed?',27)
        if(.not.ma_pop_stack(h_gradp)) call errquit
     &      ('gen_apphess_at: ma_pop_stack(h_gradp) failed?',27)
        if(.not.ma_pop_stack(h_grad0)) call errquit
     &      ('gen_apphess_at: ma_pop_stack(h_grad0) failed?',27)

      endif
      call ga_sync()
      ignore_status = rtdb_parallel(.true.)
      gen_apphess_at = .true.
      end
      subroutine wrt_fd_from_sq(hess,rank_hess,filename)
      implicit none
c
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
      character*(*) filename
c
      integer i, j, lu
c
      lu = 69
      open(unit=lu,file=filename,
     &    form='formatted',
     &    access='sequential',
     &    status='unknown')
c
      do 00100 i = 1,rank_hess
        do 00200 j = 1,i
          write(lu,10000)hess(i,j)
00200   continue
00100 continue
c
10000 format(1x,1pd20.10)
c
      end
      subroutine fd_upd_hess(hess,gradm,gradp,s_delta,delta,nat,
     &    iatom_t,ixyz_t)
      implicit none
c::passed
      integer nat
      integer iatom_t, ixyz_t
      double precision hess(3,nat,3,nat)
      double precision gradm(3,nat)
      double precision gradp(3,nat)
      double precision delta, s_delta
c::local
      integer iatom,ixyz
      double precision rdelta, value
c
c 
c finite difference  [g(x+delta) - g(x-delta)]/(s_delta*delta) (s_delta = 2.0)
c central difference [g(x+delta) - g(x)]/(s_delta*delta) (s_delta = 1.0)
c 
c
      rdelta = 1.0d00/(s_delta*delta)
      do 00100 iatom = 1,nat
        do 00200 ixyz = 1,3
          value = rdelta*(gradp(ixyz,iatom)-gradm(ixyz,iatom))
          hess(ixyz_t,iatom_t,ixyz,iatom) = value
          hess(ixyz,iatom,ixyz_t,iatom_t) = value
00200   continue
00100 continue
c
      end
      logical function gen_apphess_opt(rtdb,delta)
      implicit none
c
c function to generate approximate hessian
c
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
      integer rtdb
      double precision delta
c
      integer whoiam, master
      logical ignore_status
      integer geom,nat,basis
      character*40 mo_basis,mo_basis_ref,fd_basis,new_geom_name
      character*40 movecs_in_ref
c
      gen_apphess_opt=.false.
      whoiam = ga_nodeid()
      master = 0
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
c
c.. create/load reference geometry
        if (.not.geom_create(geom,'geometry')) call errquit
     &      ('gen_apphess_at:geom_create failed?',1)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &      ('gen_apphess_at:geom_rtdb_load failed?',2)
c..  get the number of atoms
        if (.not. geom_ncent(geom,nat)) call errquit
     &      ('gen_apphess_at:geom_ncent failed?',3)
c
c.. copy reference geometry and store on rtdb
        if (.not.
     &      geom_rtdb_store(rtdb,geom,'reference'))
     &      call errquit('gen_apphess_at: geom_rtdb_store failed',5)
c
c.. determine what fd basis name is
c
        if (.not.context_rtdb_match(rtdb,'mo basis',mo_basis))
     &      mo_basis = 'mo basis'
        mo_basis_ref = mo_basis
        if (.not.context_rtdb_match(rtdb,'fd basis',fd_basis))
     &      fd_basis = 'fd basis'
c
c.. try to load "fd basis"
c.. create basis handle
        if(.not.bas_create(basis,fd_basis)) call errquit
     &      ('gen_apphess_at: bas_create failed?',7)
        if(bas_rtdb_load(rtdb,geom,basis,fd_basis)) then
          continue
        else
c...  if fd basis not defined use current mo_basis
          fd_basis = mo_basis
          write(6,*) ' fd basis is the mo basis This will be expensive'
          call util_flush(6)
        endif
c
c.. set "mo basis" to point to "fd basis"
        if (.not.rtdb_cput(rtdb,'mo basis',1,fd_basis)) call errquit
     &      ('gen_apphess_at: rtdb_cpu failed?',8)
c... destroty basis and geom handles
        if (.not.(bas_destroy(basis))) call errquit
     &      ('gen_apphess_at: bas_destroy failed?',9)
        if (.not.(geom_destroy(geom))) call errquit
     &      ('gen_apphess_at: geom_destroy failed?',9)
c set "geometry" to "reference"
        new_geom_name='reference'
        if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &      call errquit
     &        ('gen_apphess_at: rtdb_cput failed ?',911)
      endif
c... store movecs initial state
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in_ref))
     $     movecs_in_ref = 'atomic'
c... set initial guess for fd basis to atomic
      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, 'atomic'))
     &     call errquit ('gen_apphess_opt: rtdb_cput failed?',911)
c
      ignore_status = rtdb_parallel(.true.)
      call gen_apphess_at(rtdb,delta)
      call scf_opt(rtdb)
c
c.. reset stepper so it knows to restart the calculation
      call stpr_walk_reset
      ignore_status = rtdb_parallel(.false.)
c.. set "mo basis" to point to %val<mo_basis_ref
      if (.not.rtdb_cput(rtdb,'mo basis',1,mo_basis_ref))
     &    call errquit ('gen_apphess_opt: rtdb_cput failed?',8)
c
c... reset input vectors option here to original guess if appropriate
      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, movecs_in_ref))
     &     call errquit ('gen_apphess_opt: rtdb_cput failed?',911)
c
      call scf_opt(rtdb)
c
      gen_apphess_opt = .true.
      end
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
*rak:c----#include <dimensions.fh>
*rak:c----#include <atoms.fh>
*rak:c----#include <AO_basis_set.fh>
*rak:c----#include <timer.fh>
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
      subroutine raktest3(rtdb)
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
     &         call errquit('raktest3: geom_create?', 911)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &         call errquit('raktest3: no geometry ', 911)
c get number of atoms = nat
        if (.not. geom_ncent(geom,nat))
     &      call errquit('raktest3: geom_ncent?',911)
        if (.not. geom_destroy(geom))
     &      call errquit('raktest3: geom_destroy?',911)
        if (.not.
     &         MA_Push_Get(MT_DBL,(3*nat),'stpr fake gradient',
     &         h_grad,k_grad))
     &         call errquit
     &         ('raktest3: allocation for gradient failed?',911)
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
     &      ('raktest3: failed to read converged in rtdb', 911)
c put scf:energy   =  real value
        if (.not. rtdb_put(rtdb,'scf:energy', MT_DBL, 1, energy))
     &      call errquit
     &      ('raktest3: failed to read energy in rtdb', 911)
c put scf:gradients = 3*nat reals
        if (.not. rtdb_put(rtdb, 'scf:gradients', MT_DBL, 
     &       (3*nat),dbl_mb(k_grad)))
     &      call errquit
     &        ('raktest3: reading gradients failed',911)
c free memory
        if (.not. ma_pop_stack(h_grad))
     &      call errquit('raktest3: pop failed')
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
      subroutine raktest4(rtdb)
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "basP.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer geom
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
      end

