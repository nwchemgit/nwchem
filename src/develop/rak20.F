      logical function task_pderiv(rtdb)
*
* $Id: rak20.F,v 1.1 1998-05-06 00:55:09 d3e129 Exp $
*
***********************************************************************
* This is about the nuclear coordinate derivatives of multicenter 
* integrals that we need.  Notation:
*  i,j,p  compound basis function (in practice, shell) indices, 
*       including, specification of which center the function is on
*
*  Xi, Xj   orbital basis functions (chi), position of center not shown        
*           explicitly here
*
*  Gp      density basis function, ditto
*
*  R, R2   crystal lattice translation vectors
*
*  rnuc    position of some nucleus
*
*   1,2     electronic coordinates to be integrated (vol. elements d1,d2)
*
* I will just display the integrals whose derivatives we need, not the
*  differentiation itself.
*
* 1. orbital overlap   Int d1 Xi(1) Xj(1-R)   (i think we have this)  eq. 78
*
* 2. Kinetic energy    Int d1 Xi(1) [-del^2/2] Xj(1-R)     (ditto)      eq. 30
*
* 3. 1-center attraction  -Sum(nuc)Z(nuc) Int d1 Gp(1-R) /|1-rnuc|      eq. 60
*
* 4. 2-center repulsion  Int d1 d2 Gp(1) (1/r12) Gp(2-R)  [r12 obvious, 
*     p = some other p index] eq. 52
*
* 5. 2-center attraction -Sum(nuc)Z(nuc) Int d1 Xi(1+R2)Xj(1+R2-R)/|1-rnuc|  eq.42
*
* 6. 3-center repulsion Int d1 d2 Xi(1)Xj(1-R)(1/r12)Gp(2-R2)         eq. 41
* 
*  All equations refer to our gross paper JCP 105, 10983 (1995) which I 
*  believe you have a copy of, if not Ill bring one.
* 
* Probably you will need some more details on how we want this, or we can
* talk to you to understand how we get it, then we can reshuffle things as
* necessary.  For example do we get the derivative of any integral with
* repect to any arbitrary nucleus that we specify, or with respect to all
* nuclei in a single block of data, or only for those nuclei which are
* relevant (nonzero deriv.) for a given integral?  Call or message me as
* needed and Ill come over.)
* 
* Later
* John
***********************************************************************
* 

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
c
      logical int_normalize
      external int_normalize
      logical pderiv_compute_1e1cpe
      external pderiv_compute_1e1cpe
      logical pderiv_compute_1epe
      external pderiv_compute_1epe
      logical pderiv_compute_1eov
      external pderiv_compute_1eov
      logical pderiv_compute_1eke
      external pderiv_compute_1eke
      logical pderiv_compute_2e2c
      external pderiv_compute_2e2c
      logical pderiv_compute_2e3c
      external pderiv_compute_2e3c
c
      integer rtdb
c
      logical status
      integer basis, geom
      integer nat, nat3, size, maxg1, maxs1, maxg2, maxs2, maxg, maxs
      integer nbf, nbfsq
      integer hbuf, hbufp, hbufm, hscr, h1ec, h1efd, hxyz
      integer kbuf, kbufp, kbufm, kscr, k1ec, k1efd, kxyz
c
      task_pderiv = .false.
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &    ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &    ('bas_rtdb_load failed',911)
c
      write(6,*)' geom/basis loaded'
c
      if (.not.int_normalize(basis)) stop ' norm error 1'
c
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0)
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &    ('numbf failed',911)
c
      nbfsq = nbf*nbf
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      nat3 = 3*nat
      size = nat3*nbfsq*nbf
c
      call intd_init(rtdb,1,basis)
      call int_mem_print()
      call int_mem_1e(maxg1,maxs1)
      call int_mem_2e3c(maxg2,maxs2)
      maxg2 = maxg2*3*4
      maxg = max(maxg1,maxg2)
      maxs = max(maxs1,maxs2)
      maxg = maxg + maxg/10
      maxs = maxs + maxs/10
      write(6,*)' normal maxg1 ',maxg1
      write(6,*)' normal maxs1 ',maxs1
      write(6,*)' normal maxg2 ',maxg2
      write(6,*)' normal maxs2 ',maxs2
      write(6,*)' normal maxg  ',maxg
      write(6,*)' normal maxs  ',maxs
      status = ma_alloc_get(mt_dbl,maxg,"int buffer" ,hbuf,kbuf)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg,"int buffer plus",hbufp,kbufp)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg,"int buffer minus",hbufm,kbufm)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxs,"int scratch",hscr,kscr)
      status = status .and.
     &    ma_alloc_get(mt_dbl,size,"block c",h1ec,k1ec)
      status = status .and.
     &    ma_alloc_get(mt_dbl,size,"block fd",h1efd,k1efd)
      status = status .and.
     &    ma_alloc_get(mt_dbl,3*nat,"my coords",hxyz,kxyz)
      if (.not.status) stop ' memory alloc failed rak20 (1)'
c
      task_pderiv =  pderiv_compute_1e1cpe(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
c
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_1eov(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
c
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_1eke(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
c
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_1epe(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
c
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_2e2c(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
c
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_2e3c(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
      
      call intd_terminate()
      if (.not.bas_destroy(basis)) stop ' bas destroy fail'
      if (.not.geom_destroy(geom)) stop ' geom destroy fail'
      status = .true.
      status = status.and.ma_free_heap(hbuf)
      status = status.and.ma_free_heap(hbufp)
      status = status.and.ma_free_heap(hbufm)
      status = status.and.ma_free_heap(hscr)
      status = status.and.ma_free_heap(h1ec)
      status = status.and.ma_free_heap(h1efd)
      status = status.and.ma_free_heap(hxyz)
      task_pderiv = status.and.task_pderiv
c
      end
      logical function pderiv_compute_1e1cpe(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,p1c,p1fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision p1c(nbf,3,nat)
      double precision p1fd(nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3
      integer atom, ixyz
      double precision delta, factor, thresh, norm
      double precision R(3)
      integer nzero1, nzero2
      integer ishell, nshell, ilo, ihi, nbfsh, cnt, i
      integer IR
c
      call dfill(3,0.0d00,R,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*3*nat),0.0d00,p1c,1)
      call dfill((nbf*3*nat),0.0d00,p1fd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,2
        if (IR.eq.1) call dfill(3,0.0d00,R,1)
        if (IR.eq.2) then
          R(1) = 1.0d00
          R(2) = 2.0d00
          R(3) = 3.0d00
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 1'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfsh = ihi - ilo + 1
          do atom = 1,nat
            do ixyz = 1,3
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call dfill(lbuf,0.0d00,bufp,1)
              call intp_1e1cpe(basis,ishell,R,lscr,scr,lbuf,bufp)
*
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,bufm,1)
              call intp_1e1cpe(basis,ishell,R,lscr,scr,lbuf,bufm)
*
              call dcopy(nbfsh,bufp,1,buf,1)
              call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
              call dscal(nbfsh,factor,buf,1)
              cnt = 1
              do i = ilo,ihi
                p1fd(i,ixyz,atom) = buf(cnt)
                cnt = cnt + 1
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list '
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              if (abs(p1fd(i,ixyz,atom)).le.thresh) nzero1 = nzero1 + 1
*              write(6,10000)i,ixyz,atom,p1fd(i,ixyz,atom)
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              if (abs(p1fd(i,ixyz,atom)).gt.thresh) then
*               write(6,10000)i,ixyz,atom,p1fd(i,ixyz,atom)
                continue
              else
                nzero2 = nzero2 + 1
              endif
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*3*nat)
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          call dfill(lscr,0.0d00,scr,1)
          call dfill(lbuf,0.0d00,buf,1)
          call intpd_1e1cpe(basis,ishell,R,lscr,scr,lbuf,buf)
          cnt = 1
          do atom=1,nat
            do ixyz = 1,3
              do i = ilo, ihi
                p1c(i,ixyz,atom) = buf(cnt)
                cnt = cnt + 1
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: full list '
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              if (abs(p1c(i,ixyz,atom)).le.thresh) nzero1 = nzero1 + 1
*              write(6,10001)i,ixyz,atom,p1c(i,ixyz,atom)
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              if (abs(p1c(i,ixyz,atom)).gt.thresh) then
*                write(6,10001)i,ixyz,atom,p1c(i,ixyz,atom)
                continue
              else
                nzero2 = nzero2 + 1
              endif
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*3*nat)
c
10000   format(1x,'p1fd(',i3,',',i2,',',i3,') = ',pd20.10)
10001   format(1x,' p1c(',i3,',',i2,',',i3,') = ',pd20.10)
c
        call daxpy((nbf*3*nat),-1.0d00,p1fd,1,p1c,1)
        norm = ddot((nbf*3*nat),p1c,1,p1c,1)
        write(luout,*)' 1e1cpe difference norm:',ir,' ',norm
c
        pderiv_compute_1e1cpe = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_1eov(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,ovc,ovfd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision ovc(nbf,nbf,3,nat)
      double precision ovfd(nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3
      integer atom, ixyz, IR
      double precision delta, factor, thresh, norm
      double precision R(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(2)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
c
      call dfill(3,0.0d00,R,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*3*nat),0.0d00,ovc,1)
      call dfill((nbf*nbf*3*nat),0.0d00,ovfd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,2
        if (IR.eq.1) call dfill(3,0.0d00,R,1)
        if (IR.eq.2) then
          R(1) = 1.0d00
          R(2) = 2.0d00
          R(3) = 3.0d00
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 2'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            nbfsh = nbfshi*nbfshj
            do atom = 1,nat
              do ixyz = 1,3
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,buf,1)
                call dfill(lbuf,0.0d00,bufp,1)
                call intp_1eov(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufp)
*
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,bufm,1)
                call intp_1eov(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufm)
*
                call dcopy(nbfsh,bufp,1,buf,1)
                call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                call dscal(nbfsh,factor,buf,1)
                cnt = 1
                do i = ilo,ihi
                  do j = jlo, jhi
                    ovfd(i,j,ixyz,atom) = buf(cnt)
                    ovfd(j,i,ixyz,atom) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
                if (.not.ma_verify_allocator_stuff())
     &              stop ' ma broke 3'
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(ovfd(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10000)i,j,ixyz,atom,ovfd(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                if (abs(ovfd(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10000)i,j,ixyz,atom,ovfd(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*nbf*3*nat)
        nzero1 = 0
        nzero2 = 0
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          do jshell = 1,ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            call dfill(lscr,0.0d00,scr,1)
            call dfill(lbuf,0.0d00,buf,1)
            call intpd_1eov(basis,ishell,basis,jshell,R,lscr,scr,
     &          lbuf,buf,watom)
*            write(6,*)'watom 1eov',watom
            cnt = 1
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 4'
            do atom=1,2
              if (watom(atom).gt.0) then
                do ixyz = 1,3
                  do i = ilo, ihi
                    do j = jlo, jhi
                      ovc(i,j,ixyz,watom(atom)) = buf(cnt)
                      ovc(j,i,ixyz,watom(atom)) = buf(cnt)
                      cnt = cnt + 1
                    enddo
                  enddo
                enddo
              endif
            enddo
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 5'
          enddo
        enddo
        write(6,*)'  c: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(ovc(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10001)i,j,ixyz,atom,ovc(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(ovc(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10001)i,j,ixyz,atom,ovc(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*nbf*3*nat)
c
10000   format(1x,'ovfd(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
10001   format(1x,' ovc(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
c
        call daxpy((nbf*nbf*3*nat),-1.0d00,ovfd,1,ovc,1)
        norm = ddot((nbf*nbf*3*nat),ovc,1,ovc,1)
        write(luout,*)' 1eov difference norm:',ir,' ',norm
c
        pderiv_compute_1eov = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_1eke(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,kec,kefd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision kec(nbf,nbf,3,nat)
      double precision kefd(nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3
      integer atom, ixyz, IR
      double precision delta, factor, thresh, norm
      double precision R(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(2)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
c
      call dfill(3,0.0d00,R,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*3*nat),0.0d00,kec,1)
      call dfill((nbf*nbf*3*nat),0.0d00,kefd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,2
        if (IR.eq.1) call dfill(3,0.0d00,R,1)
        if (IR.eq.2) then
          R(1) = 1.0d00
          R(2) = 2.0d00
          R(3) = 3.0d00
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 6'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            nbfsh = nbfshi*nbfshj
            do atom = 1,nat
              do ixyz = 1,3
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,buf,1)
                call dfill(lbuf,0.0d00,bufp,1)
                call intp_1eke(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufp)
*
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,bufm,1)
                call intp_1eke(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufm)
*
                call dcopy(nbfsh,bufp,1,buf,1)
                call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                call dscal(nbfsh,factor,buf,1)
                cnt = 1
                do i = ilo,ihi
                  do j = jlo, jhi
                    kefd(i,j,ixyz,atom) = buf(cnt)
                    kefd(j,i,ixyz,atom) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
                if (.not.ma_verify_allocator_stuff())
     &              stop ' ma broke 7'
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(kefd(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10000)i,j,ixyz,atom,kefd(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                if (abs(kefd(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10000)i,j,ixyz,atom,kefd(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*nbf*3*nat)
        nzero1 = 0
        nzero2 = 0
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          do jshell = 1,ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            call dfill(lscr,0.0d00,scr,1)
            call dfill(lbuf,0.0d00,buf,1)
            call intpd_1eke(basis,ishell,basis,jshell,R,lscr,scr,
     &          lbuf,buf,watom)
*            write(6,*)'watom 1eke',watom
            cnt = 1
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 8'
            do atom=1,2
              if (watom(atom).gt.0) then
                do ixyz = 1,3
                  do i = ilo, ihi
                    do j = jlo, jhi
                      kec(i,j,ixyz,watom(atom)) = buf(cnt)
                      kec(j,i,ixyz,watom(atom)) = buf(cnt)
                      cnt = cnt + 1
                    enddo
                  enddo
                enddo
              endif
            enddo
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 9'
          enddo
        enddo
        write(6,*)'  c: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(kec(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10001)i,j,ixyz,atom,kec(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(kec(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10001)i,j,ixyz,atom,kec(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*nbf*3*nat)
c
10000   format(1x,'kefd(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
10001   format(1x,' kec(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
c
        call daxpy((nbf*nbf*3*nat),-1.0d00,kefd,1,kec,1)
        norm = ddot((nbf*nbf*3*nat),kec,1,kec,1)
        write(luout,*)' 1eke difference norm:',ir,' ',norm
c
        pderiv_compute_1eke = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_1epe(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,pec,pefd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision pec(nbf,nbf,3,nat)
      double precision pefd(nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3
      integer atom, ixyz, IR
      double precision delta, factor, thresh, norm
      double precision R1(3), R2(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
c
      call dfill(3,0.0d00,R1,1)
      call dfill(3,0.0d00,R2,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*3*nat),0.0d00,pec,1)
      call dfill((nbf*nbf*3*nat),0.0d00,pefd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,3
        if (IR.eq.1) then
          call dfill(3,0.0d00,R1,1)
          call dfill(3,0.0d00,R2,1)
        elseif (IR.eq.2) then
          call dfill(3,0.0d00,R2,1)
          R1(1) = 1.0d00
          R1(2) = 2.0d00
          R1(3) = 3.0d00
        elseif (IR.eq.3) then
          R1(1) = 1.0d00
          R1(2) = 2.0d00
          R1(3) = 3.0d00
          R2(1) = 3.0d00
          R2(2) = 4.0d00
          R2(3) = 5.0d00
        else
          stop ' how did IR get here'
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 10'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            nbfsh = nbfshi*nbfshj
            do atom = 1,nat
              do ixyz = 1,3
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,buf,1)
                call dfill(lbuf,0.0d00,bufp,1)
                call intp_1epe(basis,ishell,R1,basis,jshell,
     &              R2,lscr,scr,lbuf,bufp)
*
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,bufm,1)
                call intp_1epe(basis,ishell,R1,basis,jshell,
     &              R2,lscr,scr,lbuf,bufm)
*
                call dcopy(nbfsh,bufp,1,buf,1)
                call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                call dscal(nbfsh,factor,buf,1)
                cnt = 1
                do i = ilo,ihi
                  do j = jlo, jhi
                    pefd(i,j,ixyz,atom) = buf(cnt)
                    pefd(j,i,ixyz,atom) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
                if (.not.ma_verify_allocator_stuff())
     &              stop ' ma broke 11'
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(pefd(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10000)i,j,ixyz,atom,pefd(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                if (abs(pefd(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10000)i,j,ixyz,atom,pefd(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*nbf*3*nat)
        nzero1 = 0
        nzero2 = 0
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          do jshell = 1,ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            call dfill(lscr,0.0d00,scr,1)
            call dfill(lbuf,0.0d00,buf,1)
            call intpd_1epe(basis,ishell,R1,basis,jshell,R2,
     &          lscr,scr,lbuf,buf)
            cnt = 1
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 12'
            do atom=1,nat
              do ixyz = 1,3
                do i = ilo, ihi
                  do j = jlo, jhi
                    pec(i,j,ixyz,atom) = buf(cnt)
                    pec(j,i,ixyz,atom) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
              enddo
            enddo
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 13'
          enddo
        enddo
        write(6,*)'  c: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(pec(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10001)i,j,ixyz,atom,pec(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(pec(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10001)i,j,ixyz,atom,pec(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*nbf*3*nat)
c
10000   format(1x,'pefd(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
10001   format(1x,' pec(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
c
        call daxpy((nbf*nbf*3*nat),-1.0d00,pefd,1,pec,1)
        norm = ddot((nbf*nbf*3*nat),pec,1,pec,1)
        write(luout,*)' 1epe difference norm:',IR,' ',norm
c
        pderiv_compute_1epe = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_2e2c(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,eri2c,eri2fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision eri2c(nbf,nbf,3,nat)
      double precision eri2fd(nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3
      integer atom, ixyz, IR
      double precision delta, factor, thresh, norm
      double precision R(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(2)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
c
      call dfill(3,0.0d00,R,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*3*nat),0.0d00,eri2c,1)
      call dfill((nbf*nbf*3*nat),0.0d00,eri2fd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,2
        if (IR.eq.1) call dfill(3,0.0d00,R,1)
        if (IR.eq.2) then
          R(1) = 1.0d00
          R(2) = 2.0d00
          R(3) = 3.0d00
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 14'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            nbfsh = nbfshi*nbfshj
            do atom = 1,nat
              do ixyz = 1,3
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,buf,1)
                call dfill(lbuf,0.0d00,bufp,1)
                call intp_2e2c(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufp)
*
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,bufm,1)
                call intp_2e2c(basis,ishell,basis,jshell,
     &              R,lscr,scr,lbuf,bufm)
*
                call dcopy(nbfsh,bufp,1,buf,1)
                call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                call dscal(nbfsh,factor,buf,1)
                cnt = 1
                do i = ilo,ihi
                  do j = jlo, jhi
                    eri2fd(i,j,ixyz,atom) = buf(cnt)
                    eri2fd(j,i,ixyz,atom) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
                if (.not.ma_verify_allocator_stuff())
     &              stop ' ma broke 15'
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(eri2fd(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10000)i,j,ixyz,atom,eri2fd(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                if (abs(eri2fd(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10000)i,j,ixyz,atom,eri2fd(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*nbf*3*nat)
        nzero1 = 0
        nzero2 = 0
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          do jshell = 1,ishell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            call dfill(lscr,0.0d00,scr,1)
            call dfill(lbuf,0.0d00,buf,1)
*            write(6,*)'rak20: ishell, jshell ',ishell,jshell
            call intpd_2e2c(basis,ishell,basis,jshell,R,lscr,scr,
     &          lbuf,buf,watom)
*            write(6,*)'watom 2e2c',watom
            cnt = 1
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 16'
            do atom=1,2
              if (watom(atom).gt.0) then
                do ixyz = 1,3
                  do i = ilo, ihi
                    do j = jlo, jhi
                      eri2c(i,j,ixyz,watom(atom)) = buf(cnt)
                      eri2c(j,i,ixyz,watom(atom)) = buf(cnt)
                      cnt = cnt + 1
                    enddo
                  enddo
                enddo
              endif
            enddo
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 17'
          enddo
        enddo
        write(6,*)'  c: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(eri2c(i,j,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                write(6,10001)i,j,ixyz,atom,eri2c(i,j,ixyz,atom)
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                if (abs(eri2c(i,j,ixyz,atom)).gt.thresh) then
*                  write(6,10001)i,j,ixyz,atom,eri2c(i,j,ixyz,atom)
                  continue
                else
                  nzero2 = nzero2 + 1
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*nbf*3*nat)
c
10000   format(1x,'eri2fd(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
10001   format(1x,' eri2c(',i3,',',i3,',',i2,',',i3,') = ',pd20.10)
c
        call daxpy((nbf*nbf*3*nat),-1.0d00,eri2fd,1,eri2c,1)
        norm = ddot((nbf*nbf*3*nat),eri2c,1,eri2c,1)
        write(luout,*)' 2e2c difference norm:',ir,' ',norm
c
        pderiv_compute_2e2c = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_2e3c(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,eri3c,eri3fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "geomP.fh"
#include "bas.fh"
#include "stdio.fh"
c
      double precision ddot
      external ddot
c
      integer geom
      integer basis
      integer nbf
      integer nat
      integer lbuf
      integer lscr
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision eri3c(nbf,nbf,nbf,3,nat)
      double precision eri3fd(nbf,nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3, nintz
      integer atom, ixyz, IR
      double precision delta, factor, thresh, norm
      double precision R1(3), R2(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(4)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
      integer kshell, klo, khi, nbfshk, k
c
      call dfill(3,0.0d00,R1,1)
      call dfill(3,0.0d00,R2,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*nbf*3*nat),0.0d00,eri3c,1)
      call dfill((nbf*nbf*nbf*3*nat),0.0d00,eri3fd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,3
        if (IR.eq.1) then
          call dfill(3,0.0d00,R1,1)
          call dfill(3,0.0d00,R2,1)
        elseif (IR.eq.2) then
          call dfill(3,0.0d00,R2,1)
          R1(1) = 1.0d00
          R1(2) = 2.0d00
          R1(3) = 3.0d00
        elseif (IR.eq.3) then
          R1(1) = 1.0d00
          R1(2) = 2.0d00
          R1(3) = 3.0d00
          R2(1) = 3.0d00
          R2(2) = 4.0d00
          R2(3) = 5.0d00
        else
          stop ' how did IR get here'
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 18'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, nshell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            do kshell = 1, jshell
              if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &            stop 'cn2bfr error k'
              nbfshk = khi - klo + 1
              nbfsh = nbfshi*nbfshj*nbfshk
              do atom = 1,nat
                do ixyz = 1,3
                  call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                  coords(ixyz,atom,geom) =
     &                coords(ixyz,atom,geom) + delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,buf,1)
                  call dfill(lbuf,0.0d00,bufp,1)
                  call intp_2e3c(basis,ishell,basis,jshell,kshell,
     &              R1,R2,lscr,scr,lbuf,bufp)
*
                  call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                  coords(ixyz,atom,geom) =
     &                coords(ixyz,atom,geom) - delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,bufm,1)
                  call intp_2e3c(basis,ishell,basis,jshell,kshell,
     &              R1,R2,lscr,scr,lbuf,bufm)
*
                  call dcopy(nbfsh,bufp,1,buf,1)
                  call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                  call dscal(nbfsh,factor,buf,1)
                  cnt = 1
                  do i = ilo,ihi
                    do j = jlo, jhi
                      do k = klo, khi
                        eri3fd(i,j,k,ixyz,atom) = buf(cnt)
                        eri3fd(i,k,j,ixyz,atom) = buf(cnt)
                        cnt = cnt + 1
                      enddo
                    enddo
                  enddo
                  if (.not.ma_verify_allocator_stuff())
     &                stop ' ma broke 19'
                enddo
              enddo
            enddo
          enddo
        enddo
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,nbf
                  if (abs(eri3fd(i,j,k,ixyz,atom)).le.thresh)
     &                nzero1 = nzero1 + 1
*                  write(6,10000)i,j,k,ixyz,atom,eri3fd(i,j,k,ixyz,atom)
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                do k = 1,nbf
                  if (abs(eri3fd(i,j,k,ixyz,atom)).gt.thresh) then
*                    write(6,10000)i,j,k,ixyz,atom,
*     &                  eri3fd(i,j,k,ixyz,atom)
                    continue
                  else
                    nzero2 = nzero2 + 1
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: num zeros :1:', nzero1 
        write(6,*)' fd: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)' fd: non-zero  :1:', ((nbf*nbf*nbf*3*nat)-nzero1)
        write(6,*)' fd: non-zero  :2:', ((nbf*nbf*nbf*3*nat)-nzero2)
        write(6,*)' fd: possible  : :', (nbf*nbf*nbf*3*nat)
        nzero1 = 0
        nzero2 = 0
c
        call dcopy(nat3,xyz,1,coords(1,1,geom),1)
        do ishell = 1,nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          do jshell = 1,nshell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            do kshell = 1, jshell
              if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &            stop 'cn2bfr error k'
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
*              write(6,*)'rak20: ishell, jshell, kshell ',
*     &            ishell,jshell,kshell
              call intpd_2e3c(basis,ishell,basis,jshell,kshell,R1,R2,
     &            lscr,scr,lbuf,buf,watom)
*              write(6,*)'watom 2e3c',watom
              if (.not.ma_verify_allocator_stuff()) stop ' ma broke 20'
              nintz = (ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)
              do atom=1,4
                if (watom(atom).gt.0) then
                  cnt = ((atom-1)*nintz*3) + 1
                  do ixyz = 1,3
                    do i = ilo, ihi
                      do j = jlo, jhi
                        do k = klo, khi
                          eri3c(i,j,k,ixyz,watom(atom)) = buf(cnt)
                          eri3c(i,k,j,ixyz,watom(atom)) = buf(cnt)
                          cnt = cnt + 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
              enddo
              if (.not.ma_verify_allocator_stuff()) stop ' ma broke 21'
            enddo
          enddo
        enddo
        write(6,*)'  c: full list i,j'
        nzero1 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,nbf
                  if (abs(eri3c(i,j,k,ixyz,atom)).le.thresh)
     &                nzero1 = nzero1 + 1
*                  write(6,10001)i,j,k,ixyz,atom,eri3c(i,j,k,ixyz,atom)
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,nbf
                  if (abs(eri3c(i,j,k,ixyz,atom)).gt.thresh) then
*                    write(6,10001)i,j,k,ixyz,atom,
*     &                  eri3c(i,j,k,ixyz,atom)
                    continue
                  else
                    nzero2 = nzero2 + 1
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: num zeros :1:', nzero1 
        write(6,*)'  c: num zeros :2:', nzero2, ':delta:',
     &      (nzero2-nzero1)
        write(6,*)'  c: non-zero  :1:', ((nbf*nbf*nbf*3*nat)-nzero1)
        write(6,*)'  c: non-zero  :2:', ((nbf*nbf*nbf*3*nat)-nzero2)
        write(6,*)'  c: possible  : :', (nbf*nbf*nbf*3*nat)
c
10000   format(1x,'eri3fd(',i3,',',i3,',',i3,',',i2,',',i3,
     &      ') = ',pd20.10)
10001   format(1x,' eri3c(',i3,',',i3,',',i3,',',i2,',',i3,
     &      ') = ',pd20.10)
c
        call daxpy((nbf*nbf*nbf*3*nat),-1.0d00,eri3fd,1,eri3c,1)
        norm = ddot((nbf*nbf*nbf*3*nat),eri3c,1,eri3c,1)
        write(luout,*)' 2e3c difference norm:',ir,' ',norm
c
        pderiv_compute_2e3c = norm.lt.thresh
      enddo
      end
