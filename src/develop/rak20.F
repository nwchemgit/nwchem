C> \ingroup task
C> @{
      logical function task_pderiv(rtdb)
*
* $Id$
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
* All equations refer to our gross paper JCP 105, 10983 (1995),
* doi:10.1063/1.472866 of which I believe you have a copy of, if not I will
* bring one.
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
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
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
      logical pderiv_compute_mpole
      external pderiv_compute_mpole
      logical pderiv_compute_3ov
      external pderiv_compute_3ov
      logical pderiv_compute_2e3ct, pderiv_compute_d2e3ct
      external pderiv_compute_2e3ct, pderiv_compute_d2e3ct
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
     &    ('geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911, RTDB_ERR)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &    ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &    ('bas_rtdb_load failed',911, RTDB_ERR)
c
      write(6,*)' geom/basis loaded'
c
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
c
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0, BASIS_ERR)
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &    ('numbf failed',911, BASIS_ERR)
c
      nbfsq = nbf*nbf
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      write(6,*) 'number of atoms ', nat
      nat3 = 3*nat
      size = max(nbf,56)  ! slmax = 35 for mpole
      size = nat3*nbfsq*size
      size = 2*size
c
      call intd_init(rtdb,1,basis)
      call int_mem_print()
      call int_mem_1e(maxg1,maxs1)
      call int_mem_2e3c(maxg2,maxs2)
      maxg2 = maxg2*3*4
      maxg1 = maxg1*35  ! mpole
      maxg = max(maxg1,maxg2)
      maxs = max(maxs1,maxs2)
      maxg = maxg + maxg/10
      maxs = maxs + maxs/10
      maxs = 2*maxs
      maxg = 2*maxg 
      write(6,*)' normal maxg1 ',maxg1
      write(6,*)' normal maxs1 ',maxs1
      write(6,*)' normal maxg2 ',maxg2
      write(6,*)' normal maxs2 ',maxs2
      write(6,*)' normal maxg  ',maxg
      write(6,*)' normal maxs  ',maxs
      status = ma_alloc_get(mt_dbl,maxg,'int buffer' ,hbuf,kbuf)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg,'int buffer plus',hbufp,kbufp)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg,'int buffer minus',hbufm,kbufm)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxs,'int scratch',hscr,kscr)
      status = status .and.
     &    ma_alloc_get(mt_dbl,size,'block c',h1ec,k1ec)
      status = status .and.
     &    ma_alloc_get(mt_dbl,size,'block fd',h1efd,k1efd)
      status = status .and.
     &    ma_alloc_get(mt_dbl,3*nat,'my coords',hxyz,kxyz)
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
      
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_3ov(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))
      
      task_pderiv = task_pderiv .and.
     &    pderiv_compute_mpole(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))

      task_pderiv = task_pderiv .and.
     &    pderiv_compute_2e3ct(
     &    geom,basis,nbf,nat,
     &    maxg,maxs,
     &    dbl_mb(kbuf),dbl_mb(kbufp),dbl_mb(kbufm),
     &    dbl_mb(kscr),dbl_mb(k1ec),dbl_mb(k1efd),
     &    dbl_mb(kxyz))

      task_pderiv = task_pderiv .and.
     &    pderiv_compute_d2e3ct(
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
C> @}
      logical function pderiv_compute_1e1cpe(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,p1c,p1fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
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
10000   format(1x,'p1fd(',i3,',',i2,',',i3,') = ',1pd20.10)
10001   format(1x,' p1c(',i3,',',i2,',',i3,') = ',1pd20.10)
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
#include "nwc_const.fh"
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
*buffer is <jlo:jhi, ilo:ihi, 1:3, 2>
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
10000   format(1x,'ovfd(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
10001   format(1x,' ovc(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
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
#include "nwc_const.fh"
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
10000   format(1x,'kefd(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
10001   format(1x,' kec(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
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
#include "nwc_const.fh"
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
10000   format(1x,'pefd(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
10001   format(1x,' pec(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
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
#include "nwc_const.fh"
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
10000   format(1x,'eri2fd(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
10001   format(1x,' eri2c(',i3,',',i3,',',i2,',',i3,') = ',1pd20.10)
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
#include "nwc_const.fh"
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
     &      ') = ',1pd20.10)
10001   format(1x,' eri3c(',i3,',',i3,',',i3,',',i2,',',i3,
     &      ') = ',1pd20.10)
c
        call daxpy((nbf*nbf*nbf*3*nat),-1.0d00,eri3fd,1,eri3c,1)
        norm = ddot((nbf*nbf*nbf*3*nat),eri3c,1,eri3c,1)
        write(luout,*)' 2e3c difference norm:',ir,' ',norm
c
        pderiv_compute_2e3c = norm.lt.thresh
      enddo
      end
      logical function pderiv_compute_3ov(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,ov3c,ov3fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
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
      double precision ov3c(nbf,nbf,nbf,3,nat)
      double precision ov3fd(nbf,nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3, nintz
      integer atom, ixyz
      double precision delta, factor, thresh, norm
*rak:      double precision zeta(50)
*rak:      integer ztype, znp, zng, zsph
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(3)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
      integer kshell, klo, khi, nbfshk, k
c
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*nbf*3*nat),0.0d00,ov3c,1)
      call dfill((nbf*nbf*nbf*3*nat),0.0d00,ov3fd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      if (.not.ma_verify_allocator_stuff()) stop ' ma broke 18'
      thresh = 1.0d-10
      delta = 0.0001d00
      factor = 1.0d00/(2.0d00*delta)
      nat3 = 3*nat
* store original coordintates
      call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
      if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
*
*rak:      if (.not.bas_continfo(basis,1,ztype,znp,zng,zsph)) stop 'e1'
*rak:      if (znp.gt.50) stop ' zeta dimension too small '
*rak:      call dfill(znp,0.0d00,zeta,1)
*rak:      if (.not.bas_get_exponent(basis,1,zeta)) stop 'bas_get_e failed'
*rak:      zeta(1) = 1.0d-50
*rak:      if (.not.bas_set_exponent(basis,1,zeta,znp))
*rak:     &    stop 'bas_set_e failed'
*
      do ishell = 1, nshell
        if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &      stop 'cn2bfr error i'
        nbfshi = ihi - ilo + 1
        do jshell = 1, nshell
          if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &        stop 'cn2bfr error j'
          nbfshj = jhi - jlo + 1
          do kshell = 1, nshell
            if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &          stop 'cn2bfr error k'
            nbfshk = khi - klo + 1
            nbfsh = nbfshi*nbfshj*nbfshk
            do atom = 1,nat
              do ixyz = 1,3
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) =
     &              coords(ixyz,atom,geom) + delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,buf,1)
                call dfill(lbuf,0.0d00,bufp,1)

                call int_1e3ov(
     &              basis,ishell,basis,jshell,basis,kshell,
     &              lscr,scr,lbuf,bufp)
*
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                coords(ixyz,atom,geom) =
     &              coords(ixyz,atom,geom) - delta
                call dfill(lscr,0.0d00,scr,1)
                call dfill(lbuf,0.0d00,bufm,1)
                call int_1e3ov(
     &              basis,ishell,basis,jshell,basis,kshell,
     &              lscr,scr,lbuf,bufm)
*
                call dcopy(nbfsh,bufp,1,buf,1)
                call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
                call dscal(nbfsh,factor,buf,1)
                cnt = 1
                do i = ilo,ihi
                  do j = jlo, jhi
                    do k = klo, khi
                      ov3fd(i,j,k,ixyz,atom) = buf(cnt)
                      cnt = cnt + 1
                    enddo
                  enddo
                enddo
                if (.not.ma_verify_allocator_stuff())
     &              stop ' ma broke 19'
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
                if (abs(ov3fd(i,j,k,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                  write(6,10000)i,j,k,ixyz,atom,ov3fd(i,j,k,ixyz,atom)
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
                if (abs(ov3fd(i,j,k,ixyz,atom)).gt.thresh) then
*                    write(6,10000)i,j,k,ixyz,atom,
*     &                  ov3fd(i,j,k,ixyz,atom)
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
     &    (nzero2-nzero1)
      write(6,*)' fd: non-zero  :1:', ((nbf*nbf*nbf*3*nat)-nzero1)
      write(6,*)' fd: non-zero  :2:', ((nbf*nbf*nbf*3*nat)-nzero2)
      write(6,*)' fd: possible  : :', (nbf*nbf*nbf*3*nat)
      nzero1 = 0
      nzero2 = 0
c
      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
      do ishell = 1,nshell
        if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &      stop 'cn2bfr error i'
        do jshell = 1,nshell
          if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &        stop 'cn2bfr error j'
          do kshell = 1, jshell
            if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &          stop 'cn2bfr error k'
            call dfill(lscr,0.0d00,scr,1)
            call dfill(lbuf,0.0d00,buf,1)
*              write(6,*)'rak20: ishell, jshell, kshell ',
*     &            ishell,jshell,kshell
            call intd_1e3ov(
     &          basis,ishell,basis,jshell,basis,kshell,
     &          lscr,scr,lbuf,buf,watom)
*              write(6,*)'watom 3ov',watom
            if (.not.ma_verify_allocator_stuff()) stop ' ma broke 20'
            nintz = (ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)
            do atom=1,3
              if (watom(atom).gt.0) then
                cnt = ((atom-1)*nintz*3) + 1
                do ixyz = 1,3
                  do i = ilo, ihi
                    do j = jlo, jhi
                      do k = klo, khi
                        ov3c(i,j,k,ixyz,watom(atom)) = buf(cnt)
                        ov3c(i,k,j,ixyz,watom(atom)) = buf(cnt)
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
                if (abs(ov3c(i,j,k,ixyz,atom)).le.thresh)
     &              nzero1 = nzero1 + 1
*                  write(6,10001)i,j,k,ixyz,atom,ov3c(i,j,k,ixyz,atom)
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
                if (abs(ov3c(i,j,k,ixyz,atom)).gt.thresh) then
*                    write(6,10001)i,j,k,ixyz,atom,
*     &                  ov3c(i,j,k,ixyz,atom)
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
     &    (nzero2-nzero1)
      write(6,*)'  c: non-zero  :1:', ((nbf*nbf*nbf*3*nat)-nzero1)
      write(6,*)'  c: non-zero  :2:', ((nbf*nbf*nbf*3*nat)-nzero2)
      write(6,*)'  c: possible  : :', (nbf*nbf*nbf*3*nat)
c
10000 format(1x,'ov3fd(',i3,',',i3,',',i3,',',i2,',',i3,
     &    ') = ',1pd20.10)
10001 format(1x,' ov3c(',i3,',',i3,',',i3,',',i2,',',i3,
     &    ') = ',1pd20.10)
c
      call daxpy((nbf*nbf*nbf*3*nat),-1.0d00,ov3fd,1,ov3c,1)
      norm = ddot((nbf*nbf*nbf*3*nat),ov3c,1,ov3c,1)
      write(luout,*)' 3ov difference norm: ',norm
c
      pderiv_compute_3ov = norm.lt.thresh
      end
      logical function pderiv_compute_mpole(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,mp3c,mp3fd,xyz)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
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
      integer lval_max
      parameter (lval_max = 5)
      integer slmax 
*      parameter (slmax = 1)  ! for lval_max = 0
*      parameter (slmax = 4)  ! for lval_max = 1
*      parameter (slmax = 35) ! ! for lval_max = 4
      parameter (slmax = 56) ! 1 +3+6+10+15+21 = 35+21 = 56
      double precision buf(lbuf), bufp(lbuf), bufm(lbuf)
      double precision scr(lscr)
      double precision mp3c(slmax,nbf,nbf,3,(nat+1))   ! add one for the multipole center 
      double precision mp3fd(slmax,nbf,nbf,3,(nat+1))  ! ditto
      double precision xyz(3,nat)
c
      integer nat3, nintz
      integer atom, ixyz, IR, lval, nbflval
      double precision delta, factor, thresh, norm
      double precision R1(3), l_center(3), l_orgcenter(3)
      integer nzero1, nzero2
      integer nshell, cnt, nbfsh, watom(3), nshell_bas
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
      integer k, klo, khi, num_ret
      integer h_diff, k_diff
      integer imp_cent
c
      integer pole_index
      pole_index(lval) = lval*((lval+3)*(lval+3)+2)/6 + 1
c
      call dfill(3,0.0d00,R1,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((slmax*nbf*nbf*3*(nat+1)),0.0d00,mp3c,1)
      call dfill((slmax*nbf*nbf*3*(nat+1)),0.0d00,mp3fd,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do imp_cent = 1, 2
        if (imp_cent.eq.1) then
          call dfill(3,0.0d00,l_center,1) ! start with multipole center at origin
        else if (imp_cent.eq.2) then
          call dcopy(3,coords(1,1,geom),1,l_center,1)
        endif
        call dcopy(3,l_center,1,l_orgcenter,1)

      do IR = 1, 3
        if (IR.eq.1) then
          call dfill(3,0.0d00,R1,1)
        elseif (IR.eq.2) then
          R1(1) = 1.0d00
          R1(2) = 2.0d00
          R1(3) = 3.0d00
        elseif (IR.eq.3) then
          R1(1) = 1.0d00
          R1(2) = 1.0d00
          R1(3) = 1.0d00
        else
          stop ' how did IR get here'
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 18'
        thresh = 1.0d-10
        delta = 0.0001d00
        factor = 1.0d00/(2.0d00*delta)
*rak:        write(6,*)' factor = ',factor
        nat3 = 3*nat
* store original coordintates
        call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*rak:        write(6,*)' original coords '
*rak:        call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:        write(6,*)' l_center',l_center
*rak:        write(6,*)' l_orgcenter',l_orgcenter
*
        if (.not.bas_numcont(basis,nshell_bas)) stop 'bas_numcont'
        nshell = nshell_bas
        nshell = 1 ! tmp change
*rak:        write(6,*)' nshell = ',nshell
        
        do ishell = 1, nshell
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, nshell
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            do lval = 0,lval_max
*rak:              write(6,*)'................................1',
*rak:     &            ishell,jshell,lval
              nbflval = (lval+1)*(lval+2)/2
              nbfsh = nbfshj*nbfshi*nbflval
              klo = pole_index((lval-1)) + 1
              klo = max(klo,1)  ! minumum value of klo is 1
              khi = pole_index(lval)
*rak:              write(6,*)' klo, khi',klo,khi
              do atom = 1,nat
                do ixyz = 1,3
                  call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                  coords(ixyz,atom,geom) =
     &                coords(ixyz,atom,geom) + delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,buf,1)
                  call dfill(lbuf,0.0d00,bufp,1)
*rak:                  write(6,*)' coords intp_mpolel(1)',atom,ixyz
*rak:                  call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:                  write(6,*)' l_center',l_center
*rak:                  write(6,*)' l_orgcenter',l_orgcenter,' -*'
                  call intp_mpolel(basis,ishell,basis,jshell,
     &              R1,lval,l_orgcenter,lscr,scr,lbuf,bufp,num_ret)
*                  call intp_1eov(basis,ishell,basis,jshell,
*     &                R1,lscr,scr,lbuf,buf)
*                  cnt = 0
*                  write(6,*)
*     &                ' i   j  cnt  mpole  ovl   diff + atom ixyz'
*                  do i = ilo,ihi
*                    do j = jlo,jhi
*                      cnt = cnt + 1
*                      norm = abs(bufp(cnt)-buf(cnt))
*                      write(6,*)i,j,cnt,bufp(cnt),buf(cnt),norm,
*     &                    atom, ixyz
*                    enddo
*                  enddo
*                  norm = 0.0d00
*                  call dfill(lbuf,0.0d00,buf,1)
*
                  call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                  coords(ixyz,atom,geom) =
     &                coords(ixyz,atom,geom) - delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,bufm,1)
*rak:                  write(6,*)' coords intp_mpolel(2)',atom,ixyz
*rak:                  call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:                  write(6,*)' l_center',l_center
*rak:                  write(6,*)' l_orgcenter',l_orgcenter,' -*'
                  call intp_mpolel(basis,ishell,basis,jshell,
     &              R1,lval,l_orgcenter,lscr,scr,lbuf,bufm,num_ret)
*rak:                  call intp_1eov(basis,ishell,basis,jshell,
*rak:     &                R1,lscr,scr,lbuf,buf)
*rak:                  cnt = 0
*rak:                  write(6,*)
*rak:     &                ' i   j  cnt  mpole  ovl   diff - atom ixyz'
*rak:                  do i = ilo,ihi
*rak:                    do j = jlo,jhi
*rak:                      cnt = cnt + 1
*rak:                      norm = abs(bufp(cnt)-buf(cnt))
*rak:                      write(6,*)i,j,cnt,bufp(cnt),buf(cnt),norm,
*rak:     &                    atom,ixyz
*rak:                    enddo
*rak:                  enddo
*rak:                  norm = 0.0d00
*rak:                  call dfill(lbuf,0.0d00,buf,1)
*
*rak:                  call prak3mat(bufp,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'bufp in atom move',atom,ixyz)
*rak:                  call prak3mat(bufm,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'bufm in atom move',atom,ixyz)
                  call dcopy(nbfsh,bufp,1,buf,1)
                  call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
*rak:                  call prak3mat(buf,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'buf before scale in atom move',atom,ixyz)
                  call dscal(nbfsh,factor,buf,1)
*rak:                  call prak3mat(buf,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'buf in atom move',atom,ixyz)
                  cnt = 1
                  do i = ilo,ihi
                    do j = jlo, jhi
                      do k = klo, khi
                        mp3fd(k,j,i,ixyz,atom) = buf(cnt)
*rak:                        write(70,*)
*rak:     &                      'fd<',k,j,i,ixyz,atom,'>=',buf(cnt)
                        cnt = cnt + 1
                      enddo
                    enddo
                  enddo
*rak:                  write(6,*)' nbfsh ',nbfsh,' cnt ',cnt
                  if (.not.ma_verify_allocator_stuff())
     &                stop ' ma broke 19'
                enddo
              enddo
*             now do multipole center move
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)  ! restore original atom coords
              do ixyz = 1,3
                  call dcopy(3,l_orgcenter,1,l_center,1)
                  l_center(ixyz) = l_center(ixyz) + delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,buf,1)
                  call dfill(lbuf,0.0d00,bufp,1)
*rak:                  write(6,*)' coords intp_mpolel(3)',atom,ixyz
*rak:                  call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:                  write(6,*)' l_center',l_center,' -*'
*rak:                  write(6,*)' l_orgcenter',l_orgcenter
                  call intp_mpolel(basis,ishell,basis,jshell,
     &              R1,lval,l_center,lscr,scr,lbuf,bufp,num_ret)
*
                  call dcopy(3,l_orgcenter,1,l_center,1)
                  l_center(ixyz) = l_center(ixyz) - delta
                  call dfill(lscr,0.0d00,scr,1)
                  call dfill(lbuf,0.0d00,bufm,1)
*rak:                  write(6,*)' coords intp_mpolel(4)',atom,ixyz
*rak:                  call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:                  write(6,*)' l_center',l_center,' -*'
*rak:                  write(6,*)' l_orgcenter',l_orgcenter
                  call intp_mpolel(basis,ishell,basis,jshell,
     &              R1,lval,l_center,lscr,scr,lbuf,bufm,num_ret)
*
*rak:                  call prak3mat(bufp,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,                  
*rak:     &                'bufp',(nat+1),ixyz)
*rak:                  call prak3mat(bufm,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,                  
*rak:     &                'bufm',(nat+1),ixyz)
                  call dcopy(nbfsh,bufp,1,buf,1)
                  call daxpy(nbfsh,-1.0d00,bufm,1,buf,1)
*rak:                  call prak3mat(buf,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'buf b4 scale in mp center move',(nat+1),ixyz)
                  call dscal(nbfsh,factor,buf,1)
*rak:                  call prak3mat(buf,nbfsh,nbfshi,nbfshj,nbflval,
*rak:     &                ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                'buf in mp center move',(nat+1),ixyz)
                  cnt = 1
                  do i = ilo,ihi
                    do j = jlo, jhi
                      do k = klo, khi
                        mp3fd(k,j,i,ixyz,nat+1) = buf(cnt)
*rak:                        write(70,*)
*rak:     &                      'fd<',k,j,i,ixyz,(nat+1),'>=',buf(cnt)
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
c
        write(6,*)' fd: full list i,j'
        nzero1 = 0
        do atom = 1,nat+1
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,slmax
                  if (abs(mp3fd(k,j,i,ixyz,atom)).le.thresh)
     &                nzero1 = nzero1 + 1
*                  write(6,10000)k,j,i,ixyz,atom,mp3fd(k,j,i,ixyz,atom)
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)' fd: non-zero list '
        nzero2 = 0
        do atom = 1,nat+1
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf 
                do k = 1,slmax
                  if (abs(mp3fd(k,j,i,ixyz,atom)).gt.thresh) then
*                    write(6,10000)k,j,i,ixyz,atom,
*     &                  mp3fd(k,j,i,ixyz,atom)
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
        write(6,*)' fd: non-zero  :1:',
     &      ((nbf*nbf*slmax*3*(nat+1))-nzero1)
        write(6,*)' fd: non-zero  :2:',
     &      ((nbf*nbf*slmax*3*(nat+1))-nzero2)
        write(6,*)' fd: possible  : :',
     &      (nbf*nbf*slmax*3*(nat+1))
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
            do lval = 0, lval_max
*rak:              write(6,*)'................................2',
*rak:     &            ishell,jshell,lval
              klo = pole_index((lval-1)) + 1
              klo = max(klo,1)  ! minumum value of klo is 1
              khi = pole_index(lval)
*rak:              write(6,*)' klo, khi',klo,khi
c
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call dfill(lbuf,0.0d00,bufp,1)
*rak:              write(6,*)'rak20: ishell, jshell,lval ',
*rak:     &            ishell,jshell,lval
*rak:              write(6,*)' lscr ', lscr
*rak:              call intpd_1eov(basis,ishell,basis,jshell,R1,lscr,scr,
*rak:     &            lbuf,bufp,watom)
*rak:              write(6,*)' coords intpd_mpolel'
*rak:              call output(coords(1,1,geom),1,3,1,nat,3,nat,1)
*rak:              write(6,*)' l_center',l_center
*rak:              write(6,*)' l_orgcenter',l_orgcenter,' -*'
              call intpd_mpolel(basis,ishell,basis,jshell,R1,
     &            lval,l_orgcenter,
     &            lscr,scr,lbuf,buf,num_ret,watom)
*rak:              cnt = 0
*rak:              write(6,*) ' atom ixyz  i   j  cnt  mpole  ovl   diff '
*rak:              do atom = 1,nat+1
*rak:                do ixyz = 1,3
*rak:                  do i = ilo,ihi
*rak:                    do j = jlo,jhi
*rak:                      cnt = cnt + 1
*rak:                      norm = abs(bufp(cnt)-buf(cnt))
*rak:                      write(6,*)
*rak:     &                    atom,ixyz,i,j,cnt,buf(cnt),bufp(cnt),norm
*rak:                    enddo
*rak:                  enddo
*rak:                enddo
*rak:              enddo
*rak:              call dfill(lbuf,0.0d00,bufp,1)
*rak:              write(6,*)'watom mpole',watom
              if (.not.ma_verify_allocator_stuff()) stop ' ma broke 20'
              nintz = (ihi-ilo+1)*(jhi-jlo+1)*((lval+1)*(lval+2)/2)
*rak:              write(6,*)'nintz ',nintz
*rak:              write(6,*)'ilo/hi',ilo,ihi
*rak:              write(6,*)'jlo/hi',jlo,jhi
*rak:              write(6,*)'klo/hi',klo,khi
              do atom=1,3
*rak:                write(6,*)' atom,watom(atom)',atom,watom(atom)
                if (watom(atom).gt.0) then
                  cnt = ((atom-1)*nintz*3) + 1
                  do ixyz = 1,3
*rak:                    k_diff = cnt + (ixyz-1)*nintz
*rak:                    call prak3mat(buf(k_diff),nbfsh,
*rak:     &                  nbfshi,nbfshj,nbflval,
*rak:     &                  ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                  'buf in derivs call(1)',watom(atom),ixyz)
                    do i = ilo, ihi
                      do j = jlo, jhi
                        do k = klo, khi
                          mp3c(k,j,i,ixyz,watom(atom)) = buf(cnt)
*rak:                          write(71,*)
*rak:     &                        ' c<',k,j,i,ixyz,watom(atom),'>=',buf(cnt)
                          cnt = cnt + 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
                if (watom(atom).eq.-3) then
                  cnt = ((atom-1)*nintz*3) + 1
                  do ixyz = 1,3
*rak:                    k_diff = cnt + (ixyz-1)*nintz
*rak:                    call prak3mat(buf(k_diff),nbfsh,
*rak:     &                  nbfshi,nbfshj,nbflval,
*rak:     &                  ilo,ihi,jlo,jhi,klo,khi,
*rak:     &                  'buf in derivs call(2)',watom(atom),ixyz)
                    do i = ilo, ihi
                      do j = jlo, jhi
                        do k = klo, khi
                          mp3c(k,j,i,ixyz,nat+1) = buf(cnt)
*rak:                          write(71,*)
*rak:     &                        ' c<',k,j,i,ixyz,(nat+1),'>=',buf(cnt)
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
        do atom = 1,nat+1
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,slmax
                  if (abs(mp3c(k,j,i,ixyz,atom)).le.thresh)
     &                nzero1 = nzero1 + 1
*                  write(6,10001)k,j,i,ixyz,atom,mp3c(k,j,i,ixyz,atom)
                enddo
              enddo
            enddo
          enddo
        enddo
        write(6,*)'  c: non-zero list '
        nzero2 = 0
        do atom = 1,nat+1
          do ixyz = 1,3
            do i = 1,nbf
              do j = 1,nbf
                do k = 1,slmax
                  if (abs(mp3c(k,j,i,ixyz,atom)).gt.thresh) then
*                    write(6,10001)k,j,i,ixyz,atom,
*     &                  mp3c(k,j,i,ixyz,atom)
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
        write(6,*)'  c: non-zero  :1:',
     &      ((nbf*nbf*slmax*3*(nat+1))-nzero1)
        write(6,*)'  c: non-zero  :2:',
     &      ((nbf*nbf*slmax*3*(nat+1))-nzero2)
        write(6,*)'  c: possible  : :',
     &      (nbf*nbf*slmax*3*(nat+1))
c
10000   format(1x,'mp3fd(',i3,',',i3,',',i3,',',i2,',',i3,
     &      ') = ',1pd20.10)
10001   format(1x,' mp3c(',i3,',',i3,',',i3,',',i2,',',i3,
     &      ') = ',1pd20.10)
c
        if (.not.ma_alloc_get(mt_dbl,
     &      (nbf*nbf*slmax*3*(nat+1)),
     &      'diff buffer',
     &      h_diff,
     &      k_diff)) stop 'mpole: ma alloc failed'
        call dcopy ((nbf*nbf*slmax*3*(nat+1)),mp3c,1,dbl_mb(k_diff),1)
        call  daxpy((nbf*nbf*slmax*3*(nat+1)),
     &      -1.0d00,mp3fd,1,dbl_mb(k_diff),1)
        norm = ddot((nbf*nbf*slmax*3*(nat+1)),
     &      dbl_mb(k_diff),1,dbl_mb(k_diff),1)
        write(luout,*)' mpole difference norm:ir=',ir,
     &      ':imp_cent= ',imp_cent,': ',norm
c
        pderiv_compute_mpole = norm.lt.thresh
        if (.not.pderiv_compute_mpole) then
          do atom = 1,(nat+1)
            do ixyz = 1,3
              nintz = (atom-1)*3*slmax*nbf*nbf
              nintz = nintz + (ixyz-1)*slmax*nbf*nbf
              call prak3mat(dbl_mb(k_diff + nintz),
     &            (nbf*nbf*slmax),
     &            nbf,nbf,slmax,
     &            1,nbf,1,nbf,1,slmax,
     &            'diff matrix',atom,ixyz)
            enddo
          enddo
        endif
        if (.not.ma_free_heap(h_diff)) stop 'mpole: ma free failed '
      enddo
      enddo
      end
      subroutine prak3mat(buf,nbf,nbfi,nbfj,nbfk,
     &    ilo,ihi,jlo,jhi,klo,khi,msg,atom,ixyz)
      implicit none
      integer atom
      integer ixyz
      integer nbf, nbfi, nbfj, nbfk
      integer ilo,ihi,jlo,jhi,klo,khi
      character*(*) msg
      double precision buf(nbfk,nbfj,nbfi)
c
      integer i, j, k
c
      write(6,*)' '
      write(6,*)' prak3mat<',msg,'> atom=',atom,' ixyz=',ixyz
c
      if (nbf.ne.(nbfk*nbfj*nbfi)) then
        write(6,*)' nbf error '
        write(6,*)' nbf  = ',nbf
        write(6,*)' nbfi = ',nbfi
        write(6,*)' nbfj = ',nbfj
        write(6,*)' nbfk = ',nbfk
      endif
      do i = ilo,ihi
        do j = jlo,jhi
          do k = klo,khi
            if (abs(buf(k,j,i)).gt.1.0d-07) then
              write(6,10000)k,j,i,buf(k,j,i)
            endif
          enddo
        enddo
      enddo
      write(6,*)' '
10000 format(1x,'buf(',i3,',',i3,',',i3,')=',1pd20.10)
      end
      logical function pderiv_compute_2e3ct(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,eri3,eri3t,xyz)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
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
      double precision eri3t(nbf,nbf,nbf)
      double precision eri3(nbf,nbf,nbf)
      double precision xyz(3,nat)
c
      integer IR
      double precision thresh, norm
      double precision RJ(3), RK(3)
      integer nshell, cnt, nbfsh
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
      integer kshell, klo, khi, nbfshk, k
      integer nzero, nzerot, npos, npost
      integer h_diff, k_diff
c
      call dfill(3,0.0d00,RJ,1)
      call dfill(3,0.0d00,RK,1)
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*nbf),0.0d00,eri3t,1)
      call dfill((nbf*nbf*nbf),0.0d00,eri3,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,3
        if (IR.eq.1) then
          call dfill(3,0.0d00,RJ,1)
          call dfill(3,0.0d00,RK,1)
        elseif (IR.eq.2) then
          call dfill(3,0.0d00,RK,1)
          RJ(1) = 1.0d00
          RJ(2) = 2.0d00
          RJ(3) = 3.0d00
        elseif (IR.eq.3) then
          RJ(1) = 1.0d00
          RJ(2) = 2.0d00
          RJ(3) = 3.0d00
          RK(1) = 3.0d00
          RK(2) = 4.0d00
          RK(3) = 5.0d00
        else
          stop ' how did IR get here'
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 18'
        thresh = 1.0d-10
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        npos   = 0
        nzero  = 0
        nzerot = 0
        npost  = 0
        do ishell = 1, nshell ! basis functon NOT translated
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, nshell ! basis functon translated
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            do kshell = 1, nshell   ! fitting function
              if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &            stop 'cn2bfr error k'
              nbfshk = khi - klo + 1
              nbfsh = nbfshi*nbfshj*nbfshk
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call intp_2e3c(basis,kshell,basis,jshell,ishell,
     &              Rk,Rj,lscr,scr,lbuf,buf)

              cnt = 1
              do k = klo, khi
                do j = jlo, jhi
                  do i = ilo, ihi
*                    write(70,*)i,j,k,buf(cnt)
                    npos = npos + 1
                    if (abs(buf(cnt)).lt.thresh) nzero = nzero+1
                    eri3(i,j,k) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
              enddo

              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call intp_2e3ct(basis,ishell,jshell,basis,kshell,
     &            RJ,RK,lscr,scr,lbuf,buf)
*
              cnt = 1
              do i = ilo, ihi
                do j = jlo, jhi
                  do k = klo, khi
*                    write(71,*)i,j,k,buf(cnt)
                    npost = npost + 1
                    if (abs(buf(cnt)).lt.thresh) nzerot=nzerot+1
                    eri3t(i,j,k) = buf(cnt)
                    cnt = cnt + 1
                  enddo
                enddo
              enddo
c
            enddo
          enddo
        enddo
c
        if (.not.ma_alloc_get(mt_dbl,
     &      (nbf*nbf*nbf),
     &      'diff buffer',
     &      h_diff,
     &      k_diff)) stop '2e3ct: ma alloc failed'
        call dcopy((nbf*nbf*nbf),eri3,1,dbl_mb(k_diff),1)
        call daxpy((nbf*nbf*nbf),-1.0d00,eri3t,1,dbl_mb(k_diff),1)
        norm = ddot((nbf*nbf*nbf),dbl_mb(k_diff),1,dbl_mb(k_diff),1)
        write(luout,*)' 2e3ct info:'
        write(luout,*)' number possible in regular:',npos
        write(luout,*)' number zero     in regular:',nzero
        write(luout,*)' number non-zero in regular:',(npos-nzero)
        write(luout,*)' number possible in transpo:',npost
        write(luout,*)' number zero     in transpo:',nzerot
        write(luout,*)' number non-zero in transpo:',(npost-nzerot)
        write(luout,*)' 2e3ct difference norm:',ir,' ',norm
c
        pderiv_compute_2e3ct = norm.lt.thresh
        if (.not.pderiv_compute_2e3ct) then
          call prak2e3c(nbf,eri3,eri3t,dbl_mb(k_diff))
        endif
        if (.not.ma_free_heap(h_diff)) stop '2e3ct: ma free failed '
      enddo
      end
      subroutine prak2e3c(nbf,eri3,eri3t,diff)
      implicit none
      integer nbf
      double precision eri3 (nbf,nbf,nbf)
      double precision eri3t(nbf,nbf,nbf)
      double precision diff (nbf,nbf,nbf)
c
      integer i,j,k
c      
      do i = 1,nbf
        do j = 1,nbf
          do k = 1,nbf
            if (abs(diff(i,j,k)).gt.1.0d-07) then
              write(6,10000)i,j,k,
     &            eri3(i,j,k),
     &            eri3t(i,j,k),
     &            diff(i,j,k)
            endif
          enddo
        enddo
      enddo
10000 format(1x,'<',i3,',',i3,',',i3,'> =',3(1pd20.10))
      end
      logical function pderiv_compute_d2e3ct(
     &    geom,basis,nbf,nat,lbuf,lscr,
     &    buf,bufp,bufm,scr,eri3,eri3t,xyz)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
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
      double precision eri3t(nbf,nbf,nbf,3,nat)
      double precision eri3(nbf,nbf,nbf,3,nat)
      double precision xyz(3,nat)
c
      integer nat3, nintz
      integer atom, ixyz, IR
      double precision thresh, norm
      double precision RJ(3), RK(3)
      integer nshell, cnt, nbfsh, watom(4)
      integer ishell, ilo, ihi, nbfshi, i
      integer jshell, jlo, jhi, nbfshj, j
      integer kshell, klo, khi, nbfshk, k
      integer nzero, nzerot, npos, npost
      integer h_diff, k_diff
c
      nat3 = 3*nat
      call dfill(3,0.0d00,RJ,1)
      call dfill(3,0.0d00,RK,1)
      call dfill((nat3),0.0d00,xyz,1)
      call dfill((nbf*nbf*nbf*nat3),0.0d00,eri3t,1)
      call dfill((nbf*nbf*nbf*nat3),0.0d00,eri3,1)
      call dfill(lbuf,0.0d00,buf,1)
      call dfill(lbuf,0.0d00,bufp,1)
      call dfill(lbuf,0.0d00,bufm,1)
      call dfill(lscr,0.0d00,scr,1)
      do IR = 1,3
        if (IR.eq.1) then
          call dfill(3,0.0d00,RJ,1)
          call dfill(3,0.0d00,RK,1)
        elseif (IR.eq.2) then
          call dfill(3,0.0d00,RK,1)
          RJ(1) = 1.0d00
          RJ(2) = 2.0d00
          RJ(3) = 3.0d00
        elseif (IR.eq.3) then
          RJ(1) = 1.0d00
          RJ(2) = 2.0d00
          RJ(3) = 3.0d00
          RK(1) = 3.0d00
          RK(2) = 4.0d00
          RK(3) = 5.0d00
        else
          stop ' how did IR get here'
        endif
        if (.not.ma_verify_allocator_stuff()) stop ' ma broke 18'
        thresh = 1.0d-10
*
        if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont'
        
        npos   = 0
        nzero  = 0
        nzerot = 0
        npost  = 0
        do ishell = 1, nshell ! basis functon NOT translated
          if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &        stop 'cn2bfr error i'
          nbfshi = ihi - ilo + 1
          do jshell = 1, nshell ! basis functon translated
            if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &          stop 'cn2bfr error j'
            nbfshj = jhi - jlo + 1
            do kshell = 1, nshell   ! fitting function
              if (.not.bas_cn2bfr(basis,kshell,klo,khi))
     &            stop 'cn2bfr error k'
              nbfshk = khi - klo + 1
              nbfsh = nbfshi*nbfshj*nbfshk
              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call intpd_2e3c(basis,kshell,basis,jshell,ishell,
     &              Rk,Rj,lscr,scr,lbuf,buf,watom)

              nintz = nbfsh
              do atom = 1,4
                if (watom(atom).gt.0) then
                  cnt = ((atom-1)*nintz*3) + 1
                  do ixyz = 1,3
                    do k = klo, khi
                      do j = jlo, jhi
                        do i = ilo, ihi
                          npos = npos + 1
                          if (abs(buf(cnt)).lt.thresh) nzero = nzero+1
*                          write(70,*)i,j,k,ixyz,watom(atom),buf(cnt)
                          eri3(i,j,k,ixyz,watom(atom)) = buf(cnt)
                          cnt = cnt + 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
              enddo

              call dfill(lscr,0.0d00,scr,1)
              call dfill(lbuf,0.0d00,buf,1)
              call intpd_2e3ct(basis,ishell,jshell,basis,kshell,
     &            RJ,RK,lscr,scr,lbuf,buf,watom)
*
              do atom = 1,4
                if (watom(atom).gt.0) then
                  cnt = ((atom-1)*nintz*3) + 1
                  do ixyz = 1,3
                    do i = ilo, ihi
                      do j = jlo, jhi
                        do k = klo, khi
                          npost = npost + 1
                          if (abs(buf(cnt)).lt.thresh) nzerot=nzerot+1
*                          write(71,*)i,j,k,ixyz,watom(atom),buf(cnt)
                          eri3t(i,j,k,ixyz,watom(atom)) = buf(cnt)
                          cnt = cnt + 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
              enddo
c
            enddo
          enddo
        enddo
c
        if (.not.ma_alloc_get(mt_dbl,
     &      (nbf*nbf*nbf*nat3),
     &      'diff buffer',
     &      h_diff,
     &      k_diff)) stop 'd2e3ct: ma alloc failed'
        call dcopy((nbf*nbf*nbf*nat3),eri3,1,dbl_mb(k_diff),1)
        call daxpy((nbf*nbf*nbf*nat3),-1.0d00,eri3t,1,dbl_mb(k_diff),1)
        norm =
     &      ddot((nbf*nbf*nbf*nat3),dbl_mb(k_diff),1,dbl_mb(k_diff),1)
        write(luout,*)' d2e3ct info:'
        write(luout,*)' number possible in regular:',npos
        write(luout,*)' number zero     in regular:',nzero
        write(luout,*)' number non-zero in regular:',(npos-nzero)
        write(luout,*)' number possible in transpo:',npost
        write(luout,*)' number zero     in transpo:',nzerot
        write(luout,*)' number non-zero in transpo:',(npost-nzerot)
        write(luout,*)' d2e3ct difference norm:',ir,' ',norm
c
        pderiv_compute_d2e3ct = norm.lt.thresh
        if (.not.pderiv_compute_d2e3ct) then
          call prakd2e3c(nbf,nat,eri3,eri3t,dbl_mb(k_diff))
        endif
        if (.not.ma_free_heap(h_diff)) stop 'd2e3ct: ma free failed '
      enddo
      end
      subroutine prakd2e3c(nbf,nat,eri3,eri3t,diff)
      implicit none
      integer nbf
      integer nat
      double precision eri3(nbf,nbf,nbf,3,nat)
      double precision eri3t(nbf,nbf,nbf,3,nat)
      double precision diff(nbf,nbf,nbf,3,nat)
c
      integer atom,ixyz,i,j,k
c
      do atom = 1,nat
        do ixyz = 1,3
          do i = 1,nbf
            do j = 1,nbf
              do k = 1,nbf
                if (abs(diff(i,j,k,ixyz,atom)).gt.1.0d-07) then
                  write(6,10000)i,j,k,ixyz,atom,
     &            eri3(i,j,k,ixyz,atom),
     &            eri3t(i,j,k,ixyz,atom),
     &            diff(i,j,k,ixyz,atom)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'<',i3,',',i3,',',i3,',',i3,',',i3,'> =',3(1pd20.10))
      end
