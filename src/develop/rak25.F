C> \ingroup task
C> @{
      logical function task_ecp_deriv_check(rtdb)
      implicit none
#include "errquit.fh"
* $Id$
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
c
      logical int_normalize, raktask25_a
      external int_normalize, raktask25_a
c
      integer rtdb
c
      logical status
      integer basis, geom
      integer nbf, nat, nshell
      integer maxg1, maxs1
      integer hbuf, hbufp, hbufm, hscr, hg, hfd, hxyz
      integer kbuf, kbufp, kbufm, kscr, kg, kfd, kxyz
      integer hbufpp, hbufmm, hbuf3
      integer kbufpp, kbufmm, kbuf3
c
      task_ecp_deriv_check = .false.
c
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
      if (.not.bas_numcont(basis,nshell)) call errquit
     &    ('numcont failed',911, BASIS_ERR)
c
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      write(6,*) 'number of atoms ', nat
c
      call intd_init(rtdb, 1, basis)
      call int_mem_1e(maxg1, maxs1)
      write(luout,*)' maxg1 = ',maxg1
      write(luout,*)' maxs1 = ',maxs1
      maxs1 = max(maxs1,(nbf*nbf*3*nat))
      write(luout,*)' maxs1 = ',maxs1, ' after max for copy '
      status = .true.
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer' ,hbuf,kbuf)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer' ,hbuf3,kbuf3)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer(+)' ,hbufp,kbufp)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer(2+)' ,hbufpp,kbufpp)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer(-)' ,hbufm,kbufm)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer(2-)' ,hbufmm,kbufmm)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxs1,'scr buffer' ,hscr,kscr)
      status = status .and.
     &    ma_alloc_get(mt_dbl,(nbf*nbf*3*nat),'grad' ,hg,kg)
      status = status .and.
     &    ma_alloc_get(mt_dbl,(nbf*nbf*3*nat),
     &    'finite difference grad' ,hfd,kfd)
      status = status .and.
     &    ma_alloc_get(mt_dbl,3*nat,'my coords',hxyz,kxyz)
      if (.not.status) stop ' memory failed rak25 '
      call intd_terminate()
      task_ecp_deriv_check = raktask25_a(rtdb,
     &    geom, basis, nbf, nat, nshell, maxg1, maxs1,
     &    dbl_mb(kg), dbl_mb(kfd), 
     &    dbl_mb(kbuf),
     &    dbl_mb(kbuf3),
     &    dbl_mb(kscr),
     &    dbl_mb(kbufp),
     &    dbl_mb(kbufm),
     &    dbl_mb(kbufpp),
     &    dbl_mb(kbufmm),
     &    dbl_mb(kxyz))
      status = .true.
      status = status.and.ma_free_heap(hg)
      status = status.and.ma_free_heap(hfd)
      status = status.and.ma_free_heap(hbuf)
      status = status.and.ma_free_heap(hbuf3)
      status = status.and.ma_free_heap(hscr)
      status = status.and.ma_free_heap(hbufp)
      status = status.and.ma_free_heap(hbufm)
      status = status.and.ma_free_heap(hbufpp)
      status = status.and.ma_free_heap(hbufmm)
      status = status.and.ma_free_heap(hxyz)
      status = status.and.bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      task_ecp_deriv_check = task_ecp_deriv_check.and.status
      end
C> @}
      logical function raktask25_a(rtdb,
     &    geom, basis, nbf, nat, nshell,
     &    sizeg, sizescr,
     &    grad, fdgrad, buf, buf3, scr, bufp, bufm,
     &    bufpp, bufmm, xyz)
      implicit none
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "geomP.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "bas.fh"
#include "stdio.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      logical int_ecp_init
      external int_ecp_init
      double precision ddot
      external ddot
c
      logical ignore
      integer rtdb, geom, basis, nbf, nat, nshell
      integer sizeg, sizescr, ecpid, soid
      double precision xyz(3,nat)
      double precision grad(nbf,nbf,3,nat)
      double precision fdgrad(nbf,nbf,3,nat)
      double precision buf(sizeg), buf3(sizeg)
      double precision bufp(sizeg), bufm(sizeg)
      double precision bufpp(sizeg), bufmm(sizeg)
      double precision scr(sizescr)
c
      double precision thresh, delta, factor, norm
      double precision asum, fdsum
      integer xbas, xsize
      integer ii_np, ii_gen, ii_exp, ii_cf, ii_type, ii_atom
      integer jj_np, jj_gen, jj_exp, jj_cf, jj_type, jj_atom
      integer ishell, ilo, ihi, nbfshi
      integer jshell, jlo, jhi, nbfshj
      integer nbfsh, ucont
      integer atom, ixyz, cnt, i, j
      integer nat3, nshell_use
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      call intd_init(rtdb,1,basis)
c
      nat3 = 3*nat
      call dfill((3*nat),0.0d00,xyz,1)
      call dfill((nbf*nbf*3*nat),0.0d00,grad,1)
      call dfill((nbf*nbf*3*nat),0.0d00,fdgrad,1)
      call dfill(sizeg,0.0d00,buf,1)
      call dfill(sizeg,0.0d00,bufp,1)
      call dfill(sizeg,0.0d00,bufm,1)
      call dfill(sizeg,0.0d00,bufpp,1)
      call dfill(sizeg,0.0d00,bufmm,1)
      call dfill(sizescr,0.0d00,scr,1)
* store original coordintates
*      write(6,*)' original coordinates '
*      if (.not.geom_print(geom)) stop ' gp error'
      call dcopy(nat3,coords(1,1,geom),1,xyz,1)
c
      ignore = bas_get_ecp_handle(basis,ecpid)
      write(6,*)' ecp id ',ecpid
      soid = 0
      thresh = 1.0d-09
      delta = 0.0025d00
      write(6,*)' delta =',delta
      xbas = basis + BASIS_HANDLE_OFFSET
c
      nshell_use = nshell
*      nshell_use = 1
      do ishell = 1, nshell_use
        write(6,*)' fd: ishell = ',ishell,' of ',nshell_use
        call util_flush(6)
        if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &      stop 'cn2bfr error i'
        nbfshi = ihi - ilo + 1
        ucont = (sf_ibs_cn2ucn(ishell,xbas))
        ii_np   = infbs_cont(CONT_NPRIM,ucont,xbas)
        ii_gen  = infbs_cont(CONT_NGEN,ucont,xbas)
        ii_exp  = infbs_cont(CONT_IEXP,ucont,xbas)
        ii_cf   = infbs_cont(CONT_ICFP,ucont,xbas)
        ii_type = infbs_cont(CONT_TYPE,ucont,xbas)
        ii_atom = (sf_ibs_cn2ce(ishell,xbas))
        do jshell = 1, ishell
          if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &        stop 'cn2bfr error j'
          nbfshj = jhi - jlo + 1
          nbfsh = nbfshi*nbfshj
*          write(6,*)' fd:     jshell = ',jshell,' size =',nbfsh
          ucont = (sf_ibs_cn2ucn(jshell,xbas))
          jj_np   = infbs_cont(CONT_NPRIM,ucont,xbas)
          jj_gen  = infbs_cont(CONT_NGEN,ucont,xbas)
          jj_exp  = infbs_cont(CONT_IEXP,ucont,xbas)
          jj_cf   = infbs_cont(CONT_ICFP,ucont,xbas)
          jj_type = infbs_cont(CONT_TYPE,ucont,xbas)
          jj_atom = (sf_ibs_cn2ce(jshell,xbas))
          do atom = 1,nat
            do ixyz = 1,3
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) =
     &              coords(ixyz,atom,geom) + delta + delta
              call dfill(sizescr,0.0d00,scr,1)
              call dfill(sizeg,0.0d00,buf,1)
              call dfill(sizeg,0.0d00,bufpp,1)
c
*             write(6,*)' plus ', atom, ixyz
*              if (.not.geom_print(geom)) stop ' gp error'
              call int_ecp_terminate()
              ignore = int_ecp_init(ecpid, soid,1)
              call int_ecp_hf1(
     &            coords(1,ii_atom,geom),
     &            dbl_mb(mb_exndcf(ii_exp,xbas)),
     &            dbl_mb(mb_exndcf(ii_cf,xbas)),
     &            ii_np, ii_gen, ii_type, ii_atom,

     &            coords(1,jj_atom,geom),
     &            dbl_mb(mb_exndcf(jj_exp,xbas)),
     &            dbl_mb(mb_exndcf(jj_cf,xbas)),
     &            jj_np, jj_gen, jj_type, jj_atom,

     &            bufpp,nbfsh,scr,sizescr,.false.)
*              if (nbfsh.eq.1) write(6,*)
*     &            ' ',atom,' X=',ixyz,' +', bufp(1)
*
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) = coords(ixyz,atom,geom) + delta
              call dfill(sizescr,0.0d00,scr,1)
              call dfill(sizeg,0.0d00,bufp,1)
c
*             write(6,*)' plus ', atom, ixyz
*              if (.not.geom_print(geom)) stop ' gp error'
              call int_ecp_terminate()
              ignore = int_ecp_init(ecpid, soid,1)
              call int_ecp_hf1(
     &            coords(1,ii_atom,geom),
     &            dbl_mb(mb_exndcf(ii_exp,xbas)),
     &            dbl_mb(mb_exndcf(ii_cf,xbas)),
     &            ii_np, ii_gen, ii_type, ii_atom,

     &            coords(1,jj_atom,geom),
     &            dbl_mb(mb_exndcf(jj_exp,xbas)),
     &            dbl_mb(mb_exndcf(jj_cf,xbas)),
     &            jj_np, jj_gen, jj_type, jj_atom,

     &            bufp,nbfsh,scr,sizescr,.false.)
*              if (nbfsh.eq.1) write(6,*)
*     &            ' ',atom,' X=',ixyz,' +', bufp(1)
*
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) = coords(ixyz,atom,geom) - delta
              call dfill(sizescr,0.0d00,scr,1)
              call dfill(sizeg,0.0d00,bufm,1)
c
*              write(6,*)' minus ', atom, ixyz
*              if (.not.geom_print(geom)) stop ' gp error'
              call int_ecp_terminate()
              ignore = int_ecp_init(ecpid, soid,1)
              call int_ecp_hf1(
     &            coords(1,ii_atom,geom),
     &            dbl_mb(mb_exndcf(ii_exp,xbas)),
     &            dbl_mb(mb_exndcf(ii_cf,xbas)),
     &            ii_np, ii_gen, ii_type, ii_atom,

     &            coords(1,jj_atom,geom),
     &            dbl_mb(mb_exndcf(jj_exp,xbas)),
     &            dbl_mb(mb_exndcf(jj_cf,xbas)),
     &            jj_np, jj_gen, jj_type, jj_atom,

     &            bufm,nbfsh,scr,sizescr,.false.)
*
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              coords(ixyz,atom,geom) =
     &              coords(ixyz,atom,geom) - delta - delta 
              call dfill(sizescr,0.0d00,scr,1)
              call dfill(sizeg,0.0d00,bufmm,1)
c
*              write(6,*)' minus ', atom, ixyz
*              if (.not.geom_print(geom)) stop ' gp error'
              call int_ecp_terminate()
              ignore = int_ecp_init(ecpid, soid,1)
              call int_ecp_hf1(
     &            coords(1,ii_atom,geom),
     &            dbl_mb(mb_exndcf(ii_exp,xbas)),
     &            dbl_mb(mb_exndcf(ii_cf,xbas)),
     &            ii_np, ii_gen, ii_type, ii_atom,

     &            coords(1,jj_atom,geom),
     &            dbl_mb(mb_exndcf(jj_exp,xbas)),
     &            dbl_mb(mb_exndcf(jj_cf,xbas)),
     &            jj_np, jj_gen, jj_type, jj_atom,

     &            bufmm,nbfsh,scr,sizescr,.false.)
*              if (nbfsh.eq.1) write(6,*)
*     &            ' ',atom,' X=',ixyz,' -', bufm(1)
*
*.............. 3 point
              factor = 1.0d00/(2.0d00*delta)
              call dcopy(nbfsh,bufp,1,buf3,1)
              call daxpy(nbfsh,-1.0d00,bufm,1,buf3,1)
              call dscal(nbfsh,factor,buf3,1)
*.............. 5 point
              factor = 2.0d0/(3.0d0*delta)
              call dscal(nbfsh,factor,bufp,1)
              call dscal(nbfsh,factor,bufm,1)
              factor = 1.0d00/(12.0d00*delta)
              call dscal(nbfsh,factor,bufpp,1)
              call dscal(nbfsh,factor,bufmm,1)
              call dcopy(nbfsh,          bufp, 1,buf,1)
              call daxpy(nbfsh,-1.0d00,  bufm, 1,buf,1)
              call daxpy(nbfsh,-1.0d00, bufpp, 1,buf,1)
              call daxpy(nbfsh, 1.0d00, bufmm, 1,buf,1)
              
*              if (nbfsh.eq.1) write(6,*)
*     &            ' ',atom,' X=',ixyz,' g', buf(1)
              cnt = 1
              do i = ilo,ihi
                do j = jlo, jhi
                  fdgrad(i,j,ixyz,atom) = buf(cnt)
                  fdgrad(j,i,ixyz,atom) = buf(cnt)
                    grad(i,j,ixyz,atom) = buf3(cnt)   ! grad has 3 point
                    grad(j,i,ixyz,atom) = buf3(cnt)
                  cnt = cnt + 1
                enddo
              enddo
              if (.not.ma_verify_allocator_stuff())
     &            stop ' ma broke 3'
            enddo
          enddo
        enddo
      enddo
      write(6,*)' differences between 3 point and 5 point '
      call rak25_diff(fdgrad,grad,nbf,nat)
c
*      call prt25(fdgrad,nbf,nat,'fd gradient')
c
      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
      call dfill((nbf*nbf*3*nat),0.0d00,grad,1)
      call int_ecp_terminate()
      ignore = int_ecp_init(ecpid, soid,1)
*
      nshell_use = nshell
*      nshell_use = 1
      do ishell = 1, nshell_use
        write(6,*)'  a: ishell = ',ishell,' of ',nshell_use
        call util_flush(6)
        if (.not.bas_cn2bfr(basis,ishell,ilo,ihi))
     &      stop 'cn2bfr error i'
        nbfshi = ihi - ilo + 1
        ucont = (sf_ibs_cn2ucn(ishell,xbas))
        ii_np   = infbs_cont(CONT_NPRIM,ucont,xbas)
        ii_gen  = infbs_cont(CONT_NGEN,ucont,xbas)
        ii_exp  = infbs_cont(CONT_IEXP,ucont,xbas)
        ii_cf   = infbs_cont(CONT_ICFP,ucont,xbas)
        ii_type = infbs_cont(CONT_TYPE,ucont,xbas)
        ii_atom = (sf_ibs_cn2ce(ishell,xbas))
        do jshell = 1, ishell
          if (.not.bas_cn2bfr(basis,jshell,jlo,jhi))
     &        stop 'cn2bfr error j'
          nbfshj = jhi - jlo + 1
          nbfsh = nbfshi*nbfshj
*          write(6,*)'  a:     jshell = ',jshell,' size =',nbfsh
          ucont = (sf_ibs_cn2ucn(jshell,xbas))
          jj_np   = infbs_cont(CONT_NPRIM,ucont,xbas)
          jj_gen  = infbs_cont(CONT_NGEN,ucont,xbas)
          jj_exp  = infbs_cont(CONT_IEXP,ucont,xbas)
          jj_cf   = infbs_cont(CONT_ICFP,ucont,xbas)
          jj_type = infbs_cont(CONT_TYPE,ucont,xbas)
          jj_atom = (sf_ibs_cn2ce(jshell,xbas))
c
          call dfill(sizescr,0.0d00,scr,1)
          call dfill(sizeg,0.0d00,buf,1)
c
          call intd_ecp_hf1(
     &        coords(1,ii_atom,geom),
     &        dbl_mb(mb_exndcf(ii_exp,xbas)),
     &        dbl_mb(mb_exndcf(ii_cf,xbas)),
     &        ii_np, ii_gen, ii_type, ii_atom,
          
     &        coords(1,jj_atom,geom),
     &        dbl_mb(mb_exndcf(jj_exp,xbas)),
     &        dbl_mb(mb_exndcf(jj_cf,xbas)),
     &        jj_np, jj_gen, jj_type, jj_atom,
          
     &        buf,nbfsh,nat,scr,sizescr,.false.)
          xsize = nbfsh*3*nat
*          if (nbfsh.eq.1) then
*            do i = 1, 3*nat
*              write(6,*)' grad buff ',i,' ',buf(i)
*            enddo
*          endif
          cnt = 1
          do atom=1,nat
            do ixyz = 1,3
              do i = ilo, ihi
                do j = jlo, jhi
                  grad(i,j,ixyz,atom) = buf(cnt)
                  grad(j,i,ixyz,atom) = buf(cnt)
                  cnt = cnt + 1
                enddo
              enddo
            enddo
          enddo
c
        enddo
      enddo

c
*      call prt25(grad,nbf,nat,'   gradient')
c
      xsize = nbf*nbf*3*nat
      do i = 1,nbf
        do j = 1,i
          asum = 0.0d00
          fdsum = 0.0d00
          do atom = 1,nat
            do ixyz = 1,3
              asum = asum + grad(i,j,ixyz,atom)
              fdsum = fdsum + fdgrad(i,j,ixyz,atom)
            enddo
          enddo
          if ((abs(asum).gt.1.0d-06).or.(abs(fdsum).gt.1.0d-06))
     &        write(6,*)i,j,asum,fdsum
        enddo
      enddo
c
      xsize = nbf*nbf*3*nat
      call dcopy(xsize,grad,1,scr,1)
      call daxpy(xsize,-1.0d00,fdgrad,1,scr,1)
      norm = ddot(xsize,scr,1,scr,1)
      write(luout,*)' difference norm for ecp derivs: ',norm
      call util_flush(luout)
c
      call intd_terminate()
      raktask25_a = norm .lt. thresh
      if (raktask25_a) return
c
      call rak25_diff(fdgrad,grad,nbf,nat)
c
      end
      subroutine prt25(g,nbf,nat,msg)
      implicit none
#include "stdio.fh"
      integer nbf
      integer nat
      double precision g(nbf,nbf,3,nat)
      character*(*) msg
c
      double precision thresh, val
      integer i,j,x,a
c
      thresh = 1.0d-06
      write(luout,*)'........................................   ',msg
      do i = 1,nbf
        do j = 1,i
          do x = 1,3
            do a = 1,nat
              val = g(i,j,x,a)
              if (abs(val).gt.thresh) then
                write(luout,10000)i,j,x,a,val
              endif
            enddo
          enddo
        enddo
      enddo
      write(luout,*)
      write(luout,*)
10000 format(1x,'g(',4i5,') =',1pd20.10)
      end
      subroutine rak25_diff(fd,g,nbf,nat)
      implicit none
#include "stdio.fh"
      integer nbf, nat
      double precision fd(nbf,nbf,3,nat)
      double precision  g(nbf,nbf,3,nat)
c
      double precision aval, val, thresh
      integer i, j, x, a, t
      integer hgram(10)
      double precision tt(9)
      logical pheader
c
      call ifill(10,0,hgram,1)
      do i = 1,9
        tt(i) = 10.0d00**(-(12-i))
      enddo
*      write(6,*)' tt is ', tt
      write(luout,*)
      write(luout,*)
      write(luout,*)
      pheader = .false.
      thresh = 1.0d-05
      do i = 1,nbf
        do j = 1,nbf
          do x = 1,3
            do a = 1,nat
              val = g(i,j,x,a) - fd(i,j,x,a) 
              aval = abs(val)
              if (aval.lt.tt(1)) then
                hgram(1) = hgram(1) + 1
              else if (aval.gt.tt(9)) then
                hgram(10) = hgram(10) + 1
              else
                do t = 2,9
                  if (aval.le.tt(t).and.aval.gt.tt((t-1)))
     &                hgram(t) = hgram(t) + 1
                enddo
              endif
              if (aval.gt.thresh)then
                if (.not.pheader) then
                  write(luout,10000)
                  pheader = .true.
                endif
                write(luout,10001)
     &              i,j,x,a,g(i,j,x,a),fd(i,j,x,a),val
              endif
            enddo
          enddo
        enddo
      enddo
c
10010 format(1x,'Difference count <',1x,1pd8.2,18x,i10,' values')
10011 format(1x,'Difference count <',1x,
     &    1pd8.2,' and >',1x,1pd8.2,3x,i10,' values')
10012 format(1x,'Difference count ',15x,'>',1x,1pd8.2,3x,i10,' values')
      write(luout,10010)tt(1),hgram(1) 
      do i = 2,9
        write(luout,10011) tt(i),tt(i-1),hgram(i) 
      enddo
      write(luout,10012)tt(9),hgram(10)
*                12345
10000 format(1x,'   i ','   j ','   x ','   a ',3x,
*          12345678901234567890123
     &    '   ---- analytic ------',3x,
     &    '   - FiniteDifference -',3x,
     &    '   --- difference -----')
10001 format(1x,4i5,3(3x,1pd20.10))
      end

