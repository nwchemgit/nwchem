C> \ingroup task
C> @{
      logical function task_ecp_print_integrals(rtdb)
      implicit none
#include "errquit.fh"
* $Id$
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
c
      logical int_normalize, raktask26_a
      external int_normalize, raktask26_a
c
c
      integer rtdb
c
      logical status
      integer basis, geom
      integer nbf, nat, nshell
      integer maxg1, maxs1
      integer hbuf, hscr, hint
      integer kbuf, kscr, kint
c
      task_ecp_print_integrals = .false.
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
      call int_init(rtdb, 1, basis)
      call int_mem_1e(maxg1, maxs1)
      write(luout,*)' maxg1 = ',maxg1
      write(luout,*)' maxs1 = ',maxs1
      maxs1 = max(maxs1,(nbf*nbf))
      write(luout,*)' maxs1 = ',maxs1, ' after max for copy '
      status = .true.
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxg1,'int buffer' ,hbuf,kbuf)
      status = status .and.
     &    ma_alloc_get(mt_dbl,maxs1,'scr buffer' ,hscr,kscr)
      status = status .and.
     &    ma_alloc_get(mt_dbl,(nbf*nbf),'ints' ,hint,kint)
c
      task_ecp_print_integrals = raktask26_a(rtdb,
     &    geom, basis, nbf, nat, nshell, maxg1, maxs1,
     &    dbl_mb(kint), 
     &    dbl_mb(kbuf),
     &    dbl_mb(kscr))
      status = .true.
      status = status.and.ma_free_heap(hint)
      status = status.and.ma_free_heap(hbuf)
      status = status.and.ma_free_heap(hscr)
      status = status.and.bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      task_ecp_print_integrals = task_ecp_print_integrals.and.status
      call int_terminate()
      end
C> @}
      logical function raktask26_a(rtdb,
     &    geom, basis, nbf, nat, nshell, szb, szs,
     &    zint,buf,scr)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "geomP.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "bas.fh"
#include "stdio.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      integer rtdb, geom, basis, nbf, nat, nshell, szb, szs
      double precision zint(nbf,nbf)
      double precision buf(szb)
      double precision scr(szs)
c
      double precision val
      integer nshell_use
      integer ishell, ilo, ihi, nbfshi
      integer jshell, jlo, jhi, nbfshj
      integer ii_np, ii_gen, ii_exp, ii_cf, ii_type, ii_atom
      integer jj_np, jj_gen, jj_exp, jj_cf, jj_type, jj_atom
      integer nbfsh, ucont, xbas, cnt, i, j, lu
      logical does_it_exist
c
      character*255 filename
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      call dfill (nbf*nbf,0.0d00,zint,1)
      call dfill (szb,0.0d00,buf,1)
      call dfill (szs,0.0d00,scr,1)
c
      xbas = basis + BASIS_HANDLE_OFFSET
      nshell_use = nshell
c
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
*          
          call dfill (szb,0.0d00,buf,1)
          call dfill (szs,0.0d00,scr,1)
          call int_ecp_hf1(
     &        coords(1,ii_atom,geom),
     &        dbl_mb(mb_exndcf(ii_exp,xbas)),
     &        dbl_mb(mb_exndcf(ii_cf,xbas)),
     &        ii_np, ii_gen, ii_type, ii_atom,
          
     &        coords(1,jj_atom,geom),
     &        dbl_mb(mb_exndcf(jj_exp,xbas)),
     &        dbl_mb(mb_exndcf(jj_cf,xbas)),
     &        jj_np, jj_gen, jj_type, jj_atom,
          
     &        buf,nbfsh,scr,szs,.false.)
*--------          
          cnt = 1
          do i = ilo,ihi
            do j = jlo, jhi
              zint(i,j) = buf(cnt)
              zint(j,i) = buf(cnt)
            enddo
          enddo
*--------          
        enddo
      enddo
* now have full ecp integral matrix print it out      
c
      call util_file_name('ecp_integrals',.false.,.false.,filename)
c
      lu = 69
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if (does_it_exist) then
        write(luout,*)
     &      'rak26: overwrite of existing file',
     &      filename(1:inp_strlen(filename))
        call util_file_unlink(filename)
      endif
*
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if (does_it_exist) then
        write(luout,*)
     &      'rak26: file not unlinked:  ',
     &      filename(1:inp_strlen(filename))
        call errquit('rak26: fatal error ',911, DISK_ERR)
      endif
*
      open(unit=lu,file=filename,
     &    form='formatted',
     &    access='sequential',
     &    status='new')
      cnt = 0
      do i = 1,nbf
        do j = 1,i
          cnt = cnt + 1
          val = zint(i,j)
          if (abs(val).gt.1.0d-8) then
            write(lu,10000)cnt,i,j,val
          endif
        enddo
      enddo
      close(unit=lu,status='keep')
      raktask26_a = .true.
10000 format(1x,i10,i5,i5,1pd20.10)
      end


