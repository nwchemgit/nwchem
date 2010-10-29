*
* $Id$
*
      logical function raktask_intdd_3c(rtdb)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*::functions
      logical int_normalize 
      external int_normalize
*::passed
      integer rtdb
*::local
      integer basis, geom, nbf, cn_nbf_max, nshell
      integer nat, nat3
      integer size, sizesq, sizeg, scr_size
      integer maxg1, maxg2, maxs1, maxs2
      logical status
      integer hbuf, kbuf, hscr, kscr
      integer hfd, kfd, hfdsq, kfdsq, hxyz, kxyz
      integer hgradp, kgradp
      integer hgradm, kgradm
      integer hbufsum, kbufsum
*
      raktask_intdd_3c = .false.
*
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
      write(luout,*)' geom/basis loaded'
      call util_flush(luout)
c
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
c
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0, BASIS_ERR)
      if (.not. gbs_map_print(basis))
     $    call errquit(' gbs_map_print failed', 0, UNKNOWN_ERR)
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &    ('numbf failed',911, BASIS_ERR)
c
      if (.not.bas_numcont(basis,nshell)) call errquit
     &    ('numbf failed',911, BASIS_ERR)
c
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      write(luout,*) 'number of atoms ', nat
      nat3 = 3*nat
c
      if (.not.bas_nbf_cn_max(basis,cn_nbf_max))
     &    stop 'bas_nbf_cn_max'
c
      size   = (cn_nbf_max**4)*78
      sizesq = (cn_nbf_max**4)*12*12
      sizeg  = (cn_nbf_max**4)*12
      size   = size   + size  /10
      sizesq = sizesq + sizesq/10
      sizeg  = sizeg  + sizeg /10
c
      call intdd_init(rtdb,1,basis)
      call int_mem_print()
      call int_mem_1e(maxg1,maxs1)
      call int_mem_2e4c(maxg2,maxs2)
*
      write(luout,*)' maxg1 :',maxg1
      write(luout,*)' maxs1 :',maxs1
      write(luout,*)' maxg2 :',maxg2
      write(luout,*)' size  :',size
      write(luout,*)' sizesq:',sizesq
      write(luout,*)' sizeg :',sizeg
      write(luout,*)' maxs2 :',maxs2
      call util_flush(luout)
*
      scr_size = 2*maxs2
      status =
     &    ma_alloc_get(mt_dbl,size,'intdd buffer',hbuf,kbuf)
      status = status.and.
     &    ma_alloc_get(mt_dbl,sizesq,'intdd buffer summed',
     &    hbufsum,kbufsum)      
      status = status.and.
     &    ma_alloc_get(mt_dbl,scr_size,'scr buffer',hscr,kscr)      
      status = status.and.
     &    ma_alloc_get(mt_dbl,size,'intdd fd buffer',hfd,kfd)
      status = status.and.
     &    ma_alloc_get(mt_dbl,sizesq,'intdd fd sq buffer',hfdsq,kfdsq)
      status = status.and.
     &    ma_alloc_get(mt_dbl,3*nat,'coords',hxyz,kxyz)
      status = status.and.
     &    ma_alloc_get(mt_dbl,sizeg,'grad +',hgradp,kgradp)
      status = status.and.
     &    ma_alloc_get(mt_dbl,sizeg,'grad -',hgradm,kgradm)
      if (.not.status) stop ' memory alloc failed rak23 (1)'
*
      call raktask_intdd_3ca(geom,basis,nbf,nshell,cn_nbf_max,
     &    size,dbl_mb(kbuf),
     &    scr_size,dbl_mb(kscr),
     &    dbl_mb(kgradp),
     &    dbl_mb(kgradm),
     &    dbl_mb(kfd),
     &    dbl_mb(kfdsq),
     &    nat,dbl_mb(kxyz),
     &    size,sizesq,sizeg,
     &    dbl_mb(kbufsum))
*
      call intdd_terminate()
      raktask_intdd_3c = bas_destroy(basis)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    geom_destroy(geom)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hscr)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hbuf)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hfd)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hfdsq)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hxyz)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hgradp)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hgradm)
      raktask_intdd_3c = raktask_intdd_3c.and.
     &    ma_free_heap(hbufsum)
      end
      subroutine raktask_intdd_3ca(geom,basis,nbf,nshell,cn_nbf_max,
     &    lbuf,buf,lscr,scr,gradp, gradm,buffd,buffdsq,nat,xyz,
     &    lbuffd, lbuffdsq, lgrad, bufsum)
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "inp.fh"
c::functions
      logical rakdd_checktrans
      external rakdd_checktrans
c::passed
      integer geom
      integer basis
      integer nbf
      integer nshell
      integer cn_nbf_max
      integer lbuf
      integer lscr
      integer nat
      integer lbuffd, lbuffdsq, lgrad
      double precision buffd(lbuffd), buffdsq(lbuffdsq)
      double precision gradp(lgrad), gradm(lgrad)
      double precision buf(lbuf)
      double precision bufsum(lbuffdsq)
      double precision scr(lscr)
      double precision xyz(3,nat)
*
      integer hbufcp, kbufcp
      integer nzero,ncount,count
      integer iatom, jatom, katom, latom
      integer ish, jsh, ksh
      integer ilo, ihi, inbf
      integer jlo, jhi, jnbf
      integer klo, khi, knbf
      integer nint, ninth, nintg
      integer idatom(4)
      integer idatoms(4)
      integer idatomp(4)
      integer idatomm(4)
      integer atoms2move(4)
      integer num_atoms2move, atom1, atom2
      integer nat3
      integer ixyz, zatom
      double precision thresh, delta, scale, normmax, normmin, norm
      character*99 string,strings, stringe
*
      normmax = 0.0d00
      normmin = 1.2345678d05
      call int_acc_high()
* store original coordintates
      nat3 = 3*nat
      call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
      thresh = 1.0d-07
      delta  = 0.001d00
*
      write(luout,*)'  ',nshell,' total shells '
      call util_flush(luout)
*      do ish = 1, nshell
      do ish = 1, 1
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &      stop 'cn2bfr error i'
        inbf = ihi - ilo + 1
        strings = ' '
        call util_date(strings)
        if (.not.bas_cn2ce(basis,ish,iatom))
     &      stop 'bas_cn2ce error i'
*        do jsh = 1, nshell
        do jsh = 1, 1
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &          stop 'cn2bfr error j'
          jnbf = jhi - jlo + 1
          if (.not.bas_cn2ce(basis,jsh,jatom))
     &          stop 'bas_cn2ce error j'
*          do ksh =  1, nshell
          do ksh =  1, 6
            if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &            stop 'cn2bfr error k'
            knbf = khi - klo + 1
            if (.not.bas_cn2ce(basis,ksh,katom))
     &            stop 'bas_cn2ce error k'
            if (.not.(iatom.eq.jatom.and.
     &            jatom.eq.katom)) then
              nint  = inbf*jnbf*knbf
              nintg = nint*12
              ninth = nint*78
              call dfill (lbuf,0.0d00,buf,1)
              call dcopy(nat3,xyz,1,coords(1,1,geom),1)
              call intdd_2e3c(basis,ish,basis,jsh,ksh,
     &              lscr,scr,lbuf,buf,idatom)
              call rakdd_checkt_78(nint,buf)
              write(luout,*)idatom
              idatoms(1) = idatom(1)
              idatoms(2) = idatom(2)
              idatoms(3) = idatom(3)
              idatoms(4) = idatom(4)
              write(6,*)' idatom  1:', idatom 
              write(6,*)' idatoms 1:', idatoms 
              if (.not.ma_alloc_get(mt_dbl,lbuf,
     &              'copy of buf',hbufcp,kbufcp)) stop ' ma alloc fail'
              call dcopy(lbuf,buf,1,dbl_mb(kbufcp),1)
              call rakdd_fill(12,nint,dbl_mb(kbufcp),bufsum,idatoms)
              write(6,*)' idatom  2:', idatom 
              write(6,*)' idatoms 2:', idatoms 
              if (.not.ma_free_heap(hbufcp)) stop ' ma_free fail'
              call rakdd_print_dd(nint,buf,idatom)
              if (.not.rakdd_checktrans(nint,bufsum)) then
                call rakdd_print_dd(nint,buf,idatom)
                call rakdd_printsum(nint,bufsum,idatoms)
              endif
              nzero = 0
              ncount = 0
              do count = 1,ninth
                if (abs(buf(count)).gt.thresh) then
                  ncount = ncount + 1
                else
                  nzero = nzero + 1
                endif
              enddo
              write(luout,*)nzero,'+',ncount,' != ',ninth
              if ((nzero+ncount).ne.ninth)
     &              write(luout,*)nzero,'+',ncount,' != ',ninth
              atoms2move(1) = iatom
              atoms2move(2) = iatom
              atoms2move(3) = jatom
              atoms2move(4) = katom
              num_atoms2move = 0
              do atom1 = 1,4
                do atom2 = 2,4
                  if (atom1.ne.atom2) then
                    if (atoms2move(atom1).eq.atoms2move(atom2))
     &                    atoms2move(atom2) = 0
                  endif
                enddo
              enddo
              num_atoms2move = 0
              do atom1 = 1,4
                if (atoms2move(atom1).gt.0)
     &                num_atoms2move = num_atoms2move + 1
              enddo
*
              call dfill (lbuffdsq,0.0d00,buffdsq,1)
              write(6,*)' atoms2move ',atoms2move
              do zatom = 1,4
                if (atoms2move(zatom).gt.0) then
                  do ixyz = 1,3
                    nintg = nint*12
                    call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                    coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) + delta
                    call dfill(lgrad,0.0d00,gradp,1)
                    call intd_2e3c(basis,ish,basis,jsh,ksh,
     &                    lscr,scr,lgrad,gradp,idatomp)
                    write(string,*)' grad +',ixyz,zatom
*                      call rakdd_printgrad(nint,gradp,string)
                    write(6,*)'idatomp:',idatomp
                    call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                    coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) - delta
                    call dfill(lgrad,0.0d00,gradm,1)
                    call intd_2e3c(basis,ish,basis,jsh,ksh,
     &                    lscr,scr,lgrad,gradm,idatomm)
                    write(string,*)' grad -',ixyz,zatom
*                      call rakdd_printgrad(nint,gradm,string)
                    call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                    write(6,*)'idatomm:',idatomm
                    call daxpy(nintg,-1.0d00,gradm,1,gradp,1)
                    write(string,*)' grad diff',ixyz,zatom
*                      call rakdd_printgrad(nint,gradp,string)
                    scale = 1.0d00/(2.0d00*delta)
                    call dscal(nintg,scale,gradp,1)
                    write(string,*)' grad diff scaled',ixyz,zatom
*                      call rakdd_printgrad(nint,gradp,string)
                    call rakdd_grad_fill(gradp,idatomp,
     &                    idatomm,ixyz,
     &                    zatom,idatoms,
     &                    atoms2move(zatom),buffdsq,nint)
                  enddo
                endif
              enddo             ! zatom
              call rakdd_fill_b(buffdsq,nint)  ! check if symmetric
              call rakdd_printsq(buffdsq,nint)
*                call rakdd_print_both(80,buffdsq,bufsum,nint)
              call rakdd_compare3c(
     &              ish,jsh,ksh,
     &              (nint*12*12),buffdsq,bufsum,norm)
              normmax = max(normmax,norm)
              normmin = min(normmin,norm)
            endif               ! 4 atoms the same
          enddo                 ! ksh
        enddo                   ! jsh
        string  = ' '
        stringe = ' '
        call util_date(stringe)
        write(string,'(1x,i5,a3,a27,a2,a27)')
     &        ish,' ::',
     &        strings(1:25),'::',
     &        stringe(1:25)
        do jsh = 1,inp_strlen(string)
          ksh = ichar(string(jsh:jsh))
          if (ksh.eq.10) string(jsh:jsh) = ' '
        enddo
        write(luout,'(1x,a)')string(1:inp_strlen(string))
      enddo                     ! ish
      write(luout,'(1x,a,1pd20.10)')
     &      ' maximum difference norm over all quartets:',normmax
      write(luout,'(1x,a,1pd20.10)')
     &      ' minimum difference norm over all quartets:',normmin
*
10000 format(1x,a,4(i4),1x,a,4(i3),1x,a,4(1x,i5),1x,a,i6,/,
     &      47x,a,4(1x,i5),1x,a,i2)
      end
      subroutine rakdd_compare3c(
     &              ish,jsh,ksh,
     &              len,buffd,buf,norm)
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
      integer ish,jsh,ksh
      integer len
      double precision buffd(len), buf(len)
*
      double precision ddot
      external ddot
*
      double precision norm, thresh
      integer h_diff, k_diff
*
      thresh = 1.0d-07
      if (.not.ma_alloc_get(mt_dbl,len,
     &    'diff buf',h_diff,k_diff)) stop ' ma alloc fail'

      call dcopy(len,buf,1,dbl_mb(k_diff),1)
      call daxpy(len,-1.0d00,buffd,1,dbl_mb(k_diff),1)
      norm = ddot(len,dbl_mb(k_diff),1,dbl_mb(k_diff),1)
*      if (norm.gt.thresh)
*     &    
      write(luout,*)' ',ish,jsh,ksh,
     &    ' norm = ',norm
      if (.not.ma_free_heap(h_diff)) stop ' ma free fail'
      end
      subroutine rakdd_printsq(buffsq,nint)
      implicit none
      integer nint
      double precision buffsq (3,4,3,4,nint)
c
      integer ia, ja, ix, jx, int
      logical doit
c
      do int = 1,nint
        do ia = 1,4
          do ja = 1,4
            if (ia.le.ja) then
              do ix = 1,3
                do jx = 1,3
                  doit = .true.
                  if (ia.eq.ja) doit = ix.le.jx
                  if (doit) then
                    write(6,10000)ix,ia,jx,ja,int,
     &                    buffsq(ix,ia,jx,ja,int)
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
10000 format(1x,'dd(',4(i4,','),i4,') =',1pd20.10)
      end
