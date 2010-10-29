*
* $Id$
*
      logical function raktask_intdd(rtdb)
* test the second derivative integrals api
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
      raktask_intdd = .false.
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
      call raktask_intdd_a(geom,basis,nbf,nshell,cn_nbf_max,
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
      raktask_intdd = bas_destroy(basis)
      raktask_intdd = raktask_intdd.and.
     &    geom_destroy(geom)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hscr)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hbuf)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hfd)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hfdsq)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hxyz)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hgradp)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hgradm)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hbufsum)
      end
      subroutine raktask_intdd_a(geom,basis,nbf,nshell,cn_nbf_max,
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
      integer ish, jsh, ksh, lsh
      integer ilo, ihi, inbf
      integer jlo, jhi, jnbf
      integer klo, khi, knbf
      integer llo, lhi, lnbf
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
      normmin = 1.0d05
      call int_acc_high()
* store original coordintates
      nat3 = 3*nat
      call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
      thresh = 1.0d-07
      delta  = 0.001d00
*
      write(luout,*)'  ',nshell,' total shells '
      do ish = 1, nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &      stop 'cn2bfr error i'
        inbf = ihi - ilo + 1
        strings = ' '
        call util_date(strings)
        if (.not.bas_cn2ce(basis,ish,iatom))
     &      stop 'bas_cn2ce error i'
        do jsh = 1, nshell
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &        stop 'cn2bfr error j'
          jnbf = jhi - jlo + 1
          if (.not.bas_cn2ce(basis,jsh,jatom))
     &        stop 'bas_cn2ce error j'
          do ksh =  1, nshell
            if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &          stop 'cn2bfr error k'
            knbf = khi - klo + 1
            if (.not.bas_cn2ce(basis,ksh,katom))
     &          stop 'bas_cn2ce error k'
            do lsh =  1, nshell
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi))
     &            stop 'cn2bfr error l'
              lnbf = lhi - llo + 1
              if (.not.bas_cn2ce(basis,lsh,latom))
     &            stop 'bas_cn2ce error l'
              if (.not.(iatom.eq.jatom.and.
     &                           jatom.eq.katom.and.
     &                                    katom.eq.latom)) then
*                write(79,'(/)')
*                write(79,10000)
*     &              'shells:s: ',ish,jsh,ksh,lsh
                nint  = inbf*jnbf*knbf*lnbf
                nintg = nint*12
                ninth = nint*78
*rak:                write(luout,'(/)')
*rak:                write(luout,10000)
*rak:     &              'shells:s: ',ish,jsh,ksh,lsh,
*rak:     &              'nbfX: ',inbf,jnbf,knbf,lnbf,
*rak:     &              'nint:nintg:ninth:',nint,nintg,ninth
                call dfill (lbuf,0.0d00,buf,1)
                call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                call intdd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &              lscr,scr,lbuf,buf,idatom)
                call rakdd_checkt_78(nint,buf)
                idatoms(1) = idatom(1)
                idatoms(2) = idatom(2)
                idatoms(3) = idatom(3)
                idatoms(4) = idatom(4)
                if (.not.ma_alloc_get(mt_dbl,lbuf,
     &              'copy of buf',hbufcp,kbufcp)) stop ' ma alloc fail'
                call dcopy(lbuf,buf,1,dbl_mb(kbufcp),1)
                call rakdd_fill(12,nint,dbl_mb(kbufcp),bufsum,idatoms)
                if (.not.ma_free_heap(hbufcp)) stop ' ma_free fail'
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
                if ((nzero+ncount).ne.ninth)
     &              write(luout,*)nzero,'+',ncount,' != ',ninth
                atoms2move(1) = iatom
                atoms2move(2) = jatom
                atoms2move(3) = katom
                atoms2move(4) = latom
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
*rak:00001           continue
*rak:                do atom1 = 1,4
*rak:                  if (atoms2move(atom1).eq.0) then
*rak:                    do atom2 = atom1,3
*rak:                      atoms2move(atom2) = atoms2move(atom2+1)
*rak:                    enddo
*rak:                    atoms2move(4)=0
*rak:                  endif
*rak:                enddo
*rak:                do atom1 = 1,num_atoms2move
*rak:                  if (atoms2move(atom1).eq.0) goto 00001
*rak:                enddo
*
*rak:                write(luout,10000)
*rak:     &              'shells : :',ish,jsh,ksh,lsh,
*rak:     &              'idatom : ',idatom,
*rak:     &              'idatoms: ',idatoms,
*rak:     &              'ncount : ',ncount,
*rak:     &              'atoms2move: ',atoms2move,'num',num_atoms2move
*
                call dfill (lbuffdsq,0.0d00,buffdsq,1)
                do zatom = 1,4
                  if (atoms2move(zatom).gt.0) then
                    do ixyz = 1,3
                      nintg = nint*12
                      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                      coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) + delta
                      call dfill(lgrad,0.0d00,gradp,1)
                      call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &                    lscr,scr,lgrad,gradp,idatomp)
                      write(string,*)' grad +',ixyz,zatom
*                      call rakdd_printgrad(nint,gradp,string)
*                      write(6,*)'idatomp:',idatomp
                      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                      coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) - delta
                      call dfill(lgrad,0.0d00,gradm,1)
                      call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &                    lscr,scr,lgrad,gradm,idatomm)
                      write(string,*)' grad -',ixyz,zatom
*                      call rakdd_printgrad(nint,gradm,string)
                      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
*                      write(6,*)'idatomm:',idatomm
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
                enddo ! zatom
                call rakdd_fill_b(buffdsq,nint)
*                call rakdd_print_both(80,buffdsq,bufsum,nint)
                call rakdd_compare(
     &              ish,jsh,ksh,lsh,
     &              (nint*12*12),buffdsq,bufsum,norm)
                normmax = max(normmax,norm)
                normmin = min(normmin,norm)
              endif ! 4 atoms the same
            enddo   ! lsh
          enddo     ! ksh
        enddo       ! jsh
        string  = ' '
        stringe = ' '
        call util_date(stringe)
        write(string,'(1x,i5,a3,a27,a2,a27)')
     &      ish,' ::',
     &      strings(1:25),'::',
     &      stringe(1:25)
        do jsh = 1,inp_strlen(string)
          ksh = ichar(string(jsh:jsh))
          if (ksh.eq.10) string(jsh:jsh) = ' '
        enddo
        write(luout,'(1x,a)')string(1:inp_strlen(string))
      enddo         ! ish
      write(luout,'(1x,a,1pd20.10)')
     &    ' maximum difference norm over all quartets:',normmax
      write(luout,'(1x,a,1pd20.10)')
     &    ' minimum difference norm over all quartets:',normmin
*
10000 format(1x,a,4(i4),1x,a,4(i3),1x,a,4(1x,i5),1x,a,i6,/,
     &    47x,a,4(1x,i5),1x,a,i2)
      end
      subroutine rakdd_fill(at4xyz,nint,buf,bufsum,idatom)
      implicit none
#include "mafdecls.fh"
      integer at4xyz, nint
      integer idatom(4)
      double precision buf(78,nint)
      double precision bufsum(3,4,3,4,nint)
c
      double precision factor, val
      integer itri
      integer ia,ja,ixyz,jxyz,iint
      logical doit
      integer hbufsumcp, kbufsumcp
      integer idatomB4(4)
c
      idatomB4(1) = idatom(1)
      idatomB4(2) = idatom(2)
      idatomB4(3) = idatom(3)
      idatomB4(4) = idatom(4)
      itri = 0
      do ia = 1,4
        do ja = 1,4
          if (ia.le.ja) then
            do ixyz = 1,3
              do jxyz = 1,3
                doit = .true.
                if (ia.eq.ja) doit = ixyz.le.jxyz
                if (doit) then
                  factor = 1.0d00
                  itri = itri + 1
                  do iint = 1,nint
                    val = factor*buf(itri,iint)
                    bufsum(ixyz,ia,jxyz,ja,iint) = val
                    bufsum(jxyz,ja,ixyz,ia,iint) = val
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo
c
      if (.not.ma_alloc_get(mt_dbl,(12*12*nint),
     &    'copy of bufsum',hbufsumcp,kbufsumcp)) stop ' ma alloc fail'
      call dcopy((12*12*nint),bufsum,1,dbl_mb(kbufsumcp),1)
      call rakdd_fill_a(bufsum,idatom,nint) ! sum atom parts together
*      call rakdd_printsum_b4a(79,bufsum,dbl_mb(kbufsumcp),nint,
*     &    idatom,idatomB4)
      if (.not.ma_free_heap(hbufsumcp)) stop ' ma_free fail'
      call rakdd_fill_b(bufsum,nint) ! check if symmetric
      end
      subroutine rakdd_fill_b(bufsum,nint)
      implicit none
#include "stdio.fh"
      integer nint
      double precision bufsum(12,12,nint)
*
      double precision diff, thresh
      integer i,j,iint
*
      thresh = 1.0d-12
      do iint = 1,nint
        do i=1,12
          do j = 1,i
            diff = bufsum(i,j,iint)-bufsum(j,i,iint)
            diff = abs(diff)
            if (diff.gt.thresh) then
              write(luout,*)' bufsum not symmetric i,j,iint',i,j,iint
              write(luout,*)' bufsum(i,j,iint) :',bufsum(i,j,iint)
              write(luout,*)' bufsum(j,i,iint) :',bufsum(j,i,iint)
              write(luout,*)' diff            :',diff
            endif
          enddo
        enddo
      enddo
      end
      subroutine rakdd_fill_a(buf,idatom,nint)
      implicit none
#include "stdio.fh"
      integer idatom(4)
      integer nint
      double precision buf(3,4,3,4,nint)
*
      logical sum12, sum13, sum14, sum23, sum24, sum34
      integer ixyz, jxyz, iatom, jatom, iint
      integer num_sum
*
      num_sum = 0
      if (idatom(1).gt.0.and.idatom(2).gt.0) then
        sum12 = idatom(1).eq.idatom(2)
        if (sum12) then
          num_sum = num_sum + 1
          idatom(2) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,1,jxyz,jatom,iint) =
     &                buf(ixyz,1,jxyz,jatom,iint) +
     &                buf(ixyz,2,jxyz,jatom,iint) 
                  buf(ixyz,2,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,1,iint) =
     &                buf(ixyz,iatom,jxyz,1,iint) +
     &                buf(ixyz,iatom,jxyz,2,iint) 
                  buf(ixyz,iatom,jxyz,2,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(1).gt.0.and.idatom(3).gt.0) then
        sum13 = idatom(1).eq.idatom(3)
        if (sum13) then
          num_sum = num_sum + 1
          idatom(3) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,1,jxyz,jatom,iint) =
     &                buf(ixyz,1,jxyz,jatom,iint) +
     &                buf(ixyz,3,jxyz,jatom,iint) 
                  buf(ixyz,3,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,1,iint) =
     &                buf(ixyz,iatom,jxyz,1,iint) +
     &                buf(ixyz,iatom,jxyz,3,iint) 
                  buf(ixyz,iatom,jxyz,3,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(1).gt.0.and.idatom(4).gt.0) then
        sum14 = idatom(1).eq.idatom(4)
        if (sum14) then
          num_sum = num_sum + 1
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,1,jxyz,jatom,iint) =
     &                buf(ixyz,1,jxyz,jatom,iint) +
     &                buf(ixyz,4,jxyz,jatom,iint) 
                  buf(ixyz,4,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,1,iint) =
     &                buf(ixyz,iatom,jxyz,1,iint) +
     &                buf(ixyz,iatom,jxyz,4,iint) 
                  buf(ixyz,iatom,jxyz,4,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(2).gt.0.and.idatom(3).gt.0) then
        sum23 = idatom(2).eq.idatom(3)
        if (sum23) then
          num_sum = num_sum + 1
          idatom(3) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,2,jxyz,jatom,iint) =
     &                buf(ixyz,2,jxyz,jatom,iint) +
     &                buf(ixyz,3,jxyz,jatom,iint) 
                  buf(ixyz,3,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,2,iint) =
     &                buf(ixyz,iatom,jxyz,2,iint) +
     &                buf(ixyz,iatom,jxyz,3,iint) 
                  buf(ixyz,iatom,jxyz,3,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(2).gt.0.and.idatom(4).gt.0) then
        sum24 = idatom(2).eq.idatom(4)
        if (sum24) then
          num_sum = num_sum + 1
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,2,jxyz,jatom,iint) =
     &                buf(ixyz,2,jxyz,jatom,iint) +
     &                buf(ixyz,4,jxyz,jatom,iint) 
                  buf(ixyz,4,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,2,iint) =
     &                buf(ixyz,iatom,jxyz,2,iint) +
     &                buf(ixyz,iatom,jxyz,4,iint) 
                  buf(ixyz,iatom,jxyz,4,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(3).gt.0.and.idatom(4).gt.0) then
        sum34 = idatom(3).eq.idatom(4)
        if (sum34) then
          num_sum = num_sum + 1
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do iint = 1,nint
                  buf(ixyz,3,jxyz,jatom,iint) =
     &                buf(ixyz,3,jxyz,jatom,iint) +
     &                buf(ixyz,4,jxyz,jatom,iint) 
                  buf(ixyz,4,jxyz,jatom,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do iint = 1,nint
                  buf(ixyz,iatom,jxyz,3,iint) =
     &                buf(ixyz,iatom,jxyz,3,iint) +
     &                buf(ixyz,iatom,jxyz,4,iint) 
                  buf(ixyz,iatom,jxyz,4,iint) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (num_sum.gt.2)
     &    write(luout,*)'number of sums done: ',num_sum
      end
      subroutine rakdd_printsum(nint,bufsum,idatom)
      implicit none
#include "stdio.fh"
      integer nint
      integer idatom(4)
      double precision bufsum(3,4,3,4,nint)
c
      integer ixyz,iatom,jxyz,jatom,iint
      double precision thresh, val
      integer count,itri
      character*1 atoms(4),sxyz(3)
c
      atoms(1)='a'
      atoms(2)='b'
      atoms(3)='c'
      atoms(4)='d'
      sxyz(1)='x'
      sxyz(2)='y'
      sxyz(3)='z'
      thresh = 1.0d-07
      itri = 0
      count = 0
      do ixyz = 1,3
        do iatom = 1,4
          do jxyz = 1,3
            do jatom = 1,4
              if (iatom.eq.jatom) then
                if (ixyz.le.jxyz) itri = itri + 1
              else if (iatom.le.jatom) then
                itri = itri + 1
              endif
              do iint = 1,nint
                val = bufsum(ixyz,iatom,jxyz,jatom,iint)
                if (abs(val).gt.thresh) then
                  count = count + 1
                  write(luout,10000)
     &                atoms(iatom),sxyz(ixyz),
     &                atoms(jatom),sxyz(jxyz),
     &                ixyz,iatom,jxyz,jatom,iint,
     &                itri,val,
     &                idatom(iatom),idatom(jatom), count
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'bufsum(',4(a1),',',
     &    i4,4(',',i4),')(',i5,') =',
     &    1pd20.10,2x,2(i5),1x,i8)
      end
      subroutine rakdd_print_dd(nint,buf,idatom)
      implicit none
#include "stdio.fh"
      integer nint
      integer idatom(4)
      double precision buf(78,nint)
c
      double precision thresh, val
      integer ider,iint, count
      character*4 sder(78)
c      
      sder( 1) = 'axax'
      sder( 2) = 'axay'
      sder( 3) = 'axaz'
      sder( 4) = 'ayay'
      sder( 5) = 'ayaz'
      sder( 6) = 'azaz'
      sder( 7) = 'axbx'
      sder( 8) = 'axby'
      sder( 9) = 'axbz'
      sder(10) = 'aybx'
      sder(11) = 'ayby'
      sder(12) = 'aybz'
      sder(13) = 'azbx'
      sder(14) = 'azby'
      sder(15) = 'azbz'
      sder(16) = 'axcx'
      sder(17) = 'axcy'
      sder(18) = 'axcz'
      sder(19) = 'aycx'
      sder(20) = 'aycy'
      sder(21) = 'aycz'
      sder(22) = 'azcx'
      sder(23) = 'azcy'
      sder(24) = 'azcz'
      sder(25) = 'axdx'
      sder(26) = 'axdy'
      sder(27) = 'axdz'
      sder(28) = 'aydx'
      sder(29) = 'aydy'
      sder(30) = 'aydz'
      sder(31) = 'azdx'
      sder(32) = 'azdy'
      sder(33) = 'azdz'
      sder(34) = 'bxbx'
      sder(35) = 'bxby'
      sder(36) = 'bxbz'
      sder(37) = 'byby'
      sder(38) = 'bybz'
      sder(39) = 'bzbz'
      sder(40) = 'bxcx'
      sder(41) = 'bxcy'
      sder(42) = 'bxcz'
      sder(43) = 'bycx'
      sder(44) = 'bycy'
      sder(45) = 'bycz'
      sder(46) = 'bzcx'
      sder(47) = 'bzcy'
      sder(48) = 'bzcz'
      sder(49) = 'bxdx'
      sder(50) = 'bxdy'
      sder(51) = 'bxdz'
      sder(52) = 'bydx'
      sder(53) = 'bydy'
      sder(54) = 'bydz'
      sder(55) = 'bzdx'
      sder(56) = 'bzdy'
      sder(57) = 'bzdz'
      sder(58) = 'cxcx'
      sder(59) = 'cxcy'
      sder(60) = 'cxcz'
      sder(61) = 'cycy'
      sder(62) = 'cycz'
      sder(63) = 'czcz'
      sder(64) = 'cxdx'
      sder(65) = 'cxdy'
      sder(66) = 'cxdz'
      sder(67) = 'cydx'
      sder(68) = 'cydy'
      sder(69) = 'cydz'
      sder(70) = 'czdx'
      sder(71) = 'czdy'
      sder(72) = 'czdz'
      sder(73) = 'dxdx'
      sder(74) = 'dxdy'
      sder(75) = 'dxdz'
      sder(76) = 'dydy'
      sder(77) = 'dydz'
      sder(78) = 'dzdz'
      thresh = 1.0d-07
      count = 0
      do ider = 1,78
        do iint = 1,nint
          val = buf(ider,iint)
          if (abs(val).gt.thresh) then
            count = count + 1
            write(luout,10000)sder(ider),
     &          ider,iint,val,idatom,count
          endif
        enddo
      enddo
10000 format(1x,'buf(',a4,',',i4,',',i4,') =',1pd20.10,
     &    2x,4(i5),1x,i8)
      end
      logical function rakdd_checktrans(nint,bufsum)
* true is okay
* false is a problem
      implicit none
#include "stdio.fh"
      integer nint
      double precision bufsum(3,4,3,4,nint)
*
      double precision sum, summax, summin
      integer ia,ja,ic,jc,iint
      double precision thresh
      logical first_pass
*
*----------------------------------------------------
* trans invariance of second derivative int Ikaja
* sum(k:k=1,N)[Ikajb] = 0 a,b = xyz;j=1,N
*----------------------------------------------------
* sum(k:k=1,N;a:a=xyz)
* AxAx + BxAx + CxAx + DxAx = 0
* AxBx + BxBx + CxBx + DxBx = 0
* AxCx + BxCx + CxCx + DxCx = 0
* AxDx + BxDx + CxDx + DxDx = 0
*
      rakdd_checktrans = .true.
      thresh = 1.0d-07
      first_pass = .true.
00001 continue
      summax = -565.000d00
      summin =  565.000d00
      do iint = 1,nint
        sum = 0.0d00
        do ja = 1,4
          do jc = 1,3
            do ic = 1,3
              do ia = 1,4
                sum = sum + bufsum(ic,ia,jc,ja,iint)
              enddo
            enddo
          enddo
        enddo
        if (.not.(abs(sum).lt.thresh)) then
          write(luout,*)
     &        ' sum of 144 NOT OKAY iint:',iint,' sum:',sum
          if (rakdd_checktrans)rakdd_checktrans=.false.
        endif
      enddo
      end
      subroutine rakdd_checkt_78(nint,buf)
      implicit none
#include "stdio.fh"
*      
      integer nint
      double precision buf(78,nint)
*
      double precision sum, f, thresh, sum78
      double precision sumjim(6)
      integer ider, iint, i
*
      thresh = 1.0d-07
      do iint = 1,nint
        sum78 = 0.0d00
        do ider = 1,78
          sum78 = sum78 + buf(ider,iint)
        enddo
        sumjim(1) =
     &        buf( 2,iint)  + ! IAxAy   
     &        buf( 8,iint)  + ! IAxBy   
     &        buf(17,iint)  + ! IAxCy
     &        buf(26,iint)  + ! IAxDy
     &        buf(10,iint)  + ! IBxAy aybx
     &        buf(35,iint)  + ! IBxBy
     &        buf(41,iint)  + ! IBxCy
     &        buf(50,iint)  + ! IBxDy
     &        buf(19,iint)  + ! ICxAy aycx
     &        buf(43,iint)  + ! ICxBy bycx
     &        buf(59,iint)  + ! ICxCy
     &        buf(65,iint)  + ! ICxDy
     &        buf(28,iint)  + ! IDxAy aydx
     &        buf(52,iint)  + ! IDxBy bydx
     &        buf(67,iint)  + ! IDxCy cydx
     &        buf(74,iint)    ! IDxDy
*        write(luout,*)' sum of 78(1) :',sumjim(1)
        sumjim(2) =
     &        buf( 5,iint)  + ! IAyAz
     &        buf(12,iint)  + ! IAyBz
     &        buf(21,iint)  + ! IAyCz
     &        buf(30,iint)  + ! IAyDz
     &        buf(14,iint)  + ! IByAz azby
     &        buf(38,iint)  + ! IByBz
     &        buf(45,iint)  + ! IByCz
     &        buf(54,iint)  + ! IByDz
     &        buf(23,iint)  + ! ICyAz azcy
     &        buf(47,iint)  + ! ICyBz bzcy
     &        buf(62,iint)  + ! ICyCz 
     &        buf(69,iint)  + ! ICyDz
     &        buf(32,iint)  + ! IDyAz azdy
     &        buf(56,iint)  + ! IDyBz bzdy
     &        buf(71,iint)  + ! IDyCz czdy
     &        buf(77,iint)    ! IDyDz
*        write(luout,*)' sum of 78(2) :',sumjim(2)
        sumjim(3) =
     &        buf( 3,iint)  + ! IAzAx
     &        buf(13,iint)  + ! IAzBx
     &        buf(22,iint)  + ! IAzCx
     &        buf(31,iint)  + ! IAzDx
     &        buf( 9,iint)  + ! IBzAx axbz
     &        buf(36,iint)  + ! IBzBx bxbz
     &        buf(46,iint)  + ! IBzCx 
     &        buf(55,iint)  + ! IBzDx
     &        buf(18,iint)  + ! ICzAx axcz
     &        buf(42,iint)  + ! ICzBx bxcz
     &        buf(60,iint)  + ! ICzCx cxcz
     &        buf(70,iint)  + ! ICzDx 
     &        buf(27,iint)  + ! IDzAx axdz
     &        buf(51,iint)  + ! IDzBx bxdz
     &        buf(66,iint)  + ! IDzCx cxdz
     &        buf(75,iint)    ! IDzDx dxdz
*        write(luout,*)' sum of 78(3) :',sumjim(3)
        f = 2.0d00
        sumjim(4) =
     &        buf( 1,iint)  + ! IAxAx      
     &      f*buf( 7,iint)  + !  2 * IBxAx  axbx
     &        buf(34,iint)  + ! IBxBx     
     &      f*buf(16,iint)  + !  2 * ICxAx  axcx
     &      f*buf(40,iint)  + !  2 * ICxBx  bxcx
     &        buf(58,iint)  + ! ICxCx     
     &      f*buf(25,iint)  + !  2 * IDxAx  axdx
     &      f*buf(49,iint)  + !  2 * IDxBx  bxdx
     &      f*buf(64,iint)  + !  2 * IDxCx  cxdx
     &        buf(73,iint)    ! IDxDx     
*        write(luout,*)' sum of 78(4) :',sumjim(4)
        f = 2.0d00
        sumjim(5) =
     &        buf( 4,iint)  + ! IAyAy     
     &      f*buf(11,iint)  + !  2 * IByAy ayby
     &        buf(37,iint)  + ! IByBy     
     &      f*buf(20,iint)  + !  2 * ICyAy aycy
     &      f*buf(44,iint)  + !  2 * ICyBy bycy
     &        buf(61,iint)  + ! ICyCy     
     &      f*buf(29,iint)  + !  2 * IDyAy aydy
     &      f*buf(53,iint)  + !  2 * IDyBy bydy
     &      f*buf(68,iint)  + !  2 * IDyCy cydy
     &        buf(76,iint)    ! IDyDy     
*        write(luout,*)' sum of 78(5) :',sumjim(5)
        f = 2.0d00
        sumjim(6) =
     &        buf( 6,iint)  + ! IAzAz      
     &      f*buf(15,iint)  + !  2 * IBzAz azbz
     &        buf(39,iint)  + ! IBzBz      
     &      f*buf(24,iint)  + !  2 * ICzAz azcz
     &      f*buf(48,iint)  + !  2 * ICzBz bzcz
     &        buf(63,iint)  + ! ICzCz      
     &      f*buf(33,iint)  + !  2 * IDzAz azdz
     &      f*buf(57,iint)  + !  2 * IDzBz bzdz
     &      f*buf(72,iint)  + !  2 * IDzCz czdz
     &        buf(78,iint)    ! IDzDz      
*        write(luout,*)' sum of 78(6) :',sumjim(6)
        sum = 0.0d00
        do i = 1,6
          sum = sum + sumjim(i)
        enddo
        if (.not.(abs(sum).lt.thresh)) then
          write(luout,'(1x,a,i7,a,2x,6(1pd14.6))')
     &        ' sum of 78(1-6):',iint,':',(sumjim(i),i=1,6)
          write(luout,*)' sum of 78 :',iint,':',sum78
          write(luout,'(1x,a,i7,a,2x,1pd14.6)')
     &        ' sum of 78(*) NOT OKAY  :',iint,':',sum
        endif
      enddo
      end
      subroutine rakdd_printsum_b4a(lu,bufsumAfter,bufsumB4,nint,
     &    idA,idB4)
      implicit none
      integer lu ! unit to write to
      integer nint
      integer idA(4), idB4(4)
      double precision bufsumAfter(3,4,3,4,nint)
      double precision bufsumB4   (3,4,3,4,nint)
c
      integer iatom, ixyz, jatom, jxyz, iint
      integer i, j, num_sum
      double precision before, after, thresh
      logical sum12, sum13, sum14, sum23, sum24, sum34
c
      j = 0
      do i=1,4
        if (idb4(i).gt.0) j = j + 1
      enddo
      if (j.ne.4) write(lu,*)' j not 4 ', j
      sum12 = idb4(1).eq.idb4(2)
      sum13 = idb4(1).eq.idb4(3)
      sum14 = idb4(1).eq.idb4(4)
      sum23 = idb4(2).eq.idb4(3)
      sum24 = idb4(2).eq.idb4(4)
      sum34 = idb4(3).eq.idb4(4)
      if (sum12.and.sum13)sum23=.false.
      if (sum12.and.sum14)sum24=.false.
      if (sum13.and.sum14)sum34=.false.
      if (sum23.and.sum24)sum34=.false.
c
      num_sum = 0
      if (sum12)num_sum = num_sum + 1
      if (sum13)num_sum = num_sum + 1
      if (sum14)num_sum = num_sum + 1
      if (sum23)num_sum = num_sum + 1
      if (sum24)num_sum = num_sum + 1
      if (sum34)num_sum = num_sum + 1
      write(lu,*)'num_sum=',num_sum,'  ',
     &    '(12',sum12,') ',
     &    '(13',sum13,') ',
     &    '(14',sum14,') ',
     &    '(23',sum23,') ',
     &    '(24',sum24,') ',
     &    '(34',sum34,') '
c
      thresh = 1.0d-07
      do ixyz = 1,3
        do jxyz = 1,3
          do iint = 1,nint
            write(lu,*)'idb4',idb4,' ida',ida
            do iatom = 1,4
              do jatom = 1,4
                before = bufsumB4(ixyz,iatom,jxyz,jatom,iint)
                after  = bufsumAfter(ixyz,iatom,jxyz,jatom,iint)
                write(lu,10000)ixyz,iatom,jxyz,jatom,iint,before,after
              enddo
            enddo
          enddo
        enddo
      enddo
10000 format(1x,'<',i1,',',i2,',',i1,',',i2,',',i5,'>=',2(1pd12.4))
      end
      subroutine rakdd_grad_fill(grad,idatomp,idatomm,
     &    ixyz,zatom,
     &    idatoms,matom,buffd,nint)
      implicit none
#include "stdio.fh"
      integer idatomp(4)
      integer idatomm(4)
      integer idatoms(4)
      integer ixyz, zatom, matom, nint
      double precision grad(nint,3,4)
      double precision buffd(3,4,3,4,nint)
c
      double precision val
      integer i, iint, iga, igc
      logical okay
c
      okay = .true.
      do i = 1,4
        if (idatomp(i).gt.0.and.
     &      idatomm(i).gt.0.and.
     &      idatoms(i).gt.0) then
          okay = okay.and.idatomp(i).eq.idatomm(i)
          okay = okay.and.idatomp(i).eq.idatoms(i)
        endif
      enddo
      if (.not.okay) then
        write(luout,*)' something not okay '
        write(luout,*)' idatomp ',idatomp
        write(luout,*)' idatomm ',idatomm
        write(luout,*)' idatoms ',idatoms
      endif
*rak:      write(luout,*)'grad_fill:zatom:',zatom
*rak:      write(luout,*)'grad_fill:matom:',matom
*rak:      write(luout,*)'grad_fill:ixyz :',ixyz
c
*      write(luout,*)' grad',grad
      do iga = 1,4
        if (idatomp(iga).gt.0) then
          do igc = 1,3
            do iint = 1,nint
              val = grad(iint,igc,iga)
              buffd(ixyz,zatom,igc,iga,iint) = val
              buffd(igc,iga,ixyz,zatom,iint) = val
            enddo
          enddo
        endif
      enddo
c
      end
      subroutine rakdd_compare(
     &    ish,jsh,ksh,lsh,
     &    len,buffd,buf,norm)
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
      integer ish,jsh,ksh,lsh
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
      if (norm.gt.thresh)
     &    write(luout,*)' ',ish,jsh,ksh,lsh,
     &    ' norm = ',norm
      if (.not.ma_free_heap(h_diff)) stop ' ma free fail'
      end
      subroutine rakdd_print_both(lu,buffd,buf,nint)
      implicit none
      integer lu
      integer nint
      double precision buffd(3,4,3,4,nint)
      double precision buf  (3,4,3,4,nint)
c
      integer ia, ja, ca, cj, iint
c
      write(lu,*)'                     finit diff  analytic'
      do ca = 1,3
        do ia = 1,4
          do cj = 1,3
            do ja = 1,4
              do iint = 1,nint
                write(lu,10000)ca,ia,cj,ja,iint,
     &              buffd(ca,ia,cj,ja,iint),
     &              buf  (ca,ia,cj,ja,iint)
              enddo
            enddo
          enddo
        enddo
      enddo
      write(lu,*)' '
      write(lu,*)' '
      write(lu,*)' '
10000 format(1x,'<',i1,',',i2,',',i1,',',i2,',',i5,'>=',2(1pd12.4))
      end
      subroutine rakdd_printgrad(nint,grad,msg)
      implicit none
#include "stdio.fh"
#include "inp.fh"
      integer nint
      character*(*) msg
      double precision grad(nint,3,4)
*
      integer iint,ic,ia
*
      do ia = 1,4
        do ic = 1,3
          do iint = 1, nint
            write(luout,10000)iint,ic,ia,grad(iint,ic,ia),
     &          msg(1:inp_strlen(msg))
          enddo
        enddo
      enddo
*
10000 format(1x,'grad(',i4,',',i2,',',i4,') =',1pd14.5,a)
*
      end
