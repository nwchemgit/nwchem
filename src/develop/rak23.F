      logical function raktask_intdd(rtdb)
* test the second derivative integrals api
      implicit none
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
      if (.not.bas_numcont(basis,nshell)) call errquit
     &    ('numbf failed',911)
c
      if (.not.geom_ncent(geom,nat)) stop 'geom_ncent fe'
      write(6,*) 'number of atoms ', nat
      nat3 = 3*nat
c
      if (.not.bas_nbf_cn_max(basis,cn_nbf_max))
     &    stop 'bas_nbf_cn_max'
c
      size   = (cn_nbf_max**4)*78
      sizesq = (cn_nbf_max**4)*12*12
      sizeg  = (cn_nbf_max**4)*12
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
#include "stdio.fh"
#include "bas.fh"
#include "geomP.fh"
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
      double precision thresh, delta, scale
*
* store original coordintates
      nat3 = 3*nat
      call dcopy(nat3,coords(1,1,geom),1,xyz,1)
*
      thresh = 1.0d-07
      delta  = 0.01d00
*
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi))
     &      stop 'cn2bfr error i'
        inbf = ihi - ilo + 1
        if (.not.bas_cn2ce(basis,ish,iatom))
     &      stop 'bas_cn2ce error i'
        do jsh = 1,ish
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi))
     &        stop 'cn2bfr error j'
          jnbf = jhi - jlo + 1
          if (.not.bas_cn2ce(basis,jsh,jatom))
     &        stop 'bas_cn2ce error j'
          do ksh = 1,jsh
            if (.not.bas_cn2bfr(basis,ksh,klo,khi))
     &          stop 'cn2bfr error k'
            knbf = khi - klo + 1
            if (.not.bas_cn2ce(basis,ksh,katom))
     &          stop 'bas_cn2ce error k'
            do lsh = 1,ksh
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi))
     &            stop 'cn2bfr error l'
              lnbf = lhi - llo + 1
              if (.not.bas_cn2ce(basis,lsh,latom))
     &            stop 'bas_cn2ce error l'
              if (.not.(iatom.eq.jatom.and.
     &                           jatom.eq.katom.and.
     &                                    katom.eq.latom)) then
                nint  = inbf*jnbf*knbf*lnbf
                ninth = nint*78
                call intdd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &              lscr,scr,lbuf,buf,idatom)
                call rakdd_print_dd(nint,buf,idatom)
                idatoms(1) = idatom(1)
                idatoms(2) = idatom(2)
                idatoms(3) = idatom(3)
                idatoms(4) = idatom(4)
                call rakdd_fill(12,nint,buf,bufsum,idatoms)
                call rakdd_printsum(nint,bufsum,idatoms)
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
00001           continue
                do atom1 = 1,4
                  if (atoms2move(atom1).eq.0) then
                    do atom2 = atom1,3
                      atoms2move(atom2) = atoms2move(atom2+1)
                    enddo
                    atoms2move(4)=0
                  endif
                enddo
                do atom1 = 1,num_atoms2move
                  if (atoms2move(atom1).eq.0) goto 00001
                enddo
*
                write(luout,10000)
     &              'shells: ',ish,jsh,ksh,lsh,
     &              'idatom: ',idatom,
     &              'ncount: ',ncount,
     &              'atoms2move: ',atoms2move,'num',num_atoms2move
*
                do zatom = 1,4
                  if (atoms2move(zatom).gt.0) then
                    do ixyz = 1,3
                      nintg = nint*12
                      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                      coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) + delta
                      call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &                    lscr,scr,lgrad,gradp,idatomp)
                      call dcopy(nat3,xyz,1,coords(1,1,geom),1)
                      coords(ixyz,atoms2move(zatom),geom) =
     &                    coords(ixyz,atoms2move(zatom),geom) - delta
                      call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &                    lscr,scr,lgrad,gradm,idatomm)
                      call daxpy(nintg,-1.0d00,gradm,1,gradp,1)
                      scale = 1.0d00/(2.0d00*delta)
                      call dscal(nintg,scale,gradp,1)
*                      call rakdd_gradfill(gradp,idatomp,ixyz,
*     &                    atoms2move(zatom))
                    enddo
                  endif
                enddo
              endif ! 4 atoms the same
            enddo   ! lsh
          enddo     ! ksh
        enddo       ! jsh
      enddo         ! ish
*
10000 format(1x,a,4(i4),1x,a,4(i3),1x,a,i6,1x,a,4(i3),1x,a,i2)
      end
      subroutine rakdd_fill(at4xyz,nint,buf,bufsum,idatom)
      implicit none
      integer at4xyz, nint
      integer idatom(4)
      double precision buf(78,nint)
      double precision bufsum(12,12,nint)
c
      integer itri
      integer ia,ja,ixyz,jxyz,i,j,int
      logical doit
c
      itri = 0
      do ia = 1,4
        do ja = 1,4
          if (ia.le.ja) then
            do ixyz = 1,3
              i = (ixyz-1)*4 + ia
              do jxyz = 1,3
                j = (jxyz-1)*4 + ja
                doit = .true.
                if (ia.eq.ja) doit = ixyz.le.jxyz
                if (doit) then
                  itri = itri + 1
                  do int = 1,nint
                    bufsum(i,j,int) = buf(itri,int)
                    bufsum(j,i,int) = buf(itri,int)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo
c
      call rakdd_fill_a(bufsum,idatom,nint)
      end
      subroutine rakdd_fill_a(buf,idatom,nint)
      implicit none
      integer idatom(4)
      integer nint
      double precision buf(3,4,3,4,nint)
*
      logical sum12, sum13, sum14, sum23, sum24, sum34
      integer ixyz, jxyz, iatom, jatom, int
*
      if (idatom(1).gt.0.and.idatom(2).gt.0) then
        sum12 = idatom(1).eq.idatom(2)
        if (sum12) then
          idatom(2) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,1,jxyz,jatom,int) =
     &                buf(ixyz,1,jxyz,jatom,int) +
     &                buf(ixyz,2,jxyz,jatom,int) 
                  buf(ixyz,2,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,1,int) =
     &                buf(ixyz,iatom,jxyz,1,int) +
     &                buf(ixyz,iatom,jxyz,2,int) 
                  buf(ixyz,iatom,jxyz,2,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(1).gt.0.and.idatom(3).gt.0) then
        sum13 = idatom(1).eq.idatom(3)
        if (sum13) then
          idatom(3) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,1,jxyz,jatom,int) =
     &                buf(ixyz,1,jxyz,jatom,int) +
     &                buf(ixyz,3,jxyz,jatom,int) 
                  buf(ixyz,3,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,1,int) =
     &                buf(ixyz,iatom,jxyz,1,int) +
     &                buf(ixyz,iatom,jxyz,3,int) 
                  buf(ixyz,iatom,jxyz,3,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(1).gt.0.and.idatom(4).gt.0) then
        sum14 = idatom(1).eq.idatom(4)
        if (sum14) then
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,1,jxyz,jatom,int) =
     &                buf(ixyz,1,jxyz,jatom,int) +
     &                buf(ixyz,4,jxyz,jatom,int) 
                  buf(ixyz,4,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,1,int) =
     &                buf(ixyz,iatom,jxyz,1,int) +
     &                buf(ixyz,iatom,jxyz,4,int) 
                  buf(ixyz,iatom,jxyz,4,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(2).gt.0.and.idatom(3).gt.0) then
        sum23 = idatom(2).eq.idatom(3)
        if (sum23) then
          idatom(3) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,2,jxyz,jatom,int) =
     &                buf(ixyz,2,jxyz,jatom,int) +
     &                buf(ixyz,3,jxyz,jatom,int) 
                  buf(ixyz,3,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,2,int) =
     &                buf(ixyz,iatom,jxyz,2,int) +
     &                buf(ixyz,iatom,jxyz,3,int) 
                  buf(ixyz,iatom,jxyz,3,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(2).gt.0.and.idatom(4).gt.0) then
        sum24 = idatom(2).eq.idatom(4)
        if (sum24) then
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,2,jxyz,jatom,int) =
     &                buf(ixyz,2,jxyz,jatom,int) +
     &                buf(ixyz,4,jxyz,jatom,int) 
                  buf(ixyz,4,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,2,int) =
     &                buf(ixyz,iatom,jxyz,2,int) +
     &                buf(ixyz,iatom,jxyz,4,int) 
                  buf(ixyz,iatom,jxyz,4,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      if (idatom(3).gt.0.and.idatom(4).gt.0) then
        sum34 = idatom(3).eq.idatom(4)
        if (sum34) then
          idatom(4) = 0
          do ixyz = 1,3
            do jxyz = 1,3
              do jatom = 1,4
                do int = 1,nint
                  buf(ixyz,3,jxyz,jatom,int) =
     &                buf(ixyz,3,jxyz,jatom,int) +
     &                buf(ixyz,4,jxyz,jatom,int) 
                  buf(ixyz,4,jxyz,jatom,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
          do ixyz = 1,3
            do iatom = 1,4
              do jxyz = 1,3
                do int = 1,nint
                  buf(ixyz,iatom,jxyz,3,int) =
     &                buf(ixyz,iatom,jxyz,3,int) +
     &                buf(ixyz,iatom,jxyz,4,int) 
                  buf(ixyz,iatom,jxyz,4,int) = 0.0d00
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      end
      subroutine rakdd_printsum(nint,bufsum,idatom)
      implicit none
      integer nint
      integer idatom(4)
      double precision bufsum(3,4,3,4,nint)
c
      integer ixyz,iatom,jxyz,jatom,int
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
              do int = 1,nint
                val = bufsum(ixyz,iatom,jxyz,jatom,int)
                if (abs(val).gt.thresh) then
                  count = count + 1
                  write(6,10000)
     &                atoms(iatom),sxyz(ixyz),
     &                atoms(jatom),sxyz(jxyz),
     &                ixyz,iatom,jxyz,jatom,int,
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
      integer nint
      integer idatom(4)
      double precision buf(78,nint)
c
      double precision thresh, val
      integer ider,int, count
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
        do int = 1,nint
          val = buf(ider,int)
          if (abs(val).gt.thresh) then
            count = count + 1
            write(6,10000)sder(ider),
     &          ider,int,val,idatom,count
          endif
        enddo
      enddo
10000 format(1x,'buf(',a4,','i4,',',i4,') =',1pd20.10,
     &    2x,4(i5),1x,i8)
      end
