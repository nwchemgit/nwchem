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
      integer size, scr_size, maxg1, maxg2, maxs1, maxs2
      logical status
      integer hbuf, kbuf, hscr, kscr
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
      size = (cn_nbf_max**4)*nat3*nat3
c
      call intdd_init(rtdb,1,basis)
      call int_mem_print()
      call int_mem_1e(maxg1,maxs1)
      call int_mem_2e4c(maxg2,maxs2)
*
      write(luout,*)' maxg1 :',maxg1
      write(luout,*)' maxs1 :',maxs1
      write(luout,*)' maxg2 :',maxg2
      write(luout,*)' maxs2 :',maxs2
*
      scr_size = 2*maxs2
      status =
     &    ma_alloc_get(mt_dbl,size,'int buffer',hbuf,kbuf)
      status = status.and.
     &    ma_alloc_get(mt_dbl,scr_size,'scr buffer',hscr,kscr)      
      if (.not.status) stop ' memory alloc failed rak23 (1)'
*
      call raktask_intdd_a(geom,basis,nbf,nshell,cn_nbf_max,
     &    size,dbl_mb(kbuf),
     &    scr_size,dbl_mb(kscr))
*
      call intdd_terminate()
      raktask_intdd = bas_destroy(basis)
      raktask_intdd = raktask_intdd.and.
     &    geom_destroy(geom)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hscr)
      raktask_intdd = raktask_intdd.and.
     &    ma_free_heap(hbuf)
      raktask_intdd = .true.
      end
      subroutine raktask_intdd_a(geom,basis,nbf,nshell,cn_nbf_max,
     &    lbuf,buf,lscr,scr)
      implicit none
#include "stdio.fh"
#include "bas.fh"
      integer geom
      integer basis
      integer nbf
      integer nshell
      integer cn_nbf_max
      integer lbuf
      integer lscr
      double precision buf(lbuf)
      double precision scr(lscr)
*
      integer nzero,ncount,count
      integer iatom, jatom, katom, latom
      integer ish, jsh, ksh, lsh
      integer ilo, ihi, inbf
      integer jlo, jhi, jnbf
      integer klo, khi, knbf
      integer llo, lhi, lnbf
      integer nint
      integer idatom(4)
      integer atoms2move(4)
      integer num_atoms2move, atom1, atom2
      double precision thresh
*
      thresh = 1.0d-07
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
                atoms2move(1) = iatom
                atoms2move(2) = jatom
                atoms2move(3) = katom
                atoms2move(4) = latom
                do atom1 = 1,4
                  if (atoms2move(atom1).gt.0) then
                    do atom2 = (atom1+1),4
                      if (atoms2move(atom1).eq.atoms2move(atom2)) then
                        atoms2move(atom2) = 0
                      endif
                    enddo
                  endif
                enddo
                num_atoms2move = 0
                do atom1 = 1,4
                  if (atoms2move(atom1).ne.0)
     &                num_atoms2move = num_atoms2move + 1
                enddo
                nint = inbf*jnbf*knbf*lnbf*78
                call intdd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &              lscr,scr,lbuf,buf,idatom)
                nzero = 0
                ncount = 0
                do count = 1,nint
                  if (abs(buf(count)).gt.thresh) then
                    ncount = ncount + 1
                  else
                    nzero = nzero + 1
                  endif
                enddo
                if ((nzero+ncount).ne.nint)
     &              write(luout,*)nzero,'+',ncount,' != ',nint
                write(luout,*)
     &              ' shells: ',ish,jsh,ksh,lsh,
     &              ' idatom: ',idatom,
     &              ' ncount: ',ncount,
     &              ' atoms2move: ',atoms2move,num_atoms2move
              endif
            enddo
          enddo
        enddo
      enddo
*
      end
