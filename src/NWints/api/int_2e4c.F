      subroutine int_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri)
      implicit none
c
c basic api routine to generate a block of two electron integrals
c eri = <bra_g(ish).bra_g(jsh) | ket_g(ksh).ket_g(lsh)>
c
#include "apiP.fh"
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::functions
      logical int_check_shell
      external int_check_shell
c:: passed
      integer brain, ketin, lscr, leri
      integer ish,jsh,ksh,lsh, i,j
      double precision scr(lscr), eri(leri)
c:: local
      logical status
      integer bra, ket, ijsh, klsh
      integer ab_geom, cd_geom, ucont
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_ucent
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_ucent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_ucent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_ucent
      integer isym2
c:statement function
      isym2(i,j)=(max(i,j)*(max(i,j)-1))/2 + min(i,j)
c
      if (init_int.eq.0) then
        write(6,*)' int_2e4c: ERROR '
        write(6,*)' int_init was not called '
        stop ' error '
      endif
c
c     check input shell ranges
c
      status = int_check_shell(brain,ish)
      status = status .and. int_check_shell(brain,jsh)
      status = status .and. int_check_shell(ketin,ksh)
      status = status .and. int_check_shell(ketin,lsh)
      if (.not.status) then
        write(6,*)'int_2e4c: invalid shell input '
        stop 'error'
      endif
c
c   check canonicalizations of input shells
c      
      status = ish.gt.jsh
      status = status .and. ksh.gt.lsh
      ijsh = isym2(ish,jsh)
      klsh = isym2(ksh,lsh)
      if (brain.eq.ketin) status = status .and. ijsh.gt.klsh
      if (.not.status) then
        write(6,*)'int_2e4c: shells not canonical on input '
        write(6,*)'bra basis set handle:',brain
        write(6,*)'ket basis set handle:',ketin
        write(6,*)'                 ish:',ish
        write(6,*)'                 jsh:',jsh
        write(6,*)'                 ksh:',ksh
        write(6,*)'                 lsh:',lsh
        write(6,*)'                ijsh:',ijsh
        write(6,*)'                klsh:',klsh
        stop ' error '
      endif
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
c
      ab_geom  = ibs_geom(bra)
c
      ucont = ibs_cn2ucn(ish,bra)
      La      = infbs_cont(CONT_TYPE ,ucont,bra)
      a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
      a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
      a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
      a_ucent = ibs_ce2uce(ibs_cn2ce(ish,bra),bra)
c
      ucont = ibs_cn2ucn(jsh,bra)
      Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
      b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
      b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
      b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
      b_ucent = ibs_ce2uce(ibs_cn2ce(jsh,bra),bra)
c      
      ket = ketin + BASIS_HANDLE_OFFSET
c
      cd_geom  = ibs_geom(ket)
c
      ucont   = ibs_cn2ucn(ksh,ket)
      Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
      c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
      c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
      c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
      c_ucent = ibs_ce2uce(ibs_cn2ce(ksh,ket),ket)
c
      ucont   = ibs_cn2ucn(lsh,ket)
      Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
      d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
      d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
      d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
      d_ucent = ibs_ce2uce(ibs_cn2ce(lsh,ket),ket)
c      
c
      call hf2(
     &       coords(1,a_ucent,ab_geom), exndcf(a_iexp,bra),
     &       exndcf(a_icfp,bra), a_prim, a_gen, La,
     &       coords(1,b_ucent,ab_geom), exndcf(b_iexp,bra),
     &       exndcf(b_icfp,bra), b_prim, b_gen, Lb,
     &       coords(1,c_ucent,cd_geom), exndcf(c_iexp,ket),
     &       exndcf(c_icfp,ket), c_prim, c_gen, Lc,
     &       coords(1,d_ucent,cd_geom), exndcf(d_iexp,ket),
     &       exndcf(d_icfp,ket), d_prim,d_gen,Ld,
     &       eri, leri, .false., .false., .false.,
     &       .false., scr, lscr)
c
      end
