      subroutine intp_mpole(i_basis, ish, j_basis, jsh, R,
     &    lmax, centerl,
     &    lscr, scr, lmpint, MP)
c $Id: intp_mpole.F,v 1.2 1995-04-26 00:20:31 d3e129 Exp $
c
c routine to compute multipole integrals up to level lmax with the
c jsh translated by the fractional coordinate vector R.
c The general form is <shell|pole|shell>
c
c  Integrals are returned in shell blocks of <ish|L|jsh> L=0 to lmax
c  one block for each L value.  
c  for ish = p and L = 1 and jsh = p you would get:
c      (3*1*3)+(3*3*3)=36 integrals
c  order would be 
c   <x|x>, <x|y>, <x|z>, <y|x> ... <z|z> (first nine)
c   <x|x|x>, <x|x|y>, <x|x|z>, <x|y|x>, <x|y|y> ... <x|z|z> (2nd nine)
c   <y|x|x>, <y|x|y>, <y|x|z>, <y|y|x>, <y|y|y> ... <y|z|z> (3rd nine)
c   <z|x|x>, <z|x|y>, <z|x|z>, <z|y|x>, <z|y|y> ... <z|z|z> (rth nine)
c
c
      implicit none
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "stdio.fh"
c
c::functions
      logical int_chk_init
      integer int_nint
      external int_chk_init
      external int_nint
c::passed
      integer i_basis             ! [input] basis set handle for ish
      integer ish                 ! [input] i shell/contraction
      integer j_basis             ! [input] basis set handle for jsh
      integer jsh                 ! [input] j shell/contraction
c...     translation vectors are in fractional coordinates 
      double precision R(3)       ! [input] translation vec on j cont.
      integer lmax                ! [input] maximum lvalue for 
*.......................................... multipole integrals 
*.......................................... in this batch
      double precision centerl(3) ! [input] coordinates of multipole
      integer lscr                ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      integer lmpint              ! [input/output] length of multipole 
*................................................. integrals array
      double precision MP(lmpint) ! [output] multipole integrals
c::local      
      logical shells_ok
      integer ibas, Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer jbas, Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer ucont, lval
      integer l_int, ij_int, num_int
      integer int_point, int_size
      double precision xyz_new(3)
      logical inline_chk_sh      
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization
c      
      if (.not.int_chk_init('intp_mpole'))
     &       call errquit('intp_mpole: int_init was not called' ,0)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not. shells_ok)
     &       call errquit('intp_mpole: invalid contraction/shell',0)
c
      ucont   = ibs_cn2ucn(ish,ibas)
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = ibs_cn2ce(ish,ibas)
      i_geom  = ibs_geom(ibas)
c
      ucont   = ibs_cn2ucn(jsh,jbas)
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = ibs_cn2ce(jsh,jbas)
      j_geom  = ibs_geom(jbas)
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)'intp_mpole: WARNING: possible geometry',
     &      ' inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
      call intp_txyz(j_cent,j_geom,R,xyz_new)
c
      if (i_gen.gt.1 .or. j_gen.gt.1) then
        write(luout,*)
     &      ' hf3ois does not handle general contractions yet'
        call errquit('intp_mpole: general contraction error ',911)
      endif
c
      num_int = 0
      do 00100 lval = 0,lmax
        l_int  = (lval+1)*(lval+2)/2
        ij_int = int_nint(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        int_point = num_int + 1
        int_size  = lmpint - num_int
        num_int = num_int + l_int*ij_int
        if (num_int.gt.lmpint) then
          write(luout,*)' intp_mpole: lmpint  = ',lmpint
          write(luout,*)' intp_mpole: num_int = ',num_int
          call errquit('intp_mpole: lmpint too small ',911)
        endif
        call hf3ois(
     &       coords(1,i_cent,i_geom),exndcf(i_iexp,ibas),
     &       exndcf(i_icfp,ibas),i_prim, Li,
     &       centerl,DCexp,
     &       DCcoeff, 1,  lval,
     &       xyz_new,exndcf(j_iexp,jbas),
     &       exndcf(j_icfp,jbas),j_prim, Lj,
c...................................TriDiag..DryRun
     &       MP(int_point),int_size,.false.,.false.,scr,lscr)
00100 continue
      lmpint = num_int
      end
