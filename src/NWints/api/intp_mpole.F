C> \ingroup nwint
C> @{
C>
C> \brief Compute 1-electron Gaussian periodic multipole
C>  integrals for \f$l\f$ from \f$0\f$ to \f$l_{max}\f$
C>
C> See [1] for details.
C>
C> [1] JE Jaffe, AC Hess,
C>     <i>"Gaussian basis density functional theory for systems 
C>     periodic in two or three dimensions: Energy and forces"</i>,
C>    J.Chem.Phys. <b>105</b>, 10983-10998 (1996), DOI:
C>    <a href="https://doi.org/10.1063/1.472866">
C>    10.1063/1.472866</a>
C>
      subroutine intp_mpole(i_basis, ish, j_basis, jsh, R,
     &    lmax, centerl,
     &    lscr, scr, lmpint, MP, num_mpint)
c $Id$
c
c routine to compute multipole integrals up to level lmax with the
c jsh translated by the fractional coordinate vector R.
c The general form is <shell|pole|shell>
c
c  Integrals are returned in shell blocks of <ish|L|jsh> L=0 to lmax
c  one block for each L value.  
c  for ish = p and Lmax = 1 and jsh = p you would get:
c      (3*1*3)+(3*3*3)=36 integrals
c  order would be 
c   <x|x>, <x|y>, <x|z>, <y|x> ... <z|z> (first nine)
c   <x|x|x>, <x|x|y>, <x|x|z>, <x|y|x>, <x|y|y> ... <x|z|z> (2nd nine)
c   <y|x|x>, <y|x|y>, <y|x|z>, <y|y|x>, <y|y|y> ... <y|z|z> (3rd nine)
c   <z|x|x>, <z|x|y>, <z|x|z>, <z|y|x>, <z|y|y> ... <z|z|z> (4th nine)
c
c  for ish = p and Lmax = 1 and jsh = d you would get:
c      (3*1*6) + (3*3*6) = 18 + 54 = 72 integrals
c  order would be 
c   <x|xx>, <x|xy>, <x|xz>, ... , <x|zz>, <y|xx>, ... <z|zz> (first 18)
c   <x|x|xx>, <x|x|xy>, <x|x|xz>, ... , <x|x|zz>, <x|y|xx>, ... <x|z|zz> (second 18)
c   <y|x|xx>, <y|x|xy>, <y|x|xz>, ... , <y|x|zz>, <y|y|xx>, ... <y|z|zz> (third 18)
c   <z|x|xx>, <z|x|xy>, <z|x|xz>, ... , <z|x|zz>, <z|y|xx>, ... <z|z|zz> (fourth 18)
c
c
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c
c::functions
      logical int_chk_init
      integer int_nint_cart, int_nint
      external int_chk_init
      external int_nint_cart, int_nint
c::passed
      integer i_basis             !< [Input] basis set handle for ish
      integer ish                 !< [Input] i shell/contraction
      integer j_basis             !< [Input] basis set handle for jsh
      integer jsh                 !< [Input] j shell/contraction
c...     translation vectors are in fractional coordinates 
      double precision R(3)       !< [Input] translation vec on j cont.
      integer lmax !< [Input] maximum lvalue for 
                   !< multipole integrals in this batch
      double precision centerl(3) !< [Input] coordinates of multipole
      integer lscr                !< [Input] length of scratch array
      double precision scr(lscr)  !< [Input] scratch array
      integer lmpint              !< [Input] length of multipole 
*.......................................... integrals array
      double precision MP(lmpint) !< [Output] multipole integrals
      integer num_mpint           !< [Output] number of multipole integrals
c::local      
      logical shells_ok
      integer ibas, Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer jbas, Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer i_gen_save, j_gen_save
      integer ucont, lval
      integer l_int, ij_int, ij_int_x, num_int
      integer int_point, int_size
      integer lpole
      double precision xyz_new(3)
      logical any_spherical
      logical inline_chk_sh      
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization
c      
      if (.not.int_chk_init('intp_mpole'))
     &       call errquit('intp_mpole: int_init was not called' ,0,
     &          INT_ERR)
c
c  check if gencon/sp shells
c
      call int_nogencont_check(i_basis,'intp_mpole:i_basis')
      call int_nogencont_check(j_basis,'intp_mpole:j_basis')
      call int_nospshell_check(i_basis,'intp_mpole:i_basis')
      call int_nospshell_check(j_basis,'intp_mpole:j_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not. shells_ok)
     &       call errquit('intp_mpole: invalid contraction/shell',0,
     &         BASIS_ERR)
c
***   set defNxyz such that it can handle the maximum multi-pole
c
      lpole = lmax/4 + 1
      call defNxyz(lpole)
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
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
        call errquit('intp_mpole: general contraction error ',911,
     &           BASIS_ERR)
      endif
c
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      num_int = 0
      do 00100 lval = 0,lmax
        l_int    = (lval+1)*(lval+2)/2
        ij_int   = int_nint(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        ij_int_x = int_nint_cart(i_basis, ish, j_basis, jsh, 0,0, 0,0)
        int_point = num_int + 1
        int_size  = lmpint - num_int
        if ((num_int+l_int*ij_int_x).gt.lmpint) then
          write(luout,*)' intp_mpole: lmpint  = ',lmpint
          write(luout,*)' intp_mpole: (num_int+l_int*ij_int_x) = ',
     &        (num_int+l_int*ij_int_x)
          call errquit('intp_mpole: lmpint too small ',911, INT_ERR)
        endif
        num_int = num_int + l_int*ij_int
        call hf3ois(
     &       coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &       dbl_mb(mb_exndcf(i_icfp,ibas)),i_prim, Li,
     &       centerl,DCexp,
     &       DCcoeff, 1,  lval,
     &       xyz_new,dbl_mb(mb_exndcf(j_iexp,jbas)),
     &       dbl_mb(mb_exndcf(j_icfp,jbas)),j_prim, Lj,
c...................................TriDiag..DryRun
     &       MP(int_point),int_size,.false.,.false.,scr,lscr)
        if (any_spherical) then
          if (Li.eq.-1) then
            i_gen_save = i_gen
            i_gen = 1
          endif
          if (Lj.eq.-1) then
            j_gen_save = i_gen
            j_gen = 1
          endif
          call spcart_3ctran(MP(int_point),scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(lval),int_nbf_x(lval),lval,1,.false.,
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        .false.)
          if (Li.eq.-1) i_gen = i_gen_save
          if (Lj.eq.-1) j_gen = j_gen_save
        endif
00100 continue
      num_mpint = num_int
      end
C> @}
