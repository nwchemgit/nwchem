c**********************************************************************
c     sunbroutines interfaced with and adapted from HONDO
c $Id: rhf_hessian.F,v 1.15 1998-06-05 16:55:22 d3j191 Exp $
c**********************************************************************
      subroutine hnd_der_init(geom,basis)
      implicit double precision (a-h,o-z)
c 
      parameter (zero=0.0d+00)
      parameter (mxatom=500)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      parameter (mxnbf =2048)
      parameter (mxtype=5)
      logical  status
      integer  geom 
      integer  basis
      integer  num
      logical  geom_ncent
      external geom_ncent
      logical  geom_create
      external geom_create
      logical  geom_rtdb_load
      external geom_rtdb_load
      logical  geom_cart_get
      external geom_cart_get
      logical  bas_numcont
      external bas_numcont
      logical  bas_numbf
      external bas_numbf
      logical  bas_ce2cnr
      external bas_ce2cnr
      logical  bas_cn2bfr
      external bas_cn2bfr
      character*16 atom_tags
c
      logical dbug
      logical out
      dimension atom_tags(mxatom)
      dimension mini(mxtype),maxi(mxtype)
c
c     ----- routine initializes HONDO common blocks -----
c
      common/inp_file/iread,iwrite
      common/hnd_iofile/ir,iw,ip
      common/hnd_baspar/normf,normp,itol
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_ijpair/ia(mxnbf)
      common/hnd_nshel/ex(mxprim),cs(mxprim),cp(mxprim),
     1                 cd(mxprim),cf(mxprim),cg(mxprim),
     2                 kstart(mxshel),katom(mxshel),ktype(mxshel),
     3                 kng(mxshel),kloc(mxshel),kmin(mxshel),
     4                 kmax(mxshel),nshell
      common/hnd_shlord/modshl(mxshel),invshl(mxshel)
      common/hnd_lcapid/nap,iap
      common/hnd_symtry/invt(48),nt,ntmax,ntwd,nosym
      common/hnd_symshl/mapshl(mxshel,48),mapctr(mxatom,48)
c
#include "global.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"

      data mini / 1, 2, 5,11,21/
      data maxi / 1, 4,10,20,35/

#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      dbug=.false.
      out =.false.
      ir=iread
      iw=iwrite
      ip=7
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... hnd_der_init ... '
      endif
      if(dbug) then
         write(iw,*) ' ... geom  = ',geom
         write(iw,*) ' ... basis = ',basis
      endif
c
c     ----- get basis set info + parameters -----
c
      if (.not. geom_ncent (geom , natoms)) call errquit
     $     ('hnd_der_init: geom_ncent  failed for geom', geom)
      if(dbug) then
         write(iw,*) ' ... natoms = ', natoms
      endif
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hnd_der_init: bas_numcont failed for basis', basis)

      if (.not. bas_numbf(basis, num)) call errquit
     $     ('hnd_der_init: bas_numbf failed for basis', basis)
      if(dbug) then
         write(iw,*) ' ... nshell = ', nshell
         write(iw,*) ' ... num    = ', num
      endif

      do i=1,mxnbf
         ia(i) = i*(i-1)/2 
      enddo
c
c     ----- loop over atoms -----
c
      istart=1
      do iat=1,natoms
         status=bas_ce2cnr(basis,iat,ishlo,ishhi)
         if(.not.status) call errquit
     $       ('hnd_der_init: bas_ce to cnr/bfr failed', iat)
         if(dbug) then
            write(iw,*) ' ... iat, ishlo, ishhi = ',iat,ishlo,ishhi
         endif
c 
c     ----- loop over shells on atoms -----
c
         do ish=ishlo,ishhi
            status=bas_cn2bfr(basis,ish,ishbflo,ishbfhi)
            if(.not.status) call errquit
     $          ('hnd_der_init: bas_cn to bfr failed', ish)
            ioff=basis+BASIS_HANDLE_OFFSET
            ishu=sf_ibs_cn2ucn(ish,ioff)
            ityp=infbs_cont(CONT_TYPE ,ishu,ioff)
            ngi =infbs_cont(CONT_NPRIM,ishu,ioff)
            icon=infbs_cont(CONT_ICFP ,ishu,ioff)
            iexp=infbs_cont(CONT_IEXP ,ishu,ioff)
c
            katom(ish) =iat
            kloc(ish)  =ishbflo
            ktype(ish) =ityp+1
            kmin(ish)  =mini(ityp+1)
            kmax(ish)  =maxi(ityp+1)
            kng(ish)   =ngi
            kstart(ish)=istart
c
c     ----- loop over primitives in shells -----
c
            do ig=1,ngi
               ex(ig+istart-1)=zero
               cs(ig+istart-1)=zero
               cp(ig+istart-1)=zero
               cd(ig+istart-1)=zero
               cf(ig+istart-1)=zero
               cg(ig+istart-1)=zero
               cont_expo     =dbl_mb(mb_exndcf(ig+iexp-1,ioff))
               cont_coeff    =dbl_mb(mb_exndcf(ig+icon-1,ioff))
               ex(ig+istart-1)=cont_expo
               if(ityp.eq.0) then
                  cs(ig+istart-1)=cont_coeff
               elseif(ityp.eq.1) then
                  cp(ig+istart-1)=cont_coeff
               elseif(ityp.eq.2) then
                  cd(ig+istart-1)=cont_coeff
               elseif(ityp.eq.3) then
                  cf(ig+istart-1)=cont_coeff
               elseif(ityp.eq.4) then
                  cg(ig+istart-1)=cont_coeff
               else
                  call errquit('hnd_der_init: illegal -ityp- ', ityp)
               endif
               if(dbug) then
                  write(iw,*) ' ish,ig,ityp,exp,con= ',ish,ig,ityp,
     $                         cont_expo,cont_coeff
               endif
            enddo
            istart=istart+ngi
c
         enddo
c
      enddo
c
c     ----- print /nshel/ -----
c
      if(dbug) then
         write(iw,*) ' ... kstart= ',(kstart(ish),ish=1,nshell)
         write(iw,*) ' ... katom = ',( katom(ish),ish=1,nshell)
         write(iw,*) ' ... ktype = ',( ktype(ish),ish=1,nshell)
         write(iw,*) ' ... kng   = ',(   kng(ish),ish=1,nshell)
         write(iw,*) ' ... kloc  = ',(  kloc(ish),ish=1,nshell)
         write(iw,*) ' ... kmin  = ',(  kmin(ish),ish=1,nshell)
         write(iw,*) ' ... kmax  = ',(  kmax(ish),ish=1,nshell)
      endif
c
c     ----- /baspar/ -----
c
      normf=0
      normp=0
      itol=15
      if(dbug) then
         write(iw,*) ' ... normf,normp,itol = ',normf,normp,itol
      endif
c
c     ----- /molxyz/ -----
c
      status=geom_ncent(geom,nat)
      status=geom_cart_get(geom,nat,atom_tags,c,zan)
      if(.not.status) call errquit
     $    (' int_1e_ga_hnd_init: geom_cart_get failed.',911)
      if(dbug) then
         do iat=1,nat
            write(iw,*) ' ... i,Z,c = ',iat,zan(iat),c(1,iat),c(2,iat),
     $                                               c(3,iat)
         enddo
      endif
c
c     ----- initialization for 2e integral derivatives -----
c
      do i=1,mxshel
         modshl(i)=i
         invshl(i)=i
      enddo
c
c     nap = 1
c     iap = 1
      nap = ga_nnodes()
      iap = ga_nodeid() + 1
c
      nt=1
      do i=1,mxshel
         mapshl(i,1)=i
      enddo
c
      do i=1,mxatom
         mapctr(i,1)=i
      enddo
c
c     ----- load Hermite's roots and weights -----
c
      call hnd_whermt
c
      return
      end
c
      subroutine hnd_nw_hnd
      implicit double precision (a-h,o-z)
c
c     ----- order of components in shells in HONDO -----
c
c            x    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
c    1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
c    2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
c    3              3, 1, 3, 2, 2/
c            y    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
c    1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
c    2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
c    3              1, 3, 2, 3, 2/
c            z    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
c    1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
c    2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
c    3              3, 3, 2, 2, 3/
c
c     ----- order of components in shells in NWCHEM -----
c
c            x    / 1, 2, 1, 1, 3, 2, 2, 1, 1, 1,
c    1              4, 3, 3, 2, 2, 2, 1, 1, 1, 1,
c    2              5, 4, 4, 3, 3, 3, 2, 2, 2, 2,
c    3              1, 1, 1, 1, 1/
c            y    / 1, 1, 2, 1, 1, 2, 1, 3, 2, 1,
c    1              1, 2, 1, 3, 2, 1, 4, 3, 2, 1,
c    2              1, 2, 1, 3, 2, 1, 4, 3, 2, 1,
c    3              5, 4, 3, 2, 1/
c            z    / 1, 1, 1, 2, 1, 1, 2, 1, 2, 3,
c    1              1, 1, 2, 1, 2, 3, 1, 2, 3, 4,
c    2              1, 1, 2, 1, 2, 3, 1, 2, 3, 4,
c    3              1, 2, 3, 4, 5/
c
c     ----- order of components in shells in HONDO -----
c
c                    1      2      3      4
c                '   S','   X','   Y','   Z',
c                    5      6      7      8      9     10
c    1           '  XX','  YY','  ZZ','  XY','  XZ','  YZ',
c                   11     12     13     14     15  
c    2           ' XXX',' YYY',' ZZZ',' XXY',' XXZ',
c                   16     17     18     19     20
c    3           ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',
c                   21     22     23     24     25
c    4           'XXXX','YYYY','ZZZZ','XXXY','XXXZ',
c                   26     27     28     29     30
c    5           'YYYX','YYYZ','ZZZX','ZZZY','XXYY',
c                   31     32     33     34     35
c    6           'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
c
c     ----- order of components in shells in NWCHEM -----
c
c                    1      2      3      4
c                '   S','   X','   Y','   Z',
c                    5      6      7      8      9     10
c    1           '  XX','  XY','  XZ','  YY','  YZ','  ZZ',
c                   11     12     13     14     15     
c    2           ' XXX',' XXY',' XXZ',' XYY',' XYZ',
c                   16     17     18     19     20
c    3           ' XZZ',' YYY',' YYZ',' YZZ',' ZZZ',
c                   21     22     23     24     25
c    4           'XXXX','XXXY','XXXZ','XXYY','XXYZ',
c                   26     27     28     29     30
c    5           'XXZZ','XYYY','XYYZ','XYZZ','XZZZ',
c                   31     32     33     34     35
c    6           'YYYY','YYYZ','YYZZ','YZZZ','ZZZZ'/
c
      parameter (mxatom=500)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      parameter (mxnbf =2048)
      parameter (sqrt3 = 1.73205080756888d+00)
      parameter (sqrt5 = 2.23606797749979d+00)
      parameter (sqrt7 = 2.64575131106459d+00)
      parameter (sq5d3 = sqrt5/sqrt3)
      common/hnd_iofile/ir,iw,ip
      common/hnd_nshel/ex(mxprim),cs(mxprim),cp(mxprim),
     1                 cd(mxprim),cf(mxprim),cg(mxprim),
     2                 kstart(mxshel),katom(mxshel),ktype(mxshel),
     3                 kng(mxshel),kloc(mxshel),kmin(mxshel),
     4                 kmax(mxshel),nshell
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_hndtonw/ihnd_to_nw(mxnbf)
      common/hnd_facntoh/fac_nwthnd(mxnbf)
c
#include "global.fh"
c
      logical dbug
      logical out 
c
      integer map_ntoh(35)
      double precision fac_ntoh(35)
c
c     ----- mapping of shell components from nwchem to hondo -----
c
      data map_ntoh /1,     2,     3,     4,
     1               5,     8,     9,     6,    10,    7,
     2              11,    14,    15,    16,    20,
     3              18,    12,    17,    19,    13,  
     4              21,    24,    25,    30,    33,
     5              31,    26,    34,    35,    28,
     6              22,    27,    32,    29,    23 /

      data fac_ntoh /1,     1,     1,     1,
     1               1, sqrt3, sqrt3,     1, sqrt3,    1,
     2               1, sqrt5, sqrt5, sqrt5, sqrt3,
     3           sqrt5,     1, sqrt5, sqrt5,     1,  
     4               1, sqrt7, sqrt7, sq5d3, sqrt3,
     5           sq5d3, sqrt7, sqrt3, sqrt3, sqrt7,
     6               1, sqrt7, sq5d3, sqrt7,     1/
c
      dbug=.false.
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nw_hnd ... '
      endif
c
c     ----- loop over shells -----
c
     
      do ish = 1,nshell
         mini =  kmin(ish)
         maxi =  kmax(ish)
         loci =  kloc(ish)
         do im = mini,maxi

c     ----- pointer of component of shell in nwchem ----

            inw  = loci + im - mini

c     ----- pointer of component of shell in hondo -----

            ihnd = loci + map_ntoh(im) -mini

c     ----- switch from nwchem to hondo -----

            inw_to_hnd(inw ) = ihnd
            fac_nwthnd(inw)  = fac_ntoh(im)

c     ----- switch from hondo to nwchem -----

            ihnd_to_nw(ihnd) = inw 

         enddo
      enddo
c
      return
      end

      subroutine hnd_nwhnd_tran(dnwc,dhnd,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_iofile/ir,iw,ip
      common/hnd_facntoh/fac_nwthnd(mxnbf)
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nwhnd_tran ... '
      endif
c
c     ----- matrices tran. from nwchem to hondo shell order-----
c
      do j=1,ndim
        do i=1,ndim

           dhnd(inw_to_hnd(i),inw_to_hnd(j))=
     &     dnwc(i,j)/(fac_nwthnd(i)*fac_nwthnd(j))

        enddo
      enddo
c
      return 
      end
c
      subroutine hnd_hndnw_tran(dhnd,dnwc,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_hndtonw/ihnd_to_nw(mxnbf)
      common/hnd_iofile/ir,iw,ip
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
      common/hnd_facntoh/fac_nwthnd(mxnbf)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... hndnw_tran ... '
      endif
c
c     ----- matrices tran. from hondo to nwchem shell order-----
c
      do j=1,ndim
        do i=1,ndim
           dnwc(ihnd_to_nw(i),ihnd_to_nw(j))=
     &     dhnd(i,j)*(
     &     1.0d+00/(fac_nwthnd(ihnd_to_nw(i))*
     &     fac_nwthnd(ihnd_to_nw(j) )) )
        enddo
      enddo
c
      return 
      end
c
      subroutine hnd_int_1e_stvdd(rtdb)
      implicit double precision (a-h,o-z)
#include "global.fh"
      integer rtdb
c
c     ----- routine adapted from HONDO's -STVDD- -----
c
      logical dbug
      logical out 
      dbug=.false.
      out =.false.
      if(out.and.ga_nodeid().eq.0) then
         write(6,*) ' in routine ... hnd_int_1e_stvdd ... '
      endif
c
c     ----- routine -stvdd- -----
c
      call hnd_stvdd(rtdb)
c
      return
      end
c
      subroutine hnd_int_2e_jkdd(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "global.fh"
c
c     ----- routine adapted from HONDO's -JKDD- -----
c
      integer rtdb,geom,basis
      logical dbug
      logical out 
      dbug=.false.
      out =.false.
      if(out.and.ga_nodeid().eq.0) then
         write(6,*) ' in routine ... hnd_int_2e_jkdd ... '
      endif
c
c     ----- routine -jkdd- -----
c
      call hnd_jkdd(rtdb,geom,basis)
c
      return
      end

      subroutine scf_hnd_dder(rtdb,geom,basis)
      implicit none
#include "bas.fh"
#include "util.fh"
#include "global.fh"
c
      integer rtdb, basis, geom
      double precision total_wall,total_cpu
      logical ohondo
c
      ohondo =.true.
c
      if(ohondo) then

c     ----- fill appropriate arrays for HONDO integrals -----

         call hnd_der_init(geom,basis)
   
         call hnd_nw_hnd

c     ----- fill the common block of hnd_molabl

         call hnd_atlabl(rtdb,geom,basis)
       
c     ----- calculate 1e 2nd derivatives -----

         total_wall = -util_wallsec()
         total_cpu  = -util_cpusec()

         call hnd_int_1e_stvdd(rtdb)
       
         total_wall = total_wall + util_wallsec()
         total_cpu  = total_cpu  + util_cpusec()
         if (ga_nodeid().eq.0) write(6,1) total_cpu, total_wall
 1       format(/' Times in 1e_stvdd  cpu:',
     &                f11.1,'s     wall:',f11.1,'s')
         call util_flush(6)

c      ----- calculate 2e 2nd derivatives ----- 

         total_wall = -util_wallsec()
         total_cpu  = -util_cpusec()

         call hnd_int_2e_jkdd(rtdb,geom,basis)

         total_wall = total_wall + util_wallsec()
         total_cpu  = total_cpu  + util_cpusec()
         if (ga_nodeid().eq.0) write(6,2) total_cpu, total_wall
 2       format(/' Times in 2e_jkdd   cpu:',
     &                f11.1,'s     wall:',f11.1,'s')
         call util_flush(6)

c      ------ cphf part -----

         total_wall = -util_wallsec()
         total_cpu  = -util_cpusec()

         call hnd_cphfx(rtdb,geom,basis)

         total_wall = total_wall + util_wallsec()
         total_cpu  = total_cpu  + util_cpusec()
         if (ga_nodeid().eq.0) write(6,3) total_cpu, total_wall
 3       format(/' Times in cphfx     cpu:',
     &                f11.1,'s     wall:',f11.1,'s')
         call util_flush(6)

c
         call hnd_hessian_tidy
       
      else
       
         print*, 'skip HONDO routines'

      endif
c
      return
      end
c
      subroutine hnd_hessian_tidy
      implicit none
      integer nft17
      integer nft18
#include "global.fh"
      logical out 

      out =.false.
      if(out.and.ga_nodeid().eq.0) then
         write(6,*) ' in routine ... hnd_hessian_tidy ... '
      endif
c
      nft17 = 17
      nft18 = 18
c
      close(unit=nft17,status='delete')
      close(unit=nft18,status='delete')
c
      return
      end
C
      SUBROUTINE HND_WHERMT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_HERMIT/H(28)
      COMMON/HND_WERMIT/W(28)
      DATA H11
     1 / 0.000000000000000D+00/
      DATA H21,H22
     1 /-0.707106781186548D+00, 0.707106781186548D+00/
      DATA H31,H32,H33
     1 /-1.224744871391589D+00, 0.000000000000000D+00,
     2   1.224744871391589D+00/
      DATA H41,H42,H43,H44
     1 /-1.650680123885785D+00,-0.524647623275290D+00,
     2   0.524647623275290D+00, 1.650680123885785D+00/
      DATA H51,H52,H53,H54,H55
     1 /-2.020182870456086D+00,-0.958572464613819D+00,
     2   0.000000000000000D+00, 0.958572464613819D+00,
     3   2.020182870456086D+00/
      DATA H61,H62,H63,H64,H65,H66
     1 /-2.350604973674492D+00,-1.335849074013697D+00,
     2  -0.436077411927617D+00, 0.436077411927617D+00,
     3   1.335849074013697D+00, 2.350604973674492D+00/
      DATA H71,H72,H73,H74,H75,H76,H77
     1 /-2.651961356835233D+00,-1.673551628767471D+00,
     2  -0.816287882858965D+00, 0.000000000000000D+00,
     3   0.816287882858965D+00, 1.673551628767471D+00,
     4   2.651961356835233D+00/
      DATA W11
     1 /1.772453850906D+00/
      DATA W21,W22
     1 /8.862269254528D-01,8.862269254528D-01/
      DATA W31,W32,W33
     1 /2.954089751509D-01,1.181635900604D+00,
     2  2.954089751509D-01                   /
      DATA W41,W42,W43,W44
     1 /8.131283544725D-02,8.049140900055D-01,
     2  8.049140900055D-01,8.131283544725D-02/
      DATA W51,W52,W53,W54,W55
     1 /1.995324205905D-02,3.936193231522D-01,
     2  9.453087204829D-01,3.936193231522D-01,
     3  1.995324205905D-02                   /
      DATA W61,W62,W63,W64,W65,W66
     1 /4.530009905509D-03,1.570673203229D-01,
     2  7.246295952244D-01,7.246295952244D-01,
     3  1.570673203229D-01,4.530009905509D-03/
      DATA W71,W72,W73,W74,W75,W76,W77
     1 /9.717812450995D-04,5.451558281913D-02,
     2  4.256072526101D-01,8.102646175568D-01,
     3  4.256072526101D-01,5.451558281913D-02,
     4  9.717812450995D-04/
      H( 1)=H11
      H( 2)=H21
      H( 3)=H22
      H( 4)=H31
      H( 5)=H32
      H( 6)=H33
      H( 7)=H41
      H( 8)=H42
      H( 9)=H43
      H(10)=H44
      H(11)=H51
      H(12)=H52
      H(13)=H53
      H(14)=H54
      H(15)=H55
      H(16)=H61
      H(17)=H62
      H(18)=H63
      H(19)=H64
      H(20)=H65
      H(21)=H66
      H(22)=H71
      H(23)=H72
      H(24)=H73
      H(25)=H74
      H(26)=H75
      H(27)=H76
      H(28)=H77
      W( 1)=W11
      W( 2)=W21
      W( 3)=W22
      W( 4)=W31
      W( 5)=W32
      W( 6)=W33
      W( 7)=W41
      W( 8)=W42
      W( 9)=W43
      W(10)=W44
      W(11)=W51
      W(12)=W52
      W(13)=W53
      W(14)=W54
      W(15)=W55
      W(16)=W61
      W(17)=W62
      W(18)=W63
      W(19)=W64
      W(20)=W65
      W(21)=W66
      W(22)=W71
      W(23)=W72
      W(24)=W73
      W(25)=W74
      W(26)=W75
      W(27)=W76
      W(28)=W77
      RETURN
      END
C
      SUBROUTINE HND_RT123
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C             -----   VERSION FEBRUARY 13,1975   -----
      COMMON/HND_RYS/X,U(5),W(5),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/
    1 IF(X.GT.5.0D+00) GO TO 50
      IF(X.GT.1.0D+00) GO TO 30
      IF(X.GT.3.0D-07) GO TO 20
C     X IS APPROXIMATELY ZERO.         NROOTS=1,2, OR 3
      IF(NROOTS-2) 11,12,13
   11 RT1= 0.5D+00 -X/5.0D+00
      WW1= 1.0D+00 -X/3.0D+00
      RETURN
   12 RT1=1.30693606237085D-01    -2.90430236082028D-02 *X
      RT2=2.86930639376291D+00    -6.37623643058102D-01 *X
      WW1=6.52145154862545D-01    -1.22713621927067D-01 *X
      WW2=3.47854845137453D-01    -2.10619711404725D-01 *X
      RETURN
   13 RT1=6.03769246832797D-02    -9.28875764357368D-03 *X
      RT2=7.76823355931043D-01    -1.19511285527878D-01 *X
      RT3=6.66279971938567D+00    -1.02504611068957D+00 *X
      WW1=4.67913934572691D-01    -5.64876917232519D-02 *X
      WW2=3.60761573048137D-01    -1.49077186455208D-01 *X
      WW3=1.71324492379169D-01    -1.27768455150979D-01 *X
      RETURN
C     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
   20 IF(NROOTS.EQ.3) GO TO 23
      F1=          ((((((((-8.36313918003957D-08*X+1.21222603512827D-06
     1)*X-1.15662609053481D-05 )*X+9.25197374512647D-05
     2)*X-6.40994113129432D-04 )*X+3.78787044215009D-03
     3)*X-1.85185172458485D-02 )*X+7.14285713298222D-02
     4)*X-1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1=(X+X)*F1+DEXP(-X)
      IF(NROOTS.EQ.2) GO TO 22
      RT1=F1/(WW1-F1)
      RETURN
   22 RT1=           (((((((-2.35234358048491D-09*X+2.49173650389842D-08
     1)*X-4.558315364581D-08)*X-2.447252174587D-06)*X+4.743292959463D-05
     2)*X-5.33184749432408D-04 )*X+4.44654947116579D-03
     3)*X-2.90430236084697D-02 )*X+1.30693606237085D-01
      RT2=           (((((((-2.47404902329170D-08*X+2.36809910635906D-07
     1)*X+1.835367736310D-06)*X-2.066168802076D-05)*X-1.345693393936D-04
     2)*X-5.88154362858038D-05 )*X+5.32735082098139D-02
     3)*X-6.37623643056745D-01 )*X+2.86930639376289D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   23 RT1=            ((((((-5.10186691538870D-10*X+2.40134415703450D-08
     1)*X-5.01081057744427D-07 )*X+7.58291285499256D-06
     2)*X-9.55085533670919D-05 )*X+1.02893039315878D-03
     3)*X-9.28875764374337D-03 )*X+6.03769246832810D-02
      RT2=            ((((((-1.29646524960555D-08*X+7.74602292865683D-08
     1)*X+1.56022811158727D-06 )*X-1.58051990661661D-05
     2)*X-3.30447806384059D-04 )*X+9.74266885190267D-03
     3)*X-1.19511285526388D-01 )*X+7.76823355931033D-01
      RT3=            ((((((-9.28536484109606D-09*X-3.02786290067014D-07
     1)*X-2.50734477064200D-06 )*X-7.32728109752881D-06
     2)*X+2.44217481700129D-04 )*X+4.94758452357327D-02
     3)*X-1.02504611065774D+00 )*X+6.66279971938553D+00
      F2=          ((((((((-7.60911486098850D-08*X+1.09552870123182D-06
     1)*X-1.03463270693454D-05 )*X+8.16324851790106D-05
     2)*X-5.55526624875562D-04 )*X+3.20512054753924D-03
     3)*X-1.51515139838540D-02 )*X+5.55555554649585D-02
     4)*X-1.42857142854412D-01 )*X+1.99999999999986D-01
  300 G=DEXP(-X)
      F1=((X+X)*F2+G)/3.0D+00
      WW1=(X+X)*F1+G
  301 T1=RT1/(RT1+1.0D+00)
      T2=RT2/(RT2+1.0D+00)
      T3=RT3/(RT3+1.0D+00)
      A2=F2-T1*F1
      A1=F1-T1*WW1
      WW3=(A2-T2*A1)/((T3-T2)*(T3-T1))
      WW2=(T3*A1-A2)/((T3-T2)*(T2-T1))
      WW1=WW1-WW2-WW3
      RETURN
   30 IF(X.GT.3.0D+00) GO TO 40
C     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3
      Y=X-2.0D+00
      IF(NROOTS.EQ.3) GO TO 33
      F1=        ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09
     1)*Y-2.14234468198419D-08 )*Y+2.17216556336318D-07
     2)*Y-1.98850171329371D-06 )*Y+1.62429321438911D-05
     3)*Y-1.16740298039895D-04 )*Y+7.24888732052332D-04
     4)*Y-3.79490003707156D-03 )*Y+1.61723488664661D-02
     5)*Y-5.29428148329736D-02 )*Y+1.15702180856167D-01
      WW1=(X+X)*F1+DEXP(-X)
      IF(NROOTS.EQ.2) GO TO 32
      RT1=F1/(WW1-F1)
      RETURN
   32 RT1=         (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11
     1)*Y-5.152207846962D-10)*Y-3.846389873308D-10)*Y+8.472253388380D-08
     2)*Y-1.85306035634293D-06 )*Y+2.47191693238413D-05
     3)*Y-2.49018321709815D-04 )*Y+2.19173220020161D-03
     4)*Y-1.63329339286794D-02 )*Y+8.68085688285261D-02
      RT2=         ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09
     1)*Y-1.878920917404D-08)*Y-1.725838516261D-07)*Y+2.247389642339D-06
     2)*Y+9.76783813082564D-06 )*Y-1.93160765581969D-04
     3)*Y-1.58064140671893D-03 )*Y+4.85928174507904D-02
     4)*Y-4.30761584997596D-01 )*Y+1.80400974537950D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   33 RT1=          (((((((( 1.44687969563318D-12*Y+4.85300143926755D-12
     1)*Y-6.55098264095516D-10 )*Y+1.56592951656828D-08
     2)*Y-2.60122498274734D-07 )*Y+3.86118485517386D-06
     3)*Y-5.13430986707889D-05 )*Y+6.03194524398109D-04
     4)*Y-6.11219349825090D-03 )*Y+4.52578254679079D-02
      RT2=           ((((((( 6.95964248788138D-10*Y-5.35281831445517D-09
     1)*Y-6.745205954533D-08)*Y+1.502366784525D-06)*Y+9.923326947376D-07
     2)*Y-3.89147469249594D-04 )*Y+7.51549330892401D-03
     3)*Y-8.48778120363400D-02 )*Y+5.73928229597613D-01
      RT3=          ((((((((-2.81496588401439D-10*Y+3.61058041895031D-09
     1)*Y+4.53631789436255D-08 )*Y-1.40971837780847D-07
     2)*Y-6.05865557561067D-06 )*Y-5.15964042227127D-05
     3)*Y+3.34761560498171D-05 )*Y+5.04871005319119D-02
     4)*Y-8.24708946991557D-01 )*Y+4.81234667357205D+00
      F2=        ((((((((((-1.48044231072140D-10*Y+1.78157031325097D-09
     1)*Y-1.92514145088973D-08 )*Y+1.92804632038796D-07
     2)*Y-1.73806555021045D-06 )*Y+1.39195169625425D-05
     3)*Y-9.74574633246452D-05 )*Y+5.83701488646511D-04
     4)*Y-2.89955494844975D-03 )*Y+1.13847001113810D-02
     5)*Y-3.23446977320647D-02 )*Y+5.29428148329709D-02
      GO TO 300
C     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3
   40 Y=X-4.0D+00
      IF(NROOTS.EQ.3) GO TO 43
      F1=        ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10
     1)*Y-3.614965656163D-09)*Y+3.760256799971D-08)*Y-3.553558319675D-07
     2)*Y+3.022556449731D-06)*Y-2.290098979647D-05)*Y+1.526537461148D-04
     3)*Y-8.81947375894379D-04 )*Y+4.33207949514611D-03
     4)*Y-1.75257821619926D-02 )*Y+5.28406320615584D-02
      WW1=(X+X)*F1+DEXP(-X)
      IF(NROOTS.EQ.2) GO TO 42
      RT1=F1/(WW1-F1)
      RETURN
   42 RT1=          ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11
     1)*Y-1.73508862390291D-09 )*Y+5.93066856324744D-08
     2)*Y-9.76085576741771D-07 )*Y+1.08484384385679D-05
     3)*Y-1.12608004981982D-04 )*Y+1.16210907653515D-03
     4)*Y-9.89572595720351D-03 )*Y+6.12589701086408D-02
      RT2=         (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10
     1)*Y+1.603498045240D-08)*Y-1.497986283037D-07)*Y-7.017002532106D-07
     2)*Y+1.85882653064034D-05 )*Y-2.04685420150802D-05
     3)*Y-2.49327728643089D-03 )*Y+3.56550690684281D-02
     4)*Y-2.60417417692375D-01 )*Y+1.12155283108289D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   43 RT1=           ((((((( 1.44265709189601D-11*Y-4.66622033006074D-10
     1)*Y+7.649155832025D-09)*Y-1.229940017368D-07)*Y+2.026002142457D-06
     2)*Y-2.87048671521677D-05 )*Y+3.70326938096287D-04
     3)*Y-4.21006346373634D-03 )*Y+3.50898470729044D-02
      RT2=          ((((((((-2.65526039155651D-11*Y+1.97549041402552D-10
     1)*Y+2.15971131403034D-09 )*Y-7.95045680685193D-08
     2)*Y+5.15021914287057D-07 )*Y+1.11788717230514D-05
     3)*Y-3.33739312603632D-04 )*Y+5.30601428208358D-03
     4)*Y-5.93483267268959D-02 )*Y+4.31180523260239D-01
      RT3=          ((((((((-3.92833750584041D-10*Y-4.16423229782280D-09
     1)*Y+4.42413039572867D-08 )*Y+6.40574545989551D-07
     2)*Y-3.05512456576552D-06 )*Y-1.05296443527943D-04
     3)*Y-6.14120969315617D-04 )*Y+4.89665802767005D-02
     4)*Y-6.24498381002855D-01 )*Y+3.36412312243724D+00
      F2=        ((((((((((-2.36788772599074D-11*Y+2.89147476459092D-10
     1)*Y-3.18111322308846D-09 )*Y+3.25336816562485D-08
     2)*Y-3.00873821471489D-07 )*Y+2.48749160874431D-06
     3)*Y-1.81353179793672D-05 )*Y+1.14504948737066D-04
     4)*Y-6.10614987696677D-04 )*Y+2.64584212770942D-03
     5)*Y-8.66415899015349D-03 )*Y+1.75257821619922D-02
      GO TO 300
   50 IF(X.GT.15.0D+00) GO TO 70
      G=DEXP(-X)
      IF(X.GT.10.0D+00) GO TO 60
C     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3
      WW1=    (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X
     1+5.3689283271887D-01)/X-3.2883030418398D-01)/X
     2+2.4645596956002D-01)/X-4.9984072848436D-01)/X
     3-3.1501078774085D-06)*G + DSQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 51,52,53
   51 RT1=F1/(WW1-F1)
      RETURN
   52 Y=X-7.5D+00
      RT1=     (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16
     1)*Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-7.719300212748D-13
     2)*Y+7.802544789997D-12)*Y+6.628721099436D-11)*Y-1.775564159743D-09
     3)*Y+1.713828823990D-08)*Y-1.497500187053D-07)*Y+2.283485114279D-06
     4)*Y-3.76953869614706D-05 )*Y+4.74791204651451D-04
     5)*Y-4.60448960876139D-03 )*Y+3.72458587837249D-02
      RT2=      (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13
     1)*Y-2.224334349799D-12)*Y+4.190559455515D-11)*Y-2.222722579924D-10
     2)*Y-2.624183464275D-09)*Y+6.128153450169D-08)*Y-4.383376014528D-07
     3)*Y-2.49952200232910D-06 )*Y+1.03236647888320D-04
     4)*Y-1.44614664924989D-03 )*Y+1.35094294917224D-02
     5)*Y-9.53478510453887D-02 )*Y+5.44765245686790D-01
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   53 F2=(F1+F1+F1-G)/(X+X)
      Y=X-7.5D+00
      RT1=       ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16
     1)*Y-6.736701449826D-14)*Y-6.264613873998D-13)*Y+1.315418927040D-11
     2)*Y-4.23879635610964D-11 )*Y+1.39032379769474D-09
     3)*Y-4.65449552856856D-08 )*Y+7.34609900170759D-07
     4)*Y-1.08656008854077D-05 )*Y+1.77930381549953D-04
     5)*Y-2.39864911618015D-03 )*Y+2.39112249488821D-02
      RT2=       ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15
     1)*Y-8.595618132088D-13)*Y-5.293620408757D-12)*Y-2.492175211635D-11
     2)*Y+2.73681574882729D-09 )*Y-1.06656985608482D-08
     3)*Y-4.40252529648056D-07 )*Y+9.68100917793911D-06
     4)*Y-1.68211091755327D-04 )*Y+2.69443611274173D-03
     5)*Y-3.23845035189063D-02 )*Y+2.75969447451882D-01
      RT3=      (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13
     1)*Y-1.985141104444D-12)*Y-2.309293727603D-11)*Y+3.917984522103D-10
     2)*Y+1.663165279876D-09)*Y-6.205591993923D-08)*Y+8.769581622041D-09
     3)*Y+8.97224398620038D-06 )*Y-3.14232666170796D-05
     4)*Y-1.83917335649633D-03 )*Y+3.51246831672571D-02
     5)*Y-3.22335051270860D-01 )*Y+1.73582831755430D+00
      GO TO 301
C     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3
   60 WW1=       (((-1.8784686463512D-01/X+2.2991849164985D-01)/X
     1-4.9893752514047D-01)/X-2.1916512131607D-05)*G + DSQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 61,62,63
   61 RT1=F1/(WW1-F1)
      RETURN
   62 RT1=      ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X
     1 -6.73760231824074D-02)*X+1.25705571069895D+00)*X
     2 +     (((-8.57609422987199D+03/X+5.91005939591842D+03)/X
     3 -1.70807677109425D+03)/X+2.64536689959503D+02)/X
     4 -2.38570496490846D+01)*G + R12/(X-R12)
      RT2=       ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X
     1 -4.22216483306320D+00)*X +       (((-2.08457050986847D+03/X
     2 -1.04999071905664D+03)/X+3.39891508992661D+02)/X
     3 -1.56184800325063D+02)/X+8.00839033297501D+00)*G + R22/(X-R22)
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   63 F2=(F1+F1+F1-G)/(X+X)
      Y=X-12.5D+00
      RT1=       ((((((((((( 4.42133001283090D-16*Y-2.77189767070441D-15
     1)*Y-4.084026087887D-14)*Y+5.379885121517D-13)*Y+1.882093066702D-12
     2)*Y-8.67286219861085D-11 )*Y+7.11372337079797D-10
     3)*Y-3.55578027040563D-09 )*Y+1.29454702851936D-07
     4)*Y-4.14222202791434D-06 )*Y+8.04427643593792D-05
     5)*Y-1.18587782909876D-03 )*Y+1.53435577063174D-02
      RT2=       ((((((((((( 6.85146742119357D-15*Y-1.08257654410279D-14
     1)*Y-8.579165965128D-13)*Y+6.642452485783D-12)*Y+4.798806828724D-11
     2)*Y-1.13413908163831D-09 )*Y+7.08558457182751D-09
     3)*Y-5.59678576054633D-08 )*Y+2.51020389884249D-06
     4)*Y-6.63678914608681D-05 )*Y+1.11888323089714D-03
     5)*Y-1.45361636398178D-02 )*Y+1.65077877454402D-01
      RT3=      (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14
     1)*Y-3.157134329361D-13)*Y+8.654129268056D-12)*Y-5.625235879301D-11
     2)*Y-7.718080513708D-10)*Y+2.064664199164D-08)*Y-1.567725007761D-07
     3)*Y-1.57938204115055D-06 )*Y+6.27436306915967D-05
     4)*Y-1.01308723606946D-03 )*Y+1.13901881430697D-02
     5)*Y-1.01449652899450D-01 )*Y+7.77203937334739D-01
      GO TO 301
   70 IF(X.GT.33.0D+00) GO TO 90
C     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3
      G=DEXP(-X)
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*G + DSQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 71,72,73
   71 RT1=F1/(WW1-F1)
      RETURN
   72 RT1=      ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X
     1 -1.71984023644904D-02)*X-1.37292644149838D-01)*X
     2 +       (-4.75742064274859D+01/X+9.21005186542857D+00)/X
     3 -2.31080873898939D-02)*G + R12/(X-R12)
      RT2=       ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X
     1 -4.02886174850252D+00)*X +         (-1.35831002139173D+02/X
     2 -8.66891724287962D+01)/X+2.98011277766958D+00)*G + R22/(X-R22)
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   73 F2=(F1+F1+F1-G)/(X+X)
      IF(X.GT.20.0D+00) GO TO 83
      RT1=    ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X
     1 -2.34112415981143D-02)*X+7.81425144913975D-01)*X
     2 -1.73209218219175D+01)*X+2.43517435690398D+02)*X
     3 +       (-1.97611541576986D+04/X+9.82441363463929D+03)/X
     4 -2.07970687843258D+03)*G + R13/(X-R13)
      RT2=     (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X
     1 -3.09337618731880D+00)*X+1.07037141010778D+02)*X
     2 -2.36659637247087D+03)*X +        ((-2.91669113681020D+06/X
     3 +1.41129505262758D+06)/X-2.91532335433779D+05)/X
     4 +3.35202872835409D+04)*G + R23/(X-R23)
      RT3=     ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X
     1 -7.83503697918455D-01)*X-1.84338896480695D+01)*X
     2 +4.04996712650414D+02)*X +         (-1.89829509315154D+05/X
     3 +5.11498390849158D+04)/X-6.88145821789955D+03)*G + R33/(X-R33)
      GO TO 301
   83 RT1=      ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X
     1 +1.31099142238996D+00)*X-1.88336409225481D+01)*X
     2 -6.60344754467191D+02 /X+1.64931462413877D+02)*G + R13/(X-R13)
      RT2=      ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X
     1 +1.13691058739678D+01)*X-1.65426392885291D+02)*X
     2 -6.30909125686731D+03 /X+1.52231757709236D+03)*G + R23/(X-R23)
      RT3=      ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X
     1 +1.73639054044562D+01)*X-3.57615122086961D+02)*X
     2 -1.45734701095912D+04 /X+2.69831813951849D+03)*G + R33/(X-R33)
      GO TO 301
C     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3
   90 WW1= DSQRT(PIE4/X)
      IF(NROOTS-2) 91,92,93
   91 RT1=0.5D+00/(X-0.5D+00)
      RETURN
   92 IF(X.GT.40.0D+00) GO TO 102
      G=DEXP(-X)
      RT1=(-8.78947307498880D-01*X+1.09243702330261D+01)*G + R12/(X-R12)
      RT2=(-9.28903924275977D+00*X+8.10642367843811D+01)*G + R22/(X-R22)
      WW2=( 4.46857389308400D+00*X-7.79250653461045D+01)*G + W22*WW1
      WW1=WW1-WW2
      RETURN
   93 IF(X.GT.47.0D+00) GO TO 103
      G=DEXP(-X)
      RT1=        ((-7.39058467995275D+00*X+3.21318352526305D+02)*X
     1 -3.99433696473658D+03)*G + R13/(X-R13)
      RT2=        ((-7.38726243906513D+01*X+3.13569966333873D+03)*X
     1 -3.86862867311321D+04)*G + R23/(X-R23)
      RT3=        ((-2.63750565461336D+02*X+1.04412168692352D+04)*X
     1 -1.28094577915394D+05)*G + R33/(X-R33)
      WW3=       ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X
     1 +1.92977367967984D+02)*X-1.67787926005344D+03)*G + W33*WW1
      WW2=        (( 6.15072615497811D+01*X-2.91980647450269D+03)*X
     1 +3.80794303087338D+04)*G + W23*WW1
      WW1=WW1-WW2-WW3
      RETURN
  102 RT1=R12/(X-R12)
      RT2=R22/(X-R22)
      WW2=W22*WW1
      WW1=WW1-WW2
      RETURN
  103 RT1=R13/(X-R13)
      RT2=R23/(X-R23)
      RT3=R33/(X-R33)
      WW2=W23*WW1
      WW3=W33*WW1
      WW1=WW1-WW2-WW3
      RETURN
      END

      SUBROUTINE HND_ROOT4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C          -----   VERSION FEBRUARY 16,1975   -----
      COMMON/HND_RYS/X,U(5),W(5),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R14,PIE4/1.45303521503316D-01, 7.85398163397448D-01/
      DATA R24,W24/ 1.33909728812636D+00, 2.34479815323517D-01/
      DATA R34,W34/ 3.92696350135829D+00, 1.92704402415764D-02/
      DATA R44,W44/ 8.58863568901199D+00, 2.25229076750736D-04/
      IF(X.GT.15.0D+00) GO TO 470
      IF(X.GT.5.0D+00) GO TO 450
      IF(X.GT.1.0D+00) GO TO 430
      IF(X.GT.3.0D-07) GO TO 420
C     X IS APPROXIMATELY ZERO.                   NROOTS = 4
      RT1=3.48198973061471D-02    -4.09645850660395D-03 *X
      RT2=3.81567185080042D-01    -4.48902570656719D-02 *X
      RT3=1.73730726945891D+00    -2.04389090547327D-01 *X
      RT4=1.18463056481549D+01    -1.39368301742312D+00 *X
      WW1=3.62683783378362D-01    -3.13844305713928D-02 *X
      WW2=3.13706645877886D-01    -8.98046242557724D-02 *X
      WW3=2.22381034453372D-01    -1.29314370958973D-01 *X
      WW4=1.01228536290376D-01    -8.28299075414321D-02 *X
      RETURN
C     X=0.0 TO 1.0                               NROOTS = 4
  420 RT1=            ((((((-1.95309614628539D-10*X+5.19765728707592D-09
     1)*X-1.01756452250573D-07 )*X+1.72365935872131D-06
     2)*X-2.61203523522184D-05 )*X+3.52921308769880D-04
     3)*X-4.09645850658433D-03 )*X+3.48198973061469D-02
      RT2=             (((((-1.89554881382342D-08*X+3.07583114342365D-07
     1)*X+1.270981734393D-06)*X-1.417298563884D-04)*X+3.226979163176D-03
     2)*X-4.48902570678178D-02 )*X+3.81567185080039D-01
      RT3=            (((((( 1.77280535300416D-09*X+3.36524958870615D-08
     1)*X-2.58341529013893D-07 )*X-1.13644895662320D-05
     2)*X-7.91549618884063D-05 )*X+1.03825827346828D-02
     3)*X-2.04389090525137D-01 )*X+1.73730726945889D+00
      RT4=             (((((-5.61188882415248D-08*X-2.49480733072460D-07
     1)*X+3.428685057114D-06)*X+1.679007454539D-04)*X+4.722855585715D-02
     2)*X-1.39368301737828D+00 )*X+1.18463056481543D+01
      WW1=            ((((((-1.14649303201279D-08*X+1.88015570196787D-07
     1)*X-2.33305875372323D-06 )*X+2.68880044371597D-05
     2)*X-2.94268428977387D-04 )*X+3.06548909776613D-03
     3)*X-3.13844305680096D-02 )*X+3.62683783378335D-01
      WW2=          ((((((((-4.11720483772634D-09*X+6.54963481852134D-08
     1)*X-7.20045285129626D-07 )*X+6.93779646721723D-06
     2)*X-6.05367572016373D-05 )*X+4.74241566251899D-04
     3)*X-3.26956188125316D-03 )*X+1.91883866626681D-02
     4)*X-8.98046242565811D-02 )*X+3.13706645877886D-01
      WW3=          ((((((((-3.41688436990215D-08*X+5.07238960340773D-07
     1)*X-5.01675628408220D-06 )*X+4.20363420922845D-05
     2)*X-3.08040221166823D-04 )*X+1.94431864731239D-03
     3)*X-1.02477820460278D-02 )*X+4.28670143840073D-02
     4)*X-1.29314370962569D-01 )*X+2.22381034453369D-01
      WW4=         ((((((((( 4.99660550769508D-09*X-7.94585963310120D-08
     1)*X+8.359072409485D-07)*X-7.422369210610D-06)*X+5.763374308160D-05
     2)*X-3.86645606718233D-04 )*X+2.18417516259781D-03
     3)*X-9.99791027771119D-03 )*X+3.48791097377370D-02
     4)*X-8.28299075413889D-02 )*X+1.01228536290376D-01
      RETURN
C     X= 1.0 TO 5.0                              NROOTS = 4
  430 Y=X-3.0D+00
      RT1=         (((((((((-1.48570633747284D-15*Y-1.33273068108777D-13
     1)*Y+4.068543696670D-12)*Y-9.163164161821D-11)*Y+2.046819017845D-09
     2)*Y-4.03076426299031D-08 )*Y+7.29407420660149D-07
     3)*Y-1.23118059980833D-05 )*Y+1.88796581246938D-04
     4)*Y-2.53262912046853D-03 )*Y+2.51198234505021D-02
      RT2=         ((((((((( 1.35830583483312D-13*Y-2.29772605964836D-12
     1)*Y-3.821500128045D-12)*Y+6.844424214735D-10)*Y-1.048063352259D-08
     2)*Y+1.50083186233363D-08 )*Y+3.48848942324454D-06
     3)*Y-1.08694174399193D-04 )*Y+2.08048885251999D-03
     4)*Y-2.91205805373793D-02 )*Y+2.72276489515713D-01
      RT3=         ((((((((( 5.02799392850289D-13*Y+1.07461812944084D-11
     1)*Y-1.482277886411D-10)*Y-2.153585661215D-09)*Y+3.654087802817D-08
     2)*Y+5.15929575830120D-07 )*Y-9.52388379435709D-06
     3)*Y-2.16552440036426D-04 )*Y+9.03551469568320D-03
     4)*Y-1.45505469175613D-01 )*Y+1.21449092319186D+00
      RT4=         (((((((((-1.08510370291979D-12*Y+6.41492397277798D-11
     1)*Y+7.542387436125D-10)*Y-2.213111836647D-09)*Y-1.448228963549D-07
     2)*Y-1.95670833237101D-06 )*Y-1.07481314670844D-05
     3)*Y+1.49335941252765D-04 )*Y+4.87791531990593D-02
     4)*Y-1.10559909038653D+00 )*Y+8.09502028611780D+00
      WW1=        ((((((((((-4.65801912689961D-14*Y+7.58669507106800D-13
     1)*Y-1.186387548048D-11)*Y+1.862334710665D-10)*Y-2.799399389539D-09
     2)*Y+4.148972684255D-08)*Y-5.933568079600D-07)*Y+8.168349266115D-06
     3)*Y-1.08989176177409D-04 )*Y+1.41357961729531D-03
     4)*Y-1.87588361833659D-02 )*Y+2.89898651436026D-01
      WW2=      ((((((((((((-1.46345073267549D-14*Y+2.25644205432182D-13
     1)*Y-3.116258693847D-12)*Y+4.321908756610D-11)*Y-5.673270062669D-10
     2)*Y+7.006295962960D-09)*Y-8.120186517000D-08)*Y+8.775294645770D-07
     3)*Y-8.77829235749024D-06 )*Y+8.04372147732379D-05
     4)*Y-6.64149238804153D-04 )*Y+4.81181506827225D-03
     5)*Y-2.88982669486183D-02 )*Y+1.56247249979288D-01
      WW3=     ((((((((((((( 9.06812118895365D-15*Y-1.40541322766087D-13
     1)*Y+1.919270015269D-12)*Y-2.605135739010D-11)*Y+3.299685839012D-10
     2)*Y-3.86354139348735D-09 )*Y+4.16265847927498D-08
     3)*Y-4.09462835471470D-07 )*Y+3.64018881086111D-06
     4)*Y-2.88665153269386D-05 )*Y+2.00515819789028D-04
     5)*Y-1.18791896897934D-03 )*Y+5.75223633388589D-03
     6)*Y-2.09400418772687D-02 )*Y+4.85368861938873D-02
      WW4=    ((((((((((((((-9.74835552342257D-16*Y+1.57857099317175D-14
     1)*Y-2.249993780112D-13)*Y+3.173422008953D-12)*Y-4.161159459680D-11
     2)*Y+5.021343560166D-10)*Y-5.545047534808D-09)*Y+5.554146993491D-08
     3)*Y-4.99048696190133D-07 )*Y+3.96650392371311D-06
     4)*Y-2.73816413291214D-05 )*Y+1.60106988333186D-04
     5)*Y-7.64560567879592D-04 )*Y+2.81330044426892D-03
     6)*Y-7.16227030134947D-03 )*Y+9.66077262223353D-03
      RETURN
  450 IF(X.GT.10.0D+00) GO TO 460
C     X=5.0 TO 10.0                              NROOTS = 4
      Y=X-7.5D+00
      RT1=         ((((((((( 4.64217329776215D-15*Y-6.27892383644164D-15
     1)*Y+3.462236347446D-13)*Y-2.927229355350D-11)*Y+5.090355371676D-10
     2)*Y-9.97272656345253D-09 )*Y+2.37835295639281D-07
     3)*Y-4.60301761310921D-06 )*Y+8.42824204233222D-05
     4)*Y-1.37983082233081D-03 )*Y+1.66630865869375D-02
      RT2=         ((((((((( 2.93981127919047D-14*Y+8.47635639065744D-13
     1)*Y-1.446314544774D-11)*Y-6.149155555753D-12)*Y+8.484275604612D-10
     2)*Y-6.10898827887652D-08 )*Y+2.39156093611106D-06
     3)*Y-5.35837089462592D-05 )*Y+1.00967602595557D-03
     4)*Y-1.57769317127372D-02 )*Y+1.74853819464285D-01
      RT3=        (((((((((( 2.93523563363000D-14*Y-6.40041776667020D-14
     1)*Y-2.695740446312D-12)*Y+1.027082960169D-10)*Y-5.822038656780D-10
     2)*Y-3.159991002539D-08)*Y+4.327249251331D-07)*Y+4.856768455119D-06
     3)*Y-2.54617989427762D-04 )*Y+5.54843378106589D-03
     4)*Y-7.95013029486684D-02 )*Y+7.20206142703162D-01
      RT4=       (((((((((((-1.62212382394553D-14*Y+7.68943641360593D-13
     1)*Y+5.764015756615D-12)*Y-1.380635298784D-10)*Y-1.476849808675D-09
     2)*Y+1.84347052385605D-08 )*Y+3.34382940759405D-07
     3)*Y-1.39428366421645D-06 )*Y-7.50249313713996D-05
     4)*Y-6.26495899187507D-04 )*Y+4.69716410901162D-02
     5)*Y-6.66871297428209D-01 )*Y+4.11207530217806D+00
      WW1=        ((((((((((-1.65995045235997D-15*Y+6.91838935879598D-14
     1)*Y-9.131223418888D-13)*Y+1.403341829454D-11)*Y-3.672235069444D-10
     2)*Y+6.366962546990D-09)*Y-1.039220021671D-07)*Y+1.959098751715D-06
     3)*Y-3.33474893152939D-05 )*Y+5.72164211151013D-04
     4)*Y-1.05583210553392D-02 )*Y+2.26696066029591D-01
      WW2=      ((((((((((((-3.57248951192047D-16*Y+6.25708409149331D-15
     1)*Y-9.657033089714D-14)*Y+1.507864898748D-12)*Y-2.332522256110D-11
     2)*Y+3.428545616603D-10)*Y-4.698730937661D-09)*Y+6.219977635130D-08
     3)*Y-7.83008889613661D-07 )*Y+9.08621687041567D-06
     4)*Y-9.86368311253873D-05 )*Y+9.69632496710088D-04
     5)*Y-8.14594214284187D-03 )*Y+8.50218447733457D-02
      WW3=     ((((((((((((( 1.64742458534277D-16*Y-2.68512265928410D-15
     1)*Y+3.788890667676D-14)*Y-5.508918529823D-13)*Y+7.555896810069D-12
     2)*Y-9.69039768312637D-11 )*Y+1.16034263529672D-09
     3)*Y-1.28771698573873D-08 )*Y+1.31949431805798D-07
     4)*Y-1.23673915616005D-06 )*Y+1.04189803544936D-05
     5)*Y-7.79566003744742D-05 )*Y+5.03162624754434D-04
     6)*Y-2.55138844587555D-03 )*Y+1.13250730954014D-02
      WW4=    ((((((((((((((-1.55714130075679D-17*Y+2.57193722698891D-16
     1)*Y-3.626606654097D-15)*Y+5.234734676175D-14)*Y-7.067105402134D-13
     2)*Y+8.793512664890D-12)*Y-1.006088923498D-10)*Y+1.050565098393D-09
     3)*Y-9.91517881772662D-09 )*Y+8.35835975882941D-08
     4)*Y-6.19785782240693D-07 )*Y+3.95841149373135D-06
     5)*Y-2.11366761402403D-05 )*Y+9.00474771229507D-05
     6)*Y-2.78777909813289D-04 )*Y+5.26543779837487D-04
      RETURN
C     X=10.0 TO 15.0                             NROOTS = 4
  460 Y=X-12.5D+00
      RT1=       ((((((((((( 4.94869622744119D-17*Y+8.03568805739160D-16
     1)*Y-5.599125915431D-15)*Y-1.378685560217D-13)*Y+7.006511663249D-13
     2)*Y+1.30391406991118D-11 )*Y+8.06987313467541D-11
     3)*Y-5.20644072732933D-09 )*Y+7.72794187755457D-08
     4)*Y-1.61512612564194D-06 )*Y+4.15083811185831D-05
     5)*Y-7.87855975560199D-04 )*Y+1.14189319050009D-02
      RT2=       ((((((((((( 4.89224285522336D-16*Y+1.06390248099712D-14
     1)*Y-5.446260182933D-14)*Y-1.613630106295D-12)*Y+3.910179118937D-12
     2)*Y+1.90712434258806D-10 )*Y+8.78470199094761D-10
     3)*Y-5.97332993206797D-08 )*Y+9.25750831481589D-07
     4)*Y-2.02362185197088D-05 )*Y+4.92341968336776D-04
     5)*Y-8.68438439874703D-03 )*Y+1.15825965127958D-01
      RT3=        (((((((((( 6.12419396208408D-14*Y+1.12328861406073D-13
     1)*Y-9.051094103059D-12)*Y-4.781797525341D-11)*Y+1.660828868694D-09
     2)*Y+4.499058798868D-10)*Y-2.519549641933D-07)*Y+4.977444040180D-06
     3)*Y-1.25858350034589D-04 )*Y+2.70279176970044D-03
     4)*Y-3.99327850801083D-02 )*Y+4.33467200855434D-01
      RT4=       ((((((((((( 4.63414725924048D-14*Y-4.72757262693062D-14
     1)*Y-1.001926833832D-11)*Y+6.074107718414D-11)*Y+1.576976911942D-09
     2)*Y-2.01186401974027D-08 )*Y-1.84530195217118D-07
     3)*Y+5.02333087806827D-06 )*Y+9.66961790843006D-06
     4)*Y-1.58522208889528D-03 )*Y+2.80539673938339D-02
     5)*Y-2.78953904330072D-01 )*Y+1.82835655238235D+00
      WW4=     ((((((((((((( 2.90401781000996D-18*Y-4.63389683098251D-17
     1)*Y+6.274018198326D-16)*Y-8.936002188168D-15)*Y+1.194719074934D-13
     2)*Y-1.45501321259466D-12 )*Y+1.64090830181013D-11
     3)*Y-1.71987745310181D-10 )*Y+1.63738403295718D-09
     4)*Y-1.39237504892842D-08 )*Y+1.06527318142151D-07
     5)*Y-7.27634957230524D-07 )*Y+4.12159381310339D-06
     6)*Y-1.74648169719173D-05 )*Y+8.50290130067818D-05
      WW3=      ((((((((((((-4.19569145459480D-17*Y+5.94344180261644D-16
     1)*Y-1.148797566469D-14)*Y+1.881303962576D-13)*Y-2.413554618391D-12
     2)*Y+3.372127423047D-11)*Y-4.933988617784D-10)*Y+6.116545396281D-09
     3)*Y-6.69965691739299D-08 )*Y+7.52380085447161D-07
     4)*Y-8.08708393262321D-06 )*Y+6.88603417296672D-05
     5)*Y-4.67067112993427D-04 )*Y+5.42313365864597D-03
      WW2=        ((((((((((-6.22272689880615D-15*Y+1.04126809657554D-13
     1)*Y-6.842418230913D-13)*Y+1.576841731919D-11)*Y-4.203948834175D-10
     2)*Y+6.287255934781D-09)*Y-8.307159819228D-08)*Y+1.356478091922D-06
     3)*Y-2.08065576105639D-05 )*Y+2.52396730332340D-04
     4)*Y-2.94484050194539D-03 )*Y+6.01396183129168D-02
      WW1=       (((-1.8784686463512D-01/X+2.2991849164985D-01)/X
     1-4.9893752514047D-01)/X-2.1916512131607D-05)*DEXP(-X)
     2 + DSQRT(PIE4/X)-WW4-WW3-WW2
      RETURN
  470 WW1= DSQRT(PIE4/X)
      IF(X.GT.35.0D+00) GO TO 490
      IF(X.GT.20.0D+00) GO TO 480
C     X=15.0 TO 20.0                             NROOTS = 4
      Y=X-17.5D+00
      RT1=       ((((((((((( 4.36701759531398D-17*Y-1.12860600219889D-16
     1)*Y-6.149849164164D-15)*Y+5.820231579541D-14)*Y+4.396602872143D-13
     2)*Y-1.24330365320172D-11 )*Y+6.71083474044549D-11
     3)*Y+2.43865205376067D-10 )*Y+1.67559587099969D-08
     4)*Y-9.32738632357572D-07 )*Y+2.39030487004977D-05
     5)*Y-4.68648206591515D-04 )*Y+8.34977776583956D-03
      RT2=       ((((((((((( 4.98913142288158D-16*Y-2.60732537093612D-16
     1)*Y-7.775156445127D-14)*Y+5.766105220086D-13)*Y+6.432696729600D-12
     2)*Y-1.39571683725792D-10 )*Y+5.95451479522191D-10
     3)*Y+2.42471442836205D-09 )*Y+2.47485710143120D-07
     4)*Y-1.14710398652091D-05 )*Y+2.71252453754519D-04
     5)*Y-4.96812745851408D-03 )*Y+8.26020602026780D-02
      RT3=       ((((((((((( 1.91498302509009D-15*Y+1.48840394311115D-14
     1)*Y-4.316925145767D-13)*Y+1.186495793471D-12)*Y+4.615806713055D-11
     2)*Y-5.54336148667141D-10 )*Y+3.48789978951367D-10
     3)*Y-2.79188977451042D-09 )*Y+2.09563208958551D-06
     4)*Y-6.76512715080324D-05 )*Y+1.32129867629062D-03
     5)*Y-2.05062147771513D-02 )*Y+2.88068671894324D-01
      RT4=       (((((((((((-5.43697691672942D-15*Y-1.12483395714468D-13
     1)*Y+2.826607936174D-12)*Y-1.266734493280D-11)*Y-4.258722866437D-10
     2)*Y+9.45486578503261D-09 )*Y-5.86635622821309D-08
     3)*Y-1.28835028104639D-06 )*Y+4.41413815691885D-05
     4)*Y-7.61738385590776D-04 )*Y+9.66090902985550D-03
     5)*Y-1.01410568057649D-01 )*Y+9.54714798156712D-01
      WW4=      ((((((((((((-7.56882223582704D-19*Y+7.53541779268175D-18
     1)*Y-1.157318032236D-16)*Y+2.411195002314D-15)*Y-3.601794386996D-14
     2)*Y+4.082150659615D-13)*Y-4.289542980767D-12)*Y+5.086829642731D-11
     3)*Y-6.35435561050807D-10 )*Y+6.82309323251123D-09
     4)*Y-5.63374555753167D-08 )*Y+3.57005361100431D-07
     5)*Y-2.40050045173721D-06 )*Y+4.94171300536397D-05
      WW3=       (((((((((((-5.54451040921657D-17*Y+2.68748367250999D-16
     1)*Y+1.349020069254D-14)*Y-2.507452792892D-13)*Y+1.944339743818D-12
     2)*Y-1.29816917658823D-11 )*Y+3.49977768819641D-10
     3)*Y-8.67270669346398D-09 )*Y+1.31381116840118D-07
     4)*Y-1.36790720600822D-06 )*Y+1.19210697673160D-05
     5)*Y-1.42181943986587D-04 )*Y+4.12615396191829D-03
      WW2=       (((((((((((-1.86506057729700D-16*Y+1.16661114435809D-15
     1)*Y+2.563712856363D-14)*Y-4.498350984631D-13)*Y+1.765194089338D-12
     2)*Y+9.04483676345625D-12 )*Y+4.98930345609785D-10
     3)*Y-2.11964170928181D-08 )*Y+3.98295476005614D-07
     4)*Y-5.49390160829409D-06 )*Y+7.74065155353262D-05
     5)*Y-1.48201933009105D-03 )*Y+4.97836392625268D-02
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*DEXP(-X)+WW1-WW2-WW3-WW4
      RETURN
C     X=20.0 TO 35.0                             NROOTS = 4
  480 G=DEXP(-X)
      RT1=    ((((((-4.45711399441838D-05*X+1.27267770241379D-03)*X
     1 -2.36954961381262D-01)*X+1.54330657903756D+01)*X
     2 -5.22799159267808D+02)*X+1.05951216669313D+04)*X
     3 +       (-2.51177235556236D+06/X+8.72975373557709D+05)/X
     4 -1.29194382386499D+05)*G + R14/(X-R14)
      RT2=     (((((-7.85617372254488D-02*X+6.35653573484868D+00)*X
     1 -3.38296938763990D+02)*X+1.25120495802096D+04)*X
     2 -3.16847570511637D+05)*X +        ((-1.02427466127427D+09/X
     3 +3.70104713293016D+08)/X-5.87119005093822D+07)/X
     4 +5.38614211391604D+06)*G + R24/(X-R24)
      RT3=     (((((-2.37900485051067D-01*X+1.84122184400896D+01)*X
     1 -1.00200731304146D+03)*X+3.75151841595736D+04)*X
     2 -9.50626663390130D+05)*X +        ((-2.88139014651985D+09/X
     3 +1.06625915044526D+09)/X-1.72465289687396D+08)/X
     4 +1.60419390230055D+07)*G + R34/(X-R34)
      RT4=    ((((((-6.00691586407385D-04*X-3.64479545338439D-01)*X
     1 +1.57496131755179D+01)*X-6.54944248734901D+02)*X
     2 +1.70830039597097D+04)*X-2.90517939780207D+05)*X
     3 +       (+3.49059698304732D+07/X-1.64944522586065D+07)/X
     4 +2.96817940164703D+06)*G + R44/(X-R44)
      IF(X.LE.25.0D+00)
     1WW4=   ((((((( 2.33766206773151D-07*X-3.81542906607063D-05)*X
     1 +3.51416601267000D-03)*X-1.66538571864728D-01)*X
     2 +4.80006136831847D+00)*X-8.73165934223603D+01)*X
     3 +9.77683627474638D+02)*X +           1.66000945117640D+04/X
     4 -6.14479071209961D+03)*G + W44*WW1
      IF(X.GT.25.0D+00)
     1WW4=    (((((( 5.74245945342286D-06*X-7.58735928102351D-05)*X
     1 +2.35072857922892D-04)*X-3.78812134013125D-03)*X
     2 +3.09871652785805D-01)*X-7.11108633061306D+00)*X
     3 +5.55297573149528D+01)*G + W44*WW1
      WW3=    (((((( 2.36392855180768D-04*X-9.16785337967013D-03)*X
     1 +4.62186525041313D-01)*X-1.96943786006540D+01)*X
     2 +4.99169195295559D+02)*X-6.21419845845090D+03)*X
     3 +      ((+5.21445053212414D+07/X-1.34113464389309D+07)/X
     4 +1.13673298305631D+06)/X-2.81501182042707D+03)*G + W34*WW1
      WW2=    (((((( 7.29841848989391D-04*X-3.53899555749875D-02)*X
     1 +2.07797425718513D+00)*X-1.00464709786287D+02)*X
     2 +3.15206108877819D+03)*X-6.27054715090012D+04)*X
     3 +       (+1.54721246264919D+07/X-5.26074391316381D+06)/X
     4 +7.67135400969617D+05)*G + W24*WW1
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*G + WW1-WW2-WW3-WW4
      RETURN
  490 IF(X.GT.53.0D+00) GO TO 495
C     X=35.0 TO 53.0                             NROOTS = 4
      G=DEXP(-X)*(X*X)**2
      RT4=        ((-2.19135070169653D-03*X-1.19108256987623D-01)*X
     1 -7.50238795695573D-01)*G + R44/(X-R44)
      RT3=        ((-9.65842534508637D-04*X-4.49822013469279D-02)*X
     1 +6.08784033347757D-01)*G + R34/(X-R34)
      RT2=        ((-3.62569791162153D-04*X-9.09231717268466D-03)*X
     1 +1.84336760556262D-01)*G + R24/(X-R24)
      RT1=        ((-4.07557525914600D-05*X-6.88846864931685D-04)*X
     1 +1.74725309199384D-02)*G + R14/(X-R14)
      WW4=        (( 5.76631982000990D-06*X-7.89187283804890D-05)*X
     1 +3.28297971853126D-04)*G + W44*WW1
      WW3=        (( 2.08294969857230D-04*X-3.77489954837361D-03)*X
     1 +2.09857151617436D-02)*G + W34*WW1
      WW2=        (( 6.16374517326469D-04*X-1.26711744680092D-02)*X
     1 +8.14504890732155D-02)*G + W24*WW1
      WW1=WW1-WW2-WW3-WW4
      RETURN
C     X=47.0 TO INFINITY                         NROOTS = 4
  495 RT1=R14/(X-R14)
      RT2=R24/(X-R24)
      RT3=R34/(X-R34)
      RT4=R44/(X-R44)
      WW4=W44*WW1
      WW3=W34*WW1
      WW2=W24*WW1
      WW1=WW1-WW2-WW3-WW4
      RETURN
      END

      SUBROUTINE HND_ROOT5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C          -----   VERSION  FEBRUARY 27,1975   -----
      COMMON/HND_RYS/X,U(5),W(5),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
      IF(X.GT.15.0D+00) GO TO 570
      IF(X.GT. 5.0D+00) GO TO 550
      IF(X.GT. 1.0D+00) GO TO 530
      IF(X.GT. 3.0D-07) GO TO 520
C     X IS APPROXIMATELY ZERO.                   NROOTS = 5
      RT1=2.26659266316985D-02    -2.15865967920897D-03 *X
      RT2=2.31271692140903D-01    -2.20258754389745D-02 *X
      RT3=8.57346024118836D-01    -8.16520023025515D-02 *X
      RT4=2.97353038120346D+00    -2.83193369647137D-01 *X
      RT5=1.84151859759051D+01    -1.75382723579439D+00 *X
      WW1=2.95524224714752D-01    -1.96867576909777D-02 *X
      WW2=2.69266719309995D-01    -5.61737590184721D-02 *X
      WW3=2.19086362515981D-01    -9.71152726793658D-02 *X
      WW4=1.49451349150580D-01    -1.02979262193565D-01 *X
      WW5=6.66713443086877D-02    -5.73782817488315D-02 *X
      RETURN
C     X=0.0 TO 1.0                               NROOTS = 5
  520 RT1=            ((((((-4.46679165328413D-11*X+1.21879111988031D-09
     1)*X-2.62975022612104D-08 )*X+5.15106194905897D-07
     2)*X-9.27933625824749D-06 )*X+1.51794097682482D-04
     3)*X-2.15865967920301D-03 )*X+2.26659266316985D-02
      RT2=            (((((( 1.93117331714174D-10*X-4.57267589660699D-09
     1)*X+2.48339908218932D-08 )*X+1.50716729438474D-06
     2)*X-6.07268757707381D-05 )*X+1.37506939145643D-03
     3)*X-2.20258754419939D-02 )*X+2.31271692140905D-01
      RT3=             ((((( 4.84989776180094D-09*X+1.31538893944284D-07
     1)*X-2.766753852879D-06)*X-7.651163510626D-05)*X+4.033058545972D-03
     2)*X-8.16520022916145D-02 )*X+8.57346024118779D-01
      RT4=              ((((-2.48581772214623D-07*X-4.34482635782585D-06
     1)*X-7.46018257987630D-07 )*X+1.01210776517279D-02
     2)*X-2.83193369640005D-01 )*X+2.97353038120345D+00
      RT5=             (((((-8.92432153868554D-09*X+1.77288899268988D-08
     1)*X+3.040754680666D-06)*X+1.058229325071D-04)*X+4.596379534985D-02
     2)*X-1.75382723579114D+00 )*X+1.84151859759049D+01
      WW1=            ((((((-2.03822632771791D-09*X+3.89110229133810D-08
     1)*X-5.84914787904823D-07 )*X+8.30316168666696D-06
     2)*X-1.13218402310546D-04 )*X+1.49128888586790D-03
     3)*X-1.96867576904816D-02 )*X+2.95524224714749D-01
      WW2=           ((((((( 8.62848118397570D-09*X-1.38975551148989D-07
     1)*X+1.602894068228D-06)*X-1.646364300836D-05)*X+1.538445806778D-04
     2)*X-1.28848868034502D-03 )*X+9.38866933338584D-03
     3)*X-5.61737590178812D-02 )*X+2.69266719309991D-01
      WW3=          ((((((((-9.41953204205665D-09*X+1.47452251067755D-07
     1)*X-1.57456991199322D-06 )*X+1.45098401798393D-05
     2)*X-1.18858834181513D-04 )*X+8.53697675984210D-04
     3)*X-5.22877807397165D-03 )*X+2.60854524809786D-02
     4)*X-9.71152726809059D-02 )*X+2.19086362515979D-01
      WW4=          ((((((((-3.84961617022042D-08*X+5.66595396544470D-07
     1)*X-5.52351805403748D-06 )*X+4.53160377546073D-05
     2)*X-3.22542784865557D-04 )*X+1.95682017370967D-03
     3)*X-9.77232537679229D-03 )*X+3.79455945268632D-02
     4)*X-1.02979262192227D-01 )*X+1.49451349150573D-01
      WW5=         ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08
     1)*X+6.743541482689D-07)*X-5.917993920224D-06)*X+4.531969237381D-05
     2)*X-2.99102856679638D-04 )*X+1.65695765202643D-03
     3)*X-7.40671222520653D-03 )*X+2.50889946832192D-02
     4)*X-5.73782817487958D-02 )*X+6.66713443086877D-02
      RETURN
C     X=1.0 TO 5.0                               NROOTS = 5
  530 Y=X-3.0D+00
      RT1=          ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13
     1)*Y-2.11414838976129D-11 )*Y+5.09822003260014D-10
     2)*Y-1.16002134438663D-08 )*Y+2.46810694414540D-07
     3)*Y-4.92556826124502D-06 )*Y+9.02580687971053D-05
     4)*Y-1.45190025120726D-03 )*Y+1.73416786387475D-02
      RT2=         ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14
     1)*Y-4.831059411392D-12)*Y+1.136643908832D-10)*Y-1.104373076913D-09
     2)*Y-2.35346740649916D-08 )*Y+1.43772622028764D-06
     3)*Y-4.23405023015273D-05 )*Y+9.12034574793379D-04
     4)*Y-1.52479441718739D-02 )*Y+1.76055265928744D-01
      RT3=         (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13
     1)*Y+1.847170956043D-11)*Y-3.390752744265D-10)*Y-2.995532064116D-09
     2)*Y+1.57456141058535D-07 )*Y-3.95859409711346D-07
     3)*Y-9.58924580919747D-05 )*Y+3.23551502557785D-03
     4)*Y-5.97587007636479D-02 )*Y+6.46432853383057D-01
      RT4=          ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11
     1)*Y+8.83758848468769D-10 )*Y+1.59166632851267D-08
     2)*Y-1.32581997983422D-07 )*Y-7.60223407443995D-06
     3)*Y-7.41019244900952D-05 )*Y+9.81432631743423D-03
     4)*Y-2.23055570487771D-01 )*Y+2.21460798080643D+00
      RT5=         ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12
     1)*Y-8.776668218053D-11)*Y-2.342817613343D-09)*Y-3.496962018025D-08
     2)*Y-3.03172870136802D-07 )*Y+1.50511293969805D-06
     3)*Y+1.37704919387696D-04 )*Y+4.70723869619745D-02
     4)*Y-1.47486623003693D+00 )*Y+1.35704792175847D+01
      WW1=         ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12
     1)*Y+3.485512360993D-11)*Y-6.277497362235D-10)*Y+1.100758247388D-08
     2)*Y-1.88329804969573D-07 )*Y+3.12338120839468D-06
     3)*Y-5.04404167403568D-05 )*Y+8.00338056610995D-04
     4)*Y-1.30892406559521D-02 )*Y+2.47383140241103D-01
      WW2=       ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13
     1)*Y+7.743219385056D-12)*Y-1.146022750992D-10)*Y+1.615238462197D-09
     2)*Y-2.15479017572233D-08 )*Y+2.70933462557631D-07
     3)*Y-3.18750295288531D-06 )*Y+3.47425221210099D-05
     4)*Y-3.45558237388223D-04 )*Y+3.05779768191621D-03
     5)*Y-2.29118251223003D-02 )*Y+1.59834227924213D-01
      WW3=      ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13
     1)*Y-7.184330797139D-12)*Y+9.763932908544D-11)*Y-1.244014559219D-09
     2)*Y+1.472744068942D-08)*Y-1.611749975234D-07)*Y+1.616487851917D-06
     3)*Y-1.46852359124154D-05 )*Y+1.18900349101069D-04
     4)*Y-8.37562373221756D-04 )*Y+4.93752683045845D-03
     5)*Y-2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4=     ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13
     1)*Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+3.679254029085D-10
     2)*Y-4.23775673047899D-09 )*Y+4.46559231067006D-08
     3)*Y-4.26488836563267D-07 )*Y+3.64721335274973D-06
     4)*Y-2.74868382777722D-05 )*Y+1.78586118867488D-04
     5)*Y-9.68428981886534D-04 )*Y+4.16002324339929D-03
     6)*Y-1.28290192663141D-02 )*Y+2.22353727685016D-02
      WW5=    ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14
     1)*Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y-3.372639523006D-11
     2)*Y+4.025371849467D-10)*Y-4.389453269417D-09)*Y+4.332753856271D-08
     3)*Y-3.82673275931962D-07 )*Y+2.98006900751543D-06
     4)*Y-2.00718990300052D-05 )*Y+1.13876001386361D-04
     5)*Y-5.23627942443563D-04 )*Y+1.83524565118203D-03
     6)*Y-4.37785737450783D-03 )*Y+5.36963805223095D-03
      RETURN
  550 IF(X.GT.10.0D+00) GO TO 560
C     X=5.0 TO 10.0                              NROOTS = 5
      Y=X-7.5D+00
      RT1=          ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13
     1)*Y-4.81561201185876D-12 )*Y+1.56666512163407D-10
     2)*Y-3.73782213255083D-09 )*Y+9.15858355075147D-08
     3)*Y-2.13775073585629D-06 )*Y+4.56547356365536D-05
     4)*Y-8.68003909323740D-04 )*Y+1.22703754069176D-02
      RT2=         (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14
     1)*Y-1.296476623788D-12)*Y+1.477175434354D-11)*Y+5.464102147892D-10
     2)*Y-2.42538340602723D-08 )*Y+8.20460740637617D-07
     3)*Y-2.20379304598661D-05 )*Y+4.90295372978785D-04
     4)*Y-9.14294111576119D-03 )*Y+1.22590403403690D-01
      RT3=         ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13
     1)*Y+1.176946020731D-12)*Y+1.725627235645D-10)*Y-3.686383856300D-09
     2)*Y+2.87495324207095D-08 )*Y+1.71307311000282D-06
     3)*Y-7.94273603184629D-05 )*Y+2.00938064965897D-03
     4)*Y-3.63329491677178D-02 )*Y+4.34393683888443D-01
      RT4=        ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14
     1)*Y+3.753542914426D-12)*Y-2.708018219579D-11)*Y-1.190574776587D-09
     2)*Y+1.106696436509D-08)*Y+3.954955671326D-07)*Y-4.398596059588D-06
     3)*Y-2.01087998907735D-04 )*Y+7.89092425542937D-03
     4)*Y-1.42056749162695D-01 )*Y+1.39964149420683D+00
      RT5=        ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12
     1)*Y+6.861649627426D-12)*Y+6.082671496226D-10)*Y+5.381160105420D-09
     2)*Y-6.253297138700D-08)*Y-2.135966835050D-06)*Y-2.373394341886D-05
     3)*Y+2.88711171412814D-06 )*Y+4.85221195290753D-02
     4)*Y-1.04346091985269D+00 )*Y+7.89901551676692D+00
      WW1=         ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13
     1)*Y+4.761246208720D-12)*Y-9.535763686605D-11)*Y+2.225273630974D-09
     2)*Y-4.49796778054865D-08 )*Y+9.17812870287386D-07
     3)*Y-1.86764236490502D-05 )*Y+3.76807779068053D-04
     4)*Y-8.10456360143408D-03 )*Y+2.01097936411496D-01
      WW2=       ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14
     1)*Y+3.973252415832D-13)*Y-6.830539401049D-12)*Y+1.140771033372D-10
     2)*Y-1.82546185762009D-09 )*Y+2.77209637550134D-08
     3)*Y-4.01726946190383D-07 )*Y+5.48227244014763D-06
     4)*Y-6.95676245982121D-05 )*Y+8.05193921815776D-04
     5)*Y-8.15528438784469D-03 )*Y+9.71769901268114D-02
      WW3=      ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14
     1)*Y-2.022863065220D-13)*Y+3.058055403795D-12)*Y-4.387890955243D-11
     2)*Y+5.923946274445D-10)*Y-7.503659964159D-09)*Y+8.851599803902D-08
     3)*Y-9.65561998415038D-07 )*Y+9.60884622778092D-06
     4)*Y-8.56551787594404D-05 )*Y+6.66057194311179D-04
     5)*Y-4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4=    ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16
     1)*Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y-5.881379088074D-13
     2)*Y+7.805245193391D-12)*Y-9.632707991704D-11)*Y+1.099047050624D-09
     3)*Y-1.15042731790748D-08 )*Y+1.09415155268932D-07
     4)*Y-9.33687124875935D-07 )*Y+7.02338477986218D-06
     5)*Y-4.53759748787756D-05 )*Y+2.41722511389146D-04
     6)*Y-9.75935943447037D-04 )*Y+2.57520532789644D-03
      WW5=   ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17
     1)*Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+4.114588668138D-14
     2)*Y-5.44436631413933D-13 )*Y+6.64976446790959D-12
     3)*Y-7.44560069974940D-11 )*Y+7.57553198166848D-10
     4)*Y-6.92956101109829D-09 )*Y+5.62222859033624D-08
     5)*Y-3.97500114084351D-07 )*Y+2.39039126138140D-06
     6)*Y-1.18023950002105D-05 )*Y+4.52254031046244D-05
     7)*Y-1.21113782150370D-04 )*Y+1.75013126731224D-04
      RETURN
C     X=10.0 TO 15.0                             NROOTS = 5
  560 Y=X-12.5D+00
      RT1=        ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16
     1)*Y+1.395993802064D-14)*Y+3.660484641252D-14)*Y-4.154857548139D-12
     2)*Y+2.301379846544D-11)*Y-1.033307012866D-09)*Y+3.997777641049D-08
     3)*Y-9.35118186333939D-07 )*Y+2.38589932752937D-05
     4)*Y-5.35185183652937D-04 )*Y+8.85218988709735D-03
      RT2=        ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15
     1)*Y+1.737896339191D-13)*Y+8.964205979517D-14)*Y-3.538906780633D-11
     2)*Y+9.561341254948D-11)*Y-9.772831891310D-09)*Y+4.240340194620D-07
     3)*Y-1.02384302866534D-05 )*Y+2.57987709704822D-04
     4)*Y-5.54735977651677D-03 )*Y+8.68245143991948D-02
      RT3=        ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14
     1)*Y+7.884307667104D-13)*Y-9.023398159510D-13)*Y-5.814101544957D-11
     2)*Y-1.333480437968D-09)*Y-2.217064940373D-08)*Y+1.643290788086D-06
     3)*Y-4.39602147345028D-05 )*Y+1.08648982748911D-03
     4)*Y-2.13014521653498D-02 )*Y+2.94150684465425D-01
      RT4=        ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15
     1)*Y+8.960828117859D-13)*Y+5.214153461337D-11)*Y-1.106601744067D-10
     2)*Y-2.007890743962D-08)*Y+1.543764346501D-07)*Y+4.520749076914D-06
     3)*Y-1.88893338587047D-04 )*Y+4.73264487389288D-03
     4)*Y-7.91197893350253D-02 )*Y+8.60057928514554D-01
      RT5=       (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14
     1)*Y+5.587369053655D-12)*Y-3.045253104617D-12)*Y-1.223983883080D-09
     2)*Y-2.05603889396319D-09 )*Y+2.58604071603561D-07
     3)*Y+1.34240904266268D-06 )*Y-5.72877569731162D-05
     4)*Y-9.56275105032191D-04 )*Y+4.23367010370921D-02
     5)*Y-5.76800927133412D-01 )*Y+3.87328263873381D+00
      WW1=         ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14
     1)*Y+5.851494250405D-14)*Y-4.234204823846D-11)*Y+3.911507312679D-10
     2)*Y-9.65094802088511D-09 )*Y+3.42197444235714D-07
     3)*Y-7.51821178144509D-06 )*Y+1.94218051498662D-04
     4)*Y-5.38533819142287D-03 )*Y+1.68122596736809D-01
      WW2=        ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14
     1)*Y-5.500330153548D-13)*Y+1.003849567976D-11)*Y-1.720997242621D-10
     2)*Y+3.533277061402D-09)*Y-6.389171736029D-08)*Y+1.046236652393D-06
     3)*Y-1.73148206795827D-05 )*Y+2.57820531617185D-04
     4)*Y-3.46188265338350D-03 )*Y+7.03302497508176D-02
      WW3=       ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15
     1)*Y+9.945311467434D-14)*Y-1.749051512721D-12)*Y+2.768503957853D-11
     2)*Y-4.08688551136506D-10 )*Y+6.04189063303610D-09
     3)*Y-8.23540111024147D-08 )*Y+1.01503783870262D-06
     4)*Y-1.20490761741576D-05 )*Y+1.26928442448148D-04
     5)*Y-1.05539461930597D-03 )*Y+1.15543698537013D-02
      WW4=     ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17
     1)*Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+1.491302084832D-13
     2)*Y-2.06638666222265D-12 )*Y+2.67958697789258D-11
     3)*Y-3.23322654638336D-10 )*Y+3.63722952167779D-09
     4)*Y-3.75484943783021D-08 )*Y+3.49164261987184D-07
     5)*Y-2.92658670674908D-06 )*Y+2.12937256719543D-05
     6)*Y-1.19434130620929D-04 )*Y+6.45524336158384D-04
      WW5=    ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18
     1)*Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y-6.301348858104D-15
     2)*Y+8.031304476153D-14)*Y-9.446196472547D-13)*Y+1.018245804339D-11
     3)*Y-9.96995451348129D-11 )*Y+8.77489010276305D-10
     4)*Y-6.84655877575364D-09 )*Y+4.64460857084983D-08
     5)*Y-2.66924538268397D-07 )*Y+1.24621276265907D-06
     6)*Y-4.30868944351523D-06 )*Y+9.94307982432868D-06
      RETURN
  570 IF(X.GT.25.0D+00) GO TO 590
      IF(X.GT.20.0D+00) GO TO 580
C     X=15.0 TO 20.0                             NROOTS = 5
      Y=X-17.5D+00
      RT1=        (((((((((( 1.91875764545740D-16*Y+7.8357401095707 D-16
     1)*Y-3.260875931644D-14)*Y-1.186752035569D-13)*Y+4.275180095653D-12
     2)*Y+3.357056136731D-11)*Y-1.123776903884D-09)*Y+1.231203269887D-08
     3)*Y-3.99851421361031D-07 )*Y+1.45418822817771D-05
     4)*Y-3.49912254976317D-04 )*Y+6.67768703938812D-03
      RT2=        (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14
     1)*Y-3.385363492036D-13)*Y-1.615655871159D-12)*Y+4.527419140333D-11
     2)*Y+3.853670706486D-10)*Y-1.184607130107D-08)*Y+1.347873288827D-07
     3)*Y-4.47788241748377D-06 )*Y+1.54942754358273D-04
     4)*Y-3.55524254280266D-03 )*Y+6.44912219301603D-02
      RT3=        (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14
     1)*Y-1.249779730869D-12)*Y-1.020720636353D-11)*Y+1.814709816693D-10
     2)*Y+1.766397336977D-09)*Y-4.603559449010D-08)*Y+5.863956443581D-07
     3)*Y-2.03797212506691D-05 )*Y+6.31405161185185D-04
     4)*Y-1.30102750145071D-02 )*Y+2.10244289044705D-01
      RT4=       (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14
     1)*Y+4.859447665850D-13)*Y-3.217227223463D-12)*Y-7.484522135512D-11
     2)*Y+7.19101516047753D-10 )*Y+6.88409355245582D-09
     3)*Y-1.44374545515769D-07 )*Y+2.74941013315834D-06
     4)*Y-1.02790452049013D-04 )*Y+2.59924221372643D-03
     5)*Y-4.35712368303551D-02 )*Y+5.62170709585029D-01
      RT5=       ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13
     1)*Y-3.892741622280D-12)*Y-7.755793199043D-12)*Y+9.492190032313D-10
     2)*Y-4.98680128123353D-09 )*Y-1.81502268782664D-07
     3)*Y+2.69463269394888D-06 )*Y+2.50032154421640D-05
     4)*Y-1.33684303917681D-03 )*Y+2.29121951862538D-02
     5)*Y-2.45653725061323D-01 )*Y+1.89999883453047D+00
      WW1=        (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16
     1)*Y-3.000659497257D-13)*Y+2.021279817961D-13)*Y+3.853596935400D-11
     2)*Y+1.461418533652D-10)*Y-1.014517563435D-08)*Y+1.132736008979D-07
     3)*Y-2.86605475073259D-06 )*Y+1.21958354908768D-04
     4)*Y-3.86293751153466D-03 )*Y+1.45298342081522D-01
      WW2=        ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15
     1)*Y+1.220613939709D-13)*Y+1.275068098526D-12)*Y-5.341838883262D-11
     2)*Y+6.161037256669D-10)*Y-1.009147879750D-08)*Y+2.907862965346D-07
     3)*Y-6.12300038720919D-06 )*Y+1.00104454489518D-04
     4)*Y-1.80677298502757D-03 )*Y+5.78009914536630D-02
      WW3=        ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15
     1)*Y+2.606163540537D-15)*Y+1.983799950150D-12)*Y-5.400548574357D-11
     2)*Y+6.638043374114D-10)*Y-8.799518866802D-09)*Y+1.791418482685D-07
     3)*Y-2.96075397351101D-06 )*Y+3.38028206156144D-05
     4)*Y-3.58426847857878D-04 )*Y+8.39213709428516D-03
      WW4=       ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16
     1)*Y+4.745009557656D-15)*Y-6.033814209875D-14)*Y+1.049256040808D-12
     2)*Y-1.70859789556117D-11 )*Y+2.15219425727959D-10
     3)*Y-2.52746574206884D-09 )*Y+3.27761714422960D-08
     4)*Y-3.90387662925193D-07 )*Y+3.46340204593870D-06
     5)*Y-2.43236345136782D-05 )*Y+3.54846978585226D-04
      WW5=     ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19
     1)*Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+1.246599177672D-15
     2)*Y-1.56872999797549D-14 )*Y+1.87305099552692D-13
     3)*Y-2.09498886675861D-12 )*Y+2.11630022068394D-11
     4)*Y-1.92566242323525D-10 )*Y+1.62012436344069D-09
     5)*Y-1.23621614171556D-08 )*Y+7.72165684563049D-08
     6)*Y-3.59858901591047D-07 )*Y+2.43682618601000D-06
      RETURN
C     X=20.0 TO 25.0                             NROOTS = 5
  580 Y=X-22.5D+00
      RT1=         (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15
     1)*Y+1.445982921243D-13)*Y-2.676703245252D-12)*Y+5.823521627177D-12
     2)*Y+2.17264723874381D-10 )*Y+3.56242145897468D-09
     3)*Y-3.03763737404491D-07 )*Y+9.46859114120901D-06
     4)*Y-2.30896753853196D-04 )*Y+5.24663913001114D-03
      RT2=        (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14
     1)*Y+6.184065097200D-14)*Y+1.649846591230D-12)*Y-2.729713905266D-11
     2)*Y+3.709913790650D-11)*Y+2.216486288382D-09)*Y+4.616160236414D-08
     3)*Y-3.32380270861364D-06 )*Y+9.84635072633776D-05
     4)*Y-2.30092118015697D-03 )*Y+5.00845183695073D-02
      RT3=        (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14
     1)*Y+1.136466605916D-13)*Y+7.546203883874D-12)*Y-9.635646767455D-11
     2)*Y-8.295965491209D-11)*Y+7.534109114453D-09)*Y+2.699970652707D-07
     3)*Y-1.42982334217081D-05 )*Y+3.78290946669264D-04
     4)*Y-8.03133015084373D-03 )*Y+1.58689469640791D-01
      RT4=        (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13
     1)*Y-7.522712577474D-13)*Y+3.209520801187D-11)*Y-2.075594313618D-10
     2)*Y-2.070575894402D-09)*Y+7.323046997451D-09)*Y+1.851491550417D-06
     3)*Y-6.37524802411383D-05 )*Y+1.36795464918785D-03
     4)*Y-2.42051126993146D-02 )*Y+3.97847167557815D-01
      RT5=        ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12
     1)*Y-4.286617818951D-13)*Y-2.632066100073D-10)*Y+4.804518986559D-09
     2)*Y-1.835675889421D-08)*Y-1.068175391334D-06)*Y+3.292234974141D-05
     3)*Y-5.94805357558251D-04 )*Y+8.29382168612791D-03
     4)*Y-9.93122509049447D-02 )*Y+1.09857804755042D+00
      WW1=         (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13
     1)*Y+7.817349237071D-13)*Y-2.547619474232D-11)*Y+1.479321506529D-10
     2)*Y+1.52314028857627D-09 )*Y+9.20072040917242D-09
     3)*Y-2.19427111221848D-06 )*Y+8.65797782880311D-05
     4)*Y-2.82718629312875D-03 )*Y+1.28718310443295D-01
      WW2=         ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14
     1)*Y-2.358734508092D-13)*Y+8.261326648131D-12)*Y+9.229645304956D-11
     2)*Y-5.68108973828949D-09 )*Y+1.22477891136278D-07
     3)*Y-2.11919643127927D-06 )*Y+4.23605032368922D-05
     4)*Y-1.14423444576221D-03 )*Y+5.06607252890186D-02
      WW3=         ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14
     1)*Y-4.157593182747D-14)*Y+4.214670817758D-12)*Y+6.705582751532D-11
     2)*Y-3.36086411698418D-09 )*Y+6.07453633298986D-08
     3)*Y-7.40736211041247D-07 )*Y+8.84176371665149D-06
     4)*Y-1.72559275066834D-04 )*Y+7.16639814253567D-03
      WW4=       (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18
     1)*Y+6.126212599772D-16)*Y-8.526651626939D-15)*Y+4.826636065733D-14
     2)*Y-3.39554163649740D-13 )*Y+1.67070784862985D-11
     3)*Y-4.42671979311163D-10 )*Y+6.77368055908400D-09
     4)*Y-7.03520999708859D-08 )*Y+6.04993294708874D-07
     5)*Y-7.80555094280483D-06 )*Y+2.85954806605017D-04
      WW5=      ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20
     1)*Y-1.586937691507D-18)*Y+3.357639744582D-17)*Y-4.810285046442D-16
     2)*Y+5.386312669975D-15)*Y-6.117895297439D-14)*Y+8.441808227634D-13
     3)*Y-1.18527596836592D-11 )*Y+1.36296870441445D-10
     4)*Y-1.17842611094141D-09 )*Y+7.80430641995926D-09
     5)*Y-5.97767417400540D-08 )*Y+1.65186146094969D-06
      RETURN
  590 WW1= DSQRT(PIE4/X)
      IF(X.GT.40.0D+00) GO TO 595
C     X=25.0 TO 40.0                             NROOTS = 5
      G=DEXP(-X)
      RT1=  ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X
     1 -1.59437614121125D-02)*X+1.13467897349442D+00)*X
     2 -4.47216460864586D+01)*X+1.06251216612604D+03)*X
     3 -1.52073917378512D+04)*X+1.20662887111273D+05)*X
     4 -4.07186366852475D+05)*G + R15/(X-R15)
      RT2=  ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X
     1 -1.50043662589017D-01)*X+1.05563640866077D+01)*X
     2 -4.10468817024806D+02)*X+9.62604416506819D+03)*X
     3 -1.35888069838270D+05)*X+1.06107577038340D+06)*X
     4 -3.51190792816119D+06)*G + R25/(X-R25)
      RT3=  ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X
     1 -4.01048115525954D-01)*X+2.78360021977405D+01)*X
     2 -1.04891729356965D+03)*X+2.36985942687423D+04)*X
     3 -3.19504627257548D+05)*X+2.34879693563358D+06)*X
     4 -7.16341568174085D+06)*G + R35/(X-R35)
      RT4=  ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X
     1 -7.65735935499627D-02)*X+9.12692349152792D+00)*X
     2 -2.32077034386717D+02)*X+2.81839578728845D+02)*X
     3 +9.59529683876419D+04)*X-1.77638956809518D+06)*X
     4 +1.02489759645410D+07)*G + R45/(X-R45)
      RT5=  ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X
     1 +1.12594870794668D+00)*X-4.56752462103909D+01)*X
     2 +1.05804526830637D+03)*X-1.16003199605875D+04)*X
     3 -4.07297627297272D+04)*X+2.22215528319857D+06)*X
     4 -1.61196455032613D+07)*G + R55/(X-R55)
      WW5= (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X
     1 -1.63960915431080D-05)*X+1.15791154612838D-03)*X
     2 -5.30573476742071D-02)*X+1.61156533367153D+00)*X
     3 -3.23248143316007D+01)*X+4.12007318109157D+02)*X
     4 -3.02260070158372D+03)*X+9.71575094154768D+03)*G + W55*WW1
      WW4= (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X
     1 -9.04491430884113D-04)*X+6.37686375770059D-02)*X
     2 -2.96135703135647D+00)*X+9.15142356996330D+01)*X
     3 -1.86971865249111D+03)*X+2.42945528916947D+04)*X
     4 -1.81852473229081D+05)*X+5.96854758661427D+05)*G + W45*WW1
      WW3=  (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X
     1 +1.18520453711586D-01)*X-6.69649981309161D+00)*X
     2 +2.44789386487321D+02)*X-5.68832664556359D+03)*X
     3 +8.14507604229357D+04)*X-6.55181056671474D+05)*X
     4 +2.26410896607237D+06)*G + W35*WW1
      WW2=  (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X
     1 +1.61077633475573D-01)*X-8.96743743396132D+00)*X
     2 +3.28062687293374D+02)*X-7.65722701219557D+03)*X
     3 +1.10255055017664D+05)*X-8.92528122219324D+05)*X
     4 +3.10638627744347D+06)*G + W25*WW1
      WW1=WW1-0.01962D+00*G -WW2-WW3-WW4-WW5
      RETURN
  595 IF(X.GT.59.0D+00) GO TO 599
C     X=40.0 TO 59.0                             NROOTS = 5
      XXX=X**3
      G=XXX*DEXP(-X)
      RT1=       (((-2.43758528330205D-02*X+2.07301567989771D+00)*X
     1 -6.45964225381113D+01)*X+7.14160088655470D+02)*G + R15/(X-R15)
      RT2=       (((-2.28861955413636D-01*X+1.93190784733691D+01)*X
     1 -5.99774730340912D+02)*X+6.61844165304871D+03)*G + R25/(X-R25)
      RT3=       (((-6.95053039285586D-01*X+5.76874090316016D+01)*X
     1 -1.77704143225520D+03)*X+1.95366082947811D+04)*G + R35/(X-R35)
      RT4=       (((-1.58072809087018D+00*X+1.27050801091948D+02)*X
     1 -3.86687350914280D+03)*X+4.23024828121420D+04)*G + R45/(X-R45)
      RT5=       (((-3.33963830405396D+00*X+2.51830424600204D+02)*X
     1 -7.57728527654961D+03)*X+8.21966816595690D+04)*G + R55/(X-R55)
      G=XXX*G
      WW5=        (( 1.35482430510942D-08*X-3.27722199212781D-07)*X
     1 +2.41522703684296D-06)*G + W55*WW1
      WW4=        (( 1.23464092261605D-06*X-3.55224564275590D-05)*X
     1 +3.03274662192286D-04)*G + W45*WW1
      WW3=        (( 1.34547929260279D-05*X-4.19389884772726D-04)*X
     1 +3.87706687610809D-03)*G + W35*WW1
      WW2=        (( 2.09539509123135D-05*X-6.87646614786982D-04)*X
     1 +6.68743788585688D-03)*G + W25*WW1
      WW1=WW1-WW2-WW3-WW4-WW5
      RETURN
C     X=59.0 TO INFINITY                         NROOTS = 5
  599 RT1=R15/(X-R15)
      RT2=R25/(X-R25)
      RT3=R35/(X-R35)
      RT4=R45/(X-R45)
      RT5=R55/(X-R55)
      WW2=W25*WW1
      WW3=W35*WW1
      WW4=W45*WW1
      WW5=W55*WW1
      WW1=WW1-WW2-WW3-WW4-WW5
      RETURN
      END

C
C     ------ SUBROUTINES TAKEN FROM HSS.F IN HONDO ------
C 
      SUBROUTINE HND_STVDD(RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER RTDB
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=500)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
C
C     ----- MEMORY BLOCK   -----
C
      CALL HND_DDER1_INIT1(RTDB)
C
C     ----- NUCLEAR FORCES -----
C
      CALL HND_NUCDX
C
C     ----- HELLMANN-FEYNMAN FORCES -----
C
      CALL HND_HLFX
C
C     ----- 1E-CONTRIBUTION -----
C
      CALL HND_HLFDX

      CALL HND_TVDX
C
C     ----- MEMORY BLOCK   -----
C
      CALL HND_DDER1_INIT2(RTDB)
C
C     ----- DENSITY FORCES -----
C
      CALL HND_SDX
C
C     ----- MEMORY BLOCK   -----
C
      CALL HND_DDER1_POPM
C
C     ----- SAVE 1E-GRADIENT ON -DAF- FILE AND CHECK TIME -----
C

C     IDAF=10
      CALL HND_HONDO

C
C     ----- OPEN -DA- FILE -10-
C
      NDAR=MXIODA-1
C     CALL HND_DAOPEN(IDAF,IODA,NAV,NDAR)
      NREST=4
      IST=1
      JST=1
      KST=1
      LST=1
C
C     CALL HND_WRTGRD(DE,NREST,IST,JST,KST,LST)
C
      RETURN
      END
C
      subroutine hnd_dder1_init1(rtdb)
      implicit double precision (a-h,o-z)
      parameter (mxatom=500)
#include "mafdecls.fh"
#include "global.fh"
c
c*    logical status
c*    integer g_dens, g_wdens
      integer rtdb 
c 
      common/hnd_iofile/ir,iw,ip
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_dder1_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
      dimension g_dens(3) ! max # of density matrices
      data zero    /0.0d+00/
c
      do i = 1,3
         g_dens(i) = 0
      enddo 
c
      nc=3*nat
      num2=(num*(num+1))/2
      num3= num*num
c     
      need1 = nc*nc+num2*nc+num2+num2
     &        +num2+num2+num2
c*   &        +num3+num3+num3+num3

c     ----- extra memory to hold temporary de, fc, and df -----

      need2 = need1+nc+nc*nc+num2*nc
c
      need = max(need1,need2)
c
c     ----- memory pointers -----
c
c*     if(.not.ma_push_get(mt_dbl,need , 'mem init',i_init,init))
c*    & call errquit('hnd_dder1_init, malloc of init  failed',911)
      if(.not.ma_push_get(mt_dbl,need ,'mem init',i_init,init)) then
         write(iw,*)'The totoal number of basis functions=', num 
         write(iw,*)'The stack memory for dder_1e =', need 
         call errquit('hnd_dder1_init, malloc of init  failed',911)
      endif
c
      i00=init
      i01=i00+ nc  *nc
      i10=i01+ num2*nc
      i20=i10+ num2
      i30=i20+ num2
c
      i40=i30+ num2
      i50=i40+ num2
      i60=i50+ num2
c
c*    i70=i60+ num3
c*    i80=i70+ num3
c*    i90=i80+ num3
c*    i95=i90+ num3
      i70=i60+ 0
      i80=i70+ 0
      i90=i80+ 0
      i95=i90+ 0
c
      i96=i95+ nc
      i97=i96+ nc*nc
      i98=i97+ num2*nc
c
      do 10 i=1,need
   10 x(i+init-1)=zero
c
c*    call hnd_dens_path1(rtdb,g_dens,g_wdens)
c
c*    call ga_get(g_dens(1),1,num,1,num, x(i60),num)
c
c*    call ga_get(g_wdens  ,1,num,1,num, x(i80),num)

c*    call hnd_nwhnd_tran(x(i60),x(i70),num)

c*    ij = 0
c*    do i = 1,num
c*       do j = 1,i
c*          ij = ij + 1
c*          x(i10 +ij-1) = x(i70 + (j-1)*(num)+(i-1))
c*       enddo
c*    enddo
c%%%
c*    call hnd_prtr(x(i10),num)
      call hnd_get_dens(rtdb,x(i10),num)
c%%%
c
c*    status = ga_destroy ( g_wdens )
c*    do i=1, 3
c*      status = ga_destroy ( g_dens(i) )
c*    enddo
c
      return
      end
c
      subroutine hnd_dder1_init2(rtdb)
      implicit double precision (a-h,o-z)
      parameter (mxatom=500)
c
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
c 
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_dder1_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
c
      equivalence (x(1),dbl_mb(1))
      data zero    /0.0d+00/
c
      nc=3*nat
      num2=(num*(num+1))/2
      num3= num*num
c
      i00=init
      i01=i00+ nc  *nc
      i10=i01+ num2*nc
      i20=i10+ num2
      i30=i20+ num2
c
      i40=i30+ num2
      i50=i40+ num2
      i60=i50+ num2
c
c*    i70=i60+ num3
c*    i80=i70+ num3
c*    i90=i80+ num3
c*    i95=i90+ num3
      i70=i60+ 0
      i80=i70+ 0
      i90=i80+ 0
      i95=i90+ 0
c
c*    call hnd_nwhnd_tran(x(i80),x(i90),num)
c
c*    ij = 0
c*     do i = 1,num
c*       do j = 1,i
c*          ij = ij + 1
c*          x(i10 +ij-1) = x(i90 + (j-1)*(num)+(i-1))
c*       enddo
c*    enddo
      call hnd_get_wdens(rtdb,x(i10),num)
c
      return
      end
c
c     logical function gradients(rtdb)
      subroutine hnd_dens_path1(rtdb,g_dens,g_wdens)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
c
c     gradients module.
c
c     Context is '...:gradients'
c
c     Assumes SCF has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c
      integer g_dens, g_wdens 
c
      logical int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status
      character*255 title
c
      character*255 movecs_in, movecs_out
c
      dimension g_dens(3) ! max # of density matrices
c
      status = rtdb_get(rtdb, 'scf:converged', MT_LOG, 1, converged) 
      if ( .not.( status .and. converged )) call errquit
     $     ('gradients: no converged SCF wavefunction available', 0)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('gradients: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('gradients: no ao basis ', 0)
**       call int_normalize(basis)

c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'scf:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
            movecs_out = 'movecs'
         else
            movecs_out = movecs_in
         endif
      endif
c
c
c initialize for schwarz screening
c
      call intd_init(rtdb,1,basis)
      call schwarz_init (geom, basis)
c
      call hnd_dens_path2(rtdb,basis,geom,movecs_out,g_dens,g_wdens)
c
      call schwarz_tidy()

      call intd_terminate()
c
c gradients is done destroy basis and geomtry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('gradients:error destroying geom and basis handles',911)
c
      return
      end
c
      subroutine hnd_dens_path2 ( rtdb, basis, geom, movecs,
     &                            g_dens,g_wdens)
c     subroutine grad_force ( rtdb, basis, geom, movecs )
c
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      integer rtdb, basis, geom
      character*(*) movecs
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_dens,     ! density matrices (up to 3)
     $     g_wdens,       ! energy weighted density
     $     g_eigen_diag   ! eigenvalue diagonal matrix or lagrangian
  
      integer nat, nbf, max_sh_bf, max_at_bf, nopen, nclosed, 
     $     ndens, nbf2, nocc, i
c
      integer l_evals, k_evals, l_occ, k_occ

      dimension g_dens(3) ! max # of density matrices

      character*4 scftype
      logical status

C     get information about basis set
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('path2: could not get natoms',0)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('path2: could not get nbf',0)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('path2: could not get max_sh_bf',0)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('path2: could not get max_at_bf',0)
c
C     get information about type of calculation
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen)

      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
         if (nopen .eq. 0) then
            scftype = 'RHF'
         else
            scftype = 'ROHF'
         endif
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2 = 2 * nbf
        ndens = 2
      else if ( scftype .eq. 'RHF' ) then
        nbf2 = nbf
        ndens = 2
      else    ! ROHF
        nbf2 = nbf
        ndens = 3
      endif

C  allocate and initialize global and local memory

C     energy weighted density (NxN)
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')

C     eigenvalue diagonal matrix
      nocc = nopen + nclosed
      status = ga_create ( MT_DBL, nocc, nocc,'eigen_diag', 1, nocc, 
     $                       g_eigen_diag )
      call ga_zero ( g_eigen_diag )

C     eigenvalues
      status = MA_alloc_get(MT_DBL, nbf2,'MO evals', l_evals, k_evals)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_evals',1)
      
C     occupation numbers
      status = MA_alloc_get(MT_DBL, nbf2,'occ. numbers', l_occ, k_occ)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_occ',1)

C     global density
      do i=1, ndens
        g_dens(i) = ga_create_atom_blocked ( geom, basis, 
     $              'density matrix' )
      enddo

c
      call hnd_grad_dens (geom, basis, g_dens, g_wdens, g_eigen_diag, 
     $                dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2, 
     $                nopen, nclosed, scftype, movecs)

C     free temporary arrays
      status = ga_destroy ( g_eigen_diag )
      if (.not.status) 
     $     call errquit('grad: could not destroy g_eigen_diag',1)
c
      status = MA_free_heap ( l_occ )
      status = MA_free_heap ( l_evals )
c     status = MA_free_heap ( l_act )

c include in HND_DDER1_POPM !

c     status = ga_destroy ( g_wdens )
c     do i=1, ndens
c       status = ga_destroy ( g_dens(i) )
c     enddo

      return
      end
C
      subroutine hnd_dder1_popm 
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
      common/hnd_dder1_mem/i_init,init

      if(.not.ma_pop_stack(i_init))
     & call errquit('hnd_dder1_init, ma_pop_stack of init
     & failed',911)

      return
      end
C
      SUBROUTINE HND_DDTXYZ(XT,YT,ZT,XS,YS,ZS,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XT(7,*),YT(7,*),ZT(7,*)
      DIMENSION XS(7,*),YS(7,*),ZS(7,*)
C
      DO 10 I=1,NI
      XT(I,1)=(XS(I,1  )            -XS(I,3  )*(AJ+AJ))*AJ
      YT(I,1)=(YS(I,1  )            -YS(I,3  )*(AJ+AJ))*AJ
      ZT(I,1)=(ZS(I,1  )            -ZS(I,3  )*(AJ+AJ))*AJ
   10 CONTINUE
C
      IF(NJ.EQ.1) RETURN
C
      DO 20 I=1,NI
      XT(I,2)=(XS(I,2  )*DBLE(2+2-1)-XS(I,4  )*(AJ+AJ))*AJ
      YT(I,2)=(YS(I,2  )*DBLE(2+2-1)-YS(I,4  )*(AJ+AJ))*AJ
      ZT(I,2)=(ZS(I,2  )*DBLE(2+2-1)-ZS(I,4  )*(AJ+AJ))*AJ
   20 CONTINUE
C
      IF(NJ.EQ.2) RETURN
C
      DO 30 J=3,NJ
      DO 30 I=1,NI
      XT(I,J)=(XS(I,J  )*DBLE(J+J-1)-XS(I,J+2)*(AJ+AJ))*AJ
     1        -XS(I,J-2)*DBLE(((J-1)*(J-2))/2)
      YT(I,J)=(YS(I,J  )*DBLE(J+J-1)-YS(I,J+2)*(AJ+AJ))*AJ
     1        -YS(I,J-2)*DBLE(((J-1)*(J-2))/2)
      ZT(I,J)=(ZS(I,J  )*DBLE(J+J-1)-ZS(I,J+2)*(AJ+AJ))*AJ
     1        -ZS(I,J-2)*DBLE(((J-1)*(J-2))/2)
   30 CONTINUE
C
      RETURN
      END

      SUBROUTINE HND_DDERIJ(DDXDII,DDYDII,DDZDII,DDXDIJ,DDYDIJ,DDZDIJ,
     1                  DXDI  ,DYDI  ,DZDI  ,DXDJ  ,DYDJ  ,DZDJ  ,
     2                  X,Y,Z,LIT,LJT,AI,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION      X(7,*),     Y(7,*),     Z(7,*)
      DIMENSION   DXDI(5,1),  DYDI(5,1),  DZDI(5,1)
      DIMENSION   DXDJ(5,1),  DYDJ(5,1),  DZDJ(5,1)
      DIMENSION DDXDII(5,1),DDYDII(5,1),DDZDII(5,1)
      DIMENSION DDXDIJ(5,1),DDYDIJ(5,1),DDZDIJ(5,1)
C
C     ----- DXDI ... -----
C
      DO 110 J=1,LJT
      DXDI(1,J)=X(2,J)*(AI+AI)
      DYDI(1,J)=Y(2,J)*(AI+AI)
      DZDI(1,J)=Z(2,J)*(AI+AI)
  110 CONTINUE
C
      IF(LIT.EQ.1) GO TO 130
C
      DO 120 I=2,LIT
      DO 120 J=1,LJT
      DXDI(I,J)=X(I+1,J)*(AI+AI)-X(I-1,J)*DBLE(I-1)
      DYDI(I,J)=Y(I+1,J)*(AI+AI)-Y(I-1,J)*DBLE(I-1)
      DZDI(I,J)=Z(I+1,J)*(AI+AI)-Z(I-1,J)*DBLE(I-1)
  120 CONTINUE
C
  130 CONTINUE
C
C     ----- DXDJ ... -----
C
      DO 140 I=1,LIT
      DXDJ(I,1)=X(I,2)*(AJ+AJ)
      DYDJ(I,1)=Y(I,2)*(AJ+AJ)
      DZDJ(I,1)=Z(I,2)*(AJ+AJ)
  140 CONTINUE
C
      IF(LJT.EQ.1) GO TO 200
C
      DO 150 J=2,LJT
      DO 150 I=1,LIT
      DXDJ(I,J)=X(I,J+1)*(AJ+AJ)-X(I,J-1)*DBLE(J-1)
      DYDJ(I,J)=Y(I,J+1)*(AJ+AJ)-Y(I,J-1)*DBLE(J-1)
      DZDJ(I,J)=Z(I,J+1)*(AJ+AJ)-Z(I,J-1)*DBLE(J-1)
  150 CONTINUE
C
  200 CONTINUE
C
C     ----- DDXDII ... -----
C
      DO 210 J=1,LJT
      DDXDII(1,J)=(X(3,J)*(AI+AI)-X(1,J))*(AI+AI)
      DDYDII(1,J)=(Y(3,J)*(AI+AI)-Y(1,J))*(AI+AI)
      DDZDII(1,J)=(Z(3,J)*(AI+AI)-Z(1,J))*(AI+AI)
  210 CONTINUE
C
      IF(LIT.EQ.1) GO TO 300
C
      DO 220 J=1,LJT
      DDXDII(2,J)=(X(4,J)*(AI+AI)-X(2,J)*DBLE(2+2-1))*(AI+AI)
      DDYDII(2,J)=(Y(4,J)*(AI+AI)-Y(2,J)*DBLE(2+2-1))*(AI+AI)
      DDZDII(2,J)=(Z(4,J)*(AI+AI)-Z(2,J)*DBLE(2+2-1))*(AI+AI)
  220 CONTINUE
C
      IF(LIT.EQ.2) GO TO 300
C
      DO 230 I=3,LIT
      DO 230 J=1,LJT
      DDXDII(I,J)=(X(I+2,J)*(AI+AI)-X(I,J)*DBLE(I+I-1))*(AI+AI)
     1            +X(I-2,J)*DBLE((I-1)*(I-2))
      DDYDII(I,J)=(Y(I+2,J)*(AI+AI)-Y(I,J)*DBLE(I+I-1))*(AI+AI)
     1            +Y(I-2,J)*DBLE((I-1)*(I-2))
      DDZDII(I,J)=(Z(I+2,J)*(AI+AI)-Z(I,J)*DBLE(I+I-1))*(AI+AI)
     1            +Z(I-2,J)*DBLE((I-1)*(I-2))
  230 CONTINUE
C
  300 CONTINUE
C
C     ----- DDXDIJ ... -----
C
      DDXDIJ(1,1)= X(2,2)*((AJ+AJ)*(AI+AI))
      DDYDIJ(1,1)= Y(2,2)*((AJ+AJ)*(AI+AI))
      DDZDIJ(1,1)= Z(2,2)*((AJ+AJ)*(AI+AI))
C
      IF(LJT.EQ.1) GO TO 330
      DO 320 J=2,LJT
      DDXDIJ(1,J)=(X(2,J+1)*(AJ+AJ)-X(2,J-1)*DBLE(J-1))*(AI+AI)
      DDYDIJ(1,J)=(Y(2,J+1)*(AJ+AJ)-Y(2,J-1)*DBLE(J-1))*(AI+AI)
      DDZDIJ(1,J)=(Z(2,J+1)*(AJ+AJ)-Z(2,J-1)*DBLE(J-1))*(AI+AI)
  320 CONTINUE
  330 IF(LIT.EQ.1) RETURN
      DO 340 I=2,LIT
      DDXDIJ(I,1)=(X(I+1,2)*(AI+AI)-X(I-1,2)*DBLE(I-1))*(AJ+AJ)
      DDYDIJ(I,1)=(Y(I+1,2)*(AI+AI)-Y(I-1,2)*DBLE(I-1))*(AJ+AJ)
      DDZDIJ(I,1)=(Z(I+1,2)*(AI+AI)-Z(I-1,2)*DBLE(I-1))*(AJ+AJ)
  340 CONTINUE
      IF(LJT.EQ.1) RETURN
      DO 350 I=2,LIT
      DO 350 J=2,LJT
      DDXDIJ(I,J)=
     1      (X(I+1,J+1)*(AJ+AJ)-X(I+1,J-1)*DBLE(J-1))*(AI+AI)
     2     -(X(I-1,J+1)*(AJ+AJ)-X(I-1,J-1)*DBLE(J-1))*DBLE(I-1)
      DDYDIJ(I,J)=
     1      (Y(I+1,J+1)*(AJ+AJ)-Y(I+1,J-1)*DBLE(J-1))*(AI+AI)
     2     -(Y(I-1,J+1)*(AJ+AJ)-Y(I-1,J-1)*DBLE(J-1))*DBLE(I-1)
      DDZDIJ(I,J)=
     1      (Z(I+1,J+1)*(AJ+AJ)-Z(I+1,J-1)*DBLE(J-1))*(AI+AI)
     2     -(Z(I-1,J+1)*(AJ+AJ)-Z(I-1,J-1)*DBLE(J-1))*DBLE(I-1)
  350 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_NUCDX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      LOGICAL OUT
C     LOGICAL SOME
      CHARACTER*4 KEYEFC
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      DIMENSION DDE(1)
      EQUIVALENCE (DDE(1),DBL_MB(1))
      DIMENSION DNAM(3)
C
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA THREE  /3.0D+00/
      DATA KEYEFC /' EFC'/
      DATA KEYFLD /4H FLD/
C
      OUT =.FALSE.       
C
C     ----- ZERO OUT GRADIENT+HESSIAN ARRAYS -----
C
      DO 10 IAT=1,NAT
      DO 10 IXYZ=1,3
   10 DE(IXYZ,IAT)=ZERO
      NC=3*NAT
      NC2=NC*NC
      DO 20 IXYZ=1,NC2
   20 DDE(INIT-1+IXYZ)=ZERO
C
      IF(NAT.EQ.1) GO TO 200
C
C     ----- CALCULATE FORCES -----
C
      DO 140 IAT=2,NAT
      IC=3*(IAT-1)
      DO 140 JAT=1,IAT-1
      JC=3*(JAT-1)
C
      RIJ=ZERO
      DO 110 IXYZ=1,3
  110 RIJ=RIJ+(C(IXYZ,IAT)-C(IXYZ,JAT))**2
C
C     ----- DEDI AND DEDJ -----
C
      DUM=ZAN(IAT)*ZAN(JAT)/(RIJ* SQRT(RIJ))
      DO 120 IXYZ=1,3
      DUM1=(C(IXYZ,IAT)-C(IXYZ,JAT))*DUM
      DE(IXYZ,IAT)=DE(IXYZ,IAT)-DUM1
      DE(IXYZ,JAT)=DE(IXYZ,JAT)+DUM1
  120 CONTINUE
C
C     ----- DDEDII , DDEDIJ , DDEDJJ -----
C
      DUM=DUM/RIJ
      DO 130 IXYZ=1,3
      DUM1=(THREE*(C(IXYZ,IAT)-C(IXYZ,JAT))**2-RIJ)*DUM
C     ICIC=(IC+IXYZ)+NC*(IC+IXYZ-1)
C     ICJC=(IC+IXYZ)+NC*(JC+IXYZ-1)
C     JCJC=(JC+IXYZ)+NC*(JC+IXYZ-1)
      ICIC=(IC+IXYZ)+NC*(IC+IXYZ-1) + INIT -1
      ICJC=(IC+IXYZ)+NC*(JC+IXYZ-1) + INIT -1
      JCJC=(JC+IXYZ)+NC*(JC+IXYZ-1) + INIT -1
      DDE(ICIC)=DDE(ICIC)+DUM1
      DDE(ICJC)=DDE(ICJC)-DUM1
      DDE(JCJC)=DDE(JCJC)+DUM1
  130 CONTINUE
      DUM1=(THREE*(C(1,IAT)-C(1,JAT))*(C(2,IAT)-C(2,JAT)))*DUM
C     IC2IC1=(IC+2)+NC*(IC+1-1)
C     IC2JC1=(IC+2)+NC*(JC+1-1)
C     IC1JC2=(IC+1)+NC*(JC+2-1)
C     JC2JC1=(JC+2)+NC*(JC+1-1)
      IC2IC1=(IC+2)+NC*(IC+1-1)+INIT-1
      IC2JC1=(IC+2)+NC*(JC+1-1)+INIT-1
      IC1JC2=(IC+1)+NC*(JC+2-1)+INIT-1
      JC2JC1=(JC+2)+NC*(JC+1-1)+INIT-1
      DDE(IC2IC1)=DDE(IC2IC1)+DUM1
      DDE(IC2JC1)=DDE(IC2JC1)-DUM1
      DDE(IC1JC2)=DDE(IC1JC2)-DUM1
      DDE(JC2JC1)=DDE(JC2JC1)+DUM1
      DUM1=(THREE*(C(1,IAT)-C(1,JAT))*(C(3,IAT)-C(3,JAT)))*DUM
      IC3IC1=(IC+3)+NC*(IC+1-1)+INIT-1
      IC3JC1=(IC+3)+NC*(JC+1-1)+INIT-1
      IC1JC3=(IC+1)+NC*(JC+3-1)+INIT-1
      JC3JC1=(JC+3)+NC*(JC+1-1)+INIT-1
      DDE(IC3IC1)=DDE(IC3IC1)+DUM1
      DDE(IC3JC1)=DDE(IC3JC1)-DUM1
      DDE(IC1JC3)=DDE(IC1JC3)-DUM1
      DDE(JC3JC1)=DDE(JC3JC1)+DUM1
      DUM1=(THREE*(C(2,IAT)-C(2,JAT))*(C(3,IAT)-C(3,JAT)))*DUM
      IC3IC2=(IC+3)+NC*(IC+2-1)+INIT-1
      IC3JC2=(IC+3)+NC*(JC+2-1)+INIT-1
      IC2JC3=(IC+2)+NC*(JC+3-1)+INIT-1
      JC3JC2=(JC+3)+NC*(JC+2-1)+INIT-1
      DDE(IC3IC2)=DDE(IC3IC2)+DUM1
      DDE(IC3JC2)=DDE(IC3JC2)-DUM1
      DDE(IC2JC3)=DDE(IC2JC3)-DUM1
      DDE(JC3JC2)=DDE(JC3JC2)+DUM1
C
  140 CONTINUE
C
C     ----- SYMMETRIZE -----
C
      DO 150 IC=1,NC
      DO 150 JC=1,IC
      IJC=IC+NC*(JC-1)+INIT-1
      JIC=JC+NC*(IC-1)+INIT-1
  150 DDE(JIC)=DDE(IJC)
C
C     ----- PRINT -----
C
C NOTE: OUT CLOSED; FOR PARAL. CAL. WRITE NEED CONDITIONED
  200 CONTINUE
      IF(.NOT.OUT) GO TO 300
C
      WRITE(IW,9996)
      MAX=0
  210 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9999)
      WRITE(IW,9998) (IAT,IAT=MIN,MAX)
      WRITE(IW,9999)
      DO 220 IXYZ=1,3
  220 WRITE(IW,9997) DNAM(IXYZ),(DE(IXYZ,IAT),IAT=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 210
C
      WRITE(IW,9995)
      CALL HND_HSSOUT(DDE(INIT),NC,NC,NC)
  300 CONTINUE
      RETURN
 9999 FORMAT(/)
 9998 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9997 FORMAT(6X,A3,8F15.7)
 9996 FORMAT(/,10X,32(1H-),/,10X,'NUCLEAR CONTRIBUTION TO GRADIENT',
     1       /,10X,32(1H-))
 9995 FORMAT(/,10X,31(1H-),/,10X,'NUCLEAR CONTRIBUTION TO HESSIAN',
     1       /,10X,31(1H-))
C
      END

      SUBROUTINE HND_HLFX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "global.fh"
      PARAMETER (MXATOM=500)
      CHARACTER*8 ERRMSG
      CHARACTER*8 SCF,UHF
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      DIMENSION DNAM(3)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL SOME,OUT
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN','- HLFX -'/
      DATA SCF,UHF /'SCF     ','UHF    '/
      DATA ZERO    /0.0D+00/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
C
      SOME = .FALSE.
      OUT  = .FALSE.
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
      I00=INIT
      I01=I00+ NC  *NC
      I10=I01+ NUM2*NC
      I20=I10+ NUM2
      I30=I20+ NUM2
C
      I40=I30+ NUM2
      I50=I40+ NUM2
      I60=I50+ NUM2
C
c*    I70=I60+ NUM*NUM
c*    I80=I70+ NUM*NUM
c*    I90=I80+ NUM*NUM
c*    I95=I90+ NUM*NUM
      I70=I60+ 0
      I80=I70+ 0
      I90=I80+ 0
      I95=I90+ 0
C
      I96=I95+ NC
      I97=I96+ NC*NC
      I98=I97+ NUM2*NC
C
C     ----- STORE AWAY 1E-GRADIENT -----
C
      DO 10 I=1,3*NAT
         X(I95+I-1)=DE(I,1)
   10    DE (I,1)=ZERO

C
C     ----- STORE AWAY 1E-FORCE CONSTANT -----
C
      DO 11 I=1,NC*NC
         X(I96+I-1)=X(I00+I-1)
   11    X(I00+I-1)=ZERO     
C
C
C     ----- INITIALIZE FOCK DERIVATIVE MATRICES -----
C
      DO 20 I=1,NUM2*NC
   20 X(I+I01-1)=ZERO
C
C
C     ----- CALCULATE -HELFEY- CONTRIBUTION -----
C
      CALL HND_HLF_NEW(X(I00),NC,X(I01),NUM2,X(I10))
C
C     ----- ADD 1E-GRADIENTS -----
C
      DO 31 I=1,3*NAT
   31    DE(I,1)=DE(I,1)+X(I95+I-1)
C
C     ----- ADD 1E-FORCE CONSTANT -----
C
      DO 32 I=1,NC*NC
   32    X(I00+I-1)=X(I00+I-1)+X(I96+I-1)
C
C     ----- PRINT -----
C
      IF(.NOT.SOME) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
C
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9999)
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
C
         DO 9200 I=1,3
 9200    WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      ENDIF
C
      IF(MAX.LT.NAT) GO TO 9100
C
      IF(.NOT.OUT) GO TO 9300
C     NC=3*NAT
C
      IF (GA_NODEID().EQ.0 ) THEN
C
         WRITE(IW,9991)
         CALL HND_HSSOUT(X(I00),NC,NC,NC)
C
         WRITE(IW,9990)
         DO 9250 IC=1,NC
 9250    CALL HND_PRTR(X(I01+(IC-1)*NUM2),NUM)
C
         WRITE(IW,9998)
C
      ENDIF

 9300 CONTINUE
C
 9999 FORMAT(/,10X,23(1H-),/,10X,'HELLMANN-FEYNMAN FORCES',
     1       /,10X,23(1H-))
 9998 FORMAT(' ...... END OF HELLMANN-FEYNMAN FORCES ......')
 9997 FORMAT(' IN -HLF- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS = ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8E15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,24(1H-),/,10X,'HELLMANN-FEYNMAN HESSIAN',
     1       /,10X,24(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      RETURN
      END
      SUBROUTINE HND_HLF_NEW(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL SOME,OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,*),DF(LDF,*),DAB(*)
      DIMENSION DIJ(225)
      DIMENSION DVIJ(3,225),DDVIJ(3,3,225)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HLFSPD-'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
#include "global.fh"
C
      OUT =.FALSE.
      SOME=.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
C     ----- CALCULATE -HELFEY- TERM -----
C
      IF(SOME) WRITE(IW,9999)
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
         MINI=KMIN(II)
         MAXI=KMAX(II)
         LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
         MINJ=KMIN(JJ)
         MAXJ=KMAX(JJ)
         LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
         IJSH=IJSH+1
         IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 8000
C
         IANDJ=II.EQ.JJ
C
         IJ=0
         DO I=MINI,MAXI
            JMAX=MAXJ
            IF(IANDJ) JMAX=I
            DO J=MINJ,JMAX
               IDIM=LOCI+I
               JDIM=LOCJ+J
               NN=(IDIM*(IDIM-1))/2+JDIM
               DEN=DAB(NN)
               IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
               IJ=IJ+1
               DIJ(IJ)=DEN
            ENDDO
         ENDDO
C
C     ----- HELLMANN-FEYNMAN TERM -----
C
      DO KAT=1,NAT
C
         CALL HND_HLF_IJ(II,JJ,KAT,DVIJ,DDVIJ)
C
         KC=3*(KAT-1)
C
         IJ=0
         DO I=MINI,MAXI
            JMAX=MAXJ
            IF(IANDJ) JMAX=I
            DO J=MINJ,JMAX
               IJ=IJ+1
               DUMX = DVIJ(1,IJ)
               DUMY = DVIJ(2,IJ)
               DUMZ = DVIJ(3,IJ)
               DUMXX=DDVIJ(1,1,IJ)
               DUMYY=DDVIJ(2,2,IJ)
               DUMZZ=DDVIJ(3,3,IJ)
               DUMXY=DDVIJ(2,1,IJ)
               DUMXZ=DDVIJ(3,1,IJ)
               DUMYZ=DDVIJ(3,2,IJ)
C
               IDIM=LOCI+I
               JDIM=LOCJ+J
               NN=(IDIM*(IDIM-1))/2+JDIM
               DF(NN,KC+1)=DF(NN,KC+1)+DUMX
               DF(NN,KC+2)=DF(NN,KC+2)+DUMY
               DF(NN,KC+3)=DF(NN,KC+3)+DUMZ
C     
               DE(1,KAT)=DE(1,KAT)+DUMX*DIJ(IJ)
               DE(2,KAT)=DE(2,KAT)+DUMY*DIJ(IJ)
               DE(3,KAT)=DE(3,KAT)+DUMZ*DIJ(IJ)
C
               FC(KC+1,KC+1)=FC(KC+1,KC+1)+DUMXX*DIJ(IJ)
               FC(KC+2,KC+1)=FC(KC+2,KC+1)+DUMXY*DIJ(IJ)
               FC(KC+2,KC+2)=FC(KC+2,KC+2)+DUMYY*DIJ(IJ)
               FC(KC+3,KC+1)=FC(KC+3,KC+1)+DUMXZ*DIJ(IJ)
               FC(KC+3,KC+2)=FC(KC+3,KC+2)+DUMYZ*DIJ(IJ)
               FC(KC+3,KC+3)=FC(KC+3,KC+3)+DUMZZ*DIJ(IJ)
C
            ENDDO   
         ENDDO   
C
      ENDDO   
C
      IF(OUT) THEN       
         WRITE(IW,9993) II,JJ
         MAX=0
 7500    MIN=MAX+1
         MAX=MAX+8
         IF(MAX.GT.NAT) MAX=NAT
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
         DO I=1,3
            WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
         ENDDO
         IF(MAX.LT.NAT) GO TO 7500
      ENDIF
C
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      L2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DF,(NC*L2),'+')
C
C     ----- SYMMETRIZE -FC- -----
C
      NC=3*NAT
      DO 9010 IC=1,NC
      DO 9010 JC=1,IC
 9010 FC(JC,IC)=FC(IC,JC)
C
C     ----- PRINT -----
C
      IF(.NOT.SOME) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
      ENDIF

      IF (GA_NODEID().EQ.0 ) THEN
         DO 9200 I=1,3
 9200    WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      ENDIF
      IF(MAX.LT.NAT) GO TO 9100
C
      IF(.NOT.OUT) GO TO 9300
      NC=3*NAT
      IF (GA_NODEID().EQ.0 ) THEN
         WRITE(IW,9991)
         CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
      IF (GA_NODEID().EQ.0 ) THEN
         WRITE(IW,9990)
         DO 9250 IC=1,NC
 9250    CALL HND_PRTR (DF(1,IC),NUM)
C
         WRITE(IW,9998)
      ENDIF
 9300 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,23(1H-),/,10X,'HELLMANN-FEYNMAN FORCES',
     1       /,10X,23(1H-))
 9998 FORMAT(' ...... END OF HELLMANN-FEYNMAN FORCES ......')
 9997 FORMAT(' IN -HLF- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS = ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8E15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,24(1H-),/,10X,'HELLMANN-FEYNMAN HESSIAN',
     1       /,10X,24(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      END
      SUBROUTINE HND_HLF_IJ(II,JJ,KAT,DVIJ,DDVIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL SOME,OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      DIMENSION DVIJ(3,225),DDVIJ(3,3,225)
      DIMENSION PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,5,5),  YV(5,5,5),  ZV(5,5,5)
      DIMENSION  DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
      DIMENSION DDXV(5,5,5),DDYV(5,5,5),DDZV(5,5,5)
      DIMENSION W2(5),W4(5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HLFSPD-'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA FOUR   /4.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT =.FALSE.
      SOME=.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
C     ----- CALCULATE -HELFEY- TERM -----
C
      NDER=2
C
C     ----- ISHELL -----
C
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IJ=0
      DO I=MINI,MAXI
         JMAX=MAXJ
         IF(IANDJ) JMAX=I
         DO J=MINJ,JMAX
            IJ=IJ+1
            DO K=1,3
               DVIJ(K,IJ)=ZERO
               DO L=1,3
                  DDVIJ(K,L,IJ)=ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 PIJ(IJ)=DUM2
C
C     ----- HELLMANN-FEYNMAN TERM -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
C
C     ----- KAT -----
C
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      U2=UU
      U4=UU*UU
      WW=W(IROOT)*ZNUC
      W2(IROOT)=WW*U2*TWO
      W4(IROOT)=WW*U4*FOUR
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I

      CALL HND_DSXYZ

      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT

      CALL HND_DVXYZ

      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT

      CALL HND_DDVXYZ

      DDXV(I,J,IROOT)=XINT
      DDYV(I,J,IROOT)=YINT
      DDZV(I,J,IROOT)=ZINT
  410 CONTINUE
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DUMXX=ZERO
      DUMYY=ZERO
      DUMZZ=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX =DUMX+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMY =DUMY+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMZ =DUMZ+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W2(IROOT)
      DUM  =
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMXX=DUMXX-DUM+
     1   DDXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYY=DUMYY-DUM+
     1     XV(IX,JX,IROOT)*DDYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMZZ=DUMZZ-DUM+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*DDZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXY=DUMXY+
     1    DXV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXZ=DUMXZ+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYZ=DUMYZ+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
  430 CONTINUE
      IJ=IJ+1
      DUMX =DUMX *(PI212*AA1*PIJ(IJ))
      DUMY =DUMY *(PI212*AA1*PIJ(IJ))
      DUMZ =DUMZ *(PI212*AA1*PIJ(IJ))
      DVIJ(1,IJ)=DVIJ(1,IJ)+DUMX
      DVIJ(2,IJ)=DVIJ(2,IJ)+DUMY
      DVIJ(3,IJ)=DVIJ(3,IJ)+DUMZ
      DUMXX=DUMXX*(PI212*AA1*PIJ(IJ))
      DUMYY=DUMYY*(PI212*AA1*PIJ(IJ))
      DUMZZ=DUMZZ*(PI212*AA1*PIJ(IJ))
      DUMXY=DUMXY*(PI212*AA1*PIJ(IJ))
      DUMXZ=DUMXZ*(PI212*AA1*PIJ(IJ))
      DUMYZ=DUMYZ*(PI212*AA1*PIJ(IJ))
      DDVIJ(1,1,IJ)=DDVIJ(1,1,IJ)+DUMXX    
      DDVIJ(2,1,IJ)=DDVIJ(2,1,IJ)+DUMXY
      DDVIJ(2,2,IJ)=DDVIJ(2,2,IJ)+DUMYY
      DDVIJ(3,1,IJ)=DDVIJ(3,1,IJ)+DUMXZ
      DDVIJ(3,2,IJ)=DDVIJ(3,2,IJ)+DUMYZ
      DDVIJ(3,3,IJ)=DDVIJ(3,3,IJ)+DUMZZ
      DUM=(PI212*AA1)*PIJ(IJ)
C
  440 CONTINUE
  450 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      RETURN
 9997 FORMAT(' IN -HLF- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS = ',3I3)
      END
      SUBROUTINE HND_HLF_OLD(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL SOME,OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION W2(5),W4(5)
      DIMENSION FC(LFC,1),DF(LDF,1),DAB(1)
      DIMENSION DIJ(225),PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,5,5),  YV(5,5,5),  ZV(5,5,5)
      DIMENSION  DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
      DIMENSION DDXV(5,5,5),DDYV(5,5,5),DDZV(5,5,5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HLFSPD-'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA FOUR   /4.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =.FALSE.
      SOME=.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- CALCULATE -HELFEY- TERM -----
C
      IF(SOME) WRITE(IW,9999)
      NDER=2
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
      IJSH=IJSH+1
      IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 7900
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
cd    NN=IA(LOCI+I)+(LOCJ+J)
      IDIM=LOCI+I
      JDIM=LOCJ+J
      NN=(IDIM*(IDIM-1))/2+JDIM
      DEN=DAB(NN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      IJ=IJ+1
      PIJ(IJ)=DUM2
  360 DIJ(IJ)=DUM2*DEN
C
C     ----- HELLMANN-FEYNMAN TERM -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 KAT=1,NAT
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      U2=UU
      U4=UU*UU
      WW=W(IROOT)*ZNUC
      W2(IROOT)=WW*U2*TWO
      W4(IROOT)=WW*U4*FOUR
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I

      CALL HND_DSXYZ

      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT

      CALL HND_DVXYZ

      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT

      CALL HND_DDVXYZ

      DDXV(I,J,IROOT)=XINT
      DDYV(I,J,IROOT)=YINT
      DDZV(I,J,IROOT)=ZINT
  410 CONTINUE
  420 CONTINUE
C
      KC=3*(KAT-1)
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DUMXX=ZERO
      DUMYY=ZERO
      DUMZZ=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX =DUMX+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMY =DUMY+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMZ =DUMZ+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W2(IROOT)
      DUM  =
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMXX=DUMXX-DUM+
     1   DDXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYY=DUMYY-DUM+
     1     XV(IX,JX,IROOT)*DDYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMZZ=DUMZZ-DUM+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*DDZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXY=DUMXY+
     1    DXV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXZ=DUMXZ+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYZ=DUMYZ+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
  430 CONTINUE
      IJ=IJ+1
      DUM=(PI212*AA1)*DIJ(IJ)
      DE(1,KAT)=DE(1,KAT)+DUMX*DUM
      DE(2,KAT)=DE(2,KAT)+DUMY*DUM
      DE(3,KAT)=DE(3,KAT)+DUMZ*DUM
      FC(KC+1,KC+1)=FC(KC+1,KC+1)+DUMXX*DUM
      FC(KC+2,KC+1)=FC(KC+2,KC+1)+DUMXY*DUM
      FC(KC+2,KC+2)=FC(KC+2,KC+2)+DUMYY*DUM
      FC(KC+3,KC+1)=FC(KC+3,KC+1)+DUMXZ*DUM
      FC(KC+3,KC+2)=FC(KC+3,KC+2)+DUMYZ*DUM
      FC(KC+3,KC+3)=FC(KC+3,KC+3)+DUMZZ*DUM
      DUM=(PI212*AA1)*PIJ(IJ)

      IDIM=LOCI+I
      JDIM=LOCJ+J
      NN=(IDIM*(IDIM-1))/2+JDIM

      DF(NN,KC+1)=DF(NN,KC+1)+DUMX*DUM
      DF(NN,KC+2)=DF(NN,KC+2)+DUMY*DUM
      DF(NN,KC+3)=DF(NN,KC+3)+DUMZ*DUM

  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      IF(.NOT.OUT) GO TO 8000
      WRITE(IW,9993) II,JJ
      MAX=0
 7500 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 7600 I=1,3
 7600 WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 7500
C
 7900 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      L2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DF,(NC*L2),'+')
C
C     ----- SYMMETRIZE -FC- -----
C
      NC=3*NAT
      DO 9010 IC=1,NC
      DO 9010 JC=1,IC
 9010 FC(JC,IC)=FC(IC,JC)
C
C     ----- PRINT -----
C
      IF(.NOT.SOME) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
      ENDIF

      IF (GA_NODEID().EQ.0 ) THEN
         DO 9200 I=1,3
 9200    WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      ENDIF
      IF(MAX.LT.NAT) GO TO 9100
C
      IF(.NOT.OUT) GO TO 9300
      NC=3*NAT
      IF (GA_NODEID().EQ.0 ) THEN
         WRITE(IW,9991)
         CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
      IF (GA_NODEID().EQ.0 ) THEN
         WRITE(IW,9990)
         DO 9250 IC=1,NC
 9250    CALL HND_PRTR (DF(1,IC),NUM)
C
         WRITE(IW,9998)
      ENDIF
 9300 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,23(1H-),/,10X,'HELLMANN-FEYNMAN FORCES',
     1       /,10X,23(1H-))
 9998 FORMAT(' ...... END OF HELLMANN-FEYNMAN FORCES ......')
 9997 FORMAT(' IN -HLF- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS = ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8E15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,24(1H-),/,10X,'HELLMANN-FEYNMAN HESSIAN',
     1       /,10X,24(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      END

      SUBROUTINE HND_HLFDX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
      PARAMETER (MXATOM=500)
      CHARACTER*8 ERRMSG
      CHARACTER*8 SCF,UHF
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL OUT
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','- HLFDX-'/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA ZERO   /0.0D+00/
C
      OUT = .FALSE.
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
      I00=INIT
      I01=I00+ NC  *NC
      I10=I01+ NUM2*NC
      I20=I10+ NUM2
      I30=I20+ NUM2
C
      I40=I30+ NUM2
      I50=I40+ NUM2
      I60=I50+ NUM2
C
c*    I70=I60+ NUM*NUM
c*    I80=I70+ NUM*NUM
c*    I90=I80+ NUM*NUM
c*    I95=I90+ NUM*NUM
      I70=I60+ 0
      I80=I70+ 0
      I90=I80+ 0
      I95=I90+ 0
C
      I96=I95+ NC
      I97=I96+ NC*NC
      I98=I97+ NUM2*NC
C
C     ----- STORE AWAY 1E-FORCE CONSTANT -----
C
      DO 11 I=1,NC*NC
         X(I96+I-1)=X(I00+I-1)
   11    X(I00+I-1)=ZERO     
C
C     ----- GET -HLFD- DERIVATIVES -----
C
      CALL HND_HLFD_NEW(X(I00),NC,X(I01),NUM2,X(I10))
C
C     ----- ADD 1E-FORCE CONSTANT -----
C
      DO 31 I=1,NC*NC
   31    X(I00+I-1)=X(I00+I-1)+X(I96+I-1)
C
C     ----- PRINT OUT -----
C
      IF(.NOT.OUT) GO TO 9300
C     NC=3*NAT
      IF (GA_NODEID().EQ.0 ) THEN
         WRITE(IW,9991)
         CALL HND_HSSOUT(X(I00),NC,NC,NC)
C
         WRITE(IW,9998)
      ENDIF
 9300 CONTINUE
C
      RETURN
 9998 FORMAT(' ...... END OF -HLF- DERIVATIVE TERM ......')
 9991 FORMAT(/,10X,19(1H-),/,10X,'-HLFD- HESSIAN TERM',
     1       /,10X,19(1H-))
      END
      SUBROUTINE HND_HLFD_NEW(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,1),DF(LDF,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION DIDVIJ(3,3,225)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- HLFD -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
C     ----- CALCULATE DERIVATIVES OF -HLF- TERM -----
C
      IF(OUT) WRITE(IW,9999)
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,NSHELL
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
      IJSH=IJSH+1
      IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 8000
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IDIM=MAX0(LOCI+I,LOCJ+J)
            JDIM=MIN0(LOCI+I,LOCJ+J)
            NN=(IDIM*(IDIM-1))/2+JDIM
            DEN=DAB(NN)
            DEN=DEN+DEN
            IJ=IJ+1
            DIJ(IJ)=DEN
         ENDDO    
      ENDDO
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      DO KAT=1,NAT
C
      CALL HND_HLFD_IJ(II,JJ,IAT,KAT,DIDVIJ)
C
      IF(KAT.EQ.IAT) THEN 
         FAC=TWO       
      ELSE
         FAC=ONE    
      ENDIF
      IC=3*(IAT-1)
      KC=3*(KAT-1)
      IJ=0
      DO I=MINI,MAXI
      DO J=MINJ,MAXJ
      IJ=IJ+1
      DUMXX=DIDVIJ(1,1,IJ)
      DUMXY=DIDVIJ(1,2,IJ)
      DUMXZ=DIDVIJ(1,3,IJ)
      DUMYX=DIDVIJ(2,1,IJ)
      DUMYY=DIDVIJ(2,2,IJ)
      DUMYZ=DIDVIJ(2,3,IJ)
      DUMZX=DIDVIJ(3,1,IJ)
      DUMZY=DIDVIJ(3,2,IJ)
      DUMZZ=DIDVIJ(3,3,IJ)
      IF(IC.GE.KC) THEN
         FC(IC+1,KC+1)=FC(IC+1,KC+1)+DUMXX*(FAC*DIJ(IJ))
         FC(IC+1,KC+2)=FC(IC+1,KC+2)+DUMXY*(FAC*DIJ(IJ))
         FC(IC+1,KC+3)=FC(IC+1,KC+3)+DUMXZ*(FAC*DIJ(IJ))
         FC(IC+2,KC+1)=FC(IC+2,KC+1)+DUMYX*(FAC*DIJ(IJ))
         FC(IC+2,KC+2)=FC(IC+2,KC+2)+DUMYY*(FAC*DIJ(IJ))
         FC(IC+2,KC+3)=FC(IC+2,KC+3)+DUMYZ*(FAC*DIJ(IJ))
         FC(IC+3,KC+1)=FC(IC+3,KC+1)+DUMZX*(FAC*DIJ(IJ))
         FC(IC+3,KC+2)=FC(IC+3,KC+2)+DUMZY*(FAC*DIJ(IJ))
         FC(IC+3,KC+3)=FC(IC+3,KC+3)+DUMZZ*(FAC*DIJ(IJ))
      ELSE
         FC(KC+1,IC+1)=FC(KC+1,IC+1)+DUMXX*(FAC*DIJ(IJ))
         FC(KC+1,IC+2)=FC(KC+1,IC+2)+DUMYX*(FAC*DIJ(IJ))
         FC(KC+1,IC+3)=FC(KC+1,IC+3)+DUMZX*(FAC*DIJ(IJ))
         FC(KC+2,IC+1)=FC(KC+2,IC+1)+DUMXY*(FAC*DIJ(IJ))
         FC(KC+2,IC+2)=FC(KC+2,IC+2)+DUMYY*(FAC*DIJ(IJ))
         FC(KC+2,IC+3)=FC(KC+2,IC+3)+DUMZY*(FAC*DIJ(IJ))
         FC(KC+3,IC+1)=FC(KC+3,IC+1)+DUMXZ*(FAC*DIJ(IJ))
         FC(KC+3,IC+2)=FC(KC+3,IC+2)+DUMYZ*(FAC*DIJ(IJ))
         FC(KC+3,IC+3)=FC(KC+3,IC+3)+DUMZZ*(FAC*DIJ(IJ))
      ENDIF
      ENDDO   
      ENDDO    
C
      ENDDO
C
      IF(OUT) THEN        
         WRITE(IW,9993) II,JJ
         NC=3*NAT
         CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF FC FROM ALL PROCESSORS -----
C
      NC=3*NAT
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- SYMMETRIZE -----
C
      DO 9020 IAT=1,NAT
      IC=3*(IAT-1)
      DO 9010 IXYZ=1,3
      DO 9010 JXYZ=1,IXYZ
      FCIJ=FC(IC+IXYZ,IC+JXYZ)
      FCJI=FC(IC+JXYZ,IC+IXYZ)
 9010 FC(IC+IXYZ,IC+JXYZ)=(FCIJ+FCJI)*PT5
 9020 CONTINUE
      NC=3*NAT
      DO 9030 IC=1,NC
      DO 9030 JC=1,IC
      FC(JC,IC)=FC(IC,JC)
 9030 CONTINUE
C
C     ----- PRINT OUT -----
C
      IF(OUT) THEN        
         NC=3*NAT
         WRITE(IW,9991)
         CALL HND_HSSOUT(FC,NC,NC,NC)
         WRITE(IW,9998)
      ENDIF      
C
      RETURN
 9999 FORMAT(/,10X,21(1H-),/,10X,'-HLF- DERIVATIVE TERM',
     1       /,10X,21(1H-))
 9998 FORMAT(' ...... END OF -HLF- DERIVATIVE TERM ......')
 9997 FORMAT(' IN -HLFD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8E15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,19(1H-),/,10X,'-HLFD- HESSIAN TERM',
     1       /,10X,19(1H-))
      END
      SUBROUTINE HND_HLFD_IJ(II,JJ,IAT,KAT,DIDVIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      DIMENSION PIJ(225)
      DIMENSION DIDVIJ(3,3,225)
      DIMENSION W2(5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XV(6,5,5),    YV(6,5,5),    ZV(6,5,5)
      DIMENSION    DXV(6,5,5),   DYV(6,5,5),   DZV(6,5,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION DDXVDI(5,5,5),DDYVDI(5,5,5),DDZVDI(5,5,5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- HLFD -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
C     ----- CALCULATE DERIVATIVES OF -HLF- TERM -----
C
      NDER=2
C
C     ----- ISHELL -----
C
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER-1
      IAT   =I
C
C     ----- JSHELL -----
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      JAT   =J
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            DO KXYZ=1,3
               DO IXYZ=1,3
                  DIDVIJ(IXYZ,KXYZ,IJ)=ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
      PIJ(IJ)=DUM2
  360 CONTINUE
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
C                      
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      W2(IROOT)=WW*UU*TWO
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
C
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT
      CALL HND_DVXYZ
      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT
  410 CONTINUE
C
      CALL HND_DERI( DXVDI(1,1,IROOT), DYVDI(1,1,IROOT),
     1              DZVDI(1,1,IROOT),
     1              XV(1,1,IROOT),    YV(1,1,IROOT),    ZV(1,1,IROOT),
     2          LIT,LJT,AI)
      CALL HND_DERI(DDXVDI(1,1,IROOT),DDYVDI(1,1,IROOT),
     1             DDZVDI(1,1,IROOT),
     1             DXV(1,1,IROOT),   DYV(1,1,IROOT),   DZV(1,1,IROOT),
     2          LIT,LJT,AI)
C
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 440 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMXX=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYX=ZERO
      DUMYY=ZERO
      DUMYZ=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DUMZZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMXX=DUMXX+W2(IROOT)*
     1     (DDXVDI(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMXY=DUMXY+W2(IROOT)*
     1     (DXVDI (IX,JX,IROOT)*DYV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMXZ=DUMXZ+W2(IROOT)*
     1     (DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZV   (IZ,JZ,IROOT))
      DUMYX=DUMYX+W2(IROOT)*
     1     (DXV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMYY=DUMYY+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DDYVDI(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMYZ=DUMYZ+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZV   (IZ,JZ,IROOT))
      DUMZX=DUMZX+W2(IROOT)*
     1     (DXV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT))
      DUMZY=DUMZY+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DYV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT))
      DUMZZ=DUMZZ+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DDZVDI(IZ,JZ,IROOT))
  430 CONTINUE
      IJ=IJ+1
      DIDVIJ(1,1,IJ)=DIDVIJ(1,1,IJ)+DUMXX*(AA1*PI212*PIJ(IJ))
      DIDVIJ(1,2,IJ)=DIDVIJ(1,2,IJ)+DUMXY*(AA1*PI212*PIJ(IJ))
      DIDVIJ(1,3,IJ)=DIDVIJ(1,3,IJ)+DUMXZ*(AA1*PI212*PIJ(IJ))
      DIDVIJ(2,1,IJ)=DIDVIJ(2,1,IJ)+DUMYX*(AA1*PI212*PIJ(IJ))
      DIDVIJ(2,2,IJ)=DIDVIJ(2,2,IJ)+DUMYY*(AA1*PI212*PIJ(IJ))
      DIDVIJ(2,3,IJ)=DIDVIJ(2,3,IJ)+DUMYZ*(AA1*PI212*PIJ(IJ))
      DIDVIJ(3,1,IJ)=DIDVIJ(3,1,IJ)+DUMZX*(AA1*PI212*PIJ(IJ))
      DIDVIJ(3,2,IJ)=DIDVIJ(3,2,IJ)+DUMZY*(AA1*PI212*PIJ(IJ))
      DIDVIJ(3,3,IJ)=DIDVIJ(3,3,IJ)+DUMZZ*(AA1*PI212*PIJ(IJ))
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      RETURN
 9997 FORMAT(' IN -HLFD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
      END
      SUBROUTINE HND_HLFD_OLD(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION W2(5)
      DIMENSION FC(LFC,1),DF(LDF,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XV(6,5,5),    YV(6,5,5),    ZV(6,5,5)
      DIMENSION    DXV(6,5,5),   DYV(6,5,5),   DZV(6,5,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION DDXVDI(5,5,5),DDYVDI(5,5,5),DDZVDI(5,5,5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- HLFD -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- CALCULATE DERIVATIVES OF -HLF- TERM -----
C
      IF(OUT) WRITE(IW,9999)
      NDER=2
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER-1
      IAT   =I
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
      IJSH=IJSH+1
      IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 7900
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
C     NN=IA(MAX0(LOCI+I,LOCJ+J))
C    1  +   MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
      DEN=DAB(NN)
      DEN=DEN+DEN
      IJ=IJ+1
      DIJ(IJ)=DUM2*DEN
  360 CONTINUE
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 KAT=1,NAT
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      W2(IROOT)=WW*UU*TWO
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
C
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT
      CALL HND_DVXYZ
      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT
  410 CONTINUE
C
      CALL HND_DERI( DXVDI(1,1,IROOT), DYVDI(1,1,IROOT),
     1              DZVDI(1,1,IROOT),
     1              XV(1,1,IROOT),    YV(1,1,IROOT),    ZV(1,1,IROOT),
     2          LIT,LJT,AI)
      CALL HND_DERI(DDXVDI(1,1,IROOT),DDYVDI(1,1,IROOT),
     1             DDZVDI(1,1,IROOT),
     1             DXV(1,1,IROOT),   DYV(1,1,IROOT),   DZV(1,1,IROOT),
     2          LIT,LJT,AI)
C
  420 CONTINUE
C
      FAC=PI212
      IF(KAT.EQ.IAT) FAC=FAC+FAC
      IC=3*(IAT-1)
      KC=3*(KAT-1)
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 440 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMXX=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYX=ZERO
      DUMYY=ZERO
      DUMYZ=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DUMZZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMXX=DUMXX+W2(IROOT)*
     1     (DDXVDI(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMXY=DUMXY+W2(IROOT)*
     1     (DXVDI (IX,JX,IROOT)*DYV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMXZ=DUMXZ+W2(IROOT)*
     1     (DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZV   (IZ,JZ,IROOT))
      DUMYX=DUMYX+W2(IROOT)*
     1     (DXV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMYY=DUMYY+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DDYVDI(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT))
      DUMYZ=DUMYZ+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZV   (IZ,JZ,IROOT))
      DUMZX=DUMZX+W2(IROOT)*
     1     (DXV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT))
      DUMZY=DUMZY+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)*DYV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT))
      DUMZZ=DUMZZ+W2(IROOT)*
     1     ( XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DDZVDI(IZ,JZ,IROOT))
  430 CONTINUE
      IJ=IJ+1
      IF(IC.GE.KC) THEN
      FC(IC+1,KC+1)=FC(IC+1,KC+1)+DUMXX*(AA1*FAC*DIJ(IJ))
      FC(IC+1,KC+2)=FC(IC+1,KC+2)+DUMXY*(AA1*FAC*DIJ(IJ))
      FC(IC+1,KC+3)=FC(IC+1,KC+3)+DUMXZ*(AA1*FAC*DIJ(IJ))
      FC(IC+2,KC+1)=FC(IC+2,KC+1)+DUMYX*(AA1*FAC*DIJ(IJ))
      FC(IC+2,KC+2)=FC(IC+2,KC+2)+DUMYY*(AA1*FAC*DIJ(IJ))
      FC(IC+2,KC+3)=FC(IC+2,KC+3)+DUMYZ*(AA1*FAC*DIJ(IJ))
      FC(IC+3,KC+1)=FC(IC+3,KC+1)+DUMZX*(AA1*FAC*DIJ(IJ))
      FC(IC+3,KC+2)=FC(IC+3,KC+2)+DUMZY*(AA1*FAC*DIJ(IJ))
      FC(IC+3,KC+3)=FC(IC+3,KC+3)+DUMZZ*(AA1*FAC*DIJ(IJ))
      ELSE
      FC(KC+1,IC+1)=FC(KC+1,IC+1)+DUMXX*(AA1*FAC*DIJ(IJ))
      FC(KC+1,IC+2)=FC(KC+1,IC+2)+DUMYX*(AA1*FAC*DIJ(IJ))
      FC(KC+1,IC+3)=FC(KC+1,IC+3)+DUMZX*(AA1*FAC*DIJ(IJ))
      FC(KC+2,IC+1)=FC(KC+2,IC+1)+DUMXY*(AA1*FAC*DIJ(IJ))
      FC(KC+2,IC+2)=FC(KC+2,IC+2)+DUMYY*(AA1*FAC*DIJ(IJ))
      FC(KC+2,IC+3)=FC(KC+2,IC+3)+DUMZY*(AA1*FAC*DIJ(IJ))
      FC(KC+3,IC+1)=FC(KC+3,IC+1)+DUMXZ*(AA1*FAC*DIJ(IJ))
      FC(KC+3,IC+2)=FC(KC+3,IC+2)+DUMYZ*(AA1*FAC*DIJ(IJ))
      FC(KC+3,IC+3)=FC(KC+3,IC+3)+DUMZZ*(AA1*FAC*DIJ(IJ))
      ENDIF
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      IF(OUT) THEN        
      WRITE(IW,9993) II,JJ
      NC=3*NAT
      CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
 7900 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF FC FROM ALL PROCESSORS -----
C
      NC=3*NAT
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- SYMMETRIZE -----
C
      DO 9020 IAT=1,NAT
      IC=3*(IAT-1)
      DO 9010 IXYZ=1,3
      DO 9010 JXYZ=1,IXYZ
      FCIJ=FC(IC+IXYZ,IC+JXYZ)
      FCJI=FC(IC+JXYZ,IC+IXYZ)
 9010 FC(IC+IXYZ,IC+JXYZ)=(FCIJ+FCJI)*PT5
 9020 CONTINUE
      NC=3*NAT
      DO 9030 IC=1,NC
      DO 9030 JC=1,IC
      FC(JC,IC)=FC(IC,JC)
 9030 CONTINUE
C
C     ----- PRINT OUT -----
C
      IF(OUT) THEN         
      NC=3*NAT
      WRITE(IW,9991)
      CALL HND_HSSOUT(FC,NC,NC,NC)
      WRITE(IW,9998)
      ENDIF    
C
      RETURN
 9999 FORMAT(/,10X,21(1H-),/,10X,'-HLF- DERIVATIVE TERM',
     1       /,10X,21(1H-))
 9998 FORMAT(' ...... END OF -HLF- DERIVATIVE TERM ......')
 9997 FORMAT(' IN -HLFD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8E15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,19(1H-),/,10X,'-HLFD- HESSIAN TERM',
     1       /,10X,19(1H-))
      END
      SUBROUTINE HND_TVDX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
      PARAMETER (MXATOM=500)
      CHARACTER*8 ERRMSG
      CHARACTER*8 SCF,UHF
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION X(1)
      LOGICAL OUT
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION ERRMSG(3)
      DIMENSION DNAM(3)
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ERRMSG  /'PROGRAM ','STOP IN ','- TVDX -'/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA ZERO   /0.0D+00/
      DATA NFT17  /17/
C
      OUT = .FALSE.
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
      I00=INIT
      I01=I00+ NC  *NC
      I10=I01+ NUM2*NC
      I20=I10+ NUM2
      I30=I20+ NUM2
C
      I40=I30+ NUM2
      I50=I40+ NUM2
      I60=I50+ NUM2
C
c*    I70=I60+ NUM*NUM
c*    I80=I70+ NUM*NUM
c*    I90=I80+ NUM*NUM
c*    I95=I90+ NUM*NUM
      I70=I60+ 0
      I80=I70+ 0
      I90=I80+ 0
      I95=I90+ 0
C
      I96=I95+ NC
      I97=I96+ NC*NC
      I98=I97+ NUM2*NC
C
C     ----- STORE AWAY 1E-GRADIENT -----
C
      DO 10 I=1,3*NAT
         X(I95+I-1)=DE(I,1)
   10    DE (I,1)=ZERO
C
C     ----- STORE AWAY 1E-FORCE CONSTANT -----
C
      DO 11 I=1,NC*NC
         X(I96+I-1)=X(I00+I-1)
   11    X(I00+I-1)=ZERO     
C
C     ----- STORE AWAY 1E-DERIVATIVE FOCK MATRICES -----
C
      DO 12 I=1,NC*NUM2
         X(I97+I-1)=X(I01+I-1)
   12    X (I01+I-1)=ZERO     
C
C     ----- GET -TV- DERIVATIVES -----
C
      CALL HND_TVD_NEW(X(I00),NC,X(I01),NUM2,X(I10))
C
C     ----- ADD 1E-GRADIENTS -----
C
      DO 30 I=1,3*NAT
   30    DE(I,1)=DE(I,1)+X(I95+I-1)
C
C     ----- ADD 1E-FORCE CONSTANT -----
C
      DO 31 I=1,NC*NC
   31    X(I00+I-1)=X(I00+I-1)+X(I96+I-1)
C
C     ----- ADD 1E-DERIVATIVE FOCK MATRICES -----
C
      DO 32 I=1,NC*NUM2
   32    X(I01+I-1)=X(I01+I-1)+X(I97+I-1)
C
C     ----- PRINT OUT -----
C
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9999)
      IF(.NOT.OUT) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
         DO 9200 I=1,3
 9200    WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      ENDIF
      IF(MAX.LT.NAT) GO TO 9100
C
      NC=3*NAT
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9991)
         CALL HND_HSSOUT(FC,NC,NC,NC)
         WRITE(IW,9990)
         DO 9250 IC=1,NC
 9250    CALL HND_PRTR(X(I01+(IC-1)*NUM2),NUM)
C
         WRITE(IW,9998)
      ENDIF
C
 9300 CONTINUE
C
C     ----- WRITE -DF- ON -FT11- -----
C
      CALL HND_REWFIL(NFT17)
      DO 9400 IC=1,3*NAT
      CALL HND_HSSWRT(NFT17,X(I01+(IC-1)*NUM2),NUM2)
 9400 CONTINUE
C
      RETURN
C
 9999 FORMAT(/,10X,20(1H-),/,10X,'-TV- DERIVATIVE TERM',
     1       /,10X,20(1H-))
 9998 FORMAT(' ...... END OF -TV- DERIVATIVE TERM ......')
 9997 FORMAT(' IN -TVD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,17(1H-),/,10X,'-TV- HESSIAN TERM',
     1       /,10X,17(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      END
      SUBROUTINE HND_TVD_NEW(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,*),DF(LDF,*),DAB(*)
      DIMENSION DITIJ(3,225),DIVIJ(3,225)
      DIMENSION DIDITIJ(3,3,225),DIDJTIJ(3,3,225)
      DIMENSION DIDIVIJ(3,3,225),DIDJVIJ(3,3,225)
      DIMENSION DIJ(225)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA NFT17  /17/
      DATA ERRMSG /'PROGRAM ','STOP IN ','- TVD  -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      OUT =.FALSE.
C
C     ----- CALCULATE -TV- DERIVATIVE TERM -----
C
      IF(OUT) WRITE(IW,9999)
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
         MINI=KMIN(II)
         MAXI=KMAX(II)
         LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,NSHELL
         MINJ=KMIN(JJ)
         MAXJ=KMAX(JJ)
         LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
         IJSH=IJSH+1
         IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 8000
C
         IANDJ=II.EQ.JJ
C
C     ----- GET DENSITY BLOCK -----
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IDIM=MAX0(LOCI+I,LOCJ+J)
            JDIM=MIN0(LOCI+I,LOCJ+J)
            NN=(IDIM*(IDIM-1))/2+JDIM
            DEN=DAB(NN)
            IF(.NOT.IANDJ.OR.I.NE.J) THEN
               DEN=DEN+DEN
            ENDIF
            IJ=IJ+1
            DIJ(IJ)=DEN
         ENDDO    
      ENDDO    
C
C     ----- (II,JJ) BLOCK -----   
C
      CALL HND_TVD_IJ(II,JJ,IAT,JAT,
     1                DITIJ,DIDITIJ,DIDJTIJ,
     2                DIVIJ,DIDIVIJ,DIDJVIJ)
C
C     ----- KINETIC ENERGY DERIVATIVES -----
C
      IC=3*(IAT-1)
      JC=3*(JAT-1)
      IJ=0
      DO 390 I=MINI,MAXI
      DO 380 J=MINJ,MAXJ
      IJ=IJ+1
      DUMX =DITIJ(1,IJ)
      DUMY =DITIJ(2,IJ)
      DUMZ =DITIJ(3,IJ)
      DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
      DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
      DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
C
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DF(NN,IC+1)=DF(NN,IC+1)+ DUMX
      DF(NN,IC+2)=DF(NN,IC+2)+ DUMY
      DF(NN,IC+3)=DF(NN,IC+3)+ DUMZ
C
      DUMXX=DIDITIJ(1,1,IJ)
      DUMYY=DIDITIJ(2,2,IJ)
      DUMZZ=DIDITIJ(3,3,IJ)
      DUMYX=DIDITIJ(2,1,IJ)
      DUMZX=DIDITIJ(3,1,IJ)
      DUMZY=DIDITIJ(3,2,IJ)
      FC(IC+1,IC+1)=FC(IC+1,IC+1)+ DUMXX*DIJ(IJ)
      FC(IC+1,IC+2)=FC(IC+1,IC+2)+ DUMYX*DIJ(IJ)
      FC(IC+1,IC+3)=FC(IC+1,IC+3)+ DUMZX*DIJ(IJ)
      FC(IC+2,IC+1)=FC(IC+2,IC+1)+ DUMYX*DIJ(IJ)
      FC(IC+2,IC+2)=FC(IC+2,IC+2)+ DUMYY*DIJ(IJ)
      FC(IC+2,IC+3)=FC(IC+2,IC+3)+ DUMZY*DIJ(IJ)
      FC(IC+3,IC+1)=FC(IC+3,IC+1)+ DUMZX*DIJ(IJ)
      FC(IC+3,IC+2)=FC(IC+3,IC+2)+ DUMZY*DIJ(IJ)
      FC(IC+3,IC+3)=FC(IC+3,IC+3)+ DUMZZ*DIJ(IJ)
      DUMXX=DIDJTIJ(1,1,IJ)
      DUMXY=DIDJTIJ(1,2,IJ)
      DUMXZ=DIDJTIJ(1,3,IJ)
      DUMYX=DIDJTIJ(2,1,IJ)
      DUMYY=DIDJTIJ(2,2,IJ)
      DUMYZ=DIDJTIJ(2,3,IJ)
      DUMZX=DIDJTIJ(3,1,IJ)
      DUMZY=DIDJTIJ(3,2,IJ)
      DUMZZ=DIDJTIJ(3,3,IJ)
      FC(IC+1,JC+1)=FC(IC+1,JC+1)+ DUMXX*DIJ(IJ)
      FC(IC+1,JC+2)=FC(IC+1,JC+2)+ DUMXY*DIJ(IJ)
      FC(IC+1,JC+3)=FC(IC+1,JC+3)+ DUMXZ*DIJ(IJ)
      FC(IC+2,JC+1)=FC(IC+2,JC+1)+ DUMYX*DIJ(IJ)
      FC(IC+2,JC+2)=FC(IC+2,JC+2)+ DUMYY*DIJ(IJ)
      FC(IC+2,JC+3)=FC(IC+2,JC+3)+ DUMYZ*DIJ(IJ)
      FC(IC+3,JC+1)=FC(IC+3,JC+1)+ DUMZX*DIJ(IJ)
      FC(IC+3,JC+2)=FC(IC+3,JC+2)+ DUMZY*DIJ(IJ)
      FC(IC+3,JC+3)=FC(IC+3,JC+3)+ DUMZZ*DIJ(IJ)
  380 CONTINUE
  390 CONTINUE
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      IC=3*(IAT-1)
      JC=3*(JAT-1)
      IJ=0
      DO 450 I=MINI,MAXI
      DO 440 J=MINJ,MAXJ
      IJ=IJ+1
      DUMX =DIVIJ(1,IJ)
      DUMY =DIVIJ(2,IJ)
      DUMZ =DIVIJ(3,IJ)
      DE(1,IAT)=DE(1,IAT)+DUMX*DIJ(IJ)
      DE(2,IAT)=DE(2,IAT)+DUMY*DIJ(IJ)
      DE(3,IAT)=DE(3,IAT)+DUMZ*DIJ(IJ)
C
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DF(NN,IC+1)=DF(NN,IC+1)+DUMX
      DF(NN,IC+2)=DF(NN,IC+2)+DUMY
      DF(NN,IC+3)=DF(NN,IC+3)+DUMZ
C
      DUMXX=DIDIVIJ(1,1,IJ)
      DUMYY=DIDIVIJ(2,2,IJ)
      DUMZZ=DIDIVIJ(3,3,IJ)
      DUMYX=DIDIVIJ(2,1,IJ)
      DUMZX=DIDIVIJ(3,1,IJ)
      DUMZY=DIDIVIJ(3,2,IJ)
      FC(IC+1,IC+1)=FC(IC+1,IC+1)+DUMXX*DIJ(IJ)
      FC(IC+1,IC+2)=FC(IC+1,IC+2)+DUMYX*DIJ(IJ)
      FC(IC+1,IC+3)=FC(IC+1,IC+3)+DUMZX*DIJ(IJ)
      FC(IC+2,IC+1)=FC(IC+2,IC+1)+DUMYX*DIJ(IJ)
      FC(IC+2,IC+2)=FC(IC+2,IC+2)+DUMYY*DIJ(IJ)
      FC(IC+2,IC+3)=FC(IC+2,IC+3)+DUMZY*DIJ(IJ)
      FC(IC+3,IC+1)=FC(IC+3,IC+1)+DUMZX*DIJ(IJ)
      FC(IC+3,IC+2)=FC(IC+3,IC+2)+DUMZY*DIJ(IJ)
      FC(IC+3,IC+3)=FC(IC+3,IC+3)+DUMZZ*DIJ(IJ)
C
      DUMXX=DIDJVIJ(1,1,IJ)
      DUMXY=DIDJVIJ(1,2,IJ)
      DUMXZ=DIDJVIJ(1,3,IJ)
      DUMYX=DIDJVIJ(2,1,IJ)
      DUMYY=DIDJVIJ(2,2,IJ)
      DUMYZ=DIDJVIJ(2,3,IJ)
      DUMZX=DIDJVIJ(3,1,IJ)
      DUMZY=DIDJVIJ(3,2,IJ)
      DUMZZ=DIDJVIJ(3,3,IJ)
      FC(IC+1,JC+1)=FC(IC+1,JC+1)+DUMXX*DIJ(IJ)
      FC(IC+1,JC+2)=FC(IC+1,JC+2)+DUMXY*DIJ(IJ)
      FC(IC+1,JC+3)=FC(IC+1,JC+3)+DUMXZ*DIJ(IJ)
      FC(IC+2,JC+1)=FC(IC+2,JC+1)+DUMYX*DIJ(IJ)
      FC(IC+2,JC+2)=FC(IC+2,JC+2)+DUMYY*DIJ(IJ)
      FC(IC+2,JC+3)=FC(IC+2,JC+3)+DUMYZ*DIJ(IJ)
      FC(IC+3,JC+1)=FC(IC+3,JC+1)+DUMZX*DIJ(IJ)
      FC(IC+3,JC+2)=FC(IC+3,JC+2)+DUMZY*DIJ(IJ)
      FC(IC+3,JC+3)=FC(IC+3,JC+3)+DUMZZ*DIJ(IJ)
  440 CONTINUE
  450 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      IF(OUT) THEN       
      WRITE(IW,9993) II,JJ
      MAX=0
 7500 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 7600 I=1,3
 7600 WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 7500
      NC=3*NAT
      CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DF,(NC*NUM2),'+')
C
C     ----- PRINT OUT -----
C
      IF(.NOT.OUT) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 9200 I=1,3
 9200 WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 9100
C
      NC=3*NAT
      WRITE(IW,9991)
      CALL HND_HSSOUT(FC,NC,NC,NC)
      WRITE(IW,9990)
      DO 9250 IC=1,NC

 9250  CALL HND_PRTR(DF(1,IC),NUM)
C
      WRITE(IW,9998)
 9300 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,20(1H-),/,10X,'-TV- DERIVATIVE TERM',
     1       /,10X,20(1H-))
 9998 FORMAT(' ...... END OF -TV- DERIVATIVE TERM ......')
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,17(1H-),/,10X,'-TV- HESSIAN TERM',
     1       /,10X,17(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      END
      SUBROUTINE HND_TVD_IJ(II,JJ,IAT,JAT,
     1                      DITIJ,DIDITIJ,DIDJTIJ,
     2                      DIVIJ,DIDIVIJ,DIDJVIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      DIMENSION DITIJ(3,225),DIVIJ(3,225)
      DIMENSION DIDITIJ(3,3,225),DIDJTIJ(3,3,225)
      DIMENSION DIDIVIJ(3,3,225),DIDJVIJ(3,3,225)
      DIMENSION PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XS(7,8)  ,    YS(7,8)  ,    ZS(7,8)
      DIMENSION  DXSDI(5,5)  , DYSDI(5,5)  , DZSDI(5,5)
      DIMENSION  DXSDJ(5,5)  , DYSDJ(5,5)  , DZSDJ(5,5)
      DIMENSION DXSDII(5,5)  ,DYSDII(5,5)  ,DZSDII(5,5)
      DIMENSION DXSDIJ(5,5)  ,DYSDIJ(5,5)  ,DZSDIJ(5,5)
      DIMENSION     XT(7,6)  ,    YT(7,6)  ,    ZT(7,6)
      DIMENSION  DXTDI(5,5)  , DYTDI(5,5)  , DZTDI(5,5)
      DIMENSION  DXTDJ(5,5)  , DYTDJ(5,5)  , DZTDJ(5,5)
      DIMENSION DXTDII(5,5)  ,DYTDII(5,5)  ,DZTDII(5,5)
      DIMENSION DXTDIJ(5,5)  ,DYTDIJ(5,5)  ,DZTDIJ(5,5)
      DIMENSION     XV(7,6,5),    YV(7,6,5),    ZV(7,6,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION  DXVDJ(5,5,5), DYVDJ(5,5,5), DZVDJ(5,5,5)
      DIMENSION DXVDII(5,5,5),DYVDII(5,5,5),DZVDII(5,5,5)
      DIMENSION DXVDIJ(5,5,5),DYVDIJ(5,5,5),DZVDIJ(5,5,5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- TVD  -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
      NDER=2
C
C     ----- ISHELL -----
C
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER
      IAT   =I
C
C     ----- JSHELL -----
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTDER=LJT+NDER-1
      LJTMOD=LJTDER+2
      JAT   =J
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            DO IXYZ=1,3
               DITIJ(IXYZ,IJ)=ZERO
               DIVIJ(IXYZ,IJ)=ZERO
               DO JXYZ=1,3
                  DIDITIJ(IXYZ,JXYZ,IJ)=ZERO
                  DIDJTIJ(IXYZ,JXYZ,IJ)=ZERO
                  DIDIVIJ(IXYZ,JXYZ,IJ)=ZERO
                  DIDJVIJ(IXYZ,JXYZ,IJ)=ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
      PIJ(IJ)=DUM2
      IF(IANDJ.AND.I.EQ.J) PIJ(IJ)=PIJ(IJ)+PIJ(IJ)
  360 CONTINUE
C
C     ----- KINETIC ENERGY DERIVATIVES -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTMOD
      NJ=J
      DO 370 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XS(I,J)=XINT*T
      YS(I,J)=YINT*T
      ZS(I,J)=ZINT*T
  370 CONTINUE
C
      CALL HND_DDTXYZ(XT,YT,ZT,XS,YS,ZS,LITDER,LJTDER,AJ)
C
      CALL HND_DDERIJ(DXSDII,DYSDII,DZSDII,DXSDIJ,DYSDIJ,
     1            DZSDIJ,
     1            DXSDI ,DYSDI ,DZSDI ,DXSDJ ,DYSDJ ,DZSDJ ,
     1            XS,YS,ZS,LIT,LJT,AI,AJ)
      CALL HND_DDERIJ(DXTDII,DYTDII,DZTDII,DXTDIJ,DYTDIJ,
     1            DZTDIJ,
     1            DXTDI ,DYTDI ,DZTDI ,DXTDJ ,DYTDJ ,DZTDJ ,
     1            XT,YT,ZT,LIT,LJT,AI,AJ)
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 380 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX =DXTDI (IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMY = XT   (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMZ = XT   (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      IJ=IJ+1
      DITIJ(1,IJ)=DITIJ(1,IJ)+ DUMX*PIJ(IJ)
      DITIJ(2,IJ)=DITIJ(2,IJ)+ DUMY*PIJ(IJ)
      DITIJ(3,IJ)=DITIJ(3,IJ)+ DUMZ*PIJ(IJ)
C
      DUMXX=DXTDII(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDII(IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDII(IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMYY= XT   (IX,JX)*DYSDII(IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDII(IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDII(IY,JY)* ZT   (IZ,JZ)
      DUMZZ= XT   (IX,JX)* YS   (IY,JY)*DZSDII(IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDII(IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDII(IZ,JZ)
      DUMYX=DXTDI (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMZX=DXTDI (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      DUMZY= XT   (IX,JX)*DYSDI (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)*DZTDI (IZ,JZ)
      DIDITIJ(1,1,IJ)=DIDITIJ(1,1,IJ)+ DUMXX*PIJ(IJ)
      DIDITIJ(1,2,IJ)=DIDITIJ(1,2,IJ)+ DUMYX*PIJ(IJ)
      DIDITIJ(1,3,IJ)=DIDITIJ(1,3,IJ)+ DUMZX*PIJ(IJ)
      DIDITIJ(2,1,IJ)=DIDITIJ(2,1,IJ)+ DUMYX*PIJ(IJ)
      DIDITIJ(2,2,IJ)=DIDITIJ(2,2,IJ)+ DUMYY*PIJ(IJ)
      DIDITIJ(2,3,IJ)=DIDITIJ(2,3,IJ)+ DUMZY*PIJ(IJ)
      DIDITIJ(3,1,IJ)=DIDITIJ(3,1,IJ)+ DUMZX*PIJ(IJ)
      DIDITIJ(3,2,IJ)=DIDITIJ(3,2,IJ)+ DUMZY*PIJ(IJ)
      DIDITIJ(3,3,IJ)=DIDITIJ(3,3,IJ)+ DUMZZ*PIJ(IJ)
      DUMXX=DXTDIJ(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDIJ(IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDIJ(IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMXY=DXTDI (IX,JX)*DYSDJ (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)*DYTDJ (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)*DYSDJ (IY,JY)* ZT   (IZ,JZ)
      DUMXZ=DXTDI (IX,JX)* YS   (IY,JY)*DZSDJ (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)*DZSDJ (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)*DZTDJ (IZ,JZ)
      DUMYX=DXTDJ (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDJ (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDJ (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMYY= XT   (IX,JX)*DYSDIJ(IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDIJ(IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDIJ(IY,JY)* ZT   (IZ,JZ)
      DUMYZ= XT   (IX,JX)*DYSDI (IY,JY)*DZSDJ (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)*DZSDJ (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)*DZTDJ (IZ,JZ)
      DUMZX=DXTDJ (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     +DXSDJ (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     +DXSDJ (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      DUMZY= XT   (IX,JX)*DYSDJ (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)*DYTDJ (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)*DYSDJ (IY,JY)*DZTDI (IZ,JZ)
      DUMZZ= XT   (IX,JX)* YS   (IY,JY)*DZSDIJ(IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDIJ(IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDIJ(IZ,JZ)
      DIDJTIJ(1,1,IJ)=DIDJTIJ(1,1,IJ)+ DUMXX*PIJ(IJ)
      DIDJTIJ(1,2,IJ)=DIDJTIJ(1,2,IJ)+ DUMXY*PIJ(IJ)
      DIDJTIJ(1,3,IJ)=DIDJTIJ(1,3,IJ)+ DUMXZ*PIJ(IJ)
      DIDJTIJ(2,1,IJ)=DIDJTIJ(2,1,IJ)+ DUMYX*PIJ(IJ)
      DIDJTIJ(2,2,IJ)=DIDJTIJ(2,2,IJ)+ DUMYY*PIJ(IJ)
      DIDJTIJ(2,3,IJ)=DIDJTIJ(2,3,IJ)+ DUMYZ*PIJ(IJ)
      DIDJTIJ(3,1,IJ)=DIDJTIJ(3,1,IJ)+ DUMZX*PIJ(IJ)
      DIDJTIJ(3,2,IJ)=DIDJTIJ(3,2,IJ)+ DUMZY*PIJ(IJ)
      DIDJTIJ(3,3,IJ)=DIDJTIJ(3,3,IJ)+ DUMZZ*PIJ(IJ)
  380 CONTINUE
  390 CONTINUE
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 KAT=1,NAT
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
C
      DO 410 J=1,LJTDER
      NJ=J
      DO 410 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT*WW
  410 CONTINUE
C
      CALL HND_DDERIJ(DXVDII(1,1,IROOT),DYVDII(1,1,IROOT),
     1            DZVDII(1,1,IROOT),
     1            DXVDIJ(1,1,IROOT),DYVDIJ(1,1,IROOT),DZVDIJ(1,1,IROOT),
     2             DXVDI(1,1,IROOT), DYVDI(1,1,IROOT), DZVDI(1,1,IROOT),
     3             DXVDJ(1,1,IROOT), DYVDJ(1,1,IROOT), DZVDJ(1,1,IROOT),
     4                XV(1,1,IROOT),    YV(1,1,IROOT),    ZV(1,1,IROOT),
     5            LIT,LJT,AI,AJ)
C
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 440 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
C
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX=DUMX+DXVDI(IX,JX,IROOT)* YV  (IY,JY,IROOT)* ZV  (IZ,JZ,IROOT)
      DUMY=DUMY+ XV  (IX,JX,IROOT)*DYVDI(IY,JY,IROOT)* ZV  (IZ,JZ,IROOT)
      DUMZ=DUMZ+ XV  (IX,JX,IROOT)* YV  (IY,JY,IROOT)*DZVDI(IZ,JZ,IROOT)
  430 CONTINUE
      DUMX=DUMX*AA1*PI212
      DUMY=DUMY*AA1*PI212
      DUMZ=DUMZ*AA1*PI212
      IJ=IJ+1
      DIVIJ(1,IJ)=DIVIJ(1,IJ)+DUMX*PIJ(IJ)
      DIVIJ(2,IJ)=DIVIJ(2,IJ)+DUMY*PIJ(IJ)
      DIVIJ(3,IJ)=DIVIJ(3,IJ)+DUMZ*PIJ(IJ)
C
      DUMXX=ZERO
      DUMYY=ZERO
      DUMZZ=ZERO
      DUMYX=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DO 432 IROOT=1,NROOTS
      DUMXX=DUMXX
     1     +DXVDII(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYY=DUMYY
     1     + XV   (IX,JX,IROOT)*DYVDII(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMZZ=DUMZZ
     1     + XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDII(IZ,JZ,IROOT)
      DUMYX=DUMYX
     1     +DXVDI (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMZX=DUMZX
     1     +DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZY=DUMZY
     1     + XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
  432 CONTINUE
      DUMXX=DUMXX*AA1*PI212
      DUMYY=DUMYY*AA1*PI212
      DUMZZ=DUMZZ*AA1*PI212
      DUMYX=DUMYX*AA1*PI212
      DUMZX=DUMZX*AA1*PI212
      DUMZY=DUMZY*AA1*PI212
      DIDIVIJ(1,1,IJ)=DIDIVIJ(1,1,IJ)+DUMXX*PIJ(IJ)
      DIDIVIJ(1,2,IJ)=DIDIVIJ(1,2,IJ)+DUMYX*PIJ(IJ)
      DIDIVIJ(1,3,IJ)=DIDIVIJ(1,3,IJ)+DUMZX*PIJ(IJ)
      DIDIVIJ(2,1,IJ)=DIDIVIJ(2,1,IJ)+DUMYX*PIJ(IJ)
      DIDIVIJ(2,2,IJ)=DIDIVIJ(2,2,IJ)+DUMYY*PIJ(IJ)
      DIDIVIJ(2,3,IJ)=DIDIVIJ(2,3,IJ)+DUMZY*PIJ(IJ)
      DIDIVIJ(3,1,IJ)=DIDIVIJ(3,1,IJ)+DUMZX*PIJ(IJ)
      DIDIVIJ(3,2,IJ)=DIDIVIJ(3,2,IJ)+DUMZY*PIJ(IJ)
      DIDIVIJ(3,3,IJ)=DIDIVIJ(3,3,IJ)+DUMZZ*PIJ(IJ)
C
      DUMXX=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYX=ZERO
      DUMYY=ZERO
      DUMYZ=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DUMZZ=ZERO
      DO 434 IROOT=1,NROOTS
      DUMXX=DUMXX
     1     +DXVDIJ(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMXY=DUMXY
     1     +DXVDI (IX,JX,IROOT)*DYVDJ (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMXZ=DUMXZ
     1     +DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDJ (IZ,JZ,IROOT)
      DUMYX=DUMYX
     1     +DXVDJ (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYY=DUMYY
     1     + XV   (IX,JX,IROOT)*DYVDIJ(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYZ=DUMYZ
     1     + XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZVDJ (IZ,JZ,IROOT)
      DUMZX=DUMZX
     1     +DXVDJ (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZY=DUMZY
     1     + XV   (IX,JX,IROOT)*DYVDJ (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZZ=DUMZZ
     1     + XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDIJ(IZ,JZ,IROOT)
  434 CONTINUE
      DUMXX=DUMXX*AA1*PI212
      DUMXY=DUMXY*AA1*PI212
      DUMXZ=DUMXZ*AA1*PI212
      DUMYX=DUMYX*AA1*PI212
      DUMYY=DUMYY*AA1*PI212
      DUMYZ=DUMYZ*AA1*PI212
      DUMZX=DUMZX*AA1*PI212
      DUMZY=DUMZY*AA1*PI212
      DUMZZ=DUMZZ*AA1*PI212
      DIDJVIJ(1,1,IJ)=DIDJVIJ(1,1,IJ)+DUMXX*PIJ(IJ) 
      DIDJVIJ(1,2,IJ)=DIDJVIJ(1,2,IJ)+DUMXY*PIJ(IJ) 
      DIDJVIJ(1,3,IJ)=DIDJVIJ(1,3,IJ)+DUMXZ*PIJ(IJ) 
      DIDJVIJ(2,1,IJ)=DIDJVIJ(2,1,IJ)+DUMYX*PIJ(IJ) 
      DIDJVIJ(2,2,IJ)=DIDJVIJ(2,2,IJ)+DUMYY*PIJ(IJ) 
      DIDJVIJ(2,3,IJ)=DIDJVIJ(2,3,IJ)+DUMYZ*PIJ(IJ) 
      DIDJVIJ(3,1,IJ)=DIDJVIJ(3,1,IJ)+DUMZX*PIJ(IJ) 
      DIDJVIJ(3,2,IJ)=DIDJVIJ(3,2,IJ)+DUMZY*PIJ(IJ) 
      DIDJVIJ(3,3,IJ)=DIDJVIJ(3,3,IJ)+DUMZZ*PIJ(IJ) 

  440 CONTINUE
  450 CONTINUE
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      RETURN
 9997 FORMAT(' IN -TVD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
      END
      SUBROUTINE HND_TVD_OLD(FC,LFC,DF,LDF,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,1),DF(LDF,1),DAB(1)
      DIMENSION DIJ(225),PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XS(7,8)  ,    YS(7,8)  ,    ZS(7,8)
      DIMENSION  DXSDI(5,5)  , DYSDI(5,5)  , DZSDI(5,5)
      DIMENSION  DXSDJ(5,5)  , DYSDJ(5,5)  , DZSDJ(5,5)
      DIMENSION DXSDII(5,5)  ,DYSDII(5,5)  ,DZSDII(5,5)
      DIMENSION DXSDIJ(5,5)  ,DYSDIJ(5,5)  ,DZSDIJ(5,5)
      DIMENSION     XT(7,6)  ,    YT(7,6)  ,    ZT(7,6)
      DIMENSION  DXTDI(5,5)  , DYTDI(5,5)  , DZTDI(5,5)
      DIMENSION  DXTDJ(5,5)  , DYTDJ(5,5)  , DZTDJ(5,5)
      DIMENSION DXTDII(5,5)  ,DYTDII(5,5)  ,DZTDII(5,5)
      DIMENSION DXTDIJ(5,5)  ,DYTDIJ(5,5)  ,DZTDIJ(5,5)
      DIMENSION     XV(7,6,5),    YV(7,6,5),    ZV(7,6,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION  DXVDJ(5,5,5), DYVDJ(5,5,5), DZVDJ(5,5,5)
      DIMENSION DXVDII(5,5,5),DYVDII(5,5,5),DZVDII(5,5,5)
      DIMENSION DXVDIJ(5,5,5),DYVDIJ(5,5,5),DZVDIJ(5,5,5)
      DIMENSION DNAM(3)
      DIMENSION ERRMSG(3)
      DATA NFT17  /17/
      DATA ERRMSG /'PROGRAM ','STOP IN ','- TVD  -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- CALCULATE -TV- DERIVATIVE TERM -----
C
      IF(OUT) WRITE(IW,9999)
      NDER=2
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER
      IAT   =I
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
      IJSH=IJSH+1
      IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 7900
C
      LJTDER=LJT+NDER-1
      LJTMOD=LJTDER+2
      JAT   =J
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
cd    NN=IA(MAX0(LOCI+I,LOCJ+J))
cd   1  +   MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
      DEN=DAB(NN)
      DEN=DEN+DEN
      IJ=IJ+1
      DIJ(IJ)=DUM2*DEN
      PIJ(IJ)=DUM2
      IF(IANDJ.AND.I.EQ.J) PIJ(IJ)=PIJ(IJ)+PIJ(IJ)

  360 CONTINUE
C
C     ----- KINETIC ENERGY DERIVATIVES -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTMOD
      NJ=J
      DO 370 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XS(I,J)=XINT*T
      YS(I,J)=YINT*T
      ZS(I,J)=ZINT*T
  370 CONTINUE
C
      CALL HND_DDTXYZ(XT,YT,ZT,XS,YS,ZS,LITDER,LJTDER,AJ)
C
      CALL HND_DDERIJ(DXSDII,DYSDII,DZSDII,DXSDIJ,DYSDIJ,
     1            DZSDIJ,
     1            DXSDI ,DYSDI ,DZSDI ,DXSDJ ,DYSDJ ,DZSDJ ,
     1            XS,YS,ZS,LIT,LJT,AI,AJ)
      CALL HND_DDERIJ(DXTDII,DYTDII,DZTDII,DXTDIJ,DYTDIJ,
     1            DZTDIJ,
     1            DXTDI ,DYTDI ,DZTDI ,DXTDJ ,DYTDJ ,DZTDJ ,
     1            XT,YT,ZT,LIT,LJT,AI,AJ)
C
      IC=3*(IAT-1)
      JC=3*(JAT-1)
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 380 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX =DXTDI (IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMY = XT   (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMZ = XT   (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      IJ=IJ+1
      DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
      DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
      DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
C
cd    NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DF(NN,IC+1)=DF(NN,IC+1)+ DUMX*PIJ(IJ)
      DF(NN,IC+2)=DF(NN,IC+2)+ DUMY*PIJ(IJ)
      DF(NN,IC+3)=DF(NN,IC+3)+ DUMZ*PIJ(IJ)
      DUMXX=DXTDII(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDII(IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDII(IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMYY= XT   (IX,JX)*DYSDII(IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDII(IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDII(IY,JY)* ZT   (IZ,JZ)
      DUMZZ= XT   (IX,JX)* YS   (IY,JY)*DZSDII(IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDII(IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDII(IZ,JZ)
      DUMYX=DXTDI (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMZX=DXTDI (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      DUMZY= XT   (IX,JX)*DYSDI (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)*DZTDI (IZ,JZ)
      FC(IC+1,IC+1)=FC(IC+1,IC+1)+ DUMXX*DIJ(IJ)
      FC(IC+1,IC+2)=FC(IC+1,IC+2)+ DUMYX*DIJ(IJ)
      FC(IC+1,IC+3)=FC(IC+1,IC+3)+ DUMZX*DIJ(IJ)
      FC(IC+2,IC+1)=FC(IC+2,IC+1)+ DUMYX*DIJ(IJ)
      FC(IC+2,IC+2)=FC(IC+2,IC+2)+ DUMYY*DIJ(IJ)
      FC(IC+2,IC+3)=FC(IC+2,IC+3)+ DUMZY*DIJ(IJ)
      FC(IC+3,IC+1)=FC(IC+3,IC+1)+ DUMZX*DIJ(IJ)
      FC(IC+3,IC+2)=FC(IC+3,IC+2)+ DUMZY*DIJ(IJ)
      FC(IC+3,IC+3)=FC(IC+3,IC+3)+ DUMZZ*DIJ(IJ)
      DUMXX=DXTDIJ(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDIJ(IX,JX)* YT   (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDIJ(IX,JX)* YS   (IY,JY)* ZT   (IZ,JZ)
      DUMXY=DXTDI (IX,JX)*DYSDJ (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDI (IX,JX)*DYTDJ (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDI (IX,JX)*DYSDJ (IY,JY)* ZT   (IZ,JZ)
      DUMXZ=DXTDI (IX,JX)* YS   (IY,JY)*DZSDJ (IZ,JZ)
     1     +DXSDI (IX,JX)* YT   (IY,JY)*DZSDJ (IZ,JZ)
     2     +DXSDI (IX,JX)* YS   (IY,JY)*DZTDJ (IZ,JZ)
      DUMYX=DXTDJ (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
     1     +DXSDJ (IX,JX)*DYTDI (IY,JY)* ZS   (IZ,JZ)
     2     +DXSDJ (IX,JX)*DYSDI (IY,JY)* ZT   (IZ,JZ)
      DUMYY= XT   (IX,JX)*DYSDIJ(IY,JY)* ZS   (IZ,JZ)
     1     + XS   (IX,JX)*DYTDIJ(IY,JY)* ZS   (IZ,JZ)
     2     + XS   (IX,JX)*DYSDIJ(IY,JY)* ZT   (IZ,JZ)
      DUMYZ= XT   (IX,JX)*DYSDI (IY,JY)*DZSDJ (IZ,JZ)
     1     + XS   (IX,JX)*DYTDI (IY,JY)*DZSDJ (IZ,JZ)
     2     + XS   (IX,JX)*DYSDI (IY,JY)*DZTDJ (IZ,JZ)
      DUMZX=DXTDJ (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
     1     +DXSDJ (IX,JX)* YT   (IY,JY)*DZSDI (IZ,JZ)
     2     +DXSDJ (IX,JX)* YS   (IY,JY)*DZTDI (IZ,JZ)
      DUMZY= XT   (IX,JX)*DYSDJ (IY,JY)*DZSDI (IZ,JZ)
     1     + XS   (IX,JX)*DYTDJ (IY,JY)*DZSDI (IZ,JZ)
     2     + XS   (IX,JX)*DYSDJ (IY,JY)*DZTDI (IZ,JZ)
      DUMZZ= XT   (IX,JX)* YS   (IY,JY)*DZSDIJ(IZ,JZ)
     1     + XS   (IX,JX)* YT   (IY,JY)*DZSDIJ(IZ,JZ)
     2     + XS   (IX,JX)* YS   (IY,JY)*DZTDIJ(IZ,JZ)
      FC(IC+1,JC+1)=FC(IC+1,JC+1)+ DUMXX*DIJ(IJ)
      FC(IC+1,JC+2)=FC(IC+1,JC+2)+ DUMXY*DIJ(IJ)
      FC(IC+1,JC+3)=FC(IC+1,JC+3)+ DUMXZ*DIJ(IJ)
      FC(IC+2,JC+1)=FC(IC+2,JC+1)+ DUMYX*DIJ(IJ)
      FC(IC+2,JC+2)=FC(IC+2,JC+2)+ DUMYY*DIJ(IJ)
      FC(IC+2,JC+3)=FC(IC+2,JC+3)+ DUMYZ*DIJ(IJ)
      FC(IC+3,JC+1)=FC(IC+3,JC+1)+ DUMZX*DIJ(IJ)
      FC(IC+3,JC+2)=FC(IC+3,JC+2)+ DUMZY*DIJ(IJ)
      FC(IC+3,JC+3)=FC(IC+3,JC+3)+ DUMZZ*DIJ(IJ)
  380 CONTINUE
  390 CONTINUE
C
C     ----- NUCLEAR ATTRACTION DERIVATIVES -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 KAT=1,NAT
      ZNUC=-ZAN(KAT)
      CX=C(1,KAT)
      CY=C(2,KAT)
      CZ=C(3,KAT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
C
      DO 410 J=1,LJTDER
      NJ=J
      DO 410 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT*WW
  410 CONTINUE
C
      CALL HND_DDERIJ(DXVDII(1,1,IROOT),DYVDII(1,1,IROOT),
     1            DZVDII(1,1,IROOT),
     1            DXVDIJ(1,1,IROOT),DYVDIJ(1,1,IROOT),DZVDIJ(1,1,IROOT),
     2             DXVDI(1,1,IROOT), DYVDI(1,1,IROOT), DZVDI(1,1,IROOT),
     3             DXVDJ(1,1,IROOT), DYVDJ(1,1,IROOT), DZVDJ(1,1,IROOT),
     4                XV(1,1,IROOT),    YV(1,1,IROOT),    ZV(1,1,IROOT),
     5            LIT,LJT,AI,AJ)
C
  420 CONTINUE
C
      IC=3*(IAT-1)
      JC=3*(JAT-1)
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 440 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX=DUMX+DXVDI(IX,JX,IROOT)* YV  (IY,JY,IROOT)* ZV  (IZ,JZ,IROOT)
      DUMY=DUMY+ XV  (IX,JX,IROOT)*DYVDI(IY,JY,IROOT)* ZV  (IZ,JZ,IROOT)
      DUMZ=DUMZ+ XV  (IX,JX,IROOT)* YV  (IY,JY,IROOT)*DZVDI(IZ,JZ,IROOT)
  430 CONTINUE
      IJ=IJ+1
      DE(1,IAT)=DE(1,IAT)+DUMX*(AA1*PI212*DIJ(IJ))
      DE(2,IAT)=DE(2,IAT)+DUMY*(AA1*PI212*DIJ(IJ))
      DE(3,IAT)=DE(3,IAT)+DUMZ*(AA1*PI212*DIJ(IJ))
C
cd    NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DF(NN,IC+1)=DF(NN,IC+1)+DUMX*(AA1*PI212*PIJ(IJ))
      DF(NN,IC+2)=DF(NN,IC+2)+DUMY*(AA1*PI212*PIJ(IJ))
      DF(NN,IC+3)=DF(NN,IC+3)+DUMZ*(AA1*PI212*PIJ(IJ))
C
      DUMXX=ZERO
      DUMYY=ZERO
      DUMZZ=ZERO
      DUMYX=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DO 432 IROOT=1,NROOTS
      DUMXX=DUMXX
     1     +DXVDII(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYY=DUMYY
     1     + XV   (IX,JX,IROOT)*DYVDII(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMZZ=DUMZZ
     1     + XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDII(IZ,JZ,IROOT)
      DUMYX=DUMYX
     1     +DXVDI (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMZX=DUMZX
     1     +DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZY=DUMZY
     1     + XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
  432 CONTINUE
      FC(IC+1,IC+1)=FC(IC+1,IC+1)+DUMXX*(AA1*PI212*DIJ(IJ))
      FC(IC+1,IC+2)=FC(IC+1,IC+2)+DUMYX*(AA1*PI212*DIJ(IJ))
      FC(IC+1,IC+3)=FC(IC+1,IC+3)+DUMZX*(AA1*PI212*DIJ(IJ))
      FC(IC+2,IC+1)=FC(IC+2,IC+1)+DUMYX*(AA1*PI212*DIJ(IJ))
      FC(IC+2,IC+2)=FC(IC+2,IC+2)+DUMYY*(AA1*PI212*DIJ(IJ))
      FC(IC+2,IC+3)=FC(IC+2,IC+3)+DUMZY*(AA1*PI212*DIJ(IJ))
      FC(IC+3,IC+1)=FC(IC+3,IC+1)+DUMZX*(AA1*PI212*DIJ(IJ))
      FC(IC+3,IC+2)=FC(IC+3,IC+2)+DUMZY*(AA1*PI212*DIJ(IJ))
      FC(IC+3,IC+3)=FC(IC+3,IC+3)+DUMZZ*(AA1*PI212*DIJ(IJ))
C
      DUMXX=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYX=ZERO
      DUMYY=ZERO
      DUMYZ=ZERO
      DUMZX=ZERO
      DUMZY=ZERO
      DUMZZ=ZERO
      DO 434 IROOT=1,NROOTS
      DUMXX=DUMXX
     1     +DXVDIJ(IX,JX,IROOT)* YV   (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMXY=DUMXY
     1     +DXVDI (IX,JX,IROOT)*DYVDJ (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMXZ=DUMXZ
     1     +DXVDI (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDJ (IZ,JZ,IROOT)
      DUMYX=DUMYX
     1     +DXVDJ (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYY=DUMYY
     1     + XV   (IX,JX,IROOT)*DYVDIJ(IY,JY,IROOT)* ZV   (IZ,JZ,IROOT)
      DUMYZ=DUMYZ
     1     + XV   (IX,JX,IROOT)*DYVDI (IY,JY,IROOT)*DZVDJ (IZ,JZ,IROOT)
      DUMZX=DUMZX
     1     +DXVDJ (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZY=DUMZY
     1     + XV   (IX,JX,IROOT)*DYVDJ (IY,JY,IROOT)*DZVDI (IZ,JZ,IROOT)
      DUMZZ=DUMZZ
     1     + XV   (IX,JX,IROOT)* YV   (IY,JY,IROOT)*DZVDIJ(IZ,JZ,IROOT)
  434 CONTINUE
      FC(IC+1,JC+1)=FC(IC+1,JC+1)+DUMXX*(AA1*PI212*DIJ(IJ))
      FC(IC+1,JC+2)=FC(IC+1,JC+2)+DUMXY*(AA1*PI212*DIJ(IJ))
      FC(IC+1,JC+3)=FC(IC+1,JC+3)+DUMXZ*(AA1*PI212*DIJ(IJ))
      FC(IC+2,JC+1)=FC(IC+2,JC+1)+DUMYX*(AA1*PI212*DIJ(IJ))
      FC(IC+2,JC+2)=FC(IC+2,JC+2)+DUMYY*(AA1*PI212*DIJ(IJ))
      FC(IC+2,JC+3)=FC(IC+2,JC+3)+DUMYZ*(AA1*PI212*DIJ(IJ))
      FC(IC+3,JC+1)=FC(IC+3,JC+1)+DUMZX*(AA1*PI212*DIJ(IJ))
      FC(IC+3,JC+2)=FC(IC+3,JC+2)+DUMZY*(AA1*PI212*DIJ(IJ))
      FC(IC+3,JC+3)=FC(IC+3,JC+3)+DUMZZ*(AA1*PI212*DIJ(IJ))

  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      IF(OUT) THEN        
      WRITE(IW,9993) II,JJ
      MAX=0
 7500 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 7600 I=1,3
 7600 WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 7500
      NC=3*NAT
      CALL HND_HSSOUT(FC,NC,NC,NC)
      ENDIF
C
 7900 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DF,(NC*NUM2),'+')
C
C     ----- PRINT OUT -----
C
*     IF(.NOT.OUT) GO TO 9300
*     MAX=0
*9100 MIN=MAX+1
*     MAX=MAX+8
*     IF(MAX.GT.NAT) MAX=NAT
*     WRITE(IW,9996)
*     WRITE(IW,9995) (I,I=MIN,MAX)
*     WRITE(IW,9996)
*     DO 9200 I=1,3
*9200 WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
*     IF(MAX.LT.NAT) GO TO 9100
C
*     NC=3*NAT
*     WRITE(IW,9991)
*     CALL HND_HSSOUT(FC,NC,NC,NC)
*     WRITE(IW,9990)
*     DO 9250 IC=1,NC

*9250  CALL HND_PRTR(DF(1,IC),NUM)
C9250  CALL HND_PRTRL(DF(1,IC),NUM)
C
*     WRITE(IW,9998)
*9300 CONTINUE
C
C     ----- WRITE -DF- ON -FT11- -----
C
C     CALL HND_REWFIL(NFT17)
C     DO 9400 IC=1,3*NAT
C     CALL HND_HSSWRT(NFT17,DF(1,IC),LDF)
C9400 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,20(1H-),/,10X,'-TV- DERIVATIVE TERM',
     1       /,10X,20(1H-))
 9998 FORMAT(' ...... END OF -TV- DERIVATIVE TERM ......')
 9997 FORMAT(' IN -TVD- , THE RYS QUADRATURE IS NOT IMPLEMEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,17(1H-),/,10X,'-TV- HESSIAN TERM',
     1       /,10X,17(1H-))
 9990 FORMAT(/,10X,22(1H-),/,10X,'DERIVATIVE FOCK MATRIX',
     1       /,10X,22(1H-))
      END
C
      SUBROUTINE HND_SDX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
      PARAMETER (MXATOM=500)
      CHARACTER*8 ERRMSG
      CHARACTER*8 SCF,RHF,GVB
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL OUT,SOME
      DIMENSION ERRMSG(3)
      DIMENSION DNAM(3)
      DATA ERRMSG      /'PROGRAM ','STOP IN ','- SDX  -'/
      DATA SCF,RHF,GVB /'SCF     ','RHF     ','GVB     '/
      DATA ZERO        /0.0D+00/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA NFT17  /17/
      DATA NFT18  /18/
C
      OUT = .FALSE.
      SOME= .FALSE.
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
      I00=INIT
      I01=I00+ NC  *NC
      I10=I01+ NUM2*NC
      I20=I10+ NUM2
      I30=I20+ NUM2
C
      I40=I30+ NUM2
      I50=I40+ NUM2
      I60=I50+ NUM2
C     I31=I20+ NUM* NUM
C
c*    I70=I60+ NUM*NUM
c*    I80=I70+ NUM*NUM
c*    I90=I80+ NUM*NUM
c*    I95=I90+ NUM*NUM
      I70=I60+ 0
      I80=I70+ 0
      I90=I80+ 0
      I95=I90+ 0
C
      I96=I95+ NC
      I97=I96+ NC*NC
      I98=I97+ NUM2*NC
C
C     ----- STORE AWAY 1E-GRADIENT -----
C
      DO 10 I=1,3*NAT
         X(I95+I-1)=DE(I,1)
   10    DE (I,1)=ZERO
C
C     ----- STORE AWAY 1E-FORCE CONSTANT -----
C
      DO 11 I=1,NC*NC
         X(I96+I-1)=X(I00+I-1)
   11    X(I00+I-1)=ZERO     
C
C
C     ----- INITIALIZE OVERLAP DERIVATIVE MATRICES -----
C
      DO 20 I=1,NUM2*NC
   20 X(I+I01-1)=ZERO
C
C     ----- GET -S- DERIVATIVES -----
C
      CALL HND_SD_NEW(X(I00),NC,X(I01),NUM2,X(I10))
C
C
C     ----- ADD 1E-GRADIENTS -----
C
      DO 30 I=1,3*NAT
   30    DE(I,1)=DE(I,1)+X(I95+I-1)
C
C     ----- ADD 1E-FORCE CONSTANT -----
C
      DO 31 I=1,NC*NC
   31    X(I00+I-1)=X(I00+I-1)+X(I96+I-1)
C
C     ----- WRITE -DS- ON -NFT18- AND -FC- ON -NFT17- -----
C
      CALL HND_REWFIL(NFT18)
      DO 9400 IC=1,3*NAT
      CALL HND_HSSWRT(NFT18,X(I01+(IC-1)*NUM2),NUM2)
 9400 CONTINUE
      CALL HND_HSSWRT(NFT17,X(I00),NC*NC)
C
      IF(.NOT.OUT) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      IF (GA_NODEID().EQ.0) THEN
         WRITE(IW,9989)
         WRITE(IW,9996)
         WRITE(IW,9995) (I,I=MIN,MAX)
         WRITE(IW,9996)
         DO 9200 I=1,3
 9200    WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      ENDIF

      IF(MAX.LT.NAT) GO TO 9100
C
      IF(.NOT.SOME) GO TO 9300
C     NC=3*NAT
      WRITE(IW,9991)
      CALL HND_HSSOUT(X(I00),NC,NC,NC)
      WRITE(IW,9990)
      DO 9250 IC=1,3*NAT
         CALL HND_PRTR(X(I01+(IC-1)*NUM2),NUM)
C        CALL HND_PRTR(DS(1,IC),NUM)
C        CALL HND_PRTRL(DS(1,IC),NUM)
 9250 CONTINUE
C
      WRITE(IW,9998)
 9300 CONTINUE
C
 9999 FORMAT(/,10X,19(1H-),/,10X,'-S- DERIVATIVE TERM',
     1       /,10X,19(1H-))
 9998 FORMAT(' ...... END OF -S- DERIVATIVE TERM ......')
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,16(1H-),/,10X,'-S- HESSIAN TERM',
     1       /,10X,16(1H-))
 9990 FORMAT(/,10X,21(1H-),/,10X,'-S- DERIVATIVE MATRIX',
     1       /,10X,21(1H-))
 9989 FORMAT(/,' 1E-GRADIENT MATRIX ',/,1X,24(1H-))
C
      RETURN
      END
      SUBROUTINE HND_SD_NEW(FC,LFC,DS,LDS,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,*),DS(LDS,*),DAB(*)
      DIMENSION DIDISIJ(3,3,225),DIDJSIJ(3,3,225),DISIJ(3,225)
      DIMENSION DIJ(225)
      DIMENSION DNAM(3)
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      OUT =.FALSE.
C
C     ----- CALCULATE -S- DERIVATIVE TERM -----
C
      IF(OUT) WRITE(IW,9999)
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
         MINI=KMIN(II)
         MAXI=KMAX(II)
         LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
         DO 8000 JJ=1,NSHELL
            MINJ=KMIN(JJ)
            MAXJ=KMAX(JJ)
            LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
            IJSH=IJSH+1
            IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 8000
C
            IANDJ=II.EQ.JJ
            IF(IANDJ) GO TO 8000
C
C     ----- GET DENSITY BLOCK ; NOTE: DEN FROM NWCHEM = - DEN FROM HONDO
C
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ
                  IDIM=MAX0(LOCI+I,LOCJ+J)
                  JDIM=MIN0(LOCI+I,LOCJ+J)
                  NN  =(IDIM*(IDIM-1))/2+JDIM
                  DEN = DAB(NN)
                  DEN = DEN+DEN
                  DEN =-DEN
                  IJ=IJ+1
                  DIJ(IJ)=DEN
               ENDDO
            ENDDO
C
C     ----- (II,JJ) BLOCK -----
C
            CALL HND_SD_IJ(II,JJ,IAT,JAT,DISIJ,DIDISIJ,DIDJSIJ)
C
            IC=3*(IAT-1)
            JC=3*(JAT-1)
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ
                  IJ=IJ+1
                  DUMX =DISIJ(1,IJ)
                  DUMY =DISIJ(2,IJ)
                  DUMZ =DISIJ(3,IJ)
                  DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
                  DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
                  DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
C
                  IDIM=MAX0(LOCI+I,LOCJ+J)
                  JDIM=MIN0(LOCI+I,LOCJ+J)
                  NN=(IDIM*(IDIM-1))/2+JDIM
C
                  DS(NN,IC+1)=DS(NN,IC+1)+ DUMX
                  DS(NN,IC+2)=DS(NN,IC+2)+ DUMY
                  DS(NN,IC+3)=DS(NN,IC+3)+ DUMZ
C
                  DUMXX=DIDISIJ(1,1,IJ)
                  DUMYY=DIDISIJ(2,2,IJ)
                  DUMZZ=DIDISIJ(3,3,IJ)
                  DUMXY=DIDISIJ(1,2,IJ)
                  DUMXZ=DIDISIJ(1,3,IJ)
                  DUMYZ=DIDISIJ(2,3,IJ)
                  FC(IC+1,IC+1)=FC(IC+1,IC+1)+ DUMXX*DIJ(IJ) 
                  FC(IC+1,IC+2)=FC(IC+1,IC+2)+ DUMXY*DIJ(IJ) 
                  FC(IC+1,IC+3)=FC(IC+1,IC+3)+ DUMXZ*DIJ(IJ) 
                  FC(IC+2,IC+1)=FC(IC+2,IC+1)+ DUMXY*DIJ(IJ) 
                  FC(IC+2,IC+2)=FC(IC+2,IC+2)+ DUMYY*DIJ(IJ) 
                  FC(IC+2,IC+3)=FC(IC+2,IC+3)+ DUMYZ*DIJ(IJ) 
                  FC(IC+3,IC+1)=FC(IC+3,IC+1)+ DUMXZ*DIJ(IJ) 
                  FC(IC+3,IC+2)=FC(IC+3,IC+2)+ DUMYZ*DIJ(IJ) 
                  FC(IC+3,IC+3)=FC(IC+3,IC+3)+ DUMZZ*DIJ(IJ) 
C
                  DUMXX=DIDJSIJ(1,1,IJ)
                  DUMXY=DIDJSIJ(1,2,IJ)
                  DUMXZ=DIDJSIJ(1,3,IJ)
                  DUMYX=DIDJSIJ(2,1,IJ)
                  DUMYY=DIDJSIJ(2,2,IJ)
                  DUMYZ=DIDJSIJ(2,3,IJ)
                  DUMZX=DIDJSIJ(3,1,IJ)
                  DUMZY=DIDJSIJ(3,2,IJ)
                  DUMZZ=DIDJSIJ(3,3,IJ)
                  FC(IC+1,JC+1)=FC(IC+1,JC+1)+ DUMXX*DIJ(IJ) 
                  FC(IC+1,JC+2)=FC(IC+1,JC+2)+ DUMXY*DIJ(IJ) 
                  FC(IC+1,JC+3)=FC(IC+1,JC+3)+ DUMXZ*DIJ(IJ) 
                  FC(IC+2,JC+1)=FC(IC+2,JC+1)+ DUMYX*DIJ(IJ) 
                  FC(IC+2,JC+2)=FC(IC+2,JC+2)+ DUMYY*DIJ(IJ) 
                  FC(IC+2,JC+3)=FC(IC+2,JC+3)+ DUMYZ*DIJ(IJ) 
                  FC(IC+3,JC+1)=FC(IC+3,JC+1)+ DUMZX*DIJ(IJ) 
                  FC(IC+3,JC+2)=FC(IC+3,JC+2)+ DUMZY*DIJ(IJ) 
                  FC(IC+3,JC+3)=FC(IC+3,JC+3)+ DUMZZ*DIJ(IJ) 
               ENDDO
            ENDDO
C
            IF(OUT) THEN           
               WRITE(IW,9993) II,JJ
               MAX=0
 7500          MIN=MAX+1
               MAX=MAX+8
               IF(MAX.GT.NAT) MAX=NAT
               WRITE(IW,9996)
               WRITE(IW,9995) (I,I=MIN,MAX)
               WRITE(IW,9996)
               DO I=1,3
                  WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
               ENDDO
               IF(MAX.LT.NAT) GO TO 7500
            ENDIF
 8000    CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DS DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DS,(NC*NUM2),'+')
C
      IF(.NOT.OUT) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 9200 I=1,3
 9200 WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 9100
C
      NC=3*NAT
      WRITE(IW,9991)
      CALL HND_HSSOUT(FC,NC,NC,NC)
      WRITE(IW,9990)
      DO 9250 IC=1,3*NAT
         CALL HND_PRTR(DS(1,IC),NUM)
 9250 CONTINUE
C
      WRITE(IW,9998)
 9300 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,19(1H-),/,10X,'-S- DERIVATIVE TERM',
     1       /,10X,19(1H-))
 9998 FORMAT(' ...... END OF -S- DERIVATIVE TERM ......')
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,16(1H-),/,10X,'-S- HESSIAN TERM',
     1       /,10X,16(1H-))
 9990 FORMAT(/,10X,21(1H-),/,10X,'-S- DERIVATIVE MATRIX',
     1       /,10X,21(1H-))
      END
      SUBROUTINE HND_SD_IJ(II,JJ,IAT,JAT,DISIJ,DIDISIJ,DIDJSIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      DIMENSION DIDISIJ(3,3,225),DIDJSIJ(3,3,225),DISIJ(3,225),PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION  XS   (7,6), YS   (7,6), ZS   (7,6)
      DIMENSION DXSDI (5,5),DYSDI (5,5),DZSDI (5,5)
      DIMENSION DXSDJ (5,5),DYSDJ (5,5),DZSDJ (5,5)
      DIMENSION DXSDII(5,5),DYSDII(5,5),DZSDII(5,5)
      DIMENSION DXSDIJ(5,5),DYSDIJ(5,5),DZSDIJ(5,5)
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT =.FALSE.
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      TOL =RLN10*ITOL
C
      NDER=2
C
C     ----- ISHELL -----
C
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER
      IAT   =I
C
C     ----- JSHELL -----
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTDER=LJT+NDER-1
      JAT=J
C
      IANDJ=II.EQ.JJ
      IF(.NOT.IANDJ) THEN       
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            DO IXYZ=1,3
               DISIJ(IXYZ,IJ)=ZERO
               DO JXYZ=1,3
                  DIDISIJ(IXYZ,JXYZ,IJ)=ZERO
                  DIDJSIJ(IXYZ,JXYZ,IJ)=ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
      PIJ(IJ)=DUM2
      IF(IANDJ.AND.I.EQ.J) PIJ(IJ)=PIJ(IJ)+PIJ(IJ)
  360 CONTINUE
C
C     ----- OVERLAP DERIVATIVES -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTDER
      NJ=J
      DO 370 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XS(I,J)=XINT*T
      YS(I,J)=YINT*T
      ZS(I,J)=ZINT*T
  370 CONTINUE
C
      CALL HND_DDERIJ(DXSDII,DYSDII,DZSDII,DXSDIJ,DYSDIJ,DZSDIJ,
     1            DXSDI ,DYSDI ,DZSDI ,DXSDJ ,DYSDJ ,DZSDJ ,
     2            XS,YS,ZS,LIT,LJT,AI,AJ)
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 380 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
C
      IJ=IJ+1
C
      DUMX =DXSDI (IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMY = XS   (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMZ = XS   (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      DISIJ(1,IJ)=DISIJ(1,IJ)+ DUMX*PIJ(IJ)
      DISIJ(2,IJ)=DISIJ(2,IJ)+ DUMY*PIJ(IJ)
      DISIJ(3,IJ)=DISIJ(3,IJ)+ DUMZ*PIJ(IJ)
C
      DUMXX=DXSDII(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMYY= XS   (IX,JX)*DYSDII(IY,JY)* ZS   (IZ,JZ)
      DUMZZ= XS   (IX,JX)* YS   (IY,JY)*DZSDII(IZ,JZ)
      DUMXY=DXSDI (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMXZ=DXSDI (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      DUMYZ= XS   (IX,JX)*DYSDI (IY,JY)*DZSDI (IZ,JZ)
      DIDISIJ(1,1,IJ)=DIDISIJ(1,1,IJ)+ DUMXX*PIJ(IJ) 
      DIDISIJ(1,2,IJ)=DIDISIJ(1,2,IJ)+ DUMXY*PIJ(IJ) 
      DIDISIJ(1,3,IJ)=DIDISIJ(1,3,IJ)+ DUMXZ*PIJ(IJ) 
      DIDISIJ(2,1,IJ)=DIDISIJ(2,1,IJ)+ DUMXY*PIJ(IJ) 
      DIDISIJ(2,2,IJ)=DIDISIJ(2,2,IJ)+ DUMYY*PIJ(IJ) 
      DIDISIJ(2,3,IJ)=DIDISIJ(2,3,IJ)+ DUMYZ*PIJ(IJ) 
      DIDISIJ(3,1,IJ)=DIDISIJ(3,1,IJ)+ DUMXZ*PIJ(IJ) 
      DIDISIJ(3,2,IJ)=DIDISIJ(3,2,IJ)+ DUMYZ*PIJ(IJ) 
      DIDISIJ(3,3,IJ)=DIDISIJ(3,3,IJ)+ DUMZZ*PIJ(IJ) 
C
      DUMXX=DXSDIJ(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMXY=DXSDI (IX,JX)*DYSDJ (IY,JY)* ZS   (IZ,JZ)
      DUMXZ=DXSDI (IX,JX)* YS   (IY,JY)*DZSDJ (IZ,JZ)
      DUMYX=DXSDJ (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMYY= XS   (IX,JX)*DYSDIJ(IY,JY)* ZS   (IZ,JZ)
      DUMYZ= XS   (IX,JX)*DYSDI (IY,JY)*DZSDJ (IZ,JZ)
      DUMZX=DXSDJ (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      DUMZY= XS   (IX,JX)*DYSDJ (IY,JY)*DZSDI (IZ,JZ)
      DUMZZ= XS   (IX,JX)* YS   (IY,JY)*DZSDIJ(IZ,JZ)
      DIDJSIJ(1,1,IJ)=DIDJSIJ(1,1,IJ)+ DUMXX*PIJ(IJ) 
      DIDJSIJ(1,2,IJ)=DIDJSIJ(1,2,IJ)+ DUMXY*PIJ(IJ) 
      DIDJSIJ(1,3,IJ)=DIDJSIJ(1,3,IJ)+ DUMXZ*PIJ(IJ) 
      DIDJSIJ(2,1,IJ)=DIDJSIJ(2,1,IJ)+ DUMYX*PIJ(IJ) 
      DIDJSIJ(2,2,IJ)=DIDJSIJ(2,2,IJ)+ DUMYY*PIJ(IJ) 
      DIDJSIJ(2,3,IJ)=DIDJSIJ(2,3,IJ)+ DUMYZ*PIJ(IJ) 
      DIDJSIJ(3,1,IJ)=DIDJSIJ(3,1,IJ)+ DUMZX*PIJ(IJ) 
      DIDJSIJ(3,2,IJ)=DIDJSIJ(3,2,IJ)+ DUMZY*PIJ(IJ) 
      DIDJSIJ(3,3,IJ)=DIDJSIJ(3,3,IJ)+ DUMZZ*PIJ(IJ) 

  380 CONTINUE
  390 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_SD_OLD(FC,LFC,DS,LDS,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION FC(LFC,1),DS(LDS,1),DAB(1)
      DIMENSION DIJ(225),PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION  XS   (7,6), YS   (7,6), ZS   (7,6)
      DIMENSION DXSDI (5,5),DYSDI (5,5),DZSDI (5,5)
      DIMENSION DXSDJ (5,5),DYSDJ (5,5),DZSDJ (5,5)
      DIMENSION DXSDII(5,5),DYSDII(5,5),DZSDII(5,5)
      DIMENSION DXSDIJ(5,5),DYSDIJ(5,5),DZSDIJ(5,5)
      DIMENSION DNAM(3)
      DATA NFT17  /17/
      DATA NFT18  /18/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =.TRUE. 
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- CALCULATE -S- DERIVATIVE TERM -----
C
      IF(OUT) WRITE(IW,9999)
      NDER=2
C
C     ----- INITIALIZE SH COUNTING FOR PARALLEL PROC. -----
C
      IJSH=0
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITDER=LIT+NDER
      IAT   =I
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C
      IJSH=IJSH+1
      IF(MOD(IJSH,NAP)+1.NE.IAP) GO TO 7900
C
      LJTDER=LJT+NDER-1
      JAT=J
C
      IANDJ=II.EQ.JJ
      IF(IANDJ) GO TO 8000
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
C     NN=IA(MAX0(LOCI+I,LOCJ+J))
C    1  +   MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DEN=DAB(NN)
      DEN=DEN+DEN
C
C     ----- NOTE: DEN FROM NWCHEM = - DEN FROM HONDO
C
      DEN = -DEN
C
      IJ=IJ+1
      DIJ(IJ)=DUM2*DEN
      PIJ(IJ)=DUM2
      IF(IANDJ.AND.I.EQ.J) PIJ(IJ)=PIJ(IJ)+PIJ(IJ)

  360 CONTINUE
C
C     ----- OVERLAP DERIVATIVES -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTDER
      NJ=J
      DO 370 I=1,LITDER
      NI=I

      CALL HND_DSXYZ
      XS(I,J)=XINT*T
      YS(I,J)=YINT*T
      ZS(I,J)=ZINT*T
  370 CONTINUE
C
      CALL HND_DDERIJ(DXSDII,DYSDII,DZSDII,DXSDIJ,DYSDIJ,DZSDIJ,
     1            DXSDI ,DYSDI ,DZSDI ,DXSDJ ,DYSDJ ,DZSDJ ,
     2            XS,YS,ZS,LIT,LJT,AI,AJ)
C
      IC=3*(IAT-1)
      JC=3*(JAT-1)
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      DO 380 J=MINJ,MAXJ
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX =DXSDI (IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMY = XS   (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMZ = XS   (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      IJ=IJ+1
      DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
      DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
      DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
C
C     NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
      IDIM=MAX0(LOCI+I,LOCJ+J)
      JDIM=MIN0(LOCI+I,LOCJ+J)
      NN=(IDIM*(IDIM-1))/2+JDIM
C
      DS(NN,IC+1)=DS(NN,IC+1)+ DUMX*PIJ(IJ)
      DS(NN,IC+2)=DS(NN,IC+2)+ DUMY*PIJ(IJ)
      DS(NN,IC+3)=DS(NN,IC+3)+ DUMZ*PIJ(IJ)
C
      DUMXX=DXSDII(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMYY= XS   (IX,JX)*DYSDII(IY,JY)* ZS   (IZ,JZ)
      DUMZZ= XS   (IX,JX)* YS   (IY,JY)*DZSDII(IZ,JZ)
      DUMXY=DXSDI (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMXZ=DXSDI (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      DUMYZ= XS   (IX,JX)*DYSDI (IY,JY)*DZSDI (IZ,JZ)
      FC(IC+1,IC+1)=FC(IC+1,IC+1)+(DUMXX*DIJ(IJ))
      FC(IC+1,IC+2)=FC(IC+1,IC+2)+(DUMXY*DIJ(IJ))
      FC(IC+1,IC+3)=FC(IC+1,IC+3)+(DUMXZ*DIJ(IJ))
      FC(IC+2,IC+1)=FC(IC+2,IC+1)+(DUMXY*DIJ(IJ))
      FC(IC+2,IC+2)=FC(IC+2,IC+2)+(DUMYY*DIJ(IJ))
      FC(IC+2,IC+3)=FC(IC+2,IC+3)+(DUMYZ*DIJ(IJ))
      FC(IC+3,IC+1)=FC(IC+3,IC+1)+(DUMXZ*DIJ(IJ))
      FC(IC+3,IC+2)=FC(IC+3,IC+2)+(DUMYZ*DIJ(IJ))
      FC(IC+3,IC+3)=FC(IC+3,IC+3)+(DUMZZ*DIJ(IJ))
C
      DUMXX=DXSDIJ(IX,JX)* YS   (IY,JY)* ZS   (IZ,JZ)
      DUMXY=DXSDI (IX,JX)*DYSDJ (IY,JY)* ZS   (IZ,JZ)
      DUMXZ=DXSDI (IX,JX)* YS   (IY,JY)*DZSDJ (IZ,JZ)
      DUMYX=DXSDJ (IX,JX)*DYSDI (IY,JY)* ZS   (IZ,JZ)
      DUMYY= XS   (IX,JX)*DYSDIJ(IY,JY)* ZS   (IZ,JZ)
      DUMYZ= XS   (IX,JX)*DYSDI (IY,JY)*DZSDJ (IZ,JZ)
      DUMZX=DXSDJ (IX,JX)* YS   (IY,JY)*DZSDI (IZ,JZ)
      DUMZY= XS   (IX,JX)*DYSDJ (IY,JY)*DZSDI (IZ,JZ)
      DUMZZ= XS   (IX,JX)* YS   (IY,JY)*DZSDIJ(IZ,JZ)
      FC(IC+1,JC+1)=FC(IC+1,JC+1)+(DUMXX*DIJ(IJ))
      FC(IC+1,JC+2)=FC(IC+1,JC+2)+(DUMXY*DIJ(IJ))
      FC(IC+1,JC+3)=FC(IC+1,JC+3)+(DUMXZ*DIJ(IJ))
      FC(IC+2,JC+1)=FC(IC+2,JC+1)+(DUMYX*DIJ(IJ))
      FC(IC+2,JC+2)=FC(IC+2,JC+2)+(DUMYY*DIJ(IJ))
      FC(IC+2,JC+3)=FC(IC+2,JC+3)+(DUMYZ*DIJ(IJ))
      FC(IC+3,JC+1)=FC(IC+3,JC+1)+(DUMZX*DIJ(IJ))
      FC(IC+3,JC+2)=FC(IC+3,JC+2)+(DUMZY*DIJ(IJ))
      FC(IC+3,JC+3)=FC(IC+3,JC+3)+(DUMZZ*DIJ(IJ))

  380 CONTINUE
  390 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
 7900 CONTINUE
C
      IF(.NOT.OUT) GO TO 8000
      WRITE(IW,9993) II,JJ
      MAX=0
 7500 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 7600 I=1,3
 7600 WRITE(IW,9994) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 7500
C
 8000 CONTINUE
 9000 CONTINUE
C
C     ----- GLOBAL SUM OF DE, FC AND DF FROM ALL PROCESSORS -----
C
      NC=3*NAT
      NUM2=(NUM*(NUM+1))/2
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD PARTIAL FC AND DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,FC,(NC*NC),'+')
C
C     ----- ADD PARTIAL DF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(227,DS,(NC*NUM2),'+')
C
      IF(.NOT.OUT) GO TO 9300
      MAX=0
 9100 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9996)
      WRITE(IW,9995) (I,I=MIN,MAX)
      WRITE(IW,9996)
      DO 9200 I=1,3
 9200 WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 9100
C
      NC=3*NAT
      WRITE(IW,9991)
      CALL HND_HSSOUT(FC,NC,NC,NC)
      WRITE(IW,9990)
      DO 9250 IC=1,3*NAT
         CALL HND_PRTR(DS(1,IC),NUM)
         CALL HND_PRTRL(DS(1,IC),NUM)
 9250 CONTINUE
C
      WRITE(IW,9998)
 9300 CONTINUE
C
C     ----- WRITE -DS- ON -NFT18- AND -FC- ON -FT17- -----
C
C     CALL HND_REWFIL(NFT18)
C     DO 9400 IC=1,3*NAT
C     CALL HND_HSSWRT(NFT18,DS(1,IC),LDS)
C9400 CONTINUE
C     CALL HND_HSSWRT(NFT17,FC,LFC*LFC)
C
      RETURN
 9999 FORMAT(/,10X,19(1H-),/,10X,'-S- DERIVATIVE TERM',
     1       /,10X,19(1H-))
 9998 FORMAT(' ...... END OF -S- DERIVATIVE TERM ......')
 9996 FORMAT(/)
 9995 FORMAT(5X,'ATOM',8(6X,I3,6X))
 9994 FORMAT(6X,A3,8F15.7)
 9993 FORMAT(' SHELLS II, JJ ',2I5)
 9992 FORMAT(6X,A3,8F15.7)
 9991 FORMAT(/,10X,16(1H-),/,10X,'-S- HESSIAN TERM',
     1       /,10X,16(1H-))
 9990 FORMAT(/,10X,21(1H-),/,10X,'-S- DERIVATIVE MATRIX',
     1       /,10X,21(1H-))
      END

      SUBROUTINE HND_HSSWRT(NFT,F,LF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(LF)
      WRITE(NFT) F
      RETURN
      END

      SUBROUTINE HND_HSSRD(NFT,F,LF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(LF)
      READ(NFT) F
      RETURN
      END

      SUBROUTINE HND_HSSPRT(FC,NC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*4 CLAB
      character*8 anam,bflab
      character*2 bnam
      PARAMETER (NUMMX=3072)
      PARAMETER (MXATOM=500)
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(NUMMX)
      DIMENSION FC(NC,1),CLAB(3)
      DATA CLAB /'   X','   Y','   Z'/
C
      WRITE(IW,9993)
      MAX=0
   20 MIN=MAX+1
      MAX=MAX+4
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,9999)
      WRITE(IW,9998) (N,N=MIN,MAX)
      WRITE(IW,9999)
      WRITE(IW,9997) (ANAM(N),BNAM(N),N=MIN,MAX)
      WRITE(IW,9999)
      WRITE(IW,9996) ((CLAB(M),M=1,3),N=MIN,MAX)
      WRITE(IW,9999)
      J0=3*(MIN-1)+1
      J1=3*MAX
      DO 30 IAT=1,NAT
      I0=3*(IAT-1)
      WRITE(IW,9995) IAT,ANAM(IAT),BNAM(IAT),CLAB(1),
     1               (FC(I0+1,J),J=J0,J1)
      WRITE(IW,9994) CLAB(2),(FC(I0+2,J),J=J0,J1)
      WRITE(IW,9994) CLAB(3),(FC(I0+3,J),J=J0,J1)
   30 CONTINUE
      IF(MAX.LT.NAT) GO TO 20
      RETURN
 9999 FORMAT(/)
 9998 FORMAT(20X,4(12X,I3,12X))
 9997 FORMAT(20X,4(10X,A8,A2,7X))
 9996 FORMAT(20X,4(2X,A4,5X,A4,5X,A4,3X))
 9995 FORMAT(I3,3X,A8,A2,A4,12F9.4)
 9994 FORMAT(16X,A4,12F9.4)
 9993 FORMAT(10X,49(1H-),
     1     /,10X,'Cartesian force constant matrix (hartree/bohr**2)',
     2     /,10X,49(1H-))
      END

      SUBROUTINE HND_HSSOUT(V,M,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION V(NDIM,1)
      MAX=7
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      WRITE (IW,9128) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
 9008 FORMAT(/)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.7)
      END

C
C     ------ SUBROUTINE CREATED ON AUG.6, 1996 ------
C
      SUBROUTINE HND_DENDDM(D,DA,NUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(*),DA(*)
      IJ=0
      DO 10 I=1,NUM
         DO 20 J=1,I
            IJ=IJ+1
            D(IJ) =1.0
            DA(IJ)=1.0
 20      CONTINUE
 10      CONTINUE
      RETURN
      END
C
C
C     ------ BELOW TAKEN FROM PRP.F IN HONDO ------
C

      SUBROUTINE HND_DSXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                               ,NI,NJ
     1                              ,CX,CY,CZ
      COMMON/HND_HERMIT/H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON/HND_WERMIT/W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
      DIMENSION H(28),W(28),MIN(7),MAX(7)
      EQUIVALENCE (H(1),H1),(W(1),W1)
      DATA MIN  /1,2,4, 7,11,16,22/
      DATA MAX  /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 16 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)*T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    7 GO TO (15,14,13,12,11,10,9,8),NJ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   13 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   14 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   15 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
   16 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_DVXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                              ,NI,NJ
     1                              ,CX,CY,CZ
      COMMON/HND_HERMIT/H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON/HND_WERMIT/W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
      DIMENSION H(28),W(28),MIN(7),MAX(7)
      EQUIVALENCE (H(1),H1),(W(1),W1)
      DATA MIN  /1,2,4, 7,11,16,22/
      DATA MAX  /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      NPTS=(NI+NJ+1-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 12 I=IMIN,IMAX
      DUM=H(I)*T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      PX=PTX-CX
      PY=PTY-CY
      PZ=PTZ-CZ
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 GO TO (11,10,9,8,7),NJ
    7 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 DUM=W(I)
      XINT=XINT+DUM*PX
      YINT=YINT+DUM*PY
      ZINT=ZINT+DUM*PZ
   12 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_DDVXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                              ,NI,NJ
     1                              ,CX,CY,CZ
      COMMON/HND_HERMIT/H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON/HND_WERMIT/W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
      DIMENSION H(28),W(28),MIN(7),MAX(7)
      EQUIVALENCE (H(1),H1),(W(1),W1)
      DATA MIN  /1,2,4, 7,11,16,22/
      DATA MAX  /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ+2-2)/2+1
      IMIN = MIN(NPTS)
      IMAX = MAX(NPTS)
      DO 11 I = IMIN,IMAX
      DUM = H(I)*T
      PTX = DUM+X0
      PTY = DUM+Y0
      PTZ = DUM+Z0
      PX = PTX-CX
      PY = PTY-CY
      PZ = PTZ-CZ
      PX = PX*PX
      PY = PY*PY
      PZ = PZ*PZ
      AX = PTX-XI
      AY = PTY-YI
      AZ = PTZ-ZI
      BX = PTX-XJ
      BY = PTY-YJ
      BZ = PTZ-ZJ
      GO TO (5,4,3,2,1),NI
    1 PX = PX*AX
      PY = PY*AY
      PZ = PZ*AZ
    2 PX = PX*AX
      PY = PY*AY
      PZ = PZ*AZ
    3 PX = PX*AX
      PY = PY*AY
      PZ = PZ*AZ
    4 PX = PX*AX
      PY = PY*AY
      PZ = PZ*AZ
    5 GO TO (10,9,8,7,6),NJ
    6 PX = PX*BX
      PY = PY*BY
      PZ = PZ*BZ
    7 PX = PX*BX
      PY = PY*BY
      PZ = PZ*BZ
    8 PX = PX*BX
      PY = PY*BY
      PZ = PZ*BZ
    9 PX = PX*BX
      PY = PY*BY
      PZ = PZ*BZ
   10 DUM = W(I)
      XINT = XINT+DUM*PX
      YINT = YINT+DUM*PY
      ZINT = ZINT+DUM*PZ
   11 CONTINUE
      RETURN
      END

C
C     -------ABOVE TAKEN FROM PRP.F IN HONDO ------
C

C
C     ------ TAKEN FROM DER.F IN HONDO -------
C
      SUBROUTINE HND_DERI(DXDI,DYDI,DZDI,X,Y,Z,LIT,LJT,AI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X  (6,1), Y  (6,1), Z  (6,1)
      DIMENSION DXDI(5,1),DYDI(5,1),DZDI(5,1)
C
C     ----- DXDI ... -----
C
      DO 110 J=1,LJT
      DXDI(1,J)=X(2,J)*(AI+AI)
      DYDI(1,J)=Y(2,J)*(AI+AI)
      DZDI(1,J)=Z(2,J)*(AI+AI)
  110 CONTINUE
C
      IF(LIT.EQ.1) RETURN
C
      DO 120 I=2,LIT
      DO 120 J=1,LJT
      DXDI(I,J)=X(I+1,J)*(AI+AI)-X(I-1,J)*DBLE(I-1)
      DYDI(I,J)=Y(I+1,J)*(AI+AI)-Y(I-1,J)*DBLE(I-1)
      DZDI(I,J)=Z(I+1,J)*(AI+AI)-Z(I-1,J)*DBLE(I-1)
  120 CONTINUE
C
      RETURN
      END
C
C
C-----------------------------------------------------------------
C
C     ------ SUBROUTINES TAKEN FROM HSS.F IN HONDO ------
C
      SUBROUTINE HND_JKDD(RTDB,GEOM,BASIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER RTDB,GEOM,BASIS
C
      PARAMETER (MXATOM=500)
      COMMON/HND_GRAD12/DE(3,MXATOM)
C 
      CALL HND_DDERJK(RTDB,GEOM,BASIS)
C
      CALL HND_DDER2_POPM
C
      RETURN
      END
C
C     ----- ASSIGN SCFTYP, WFNTYP AND ...-----
C
      BLOCK DATA HND_TYP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 WFNTYP,SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
C
      DATA WFNTYP  /'SCF     '/
      DATA SCFTYP  /'RHF     '/
C
      END
      subroutine hnd_dder2_init(rtdb)
      implicit double precision (a-h,o-z)
      parameter (mxatom=500)
#include "mafdecls.fh"
#include "global.fh"
c     logical status
      integer rtdb 
      integer g_dens, g_wdens
c 
      common/hnd_dermem/iwfn,iclb,inijg,igint,ifint,isint,iijklg,
     1 idab,ichrg,ixy,ixz,iyz,ix,iy,iz,isj,isk,isl,igijkl,ignkl,ignm,
     2 idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     3 iaai,iaaj,ibbk,ibbl,ifi,ifj,ifk,ifl,
     4 isii,isjj,iskk,isll,isij,isik,isil,isjk,isjl,iskl,
     5 idijsi,idijsj,idklsk,idklsl,iabv,icv,irw
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_intfil/nintmx
      common/ hnd_mcpar/nfzc,ncorbs,nci
      common/ hnd_needm/needm
      common/hnd_dder2_mem/i_init2,init2
      common/hnd_iofile/ir,iw,ip
c
      dimension g_dens(3) ! max # of density matrices
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      data zero    /0.0d+00/
c
      do i = 1,3
         g_dens(i) = 0
      enddo
      g_wdens = 0 
c
      nc= nat*3
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
c
      nintmx = 0
      nci    = 0
      ncorbs = 0
      nfzc   = 0

c     ----- memory blocks -----

      need =  0
c*    needa = num*num + num*num
      needa = 0
      need  = needm + needa
   
c     ----- memory block to store 1-e results fc + df -----

      needc = nc*nc +nc*l2
      need  = need + needc
c
c*     if(.not.ma_push_get(mt_dbl,need , 'mem init2',i_init2,init2))
c*   & call errquit('hnd_dder2_init, malloc of init2  failed',911)
      if(.not.ma_push_get(mt_dbl,need,'mem init2',i_init2,init2)) then
         write(iw,*)'The total number of basis functions =',num
         write(iw,*)'The stack memory for dder_2e =',need
         call errquit('hnd_dder2_init, malloc of init2  failed',911)
      endif
c
      return
      end
c
      subroutine hnd_dder2_inita(rtdb)
      implicit double precision (a-h,o-z)
      parameter (mxatom=500)
#include "mafdecls.fh"
#include "global.fh"
      logical status
      integer rtdb 
      integer g_dens, g_wdens
      common/hnd_dermem/iwfn,iclb,inijg,igint,ifint,isint,iijklg,
     1 idab,ichrg,ixy,ixz,iyz,ix,iy,iz,isj,isk,isl,igijkl,ignkl,ignm,
     2 idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     3 iaai,iaaj,ibbk,ibbl,ifi,ifj,ifk,ifl,
     4 isii,isjj,iskk,isll,isij,isik,isil,isjk,isjl,iskl,
     5 idijsi,idijsj,idklsk,idklsl,iabv,icv,irw
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_dder2_mem/i_init2,init2
      common/hnd_basnum/num
      common/hnd_intfil/nintmx
      common/ hnd_mcpar/nfzc,ncorbs,nci
      common/ hnd_needm/needm
c
      dimension g_dens(3) ! max # of density matrices
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      data zero    /0.0d+00/
c
      do i = 1,3
        g_dens(i) = 0
      enddo
      g_wdens = 0 
c
      nc= nat*3
      l1= num
      l2=(num*(num+1))/2
      l3= num*num

c     ----- memory block of density matrix -----

c*    needa = num*num + num*num
      needa = 0       + 0

c     ----- total memory block -----

      need  = needm   + needa
c
      i00=iwfn
      i01=i00+ nc* nc
      i02=i01+ nc* l2
      i10=i02
c
      inita = init2 + needm
c*    i40   = inita + num*num
c*    i50   = i40   + num*num
c
      do 10 i=1,need
   10 x(i+init2-1)=zero
c
c*    call hnd_dens_path1(rtdb,g_dens,g_wdens)

c      ----- check g_dens(1),g_dens(2),g_dens(3) -----

c*    call ga_get(g_dens(1),1,num,1,num, x(inita),num)

c     ----- shell order from nwchem to hondo -----

c*    call hnd_nwhnd_tran(x(inita),x(i40),num)
c     
c*    ij = 0
c*    do i = 1,num
c*       do j = 1,i
c*          ij = ij + 1
c*          x(i10 +ij-1) = x(i40 + (j-1)*(num)+(i-1))
c*       enddo
c*    enddo
c
c*    status = ga_destroy ( g_wdens )
c*    do i=1, 3
c*      status = ga_destroy ( g_dens(i) )
c*    enddo
      call hnd_get_dens(rtdb,x(i10),num)
c
      return
      end
c
      subroutine hnd_dder2_popm 
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
      common/hnd_dder2_mem/i_init2,init2
c
      if(.not.ma_pop_stack(i_init2))
     & call errquit('hnd_dder2_init, ma_pop_stack of init2
     & failed',911)
c
      return
      end
C
      SUBROUTINE HND_DDERJK(RTDB,GEOM,BASIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
      INTEGER  RTDB,GEOM,BASIS
c
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      CHARACTER*8 ERRMSG
      CHARACTER*8 CPMCHF,RHF,UHF,GVB,SCF,TYPMC,TYPCI
      CHARACTER*8 WFNTYP,SCFTYP
      LOGICAL SOME
      LOGICAL NOCORE,NOOPEN,NOPAIR,SUHF,ROHF,RGVB,ROGVB,MC,CI,MCCI
      COMMON/HND_RESTAR/NREST
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      COMMON/HND_JKDTIM/DG80,DHND
      COMMON/HND_IJPAIR/IA(MXNBF)
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_DEROPT/DERTYP
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      COMMON/HND_INTTYP/NPOPLE,NHONDO
      COMMON/HND_INTFIL/NINTMX
      COMMON/HND_SHLBAS/MAXTYP,MAXNUM
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SCFPAR/NOCORE,NOOPEN,NOPAIR,SUHF,ROHF,RGVB,ROGVB,
     1 MC,CI,MCCI,MP2,MP4
C    1 MC,CI,MCCI
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/ HND_MCPAR/NFZC,NCORBS,NCI
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DERMEM/IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/HND_DERARG/I10,I20,J10,J20,J30,NORB,
     1              K10,K20,K30,K40,K50,K60,K70,K80,LDM2,
     2              NBUF,I00,NC,I01,L2,L1
C
      COMMON/HND_NEEDM/NEEDM
      COMMON/HND_DDER2_MEM/I_INIT2,INIT2
      DIMENSION DE1(3,MXATOM)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
C
      DIMENSION LENSHL(5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG      /'PROGRAM ','STOP IN ','-DDERJK-'/
      DATA CPMCHF      /'CPMCHF  '/
      DATA RHF,UHF,GVB /'RHF     ','UHF     ','GVB     '/
      DATA SCF         /'SCF     '/
      DATA TYPMC,TYPCI /'MC      ','CI      '/
      DATA LENSHL      / 1, 4,10,20,35/
      DATA ZERO   /0.0D+00/
C
C     ----- THIS CODE FOR GRADIENT + HESSIAN ( -MDER- = 2 ) -----
C
      MDER=2
C
      SOME=.TRUE.
      CALL HND_SECOND(TIM0,WTIM0)
C
      CALL HND_HONDO
C
      SUHF=WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF
      ROHF=WFNTYP.EQ.SCF.AND.SCFTYP.EQ.RHF.AND.NSETO.GT.0
      RGVB=WFNTYP.EQ.SCF.AND.SCFTYP.EQ.GVB
      ROGVB=ROHF.OR.RGVB
      CI=WFNTYP.EQ.TYPCI
      MC=WFNTYP.EQ.TYPMC
      MCCI=MC.OR.CI
C
C     ----- CODE VALID ONLY FOR CLOSED-SHELL CASES -----
C
      IF(SUHF.OR.ROHF.OR.RGVB.OR.MC.OR.CI) THEN
         WRITE(IW,9997)
         CALL HND_HNDERR(3,ERRMSG)
         RETURN
      ENDIF
C
      NC= NAT*3
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
C
      MAXTYP=0
      DO 5 II=1,NSHELL
      IF(KTYPE(II).GT.MAXTYP) MAXTYP=KTYPE(II)
    5 CONTINUE
      LDM2=LENSHL(MAXTYP)
      IF(LDM2.LT.4) LDM2=4
C
C     ----- INITIALIZATION -----
C
      LENRHF=L2
      LENUHF=L2+L2
      LENGVB=L2+L3+L1
      LENMC =L2+L2+L3+(NCI*NCI)*(NCI*NCI)
      LENMC1=NINTMX+NINTMX+ (NCI*(NCI+1))/2
      LENMC2=LDM2*NCI * NCI * NCI
     1      +LDM2*LDM2* NCI * NCI
     2      +LDM2*LDM2* LDM2* NCI
     3      +LDM2*LDM2* LDM2* LDM2
      LENMC =LENMC +MAX0(LENMC1,LENMC2)
      LENCI =NINTMX+NINTMX
C
      LENGTH=LENRHF
      IF(SCFTYP.EQ.UHF) LENGTH=LENUHF
      IF(ROGVB        ) LENGTH=LENGVB
      IF(MC           ) LENGTH=LENMC
      IF(CI           ) LENGTH=LENCI
      LENGTH=LENGTH+ NC* NC
      LENGTH=LENGTH+ NC* L2
C
C
      CALL HND_JKDSET(MDER)
C
C     ----- CALCULATE THE SIZE OF MEMORY BLOCK,IADDR -----
C
      IADDR=1
      CALL HND_JKDMEM(MDER,IADDR,LENGTH)
C
C     ----- GET MEMORY BLOCK -----
C
      NEEDM = IADDR
      CALL HND_DDER2_INIT(RTDB)
C
C     ----- RESCALE MEMORY POINTER ----- 
C
      IADDR=INIT2
      CALL HND_JKDMEM(MDER,IADDR,LENGTH)
C
C     ----- GET DENSITY MATRIX ----
C
      CALL HND_DDER2_INITA(RTDB)
C
      IF(MCCI) GO TO 100
C
C     ----- SCF WAVEFUNCTIONS -----
C

      I00=IWFN
      I01=I00+ NC* NC
      I02=I01+ NC* L2
      I10=I02
      I20=I10+L2
      I30=I20+L2
      J10=I02
      J20=J10+L2
      J30=J20+L3
      J40=J30+L1
C
C     ----- READ WAVEFUNCTION INFORMATION -----
C
C     ----- X(I10) PASS THROUGH COMMON BLOCK -----
C     CALL HND_JKDHF(X(I10),X(I20),X(J10),X(J20),X(J30),
C    1 IA,NUM,NORB,L1)
C 
      CALL HND_JDDHF(X(I00),NC,X(I01),L2,L1)
C
C     ----- STORE AWAY 1E-GRADIENT -----
C
      DO 20 I=1,3*NAT
         DE1(I,1)=DE(I,1)
   20    DE (I,1)=ZERO
C
C
C     ----- POINTERS OF MEMORY BLOCK -----

      INITA = INIT2 + NEEDM

C     ----- MEMORY BLOCK FOR DENSITY MATRIX -----

c*    I40   = INITA + NUM*NUM
c*    I50   = I40   + NUM*NUM
      I40   = INITA + 0
      I50   = I40   + 0
C
C     ----- STORE AWAY 1E-FORCE CONSTANT -----
C
      M00 = I50
      M01 = M00 + NC*NC

      DO 21 I=1,NC*NC
         X (M00+I-1)=X(I00+I-1)
   21    X (I00+I-1)=ZERO     
C
C     ----- STORE AWAY 1E-DERIVATIVE FOCK MATRICE -----
C
      DO 22 I=1,NC*L2
         X (M01+I-1)=X(I01+I-1)
   22    X (I01+I-1)=ZERO 
C
C     ----- MODIFY DENSITY MATRICES AND VECTORS -----
C     ----- BY APPLYING NORMALIZATION FACTORS   -----
C
C.    IF(.NOT.ROGVB) CALL HND_JKDCLU(SUHF,IA,DDA,DDB)
C.    IF(     ROHF ) CALL HND_JKDROH(NORB,IA,DDC,V,L1)
C.    IF(     RGVB ) CALL HND_JKDGVB(NORB,IA,DDC,V,L1)
C
      GO TO 300
  100 CONTINUE
      IF(CI) GO TO 200
C
C     ----- -MC- WAVEFUNCTIONS -----
C
      I00=IWFN
      I01=I00+ NC* NC
      I02=I01+ NC* L2
      K10=I02
      K20=K10+L2
      K30=K20+L2
      K40=K30+L3
      K50=K40+(NCI*NCI)*(NCI*NCI)
      K60=K50+(NCI*(NCI+1))/2
      K70=K60+NINTMX
      K80=K70+NINTMX
C
C     ----- READ WAVEFUNCTION INFORMATION -----
C
      CALL HND_JKDMC(X(K10),X(K20),X(K30),X(K40),X(K50),
     1           X(K60),X(K70),IA,NUM,NCORBS,NCI,L1)
C
      K40=K30+L3
      K50=K40+(NCI*NCI)*(NCI*NCI)
      K60=K50+LDM2* NCI*NCI *NCI
      K70=K60+LDM2*LDM2*NCI *NCI
      K80=K70+LDM2*LDM2*LDM2*NCI
      K90=K80+LDM2*LDM2*LDM2*LDM2
C
C     ----- MODIFY DENSITY MATRICES AND VECTORS -----
C     ----- BY APPLYING NORMALIZATION FACTORS   -----
C
C.    CALL HND_JKDDMC
C
      GO TO 300
  200 CONTINUE
C
C     ----- -CI- WAVEFUNCTIONS -----
C
      NBUF=LENGTH/(2*NINTMX)
      I00=IWFN
      I01=I00+ NC* NC
      I02=I01+ NC* L2
      K10=I02
      K20=K10+NBUF*NINTMX
C
C     ----- READ WAVEFUNCTION INFORMATION -----
C
      CALL HND_JKDCI
      CALL HND_JDDCI(X(I00),NC,X(I01),L2,L1)
C
C     ----- MODIFY DENSITY MATRICES AND VECTORS -----
C     ----- BY APPLYING NORMALIZATION FACTORS   -----
C
C.    CALL HND_JKDDCI
C
C     ----- -DM2- OUT OF -TRF2DM- IS NOT IN CORRECT ORDER FOR -----
C           USING BOTH -G80- AND -HND-. SO FORCE -HND- ONLY
C           FOR POLARIZED BASIS SET.
C
      IPOPLE=NPOPLE
      IF(NPOPLE.EQ.1.AND.NHONDO.EQ.1) NPOPLE=0
C
  300 CONTINUE
C
C     ----- READ INTEGRAL BOUNDS ON -DAF- -----
C
C     
C*     NSH2=(NSHELL*(NSHELL+1))/2
C*     CALL HND_DAREAD(IDAF,IODA,X(ICLB),NSH2,100)
C
C     ----- CALCULATE 2E-DERIVATIVE TERM -----
C
      CALL HND_JDDER
C
C     ----- ADD PARTIAL GRADIENTS IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(225,DE,3*NAT,'+')
C
C     ----- ADD 1E-GRADIENTS -----
C
      DO 410 I=1,3*NAT
  410    DE(I,1)=DE(I,1)+DE1(I,1)
C
C     ----- ADD PARTIAL FC AND DF IF DISTRIBUTED PROCESSING -----
C
      CALL GA_DGOP(226,X(I00),(NC*NC+NC*L2),'+')
C
C     ----- ADD 1E-FORCE CONSTANT -----
C
      DO 420 I=1,NC*NC
  420    X(I00+I-1)=X(I00+I-1)+X(M00+I-1)
C
C     ----- ADD 1E-DERIVATIVE FOCK MATRICES -----
C
      DO 430 I=1,NC*L2
  430    X(I01+I-1)=X(I01+I-1)+X(M01+I-1)
C
C     ----- IF FINISHED, THEN SYMMETRIZE -----
C
      NREST = 0 
      IF(NREST.EQ.0) CALL HND_JKDSYM(DE)
      IF(NREST.EQ.0) CALL HND_JDDSYM(X(I00),NC,X(I01),L2)
C
C     ----- IF FINISHED, THEN PRINT AND SAVE -----
C
      IF (GA_NODEID().EQ.0) THEN
        IF(NREST.EQ.0) CALL HND_JDDSAV(X(I00),NC,X(I01),L2,L1)
        IF(NREST.EQ.0) CALL HND_JKDSAV(DE,MDER)
      ENDIF
C
C     ----- RESET OPTION IF -CI- -----
C
      IF(CI) NPOPLE=IPOPLE
C
      RETURN
C
 9999 FORMAT(' NOT ENOUGH CORE IN -DDERJK- . STOP . ')
c9998 FORMAT(10X,'ELAPSED - CPU- TIME = ',F10.3,
c    1        5X,'TOTAL - CPU- TIME = ',F10.3)
 9997 FORMAT(' SECOND-DERIVATIVES CALCULATIONS FOR',
     1       ' CLOSED-SHELL ONLY. STOP            ')
c9996 FORMAT(10X,'ELAPSED -WALL- TIME = ',F10.3,
c    1        5X,'TOTAL -WALL- TIME = ',F10.3)
C
      END
C
      SUBROUTINE HND_JDDHF(FC,NC,DF,LDF,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
C     COMMON/HND_DSHPRT/OUT
      DIMENSION FC(NC,NC),DF(LDF,NC)
      DATA NFT17 /17/
C
      ENTRY      HND_JDDCI(FC,NC,DF,LDF,NBF)
C
C     ----- READ DERIVATIVE FOCK MATRICES -----
C     ----- DIAGONAL IS MULTIPLIED BY 2.0 -----
C
      OUT = .FALSE.
C
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9999)
C
      CALL HND_REWFIL(NFT17)
      DO 20 IC=1,NC
      CALL HND_HSSRD(NFT17,DF(1,IC),LDF)
C
      IJ=0
      DO 10 I=1,NBF
      DO 10 J=1,I
      IJ=IJ+1
      IF(J.EQ.I) GO TO 10
      DF(IJ,IC)=DF(IJ,IC)+DF(IJ,IC)
   10 CONTINUE
C
c*     IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRTRL(DF(1,IC),NBF)
C*     IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRTR (DF(1,IC),NBF)
C
   20 CONTINUE
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRTR(DF(1,1 ),NBF)
C
C     ----- READ 1E-HESSIAN -----
C
      CALL HND_HSSRD(NFT17,FC,NC*NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9998)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_HSSOUT(FC,NC,NC,NC)
C
      RETURN
C9999 FORMAT(/,' DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
C9998 FORMAT(/,' HESSIAN MATRIX',/,1X,14(1H-))
 9999 FORMAT(/,' 1E-DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
 9998 FORMAT(/,' 1E-HESSIAN MATRIX',/,1X,14(1H-))
C
      END
C
      SUBROUTINE HND_JDDER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL SOME
      LOGICAL SPSKIP
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_JKDTIM/DG80,DHND
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      COMMON/HND_INTHND/SPSKIP
      COMMON/HND_DERMEM/IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/HND_DERARG/I10,I20,J10,J20,J30,NORB,
     1              K10,K20,K30,K40,K50,K60,K70,K80,LDM2,
     2              NBUF,I00,NC,I01,L2,L1
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DATA ZERO /0.0D+00/
C
      SOME=.FALSE.
C     CALL HND_SECOND(TIM0,WTIM0)
C
      SPSKIP=.FALSE.
C
      CALL HND_OEDHND(X(INIJG))

      total_wall = -util_wallsec()
      total_cpu  = -util_cpusec()

      CALL HND_JKDHND(X(ICLB ),X(INIJG),
     1    X(I10),X(I20),X(J10),X(J20),X(J30),NORB,
     2    X(K10),X(K20),X(K30),X(K40),X(K50),X(K60),X(K70),X(K80),LDM2,
     3    X(K10),X(K20),NBUF,
     4    X(I00),NC,X(I01),L2,L1)
C
      total_wall = total_wall + util_wallsec()
      total_cpu  = total_cpu  + util_cpusec()
      if (some.and.ga_nodeid().eq.0) write(6,1) total_cpu, total_wall
 1          format(/' Times in hnd_jkdhnd   cpu:',
     &                f11.1,'s     wall:',f11.1,'s')
            call util_flush(6)

      CALL HND_SECOND(TIM,WTIM)
      DHND=TIM-TIM0-DG80
      DG80=ZERO
C
      RETURN
      END
C
      SUBROUTINE HND_DDSPDV(NG,NR,IJKLG,GIJKL,FIJKL,SIJKL,
     1 XYZ,
     1 FIXYZ,FJXYZ,FKXYZ,FLXYZ,SIIXYZ,SJJXYZ,SKKXYZ,SLLXYZ,
     2 SIJXYZ,SIKXYZ,SILXYZ,SJKXYZ,SJLXYZ,SKLXYZ,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IS,JS,KS,LS
      LOGICAL IJS,IJKS,IJKLS
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,
     1 SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1 LGXYZ(4,35)
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION    XYZ(NG*NR,1)
      DIMENSION  FIXYZ(NG*NR,1), FJXYZ(NG*NR,1)
      DIMENSION  FKXYZ(NG*NR,1), FLXYZ(NG*NR,1)
      DIMENSION SIIXYZ(NG*NR,1),SJJXYZ(NG*NR,1)
      DIMENSION SKKXYZ(NG*NR,1),SLLXYZ(NG*NR,1)
      DIMENSION SIJXYZ(NG*NR,1),SIKXYZ(NG*NR,1)
      DIMENSION SILXYZ(NG*NR,1),SJKXYZ(NG*NR,1)
      DIMENSION SJLXYZ(NG*NR,1),SKLXYZ(NG*NR,1)
      DIMENSION XY(NG*NR),XZ(NG*NR),YZ(NG*NR)
      DIMENSION  X(NG*NR), Y(NG*NR), Z(NG*NR)
      DIMENSION SJ(NG*NR),SK(NG*NR),SL(NG*NR)
      DIMENSION DIJSI(NG*NR),DIJSJ(NG*NR),DKLSK(NG*NR),DKLSL(NG*NR)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
C
C     ----- FIRST DERIVATIVE INTEGRALS + GRADIENT + HESSIAN -----
C
      IF(SPIJKL) GO TO 1000
C
C     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
C
      IJKLN=0
      DO 440 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 430 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 420 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 410 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 110 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  110 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
C
      IF(SKIPI) GO TO 130
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 120 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SIIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SIIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SIIXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FIXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FIXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FIXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
  120 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
      SD( 1)=SD( 1)+DAB(NN)*DUMSXX
      SD( 3)=SD( 3)+DAB(NN)*DUMSYY
      SD( 6)=SD( 6)+DAB(NN)*DUMSZZ
      SD( 2)=SD( 2)+DAB(NN)*DUMSXY
      SD( 4)=SD( 4)+DAB(NN)*DUMSXZ
      SD( 5)=SD( 5)+DAB(NN)*DUMSYZ
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 130
C-    SIJKL( 1,NN)=SIJKL( 1,NN)+DUMSXX
C-    SIJKL( 3,NN)=SIJKL( 3,NN)+DUMSYY
C-    SIJKL( 6,NN)=SIJKL( 6,NN)+DUMSZZ
C-    SIJKL( 2,NN)=SIJKL( 2,NN)+DUMSXY
C-    SIJKL( 4,NN)=SIJKL( 4,NN)+DUMSXZ
C-    SIJKL( 5,NN)=SIJKL( 5,NN)+DUMSYZ
  130 IF(SKIPJ) GO TO 150
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 140 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SJJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SJJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SJJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
  140 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
      SD(10)=SD(10)+DAB(NN)*DUMSXX
      SD(15)=SD(15)+DAB(NN)*DUMSYY
      SD(21)=SD(21)+DAB(NN)*DUMSZZ
      SD(14)=SD(14)+DAB(NN)*DUMSXY
      SD(19)=SD(19)+DAB(NN)*DUMSXZ
      SD(20)=SD(20)+DAB(NN)*DUMSYZ
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 150
C-    SIJKL(10,NN)=SIJKL(10,NN)+DUMSXX
C-    SIJKL(15,NN)=SIJKL(15,NN)+DUMSYY
C-    SIJKL(21,NN)=SIJKL(21,NN)+DUMSZZ
C-    SIJKL(14,NN)=SIJKL(14,NN)+DUMSXY
C-    SIJKL(19,NN)=SIJKL(19,NN)+DUMSXZ
C-    SIJKL(20,NN)=SIJKL(20,NN)+DUMSYZ
  150 IF(SKIPK) GO TO 170
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 160 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SKKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SKKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SKKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FKXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FKXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FKXYZ(IGR,NZ)*XYZ(IGR,NX)
  160 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
      SD(28)=SD(28)+DAB(NN)*DUMSXX
      SD(36)=SD(36)+DAB(NN)*DUMSYY
      SD(45)=SD(45)+DAB(NN)*DUMSZZ
      SD(35)=SD(35)+DAB(NN)*DUMSXY
      SD(43)=SD(43)+DAB(NN)*DUMSXZ
      SD(44)=SD(44)+DAB(NN)*DUMSYZ
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 170
C-    SIJKL(28,NN)=SIJKL(28,NN)+DUMSXX
C-    SIJKL(36,NN)=SIJKL(36,NN)+DUMSYY
C-    SIJKL(45,NN)=SIJKL(45,NN)+DUMSZZ
C-    SIJKL(35,NN)=SIJKL(35,NN)+DUMSXY
C-    SIJKL(43,NN)=SIJKL(43,NN)+DUMSXZ
C-    SIJKL(44,NN)=SIJKL(44,NN)+DUMSYZ
  170 IF(SKIPL) GO TO 190
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 180 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SLLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SLLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SLLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FLXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FLXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FLXYZ(IGR,NZ)*XYZ(IGR,NX)
  180 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
      SD(55)=SD(55)+DAB(NN)*DUMSXX
      SD(66)=SD(66)+DAB(NN)*DUMSYY
      SD(78)=SD(78)+DAB(NN)*DUMSZZ
      SD(65)=SD(65)+DAB(NN)*DUMSXY
      SD(76)=SD(76)+DAB(NN)*DUMSXZ
      SD(77)=SD(77)+DAB(NN)*DUMSYZ
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 190
C-    SIJKL(55,NN)=SIJKL(55,NN)+DUMSXX
C-    SIJKL(66,NN)=SIJKL(66,NN)+DUMSYY
C-    SIJKL(78,NN)=SIJKL(78,NN)+DUMSZZ
C-    SIJKL(65,NN)=SIJKL(65,NN)+DUMSXY
C-    SIJKL(76,NN)=SIJKL(76,NN)+DUMSXZ
C-    SIJKL(77,NN)=SIJKL(77,NN)+DUMSYZ
  190 IF(SKIPI) GO TO 250
      IF(SKIPJ) GO TO 210
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 200 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIJXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIJXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FJXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  200 CONTINUE
      SD( 7)=SD( 7)+DAB(NN)*DUMSXX
      SD( 8)=SD( 8)+DAB(NN)*DUMSXY
      SD( 9)=SD( 9)+DAB(NN)*DUMSXZ
      SD(11)=SD(11)+DAB(NN)*DUMSYX
      SD(12)=SD(12)+DAB(NN)*DUMSYY
      SD(13)=SD(13)+DAB(NN)*DUMSYZ
      SD(16)=SD(16)+DAB(NN)*DUMSZX
      SD(17)=SD(17)+DAB(NN)*DUMSZY
      SD(18)=SD(18)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 210
C-    SIJKL( 7,NN)=SIJKL( 7,NN)+DUMSXX
C-    SIJKL( 8,NN)=SIJKL( 8,NN)+DUMSXY
C-    SIJKL( 9,NN)=SIJKL( 9,NN)+DUMSXZ
C-    SIJKL(11,NN)=SIJKL(11,NN)+DUMSYX
C-    SIJKL(12,NN)=SIJKL(12,NN)+DUMSYY
C-    SIJKL(13,NN)=SIJKL(13,NN)+DUMSYZ
C-    SIJKL(16,NN)=SIJKL(16,NN)+DUMSZX
C-    SIJKL(17,NN)=SIJKL(17,NN)+DUMSZY
C-    SIJKL(18,NN)=SIJKL(18,NN)+DUMSZZ
  210 IF(SKIPK) GO TO 230
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 220 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  220 CONTINUE
      SD(22)=SD(22)+DAB(NN)*DUMSXX
      SD(23)=SD(23)+DAB(NN)*DUMSXY
      SD(24)=SD(24)+DAB(NN)*DUMSXZ
      SD(29)=SD(29)+DAB(NN)*DUMSYX
      SD(30)=SD(30)+DAB(NN)*DUMSYY
      SD(31)=SD(31)+DAB(NN)*DUMSYZ
      SD(37)=SD(37)+DAB(NN)*DUMSZX
      SD(38)=SD(38)+DAB(NN)*DUMSZY
      SD(39)=SD(39)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 230
C-    SIJKL(22,NN)=SIJKL(22,NN)+DUMSXX
C-    SIJKL(23,NN)=SIJKL(23,NN)+DUMSXY
C-    SIJKL(24,NN)=SIJKL(24,NN)+DUMSXZ
C-    SIJKL(29,NN)=SIJKL(29,NN)+DUMSYX
C-    SIJKL(30,NN)=SIJKL(30,NN)+DUMSYY
C-    SIJKL(31,NN)=SIJKL(31,NN)+DUMSYZ
C-    SIJKL(37,NN)=SIJKL(37,NN)+DUMSZX
C-    SIJKL(38,NN)=SIJKL(38,NN)+DUMSZY
C-    SIJKL(39,NN)=SIJKL(39,NN)+DUMSZZ
  230 IF(SKIPL) GO TO 250
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 240 IGR=1,NG*NR
      DUMSXX=DUMSXX+SILXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SILXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SILXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  240 CONTINUE
      SD(46)=SD(46)+DAB(NN)*DUMSXX
      SD(47)=SD(47)+DAB(NN)*DUMSXY
      SD(48)=SD(48)+DAB(NN)*DUMSXZ
      SD(56)=SD(56)+DAB(NN)*DUMSYX
      SD(57)=SD(57)+DAB(NN)*DUMSYY
      SD(58)=SD(58)+DAB(NN)*DUMSYZ
      SD(67)=SD(67)+DAB(NN)*DUMSZX
      SD(68)=SD(68)+DAB(NN)*DUMSZY
      SD(69)=SD(69)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 250
C-    SIJKL(46,NN)=SIJKL(46,NN)+DUMSXX
C-    SIJKL(47,NN)=SIJKL(47,NN)+DUMSXY
C-    SIJKL(48,NN)=SIJKL(48,NN)+DUMSXZ
C-    SIJKL(56,NN)=SIJKL(56,NN)+DUMSYX
C-    SIJKL(57,NN)=SIJKL(57,NN)+DUMSYY
C-    SIJKL(58,NN)=SIJKL(58,NN)+DUMSYZ
C-    SIJKL(67,NN)=SIJKL(67,NN)+DUMSZX
C-    SIJKL(68,NN)=SIJKL(68,NN)+DUMSZY
C-    SIJKL(69,NN)=SIJKL(69,NN)+DUMSZZ
  250 IF(SKIPJ) GO TO 290
      IF(SKIPK) GO TO 270
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 260 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FJXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NY)*XYZ(IGR,NX)
  260 CONTINUE
      SD(25)=SD(25)+DAB(NN)*DUMSXX
      SD(26)=SD(26)+DAB(NN)*DUMSXY
      SD(27)=SD(27)+DAB(NN)*DUMSXZ
      SD(32)=SD(32)+DAB(NN)*DUMSYX
      SD(33)=SD(33)+DAB(NN)*DUMSYY
      SD(34)=SD(34)+DAB(NN)*DUMSYZ
      SD(40)=SD(40)+DAB(NN)*DUMSZX
      SD(41)=SD(41)+DAB(NN)*DUMSZY
      SD(42)=SD(42)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 270
C-    SIJKL(25,NN)=SIJKL(25,NN)+DUMSXX
C-    SIJKL(26,NN)=SIJKL(26,NN)+DUMSXY
C-    SIJKL(27,NN)=SIJKL(27,NN)+DUMSXZ
C-    SIJKL(32,NN)=SIJKL(32,NN)+DUMSYX
C-    SIJKL(33,NN)=SIJKL(33,NN)+DUMSYY
C-    SIJKL(34,NN)=SIJKL(34,NN)+DUMSYZ
C-    SIJKL(40,NN)=SIJKL(40,NN)+DUMSZX
C-    SIJKL(41,NN)=SIJKL(41,NN)+DUMSZY
C-    SIJKL(42,NN)=SIJKL(42,NN)+DUMSZZ
  270 IF(SKIPL) GO TO 290
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 280 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FJXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NY)*XYZ(IGR,NX)
  280 CONTINUE
      SD(49)=SD(49)+DAB(NN)*DUMSXX
      SD(50)=SD(50)+DAB(NN)*DUMSXY
      SD(51)=SD(51)+DAB(NN)*DUMSXZ
      SD(59)=SD(59)+DAB(NN)*DUMSYX
      SD(60)=SD(60)+DAB(NN)*DUMSYY
      SD(61)=SD(61)+DAB(NN)*DUMSYZ
      SD(70)=SD(70)+DAB(NN)*DUMSZX
      SD(71)=SD(71)+DAB(NN)*DUMSZY
      SD(72)=SD(72)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 290
C-    SIJKL(49,NN)=SIJKL(49,NN)+DUMSXX
C-    SIJKL(50,NN)=SIJKL(50,NN)+DUMSXY
C-    SIJKL(51,NN)=SIJKL(51,NN)+DUMSXZ
C-    SIJKL(59,NN)=SIJKL(59,NN)+DUMSYX
C-    SIJKL(60,NN)=SIJKL(60,NN)+DUMSYY
C-    SIJKL(61,NN)=SIJKL(61,NN)+DUMSYZ
C-    SIJKL(70,NN)=SIJKL(70,NN)+DUMSZX
C-    SIJKL(71,NN)=SIJKL(71,NN)+DUMSZY
C-    SIJKL(72,NN)=SIJKL(72,NN)+DUMSZZ
  290 IF(SKIPK) GO TO 400
      IF(SKIPL) GO TO 400
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 300 IGR=1,NG*NR
      DUMSXX=DUMSXX+SKLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SKLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SKLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FKXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FKXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FKXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FKXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NY)*XYZ(IGR,NX)
  300 CONTINUE
      SD(52)=SD(52)+DAB(NN)*DUMSXX
      SD(53)=SD(53)+DAB(NN)*DUMSXY
      SD(54)=SD(54)+DAB(NN)*DUMSXZ
      SD(62)=SD(62)+DAB(NN)*DUMSYX
      SD(63)=SD(63)+DAB(NN)*DUMSYY
      SD(64)=SD(64)+DAB(NN)*DUMSYZ
      SD(73)=SD(73)+DAB(NN)*DUMSZX
      SD(74)=SD(74)+DAB(NN)*DUMSZY
      SD(75)=SD(75)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 400
C-    SIJKL(52,NN)=SIJKL(52,NN)+DUMSXX
C-    SIJKL(53,NN)=SIJKL(53,NN)+DUMSXY
C-    SIJKL(54,NN)=SIJKL(54,NN)+DUMSXZ
C-    SIJKL(62,NN)=SIJKL(62,NN)+DUMSYX
C-    SIJKL(63,NN)=SIJKL(63,NN)+DUMSYY
C-    SIJKL(64,NN)=SIJKL(64,NN)+DUMSYZ
C-    SIJKL(73,NN)=SIJKL(73,NN)+DUMSZX
C-    SIJKL(74,NN)=SIJKL(74,NN)+DUMSZY
C-    SIJKL(75,NN)=SIJKL(75,NN)+DUMSZZ
C
  400 CONTINUE
C
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
  440 CONTINUE
C
C     ----- ZEROTH ORDER INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 940 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 930 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 920 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 910 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 710 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  710 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      DUM=ZERO
      DO 720 IGR=1,NG*NR
  720 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
  910 CONTINUE
  920 CONTINUE
  930 CONTINUE
  940 CONTINUE
C
      RETURN
C
C     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
C
 1000 CONTINUE
C
C     ----- FIRST DERIVATIVE INTEGRALS + GRADIENT + HESSIAN -----
C
      IJKLN=0
      DO 1840 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 1830 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 1110 IGR=1,NG*NR
 1110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 1120 IGR=1,NG*NR
 1120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 1130 IGR=1,NG*NR
 1130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 1820 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 1210 IGR=1,NG*NR
 1210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 1220 IGR=1,NG*NR
 1220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 1230 IGR=1,NG*NR
 1230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 1810 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 1310 IGR=1,NG*NR
 1310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 1320 IGR=1,NG*NR
 1320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 1330 IGR=1,NG*NR
 1330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 1500 IGR=1,NG*NR
          X(IGR)=XYZ(IGR,NX)*SL(IGR)
          Y(IGR)=XYZ(IGR,NY)*SL(IGR)
          Z(IGR)=XYZ(IGR,NZ)*SL(IGR)
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 1500    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 1510 IGR=1,NG*NR
          X(IGR)=XYZ(IGR,NX)
          Y(IGR)=XYZ(IGR,NY)
          Z(IGR)=XYZ(IGR,NZ)
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 1510    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      IF(SKIPI) GO TO 1530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SIIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SIIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SIIXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FIXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FIXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FIXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
 1520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
      SD( 1)=SD( 1)+DAB(NN)*DUMSXX
      SD( 3)=SD( 3)+DAB(NN)*DUMSYY
      SD( 6)=SD( 6)+DAB(NN)*DUMSZZ
      SD( 2)=SD( 2)+DAB(NN)*DUMSXY
      SD( 4)=SD( 4)+DAB(NN)*DUMSXZ
      SD( 5)=SD( 5)+DAB(NN)*DUMSYZ
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1530
C-    SIJKL( 1,NN)=SIJKL( 1,NN)+DUMSXX
C-    SIJKL( 3,NN)=SIJKL( 3,NN)+DUMSYY
C-    SIJKL( 6,NN)=SIJKL( 6,NN)+DUMSZZ
C-    SIJKL( 2,NN)=SIJKL( 2,NN)+DUMSXY
C-    SIJKL( 4,NN)=SIJKL( 4,NN)+DUMSXZ
C-    SIJKL( 5,NN)=SIJKL( 5,NN)+DUMSYZ
 1530 IF(SKIPJ) GO TO 1550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SJJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SJJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SJJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
 1540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
      SD(10)=SD(10)+DAB(NN)*DUMSXX
      SD(15)=SD(15)+DAB(NN)*DUMSYY
      SD(21)=SD(21)+DAB(NN)*DUMSZZ
      SD(14)=SD(14)+DAB(NN)*DUMSXY
      SD(19)=SD(19)+DAB(NN)*DUMSXZ
      SD(20)=SD(20)+DAB(NN)*DUMSYZ
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1550
C-    SIJKL(10,NN)=SIJKL(10,NN)+DUMSXX
C-    SIJKL(15,NN)=SIJKL(15,NN)+DUMSYY
C-    SIJKL(21,NN)=SIJKL(21,NN)+DUMSZZ
C-    SIJKL(14,NN)=SIJKL(14,NN)+DUMSXY
C-    SIJKL(19,NN)=SIJKL(19,NN)+DUMSXZ
C-    SIJKL(20,NN)=SIJKL(20,NN)+DUMSYZ
 1550 IF(SKIPK) GO TO 1570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SKKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SKKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SKKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FKXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FKXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FKXYZ(IGR,NZ)*X(IGR)
 1560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
      SD(28)=SD(28)+DAB(NN)*DUMSXX
      SD(36)=SD(36)+DAB(NN)*DUMSYY
      SD(45)=SD(45)+DAB(NN)*DUMSZZ
      SD(35)=SD(35)+DAB(NN)*DUMSXY
      SD(43)=SD(43)+DAB(NN)*DUMSXZ
      SD(44)=SD(44)+DAB(NN)*DUMSYZ
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1570
C-    SIJKL(28,NN)=SIJKL(28,NN)+DUMSXX
C-    SIJKL(36,NN)=SIJKL(36,NN)+DUMSYY
C-    SIJKL(45,NN)=SIJKL(45,NN)+DUMSZZ
C-    SIJKL(35,NN)=SIJKL(35,NN)+DUMSXY
C-    SIJKL(43,NN)=SIJKL(43,NN)+DUMSXZ
C-    SIJKL(44,NN)=SIJKL(44,NN)+DUMSYZ
 1570 IF(SKIPL) GO TO 1590
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SLLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SLLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SLLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FLXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FLXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FLXYZ(IGR,NZ)*X(IGR)
 1580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
      SD(55)=SD(55)+DAB(NN)*DUMSXX
      SD(66)=SD(66)+DAB(NN)*DUMSYY
      SD(78)=SD(78)+DAB(NN)*DUMSZZ
      SD(65)=SD(65)+DAB(NN)*DUMSXY
      SD(76)=SD(76)+DAB(NN)*DUMSXZ
      SD(77)=SD(77)+DAB(NN)*DUMSYZ
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1590
C-    SIJKL(55,NN)=SIJKL(55,NN)+DUMSXX
C-    SIJKL(66,NN)=SIJKL(66,NN)+DUMSYY
C-    SIJKL(78,NN)=SIJKL(78,NN)+DUMSZZ
C-    SIJKL(65,NN)=SIJKL(65,NN)+DUMSXY
C-    SIJKL(76,NN)=SIJKL(76,NN)+DUMSXZ
C-    SIJKL(77,NN)=SIJKL(77,NN)+DUMSYZ
 1590 IF(SKIPI) GO TO 1650
      IF(SKIPJ) GO TO 1610
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1600 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIJXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIJXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FJXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1600 CONTINUE
      SD( 7)=SD( 7)+DAB(NN)*DUMSXX
      SD( 8)=SD( 8)+DAB(NN)*DUMSXY
      SD( 9)=SD( 9)+DAB(NN)*DUMSXZ
      SD(11)=SD(11)+DAB(NN)*DUMSYX
      SD(12)=SD(12)+DAB(NN)*DUMSYY
      SD(13)=SD(13)+DAB(NN)*DUMSYZ
      SD(16)=SD(16)+DAB(NN)*DUMSZX
      SD(17)=SD(17)+DAB(NN)*DUMSZY
      SD(18)=SD(18)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1610
C-    SIJKL( 7,NN)=SIJKL( 7,NN)+DUMSXX
C-    SIJKL( 8,NN)=SIJKL( 8,NN)+DUMSXY
C-    SIJKL( 9,NN)=SIJKL( 9,NN)+DUMSXZ
C-    SIJKL(11,NN)=SIJKL(11,NN)+DUMSYX
C-    SIJKL(12,NN)=SIJKL(12,NN)+DUMSYY
C-    SIJKL(13,NN)=SIJKL(13,NN)+DUMSYZ
C-    SIJKL(16,NN)=SIJKL(16,NN)+DUMSZX
C-    SIJKL(17,NN)=SIJKL(17,NN)+DUMSZY
C-    SIJKL(18,NN)=SIJKL(18,NN)+DUMSZZ
 1610 IF(SKIPK) GO TO 1630
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1620 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1620 CONTINUE
      SD(22)=SD(22)+DAB(NN)*DUMSXX
      SD(23)=SD(23)+DAB(NN)*DUMSXY
      SD(24)=SD(24)+DAB(NN)*DUMSXZ
      SD(29)=SD(29)+DAB(NN)*DUMSYX
      SD(30)=SD(30)+DAB(NN)*DUMSYY
      SD(31)=SD(31)+DAB(NN)*DUMSYZ
      SD(37)=SD(37)+DAB(NN)*DUMSZX
      SD(38)=SD(38)+DAB(NN)*DUMSZY
      SD(39)=SD(39)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1630
C-    SIJKL(22,NN)=SIJKL(22,NN)+DUMSXX
C-    SIJKL(23,NN)=SIJKL(23,NN)+DUMSXY
C-    SIJKL(24,NN)=SIJKL(24,NN)+DUMSXZ
C-    SIJKL(29,NN)=SIJKL(29,NN)+DUMSYX
C-    SIJKL(30,NN)=SIJKL(30,NN)+DUMSYY
C-    SIJKL(31,NN)=SIJKL(31,NN)+DUMSYZ
C-    SIJKL(37,NN)=SIJKL(37,NN)+DUMSZX
C-    SIJKL(38,NN)=SIJKL(38,NN)+DUMSZY
C-    SIJKL(39,NN)=SIJKL(39,NN)+DUMSZZ
 1630 IF(SKIPL) GO TO 1650
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1640 IGR=1,NG*NR
      DUMSXX=DUMSXX+SILXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SILXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SILXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1640 CONTINUE
      SD(46)=SD(46)+DAB(NN)*DUMSXX
      SD(47)=SD(47)+DAB(NN)*DUMSXY
      SD(48)=SD(48)+DAB(NN)*DUMSXZ
      SD(56)=SD(56)+DAB(NN)*DUMSYX
      SD(57)=SD(57)+DAB(NN)*DUMSYY
      SD(58)=SD(58)+DAB(NN)*DUMSYZ
      SD(67)=SD(67)+DAB(NN)*DUMSZX
      SD(68)=SD(68)+DAB(NN)*DUMSZY
      SD(69)=SD(69)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1650
C-    SIJKL(46,NN)=SIJKL(46,NN)+DUMSXX
C-    SIJKL(47,NN)=SIJKL(47,NN)+DUMSXY
C-    SIJKL(48,NN)=SIJKL(48,NN)+DUMSXZ
C-    SIJKL(56,NN)=SIJKL(56,NN)+DUMSYX
C-    SIJKL(57,NN)=SIJKL(57,NN)+DUMSYY
C-    SIJKL(58,NN)=SIJKL(58,NN)+DUMSYZ
C-    SIJKL(67,NN)=SIJKL(67,NN)+DUMSZX
C-    SIJKL(68,NN)=SIJKL(68,NN)+DUMSZY
C-    SIJKL(69,NN)=SIJKL(69,NN)+DUMSZZ
 1650 IF(SKIPJ) GO TO 1690
      IF(SKIPK) GO TO 1670
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1660 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FJXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NY)*X(IGR)
 1660 CONTINUE
      SD(25)=SD(25)+DAB(NN)*DUMSXX
      SD(26)=SD(26)+DAB(NN)*DUMSXY
      SD(27)=SD(27)+DAB(NN)*DUMSXZ
      SD(32)=SD(32)+DAB(NN)*DUMSYX
      SD(33)=SD(33)+DAB(NN)*DUMSYY
      SD(34)=SD(34)+DAB(NN)*DUMSYZ
      SD(40)=SD(40)+DAB(NN)*DUMSZX
      SD(41)=SD(41)+DAB(NN)*DUMSZY
      SD(42)=SD(42)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1670
C-    SIJKL(25,NN)=SIJKL(25,NN)+DUMSXX
C-    SIJKL(26,NN)=SIJKL(26,NN)+DUMSXY
C-    SIJKL(27,NN)=SIJKL(27,NN)+DUMSXZ
C-    SIJKL(32,NN)=SIJKL(32,NN)+DUMSYX
C-    SIJKL(33,NN)=SIJKL(33,NN)+DUMSYY
C-    SIJKL(34,NN)=SIJKL(34,NN)+DUMSYZ
C-    SIJKL(40,NN)=SIJKL(40,NN)+DUMSZX
C-    SIJKL(41,NN)=SIJKL(41,NN)+DUMSZY
C-    SIJKL(42,NN)=SIJKL(42,NN)+DUMSZZ
 1670 IF(SKIPL) GO TO 1690
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1680 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FJXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NY)*X(IGR)
 1680 CONTINUE
      SD(49)=SD(49)+DAB(NN)*DUMSXX
      SD(50)=SD(50)+DAB(NN)*DUMSXY
      SD(51)=SD(51)+DAB(NN)*DUMSXZ
      SD(59)=SD(59)+DAB(NN)*DUMSYX
      SD(60)=SD(60)+DAB(NN)*DUMSYY
      SD(61)=SD(61)+DAB(NN)*DUMSYZ
      SD(70)=SD(70)+DAB(NN)*DUMSZX
      SD(71)=SD(71)+DAB(NN)*DUMSZY
      SD(72)=SD(72)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1690
C-    SIJKL(49,NN)=SIJKL(49,NN)+DUMSXX
C-    SIJKL(50,NN)=SIJKL(50,NN)+DUMSXY
C-    SIJKL(51,NN)=SIJKL(51,NN)+DUMSXZ
C-    SIJKL(59,NN)=SIJKL(59,NN)+DUMSYX
C-    SIJKL(60,NN)=SIJKL(60,NN)+DUMSYY
C-    SIJKL(61,NN)=SIJKL(61,NN)+DUMSYZ
C-    SIJKL(70,NN)=SIJKL(70,NN)+DUMSZX
C-    SIJKL(71,NN)=SIJKL(71,NN)+DUMSZY
C-    SIJKL(72,NN)=SIJKL(72,NN)+DUMSZZ
 1690 IF(SKIPK) GO TO 1800
      IF(SKIPL) GO TO 1800
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1700 IGR=1,NG*NR
      DUMSXX=DUMSXX+SKLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SKLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SKLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FKXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FKXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FKXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FKXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NY)*X(IGR)
 1700 CONTINUE
      SD(52)=SD(52)+DAB(NN)*DUMSXX
      SD(53)=SD(53)+DAB(NN)*DUMSXY
      SD(54)=SD(54)+DAB(NN)*DUMSXZ
      SD(62)=SD(62)+DAB(NN)*DUMSYX
      SD(63)=SD(63)+DAB(NN)*DUMSYY
      SD(64)=SD(64)+DAB(NN)*DUMSYZ
      SD(73)=SD(73)+DAB(NN)*DUMSZX
      SD(74)=SD(74)+DAB(NN)*DUMSZY
      SD(75)=SD(75)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1800
C-    SIJKL(52,NN)=SIJKL(52,NN)+DUMSXX
C-    SIJKL(53,NN)=SIJKL(53,NN)+DUMSXY
C-    SIJKL(54,NN)=SIJKL(54,NN)+DUMSXZ
C-    SIJKL(62,NN)=SIJKL(62,NN)+DUMSYX
C-    SIJKL(63,NN)=SIJKL(63,NN)+DUMSYY
C-    SIJKL(64,NN)=SIJKL(64,NN)+DUMSYZ
C-    SIJKL(73,NN)=SIJKL(73,NN)+DUMSZX
C-    SIJKL(74,NN)=SIJKL(74,NN)+DUMSZY
C-    SIJKL(75,NN)=SIJKL(75,NN)+DUMSZZ
C
 1800 CONTINUE
C
 1810 CONTINUE
 1820 CONTINUE
 1830 CONTINUE
 1840 CONTINUE
C
C     ----- ZEROTH ORDER INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 3640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 3630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 3110 IGR=1,NG*NR
 3110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 3120 IGR=1,NG*NR
 3120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 3130 IGR=1,NG*NR
 3130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 3620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 3210 IGR=1,NG*NR
 3210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 3220 IGR=1,NG*NR
 3220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 3230 IGR=1,NG*NR
 3230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 3610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 3310 IGR=1,NG*NR
 3310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 3320 IGR=1,NG*NR
 3320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 3330 IGR=1,NG*NR
 3330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 3400 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 3400    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 3410 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 3410    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      DUM=ZERO
      DO 3420 IGR=1,NG*NR
 3420 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
 3610 CONTINUE
 3620 CONTINUE
 3630 CONTINUE
 3640 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_DDSPDS(NG,NR,IJKLG,GIJKL,FIJKL,SIJKL,XYZ,
     1 FIXYZ,FJXYZ,FKXYZ,FLXYZ,SIIXYZ,SJJXYZ,SKKXYZ,SLLXYZ,
     2 SIJXYZ,SIKXYZ,SILXYZ,SJKXYZ,SJLXYZ,SKLXYZ,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IS,JS,KS,LS
      LOGICAL IJS,IJKS,IJKLS
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1 LGXYZ(4,35)
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION    XYZ(NG*NR,1)
      DIMENSION  FIXYZ(NG*NR,1), FJXYZ(NG*NR,1)
      DIMENSION  FKXYZ(NG*NR,1), FLXYZ(NG*NR,1)
      DIMENSION SIIXYZ(NG*NR,1),SJJXYZ(NG*NR,1)
      DIMENSION SKKXYZ(NG*NR,1),SLLXYZ(NG*NR,1)
      DIMENSION SIJXYZ(NG*NR,1),SIKXYZ(NG*NR,1)
      DIMENSION SILXYZ(NG*NR,1),SJKXYZ(NG*NR,1)
      DIMENSION SJLXYZ(NG*NR,1),SKLXYZ(NG*NR,1)
      DIMENSION XY(NG*NR),XZ(NG*NR),YZ(NG*NR)
      DIMENSION  X(NG*NR), Y(NG*NR), Z(NG*NR)
      DIMENSION SJ(NG*NR),SK(NG*NR),SL(NG*NR)
      DIMENSION DIJSI(NG*NR),DIJSJ(NG*NR),DKLSK(NG*NR),DKLSL(NG*NR)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
C
C     ----- FIRST DERIVATIVE INTEGRALS + GRADIENT + HESSIAN -----
C
      IF(SPIJKL) GO TO 1000
C
C     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
C
      IJKLN=0
      DO 440 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 430 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 420 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 410 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 110 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  110 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
C
      IF(SKIPI) GO TO 130
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 120 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SIIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SIIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SIIXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FIXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FIXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FIXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
  120 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
      SD( 1)=SD( 1)+DAB(NN)*DUMSXX
      SD( 3)=SD( 3)+DAB(NN)*DUMSYY
      SD( 6)=SD( 6)+DAB(NN)*DUMSZZ
      SD( 2)=SD( 2)+DAB(NN)*DUMSXY
      SD( 4)=SD( 4)+DAB(NN)*DUMSXZ
      SD( 5)=SD( 5)+DAB(NN)*DUMSYZ
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 130
C-    SIJKL( 1,NN)=SIJKL( 1,NN)+DUMSXX
C-    SIJKL( 3,NN)=SIJKL( 3,NN)+DUMSYY
C-    SIJKL( 6,NN)=SIJKL( 6,NN)+DUMSZZ
C-    SIJKL( 2,NN)=SIJKL( 2,NN)+DUMSXY
C-    SIJKL( 4,NN)=SIJKL( 4,NN)+DUMSXZ
C-    SIJKL( 5,NN)=SIJKL( 5,NN)+DUMSYZ
  130 IF(SKIPJ) GO TO 150
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 140 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SJJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SJJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SJJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
  140 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
      SD(10)=SD(10)+DAB(NN)*DUMSXX
      SD(15)=SD(15)+DAB(NN)*DUMSYY
      SD(21)=SD(21)+DAB(NN)*DUMSZZ
      SD(14)=SD(14)+DAB(NN)*DUMSXY
      SD(19)=SD(19)+DAB(NN)*DUMSXZ
      SD(20)=SD(20)+DAB(NN)*DUMSYZ
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 150
C-    SIJKL(10,NN)=SIJKL(10,NN)+DUMSXX
C-    SIJKL(15,NN)=SIJKL(15,NN)+DUMSYY
C-    SIJKL(21,NN)=SIJKL(21,NN)+DUMSZZ
C-    SIJKL(14,NN)=SIJKL(14,NN)+DUMSXY
C-    SIJKL(19,NN)=SIJKL(19,NN)+DUMSXZ
C-    SIJKL(20,NN)=SIJKL(20,NN)+DUMSYZ
  150 IF(SKIPK) GO TO 170
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 160 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SKKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SKKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SKKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FKXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FKXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FKXYZ(IGR,NZ)*XYZ(IGR,NX)
  160 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
      SD(28)=SD(28)+DAB(NN)*DUMSXX
      SD(36)=SD(36)+DAB(NN)*DUMSYY
      SD(45)=SD(45)+DAB(NN)*DUMSZZ
      SD(35)=SD(35)+DAB(NN)*DUMSXY
      SD(43)=SD(43)+DAB(NN)*DUMSXZ
      SD(44)=SD(44)+DAB(NN)*DUMSYZ
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 170
C-    SIJKL(28,NN)=SIJKL(28,NN)+DUMSXX
C-    SIJKL(36,NN)=SIJKL(36,NN)+DUMSYY
C-    SIJKL(45,NN)=SIJKL(45,NN)+DUMSZZ
C-    SIJKL(35,NN)=SIJKL(35,NN)+DUMSXY
C-    SIJKL(43,NN)=SIJKL(43,NN)+DUMSXZ
C-    SIJKL(44,NN)=SIJKL(44,NN)+DUMSYZ
  170 IF(SKIPL) GO TO 190
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 180 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SLLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SLLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SLLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FLXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FLXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FLXYZ(IGR,NZ)*XYZ(IGR,NX)
  180 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
      SD(55)=SD(55)+DAB(NN)*DUMSXX
      SD(66)=SD(66)+DAB(NN)*DUMSYY
      SD(78)=SD(78)+DAB(NN)*DUMSZZ
      SD(65)=SD(65)+DAB(NN)*DUMSXY
      SD(76)=SD(76)+DAB(NN)*DUMSXZ
      SD(77)=SD(77)+DAB(NN)*DUMSYZ
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 190
C-    SIJKL(55,NN)=SIJKL(55,NN)+DUMSXX
C-    SIJKL(66,NN)=SIJKL(66,NN)+DUMSYY
C-    SIJKL(78,NN)=SIJKL(78,NN)+DUMSZZ
C-    SIJKL(65,NN)=SIJKL(65,NN)+DUMSXY
C-    SIJKL(76,NN)=SIJKL(76,NN)+DUMSXZ
C-    SIJKL(77,NN)=SIJKL(77,NN)+DUMSYZ
  190 IF(SKIPI) GO TO 250
      IF(SKIPJ) GO TO 210
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 200 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIJXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIJXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FJXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  200 CONTINUE
      SD( 7)=SD( 7)+DAB(NN)*DUMSXX
      SD( 8)=SD( 8)+DAB(NN)*DUMSXY
      SD( 9)=SD( 9)+DAB(NN)*DUMSXZ
      SD(11)=SD(11)+DAB(NN)*DUMSYX
      SD(12)=SD(12)+DAB(NN)*DUMSYY
      SD(13)=SD(13)+DAB(NN)*DUMSYZ
      SD(16)=SD(16)+DAB(NN)*DUMSZX
      SD(17)=SD(17)+DAB(NN)*DUMSZY
      SD(18)=SD(18)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 210
C-    SIJKL( 7,NN)=SIJKL( 7,NN)+DUMSXX
C-    SIJKL( 8,NN)=SIJKL( 8,NN)+DUMSXY
C-    SIJKL( 9,NN)=SIJKL( 9,NN)+DUMSXZ
C-    SIJKL(11,NN)=SIJKL(11,NN)+DUMSYX
C-    SIJKL(12,NN)=SIJKL(12,NN)+DUMSYY
C-    SIJKL(13,NN)=SIJKL(13,NN)+DUMSYZ
C-    SIJKL(16,NN)=SIJKL(16,NN)+DUMSZX
C-    SIJKL(17,NN)=SIJKL(17,NN)+DUMSZY
C-    SIJKL(18,NN)=SIJKL(18,NN)+DUMSZZ
  210 IF(SKIPK) GO TO 230
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 220 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  220 CONTINUE
      SD(22)=SD(22)+DAB(NN)*DUMSXX
      SD(23)=SD(23)+DAB(NN)*DUMSXY
      SD(24)=SD(24)+DAB(NN)*DUMSXZ
      SD(29)=SD(29)+DAB(NN)*DUMSYX
      SD(30)=SD(30)+DAB(NN)*DUMSYY
      SD(31)=SD(31)+DAB(NN)*DUMSYZ
      SD(37)=SD(37)+DAB(NN)*DUMSZX
      SD(38)=SD(38)+DAB(NN)*DUMSZY
      SD(39)=SD(39)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 230
C-    SIJKL(22,NN)=SIJKL(22,NN)+DUMSXX
C-    SIJKL(23,NN)=SIJKL(23,NN)+DUMSXY
C-    SIJKL(24,NN)=SIJKL(24,NN)+DUMSXZ
C-    SIJKL(29,NN)=SIJKL(29,NN)+DUMSYX
C-    SIJKL(30,NN)=SIJKL(30,NN)+DUMSYY
C-    SIJKL(31,NN)=SIJKL(31,NN)+DUMSYZ
C-    SIJKL(37,NN)=SIJKL(37,NN)+DUMSZX
C-    SIJKL(38,NN)=SIJKL(38,NN)+DUMSZY
C-    SIJKL(39,NN)=SIJKL(39,NN)+DUMSZZ
  230 IF(SKIPL) GO TO 250
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 240 IGR=1,NG*NR
      DUMSXX=DUMSXX+SILXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SILXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SILXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FIXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FIXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FIXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FIXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NY)*XYZ(IGR,NX)
  240 CONTINUE
      SD(46)=SD(46)+DAB(NN)*DUMSXX
      SD(47)=SD(47)+DAB(NN)*DUMSXY
      SD(48)=SD(48)+DAB(NN)*DUMSXZ
      SD(56)=SD(56)+DAB(NN)*DUMSYX
      SD(57)=SD(57)+DAB(NN)*DUMSYY
      SD(58)=SD(58)+DAB(NN)*DUMSYZ
      SD(67)=SD(67)+DAB(NN)*DUMSZX
      SD(68)=SD(68)+DAB(NN)*DUMSZY
      SD(69)=SD(69)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 250
C-    SIJKL(46,NN)=SIJKL(46,NN)+DUMSXX
C-    SIJKL(47,NN)=SIJKL(47,NN)+DUMSXY
C-    SIJKL(48,NN)=SIJKL(48,NN)+DUMSXZ
C-    SIJKL(56,NN)=SIJKL(56,NN)+DUMSYX
C-    SIJKL(57,NN)=SIJKL(57,NN)+DUMSYY
C-    SIJKL(58,NN)=SIJKL(58,NN)+DUMSYZ
C-    SIJKL(67,NN)=SIJKL(67,NN)+DUMSZX
C-    SIJKL(68,NN)=SIJKL(68,NN)+DUMSZY
C-    SIJKL(69,NN)=SIJKL(69,NN)+DUMSZZ
  250 IF(SKIPJ) GO TO 290
      IF(SKIPK) GO TO 270
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 260 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FJXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NY)*XYZ(IGR,NX)
  260 CONTINUE
      SD(25)=SD(25)+DAB(NN)*DUMSXX
      SD(26)=SD(26)+DAB(NN)*DUMSXY
      SD(27)=SD(27)+DAB(NN)*DUMSXZ
      SD(32)=SD(32)+DAB(NN)*DUMSYX
      SD(33)=SD(33)+DAB(NN)*DUMSYY
      SD(34)=SD(34)+DAB(NN)*DUMSYZ
      SD(40)=SD(40)+DAB(NN)*DUMSZX
      SD(41)=SD(41)+DAB(NN)*DUMSZY
      SD(42)=SD(42)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 270
C-    SIJKL(25,NN)=SIJKL(25,NN)+DUMSXX
C-    SIJKL(26,NN)=SIJKL(26,NN)+DUMSXY
C-    SIJKL(27,NN)=SIJKL(27,NN)+DUMSXZ
C-    SIJKL(32,NN)=SIJKL(32,NN)+DUMSYX
C-    SIJKL(33,NN)=SIJKL(33,NN)+DUMSYY
C-    SIJKL(34,NN)=SIJKL(34,NN)+DUMSYZ
C-    SIJKL(40,NN)=SIJKL(40,NN)+DUMSZX
C-    SIJKL(41,NN)=SIJKL(41,NN)+DUMSZY
C-    SIJKL(42,NN)=SIJKL(42,NN)+DUMSZZ
  270 IF(SKIPL) GO TO 290
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 280 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FJXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FJXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FJXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FJXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NY)*XYZ(IGR,NX)
  280 CONTINUE
      SD(49)=SD(49)+DAB(NN)*DUMSXX
      SD(50)=SD(50)+DAB(NN)*DUMSXY
      SD(51)=SD(51)+DAB(NN)*DUMSXZ
      SD(59)=SD(59)+DAB(NN)*DUMSYX
      SD(60)=SD(60)+DAB(NN)*DUMSYY
      SD(61)=SD(61)+DAB(NN)*DUMSYZ
      SD(70)=SD(70)+DAB(NN)*DUMSZX
      SD(71)=SD(71)+DAB(NN)*DUMSZY
      SD(72)=SD(72)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 290
C-    SIJKL(49,NN)=SIJKL(49,NN)+DUMSXX
C-    SIJKL(50,NN)=SIJKL(50,NN)+DUMSXY
C-    SIJKL(51,NN)=SIJKL(51,NN)+DUMSXZ
C-    SIJKL(59,NN)=SIJKL(59,NN)+DUMSYX
C-    SIJKL(60,NN)=SIJKL(60,NN)+DUMSYY
C-    SIJKL(61,NN)=SIJKL(61,NN)+DUMSYZ
C-    SIJKL(70,NN)=SIJKL(70,NN)+DUMSZX
C-    SIJKL(71,NN)=SIJKL(71,NN)+DUMSZY
C-    SIJKL(72,NN)=SIJKL(72,NN)+DUMSZZ
  290 IF(SKIPK) GO TO 400
      IF(SKIPL) GO TO 400
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 300 IGR=1,NG*NR
      DUMSXX=DUMSXX+SKLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SKLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SKLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FKXYZ(IGR,NY)*XYZ(IGR,NZ)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FKXYZ(IGR,NZ)*XYZ(IGR,NY)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FKXYZ(IGR,NX)*XYZ(IGR,NZ)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FKXYZ(IGR,NZ)*XYZ(IGR,NX)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NX)*XYZ(IGR,NY)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NY)*XYZ(IGR,NX)
  300 CONTINUE
      SD(52)=SD(52)+DAB(NN)*DUMSXX
      SD(53)=SD(53)+DAB(NN)*DUMSXY
      SD(54)=SD(54)+DAB(NN)*DUMSXZ
      SD(62)=SD(62)+DAB(NN)*DUMSYX
      SD(63)=SD(63)+DAB(NN)*DUMSYY
      SD(64)=SD(64)+DAB(NN)*DUMSYZ
      SD(73)=SD(73)+DAB(NN)*DUMSZX
      SD(74)=SD(74)+DAB(NN)*DUMSZY
      SD(75)=SD(75)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 400
C-    SIJKL(52,NN)=SIJKL(52,NN)+DUMSXX
C-    SIJKL(53,NN)=SIJKL(53,NN)+DUMSXY
C-    SIJKL(54,NN)=SIJKL(54,NN)+DUMSXZ
C-    SIJKL(62,NN)=SIJKL(62,NN)+DUMSYX
C-    SIJKL(63,NN)=SIJKL(63,NN)+DUMSYY
C-    SIJKL(64,NN)=SIJKL(64,NN)+DUMSYZ
C-    SIJKL(73,NN)=SIJKL(73,NN)+DUMSZX
C-    SIJKL(74,NN)=SIJKL(74,NN)+DUMSZY
C-    SIJKL(75,NN)=SIJKL(75,NN)+DUMSZZ
C
  400 CONTINUE
C
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
  440 CONTINUE
C
C     ----- ZEROTH ORDER INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 940 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 930 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 920 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 910 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 710 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  710 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      DUM=ZERO
      DO 720 IGR=1,NG*NR
  720 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
  910 CONTINUE
  920 CONTINUE
  930 CONTINUE
  940 CONTINUE
C
      RETURN
C
C     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
C
 1000 CONTINUE
C
C     ----- FIRST DERIVATIVE INTEGRALS + GRADIENT + HESSIAN -----
C
      IJKLN=0
      DO 1840 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 1830 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 1110 IGR=1,NG*NR
 1110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 1120 IGR=1,NG*NR
 1120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 1130 IGR=1,NG*NR
 1130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 1820 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 1210 IGR=1,NG*NR
 1210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 1220 IGR=1,NG*NR
 1220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 1230 IGR=1,NG*NR
 1230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 1810 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 1310 IGR=1,NG*NR
 1310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 1320 IGR=1,NG*NR
 1320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 1330 IGR=1,NG*NR
 1330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 1500 IGR=1,NG*NR
          X(IGR)=XYZ(IGR,NX)*SL(IGR)
          Y(IGR)=XYZ(IGR,NY)*SL(IGR)
          Z(IGR)=XYZ(IGR,NZ)*SL(IGR)
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 1500    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 1510 IGR=1,NG*NR
          X(IGR)=XYZ(IGR,NX)
          Y(IGR)=XYZ(IGR,NY)
          Z(IGR)=XYZ(IGR,NZ)
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 1510    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      IF(SKIPI) GO TO 1530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SIIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SIIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SIIXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FIXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FIXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FIXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
 1520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
      SD( 1)=SD( 1)+DAB(NN)*DUMSXX
      SD( 3)=SD( 3)+DAB(NN)*DUMSYY
      SD( 6)=SD( 6)+DAB(NN)*DUMSZZ
      SD( 2)=SD( 2)+DAB(NN)*DUMSXY
      SD( 4)=SD( 4)+DAB(NN)*DUMSXZ
      SD( 5)=SD( 5)+DAB(NN)*DUMSYZ
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1530
C-    SIJKL( 1,NN)=SIJKL( 1,NN)+DUMSXX
C-    SIJKL( 3,NN)=SIJKL( 3,NN)+DUMSYY
C-    SIJKL( 6,NN)=SIJKL( 6,NN)+DUMSZZ
C-    SIJKL( 2,NN)=SIJKL( 2,NN)+DUMSXY
C-    SIJKL( 4,NN)=SIJKL( 4,NN)+DUMSXZ
C-    SIJKL( 5,NN)=SIJKL( 5,NN)+DUMSYZ
 1530 IF(SKIPJ) GO TO 1550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SJJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SJJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SJJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
 1540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
      SD(10)=SD(10)+DAB(NN)*DUMSXX
      SD(15)=SD(15)+DAB(NN)*DUMSYY
      SD(21)=SD(21)+DAB(NN)*DUMSZZ
      SD(14)=SD(14)+DAB(NN)*DUMSXY
      SD(19)=SD(19)+DAB(NN)*DUMSXZ
      SD(20)=SD(20)+DAB(NN)*DUMSYZ
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1550
C-    SIJKL(10,NN)=SIJKL(10,NN)+DUMSXX
C-    SIJKL(15,NN)=SIJKL(15,NN)+DUMSYY
C-    SIJKL(21,NN)=SIJKL(21,NN)+DUMSZZ
C-    SIJKL(14,NN)=SIJKL(14,NN)+DUMSXY
C-    SIJKL(19,NN)=SIJKL(19,NN)+DUMSXZ
C-    SIJKL(20,NN)=SIJKL(20,NN)+DUMSYZ
 1550 IF(SKIPK) GO TO 1570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SKKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SKKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SKKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FKXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FKXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FKXYZ(IGR,NZ)*X(IGR)
 1560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
      SD(28)=SD(28)+DAB(NN)*DUMSXX
      SD(36)=SD(36)+DAB(NN)*DUMSYY
      SD(45)=SD(45)+DAB(NN)*DUMSZZ
      SD(35)=SD(35)+DAB(NN)*DUMSXY
      SD(43)=SD(43)+DAB(NN)*DUMSXZ
      SD(44)=SD(44)+DAB(NN)*DUMSYZ
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1570
C-    SIJKL(28,NN)=SIJKL(28,NN)+DUMSXX
C-    SIJKL(36,NN)=SIJKL(36,NN)+DUMSYY
C-    SIJKL(45,NN)=SIJKL(45,NN)+DUMSZZ
C-    SIJKL(35,NN)=SIJKL(35,NN)+DUMSXY
C-    SIJKL(43,NN)=SIJKL(43,NN)+DUMSXZ
C-    SIJKL(44,NN)=SIJKL(44,NN)+DUMSYZ
 1570 IF(SKIPL) GO TO 1590
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYZ=ZERO
      DO 1580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMSXX=DUMSXX+SLLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMSYY=DUMSYY+SLLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
      DUMSZZ=DUMSZZ+SLLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FLXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FLXYZ(IGR,NZ)*Y(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FLXYZ(IGR,NZ)*X(IGR)
 1580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
      SD(55)=SD(55)+DAB(NN)*DUMSXX
      SD(66)=SD(66)+DAB(NN)*DUMSYY
      SD(78)=SD(78)+DAB(NN)*DUMSZZ
      SD(65)=SD(65)+DAB(NN)*DUMSXY
      SD(76)=SD(76)+DAB(NN)*DUMSXZ
      SD(77)=SD(77)+DAB(NN)*DUMSYZ
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
C-    IF(.NOT.OUT) GO TO 1590
C-    SIJKL(55,NN)=SIJKL(55,NN)+DUMSXX
C-    SIJKL(66,NN)=SIJKL(66,NN)+DUMSYY
C-    SIJKL(78,NN)=SIJKL(78,NN)+DUMSZZ
C-    SIJKL(65,NN)=SIJKL(65,NN)+DUMSXY
C-    SIJKL(76,NN)=SIJKL(76,NN)+DUMSXZ
C-    SIJKL(77,NN)=SIJKL(77,NN)+DUMSYZ
 1590 IF(SKIPI) GO TO 1650
      IF(SKIPJ) GO TO 1610
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1600 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIJXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIJXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIJXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FJXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FJXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FJXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FJXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FJXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1600 CONTINUE
      SD( 7)=SD( 7)+DAB(NN)*DUMSXX
      SD( 8)=SD( 8)+DAB(NN)*DUMSXY
      SD( 9)=SD( 9)+DAB(NN)*DUMSXZ
      SD(11)=SD(11)+DAB(NN)*DUMSYX
      SD(12)=SD(12)+DAB(NN)*DUMSYY
      SD(13)=SD(13)+DAB(NN)*DUMSYZ
      SD(16)=SD(16)+DAB(NN)*DUMSZX
      SD(17)=SD(17)+DAB(NN)*DUMSZY
      SD(18)=SD(18)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1610
C-    SIJKL( 7,NN)=SIJKL( 7,NN)+DUMSXX
C-    SIJKL( 8,NN)=SIJKL( 8,NN)+DUMSXY
C-    SIJKL( 9,NN)=SIJKL( 9,NN)+DUMSXZ
C-    SIJKL(11,NN)=SIJKL(11,NN)+DUMSYX
C-    SIJKL(12,NN)=SIJKL(12,NN)+DUMSYY
C-    SIJKL(13,NN)=SIJKL(13,NN)+DUMSYZ
C-    SIJKL(16,NN)=SIJKL(16,NN)+DUMSZX
C-    SIJKL(17,NN)=SIJKL(17,NN)+DUMSZY
C-    SIJKL(18,NN)=SIJKL(18,NN)+DUMSZZ
 1610 IF(SKIPK) GO TO 1630
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1620 IGR=1,NG*NR
      DUMSXX=DUMSXX+SIKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SIKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SIKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1620 CONTINUE
      SD(22)=SD(22)+DAB(NN)*DUMSXX
      SD(23)=SD(23)+DAB(NN)*DUMSXY
      SD(24)=SD(24)+DAB(NN)*DUMSXZ
      SD(29)=SD(29)+DAB(NN)*DUMSYX
      SD(30)=SD(30)+DAB(NN)*DUMSYY
      SD(31)=SD(31)+DAB(NN)*DUMSYZ
      SD(37)=SD(37)+DAB(NN)*DUMSZX
      SD(38)=SD(38)+DAB(NN)*DUMSZY
      SD(39)=SD(39)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1630
C-    SIJKL(22,NN)=SIJKL(22,NN)+DUMSXX
C-    SIJKL(23,NN)=SIJKL(23,NN)+DUMSXY
C-    SIJKL(24,NN)=SIJKL(24,NN)+DUMSXZ
C-    SIJKL(29,NN)=SIJKL(29,NN)+DUMSYX
C-    SIJKL(30,NN)=SIJKL(30,NN)+DUMSYY
C-    SIJKL(31,NN)=SIJKL(31,NN)+DUMSYZ
C-    SIJKL(37,NN)=SIJKL(37,NN)+DUMSZX
C-    SIJKL(38,NN)=SIJKL(38,NN)+DUMSZY
C-    SIJKL(39,NN)=SIJKL(39,NN)+DUMSZZ
 1630 IF(SKIPL) GO TO 1650
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1640 IGR=1,NG*NR
      DUMSXX=DUMSXX+SILXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SILXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SILXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FIXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FIXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FIXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FIXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FIXYZ(IGR,NY)*X(IGR)
 1640 CONTINUE
      SD(46)=SD(46)+DAB(NN)*DUMSXX
      SD(47)=SD(47)+DAB(NN)*DUMSXY
      SD(48)=SD(48)+DAB(NN)*DUMSXZ
      SD(56)=SD(56)+DAB(NN)*DUMSYX
      SD(57)=SD(57)+DAB(NN)*DUMSYY
      SD(58)=SD(58)+DAB(NN)*DUMSYZ
      SD(67)=SD(67)+DAB(NN)*DUMSZX
      SD(68)=SD(68)+DAB(NN)*DUMSZY
      SD(69)=SD(69)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1650
C-    SIJKL(46,NN)=SIJKL(46,NN)+DUMSXX
C-    SIJKL(47,NN)=SIJKL(47,NN)+DUMSXY
C-    SIJKL(48,NN)=SIJKL(48,NN)+DUMSXZ
C-    SIJKL(56,NN)=SIJKL(56,NN)+DUMSYX
C-    SIJKL(57,NN)=SIJKL(57,NN)+DUMSYY
C-    SIJKL(58,NN)=SIJKL(58,NN)+DUMSYZ
C-    SIJKL(67,NN)=SIJKL(67,NN)+DUMSZX
C-    SIJKL(68,NN)=SIJKL(68,NN)+DUMSZY
C-    SIJKL(69,NN)=SIJKL(69,NN)+DUMSZZ
 1650 IF(SKIPJ) GO TO 1690
      IF(SKIPK) GO TO 1670
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1660 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJKXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJKXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJKXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FKXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FKXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FKXYZ(IGR,NY)*FJXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FKXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FKXYZ(IGR,NZ)*FJXYZ(IGR,NY)*X(IGR)
 1660 CONTINUE
      SD(25)=SD(25)+DAB(NN)*DUMSXX
      SD(26)=SD(26)+DAB(NN)*DUMSXY
      SD(27)=SD(27)+DAB(NN)*DUMSXZ
      SD(32)=SD(32)+DAB(NN)*DUMSYX
      SD(33)=SD(33)+DAB(NN)*DUMSYY
      SD(34)=SD(34)+DAB(NN)*DUMSYZ
      SD(40)=SD(40)+DAB(NN)*DUMSZX
      SD(41)=SD(41)+DAB(NN)*DUMSZY
      SD(42)=SD(42)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1670
C-    SIJKL(25,NN)=SIJKL(25,NN)+DUMSXX
C-    SIJKL(26,NN)=SIJKL(26,NN)+DUMSXY
C-    SIJKL(27,NN)=SIJKL(27,NN)+DUMSXZ
C-    SIJKL(32,NN)=SIJKL(32,NN)+DUMSYX
C-    SIJKL(33,NN)=SIJKL(33,NN)+DUMSYY
C-    SIJKL(34,NN)=SIJKL(34,NN)+DUMSYZ
C-    SIJKL(40,NN)=SIJKL(40,NN)+DUMSZX
C-    SIJKL(41,NN)=SIJKL(41,NN)+DUMSZY
C-    SIJKL(42,NN)=SIJKL(42,NN)+DUMSZZ
 1670 IF(SKIPL) GO TO 1690
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1680 IGR=1,NG*NR
      DUMSXX=DUMSXX+SJLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SJLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SJLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FJXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FJXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FJXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FJXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FJXYZ(IGR,NY)*X(IGR)
 1680 CONTINUE
      SD(49)=SD(49)+DAB(NN)*DUMSXX
      SD(50)=SD(50)+DAB(NN)*DUMSXY
      SD(51)=SD(51)+DAB(NN)*DUMSXZ
      SD(59)=SD(59)+DAB(NN)*DUMSYX
      SD(60)=SD(60)+DAB(NN)*DUMSYY
      SD(61)=SD(61)+DAB(NN)*DUMSYZ
      SD(70)=SD(70)+DAB(NN)*DUMSZX
      SD(71)=SD(71)+DAB(NN)*DUMSZY
      SD(72)=SD(72)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1690
C-    SIJKL(49,NN)=SIJKL(49,NN)+DUMSXX
C-    SIJKL(50,NN)=SIJKL(50,NN)+DUMSXY
C-    SIJKL(51,NN)=SIJKL(51,NN)+DUMSXZ
C-    SIJKL(59,NN)=SIJKL(59,NN)+DUMSYX
C-    SIJKL(60,NN)=SIJKL(60,NN)+DUMSYY
C-    SIJKL(61,NN)=SIJKL(61,NN)+DUMSYZ
C-    SIJKL(70,NN)=SIJKL(70,NN)+DUMSZX
C-    SIJKL(71,NN)=SIJKL(71,NN)+DUMSZY
C-    SIJKL(72,NN)=SIJKL(72,NN)+DUMSZZ
 1690 IF(SKIPK) GO TO 1800
      IF(SKIPL) GO TO 1800
      DUMSXX=ZERO
      DUMSYY=ZERO
      DUMSZZ=ZERO
      DUMSXY=ZERO
      DUMSXZ=ZERO
      DUMSYX=ZERO
      DUMSYZ=ZERO
      DUMSZX=ZERO
      DUMSZY=ZERO
      DO 1700 IGR=1,NG*NR
      DUMSXX=DUMSXX+SKLXYZ(IGR,NX)*YZ(IGR)
      DUMSYY=DUMSYY+SKLXYZ(IGR,NY)*XZ(IGR)
      DUMSZZ=DUMSZZ+SKLXYZ(IGR,NZ)*XY(IGR)
      DUMSXY=DUMSXY+ FLXYZ(IGR,NX)*FKXYZ(IGR,NY)*Z(IGR)
      DUMSXZ=DUMSXZ+ FLXYZ(IGR,NX)*FKXYZ(IGR,NZ)*Y(IGR)
      DUMSYX=DUMSYX+ FLXYZ(IGR,NY)*FKXYZ(IGR,NX)*Z(IGR)
      DUMSYZ=DUMSYZ+ FLXYZ(IGR,NY)*FKXYZ(IGR,NZ)*X(IGR)
      DUMSZX=DUMSZX+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NX)*Y(IGR)
      DUMSZY=DUMSZY+ FLXYZ(IGR,NZ)*FKXYZ(IGR,NY)*X(IGR)
 1700 CONTINUE
      SD(52)=SD(52)+DAB(NN)*DUMSXX
      SD(53)=SD(53)+DAB(NN)*DUMSXY
      SD(54)=SD(54)+DAB(NN)*DUMSXZ
      SD(62)=SD(62)+DAB(NN)*DUMSYX
      SD(63)=SD(63)+DAB(NN)*DUMSYY
      SD(64)=SD(64)+DAB(NN)*DUMSYZ
      SD(73)=SD(73)+DAB(NN)*DUMSZX
      SD(74)=SD(74)+DAB(NN)*DUMSZY
      SD(75)=SD(75)+DAB(NN)*DUMSZZ
C-    IF(.NOT.OUT) GO TO 1800
C-    SIJKL(52,NN)=SIJKL(52,NN)+DUMSXX
C-    SIJKL(53,NN)=SIJKL(53,NN)+DUMSXY
C-    SIJKL(54,NN)=SIJKL(54,NN)+DUMSXZ
C-    SIJKL(62,NN)=SIJKL(62,NN)+DUMSYX
C-    SIJKL(63,NN)=SIJKL(63,NN)+DUMSYY
C-    SIJKL(64,NN)=SIJKL(64,NN)+DUMSYZ
C-    SIJKL(73,NN)=SIJKL(73,NN)+DUMSZX
C-    SIJKL(74,NN)=SIJKL(74,NN)+DUMSZY
C-    SIJKL(75,NN)=SIJKL(75,NN)+DUMSZZ
C
 1800 CONTINUE
C
 1810 CONTINUE
 1820 CONTINUE
 1830 CONTINUE
 1840 CONTINUE
C
C     ----- ZEROTH ORDER INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 3640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 3630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 3110 IGR=1,NG*NR
 3110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 3120 IGR=1,NG*NR
 3120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 3130 IGR=1,NG*NR
 3130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 3620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 3210 IGR=1,NG*NR
 3210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 3220 IGR=1,NG*NR
 3220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 3230 IGR=1,NG*NR
 3230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 3610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 3310 IGR=1,NG*NR
 3310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 3320 IGR=1,NG*NR
 3320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 3330 IGR=1,NG*NR
 3330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 3400 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 3400    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 3410 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 3410    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      DUM=ZERO
      DO 3420 IGR=1,NG*NR
 3420 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
 3610 CONTINUE
 3620 CONTINUE
 3630 CONTINUE
 3640 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JDDINV(FIJKL,SIJKL,IJKLG,FC,NC,DF,LDF,NBF,DEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXNBF =2048)
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERINV/INVTYP
      COMMON/HND_DERSHL/FD(3,4),SD(78)
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION FIJKL(3,4,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION FC(NC,1),DF(LDF,1),DEN(1)
      DIMENSION KAT(4),KC(4)
      DIMENSION DD(12,12),DS(3,4,3,4)
      EQUIVALENCE (DD(1,1),DS(1,1,1,1))
      EQUIVALENCE (KAT(1),IIAT)

C
C     ----- TRANSLATIONAL INVARIANCE FOR GRADIENT   ELEMENTS  -----
C     ----- TRANSLATIONAL INVARIANCE FOR DERIVATIVE INTEGRALS -----
C     ----- TRANSLATIONAL INVARIANCE FOR HESSIAN    ELEMENTS  -----
C
      DO 3010 K=1,4
 3010 KC(K)=3*(KAT(K)-1)
      DO 3020 IC=1,12
      DO 3020 JC=1,IC
      DD(IC,JC)=SD((IA(IC)+JC))
      DD(JC,IC)=SD((IA(IC)+JC))
C 
 3020 CONTINUE
C
      GO TO (5000,3200,3300,3400,3400,3600,3700,3800,
     1       3900,3900,4100,4200,4300,4400,4500,5000),INVTYP
C
 3200 CONTINUE
      DO 3210 IXYZ=1,3
 3210 FD(IXYZ,1)=- FD(IXYZ,4)
      DO 3220 I=1,IJKL
      N=IJKLG(1,I)
      DO 3220 IXYZ=1,3
 3220 FIJKL(IXYZ,1,N)=- FIJKL(IXYZ,4,N)
      DO 3230 JXYZ=1,3
      DO 3230 IXYZ=1,3
      DS(IXYZ,1,JXYZ,4)=- DS(IXYZ,4,JXYZ,4)
      DS(JXYZ,4,IXYZ,1)=- DS(IXYZ,4,JXYZ,4)
 3230 DS(IXYZ,1,JXYZ,1)=  DS(IXYZ,4,JXYZ,4)
      GO TO 5000
C
 3300 CONTINUE
      DO 3310 IXYZ=1,3
 3310 FD(IXYZ,1)=- FD(IXYZ,3)
      DO 3320 I=1,IJKL
      N=IJKLG(1,I)
      DO 3320 IXYZ=1,3
 3320 FIJKL(IXYZ,1,N)=- FIJKL(IXYZ,3,N)
      DO 3330 JXYZ=1,3
      DO 3330 IXYZ=1,3
      DS(IXYZ,1,JXYZ,3)=- DS(IXYZ,3,JXYZ,3)
      DS(JXYZ,3,IXYZ,1)=- DS(IXYZ,3,JXYZ,3)
 3330 DS(IXYZ,1,JXYZ,1)=  DS(IXYZ,3,JXYZ,3)
      GO TO 5000
C
 3400 CONTINUE
      DO 3410 IXYZ=1,3
 3410 FD(IXYZ,1)=-(FD(IXYZ,3)+FD(IXYZ,4))
      DO 3420 I=1,IJKL
      N=IJKLG(1,I)
      DO 3420 IXYZ=1,3
 3420 FIJKL(IXYZ,1,N)=-(FIJKL(IXYZ,3,N)+FIJKL(IXYZ,4,N))
      DO 3430 JAT =3,4
      DO 3430 JXYZ=1,3
      DO 3430 IXYZ=1,3
      DS(IXYZ,  1,JXYZ,JAT)=-(DS(IXYZ,3,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
 3430 DS(JXYZ,JAT,IXYZ,  1)=-(DS(IXYZ,3,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
      DO 3440 JXYZ=1,3
      DO 3440 IXYZ=1,3
 3440 DS(IXYZ,  1,JXYZ,  1)=-(DS(IXYZ,1,JXYZ,  3)+DS(IXYZ,1,JXYZ,  4))
      GO TO 5000
C
 3600 CONTINUE
      DO 3610 IXYZ=1,3
 3610 FD(IXYZ,1)=- FD(IXYZ,2)
      DO 3620 I=1,IJKL
      N=IJKLG(1,I)
      DO 3620 IXYZ=1,3
 3620 FIJKL(IXYZ,1,N)=- FIJKL(IXYZ,2,N)
      DO 3630 JXYZ=1,3
      DO 3630 IXYZ=1,3
      DS(IXYZ,1,JXYZ,2)=- DS(IXYZ,2,JXYZ,2)
      DS(JXYZ,2,IXYZ,1)=- DS(IXYZ,2,JXYZ,2)
 3630 DS(IXYZ,1,JXYZ,1)=  DS(IXYZ,2,JXYZ,2)
      GO TO 5000
C
 3700 CONTINUE
      DO 3710 IXYZ=1,3
 3710 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      DO 3720 I=1,IJKL
      N=IJKLG(1,I)
      DO 3720 IXYZ=1,3
 3720 FIJKL(IXYZ,1,N)=-(FIJKL(IXYZ,2,N)+FIJKL(IXYZ,4,N))
      DO 3730 JAT =2,4,2
      DO 3730 JXYZ=1,3
      DO 3730 IXYZ=1,3
      DS(IXYZ,  1,JXYZ,JAT)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
 3730 DS(JXYZ,JAT,IXYZ,  1)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
      DO 3740 JXYZ=1,3
      DO 3740 IXYZ=1,3
 3740 DS(IXYZ,  1,JXYZ,  1)=-(DS(IXYZ,2,JXYZ,  1)+DS(IXYZ,4,JXYZ,  1))
      GO TO 5000
C
 3800 CONTINUE
      DO 3810 IXYZ=1,3
 3810 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      DO 3820 I=1,IJKL
      N=IJKLG(1,I)
      DO 3820 IXYZ=1,3
 3820 FIJKL(IXYZ,1,N)=-(FIJKL(IXYZ,2,N)+FIJKL(IXYZ,4,N))
      DO 3830 JAT =2,4,2
      DO 3830 JXYZ=1,3
      DO 3830 IXYZ=1,3
      DS(IXYZ,  1,JXYZ,JAT)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
 3830 DS(JXYZ,JAT,IXYZ,  1)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
      DO 3840 JXYZ=1,3
      DO 3840 IXYZ=1,3
 3840 DS(IXYZ,  1,JXYZ,  1)=-(DS(IXYZ,2,JXYZ,  1)+DS(IXYZ,4,JXYZ,  1))
      GO TO 5000
C
 3900 CONTINUE
      DO 3910 IXYZ=1,3
 3910 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,3))
      DO 3920 I=1,IJKL
      N=IJKLG(1,I)
      DO 3920 IXYZ=1,3
 3920 FIJKL(IXYZ,1,N)=-(FIJKL(IXYZ,2,N)+FIJKL(IXYZ,3,N))
      DO 3930 JAT =2,3
      DO 3930 JXYZ=1,3
      DO 3930 IXYZ=1,3
      DS(IXYZ,  1,JXYZ,JAT)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,3,JXYZ,JAT))
 3930 DS(JXYZ,JAT,IXYZ,  1)=-(DS(IXYZ,2,JXYZ,JAT)+DS(IXYZ,3,JXYZ,JAT))
      DO 3940 JXYZ=1,3
      DO 3940 IXYZ=1,3
 3940 DS(IXYZ,  1,JXYZ,  1)=-(DS(IXYZ,2,JXYZ,  1)+DS(IXYZ,3,JXYZ,  1))
      GO TO 5000
C
 4100 CONTINUE
      DO 4110 IXYZ=1,3
 4110 FD(IXYZ,2)=- FD(IXYZ,1)
      DO 4120 I=1,IJKL
      N=IJKLG(1,I)
      DO 4120 IXYZ=1,3
 4120 FIJKL(IXYZ,2,N)=- FIJKL(IXYZ,1,N)
      DO 4130 JXYZ=1,3
      DO 4130 IXYZ=1,3
      DS(IXYZ,1,JXYZ,2)=- DS(IXYZ,1,JXYZ,1)
      DS(JXYZ,2,IXYZ,1)=- DS(IXYZ,1,JXYZ,1)
 4130 DS(IXYZ,2,JXYZ,2)=  DS(IXYZ,1,JXYZ,1)
      GO TO 5000
C
 4200 CONTINUE
      DO 4210 IXYZ=1,3
 4210 FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,4))
      DO 4220 I=1,IJKL
      N=IJKLG(1,I)
      DO 4220 IXYZ=1,3
 4220 FIJKL(IXYZ,2,N)=-(FIJKL(IXYZ,1,N)+FIJKL(IXYZ,4,N))
      DO 4230 JAT =1,4,3
      DO 4230 JXYZ=1,3
      DO 4230 IXYZ=1,3
      DS(IXYZ,  2,JXYZ,JAT)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
 4230 DS(JXYZ,JAT,IXYZ,  2)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,4,JXYZ,JAT))
      DO 4240 JXYZ=1,3
      DO 4240 IXYZ=1,3
 4240 DS(IXYZ,  2,JXYZ,  2)=-(DS(IXYZ,1,JXYZ,  2)+DS(IXYZ,4,JXYZ,  2))
      GO TO 5000
C
 4300 CONTINUE
      DO 4310 IXYZ=1,3
 4310 FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,3))
      DO 4320 I=1,IJKL
      N=IJKLG(1,I)
      DO 4320 IXYZ=1,3
 4320 FIJKL(IXYZ,2,N)=-(FIJKL(IXYZ,1,N)+FIJKL(IXYZ,3,N))
      DO 4330 JAT =1,3,2
      DO 4330 JXYZ=1,3
      DO 4330 IXYZ=1,3
      DS(IXYZ,  2,JXYZ,JAT)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,3,JXYZ,JAT))
 4330 DS(JXYZ,JAT,IXYZ,  2)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,3,JXYZ,JAT))
      DO 4340 JXYZ=1,3
      DO 4340 IXYZ=1,3
 4340 DS(IXYZ,  2,JXYZ,  2)=-(DS(IXYZ,1,JXYZ,  2)+DS(IXYZ,3,JXYZ,  2))
      GO TO 5000
C
 4400 CONTINUE
      DO 4410 IXYZ=1,3
 4410 FD(IXYZ,3)=-(FD(IXYZ,1)+FD(IXYZ,2))
      DO 4420 I=1,IJKL
      N=IJKLG(1,I)
      DO 4420 IXYZ=1,3
 4420 FIJKL(IXYZ,3,N)=-(FIJKL(IXYZ,1,N)+FIJKL(IXYZ,2,N))
      DO 4430 JAT =1,2
      DO 4430 JXYZ=1,3
      DO 4430 IXYZ=1,3
      DS(IXYZ,  3,JXYZ,JAT)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,2,JXYZ,JAT))
 4430 DS(JXYZ,JAT,IXYZ,  3)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,2,JXYZ,JAT))
      DO 4440 JXYZ=1,3
      DO 4440 IXYZ=1,3
 4440 DS(IXYZ,  3,JXYZ,  3)=-(DS(IXYZ,1,JXYZ,  3)+DS(IXYZ,2,JXYZ,  3))
      GO TO 5000
C
 4500 CONTINUE
      DO 4510 IXYZ=1,3
 4510 FD(IXYZ,4)=-(FD(IXYZ,1)+FD(IXYZ,2)+FD(IXYZ,3))
      DO 4520 I=1,IJKL
      N=IJKLG(1,I)
      DO 4520 IXYZ=1,3
 4520 FIJKL(IXYZ,4,N)=-(FIJKL(IXYZ,1,N)+FIJKL(IXYZ,2,N)+FIJKL(IXYZ,3,N))
      DO 4530 JAT =1,3
      DO 4530 JXYZ=1,3
      DO 4530 IXYZ=1,3
      DS(IXYZ,  4,JXYZ,JAT)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,2,JXYZ,JAT)
     1                                           +DS(IXYZ,3,JXYZ,JAT))
      DS(JXYZ,JAT,IXYZ,  4)=-(DS(IXYZ,1,JXYZ,JAT)+DS(IXYZ,2,JXYZ,JAT)
     1                                           +DS(IXYZ,3,JXYZ,JAT))
 4530 CONTINUE
      DO 4540 JXYZ=1,3
      DO 4540 IXYZ=1,3
      DS(IXYZ,  4,JXYZ,  4)=-(DS(IXYZ,1,JXYZ,  4)+DS(IXYZ,2,JXYZ,  4)
     1                                           +DS(IXYZ,3,JXYZ,  4))
 4540 CONTINUE
C
 5000 CONTINUE
C
C     ----- GRADIENT CONTRIBUTION -----
C
      DO 5010 I=1,4
      IAT=KAT(I)
      DO 5010 IXYZ=1,3
 5010 DE(IXYZ,IAT)=DE(IXYZ,IAT)+FD(IXYZ,I)
C
      IF(OUT) WRITE(IW,9999)(I,IIAT,DE(I,IIAT),I,JJAT,DE(I,JJAT),
     1                       I,KKAT,DE(I,KKAT),I,LLAT,DE(I,LLAT),I=1,3)
C
C     ----- HESSIAN CONTRIBUTION -----
C
      DO 5020 J=1,4
      JC=KC(J)
      DO 5020 JXYZ=1,3
      DO 5020 I=1,4
      IC=KC(I)
      DO 5020 IXYZ=1,3
 5020 FC(IC+IXYZ,JC+JXYZ)=FC(IC+IXYZ,JC+JXYZ)+DS(IXYZ,I,JXYZ,J)
C
      IF(OUT) CALL HND_HSSOUT(FC,NC,NC,NC)
C
C     ----- DERIVATIVE FOCK MATRIX CONTRIBUTION -----
C
      CALL HND_JDDFCK(FIJKL,IJKLG,DF,LDF,NBF,DEN)
C
      RETURN
 9999 FORMAT(4(' DE(',I2,',',I3,') = ',E18.11))
      END
C
      SUBROUTINE HND_JDDFCK(FIJKL,IJKLG,DF,LDF,NBF,DEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXNBF =2048)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION FIJKL(3,4,1),IJKLG(4,1)
      DIMENSION DF(LDF,1),DEN(1)
      DIMENSION KC(4)
      DATA FOUR /4.0D+00/
      DATA PT5  /0.5D+00/
      DATA TOL  /1.0D-09/
C
      KC(1)=3*(IIAT-1)
      KC(2)=3*(JJAT-1)
      KC(3)=3*(KKAT-1)
      KC(4)=3*(LLAT-1)
C
      IJKLN=0
      DO 140 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 130 J=MINJ,JMAX
      PIJ=PNRM(J)*PI
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 120 K=MINK,KMAX
      PIJK=PNRM(K)*PIJ
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 110 L=MINL,LMAX
      PIJKL=PNRM(L)*PIJK
C
      IJKLN=IJKLN+1
      IJKLF=IJKLG(1,IJKLN)
      INUM=LOCI+I
      JNUM=LOCJ+J
      KNUM=LOCK+K
      LNUM=LOCL+L
C
      IF((KNUM.LT.INUM).OR.(KNUM.EQ.INUM.AND.LNUM.LE.JNUM)) THEN
         IJ=IA(INUM)+JNUM
         KL=IA(KNUM)+LNUM
         IK=IA(INUM)+KNUM
         IL=IA(INUM)+LNUM
         JK=IA(MAX0(JNUM,KNUM))+MIN0(JNUM,KNUM)
         JL=IA(MAX0(JNUM,LNUM))+MIN0(JNUM,LNUM)
      ELSE
         IJ=IA(KNUM)+LNUM
         KL=IA(INUM)+JNUM
         IK=IA(KNUM)+INUM
         IL=IA(KNUM)+JNUM
         JK=IA(MAX0(LNUM,INUM))+MIN0(LNUM,INUM)
         JL=IA(MAX0(LNUM,JNUM))+MIN0(LNUM,JNUM)
      ENDIF
      IF(INUM.EQ.JNUM                 ) PIJKL=PIJKL*PT5
      IF(KNUM.EQ.LNUM                 ) PIJKL=PIJKL*PT5
      IF(INUM.EQ.KNUM.AND.JNUM.EQ.LNUM) PIJKL=PIJKL*PT5
C
C
      DO 20 KAT=1,4
      IC=KC(KAT)
      DO 10 IXYZ=1,3
      VAL=FIJKL(IXYZ,KAT,IJKLF)*PIJKL
      IF( ABS(VAL).LT.TOL) GO TO 5
      DF(IJ,IC+IXYZ)=DF(IJ,IC+IXYZ)+( FOUR*VAL)*DEN(KL)
      DF(KL,IC+IXYZ)=DF(KL,IC+IXYZ)+( FOUR*VAL)*DEN(IJ)
      DF(IK,IC+IXYZ)=DF(IK,IC+IXYZ)+(-     VAL)*DEN(JL)
      DF(IL,IC+IXYZ)=DF(IL,IC+IXYZ)+(-     VAL)*DEN(JK)
      DF(JK,IC+IXYZ)=DF(JK,IC+IXYZ)+(-     VAL)*DEN(IL)
      DF(JL,IC+IXYZ)=DF(JL,IC+IXYZ)+(-     VAL)*DEN(IK)
    5 CONTINUE
   10 CONTINUE
   20 CONTINUE
C
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
      IF(OUT) THEN
      WRITE(IW,9999)
c     CALL HND_PRTRL(DF(1,1),NBF)
      CALL HND_PRTR(DF(1,1),NBF)
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -JDDFCK- , FIRST DERIVATIVE FOCK MATRIX = ')
      END
C
      SUBROUTINE HND_JDDSYM(FC,NC,DF,LDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- SYMMETRIZE THE SKELETON DERIVATIVE FOCK MATRIX -----
C     ----- SYMMETRIZE THE SECOND   DERIVATIVE      MATRIX -----
C
      PARAMETER (NDIM=35)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SYMSPD/PTR(3,144),DTR(6,288),FTR(10,480),GTR(15,720)
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSHL/MAPSHL(MXSHEL,48),MAPCTR(MXATOM,48)
      COMMON/HND_ISOPAC/INDIN(48),INDOUT(24)
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION MI(48),MJ(48)
      DIMENSION TX(NDIM,NDIM),TY(NDIM,NDIM),TZ(NDIM,NDIM),
     1          UX(NDIM,NDIM),UY(NDIM,NDIM),UZ(NDIM,NDIM)
      DIMENSION T(3,3),V(3,3),U(3,3)
      DIMENSION FC(NC,1),DF(LDF,3,1)
      DIMENSION ISO(MXSHEL,24),ISOC(MXATOM,24)
      EQUIVALENCE (ISO(1,1),MAPSHL(1,1)),(ISOC(1,1),MAPCTR(1,1))
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
C
      IF(NT.EQ.1) RETURN
C
C     ----- GET SYMMETRY DATA -----
C
      CALL HND_DAREAD(IDAF,IODA,PTR,  432, 7)
      CALL HND_DAREAD(IDAF,IODA,DTR, 1728, 8)
      CALL HND_DAREAD(IDAF,IODA,FTR, 4800, 9)
      CALL HND_DAREAD(IDAF,IODA,GTR,10800,10)
C
C     NISO =512*24/ISINGL
      NISO =MXSHEL*24/ISINGL
      CALL HND_DAREAD(IDAF,IODA,ISO ,NISO ,5)
C     NISOC=128*24/ISINGL
      NISOC=MXATOM*24/ISINGL
      CALL HND_DAREAD(IDAF,IODA,ISOC,NISOC,6)
C
      DO 30 II=1,NSHELL
      DO 10 IT=1,NTWD
   10 INDOUT(IT)=ISO(II,IT)
      CALL HND_ISOOUT(NT)
      DO 20 IT=1,NT
   20 ISO(II,IT)=INDIN(IT)
   30 CONTINUE
      DO 60 IAT=1,NAT
      DO 40 IT=1,NTWD
   40 INDOUT(IT)=ISOC(IAT,IT)
      CALL HND_ISOOUT(NT)
      DO 50 IT=1,NT
   50 ISOC(IAT,IT)=INDIN(IT)
   60 CONTINUE
C
C     ----- SYMMETRIZE DERIVATIVE FOCK MATRICES -----
C     -----   LOOP OVER UNIQUE BLOCKS (II,JJ)   -----
C
      DO 3000 II=1,NSHELL
      DO 110 ITR=1,NT
      ISH=MAPSHL(II,ITR)
      IF(ISH.GT.II) GO TO 3000
  110 MI(ITR)=ISH
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      DO 2000 JJ=1,II
      DO 220 ITR=1,NT
      JSH=MAPSHL(JJ,ITR)
      MJ(ITR)=JSH
      IF(JSH.GT.II) GO TO 2000
      ISH=MI(ITR)
      IF(ISH.GE.JSH) GO TO 210
      NSH=ISH
      ISH=JSH
      JSH=NSH
  210 IF(ISH.EQ.II.AND.JSH.GT.JJ) GO TO 2000
  220 CONTINUE
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
C     ----- LOOP OVER UNIQUE ATOMS -IAT- -----
C
      DO 1000 IAT=1,NAT
C
C     ----- APPLY PROJECTION OPERATOR -----
C
      DO 410 J=MINJ,MAXJ
      DO 410 I=MINI,MAXI
      UX(I,J)=ZERO
      UY(I,J)=ZERO
  410 UZ(I,J)=ZERO
C
      DO 500 ITR=1,NT
      KK=MI(ITR)
      LKT=KTYPE(KK)
      MINK=KMIN(KK)
      MAXK=KMAX(KK)
      LOCK=KLOC(KK)-MINK
      LL=MJ(ITR)
      LLT=KTYPE(LL)
      MINL=KMIN(LL)
      MAXL=KMAX(LL)
      LOCL=KLOC(LL)-MINL
C
      KAT=MAPCTR(IAT,ITR)
      DO 420 L=MINL,MAXL
      DO 420 K=MINK,MAXK
      LCK=LOCK+K
      LCL=LOCL+L
      LCKL=IA(MAX0(LCK,LCL))+MIN0(LCK,LCL)
      TX(K,L)=DF(LCKL,1,KAT)
      TY(K,L)=DF(LCKL,2,KAT)
      TZ(K,L)=DF(LCKL,3,KAT)
  420 CONTINUE
C
      CALL HND_RFCKR(TX,MINK,MAXK,LKT,MINL,MAXL,LLT,ITR)
      CALL HND_RFCKR(TY,MINK,MAXK,LKT,MINL,MAXL,LLT,ITR)
      CALL HND_RFCKR(TZ,MINK,MAXK,LKT,MINL,MAXL,LLT,ITR)
      DO 430 L=MINL,MAXL
      DO 430 K=MINK,MAXK
      UX(K,L)=UX(K,L)+PTR(1,3*(ITR-1)+1)*TX(K,L)
     1               +PTR(2,3*(ITR-1)+1)*TY(K,L)
     2               +PTR(3,3*(ITR-1)+1)*TZ(K,L)
      UY(K,L)=UY(K,L)+PTR(1,3*(ITR-1)+2)*TX(K,L)
     1               +PTR(2,3*(ITR-1)+2)*TY(K,L)
     2               +PTR(3,3*(ITR-1)+2)*TZ(K,L)
      UZ(K,L)=UZ(K,L)+PTR(1,3*(ITR-1)+3)*TX(K,L)
     1               +PTR(2,3*(ITR-1)+3)*TY(K,L)
     2               +PTR(3,3*(ITR-1)+3)*TZ(K,L)
  430 CONTINUE
C
  500 CONTINUE
C
      DUM=ONE/DBLE(NT)
      DO 510 J=MINJ,MAXJ
      DO 510 I=MINI,MAXI
      UX(I,J)=UX(I,J)*DUM
      UY(I,J)=UY(I,J)*DUM
      UZ(I,J)=UZ(I,J)*DUM
      LCI=LOCI+I
      LCJ=LOCJ+J
      LCIJ=IA(MAX0(LCI,LCJ))+MIN0(LCI,LCJ)
      DF(LCIJ,1,IAT)=UX(I,J)
      DF(LCIJ,2,IAT)=UY(I,J)
      DF(LCIJ,3,IAT)=UZ(I,J)
  510 CONTINUE
C
C     ----- SET EQUIVALENT BLOCKS (KK,LL) FOR EQUIVALENT ATOMS -----
C
      DO 600 ITR=1,NT
      KK=MI(ITR)
      LKT=KTYPE(KK)
      MINK=KMIN(KK)
      MAXK=KMAX(KK)
      LOCK=KLOC(KK)-MINK
      LL=MJ(ITR)
      LLT=KTYPE(LL)
      MINL=KMIN(LL)
      MAXL=KMAX(LL)
      LOCL=KLOC(LL)-MINL
C
      DO 520 J=MINJ,MAXJ
      DO 520 I=MINI,MAXI
      LCI=LOCI+I
      LCJ=LOCJ+J
      LCIJ=IA(MAX0(LCI,LCJ))+MIN0(LCI,LCJ)
      TX(I,J)=DF(LCIJ,1,IAT)
      TY(I,J)=DF(LCIJ,2,IAT)
      TZ(I,J)=DF(LCIJ,3,IAT)
  520 CONTINUE
C
      CALL HND_RFCKR(TX,MINK,MAXK,LKT,MINL,MAXL,LLT,INVT(ITR))
      CALL HND_RFCKR(TY,MINK,MAXK,LKT,MINL,MAXL,LLT,INVT(ITR))
      CALL HND_RFCKR(TZ,MINK,MAXK,LKT,MINL,MAXL,LLT,INVT(ITR))
      DO 530 L=MINL,MAXL
      DO 530 K=MINK,MAXK
      UX(K,L)= PTR(1,3*(INVT(ITR)-1)+1)*TX(K,L)
     1        +PTR(2,3*(INVT(ITR)-1)+1)*TY(K,L)
     2        +PTR(3,3*(INVT(ITR)-1)+1)*TZ(K,L)
      UY(K,L)= PTR(1,3*(INVT(ITR)-1)+2)*TX(K,L)
     1        +PTR(2,3*(INVT(ITR)-1)+2)*TY(K,L)
     2        +PTR(3,3*(INVT(ITR)-1)+2)*TZ(K,L)
      UZ(K,L)= PTR(1,3*(INVT(ITR)-1)+3)*TX(K,L)
     1        +PTR(2,3*(INVT(ITR)-1)+3)*TY(K,L)
     2        +PTR(3,3*(INVT(ITR)-1)+3)*TZ(K,L)
  530 CONTINUE
C
      KAT=MAPCTR(IAT,ITR)
      DO 540 L=MINL,MAXL
      DO 540 K=MINK,MAXK
      LCK=LOCK+K
      LCL=LOCL+L
      LCKL=IA(MAX0(LCK,LCL))+MIN0(LCK,LCL)
      DF(LCKL,1,KAT)=UX(K,L)
      DF(LCKL,2,KAT)=UY(K,L)
      DF(LCKL,3,KAT)=UZ(K,L)
C
  540 CONTINUE
  600 CONTINUE
C
 1000 CONTINUE
C
 2000 CONTINUE
 3000 CONTINUE
C
C     ----- SYMMETRIZE HESSIAN MATRIX -----
C
      DO 8000 IC=1,NAT
      DO 7000 JC=1,NAT
C
C     ----- APPLY PROJECTION OPERATOR -----
C
      DO 5100 J=1,3
      DO 5100 I=1,3
 5100 V(I,J)=ZERO
C
      DO 5700 IT=1,NT
      KC=MAPCTR(IC,IT)
      LC=MAPCTR(JC,IT)
C
      LOCL=3*(LC-1)
      LOCK=3*(KC-1)
      DO 5200 L=1,3
      DO 5200 K=1,3
 5200 T(K,L)=FC(LOCK+K,LOCL+L)
      N=3*(IT-1)
      DO 5400 L=1,3
      DO 5400 K=1,3
      SUM=ZERO
      DO 5300 M=1,3
 5300 SUM=SUM+T(K,M)*PTR(M,L+N)
 5400 U(K,L)=SUM
      DO 5600 K=1,3
      DO 5600 L=1,3
      SUM=ZERO
      DO 5500 M=1,3
 5500 SUM=SUM+PTR(M,K+N)*U(M,L)
 5600 V(K,L)=V(K,L)+SUM
 5700 CONTINUE
C
      FAC=ONE/DBLE(NT)
      LOCJ=3*(JC-1)
      LOCI=3*(IC-1)
      DO 5800 J=1,3
      DO 5800 I=1,3
 5800 V(I,J)=V(I,J)*FAC
C
C     ----- REMAP PROJECTED BLOCK ONTO EQUIVALENT BLOCKS -----
C
      DO 6600 IT=1,NT
      KC=MAPCTR(IC,IT)
      LC=MAPCTR(JC,IT)
C
      N=3*(INVT(IT)-1)
      DO 6200 L=1,3
      DO 6200 K=1,3
      SUM=ZERO
      DO 6100 M=1,3
 6100 SUM=SUM+V(K,M)*PTR(M,L+N)
 6200 U(K,L)=SUM
      DO 6400 L=1,3
      DO 6400 K=1,3
      SUM=ZERO
      DO 6300 M=1,3
 6300 SUM=SUM+PTR(M,K+N)*U(M,L)
 6400 T(K,L)=SUM
C
      LOCL=3*(LC-1)
      LOCK=3*(KC-1)
      DO 6500 L=1,3
      DO 6500 K=1,3
 6500 FC(LOCK+K,LOCL+L)=T(K,L)
 6600 CONTINUE
C
 7000 CONTINUE
 8000 CONTINUE
C
C     ----- END OF SYMMETRIZATION -----
C
      RETURN
      END
C
      SUBROUTINE HND_JDDSAV(FC,NC,DF,LDF,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
C     COMMON/HND_DSHPRT/OUT
      DIMENSION FC(NC,NC),DF(LDF,NC)
      DATA NFT17 /17/
      DATA PT5   /0.5D+00/
C 
      OUT = .FALSE.
C
C     ----- WRITE DERIVATIVE FOCK MATRICES  -----
C     ----- OFF-DIAGONAL ARE DIVIDED BY 2.0 -----
C
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9999)
C
      CALL HND_REWFIL(NFT17)
      DO 20 IC=1,NC
C
      IJ=0
      DO 10 I=1,NBF
      DO 10 J=1,I
      IJ=IJ+1
      IF(J.EQ.I) GO TO 10
      DF(IJ,IC)=DF(IJ,IC)*PT5
   10 CONTINUE
C
c*    IF(OUT) CALL HND_PRTRL(DF(1,IC),NBF)
c*    IF(OUT) CALL HND_PRTR (DF(1,IC),NBF)
C
      CALL HND_HSSWRT(NFT17,DF(1,IC),LDF)
   20 CONTINUE
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRTR(DF(1,1 ),NBF)
C
C     ----- WRITE 2E-HESSIAN -----
C
      CALL HND_HSSWRT(NFT17,FC,NC*NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9998)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_HSSOUT(FC,NC,NC,NC)
      RETURN
C9999 FORMAT(/,' DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
C9998 FORMAT(/,' HESSIAN MATRIX',/,1X,14(1H-))
 9999 FORMAT(/,' 2E-DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
 9998 FORMAT(/,' 2E-HESSIAN MATRIX',/,1X,14(1H-))
C
      END
C
C     ------ SUBROUTINES ABOVE FROM HSS.F IN HONDO -------
C

C
C     ------ SUBROUTINES TAKEN FROM DER.F IN HONDO ------
C
C     ------ SUBSTITUTE BY NWCHEM ------
c$$$      SUBROUTINE LGRDDM(D,V,E,EIJ,T,IA,NDIM)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      PARAMETER (MXATOM=500)
c$$$      CHARACTER*8 ERRMSG
c$$$      CHARACTER*8 WFNTYP,SCFTYP
c$$$      CHARACTER*8 RHF,UHF,GVB,SCF,WFNMC,WFNCI,WFNMP2,WFNMP4
c$$$      LOGICAL OUT
c$$$      COMMON/IOFILE/IR,IW
c$$$      COMMON/DAFILE/IDAF,NAV,IODA(255)
c$$$      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATOM),
c$$$     1                                         C(3,MXATOM)
c$$$      COMMON/SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
c$$$     1 NO(10),NCO,NSETO,NPAIR
c$$$      COMMON/SCFOPT/SCFTYP
c$$$      COMMON/WFNOPT/WFNTYP
c$$$      COMMON/QMTOPT/QMTTOL,NQMT
c$$$      COMMON/ MCPAR/NFZC,NCORBS,NCI
c$$$      DIMENSION EIJ(NDIM,1),T(1),IA(1)
c$$$      DIMENSION V(NDIM,1),D(1),E(1)
c$$$      DIMENSION ERRMSG(3)
c$$$      DATA ERRMSG /'PROGRAM ','STOP IN ','-LGRDDM-'/
c$$$      DATA ZERO   /0.0D+00/
c$$$      DATA PT5    /0.5D+00/
c$$$      DATA ONE    /1.0D+00/
c$$$      DATA TWO    /2.0D+00/
c$$$      DATA RHF    /'RHF     '/
c$$$      DATA UHF    /'UHF     '/
c$$$      DATA GVB    /'GVB     '/
c$$$      DATA SCF    /'SCF     '/
c$$$      DATA WFNMC  /'MC      '/
c$$$      DATA WFNCI  /'CI      '/
c$$$      DATA WFNMP2 /'MP2     '/
c$$$      DATA WFNMP4 /'MP4     '/
c$$$C
c$$$      OUT =.FALSE.
c$$$      NUM2=(NUM*(NUM+1))/2
c$$$      NUM3= NUM* NUM
c$$$C
c$$$      IF(WFNTYP.NE.SCF                  ) GO TO 100
c$$$      IF(SCFTYP.EQ.GVB                  ) GO TO 100
c$$$      IF(SCFTYP.EQ.RHF.AND.NSETO.GT.0   ) GO TO 100
c$$$C
c$$$C     ----- CLOSED-SHELL OR -UHF- CASE -----
c$$$C
c$$$      OCC=ONE
c$$$      IF(SCFTYP.NE.UHF) OCC=TWO
c$$$      CALL HND_DAREAD(IDAF,IODA,V,NUM3,15)
c$$$      CALL HND_DAREAD(IDAF,IODA,E,NUM ,17)
c$$$C
c$$$      IJ=0
c$$$      DO 20 I=1,NUM
c$$$         DO 20 J=1,I
c$$$            IJ=IJ+1
c$$$            DUM=ZERO
c$$$            DO 10 K=1,NA
c$$$   10          DUM=DUM-E(K)*V(I,K)*V(J,K)
c$$$   20       D(IJ)=DUM*OCC
c$$$C
c$$$      IF(SCFTYP.EQ.UHF) THEN
c$$$         CALL HND_DAREAD(IDAF,IODA,V,NUM3,19)
c$$$         CALL HND_DAREAD(IDAF,IODA,E,NUM ,21)
c$$$         IJ=0
c$$$         DO 40 I=1,NUM
c$$$            DO 40 J=1,I
c$$$               IJ=IJ+1
c$$$               DUM=ZERO
c$$$               DO 30 K=1,NB
c$$$   30             DUM=DUM-E(K)*V(I,K)*V(J,K)
c$$$   40          D(IJ)=D(IJ)+DUM*OCC
c$$$      ENDIF
c$$$C
c$$$      IF(OUT) THEN
c$$$         WRITE(IW,9998)
c$$$         CALL HND_PRTRL(D,NUM)
c$$$      ENDIF
c$$$      RETURN
c$$$C
c$$$  100 IF(WFNTYP.EQ.WFNCI.OR.WFNTYP.EQ.WFNMP2
c$$$     1                  .OR.WFNTYP.EQ.WFNMP4) GO TO 200
c$$$C
c$$$C     -----       -ROHF- OR -GVB- OR -MCHF- CASE           -----
c$$$C     ----- IF -NFZC- .NE.0 ( ORBITALS ARE NOT OPTIMIZED ) -----
c$$$C           THE GRADIENT CANNOT BE DONE PRESENTLY.
c$$$C
c$$$      IF(NFZC.NE.0) THEN
c$$$         WRITE(IW,9999)
c$$$         CALL HND_HNDERR(3,ERRMSG)
c$$$      ENDIF
c$$$C
c$$$      CALL HND_DAREAD(IDAF,IODA,  V,NUM3,15)
c$$$      CALL HND_DAREAD(IDAF,IODA,EIJ,NUM3,22)
c$$$      NORB=NA
c$$$      IF(WFNTYP.EQ.WFNMC) NORB=NCORBS+NCI
c$$$      IF(OUT) THEN
c$$$         WRITE(IW,9997)
c$$$         CALL HND_PRSQL(V,NORB,NUM,NDIM)
c$$$         WRITE(IW,9996)
c$$$         CALL HND_PRSQ(EIJ,NORB,NQMT,NDIM)
c$$$      ENDIF
c$$$C
c$$$C     ----- CALCULATE - C * E * C(T) - -----
c$$$C
c$$$      DO 110 IJ=1,NUM2
c$$$  110    D(IJ)=ZERO
c$$$      DO 150 I=1,NUM
c$$$         DO 130 L=1,NORB
c$$$            DUM=ZERO
c$$$            DO 120 K=1,NORB
c$$$  120          DUM=DUM+V(I,K)*EIJ(K,L)
c$$$  130       T(L)=DUM
c$$$      DO 150 J=1,NUM
c$$$         II=MAX0(I,J)
c$$$         JJ=MIN0(I,J)
c$$$         IJ=IA(II)+JJ
c$$$         DUM=ZERO
c$$$         DO 140 L=1,NORB
c$$$  140       DUM=DUM+T(L)*V(J,L)
c$$$         D(IJ)=D(IJ)-DUM*PT5
c$$$  150    CONTINUE
c$$$C
c$$$C     -----             -MP- OR -CI- WAVEFUNCTIONS              -----
c$$$C     ----- TOTAL ENERGY-WEIGHTED DENSITY IS ALREADY CALCULATED -----
c$$$C
c$$$  200 IF(WFNTYP.EQ.WFNCI ) THEN
c$$$         WRITE(IW,9999)
c$$$         CALL HND_HNDERR(3,ERRMSG)
c$$$      ENDIF
c$$$      IF(WFNTYP.EQ.WFNMP4) THEN
c$$$         WRITE(IW,9999)
c$$$         CALL HND_HNDERR(3,ERRMSG)
c$$$      ENDIF
c$$$
c$$$C     ----- -MP- WAVEFUNCTIONS -----
c$$$
c$$$      IF(WFNTYP.EQ.WFNMP2) THEN
c$$$         CALL HND_DAREAD(IDAF,IODA,D,NUM2,111)
c$$$      ENDIF
c$$$C
c$$$      IF(OUT) THEN
c$$$         WRITE(IW,9998)
c$$$         CALL HND_PRTR(D,NUM)
c$$$c        CALL HND_PRTRL(D,NUM)
c$$$      ENDIF
c$$$      RETURN
c$$$ 9999 FORMAT(' THERE ARE FROZEN CORE ORBITALS IN THE MCSCF.',/,
c$$$     1       ' OR ELSE THIS IS A -CI- OR -MP4- WAVEFUNCTION.',/,
c$$$     2       ' IN ANY CASE THE GRADIENT CAN NOT BE COMPUTED PRESENTLY')
c$$$ 9998 FORMAT(/,' -D- IN -LGRDDM- ')
c$$$ 9997 FORMAT(/,' -V- IN -LGRDDM- ')
c$$$ 9996 FORMAT(/,' -EIJ- IN -LGRDDM- ')
c$$$      END
c$$$C
      SUBROUTINE HND_JKDSET(MDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL OUT
      LOGICAL SP
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      COMMON/HND_RESTAR/NREST
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_SHLTOL/RTOL,DTOL
      COMMON/HND_SHLORD/MODSHL(512),INVSHL(512)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DSHINT/IST,JST,KST,LST
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/SYMSHL/MAPSHL(MXSHEL,48),MAPCTR(MXATOM,48)
      COMMON/ISOPAC/INDIN(48),INDOUT(24)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION MI(48)
      DIMENSION ISO(512,24)
      EQUIVALENCE (ISO(1,1),MAPSHL(1,1))
      DATA ZERO  /0.0D+00/
      DATA TEN   /1.0D+01/
      DATA RLN10 /2.30258D+00/
      DATA NTYP  /6/
C
      OUT=.FALSE.
C
      IF (OUT.AND.GA_NODEID().EQ.0) THEN
          IF(MDER.EQ.1) WRITE(IW,9999)
          IF(MDER.EQ.2) WRITE(IW,9998)
      ENDIF
C
      RTOL=RLN10*ITOL
      DTOL=TEN**(-ITOL)
C
C     ----- READ IN 1E-GRADIENT -----
C

C     CALL HND_RDGRD(DE,NREST,IST,JST,KST,LST)
      NREST = 4
      IST =1
      JST =1
      KST =1
      LST =1
C
C     ----- GET SYMMETRY MAPPING OF SHELLS -----
C
C     -> ISINGL = 2
      CALL HND_HONDO
C
      NISO=512*24/ISINGL
      CALL HND_DAREAD(IDAF,IODA,ISO,NISO,5)
      DO 130 II=1,NSHELL
      DO 110 IT=1,NTWD
  110 INDOUT(IT)=ISO(II,IT)
      CALL HND_ISOOUT(NT)
      DO 120 IT=1,NT
  120 ISO(II,IT)=INDIN(IT)
  130 CONTINUE
C
C     ----- DEFINE NEW ORDER OF SHELLS ACCORDING TO TYPE -----
C
      DO 200 II=1,NSHELL
      MODSHL(II)=II
  200 INVSHL(II)=II
C
      IINEW=0
      DO 270 IITYP=1,NTYP
      DO 260 II=1,NSHELL
      LIT=KTYPE(II)
      SP =KTYPE(II).EQ.2.AND.KMIN(II).EQ.1
      IF(SP) LIT=NTYP
      IF(LIT.NE.IITYP) GO TO 260
      DO 210 IT=1,NT
      ID=MAPSHL(II,IT)
      IF( ID.GT.II   ) GO TO 260
  210 MI(IT)=ID
C
      IF(NT.EQ.1) GO TO 240
      DO 230 IT=2,NT
      MAX=IT-1
      DO 220 JT=1,MAX
      IF(MI(JT).NE.MI(IT)) GO TO 220
      MI(IT)=0
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
C
C     ----- LOOP -250- MUST BE ACTIVATED ONLY WHEN -----
C           LOOP -250- IS ACTIVATED IN -JKSET- .
C     ----- PRESENTLY IT IS NOT ACTIVATED          -----
C
      DO 250 IT=1,NT
      IIOLD=MI(NT+1-IT)
      IF(IIOLD.EQ.0) GO TO 250
      IINEW=IINEW+1
C.... MODSHL(IINEW)=IIOLD
C.... INVSHL(IIOLD)=IINEW
  250 CONTINUE
C
  260 CONTINUE
  270 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,22(1H-),/,10X,'GRADIENT OF THE ENERGY',
     1       /,10X,22(1H-))
 9998 FORMAT(/,10X,25(1H-),/,10X,'GRAD + HESS OF THE ENERGY',
     1       /,10X,25(1H-))
      END
C
      SUBROUTINE HND_JKDMEM(MDER,IADDR,LFIX0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL OUT
      LOGICAL SP
      LOGICAL NORM
      LOGICAL VECTOR
      PARAMETER (LENVEC=255)
      CHARACTER*8 ERRMSG
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MACHNV/VECTOR
      COMMON/HND_MEMORY/MAXCOR,MAXSCM
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLBAS/MAXTYP,MAXNUM
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERMEM/IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,
     1 IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/HND_INTTYP/NPOPLE,NHONDO
      COMMON/HND_INTXYZ/MAXXYZ,NUMXYZ
      COMMON/HND_INTVEC/MAXVEC,MINVEC
      DIMENSION IJKLN(5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-JKDMEM-'/
      DATA IJKLN /   1,  4, 10, 20, 35/
      DATA ONE   /1.00000000000000D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
C
      MAXVEC=LENVEC/3
      IF(     VECTOR) MINVEC=24
      IF(.NOT.VECTOR) MINVEC=MAXVEC*3+1
C
C     ----- THIS PROGRAM FOR DERIVATIVES -----
C
      NDER=MDER
      IF(NDER.NE.1.AND.NDER.NE.2) NDER=0
      NDER0=0
      NDER1=0
      NDER2=0
      IF(OUT              ) NDER0=1
      IF(OUT.AND.NDER.GE.1) NDER1=12
      IF(OUT.AND.NDER.EQ.2) NDER2=78
      IF(        NDER.EQ.2) NDER1=12
C
C     ----- CHECK MAXIMUM ANGULAR MOMENTUM -----
C
      SP=.FALSE.
      MAXTYP=0
      DO 10 I=1,NSHELL
      SP=SP.OR.(KTYPE(I).EQ.2.AND.KMIN(I).EQ.1)
      IF(KTYPE(I).GT.MAXTYP) MAXTYP=KTYPE(I)
   10 CONTINUE
      IF(MAXTYP.LT.5              ) GO TO 15
      IF(MAXTYP.EQ.5.AND.NDER.LE.1) GO TO 15
      WRITE(IW,9999)
      CALL HND_HNDERR(3,ERRMSG)
   15 CONTINUE
      MINXYZ=(4*MAXTYP -2 +NDER)/2
C
C     ----- GET NUMBER OF PRIMITIVE CHARGE DISTRIBUTIONS -----
C
      NIJG=0
      DO 20 II=1,NSHELL
      DO 20 JJ=1,II
      NIJG=NIJG+KNG(II)*KNG(JJ)
   20 CONTINUE
C
C     -----  AT THIS POINT IT IS GOOD TO REMEMBER THAT    -----
C            -MAXTYP- = HIGHEST SHELL ANGULAR MOMENTUM
C            -MAXFUN- = NUMBER OF FUNCTIONS WITH ANGULAR
C                       MOMENTUM LESS OR EQUAL TO -MAXTYP-
C            -MAXNUM- = NUMBER OF FUNCTIONS WITH ANGULAR
C                       MOMENTUM         EQUAL TO -MAXTYP-
C            -MAXXYZ- = MAXIMUM NUMBER OF PRIMITIVE INTEGRALS
C                       THAT CAN BE HANDLED IN ONE -VECTOR-
C            -NUMXYZ- = ACTUAL MAXIMUM LENGTH OF ONE -VECTOR-
C            -MAXXYZ- = IT IS NUMXYZ/3 . SINCE THE X, Y, AND Z
C                       COMPONENTS ARE TREATED AS A SINGLE VECTOR,
C                       -MAXXYZ- IS THE NUMBER OF (PRIMITIVE-ROOTS)
C                       COMBINATIONS WHICH CAN BE TREATED IN ONE
C                       VECTOR. FOR -SSSS- INTEGRALS WHICH REQUIRE
C                       ONE RYS ROOT, MAXXYZ HAPPENS TO COINCIDE WITH
C                       THE NUMBER OF PRIMITIVE INTEGRALS TREATED IN
C                       ONE VECTOR. FOR -DDDD- INTEGRALS WHICH
C                       REQUIRE FIVE RYS ROOTS, THE NUMBER OF PRIMITIVE
C                       INTEGRALS TREATED IN ONE VECTOR IS -MAXXYZ-/5 .
C
C
C     ----- SET NORMALIZATION CONSTANTS -----
C
      MAXFUN=IJKLN(MAXTYP)
      DO 100 I=1,MAXFUN
  100 PNRM(I)=ONE
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      IF(.NOT.NORM) GO TO 180
C
      SQRT53=SQRT5/SQRT3
      DO 170 I=1,MAXFUN
      GO TO(110,110,160,160,110,160,160,120,160,160,
     1      110,160,160,130,160,160,160,160,160,120,
     2      110,160,160,140,160,160,160,160,160,150,
     3      160,160,120,160,160),I
  110 FI=ONE
      GO TO 160
  120 FI=SQRT3*FI
      GO TO 160
  130 FI=SQRT5*FI
      GO TO 160
  140 FI=SQRT7*FI
      GO TO 160
  150 FI=SQRT53*FI
  160 CONTINUE
      PNRM(I)=FI
  170 CONTINUE
C
  180 CONTINUE
C
      IDUM=MAXTYP-1
      MAXNUM=((IDUM+1)*(IDUM+2))/2
      IF(MAXNUM.EQ.3.AND.SP) MAXNUM=4
      NGIJKL=(MAXNUM**4)
C
C     ----- FOR DERIVATIVES -----
C
      MODTYP=MAXTYP+NDER
C
C     ----- CALCULATE VECTOR LENGTH AND SET CORE POINTERS -----
C
      LVAR=0
      LFIX=LFIX0
      LFIX=LFIX  +(NSHELL*(NSHELL+1))/2
C
C     ----- -SP- FUNCTIONS FOR FIRST DERIVATIVES ARE SPECIAL -----
C
      IWFN =IADDR
      ICLB =IWFN + LFIX0
      ILAST=ICLB + (NSHELL*(NSHELL+1))/2
      INEED=ILAST- 1
      IADDR=ILAST
      IF(NDER.EQ.1.AND.NHONDO.EQ.0) GO TO 300
C
      LFIX=LFIX+( (NSHELL*(NSHELL+1))/2 )*2
      LFIX=LFIX+NGIJKL*(NDER0+NDER1+NDER2)
      LFIX=LFIX+NGIJKL* 4
      LFIX=LFIX+NGIJKL
      LFIX=LFIX+NIJG*15
      LVAR=     ( MODTYP**2       * MODTYP**2       )*3
      LVAR=LVAR+( MODTYP**2       *(MODTYP+MODTYP-1))*3
      LVAR=LVAR+((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*3
      LVAR=LVAR+( MODTYP**2                         )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)                  )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)* 3               )*3
      LVAR=LVAR+(  3                                )*3
      LVAR=LVAR+(  9                                )
      LVAR=LVAR+(  4                                )
      LVAR=LVAR+(  5                                )
      LVAR=LVAR+( 18                                )
      LVAR=LVAR+(  2                                )
      LVAR=LVAR+(  4                                )*3
      LVAR=LVAR+( MODTYP**2       * MODTYP**2       )*3*14
C
      MAXCOR = MAXVEC*LVAR+LFIX+(IADDR-1)
C
      MAXXYZ=(MAXCOR-LFIX-(IADDR-1))/LVAR
      IF(MAXXYZ.GE.MINXYZ) GO TO 200
      WRITE(IW,9998) LVAR,LFIX,MAXXYZ,MINXYZ,MAXCOR
      CALL HND_HNDERR(3,ERRMSG)
  200 CONTINUE
      IF(MAXXYZ.GT.MAXVEC) MAXXYZ=MAXVEC
      NUMXYZ=3*MAXXYZ
C
C     X(IWFN  ) = WAVEFUNCTION DATA
C     X(ICLB  ) = COULOMB INTEGRAL THRESHOLD
C     X(INIJG ) = CHARGE DISTRIBUTION POINTERS
C     X(IGINT ) = ELECTRON REPULSION INTEGRALS
C     X(IFINT ) = FIRST DERIVATIVE INTEGRALS
C     X(ISINT ) = SECOND DERIVATIVE INTEGRALS
C     X(IIJKLG) = INDICES
C     X(IDAB  ) = DENSITY ARRAY -DAB-
C     X(ICHRG ) = CHARGE DISTRIBUTION PARAMETERS
C     X(IXY   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IXZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IYZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IX    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IY    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IZ    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(ISJ   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(ISK   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(ISL   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(IGIJKL) = ( 2-D , 4 CENTERS ) INTEGRALS
C     X(IGNKL ) = ( 2-D , 3 CENTERS ) INTEGRALS
C     X(IGNM  ) = ( 2-D , 2 CENTERS ) INTEGRALS
C     X(IDIJ  ) = CONTRACTION DENSITY FOR -IJ- CHARGE DISTRIBUTION
C     X(IDKL  ) = CONTRACTION DENSITY FOR -KL- CHARGE DISTRIBUTION
C     X(IB00  ) = -B00-
C     X(IB01  ) = -B01-
C     X(IB10  ) = -B10-
C     X(IC00  ) = -C00-
C     X(ID00  ) = -D00-
C     X(IF00  ) = -F00-
C     X(IDIJSI) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- II SHELL
C     X(IDIJSJ) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- JJ SHELL
C     X(IDKLSK) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- KK SHELL
C     X(IDKLSL) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- LL SHELL
C     X(IABV  ) = -AB- VECTOR FOR PRIMITIVE INTEGRALS
C     X(ICV   ) = -CV- VECTOR FOR PRIMITIVE INTEGRALS
C     X(IRW   ) = -RW- VECTOR FOR RYS ROOTS AND WEIGHTS
C     X(IAAI  ) = EXPONENT FOR DERIVATIVE OF II SHELL
C     X(IAAJ  ) = EXPONENT FOR DERIVATIVE OF JJ SHELL
C     X(IBBK  ) = EXPONENT FOR DERIVATIVE OF KK SHELL
C     X(IBBL  ) = EXPONENT FOR DERIVATIVE OF LL SHELL
C     X(IFI   ) = FIRST DERIVATIVE WRT. II OF ( 2-D , 4 CENTERS ) INT.
C     X(IFJ   ) = FIRST DERIVATIVE WRT. JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(IFK   ) = FIRST DERIVATIVE WRT. KK OF ( 2-D , 4 CENTERS ) INT.
C     X(IFL   ) = FIRST DERIVATIVE WRT. LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISII  ) = SECOND DER. WRT II AND II OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJJ  ) = SECOND DER. WRT JJ AND JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(ISKK  ) = SECOND DER. WRT KK AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISLL  ) = SECOND DER. WRT LL AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIJ  ) = SECOND DER. WRT II AND JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIK  ) = SECOND DER. WRT II AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIL  ) = SECOND DER. WRT II AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJK  ) = SECOND DER. WRT JJ AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJL  ) = SECOND DER. WRT JJ AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISKL  ) = SECOND DER. WRT KK AND LL OF ( 2-D , 4 CENTERS ) INT.
C
      IWFN  =IWFN
      ICLB  =ICLB
      INIJG =ICLB  + (NSHELL*(NSHELL+1))/2
      IGINT =INIJG +((NSHELL*(NSHELL+1))/2)*2
      IFINT =IGINT +  NGIJKL*NDER0
      ISINT =IFINT +  NGIJKL*NDER1
      IIJKLG=ISINT +  NGIJKL*NDER2
      IDAB  =IIJKLG+  NGIJKL*4
      ICHRG =IDAB  +  NGIJKL
      IXY   =ICHRG +  NIJG  *15
      IXZ   =IXY   +(  1                                )*MAXXYZ
      IYZ   =IXZ   +(  1                                )*MAXXYZ
      IX    =IYZ   +(  1                                )*MAXXYZ
      IY    =IX    +(  1                                )*MAXXYZ
      IZ    =IY    +(  1                                )*MAXXYZ
      ISJ   =IZ    +(  1                                )*MAXXYZ
      ISK   =ISJ   +(  1                                )*MAXXYZ
      ISL   =ISK   +(  1                                )*MAXXYZ
      IGIJKL=ISL   +(  1                                )*MAXXYZ
      IGNKL =IGIJKL+( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IGNM  =IGNKL +( MODTYP**2       *(MODTYP+MODTYP-1))*MAXXYZ*3
      IDIJ  =IGNM  +((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*MAXXYZ*3
      IDKL  =IDIJ  +( MODTYP**2                         )*MAXXYZ*3
      IB00  =IDKL  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IB01  =IB00  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IB10  =IB01  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IC00  =IB10  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      ID00  =IC00  +(  1                                )*MAXXYZ*3
      IF00  =ID00  +(  1                                )*MAXXYZ*3
      IDIJSI=IF00  +(  1                                )*MAXXYZ*3
      IDIJSJ=IDIJSI+(  1                                )*MAXXYZ
      IDKLSK=IDIJSJ+(  1                                )*MAXXYZ
      IDKLSL=IDKLSK+(  1                                )*MAXXYZ
      IABV  =IDKLSL+(  1                                )*MAXXYZ
      ICV   =IABV  +(  5                                )*MAXXYZ
      IRW   =ICV   +( 18                                )*MAXXYZ
      IAAI  =IRW   +(  2                                )*MAXXYZ
      IAAJ  =IAAI  +(  1                                )*MAXXYZ*3
      IBBK  =IAAJ  +(  1                                )*MAXXYZ*3
      IBBL  =IBBK  +(  1                                )*MAXXYZ*3
      IFI   =IBBL  +(  1                                )*MAXXYZ*3
      IFJ   =IFI   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IFK   =IFJ   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IFL   =IFK   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISII  =IFL   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJJ  =ISII  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISKK  =ISJJ  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISLL  =ISKK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIJ  =ISLL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIK  =ISIJ  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIL  =ISIK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJK  =ISIL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJL  =ISJK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISKL  =ISJL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ILAST =ISKL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      INEED =ILAST -1
      IADDR =ILAST

C     WRITE(IW,9997) NUMXYZ,INEED,MAXCOR

      IF(OUT) WRITE(IW,9996) IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
C
      RETURN
  300 CONTINUE
      IF(OUT) WRITE(IW,9995) IWFN,ICLB,INEED,MAXCOR
C
      RETURN
C
 9999 FORMAT(' FOR -G- FUNCTIONS , -NDER- MAY NOT BE GREATER THAN 1 ',
     1     /,' AS THE PROGRAM IS DIMENSIONED TO -9- RYS ROOTS. STOP ',
     2     /,' FOR HIGHER FUNCTIONS NO DERIVATIVES ARE AVAILABLE    ')
 9998 FORMAT(/,' NOT ENOUGH CORE FOR -JKDMEM-',/,
     1         ' LVAR,LFIX,MAXXYZ,MINXYZ,MAXCOR = ',5I10)
C9997 FORMAT(' - DEW  - , MAXVEC = ',I3,' INEED = ',I10,
C    1                  ' MAXCOR = ',I10)
 9996 FORMAT(
     1 ' IWFN  ',I8,' ICLB  ',I8,' INIJG ',I8,' IGINT ',I8,' IFINT ',I8,
     1 ' ISINT ',I8,' IIJKLG',I8,/,
     1 ' IDAB  ',I8,' ICHRG ',I8,' IXY   ',I8,' IXZ   ',I8,' IYZ   ',I8,
     1 ' IX    ',I8,' IY    ',I8,' IZ    ',I8,/,
     1 ' ISJ   ',I8,' ISK   ',I8,' ISL   ',I8,
     1 ' IGIJKL',I8,' IGNKL ',I8,' IGNM  ',I8,/,
     2 ' IDIJ  ',I8,' IDKL  ',I8,' IB00  ',I8,' IB01  ',I8,' IB10  ',I8,
     2 ' IC00  ',I8,' ID00  ',I8,' IF00  ',I8,/,
     3 ' IAAI  ',I8,' IAAJ  ',I8,' IBBK  ',I8,' IBBL  ',I8,
     3 ' IFI   ',I8,' IFJ   ',I8,' IFK   ',I8,' IFL   ',I8,/,
     4 ' ISII  ',I8,' ISJJ  ',I8,' ISKK  ',I8,' ISLL  ',I8,
     4 ' ISIJ  ',I8,' ISIK  ',I8,' ISIL  ',I8,' ISJK  ',I8,/
     4 ' ISJL  ',I8,' ISKL  ',I8,/,
     5 ' IDIJSI',I8,' IDIJSJ',I8,' IDKLSK',I8,' IDKLSL',I8,/,
     5 ' IABV  ',I8,' ICV   ',I8,' IRW   ',I8)
 9995 FORMAT(' SPECIAL -SP- ROUTINES USED.',
     1       ' IWFN,ICLB,INEED,MAXCOR = ',4I8)
      END
C
      SUBROUTINE HND_JKDHF(DA,DB,DC,V,NCONF,IA,NUM,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL NOCORE,NOOPEN,NOPAIR,SUHF,ROHF,RGVB,ROGVB,MC,CI,MCCI
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SCFPAR/NOCORE,NOOPEN,NOPAIR,SUHF,ROHF,RGVB,ROGVB,
     1 MC,CI,MCCI,MP2,MP4
C    1 MC,CI,MCCI
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      DIMENSION DA(1),DB(1),DC(1),IA(1),V(NDIM,1),NCONF(1)
C
      DATA ZERO /0.0D+00/
C
      NUM2=(NUM*(NUM+1))/2
      NUM3= NUM* NUM
      IF(ROGVB) GO TO 100
C
C     ----- CLOSED-SHELL OR -UHF- CASE -----
C
C
C     ----- READ IN DENSITY MATRICES (ALPHA+BETA) -----
C
      CALL HND_DAREAD(IDAF,IODA,DA,NUM2,16)
      IF(.NOT.SUHF) GO TO 30
      CALL HND_DAREAD(IDAF,IODA,DB,NUM2,20)
      DO 20 I=1,NUM2
      DUMA=DA(I)
      DUMB=DB(I)
      DA(I)=DUMA+DUMB
      DB(I)=DUMA-DUMB
   20 CONTINUE
   30 CONTINUE
      GO TO 200
C
C     ----- -ROHF- OR -GVB- CASE -----
C
  100 CONTINUE
C
C     ----- SET UP MO'S TO FOCK OPERATORS MAPPING -----
C
      IF(NCO.EQ.0) GO TO 120
      DO 110 I=1,NCO
  110 NCONF(I)=1
  120 CONTINUE
      NOP=0
      IF(NSETO.EQ.0) GO TO 150
      DO 140 ISET=1,NSETO
      IOP=NO(ISET)
      DO 130 I=1,IOP
  130 NCONF(NCO+NOP+I)=1+ISET
  140 NOP=NOP+IOP
  150 CONTINUE
      IF(NPAIR.EQ.0) GO TO 170
      NGEM=NPAIR+NPAIR
      DO 160 IGEM=1,NGEM
  160 NCONF(NCO+NOP+IGEM)=1+NSETO+IGEM
  170 CONTINUE
C
      NORB=NCO+NOP+(NPAIR+NPAIR)
      NOCORE=NCO.EQ.0
      NOOPEN=NSETO.EQ.0
      NOPAIR=NPAIR.EQ.0

C
C     ----- GET -ROHF- OR -GVB- MO'S. FORM CORE DENSITY -----
C
      CALL HND_DAREAD(IDAF,IODA,V,NUM3,15)
      IF(NOCORE) GO TO 200
      DO 190 I=1,NUM
      DO 190 J=1,I
      DUM=ZERO
      DO 180 K=1,NCO
  180 DUM=DUM+V(I,K)*V(J,K)
      IJ=IA(I)+J
  190 DC(IJ)=DUM
C
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JKDCLU(UHFTYP,IA,DA,DB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL UHFTYP
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLNRM/PNRM(35)
      DIMENSION IA(1),DA(1),DB(1)
C
      RETURN
      END
C
      SUBROUTINE HND_JKDROH(NORB,IA,DC,V,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLNRM/PNRM(35)
      DIMENSION IA(1),DC(1),V(NDIM,1)
C
      RETURN
      END
      SUBROUTINE HND_JKDGVB(NORB,IA,DC,V,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLNRM/PNRM(35)
      DIMENSION IA(1),DC(1),V(NDIM,1)
C
      RETURN
      END
C
      SUBROUTINE HND_JKDMC(DC,DV,V,GIJKL,GIJ,IX,XX,IA,NUM,NCORB,NCI,
     1NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      INTEGER ZSHIFTR
      COMMON/HND_INTFIL/NINTMX
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      DIMENSION GIJKL(NCI,NCI,NCI,NCI),GIJ(1),DC(1),DV(1),V(NUM,1),IA(1)
      DIMENSION XX(1),IX(1)
      DIMENSION ERRMSG(3)
#include "bitops_decls.fh"
      DATA ERRMSG /'PROGRAM ','STOP IN ','- JKDMC-'/
      DATA ZERO   /0.0D+00/
      DATA TWO    /2.0D+00/
      DATA GTOL   /1.0D-09/
      DATA NFTG   /15/
c
#include "bitops_funcs.fh"
C
      ZSHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      MASK1(IARG)=2**IARG-1
      IWORD1(IARG)=LAND(ZSHIFTR(IARG,24),MASK1(8))
      IWORD2(IARG)=LAND(ZSHIFTR(IARG,16),MASK1(8))
      IWORD3(IARG)=LAND(ZSHIFTR(IARG, 8),MASK1(8))
      IWORD4(IARG)=LAND(ZSHIFTR(IARG, 0),MASK1(8))
C
C     ----- LECTURE DES GIJ -----
C
      NGIJ=(NCI*(NCI+1))/2
      CALL HND_DAREAD(IDAF,IODA,GIJ,NGIJ,23)
C
      DO 10 I=1,NCI
      DO 10 J=1,NCI
      DO 10 K=1,NCI
      DO 10 L=1,NCI
   10 GIJKL(L,K,J,I)=ZERO
C
C
C     -----  CONSTRUCT -DM2- IN CORE  FROM NFTG ------
C
      CALL HND_REWFIL(NFTG)
  100 CALL HND_PREAD(NFTG,XX,IX,NXX,NINTMX)
      IF(NXX.EQ.0) GO TO 120
      MGIJKL=IABS(NXX)
      IF(MGIJKL.GT.NINTMX) CALL HND_HNDERR(3,ERRMSG)
      DO 110 MG=1,MGIJKL
      GGIJKL=XX(MG)
      IF( ABS(GGIJKL).LT.GTOL) GO TO 110
      IJKLG=IX(MG)
      IG=IWORD1(IJKLG)
      JG=IWORD2(IJKLG)
      KG=IWORD3(IJKLG)
      LG=IWORD4(IJKLG)
      GIJKL(LG,KG,IG,JG) = GGIJKL
      GIJKL(LG,KG,JG,IG) = GGIJKL
      GIJKL(KG,LG,IG,JG) = GGIJKL
      GIJKL(KG,LG,JG,IG) = GGIJKL
      GIJKL(IG,JG,KG,LG) = GGIJKL
      GIJKL(IG,JG,LG,KG) = GGIJKL
      GIJKL(JG,IG,KG,LG) = GGIJKL
      GIJKL(JG,IG,LG,KG) = GGIJKL
  110 CONTINUE
      IF(NXX.GT.0) GO TO 100
  120 CONTINUE
C
C     -----  FORM CORE DENSITY -----
C
      NUM3= NUM* NUM
      CALL HND_DAREAD(IDAF,IODA,V,NUM3,15)
      IF(NCORB.EQ.0) GO TO 300
      DO 220 I=1,NUM
      DO 220 J=1,I
      DUM=ZERO
      DO 210 K=1,NCORB
  210 DUM=DUM+V(I,K)*V(J,K)
      IJ=IA(I)+J
  220 DC(IJ)=DUM+DUM
C
  300 CONTINUE
C
C     -----  FORM VALENCE DENSITY -----
C
      IF (NCI.EQ.0) GO TO 400
C
      DO 320 I=1,NUM
      DO 320 J=1,I
      DUM=ZERO
      DO 310 K=1,NCI
      DO 310 L=1,NCI
      KL=IA(MAX0(K,L))+MIN0(K,L)
      DUM=DUM+V(I,K+NCORB)*V(J,L+NCORB)*GIJ(KL)
  310 CONTINUE
      IJ=IA(I)+J
  320 DV(IJ)=DUM
C
  400 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKDCI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_CIDAB/NDMAX,ND
      DATA NFTI /16/
C
      ND=0
      CALL HND_REWFIL(NFTI)
      RETURN
      END
C
      SUBROUTINE HND_JKDSYM(DE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- SYMMETRYZE GRADIENT VECTOR -----
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSHL/MAPSHL(MXSHEL,48),MAPCTR(MXATOM,48)
      COMMON/HND_SYMSPD/PTR(3,144),DTR(6,288),FTR(10,480),GTR(15,720)
C     COMMON/HND_SYMSPD/PTR(3,144)
      COMMON/HND_ISOPAC/INDIN(48),INDOUT(24)
      DIMENSION DE(3,1)
      DIMENSION ISOC(MXATOM,24)
      EQUIVALENCE (ISOC(1,1),MAPCTR(1,1))
      DATA ZERO,ONE /0.0D+00,1.0D+00/
C
      IF(NT.EQ.1) RETURN
C
C     ----- READ IN TRANFORMATION MATRICES OF COORDINATES. -----
C
      CALL HND_DAREAD(IDAF,IODA,PTR,432,7)
C
C     ----- READ IN TRANSFORMATION TABLE OF ATOMS -----
C           VERSUS SYMMETRY OPERATIONS
C
      NISOC=128*24/ISINGL
      CALL HND_DAREAD(IDAF,IODA,ISOC,NISOC,6)
C
      DO 30 IAT=1,NAT
      DO 10 IT=1,NTWD
   10 INDOUT(IT)=ISOC(IAT,IT)
      CALL HND_ISOOUT(NT)
      DO 20 IT=1,NT
   20 ISOC(IAT,IT)=INDIN(IT)
   30 CONTINUE
C
C     ----- LOOP OVER UNIQUE ATOMS -----
C
      DO 200 IC=1,NAT
      DO 110 IT=1,NT
      IF(MAPCTR(IC,IT).GT.IC) GO TO 200
  110 CONTINUE
C
C     ----- APPLY PROJECTION OPERATOR -----
C
      DEDX=ZERO
      DEDY=ZERO
      DEDZ=ZERO
      DO 120 IT=1,NT
      JC   =MAPCTR(IC,IT)
      DEDXP=DE(1,JC)
      DEDYP=DE(2,JC)
      DEDZP=DE(3,JC)
      N=3*(IT-1)
      DEDX=DEDX+DEDXP*PTR(1,N+1)+DEDYP*PTR(2,N+1)+DEDZP*PTR(3,N+1)
      DEDY=DEDY+DEDXP*PTR(1,N+2)+DEDYP*PTR(2,N+2)+DEDZP*PTR(3,N+2)
      DEDZ=DEDZ+DEDXP*PTR(1,N+3)+DEDYP*PTR(2,N+3)+DEDZP*PTR(3,N+3)
  120 CONTINUE
      DUM=ONE/DBLE(NT)
      DEDX=DEDX*DUM
      DEDY=DEDY*DUM
      DEDZ=DEDZ*DUM
C
C     ----- REMAP PROJECTED BLOCK ONTO EQUIVALENT ATOMS -----
C
      DO 130 IT=1,NT
      JC=MAPCTR(IC,IT)
      N =3*(INVT(IT)-1)
      DE(1,JC)=DEDX*PTR(1,N+1)+DEDY*PTR(2,N+1)+DEDZ*PTR(3,N+1)
      DE(2,JC)=DEDX*PTR(1,N+2)+DEDY*PTR(2,N+2)+DEDZ*PTR(3,N+2)
      DE(3,JC)=DEDX*PTR(1,N+3)+DEDY*PTR(2,N+3)+DEDZ*PTR(3,N+3)
  130 CONTINUE
C
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKDSAV(DE,MDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_RESTAR/NREST
      COMMON/HND_DSHINT/IST,JST,KST,LST
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_FUNCT/EHF,EG(3,MXATOM)
      DIMENSION DE(3,1)
      LOGICAL OUT
      DIMENSION DNAM(3)
      DATA DNAM /4HE'X ,4HE'Y ,4HE'Z /
C
      OUT = .FALSE.
C
C     ----- TRANSFER GRADIENT INTO -EG- -----
C
      DO 10 IAT=1,NAT
      DO 10 IXYZ=1,3
   10 EG(IXYZ,IAT)=DE(IXYZ,IAT)
C
C     ----- SAVE FINAL-GRADIENT -----
C
      NREST=0
      IST=1
      JST=1
      KST=1
      LST=1
      CALL HND_WRTGRD(EG,NREST,IST,JST,KST,LST)
C
C     ----- PRINT FINAL GRADIENT -----
C
      MAX=0
   20 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      IF (OUT.AND.GA_NODEID().EQ.0) THEN
         WRITE(IW,9994)
         WRITE(IW,9998)
         WRITE(IW,9997) (IAT,IAT=MIN,MAX)
         WRITE(IW,9998)
         DO 30 IXYZ=1,3
   30    WRITE(IW,9996) DNAM(IXYZ),(EG(IXYZ,IAT),IAT=MIN,MAX)
      ENDIF
      IF(MAX.LT.NAT) GO TO 20
C
      IF (OUT.AND.GA_NODEID().EQ.0) THEN
         IF(MDER.EQ.1) WRITE(IW,9999)
         IF(MDER.EQ.2) WRITE(IW,9995)
      ENDIF
C
      RETURN
 9999 FORMAT(/,
     1 ' ...... END OF CALCULATION OF THE ENERGY GRADIENT ...... ',/)
 9998 FORMAT(/)
 9997 FORMAT(5X,'ATOM',8(9X,I2,4X))
 9996 FORMAT(6X,A3,8F15.7)
 9995 FORMAT(/,
     1 ' ...... END OF CALCULATION OF THE ENERGY GRAD + HESS ...... ',/)
 9994 FORMAT(/,' 1E+2E-GRADIENT MATRIX ',/,1X,24(1H-))
      END
C

      SUBROUTINE HND_OEDHND(NIJG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
      PARAMETER (MXGSH =30,MXGSH2=MXGSH*MXGSH)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      COMMON/HND_DERMEM/IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLORD/MODSHL(MXSHEL),INVSHL(MXSHEL)
      COMMON/HND_IJDCHR/ADCHRG(13,MXGSH2),DASI(MXGSH2),DASJ(MXGSH2),NIJ
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION NIJG(2,*)
C
C     ----- ONE-ELECTRON CHARGE DISTRIBUTION AND DERIVATIVES -----
C
      NIJ0=0
      DO 9000 II=1,NSHELL
      IIMOD=MODSHL(II)
      DO 8000 JJ=1,II
      JJMOD=MODSHL(JJ)
      ISH=IIMOD
      JSH=JJMOD
      CALL HND_OEDSHL(ISH,JSH)
      CALL HND_OEDCHR
      IIJJ=IA(MAX0(IIMOD,JJMOD))+MIN0(IIMOD,JJMOD)
      NIJG(1,IIJJ)=NIJ0
      NIJG(2,IIJJ)=NIJ
      IF(NIJ.EQ.0) GO TO 7000
      CALL HND_OEDWR(X(ICHRG),NIJ0,ADCHRG,DASI,DASJ,NIJ)
 7000 CONTINUE
      NIJ0=NIJ0+NIJ
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKDHND(GCLB,NIJG,
     1 DDA,DDB,DDC,V,NCONF,NORB,
     2 DCMC,DVMC,VMC,DIJKL,DMJKL,DMNKL,DMNPL,DMNPQ,LDM2,
     3 DM2,IDM2,NBUF,
     4 FC,NC,DF,LDF,NBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "global.fh"
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      PARAMETER (MXGSH =30,MXGSH2=MXGSH*MXGSH)
      LOGICAL LCAP
      LOGICAL SOME,OUT
      LOGICAL SP
      LOGICAL SPSKIP
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL NOCORE,NOOPEN,NOPAIR
      LOGICAL SUHF,ROHF,RGVB,ROGVB,MC,CI,MCCI,MP2,MP4
      LOGICAL IEQJ0,KEQL0,SAME0,FIRST
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_JKSTAT/LVEC,NVEC
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERMEM/IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/HND_INTXYZ/MAXXYZ,NUMXYZ
      COMMON/HND_INTHND/SPSKIP
      COMMON/HND_INTTYP/NPOPLE,NHONDO
      COMMON/HND_SCFPAR/NOCORE,NOOPEN,NOPAIR,SUHF,ROHF,RGVB,ROGVB,
     1 MC,CI,
     1              MCCI,MP2,MP4
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLORD/MODSHL(MXSHEL),INVSHL(MXSHEL)
      COMMON/HND_SHLFAC/QQ4
      COMMON/HND_SHLBAS/MAXTYP,MAXNUM
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLOLD/NUMI0,NUMJ0,NUMK0,NUML0,
     1              IEQJ0,KEQL0,SAME0,FIRST
      COMMON/HND_IJDCHR/ADCHRG(13,MXGSH2),DASI(MXGSH2),DASJ(MXGSH2),NIJ
      COMMON/HND_KLDCHR/BDCHRG(13,MXGSH2),DBSK(MXGSH2),DBSL(MXGSH2),NKL
      COMMON/HND_DSHINT/IST,JST,KST,LST
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DSHDAB/DABTOL,DABMAX
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSHL/MAPSHL(MXSHEL,48),MAPCTR(MXATOM,48)
      COMMON/HND_IJPAIR/IA(MXNBF)
      COMMON/ HND_MCPAR/NFZC,NCORBS,NCI
      DIMENSION FC(NC,1),DF(LDF,1)
      DIMENSION DDA(1),DDB(1),DDC(1),V(1),NCONF(1),DM2(1),IDM2(1)
      DIMENSION DCMC(1),DVMC(1),VMC(NUM,1)
      DIMENSION DIJKL(1),DMJKL(1),DMNKL(1),DMNPL(1),DMNPQ(1)
      DIMENSION MI(48),MJ(48),MK(48),M0(48)
      DIMENSION GCLB(1),NIJG(2,1)
      DATA TENM10 /1.0D-10/
      DATA TENM09 /1.0D-09/
C
      SPSKIP = .FALSE.
C
      LVEC=0
      NVEC=0
C
C     ----- TWO-ELECTRON INTEGRALS DERIVATIVES -----
C
      LCAP=NAP.GT.1
      DABTOL=TENM10
      CUTOFF=TENM09
      SOME=.TRUE.
C
      NUMI0=0
      NUMJ0=0
      NUMK0=0
      NUML0=0
      IEQJ0=.FALSE.
      KEQL0=.FALSE.
      SAME0=.FALSE.
      FIRST=.FALSE.
C
C     ----- ISHELL -----
C
      IF(IST.GT.NSHELL) GO TO 9100
      IJKLSH=0
      DO 9010 II=IST,NSHELL
      IIMOD=MODSHL(II)
      DO 20 IT=1,NT
      ID=MAPSHL(IIMOD,IT)
      ID=INVSHL(ID)
      IF(ID.GT.II) GO TO 9000
      M0(IT)=ID
   20 CONTINUE
C
      IF(NCI.GT.0.AND.MC)
     1 CALL HND_DABMCI(VMC(1,NCORBS+1),DIJKL,DMJKL,NCI,
     2             IIMOD,NUM,LDM2)
C
C     ----- JSHELL -----
C
      J0=JST
      DO 8000 JJ=J0,II
      JJMOD=MODSHL(JJ)
      JST=1
      DO 50 IT=1,NT
      JD=MAPSHL(JJMOD,IT)
      JD=INVSHL(JD)
      IF(JD.GT.II) GO TO 8000
      ID=M0(IT)
      IF(ID.GE.JD) GO TO 40
      ND=ID
      ID=JD
      JD=ND
   40 IF(ID.EQ.II.AND.JD.GT.JJ) GO TO 8000
      MI(IT)=ID
      MJ(IT)=JD
   50 CONTINUE
C
      IF(NCI.GT.0.AND.MC)
     1 CALL HND_DABMCJ(VMC(1,NCORBS+1),DMJKL,DMNKL,NCI,
     2             IIMOD,JJMOD,NUM,LDM2)
C
C     ----- GET -IJ- CHARGE DISTRIBUTION -----
C
      IIJJ=IA(MAX0(IIMOD,JJMOD))+MIN0(IIMOD,JJMOD)
      NIJ0=NIJG(1,IIJJ)
      NIJ =NIJG(2,IIJJ)
      IF(NIJ.EQ.0) GO TO 8000
      CALL HND_OEDRD(X(ICHRG),NIJ0,ADCHRG,DASI,DASJ,NIJ)
C
C     ----- KSHELL -----
C
      K0=KST
      DO 7000 KK=K0,II
      KKMOD=MODSHL(KK)
      KST=1
      DO 70 IT=1,NT
      KD=MAPSHL(KKMOD,IT)
      KD=INVSHL(KD)
      IF(KD.GT.II) GO TO 7000
   70 MK(IT)=KD
C
      IF(NCI.GT.0.AND.MC)
     1 CALL HND_DABMCK(VMC(1,NCORBS+1),DMNKL,DMNPL,NCI,
     2             IIMOD,JJMOD,KKMOD,NUM,LDM2)

C
C     ----- LSHELL ----
C
      MAXLL=KK
      IF(KK.EQ.II) MAXLL=JJ
      L0=LST
      DO 6000 LL=L0,MAXLL
      LLMOD=MODSHL(LL)
      LST=1
      N4=0
      DO 110 IT=1,NT
      LD=MAPSHL(LLMOD,IT)
      LD=INVSHL(LD)
      IF(LD.GT.II) GO TO 6000
      KD=MK(IT)
      IF(KD.GE.LD) GO TO 90
      ND=KD
      KD=LD
      LD=ND
   90 ID=MI(IT)
      JD=MJ(IT)
      IF(ID.NE.II.AND.KD.NE.II) GO TO 110
      IF(KD.LT.ID) GO TO 100
      IF(KD.EQ.ID.AND.LD.LE.JD) GO TO 100
      ND=ID
      ID=KD
      KD=ND
      ND=JD
      JD=LD
      LD=ND
  100 IF(JD.LT.JJ) GO TO 110
      IF(JD.GT.JJ) GO TO 6000
      IF(KD.LT.KK) GO TO 110
      IF(KD.GT.KK) GO TO 6000
      IF(LD.LT.LL) GO TO 110
      IF(LD.GT.LL) GO TO 6000
      N4=N4+1
  110 CONTINUE
      Q4= DBLE(NT)/ DBLE(N4)
      QQ4=Q4
C
C     ----- SELECT SHELL BLOCK IF PARALLEL PROCESSING -----
C

      IJKLSH=IJKLSH+1
      IF(MOD(IJKLSH,NAP)+1.NE.IAP) GO TO 6000
C
C     ----- CHECK -SP- TYPE -----
C
      SP=KTYPE(IIMOD).LE.2.AND.KTYPE(JJMOD).LE.2.AND.
     1   KTYPE(KKMOD).LE.2.AND.KTYPE(LLMOD).LE.2
      IF(SP.AND.SPSKIP) GO TO 6000
C
C     ----- CHECK INTEGRAL BOUNDS -----
C
      IJIJ=IA(MAX0(IIMOD,JJMOD))+MIN0(IIMOD,JJMOD)
      KLKL=IA(MAX0(KKMOD,LLMOD))+MIN0(KKMOD,LLMOD)

C     GMAX= SQRT(GCLB(IJIJ)*GCLB(KLKL))
C     IF(GMAX.LT.CUTOFF) GO TO 6000

C
C     ----- GET -KL- CHARGE DISTRIBUTION -----
C
      KKLL=IA(MAX0(KKMOD,LLMOD))+MIN0(KKMOD,LLMOD)
      NKL0=NIJG(1,KKLL)
      NKL =NIJG(2,KKLL)
      IF(NKL.EQ.0) GO TO 6000
      CALL HND_OEDRD(X(ICHRG),NKL0,BDCHRG,DBSK,DBSL,NKL)
C
C     ----- SELECT CENTERS FOR DERIVATIVES -----
C
      ISH=IIMOD
      JSH=JJMOD
      KSH=KKMOD
      LSH=LLMOD
      CALL HND_JKDATM(ISH,JSH,KSH,LSH)
      IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) GO TO 6000
C
C     ----- SET INDICES FOR SHELL BLOCK -----
C
      CALL HND_JKDSHL(ISH,JSH,KSH,LSH)
      CALL HND_JKDNDX(X(IIJKLG))
C
C     ----- FORM PRODUCTS OF DENSITY MATRIX ELEMENTS -----
C           SKIP CALCULATION IF DABMAX.LT.DABTOL
C

      IF(.NOT.MCCI.AND..NOT.ROGVB)
     1 CALL HND_DABCLU(ISH,JSH,KSH,LSH,SUHF,IA,DDA,DDB,X(IDAB),
     1 DABMAX)
C
      IF(.NOT.MCCI.AND.ROHF)
     1 CALL HND_DABROH(ISH,JSH,KSH,LSH,NORB,NOCORE,NOOPEN,
     2 IA,DDC,V,NCONF,NUM,X(IDAB),DABMAX)
C
      IF(.NOT.MCCI.AND.RGVB)
     1 CALL HND_DABGVB(ISH,JSH,KSH,LSH,NORB,NOCORE,NOOPEN,NOPAIR,
     2 IA,DDC,V,NCONF,NUM,X(IDAB),DABMAX)
C
      IF(CI ) CALL HND_DABCI(ISH,JSH,KSH,LSH,NBUF,DM2,IDM2,X(IDAB),
     1 DABMAX)

      IF(MP2) CALL HND_DABMP(ISH,JSH,KSH,LSH,DM2,V,DDA,DDB,X(IDAB),
     1 DABMAX,
     1                   LDF,NBF)
C
      IF(MC ) THEN
       IF(NCI.GT.0)
     1 CALL HND_DABMCL(VMC(1,NCORBS+1),DMNPL,DMNPQ,NCI,
     2             IIMOD,JJMOD,KKMOD,LLMOD,NUM,LDM2)
       CALL HND_DABMC(ISH,JSH,KSH,LSH,DCMC,DVMC,DMNPQ,IA,NCORBS,NCI,
     1 LDM2,
     1            X(IDAB),DABMAX)
      ENDIF
C
      IF(DABMAX.LT.DABTOL) GO TO 6000
C
C
C     ----- (II,JJ//KK,LL) -----
C
      CALL HND_JKDSPD(X(IGINT),X(IFINT),X(ISINT),X(IIJKLG),
     1 X(IGIJKL),X(IGNKL),X(IGNM),
     2 X(IXY ),X(IXZ ),X(IYZ ),X(IX  ),X(IY  ),X(IZ  ),
     3 X(ISJ ),X(ISK ),X(ISL ),
     4 X(IB00),X(IB01),X(IB10),X(IC00),X(ID00),X(IF00),
     5 X(IDIJ),X(IDKL),X(IDIJSI),X(IDIJSJ),X(IDKLSK),X(IDKLSL),
     6 X(IABV),X(ICV ),X(IRW ),
     7 X(IAAI),X(IAAJ),X(IBBK),X(IBBL),
     8 X(IFI ),X(IFJ ),X(IFK ),X(IFL ),
     9 X(ISII),X(ISJJ),X(ISKK),X(ISLL),
     1 X(ISIJ),X(ISIK),X(ISIL),X(ISJK),X(ISJL),X(ISKL),
     2 X(IDAB),MAXXYZ,FC,NC,DF,LDF,NBF,DDA)
C
C
 6000 CONTINUE
 7000 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
C
      CALL HND_JKDEND(II)
 9010 CONTINUE
 9100 CONTINUE
C
c%% both # of vector and vector length are wrong ?
      AVGVEC=DBLE(LVEC)/DBLE(NVEC)
C     IF(SOME.and.ga_nodeid().eq.0) WRITE(IW,9999) NVEC,AVGVEC
      RETURN
C9999 FORMAT(' # OF VECTORS = ',I10,' AVERAGE VECTOR LENGTH = ',F10.1)
      END
C
C     ------SUBROUTINES ABOVE TAKEN FROM DER.F IN HONDO-------
C


C
C     SUBROUTINES TAKEN FROM CTL.F IN HONDO -------
C
      SUBROUTINE HND_PREAD(IJK,XX,IX,NXX,MAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PACK2E
      COMMON/HND_PCKLAB/LABSIZ
      COMMON/HND_PCKOPT/NHEX,NTUPL,PACK2E
      COMMON/HND_INTOPT/NOPK,NOK,NOSQUR
      COMMON/HND_INTFIL/NINTMX
      DIMENSION XX(MAX),IX(LABSIZ*MAX)
      IF(PACK2E) GO TO 10
      READ(IJK) XX,IX,NX
      NXX=NX
      RETURN
   10 CALL HND_READPK(IJK,XX,XX,       NXX,NH,IERR,IEND)
      CALL HND_READPN(IJK,IX,IX,LABSIZ*NXX,NT,IERR,IEND)
      RETURN
      END
      SUBROUTINE HND_READPN(IS,IX,IY,NNX,NT,IERR,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-READPN-'/
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
      SUBROUTINE HND_WRITPK(IS,XX,YY,N,NHEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-WRITPK-'/
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
C
      SUBROUTINE HND_ISOOUT(NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXSYM =48)
      INTEGER ZSHIFTR
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_ISOPAC/INDIN(MXSYM),INDOUT(MXSYM/2)
C
#include "bitops.fh"
C
      MASK1(IARG)=2**IARG-1
      ZSHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      ISOPCK(IARG,IBIT)=LAND(ZSHIFTR(IARG,IBIT),MASK1(16))
      IBIT1=16
C
      DO 10 IT=1,NT
   10 INDIN(IT)=0
      IT  =0
      NOUT=0
      MAX =0
   20 MIN =MAX+1
      MINP=MIN
      MAX =MAX+(4/ISINGL)
      MAXP=MAX
      IF(MAX.GT.NT) MAX=NT
      NOUT=NOUT+1
      IPACK=INDOUT(NOUT)
      DO 30 I=MINP,MAXP
      IBIT=IBIT1*(MAXP-I)
      IDUM=ISOPCK(IPACK,IBIT)
      IF(IDUM.EQ.0) GO TO 30
      IT=IT+1
      INDIN(IT)=IDUM
   30 CONTINUE
C
      IF(MAX.LT.NT) GO TO 20
      RETURN
      END
C
C     ------ SUBROUTINES TAKEN FROM DER.F IN HONDO ------
C 
      SUBROUTINE HND_OEDSHL(ISH,JSH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXGSH =30,MXGSH2=MXGSH*MXGSH)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL EXPNDI,EXPNDK
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/     
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLXPN/EXPNDI,EXPNDK
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLINF/GA(MXGSH),CCA(MXGSH),CCAS(MXGSH),
     1              GB(MXGSH),CCB(MXGSH),CCBS(MXGSH),
     1              GC(MXGSH),CCC(MXGSH),CCCS(MXGSH),
     1              GD(MXGSH),CCD(MXGSH),CCDS(MXGSH),
     1              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     1              NGA,NGB,NGC,NGD
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      DIMENSION CSPDFG(2048,5)
      EQUIVALENCE (CSPDFG(1,1),CS(1))
C
      IIEQJJ=ISH.EQ.JSH
C
C     ----- ISHELL -----
C
      I=KATOM(ISH)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(ISH)
      I2=I1+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
      NGA=0
      DO 10 I=I1,I2
      NGA=NGA+1
      GA(NGA)=EX(I)
      CCA(NGA)=CSPDFG(I,LIT)
      IF(SPI) CCAS(NGA)=CSPDFG(I,1)/CSPDFG(I,2)
   10 CONTINUE
C
C     ----- JSHELL -----
C
      J=KATOM(JSH)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JSH)
      J2=J1+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      NGB=0
      DO 20 J=J1,J2
      NGB=NGB+1
      GB(NGB)=EX(J)
      CCB(NGB)=CSPDFG(J,LJT)
      IF(SPJ) CCBS(NGB)=CSPDFG(J,1)/CSPDFG(J,2)
   20 CONTINUE
      RRI=((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      SPIJ=SPI.OR.SPJ
      EXPNDI=LIT.GE.LJT
      RETURN
      END
C
      SUBROUTINE HND_OEDCHR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXGSH =30,MXGSH2=MXGSH*MXGSH)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL EXPNDI,EXPNDK
      COMMON/HND_SHLTOL/RTOL,DTOL
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLXPN/EXPNDI,EXPNDK
      COMMON/HND_SHLINF/GA(MXGSH),CCA(MXGSH),CCAS(MXGSH),
     1              GB(MXGSH),CCB(MXGSH),CCBS(MXGSH),
     1              GC(MXGSH),CCC(MXGSH),CCCS(MXGSH),
     1              GD(MXGSH),CCD(MXGSH),CCDS(MXGSH),
     1              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     1              NGA,NGB,NGC,NGD
      COMMON/HND_IJDCHR/ADCHRG(13,MXGSH2),DASI(MXGSH2),DASJ(MXGSH2),
     1 NIJ
      DATA ONE /1.0D+00/
C
C     ----- -IJ- CHARGE DISTRIBUTION -----
C
      XC=XI
      YC=YI
      ZC=ZI
      DXIJ=XI-XJ
      DYIJ=YI-YJ
      DZIJ=ZI-ZJ
      IF(EXPNDI) GO TO 10
      XC=XJ
      YC=YJ
      ZC=ZJ
      DXIJ=XJ-XI
      DYIJ=YJ-YI
      DZIJ=ZJ-ZI
   10 CONTINUE
C
C     ----- - I- PRIMITIVE           -----
C
      NIJ=0
      DO 30 IA=1,NGA
      AI=GA(IA)
      ARRI=AI*RRI
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CCI=CCA(IA)
C
C     ----- - J- PRIMITIVE           -----
C
      DO 20 JB=1,NGB
      AJ=GB(JB)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.RTOL) GO TO 20
      DAEXPA=CCI*CCB(JB)* EXP(-DUM)*AA1
      DUM=  ABS(DAEXPA)
      IF(DUM.LE.DTOL) GO TO 20
C
      NIJ=NIJ+1
      ADCHRG( 1,NIJ)= DAEXPA
      ADCHRG( 2,NIJ)= AA
      ADCHRG( 3,NIJ)=(AXI+AJ*XJ)*AA1
      ADCHRG( 4,NIJ)=(AYI+AJ*YJ)*AA1
      ADCHRG( 5,NIJ)=(AZI+AJ*ZJ)*AA1
      ADCHRG( 6,NIJ)= XC
      ADCHRG( 7,NIJ)= YC
      ADCHRG( 8,NIJ)= ZC
      ADCHRG( 9,NIJ)= DXIJ
      ADCHRG(10,NIJ)= DYIJ
      ADCHRG(11,NIJ)= DZIJ
      ADCHRG(12,NIJ)= AI+AI
      ADCHRG(13,NIJ)= AJ+AJ
      IF(SPI) DASI(NIJ)=CCAS(IA)
      IF(SPJ) DASJ(NIJ)=CCBS(JB)
C
   20 CONTINUE
   30 CONTINUE
      RETURN
      END
C
C     ------- SUBROUTINE TAKEN FROM CTL.F IN HONDO ------
C
      SUBROUTINE HND_WRTGRD(G,N1,N2,N3,N4,N5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXTEMP=MXCART+5)
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_MACHIN/ISINGL,NBITS
      DIMENSION G(MXCART),X(MXTEMP),IX(1)
      EQUIVALENCE (IX(1),X(1))
      DO I=1,MXCART
         X(I)=G(I)
      ENDDO
      I0=MXCART*ISINGL
      IX(1+I0)=N1
      IX(2+I0)=N2
      IX(3+I0)=N3
      IX(4+I0)=N4
      IX(5+I0)=N5
      NX=MXTEMP
      CALL HND_DAWRIT(IDAF,IODA,X,NX,3,NAV)
      RETURN
      END
C
      SUBROUTINE HND_DAWRIT(IDAF,IODA,IX,NX,IDAR,NAV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      PARAMETER (MXIODA=255)
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_DAFNAV/IXDDAF(2048),NAV10,NAV20
      COMMON/HND_DAFREC/LENDAR(MXIODA)
      DIMENSION ERRMSG(3)
      DIMENSION IODA(1),IX(1)
      DIMENSION IXSDAF(1024)
      EQUIVALENCE (IXSDAF(1),IXDDAF(1))
      DATA ERRMSG /'PROGRAM ','STOP IN ','-DAWRIT-'/
      DATA IXMAX  /1024/
C
      OUT=.FALSE.
C
      MAXIX=IXMAX*ISINGL
      LDAR =   NX*ISINGL
      IF(IODA(IDAR+0).NE.0) GO TO 100
C
C     ----- FIRST WRITE -----
C
      IODA(IDAR+0)=NAV
      IF(IDAF.EQ.10) NAV10=NAV
      IF(IDAF.EQ.20) NAV20=NAV
C
      IF(IDAF.EQ.10) THEN
         LENDAR(IDAR+0)=NX
         IF(OUT) THEN
            WRITE(IW,9998) IDAR,NX
         ENDIF
      ENDIF
C
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      DO 20 I=MIN,MAX
   20 IXDDAF(I-MIN+1)=IX(I)
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV20) IXDDAF
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 10
C
      IF(IDAF.EQ.10) NAV=NAV10
      IF(IDAF.EQ.20) NAV=NAV20
      RETURN
C
C     ----- REWRITE -----
C
  100 CONTINUE
      IF(IDAF.EQ.10) NAV10=IODA(IDAR+0)
      IF(IDAF.EQ.20) NAV20=IODA(IDAR+0)
C
      IF(IDAF.EQ.10) THEN
         IF(OUT) THEN
            WRITE(IW,9997) IDAR,NX
         ENDIF
         NX0=LENDAR(IDAR+0)
         IF(NX.GT.NX0) THEN
            IF(OUT) THEN
               WRITE(IW,9999) IDAR,NX,NX0
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ELSEIF(NX.LT.NX0) THEN
            IF(OUT) THEN
               WRITE(IW,9999) IDAR,NX,NX0
            ENDIF
         ENDIF
      ENDIF
C
      MAX=0
  110 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      DO 120 I=MIN,MAX
  120 IXDDAF(I-MIN+1)=IX(I)
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV20) IXDDAF
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 110
C
      RETURN
 9999 FORMAT(' INCONSISTENT RECORD LENGTH FOR -IDAR- = ',I4,/,
     1       ' NX,NX0 = ',2I10)
 9998 FORMAT(' FIRST-WRITE, IDAR = ',I5,' LENGTH = ',I10)
 9997 FORMAT('    RE-WRITE, IDAR = ',I5,' LENGTH = ',I10)
      END
C

C
C     ------ SUBROUTINES TAKEN FROM DER.F IN HONDO ------
C
      SUBROUTINE HND_DABMCI(V,GIN,GOUT,NCI,II,NUM,LDM2)      
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      DIMENSION V(NUM,1)
      DIMENSION GIN(NCI,NCI,NCI,NCI)
      DIMENSION GOUT(LDM2,NCI,NCI,NCI)
      DATA ZERO /0.0D+00/
C
C     ----- PARTIAL TRANSFORMATION OF G(IJKL) TO G(MU,J,K,L) -----
C
      IX = KLOC(II)-KMIN(II)
      DO 30 L = 1, NCI
      DO 30 K = 1, NCI
      DO 30 J = 1, NCI
        DO 20 M = KMIN(II),KMAX(II)
        IN = IX + M
        DUM  = ZERO
        DO 10 I = 1, NCI
   10     DUM = DUM + GIN(I,J,K,L)* V(IN,I)
   20   GOUT(M,J,K,L)=DUM
   30 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE HND_DABMCJ(V,GIN,GOUT,NCI,II,JJ,NUM,LDM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      DIMENSION V(NUM,1)
      DIMENSION GIN(LDM2,NCI,NCI,NCI)
      DIMENSION GOUT(LDM2,LDM2,NCI,NCI)
      DATA ZERO /0.0D+00/
C
C     ----- PARTIAL TRANSFORMATION OF G(MU,JKL) IN G(MU,NU,K,L) -----
C
      JX = KLOC(JJ)-KMIN(JJ)
      DO 30 L = 1, NCI
      DO 30 K = 1, NCI
      DO 30 M = KMIN(II),KMAX(II)
        DO 20 N = KMIN(JJ),KMAX(JJ)
        JN = JX + N
        DUM  = ZERO
        DO 10 J = 1, NCI
   10     DUM = DUM + GIN(M,J,K,L) * V(JN,J)
   20   GOUT(M,N,K,L)=DUM
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_DABMCK(V,GIN,GOUT,NCI,II,JJ,KK,NUM,LDM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      DIMENSION V(NUM,1)
      DIMENSION GIN(LDM2,LDM2,NCI,NCI)
      DIMENSION GOUT(LDM2,LDM2,LDM2,NCI)
      DATA ZERO /0.0D+00/
C
C     ----- PARTIAL TRANSFORMATION OF G(M,N,KL) IN G(M,N,P,L) -----
C
      KX = KLOC(KK)-KMIN(KK)
      DO 30 L = 1, NCI
      DO 30 N = KMIN(JJ),KMAX(JJ)
      DO 30 M = KMIN(II),KMAX(II)
        DO 20 NP = KMIN(KK),KMAX(KK)
        KN = KX + NP
        DUM  = ZERO
        DO 10 K = 1, NCI
   10     DUM = DUM + GIN(M,N,K,L) * V(KN,K)
   20   GOUT(M,N,NP,L)=DUM
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_OEDRD(DCHRG,NIJ0,ADCHRG,DASI,DASJ,NIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ADCHRG(13,1),DASI(1),DASJ(1)
      DIMENSION  DCHRG(15,1)
C
      DO 20 IJ=1,NIJ
      DO 10 I=1,13
   10 ADCHRG(I,IJ)=DCHRG( I,IJ+NIJ0)
      DASI(IJ)    =DCHRG(14,IJ+NIJ0)
      DASJ(IJ)    =DCHRG(15,IJ+NIJ0)
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKDATM(II,JJ,KK,LL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL DBUG
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IANDJ,IANDK,IANDL,JANDK,JANDL,KANDL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERINV/INVTYP
      COMMON/HND_DSHPRT/OUT
C
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/
C
      SKIPI=.TRUE.
      SKIPJ=.TRUE.
      SKIPK=.TRUE.
      SKIPL=.TRUE.
      IIAT=KATOM(II)
      JJAT=KATOM(JJ)
      KKAT=KATOM(KK)
      LLAT=KATOM(LL)
      IANDJ=IIAT.EQ.JJAT
      IANDK=IIAT.EQ.KKAT
      IANDL=IIAT.EQ.LLAT
      JANDK=JJAT.EQ.KKAT
      JANDL=JJAT.EQ.LLAT
      KANDL=KKAT.EQ.LLAT
C
      IF(.NOT.IANDJ) GO TO 500
      IF(.NOT.IANDK) GO TO 200
      IF(.NOT.IANDL) GO TO 100
C
C     ----- IAT = JAT = KAT = LAT ----- (IAT,IAT/IAT,IAT) -----
C
      INVTYP=1
      GO TO 1500
  100 CONTINUE
C
C     ----- IAT = JAT = KAT ; LAT ----- (IAT,IAT/IAT,LAT) -----
C
      SKIPL=.FALSE.
      INVTYP=2
      GO TO 1500
  200 IF(.NOT.IANDL) GO TO 300
C
C     ----- IAT = JAT = LAT ; KAT ----- (IAT,IAT/KAT,IAT) -----
C
      SKIPK=.FALSE.
      INVTYP=3
      GO TO 1500
  300 IF(.NOT.KANDL) GO TO 400
C
C     ----- IAT = JAT ; KAT = LAT ----- (IAT,IAT/KAT,KAT) -----
C
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=4
      GO TO 1500
  400 CONTINUE
C
C     ----- IAT = JAT ; KAT ; LAT ----- (IAT,IAT/KAT,LAT) -----
C
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=5
      GO TO 1500
  500 IF(.NOT.IANDK) GO TO 800
      IF(.NOT.IANDL) GO TO 600
C
C     ----- IAT = KAT = LAT ; JAT ----- (IAT,JAT/IAT,IAT) -----
C
      SKIPJ=.FALSE.
      INVTYP=6
      GO TO 1500
  600 IF(.NOT.JANDL) GO TO 700
C
C     ----- IAT = KAT ; JAT = LAT ----- (IAT,JAT/IAT,JAT) -----
C
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=7
      GO TO 1500
  700 CONTINUE
C
C     ----- IAT = KAT ; JAT ; LAT ----- (IAT,JAT/IAT,LAT) -----
C
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=8
      GO TO 1500
  800 IF(.NOT.IANDL) GO TO 1000
      IF(.NOT.JANDK) GO TO 900
C
C     ----- IAT = LAT ; JAT = KAT ----- (IAT,JAT/JAT,IAT) -----
C
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=9
      GO TO 1500
  900 CONTINUE
C
C     ----- IAT = LAT ; JAT , KAT ----- (IAT,JAT/KAT,IAT) -----
C
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=10
      GO TO 1500
 1000 IF(.NOT.JANDK) GO TO 1200
      IF(.NOT.JANDL) GO TO 1100
C
C     ----- IAT ; JAT = JAT = JAT ----- (IAT,JAT/JAT,JAT) -----
C
      SKIPI=.FALSE.
      INVTYP=11
      GO TO 1500
 1100 CONTINUE
C
C     ----- IAT ; JAT = KAT ; LAT ----- (IAT,JAT/JAT,LAT) -----
C
      SKIPI=.FALSE.
      SKIPL=.FALSE.
      INVTYP=12
      GO TO 1500
 1200 IF(.NOT.JANDL) GO TO 1300
C
C     ----- JAT = LAT ; IAT ; KAT ----- (IAT,JAT/KAT,JAT) -----
C
      SKIPI=.FALSE.
      SKIPK=.FALSE.
      INVTYP=13
      GO TO 1500
 1300 IF(.NOT.KANDL) GO TO 1400
C
C     ----- KAT = LAT ; IAT ; JAT ----- (IAT,JAT/KAT,KAT) -----
C
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      INVTYP=14
      GO TO 1500
 1400 CONTINUE
C
C     ----- IAT ; JAT ; KAT ; LAT ----- (IAT,JAT/KAT,LAT) -----
C
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=15
C
 1500 CONTINUE
C
C     ----- FOR DEBUGGING PURPOSES CALCULATE ALL TERMS -----
C
      DBUG=.FALSE.
      IF(.NOT.DBUG) GO TO 1600
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=16
C
 1600 CONTINUE
C
      IF(.NOT.OUT) RETURN
      WRITE(IW,9999) II,JJ,KK,LL,SKIPI,SKIPJ,SKIPK,SKIPL
      RETURN
 9999 FORMAT(/,' -- II,JJ,KK,LL =',4I3,' SKIP1,2,3,4 =',4L3)
      END
C
      SUBROUTINE HND_JKDSHL(ISH,JSH,KSH,LSH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL EXPNDI,EXPNDK
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLXPN/EXPNDI,EXPNDK
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
C
      IIEQJJ=ISH.EQ.JSH
      KKEQLL=KSH.EQ.LSH
      IJEQKL=ISH.EQ.KSH.AND.JSH.EQ.LSH
      IJGTKL=MAX0(ISH,JSH).GT.MAX0(KSH,LSH)
      IJLTKL=MAX0(ISH,JSH).LT.MAX0(KSH,LSH)
C
C     ----- ISHELL -----
C
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
C
C     ----- JSHELL -----
C
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      SPIJ=SPI.OR.SPJ
      EXPNDI=LIT.GE.LJT
C
C     ----- KSHELL -----
C
      LKT=KTYPE(KSH)
      MINK=KMIN(KSH)
      MAXK=KMAX(KSH)
      NUMK=MAXK-MINK+1
      LOCK=KLOC(KSH)-MINK
      SPK=LKT.EQ.2.AND.MINK.EQ.1
C
C     ----- LSHELL -----
C
      LLT=KTYPE(LSH)
      MINL=KMIN(LSH)
      MAXL=KMAX(LSH)
      NUML=MAXL-MINL+1
      LOCL=KLOC(LSH)-MINL
      SPL=LLT.EQ.2.AND.MINL.EQ.1
      SPKL=SPK.OR.SPL
      SPIJKL=SPIJ.OR.SPKL
      EXPNDK=LKT.GE.LLT
C
      RETURN
      END
C
      SUBROUTINE HND_JKDNDX(IJKLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL IEQJ0,KEQL0,SAME0,FIRST
      LOGICAL SAMTYP,SAMSYM
      COMMON/HND_HFK/XX,U(9),W(9),NROOTS
      COMMON/HND_RYS/YY,T(5),V(5),MROOTS
      COMMON/HND_SHLBAS/MAXTYP,MAXNUM
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),LGXYZ(4,35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLOLD/NUMI0,NUMJ0,NUMK0,NUML0,
     1              IEQJ0,KEQL0,SAME0,FIRST
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      DIMENSION IJKLX(35),IJKLY(35),IJKLZ(35)
      DIMENSION IJKLN(5)
      DIMENSION IJKLG(4,1)
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN HONDO -----
C
      DATA IJKLN /   1,  4, 10, 20, 35/
      DATA IJKLX /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     1               3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     2               4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     3               2,  0,  2,  1,  1/
      DATA IJKLY /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     1               0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     2               0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     3               0,  2,  1,  2,  1/
      DATA IJKLZ /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     1               0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     2               0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     3               2,  2,  1,  1,  2/
C
C     ----- ORDER OF COMPONENTS IN SHELLS IN NWCHEM -----
C
C     DATA IJKLN /   1,  4, 10, 20, 35/
C     DATA IJKLX /   0,  1,  0,  0,  2,  1,  1,  0,  0,  0,
C    1               3,  2,  2,  1,  1,  1,  0,  0,  0,  0,
C    2               4,  3,  3,  2,  2,  2,  1,  1,  1,  1,
C    3               0,  0,  0,  0,  0/
C     DATA IJKLY /   0,  0,  1,  0,  0,  1,  0,  2,  1,  0,
C    1               0,  1,  0,  2,  1,  0,  3,  2,  1,  0,
C    2               0,  1,  0,  2,  1,  0,  3,  2,  1,  0,
C    3               4,  3,  2,  1,  0/
C     DATA IJKLZ /   0,  0,  0,  1,  0,  0,  1,  0,  1,  2,
C    1               0,  0,  1,  0,  1,  2,  0,  1,  2,  3,
C    2               0,  0,  1,  0,  1,  2,  0,  1,  2,  3,
C    3               0,  1,  2,  3,  4/
C
      SAMTYP=NUMI.EQ.NUMI0.AND.NUMJ.EQ.NUMJ0.AND.
     1       NUMK.EQ.NUMK0.AND.NUML.EQ.NUML0          .AND.(.NOT.FIRST)
      SAMSYM=(IIEQJJ.EQV.IEQJ0).AND.
     1       (KKEQLL.EQV.KEQL0).AND.(IJEQKL.EQV.SAME0).AND.(.NOT.FIRST)
      NUMI0=NUMI
      NUMJ0=NUMJ
      NUMK0=NUMK
      NUML0=NUML
      IEQJ0=IIEQJJ
      KEQL0=KKEQLL
      SAME0=IJEQKL
      FIRST=.FALSE.
C.... IF(SAMTYP.AND.SAMSYM) RETURN
C
      IDER=NDER
      JDER=NDER
      KDER=NDER
      LDER=NDER
      IF(SKIPI) IDER=0
      IF(SKIPJ) JDER=0
      IF(SKIPK) KDER=0
      IF(SKIPL) LDER=0
      LITMOD=LIT + IDER
      LJTMOD=LJT + JDER
      LKTMOD=LKT + KDER
      LLTMOD=LLT + LDER
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      NI=NUML*NUMK*NUMJ
      DO 10 I=MINI,MAXI
   10 IGXYZ(1,I)=NI*(I-MINI)+1
      LLKJT=LLTMOD*LKTMOD*LJTMOD
      DO 20 I=1,IJKLN(LIT)
      IGXYZ(2,I)=IJKLX(I)*LLKJT+1
      IGXYZ(3,I)=IJKLY(I)*LLKJT+1
   20 IGXYZ(4,I)=IJKLZ(I)*LLKJT+1
C
      NJ=NUML*NUMK
      DO 30 J=MINJ,MAXJ
   30 JGXYZ(1,J)=NJ*(J-MINJ)
      LLKT=LLTMOD*LKTMOD
      DO 40 J=1,IJKLN(LJT)
      JGXYZ(2,J)=IJKLX(J)*LLKT
      JGXYZ(3,J)=IJKLY(J)*LLKT
   40 JGXYZ(4,J)=IJKLZ(J)*LLKT
C
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS -----
C
      NK=NUML
      DO 110 K=MINK,MAXK
  110 KGXYZ(1,K)=NK*(K-MINK)
      DO 120 K=1,IJKLN(LKT)
      KGXYZ(2,K)=IJKLX(K)*(LLTMOD)
      KGXYZ(3,K)=IJKLY(K)*(LLTMOD)
  120 KGXYZ(4,K)=IJKLZ(K)*(LLTMOD)
C
      NL=1
      DO 130 L=MINL,MAXL
  130 LGXYZ(1,L)=NL*(L-MINL)
      DO 140 L=1,IJKLN(LLT)
      LGXYZ(2,L)=IJKLX(L)
      LGXYZ(3,L)=IJKLY(L)
  140 LGXYZ(4,L)=IJKLZ(L)
C
C     ----- PREPARE INDICES FOR (IJ/KL) -----
C
  200 CONTINUE
C
      IJKL=0
      DO 240 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 230 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 220 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 210 L=MINL,LMAX
      IJKL=IJKL+1
      NN=((IGXYZ(1,I)+JGXYZ(1,J))+KGXYZ(1,K))+LGXYZ(1,L)
      NX=((IGXYZ(2,I)+JGXYZ(2,J))+KGXYZ(2,K))+LGXYZ(2,L)
      NY=((IGXYZ(3,I)+JGXYZ(3,J))+KGXYZ(3,K))+LGXYZ(3,L)
      NZ=((IGXYZ(4,I)+JGXYZ(4,J))+KGXYZ(4,K))+LGXYZ(4,L)
      IJKLG(1,IJKL)=   NN
      IJKLG(2,IJKL)=3*(NX-1)+1
      IJKLG(3,IJKL)=3*(NY-1)+2
      IJKLG(4,IJKL)=3*(NZ-1)+3
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
C
C     ----- SET NUMBER OF QUADRATURE POINTS -----
C
      NROOTS=(LIT+LJT+LKT+LLT-2 + NDER )/2
      MROOTS=NROOTS
C
      RETURN
      END
C
      SUBROUTINE HND_DABCLU(II,JJ,KK,LL,UHFTYP,IA,DA,DB,DAB,DABMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL UHFTYP
      LOGICAL OUT
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1LGXYZ(4,35)
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_DSHPRT/OUT
      DIMENSION IA(1),DA(1),DB(1),DAB(1)
      DATA ZERO     /0.0D+00/
      DATA PT5,FOUR /0.5D+00,4.0D+00/
C
      MINI=KMIN(II)
      MINJ=KMIN(JJ)
      MINK=KMIN(KK)
      MINL=KMIN(LL)
      MAXI=KMAX(II)
      MAXJ=KMAX(JJ)
      MAXK=KMAX(KK)
      MAXL=KMAX(LL)
      LOCI=KLOC(II)-MINI
      LOCJ=KLOC(JJ)-MINJ
      LOCK=KLOC(KK)-MINK
      LOCL=KLOC(LL)-MINL
C
      DABMAX=ZERO
      DO 500 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 400 J=MINJ,JMAX
      PIJ=PI*PNRM(J)
       MAX=MAXK
      IF(IJEQKL)  MAX=I
      DO 300 K=MINK, MAX
      PIJK=PIJ*PNRM(K)
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 200 L=MINL,LMAX
      PIJKL=PIJK*PNRM(L)
      I1=LOCI+I
      I2=LOCJ+J
      I3=LOCK+K
      I4=LOCL+L
      IF(I1.GE.I2) GO TO 10
      N=I1
      I1=I2
      I2=N
   10 IF(I3.GE.I4) GO TO 20
      N=I3
      I3=I4
      I4=N
   20 IF(I1-I3) 30,40,50
   30 N=I1
      I1=I3
      I3=N
      N=I2
      I2=I4
      I4=N
      GO TO 50
   40 IF(I2.LT.I4) GO TO 30
   50 CONTINUE
      MIJ=IA(I1)+I2
      MIK=IA(I1)+I3
      MIL=IA(I1)+I4
      MKL=IA(I3)+I4
      IF(I2.LT.I3) GO TO 60
      MJK=IA(I2)+I3
      MJL=IA(I2)+I4
      GO TO 80
   60 MJK=IA(I3)+I2
      IF(I2.LT.I4) GO TO 70
      MJL=IA(I2)+I4
      GO TO 80
   70 MJL=IA(I4)+I2
   80 CONTINUE
      DFAC= DA(MIJ)*DA(MKL)*FOUR
     1     -DA(MIK)*DA(MJL)
     2     -DA(MIL)*DA(MJK)
      IF(UHFTYP) DFAC=DFAC-DB(MIK)*DB(MJL)
     1                    -DB(MIL)*DB(MJK)
      IF(I2.EQ.I1             ) DFAC=DFAC*PT5
      IF(I4.EQ.I3             ) DFAC=DFAC*PT5
      IF(I3.EQ.I1.AND.I4.EQ.I2) DFAC=DFAC*PT5
      IF( ABS(DFAC).GT.DABMAX) DABMAX= ABS(DFAC)
      NN=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
      DAB(NN)=DFAC*PIJKL
      IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,I,J,K,L,NN,DFAC
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
      RETURN
 9999 FORMAT(' -DABCLU- ',4I4,4I3,I5,D20.12)
      END
C
      SUBROUTINE HND_DABROH(II,JJ,KK,LL,NORB,NOCORE,NOOPEN,
     1 IA,DC,V,NCONF,NDIM,DAB,DABMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL NOCORE,NOOPEN
      LOGICAL OUT
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1LGXYZ(4,35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      DIMENSION IA(1),DC(1),V(NDIM,1),NCONF(1),DAB(1)
      DATA ZERO,PT5,TWO,EIGHT /0.0D+00,0.5D+00,2.0D+00,8.0D+00/
      DATA PT125,PT25,ONE     /0.125D+00,0.25D+00,1.0D+00/
C
      MINI=KMIN(II)
      MINJ=KMIN(JJ)
      MINK=KMIN(KK)
      MINL=KMIN(LL)
      MAXI=KMAX(II)
      MAXJ=KMAX(JJ)
      MAXK=KMAX(KK)
      MAXL=KMAX(LL)
      LOCI=KLOC(II)-MINI
      LOCJ=KLOC(JJ)-MINJ
      LOCK=KLOC(KK)-MINK
      LOCL=KLOC(LL)-MINL
C
      DABMAX=ZERO
      DO 500 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 400 J=MINJ,JMAX
      PIJ=PI*PNRM(J)
       MAX=MAXK
      IF(IJEQKL)  MAX=I
      DO 300 K=MINK, MAX
      PIJK=PIJ*PNRM(K)
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 200 L=MINL,LMAX
      PIJKL=PIJK*PNRM(L)
      I1=LOCI+I
      I2=LOCJ+J
      I3=LOCK+K
      I4=LOCL+L
      IF(I1.GE.I2) GO TO 10
      N=I1
      I1=I2
      I2=N
   10 IF(I3.GE.I4) GO TO 20
      N=I3
      I3=I4
      I4=N
   20 IF(I1-I3) 30,40,50
   30 N=I1
      I1=I3
      I3=N
      N=I2
      I2=I4
      I4=N
      GO TO 50
   40 IF(I2.LT.I4) GO TO 30
   50 CONTINUE
      MIJ=IA(I1)+I2
      MIK=IA(I1)+I3
      MIL=IA(I1)+I4
      MKL=IA(I3)+I4
      IF(I2.LT.I3) GO TO 60
      MJK=IA(I2)+I3
      MJL=IA(I2)+I4
      GO TO 80
   60 MJK=IA(I3)+I2
      IF(I2.LT.I4) GO TO 70
      MJL=IA(I2)+I4
      GO TO 80
   70 MJL=IA(I4)+I2
   80 CONTINUE
      DFAC = ZERO
      IF(NOCORE) GO TO 100
      DFAC = DFAC + TWO* DC(MIJ)*DC(MKL)
     1            - PT5*(DC(MIK)*DC(MJL) + DC(MIL)*DC(MJK))
      IF(NOOPEN) GO TO 120
      NCO1 = NCO + 1
      DO 90 IO = NCO1,NORB
      DFAC = DFAC +     (DC(MIJ)*V(I3,IO)*V(I4,IO)
     1                  +DC(MKL)*V(I1,IO)*V(I2,IO))
     1            -PT25*(V(I2,IO)*(DC(MIK)*V(I4,IO)
     1                            +DC(MIL)*V(I3,IO))
     1                  +V(I1,IO)*(DC(MJL)*V(I3,IO)
     1                            +DC(MJK)*V(I4,IO)))
   90 CONTINUE
  100 CONTINUE
      IF(NOOPEN) GO TO 120
      NCO1 = NCO + 1
      DO 110  IO = NCO1,NORB
      IOF=NCONF(IO)-1
      DO 110  JO = NCO1,NORB
      JOF=NCONF(JO)-1
      IOF0=MAX0(IOF,JOF)
      JOF0=MIN0(IOF,JOF)
      IOJO=IA(IOF0)+JOF0
      DFAC = DFAC +
     1                    PT5 *V(I1,IO)*V(I2,IO)*V(I3,JO)*V(I4,JO)
     1 +PT125*(BETA(IOJO)-ONE)*V(I1,IO)*V(I2,JO)*
     1                        (V(I3,IO)*V(I4,JO)+V(I4,IO)*V(I3,JO))
  110 CONTINUE
  120 CONTINUE
C
      IF(I2.EQ.I1             ) DFAC=DFAC*PT5
      IF(I4.EQ.I3             ) DFAC=DFAC*PT5
      IF(I3.EQ.I1.AND.I4.EQ.I2) DFAC=DFAC*PT5
      DFAC=DFAC*EIGHT
      IF( ABS(DFAC).GT.DABMAX) DABMAX= ABS(DFAC)
      NN=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
      DAB(NN)=DFAC*PIJKL
      IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,I,J,K,L,NN,DFAC
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
      RETURN
 9999 FORMAT(' -DABROH- ',4I4,4I3,I5,F20.12)
      END
C
      SUBROUTINE HND_DABGVB(II,JJ,KK,LL,NORB,NOCORE,NOOPEN,NOPAIR,
     1 IA,DC,V,NCONF,NDIM,DAB,DABMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL NOCORE,NOOPEN,NOPAIR
      LOGICAL OUT
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1LGXYZ(4,35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      DIMENSION IA(1),DC(1),V(NDIM,1),NCONF(1),DAB(1)
      DATA ZERO,PT5,EIGHT /0.0D+00,0.5D+00,8.0D+00/
C
      MINI=KMIN(II)
      MINJ=KMIN(JJ)
      MINK=KMIN(KK)
      MINL=KMIN(LL)
      MAXI=KMAX(II)
      MAXJ=KMAX(JJ)
      MAXK=KMAX(KK)
      MAXL=KMAX(LL)
      LOCI=KLOC(II)-MINI
      LOCJ=KLOC(JJ)-MINJ
      LOCK=KLOC(KK)-MINK
      LOCL=KLOC(LL)-MINL
C
      DABMAX=ZERO
      DO 500 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 400 J=MINJ,JMAX
      PIJ=PI*PNRM(J)
       MAX=MAXK
      IF(IJEQKL)  MAX=I
      DO 300 K=MINK, MAX
      PIJK=PIJ*PNRM(K)
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 200 L=MINL,LMAX
      PIJKL=PIJK*PNRM(L)
      I1=LOCI+I
      I2=LOCJ+J
      I3=LOCK+K
      I4=LOCL+L
      IF(I1.GE.I2) GO TO 10
      N=I1
      I1=I2
      I2=N
   10 IF(I3.GE.I4) GO TO 20
      N=I3
      I3=I4
      I4=N
   20 IF(I1-I3) 30,40,50
   30 N=I1
      I1=I3
      I3=N
      N=I2
      I2=I4
      I4=N
      GO TO 50
   40 IF(I2.LT.I4) GO TO 30
   50 CONTINUE
      MIJ=IA(I1)+I2
      MIK=IA(I1)+I3
      MIL=IA(I1)+I4
      MKL=IA(I3)+I4
      IF(I2.LT.I3) GO TO 60
      MJK=IA(I2)+I3
      MJL=IA(I2)+I4
      GO TO 80
   60 MJK=IA(I3)+I2
      IF(I2.LT.I4) GO TO 70
      MJL=IA(I2)+I4
      GO TO 80
   70 MJL=IA(I4)+I2
   80 CONTINUE
      DFAC = ZERO
      IF(NOCORE) GO TO 100
      DFAC = DFAC + ALPHA(1)* DC(MIJ)*DC(MKL) +
     1          PT5* BETA(1)*(DC(MIK)*DC(MJL) + DC(MIL)*DC(MJK))
      IF(NOOPEN.AND.NOPAIR) GO TO 120
      NCO1 = NCO + 1
      DO 90 IO = NCO1,NORB
      IOJO = IA(NCONF(IO)) + 1
      DFAC = DFAC + ALPHA(IOJO)*(DC(MIJ)*V(I3,IO)*V(I4,IO)
     1                          +DC(MKL)*V(I1,IO)*V(I2,IO)) +
     1 PT5*BETA(IOJO)*(V(I2,IO)*(DC(MIK)*V(I4,IO)
     1                          +DC(MIL)*V(I3,IO)) +
     1                 V(I1,IO)*(DC(MJL)*V(I3,IO)
     1                          +DC(MJK)*V(I4,IO)))
   90 CONTINUE
  100 CONTINUE
      IF(NOOPEN.AND.NOPAIR) GO TO 120
      NCO1 = NCO + 1
      DO 110 IO = NCO1,NORB
      DO 110 JO = NCO1,NORB
      IOF = NCONF(IO)
      JOF = NCONF(JO)
      IOJO = IA(IOF) + JOF
      IF(JOF.GT.IOF) IOJO = IA(JOF) + IOF
      DFAC = DFAC +
     1       ALPHA(IOJO)*V(I1,IO)*V(I2,IO)*V(I3,JO)*V(I4,JO) +
     1   PT5* BETA(IOJO)*V(I1,IO)*V(I2,JO)*
     1                  (V(I3,IO)*V(I4,JO)+V(I4,IO)*V(I3,JO))
  110 CONTINUE
  120 CONTINUE
C
      IF(I2.EQ.I1             ) DFAC=DFAC*PT5
      IF(I4.EQ.I3             ) DFAC=DFAC*PT5
      IF(I3.EQ.I1.AND.I4.EQ.I2) DFAC=DFAC*PT5
      DFAC=DFAC*EIGHT
      IF( ABS(DFAC).GT.DABMAX) DABMAX= ABS(DFAC)
      NN=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
      DAB(NN)=DFAC*PIJKL
      IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,I,J,K,L,NN,DFAC
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
      RETURN
 9999 FORMAT(' -DABGVB- ',4I4,4I3,I5,F20.12)
      END
C
      SUBROUTINE HND_DABMC(II,JJ,KK,LL,DC,DV,DM2,IA,
     1                 NCORBS,NCI,LDM2,DAB,DABMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL OUT
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1LGXYZ(4,35)
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_DSHPRT/OUT
      DIMENSION IA(1),DC(1),DV(1)
      DIMENSION DM2(LDM2,LDM2,LDM2,LDM2)
      DIMENSION DAB(1)
      DATA ZERO,PT5,TWO,FOUR /0.0D+00,0.5D+00,2.0D+00,4.0D+00/
      DATA PT125,PT25,ONE    /0.125D+00,0.25D+00,1.0D+00/
C
      MINI=KMIN(II)
      MINJ=KMIN(JJ)
      MINK=KMIN(KK)
      MINL=KMIN(LL)
      MAXI=KMAX(II)
      MAXJ=KMAX(JJ)
      MAXK=KMAX(KK)
      MAXL=KMAX(LL)
      LOCI=KLOC(II)-KMIN(II)
      LOCJ=KLOC(JJ)-KMIN(JJ)
      LOCK=KLOC(KK)-KMIN(KK)
      LOCL=KLOC(LL)-KMIN(LL)
C
      DABMAX=ZERO
      DO 500 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 400 J=MINJ,JMAX
      PIJ=PI*PNRM(J)
       MAX=MAXK
      IF(IJEQKL)  MAX=I
      DO 300 K=MINK, MAX
      PIJK=PIJ*PNRM(K)
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 200 L=MINL,LMAX
      PIJKL=PIJK*PNRM(L)
C
      I1=LOCI+I
      I2=LOCJ+J
      I3=LOCK+K
      I4=LOCL+L
       IF(I1.LT.I2) THEN
         N =I1
         I1=I2
         I2=N
       ENDIF
       IF(I3.LT.I4) THEN
         N =I3
         I3=I4
         I4=N
       ENDIF
       IF(I1-I3) 10,20,30
   10    N =I1
         I1=I3
         I3=N
         N =I2
         I2=I4
         I4=N
      GO TO 30
   20  IF(I2.LT.I4) GO TO 10
   30 CONTINUE
      MIJ=IA(I1)+I2
      MIK=IA(I1)+I3
      MIL=IA(I1)+I4
      MKL=IA(I3)+I4
       IF(I2.LT.I3) GO TO 40
         MJK=IA(I2)+I3
         MJL=IA(I2)+I4
       GO TO 60
   40 MJK=IA(I3)+I2
       IF(I2.LT.I4) GO TO 50
         MJL=IA(I2)+I4
       GO TO 60
   50 MJL=IA(I4)+I2
   60 CONTINUE
      DFAC = ZERO
C
C      -----  CONTRIBUTION CORE - CORE  ---------
C
      IF(NCORBS.EQ.0) GO TO 70
C
      DFAC = DFAC +       DC(MIJ)*DC(MKL)
     1            - PT25*(DC(MIK)*DC(MJL) + DC(MIL)*DC(MJK))
C
      IF(NCI.EQ.0) GO TO 80
C
C      -----  CONTRIBUTION CORE - VALENCE ---------
C
      DFAC = DFAC +       DC(MIJ)*DV(MKL)
     1            +       DC(MKL)*DV(MIJ)
     2            - PT25*(DC(MIK)*DV(MJL)
     3                   +DC(MIL)*DV(MJK)
     4                   +DC(MJK)*DV(MIL)
     5                   +DC(MJL)*DV(MIK))
C
C      -----  CONTRIBUTION VALENCE - VALENCE ---------
C
   70 CONTINUE
C
      DFAC= DFAC + DM2(I,J,K,L)
C
   80 CONTINUE
      IF(I2.EQ.I1             ) DFAC=DFAC*PT5
      IF(I4.EQ.I3             ) DFAC=DFAC*PT5
      IF(I3.EQ.I1.AND.I4.EQ.I2) DFAC=DFAC*PT5
      DFAC=DFAC*FOUR
      IF( ABS(DFAC).GT.DABMAX) DABMAX= ABS(DFAC)
      NN=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
      DAB(NN)=DFAC*PIJKL
      IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,I,J,K,L,NN,DFAC,DM2(I,J,K,L)
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
      RETURN
 9999 FORMAT(' -DABMC - ',4I4,4I3,I5,E20.12,E20.12)
      END
C
C
      SUBROUTINE HND_DABMCL(V,GIN,GOUT,NCI,II,JJ,KK,LL,NUM,LDM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      DIMENSION V(NUM,1)
      DIMENSION GIN(LDM2,LDM2,LDM2,NCI)
      DIMENSION GOUT(LDM2,LDM2,LDM2,LDM2)
      DATA ZERO /0.0D+00/
C
C     ----- FINAL TRANSFORMATION OF G(M,N,P,L) IN G(M,N,P,Q) -----
C
      LX = KLOC(LL)-KMIN(LL)
      DO 30 NP =KMIN(KK),KMAX(KK)
      DO 30 N  =KMIN(JJ),KMAX(JJ)
      DO 30 M  =KMIN(II),KMAX(II)
        DO 20 NQ =KMIN(LL),KMAX(LL)
        LN = LX + NQ
        DUM  = ZERO
        DO 10 L = 1, NCI
   10     DUM = DUM + GIN(M,N,NP,L) * V(LN,L)
   20   GOUT(M,N,NP,NQ)=DUM
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_DABCI(II,JJ,KK,LL,NBUF,D,ID,DAB,DABMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      INTEGER ZSHIFTR
      LOGICAL OUT
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_INTFIL/NINTMX
      COMMON/NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1LGXYZ(4,35)
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_CIDAB/NDMAX,ND
      DIMENSION D(1),ID(1),DAB(1)
#include "bitops_decls.fh"
      DATA ZERO /0.0D+00/
      DATA NFTI /16/
C
#include "bitops_funcs.fh"
      ZSHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      MASK1(IARG)=2**IARG-1
      IWORD1(IARG)=LAND(ZSHIFTR(IARG,24),MASK1(8))
      IWORD2(IARG)=LAND(ZSHIFTR(IARG,16),MASK1(8))
      IWORD3(IARG)=LAND(ZSHIFTR(IARG, 8),MASK1(8))
      IWORD4(IARG)=LAND(ZSHIFTR(IARG, 0),MASK1(8))
C
      DABMAX=ZERO
      IF(ND.GT.0) GO TO 200
C
C     ----- READ IN BUFFERS OF -2DM- -----
C
  100 CONTINUE
      NDMAX=0
      IADRSS=0
      DO 120 IBUF=1,NBUF
      CALL HND_PREAD(NFTI,D(IADRSS+1),ID(IADRSS+1),MD,NINTMX)
      NDMAX=NDMAX+ ABS(MD)
      IF(MD.LE.0) GO TO 130
  120 IADRSS=IADRSS+NINTMX
  130 CONTINUE
      ND=0
C
C     ----- SET UP -ABDENS- ARRAY -----
C
  200 CONTINUE
      ND=ND+1
      IF(ND.GT.NDMAX) GO TO 100
      LABEL=ID(ND)
        VAL= D(ND)
      IF(LABEL.EQ.0) GO TO 300
      I=IWORD1(LABEL)
      J=IWORD2(LABEL)
      K=IWORD3(LABEL)
      L=IWORD4(LABEL)
      IF( ABS(VAL).GT.DABMAX) DABMAX= ABS(VAL)
      PIJKL=PNRM(I)*PNRM(J)*PNRM(K)*PNRM(L)
      NABD=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
      DAB(NABD)=VAL*PIJKL
      GO TO 200
C
  300 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JKDSPD(GINT,FINT,SINT,IJKLG,
     1 GIJKL,GNKL,GNM,XY,XZ,YZ,X,Y,Z,SJ,SK,SL,
     2 B00,B01,B10,C00,D00,F00,
     3 DIJ,DKL,DIJSI,DIJSJ,DKLSK,DKLSL,
     4 ABV,CV,RWV,AAI,AAJ,BBK,BBL,
     5 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     6 DAB,MAXXYZ,FC,NC,DF,LDF,NBF,DEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXGSH =30,MXGSH2=MXGSH*MXGSH)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      LOGICAL EXPNDI,EXPNDK
      LOGICAL FIRST
      LOGICAL LAST
      LOGICAL OUT
      COMMON/HND_JKSTAT/LVEC,NVEC
      COMMON/HND_IJDCHR/ADCHRG(13,MXGSH2),DASI(MXGSH2),DASJ(MXGSH2),
     1NIJ
      COMMON/HND_KLDCHR/BDCHRG(13,MXGSH2),DBSK(MXGSH2),DBSL(MXGSH2),
     1NKL
      COMMON/HND_SHLFAC/QQ4
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLXPN/EXPNDI,EXPNDK
      COMMON/HND_SHLTOL/RTOL,DTOL
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DSHDAB/DABTOL,DABMAX
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_HFK/XX,U(9),W(9),NROOTS
      COMMON/HND_RYS/YY,T(5),V(5),MROOTS
      COMMON/HND_INTVEC/MAXVEC,MINVEC
      DIMENSION GINT(1),FINT(12,1),SINT(78,1),IJKLG(4,1)
      DIMENSION XY(1),XZ(1),YZ(1),X(1),Y(1),Z(1),SJ(1),SK(1),SL(1)
      DIMENSION GIJKL(1),GNKL(1),GNM(1)
      DIMENSION DIJ(1),DKL(1),DIJSI(1),DIJSJ(1),DKLSK(1),DKLSL(1)
      DIMENSION FI(1),FJ(1),FK(1),FL(1)
      DIMENSION SII(1),SJJ(1),SKK(1),SLL(1)
      DIMENSION SIJ(1),SIK(1),SIL(1),SJK(1),SJL(1),SKL(1)
      DIMENSION B00(1),B01(1),B10(1),C00(1),D00(1),F00(1)
      DIMENSION AAI(1),AAJ(1),BBK(1),BBL(1)
      DIMENSION ABV(5,1),CV(18,1),RWV(2,1)
      DIMENSION DAB(1)
      DIMENSION FC(NC,1),DF(LDF,1),DEN(1)
      DATA PI252  /34.986836655250D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TENM01 /1.0D-01/
C
      Q4=PI252*QQ4
      DABCUT=DABTOL*TENM01
C
      NIMAX=LIT + IDER
      NJMAX=LJT + JDER
      NKMAX=LKT + KDER
      NLMAX=LLT + LDER
      NMAX=LIT+LJT-1 + MIN0(IDER+JDER,NDER)
      MMAX=LKT+LLT-1 + MIN0(KDER+LDER,NDER)
      NMAXS=NMAX.EQ.1
      NMAXP=NMAX.LE.2
      MMAXS=MMAX.EQ.1
      MMAXP=MMAX.LE.2
C
      MAXG=MAXXYZ/NROOTS
      DTOL2=DTOL*DTOL
C
C     ----- PAIR OF K,L PRIMITIVES -----
C
      FIRST=.TRUE.
      NG=0
      KLG=0
  100 KLG=KLG+1
      IF(KLG.GT.NKL) GO TO 300
      DB=BDCHRG( 1,KLG)
      BB=BDCHRG( 2,KLG)
      XB=BDCHRG( 3,KLG)
      YB=BDCHRG( 4,KLG)
      ZB=BDCHRG( 5,KLG)
      XD=BDCHRG( 6,KLG)
      YD=BDCHRG( 7,KLG)
      ZD=BDCHRG( 8,KLG)
      DXKL=BDCHRG( 9,KLG)
      DYKL=BDCHRG(10,KLG)
      DZKL=BDCHRG(11,KLG)
      Q4DB=Q4*DB
C
C     ----- PAIR OF I,J PRIMITIVES -----
C
      IJG=0
  200 IJG=IJG+1
      IF(IJG.GT.NIJ) GO TO 100
      DA=ADCHRG( 1,IJG)
      AA=ADCHRG( 2,IJG)
      XA=ADCHRG( 3,IJG)
      YA=ADCHRG( 4,IJG)
      ZA=ADCHRG( 5,IJG)
      AANDB1=ONE/(AA+BB)
      Q4DBDA=Q4DB*DA
      DUM   =Q4DBDA*Q4DBDA*AANDB1
      IF(DUM.LE.DTOL2) GO TO 200
      Q4DBDA=Q4DBDA* SQRT(AANDB1)
      IF( ABS(Q4DBDA*DABMAX).LT.DABCUT) GO TO 200
      RHO   =AA*BB*AANDB1
      XX=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
C
      NG=NG+1
      ABV(1,NG)=AA
      ABV(2,NG)=BB
      ABV(3,NG)=RHO
      ABV(4,NG)=Q4DBDA
      ABV(5,NG)=XX
C
      XC=ADCHRG( 6,IJG)
      YC=ADCHRG( 7,IJG)
      ZC=ADCHRG( 8,IJG)
      DXIJ=ADCHRG( 9,IJG)
      DYIJ=ADCHRG(10,IJG)
      DZIJ=ADCHRG(11,IJG)
C
      AAI(NG)=ADCHRG(12,IJG)
      AAJ(NG)=ADCHRG(13,IJG)
      BBK(NG)=BDCHRG(12,KLG)
      BBL(NG)=BDCHRG(13,KLG)
C
      IF(MMAXS) GO TO 210
      CV( 1,NG)=AA*(XA-XD)
      CV( 2,NG)=BB*(XB-XD)
      CV( 3,NG)=AA*(YA-YD)
      CV( 4,NG)=BB*(YB-YD)
      CV( 5,NG)=AA*(ZA-ZD)
      CV( 6,NG)=BB*(ZB-ZD)
  210 IF(NMAXS) GO TO 220
      CV( 7,NG)=AA*(XA-XC)
      CV( 8,NG)=BB*(XB-XC)
      CV( 9,NG)=AA*(YA-YC)
      CV(10,NG)=BB*(YB-YC)
      CV(11,NG)=AA*(ZA-ZC)
      CV(12,NG)=BB*(ZB-ZC)
  220 CONTINUE
      CV(13,NG)=DXIJ
      CV(14,NG)=DYIJ
      CV(15,NG)=DZIJ
      CV(16,NG)=DXKL
      CV(17,NG)=DYKL
      CV(18,NG)=DZKL
      IF(SPI) DIJSI(NG)=DASI(IJG)
      IF(SPJ) DIJSJ(NG)=DASJ(IJG)
      IF(SPK) DKLSK(NG)=DBSK(KLG)
      IF(SPL) DKLSL(NG)=DBSL(KLG)
C
      IF(NG.LT.MAXG) GO TO 200
      LAST=.FALSE.
      GO TO 310
C
  300 CONTINUE
      LAST=.TRUE.
  310 CONTINUE
      NUMG=NG
      IF(NUMG.EQ.0) GO TO 1000
C
      IF(NROOTS.EQ.1) GO TO 480
      IF(.NOT.SPI) GO TO 420
      DO 410 IROOT=2,NROOTS
      DO 410 IG=1,NUMG
      DIJSI(IG+NUMG*(IROOT-1))=DIJSI(IG)
  410 CONTINUE
  420 IF(.NOT.SPJ) GO TO 440
      DO 430 IROOT=2,NROOTS
      DO 430 IG=1,NUMG
      DIJSJ(IG+NUMG*(IROOT-1))=DIJSJ(IG)
  430 CONTINUE
  440 IF(.NOT.SPK) GO TO 460
      DO 450 IROOT=2,NROOTS
      DO 450 IG=1,NUMG
      DKLSK(IG+NUMG*(IROOT-1))=DKLSK(IG)
  450 CONTINUE
  460 IF(.NOT.SPL) GO TO 480
      DO 470 IROOT=2,NROOTS
      DO 470 IG=1,NUMG
      DKLSL(IG+NUMG*(IROOT-1))=DKLSL(IG)
  470 CONTINUE
  480 CONTINUE
C
      IF(SKIPI) GO TO 520
      DO 510 IRXYZ=2,NROOTS*3
      DO 510 IG=1,NUMG
      AAI(IG+NUMG*(IRXYZ-1))=AAI(IG)
  510 CONTINUE
  520 IF(SKIPJ) GO TO 540
      DO 530 IRXYZ=2,NROOTS*3
      DO 530 IG=1,NUMG
      AAJ(IG+NUMG*(IRXYZ-1))=AAJ(IG)
  530 CONTINUE
  540 IF(SKIPK) GO TO 560
      DO 550 IRXYZ=2,NROOTS*3
      DO 550 IG=1,NUMG
      BBK(IG+NUMG*(IRXYZ-1))=BBK(IG)
  550 CONTINUE
  560 IF(SKIPL) GO TO 580
      DO 570 IRXYZ=2,NROOTS*3
      DO 570 IG=1,NUMG
      BBL(IG+NUMG*(IRXYZ-1))=BBL(IG)
  570 CONTINUE
  580 CONTINUE
C
      LVEC=LVEC+NUMG*NROOTS*3
      NVEC=NVEC+1
C
C     ----- COMPUTE ROOTS AND WEIGHTS FOR QUADRATURE -----
C
      CALL HND_JKWRYS(RWV,ABV,NUMG)
C
C     ----- COMPUTE COEFFICIENTS FOR RECURSION FORMULAE -----
C
      CALL HND_JKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1            ABV,CV,RWV,NUMG,NROOTS)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 2 CENTERS, 2-D ) -----
C
      IF(NUMG*NROOTS*3.LT.MINVEC)
     1 CALL HND_JKGNMS(GNM,NUMG*NROOTS*3,NMAX,MMAX,
     2             B00,B01,B10,C00,D00,F00)
      IF(NUMG*NROOTS*3.GE.MINVEC)
     1 CALL HND_JKGNMV(GNM,NUMG*NROOTS*3,NMAX,MMAX,
     2             B00,B01,B10,C00,D00,F00)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 4 CENTERS, 2-D ) -----
C
      IF(NUMG*NROOTS*3.LT.MINVEC)
     1 CALL HND_JKXYZS(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,GNM,
     2 NUMG*NROOTS*3,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     3 EXPNDI,EXPNDK)
      IF(NUMG*NROOTS*3.GE.MINVEC)
     1 CALL HND_JKXYZV(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,GNM,
     2 NUMG*NROOTS*3,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     3 EXPNDI,EXPNDK)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS FOR DERIVATIVES -----
C
      IF(NUMG*NROOTS*3.LT.MINVEC)
     1 CALL HND_JDXYZS(GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,
     2 NUMG*NROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,LIT,LJT,LKT,LLT,
     3 AAI,AAJ,BBK,BBL,
     4 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
      IF(NUMG*NROOTS*3.GE.MINVEC)
     1 CALL HND_JDXYZV(GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,
     2 NUMG*NROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,LIT,LJT,LKT,LLT,
     3 AAI,AAJ,BBK,BBL,
     4 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
C
C     ----- ZERO OUT FIRST TIME AROUND -----
C
      IF(.NOT.FIRST) GO TO 600
      CALL HND_JKDZER(GINT,FINT,SINT,IJKLG)
      FIRST=.FALSE.
  600 CONTINUE
C
C     ----- COMPUTE DERIVATIVE INTEGRALS -----
C
      IF(NDER.LE.1) THEN
C
      IF(NUMG*NROOTS.LT.MINVEC)
     1 CALL HND_DSPDFS(NUMG,NROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IF(NUMG*NROOTS.GE.MINVEC)
     1 CALL HND_DSPDFV(NUMG,NROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
C
      ELSE
C
      IF(NUMG*NROOTS.LT.MINVEC)
     1 CALL HND_DDSPDS(NUMG,NROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IF(NUMG*NROOTS.GE.MINVEC)
     1 CALL HND_DDSPDV(NUMG,NROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
C
      ENDIF
      IF(LAST) GO TO 1000
      NG=0
      GO TO 200
 1000 CONTINUE
      IF(NUMG.EQ.0.AND.FIRST) RETURN
C
C     ----- PROCESS DERIVATIVE INTEGRALS -----
C
      IF(NDER.LE.1) THEN
      CALL HND_JKDINV(FINT,SINT,IJKLG,FC,NC,DF,LDF,NBF,DEN)
      ELSE
      CALL HND_JDDINV(FINT,SINT,IJKLG,FC,NC,DF,LDF,NBF,DEN)
      ENDIF
C
      IF(OUT) CALL HND_JKDOUT(GINT,FINT,SINT,IJKLG,DAB)
      RETURN
      END
C
C
C     ------ SUBROUTINES ABOVE TAKEN FROM DER.F IN HONDO ------
C
C
C     ------ SUBROUTINES TAKEN FROM INT.F IN HONDO ------
C
      SUBROUTINE HND_JKWRYS(RWV,ABV,NUMG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_HFK/XX,U(9),W(9),NROOTS
      COMMON/HND_RYS/YY,T(5),V(5),MROOTS
      DIMENSION RWV(2,NUMG,1),ABV(5,1)
C
      IF(NROOTS.GT.5) GO TO 100
C
C     ----- NROOTS .LE. 5 -----
C
      DO 20 NG=1,NUMG
      YY=ABV(5,NG)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 10 NR=1,NROOTS
      RWV(1,NG,NR)=T(NR)
      RWV(2,NG,NR)=V(NR)
   10 CONTINUE
   20 CONTINUE
      RETURN
C
  100 CONTINUE
C
C     ----- NROOTS .GT. 5 -----
C
      DO 120 NG=1,NUMG
      XX=ABV(5,NG)
      CALL HND_DROOT
      DO 110 NR=1,NROOTS
      RWV(1,NG,NR)=U(NR)
      RWV(2,NG,NR)=W(NR)
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1                  ABV,CV,RWV,NUMG,NROOTS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION B00(NUMG,NROOTS,3),B01(NUMG,NROOTS,3),B10(NUMG,NROOTS,3)
      DIMENSION C00(NUMG,NROOTS,3)
      DIMENSION D00(NUMG,NROOTS,3)
      DIMENSION F00(NUMG,NROOTS,3),DIJ(NUMG,NROOTS,3),DKL(NUMG,NROOTS,3)
      DIMENSION ABV(5,1),CV(18,1)
      DIMENSION RWV(2,NUMG,NROOTS)
      DATA PT5,ONE /0.5D+00,1.0D+00/
C
      DO 40 NR=1,NROOTS
      DO 30 NG=1,NUMG
      AA =ABV(1,NG)
      BB =ABV(2,NG)
      RHO=ABV(3,NG)
      QAB=ABV(4,NG)
      UU =RHO*RWV(1,NG,NR)
      WW =    RWV(2,NG,NR)
      AAUU=AA+UU
      BBUU=BB+UU
      F00(NG,NR,1)=WW*QAB
      F00(NG,NR,2)=ONE
      F00(NG,NR,3)=ONE
      DUM2=PT5/(AA*BB+UU*(AA+BB))
      AUDUM=AAUU*DUM2
      BUDUM=BBUU*DUM2
       UDUM=  UU*DUM2
      B00(NG,NR,1)= UDUM
      B00(NG,NR,2)= UDUM
      B00(NG,NR,3)= UDUM
      B01(NG,NR,1)=AUDUM
      B01(NG,NR,2)=AUDUM
      B01(NG,NR,3)=AUDUM
      B10(NG,NR,1)=BUDUM
      B10(NG,NR,2)=BUDUM
      B10(NG,NR,3)=BUDUM
       UDUM= UDUM+ UDUM
      IF(MMAXS) GO TO 10
      AUDUM=AUDUM+AUDUM
      D00(NG,NR,1)= UDUM*CV( 1,NG) + AUDUM*CV( 2,NG)
      D00(NG,NR,2)= UDUM*CV( 3,NG) + AUDUM*CV( 4,NG)
      D00(NG,NR,3)= UDUM*CV( 5,NG) + AUDUM*CV( 6,NG)
   10 IF(NMAXS) GO TO 20
      BUDUM=BUDUM+BUDUM
      C00(NG,NR,1)= UDUM*CV( 8,NG) + BUDUM*CV( 7,NG)
      C00(NG,NR,2)= UDUM*CV(10,NG) + BUDUM*CV( 9,NG)
      C00(NG,NR,3)= UDUM*CV(12,NG) + BUDUM*CV(11,NG)
   20 CONTINUE
C
   30 CONTINUE
   40 CONTINUE
C
      DO 60 NR=1,NROOTS
      DO 50 NG=1,NUMG
      DIJ(NG,NR,1)=CV(13,NG)
      DIJ(NG,NR,2)=CV(14,NG)
      DIJ(NG,NR,3)=CV(15,NG)
      DKL(NG,NR,1)=CV(16,NG)
      DKL(NG,NR,2)=CV(17,NG)
      DKL(NG,NR,3)=CV(18,NG)
   50 CONTINUE
   60 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKGNMS(GNM,NG,NMAX,MMAX,
     1 B00,B01,B10,C00,D00,F00)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION GNM(NG,NMAX,MMAX)
      DIMENSION C00(NG),D00(NG),F00(NG)
      DIMENSION B00(NG,1),B01(NG,1),B10(NG,1)
      DATA ZERO,ONE /0.0D+00,1.0D+00/
C
C     ----- G(0,0) -----
C
      DO 10 IG=1,NG
      GNM(IG,1,1)=F00(IG)
   10 CONTINUE
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
C
C     ----- G(1,0) = C00 * G(0,0) -----
C
      DO 20 IG=1,NG
      GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
   20 CONTINUE
C
   30 CONTINUE
      IF(MMAXS) GO TO 60
C
C     ----- G(0,1) = D00 * G(0,0) -----
C
      DO 40 IG=1,NG
      GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
   40 CONTINUE
      IF(NMAXS) GO TO 60
C
C     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
C
      DO 50 IG=1,NG
      GNM(IG,2,2)=B00(IG,1)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
   50 CONTINUE
C
   60 CONTINUE
      MAX=MAX0(NMAX-1,MMAX-1)
      DO 70 M=2,MAX
      DO 70 IG=1,NG
   70 B00(IG,M)=B00(IG,M-1)+B00(IG,1)
C
      IF(NMAXP) GO TO 120
C
C     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
C
      DO 80 N=2,NMAX-1
      DO 80 IG=1,NG
      B10(IG,N)=B10(IG,N-1)+B10(IG,1)
   80 CONTINUE
      DO 90 N=2,NMAX-1
      DO 90 IG=1,NG
      GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
   90 CONTINUE
      IF(MMAXS) GO TO 110
C
C     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
C
      DO 100 N=2,NMAX-1
      DO 100 IG=1,NG
      GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
  100 CONTINUE
C
  110 CONTINUE
C
  120 CONTINUE
      IF(MMAXP) GO TO 170
C
C     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
C
      DO 130 M=2,MMAX-1
      DO 130 IG=1,NG
      B01(IG,M)=B01(IG,M-1)+B01(IG,1)
  130 CONTINUE
      DO 140 M=2,MMAX-1
      DO 140 IG=1,NG
      GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
  140 CONTINUE
      IF(NMAXS) GO TO 160
C
C     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
C
      DO 150 M=2,MMAX-1
      DO 150 IG=1,NG
      GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
  150 CONTINUE
C
  160 CONTINUE
  170 IF(NMAXP.OR.MMAXP) RETURN
C
C     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
C                    +     C00 * G(N  ,M  )
C                    + M * B00 * G(N  ,M-1)
C
      DO 180 M=2,MMAX-1
      DO 180 N=2,NMAX-1
      DO 180 IG=1,NG
      GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+
     1                C00(IG    )*GNM(IG,N  ,M+1)+
     2                B00(IG,M  )*GNM(IG,N  ,M  )
  180 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKGNMV(GNM,NG,NMAX,MMAX,
     1 B00,B01,B10,C00,D00,F00)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION GNM(NG,NMAX,MMAX)
      DIMENSION C00(NG),D00(NG),F00(NG)
      DIMENSION B00(NG,1),B01(NG,1),B10(NG,1)
      DATA ZERO,ONE /0.0D+00,1.0D+00/
C
C     ----- G(0,0) -----
C
      DO 10 IG=1,NG
      GNM(IG,1,1)=F00(IG)
   10 CONTINUE
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
C
C     ----- G(1,0) = C00 * G(0,0) -----
C
      DO 20 IG=1,NG
      GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
   20 CONTINUE
C
   30 CONTINUE
      IF(MMAXS) GO TO 60
C
C     ----- G(0,1) = D00 * G(0,0) -----
C
      DO 40 IG=1,NG
      GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
   40 CONTINUE
      IF(NMAXS) GO TO 60
C
C     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
C
      DO 50 IG=1,NG
      GNM(IG,2,2)=B00(IG,1)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
   50 CONTINUE
C
   60 CONTINUE
      MAX=MAX0(NMAX-1,MMAX-1)
      DO 70 IG=1,NG
      DO 70 M=2,MAX
   70 B00(IG,M)=B00(IG,M-1)+B00(IG,1)
C
      IF(NMAXP) GO TO 120
C
C     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
C
      DO 80 IG=1,NG
      DO 80 N=2,NMAX-1
      B10(IG,N)=B10(IG,N-1)+B10(IG,1)
   80 CONTINUE
      DO 90 IG=1,NG
      DO 90 N=2,NMAX-1
      GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
   90 CONTINUE
      IF(MMAXS) GO TO 110
C
C     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
C
      DO 100 IG=1,NG
      DO 100 N=2,NMAX-1
      GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
  100 CONTINUE
C
  110 CONTINUE
C
  120 CONTINUE
      IF(MMAXP) GO TO 170
C
C     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
C
      DO 130 IG=1,NG
      DO 130 M=2,MMAX-1
      B01(IG,M)=B01(IG,M-1)+B01(IG,1)
  130 CONTINUE
      DO 140 IG=1,NG
      DO 140 M=2,MMAX-1
      GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
  140 CONTINUE
      IF(NMAXS) GO TO 160
C
C     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
C
      DO 150 IG=1,NG
      DO 150 M=2,MMAX-1
      GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
  150 CONTINUE
C
  160 CONTINUE
  170 IF(NMAXP.OR.MMAXP) RETURN
C
C     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
C                    +     C00 * G(N  ,M  )
C                    + M * B00 * G(N  ,M-1)
C
      DO 180 IG=1,NG
      DO 180 N=2,NMAX-1
      DO 180 M=2,MMAX-1
      GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+
     1                C00(IG    )*GNM(IG,N  ,M+1)+
     2                B00(IG,M  )*GNM(IG,N  ,M  )
  180 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKXYZS(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,HNM,
     1 NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     2 EXPNDI,EXPNDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EXPNDI,EXPNDK
      DIMENSION GIJKL(NG*NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION HIJKL(NG*NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION  GNKL(NG,NLMAX,NKMAX,NMAX)
      DIMENSION  HNKL(NG*NLMAX*NKMAX,NMAX)
      DIMENSION  FNKL(NG*NLMAX*NKMAX*NMAX)
      DIMENSION   GNM(NG,NMAX,MMAX)
      DIMENSION   DIJ(NG)
      DIMENSION   DKL(NG)
C
C     ----- G(N,K,L) -----
C
      IF(EXPNDK) GO TO 40
C
      DO 30 NK=1,NKMAX
      DO 10 NL=1,NLMAX
      DO 10  N=1,NMAX
      DO 10 IG=1,NG
   10 GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
      IF(NK.EQ.NKMAX) GO TO 30
      MAX=MMAX-NK
      DO 20  M=1,MAX
      DO 20  N=1,NMAX
      DO 20 IG=1,NG
   20 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   30 CONTINUE
C
      GO TO 100
   40 CONTINUE
C
      DO 70 NL=1,NLMAX
      DO 50 NK=1,NKMAX
      DO 50  N=1,NMAX
      DO 50 IG=1,NG
   50 GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
      IF(NL.EQ.NLMAX) GO TO 70
      MAX=MMAX-NL
      DO 60  M=1,MAX
      DO 60  N=1,NMAX
      DO 60 IG=1,NG
   60 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   70 CONTINUE
C
  100 CONTINUE
C
C     ----- G(I,J,K,L) -----
C
      IF(EXPNDI) GO TO 140
C
      DO 130 NI=1,NIMAX
      DO 110 IGLKJ=1,NG*NLMAX*NKMAX*NJMAX
  110 HIJKL(IGLKJ,NI)=FNKL(IGLKJ)
      IF(NI.EQ.NIMAX) GO TO 130
      MAX=NMAX-NI
      DO 120  N=1,MAX
      DO 120 NK=1,NKMAX
      DO 120 NL=1,NLMAX
      DO 120 IG=1,NG
  120 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  130 CONTINUE
C
      RETURN
  140 CONTINUE
C
      DO 170 NJ=1,NJMAX
      DO 150 NI=1,NIMAX
      DO 150 IGLK=1,NG*NLMAX*NKMAX
  150 GIJKL(IGLK,NJ,NI)=HNKL(IGLK,NI)
      IF(NJ.EQ.NJMAX) GO TO 170
      MAX=NMAX-NJ
      DO 160  N=1,MAX
      DO 160 NK=1,NKMAX
      DO 160 NL=1,NLMAX
      DO 160 IG=1,NG
  160 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  170 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKXYZV(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,HNM,
     1 NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     2 EXPNDI,EXPNDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EXPNDI,EXPNDK
      DIMENSION GIJKL(NG*NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION HIJKL(NG*NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION  GNKL(NG,NLMAX,NKMAX,NMAX)
      DIMENSION  HNKL(NG*NLMAX*NKMAX,NMAX)
      DIMENSION  FNKL(NG*NLMAX*NKMAX*NMAX)
      DIMENSION   GNM(NG,NMAX,MMAX)
      DIMENSION   DIJ(NG)
      DIMENSION   DKL(NG)
C
C     ----- G(N,K,L) -----
C
      IF(EXPNDK) GO TO 40
C
      DO 30 NK=1,NKMAX
      DO 10 IG=1,NG
      DO 10 NL=1,NLMAX
      DO 10  N=1,NMAX
   10 GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
      IF(NK.EQ.NKMAX) GO TO 30
      MAX=MMAX-NK
      DO 20 IG=1,NG
      DO 20  M=1,MAX
      DO 20  N=1,NMAX
   20 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   30 CONTINUE
C
      GO TO 100
   40 CONTINUE
C
      DO 70 NL=1,NLMAX
      DO 50 IG=1,NG
      DO 50 NK=1,NKMAX
      DO 50  N=1,NMAX
   50 GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
      IF(NL.EQ.NLMAX) GO TO 70
      MAX=MMAX-NL
      DO 60 IG=1,NG
      DO 60  N=1,NMAX
      DO 60  M=1,MAX
   60 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   70 CONTINUE
C
  100 CONTINUE
C
C     ----- G(I,J,K,L) -----
C
      IF(EXPNDI) GO TO 140
C
      DO 130 NI=1,NIMAX
      DO 110 IGLKJ=1,NG*NLMAX*NKMAX*NJMAX
  110 HIJKL(IGLKJ,NI)=FNKL(IGLKJ)
      IF(NI.EQ.NIMAX) GO TO 130
      MAX=NMAX-NI
      DO 120 IG=1,NG
      DO 120 NL=1,NLMAX
      DO 120 NK=1,NKMAX
      DO 120  N=1,MAX
  120 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  130 CONTINUE
C
      RETURN
  140 CONTINUE
C
      DO 170 NJ=1,NJMAX
      DO 150 IGLK=1,NG*NLMAX*NKMAX
      DO 150 NI=1,NIMAX
  150 GIJKL(IGLK,NJ,NI)=HNKL(IGLK,NI)
      IF(NJ.EQ.NJMAX) GO TO 170
      MAX=NMAX-NJ
      DO 160 IG=1,NG
      DO 160 NL=1,NLMAX
      DO 160 NK=1,NKMAX
      DO 160  N=1,MAX
  160 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  170 CONTINUE
C
      RETURN
      END
C

C
C     ------- SUBROUTINES TAKEN FROM DER.F IN HONDO ------
C

      SUBROUTINE HND_JDXYZS(GI,GIJ,GIJK,GIJKL,HIJK,HIJKL,FIJK,
     1 FIJKL,
     1 NG,NIMAX,NJMAX,NKMAX,NLMAX,NI,NJ,NK,NL,AAI,AAJ,AAK,AAL,
     2 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IS,JS,KS,LS,IP,JP,KP,LP
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,
     1 SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      DIMENSION GIJKL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION  GIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   GIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    GI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION HIJKL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION  HIJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION FIJKL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION  FIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   AAI(NG)
      DIMENSION   AAJ(NG)
      DIMENSION   AAK(NG)
      DIMENSION   AAL(NG)
      DIMENSION    FI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION    FJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    FK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION    FL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SII(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION   SJJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SKK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   SLL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SIK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SIL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION   SJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SJL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SKL(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DATA THREE  /3.0D+00/
C
      IS=NI.EQ.1
      JS=NJ.EQ.1
      KS=NK.EQ.1
      LS=NL.EQ.1
      IP=NI.EQ.2
      JP=NJ.EQ.2
      KP=NK.EQ.2
      LP=NL.EQ.2
C
      IF(NDER.EQ.1) GO TO 1000
C
C     ----- FIRST AND SECOND DERIVATIVES -----
C
      IF(SKIPI) GO TO 40
C
C     ----- -FI- AND -SII- -----
C
      DO  10 LKJ=1,NLMAX*NKMAX*NJMAX
      DO  10 IG =1,NG
       FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
      SII(IG,LKJ,1)=( GI(IG,LKJ,3)*AAI(IG)
     1               -GI(IG,LKJ,1)        )*AAI(IG)
   10 CONTINUE
      IF(IS) GO TO 40
      DO  20 LKJ=1,NLMAX*NKMAX*NJMAX
      DO  20 IG =1,NG
       FI(IG,LKJ,2)= GI(IG,LKJ,3)*AAI(IG)
     1              -GI(IG,LKJ,1)
      SII(IG,LKJ,2)=(GI(IG,LKJ,4)*AAI(IG)
     1              -GI(IG,LKJ,2)*THREE  )*AAI(IG)
   20 CONTINUE
      IF(IP) GO TO 40
      DO  30 I  =3,NI
      DO  30 LKJ=1,NLMAX*NKMAX*NJMAX
      DO  30 IG =1,NG
       FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)
     1              -GI(IG,LKJ,I-1)*DBLE(I-1)
      SII(IG,LKJ,I)=(GI(IG,LKJ,I+2)*AAI(IG)
     1              -GI(IG,LKJ,I  )*DBLE(I+I-1)     )*AAI(IG)
     2              +GI(IG,LKJ,I-2)*DBLE((I-1)*(I-2))
   30 CONTINUE
C
   40 CONTINUE
      IF(SKIPJ) GO TO 140
C
C     ----- -FJ- AND -SJJ- -----
C
      DO 110 I =1,NIMAX
      DO 110 LK=1,NLMAX*NKMAX
      DO 110 IG=1,NG
       FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
      SJJ(IG,LK,1,I)=( GIJ(IG,LK,3,I)*AAJ(IG)
     1                -GIJ(IG,LK,1,I)        )*AAJ(IG)
  110 CONTINUE
      IF(JS) GO TO 140
      DO 120 I =1,NIMAX
      DO 120 LK=1,NLMAX*NKMAX
      DO 120 IG=1,NG
       FJ(IG,LK,2,I)= GIJ(IG,LK,3,I)*AAJ(IG)
     1               -GIJ(IG,LK,1,I)
      SJJ(IG,LK,2,I)=(GIJ(IG,LK,4,I)*AAJ(IG)
     1               -GIJ(IG,LK,2,I)*THREE  )*AAJ(IG)
  120 CONTINUE
      IF(JP) GO TO 140
      DO 130 I =1,NIMAX
      DO 130 J =3,NJ
      DO 130 LK=1,NLMAX*NKMAX
      DO 130 IG=1,NG
       FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG)
     1               -GIJ(IG,LK,J-1,I)*DBLE(J-1)
      SJJ(IG,LK,J,I)=(GIJ(IG,LK,J+2,I)*AAJ(IG)
     1               -GIJ(IG,LK,J  ,I)*DBLE(J+J-1)     )*AAJ(IG)
     2               +GIJ(IG,LK,J-2,I)*DBLE((J-1)*(J-2))
  130 CONTINUE
C
  140 CONTINUE
      IF(SKIPK) GO TO 240
C
C     ----- -FK- AND -SKK- -----
C
      DO 210 JI=1,NJMAX*NIMAX
      DO 210 L =1,NLMAX
      DO 210 IG=1,NG
       FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
      SKK(IG,L,1,JI)=( GIJK(IG,L,3,JI)*AAK(IG)
     1                -GIJK(IG,L,1,JI)        )*AAK(IG)
  210 CONTINUE
      IF(KS) GO TO 240
      DO 220 JI=1,NJMAX*NIMAX
      DO 220 L =1,NLMAX
      DO 220 IG=1,NG
       FK(IG,L,2,JI)= GIJK(IG,L,3,JI)*AAK(IG)
     1               -GIJK(IG,L,1,JI)
      SKK(IG,L,2,JI)=(GIJK(IG,L,4,JI)*AAK(IG)
     1               -GIJK(IG,L,2,JI)*THREE  )*AAK(IG)
  220 CONTINUE
      IF(KP) GO TO 240
      DO 230 JI=1,NJMAX*NIMAX
      DO 230 K =3,NK
      DO 230 L =1,NLMAX
      DO 230 IG=1,NG
       FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)
     1               -GIJK(IG,L,K-1,JI)*DBLE(K-1)
      SKK(IG,L,K,JI)=(GIJK(IG,L,K+2,JI)*AAK(IG)
     1               -GIJK(IG,L,K  ,JI)*DBLE(K+K-1)     )*AAK(IG)
     2               +GIJK(IG,L,K-2,JI)*DBLE((K-1)*(K-2))
  230 CONTINUE
C
  240 CONTINUE
      IF(SKIPL) GO TO 340
C
C     ----- -FL- AND -SLL- -----
C
      DO 310 KJI=1,NKMAX*NJMAX*NIMAX
      DO 310 IG =1,NG
       FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
      SLL(IG,1,KJI)=( GIJKL(IG,3,KJI)*AAL(IG)
     1               -GIJKL(IG,1,KJI)        )*AAL(IG)
  310 CONTINUE
      IF(LS) GO TO 340
      DO 320 KJI=1,NKMAX*NJMAX*NIMAX
      DO 320 IG =1,NG
       FL(IG,2,KJI)= GIJKL(IG,3,KJI)*AAL(IG)
     1              -GIJKL(IG,1,KJI)
      SLL(IG,2,KJI)=(GIJKL(IG,4,KJI)*AAL(IG)
     1              -GIJKL(IG,2,KJI)*THREE  )*AAL(IG)
  320 CONTINUE
      IF(LP) GO TO 340
      DO 330 KJI=1,NKMAX*NJMAX*NIMAX
      DO 330 L  =3,NL
      DO 330 IG =1,NG
       FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)
     1              -GIJKL(IG,L-1,KJI)*DBLE(L-1)
      SLL(IG,L,KJI)=(GIJKL(IG,L+2,KJI)*AAL(IG)
     1              -GIJKL(IG,L  ,KJI)*DBLE(L+L-1)     )*AAL(IG)
     2              +GIJKL(IG,L-2,KJI)*DBLE((L-1)*(L-2))
  330 CONTINUE
C
  340 CONTINUE
      IF(SKIPI) GO TO 660
      IF(SKIPJ) GO TO 460
C
C     ----- -SIJ- -----
C
      DO 410 LK=1,NLMAX*NKMAX
      DO 410 IG=1,NG
  410 SIJ(IG,LK,1,1)=GIJ(IG,LK,2,2)*AAJ(IG)*AAI(IG)
      IF(JS) GO TO 430
      DO 420 J =2,NJ
      DO 420 LK=1,NLMAX*NKMAX
      DO 420 IG=1,NG
      SIJ(IG,LK,J,1)=( GIJ(IG,LK,J+1,2)*AAJ(IG)
     1                -GIJ(IG,LK,J-1,2)*DBLE(J-1))*AAI(IG)
  420 CONTINUE
  430 IF(IS) GO TO 460
      DO 440 I =2,NI
      DO 440 LK=1,NLMAX*NKMAX
      DO 440 IG=1,NG
      SIJ(IG,LK,1,I)=( GIJ(IG,LK,2,I+1)*AAI(IG)
     1                -GIJ(IG,LK,2,I-1)*DBLE(I-1))*AAJ(IG)
  440 CONTINUE
      IF(JS) GO TO 460
      DO 450 I =2,NI
      DO 450 J =2,NJ
      DO 450 LK=1,NLMAX*NKMAX
      DO 450 IG=1,NG
      SIJ(IG,LK,J,I)=( GIJ(IG,LK,J+1,I+1)*AAJ(IG)
     1                -GIJ(IG,LK,J-1,I+1)*DBLE(J-1))*AAI(IG)
     2              +(-GIJ(IG,LK,J+1,I-1)*AAJ(IG)
     3                +GIJ(IG,LK,J-1,I-1)*DBLE(J-1))*DBLE(I-1)
  450 CONTINUE
C
  460 CONTINUE
      IF(SKIPK) GO TO 560
C
C     ----- -SIK- -----
C
      DO 510 J =1,NJMAX
      DO 510 L =1,NLMAX
      DO 510 IG=1,NG
  510 SIK(IG,L,1,J,1)=HIJK(IG,L,2,J,2)*AAK(IG)*AAI(IG)
      IF(KS) GO TO 530
      DO 520 J =1,NJMAX
      DO 520 K =2,NK
      DO 520 L =1,NLMAX
      DO 520 IG=1,NG
      SIK(IG,L,K,J,1)=( HIJK(IG,L,K+1,J,2)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J,2)*DBLE(K-1))*AAI(IG)
  520 CONTINUE
  530 IF(IS) GO TO 560
      DO 540 I =2,NI
      DO 540 J =1,NJMAX
      DO 540 L =1,NLMAX
      DO 540 IG=1,NG
      SIK(IG,L,1,J,I)=( HIJK(IG,L,2,J,I+1)*AAI(IG)
     1                 -HIJK(IG,L,2,J,I-1)*DBLE(I-1))*AAK(IG)
  540 CONTINUE
      IF(KS) GO TO 560
      DO 550 I =2,NI
      DO 550 J =1,NJMAX
      DO 550 K =2,NK
      DO 550 L =1,NLMAX
      DO 550 IG=1,NG
      SIK(IG,L,K,J,I)=( HIJK(IG,L,K+1,J,I+1)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J,I+1)*DBLE(K-1))*AAI(IG)
     2               +(-HIJK(IG,L,K+1,J,I-1)*AAK(IG)
     3                 +HIJK(IG,L,K-1,J,I-1)*DBLE(K-1))*DBLE(I-1)
  550 CONTINUE
C
  560 CONTINUE
      IF(SKIPL) GO TO 660
C
C     ----- -SIL- -----
C
      DO 610 KJ=1,NKMAX*NJMAX
      DO 610 IG=1,NG
  610 SIL(IG,1,KJ,1)=HIJKL(IG,2,KJ,2)*AAL(IG)*AAI(IG)
      IF(LS) GO TO 630
      DO 620 KJ=1,NKMAX*NJMAX
      DO 620 L =2,NL
      DO 620 IG=1,NG
      SIL(IG,L,KJ,1)=( HIJKL(IG,L+1,KJ,2)*AAL(IG)
     1                -HIJKL(IG,L-1,KJ,2)*DBLE(L-1))*AAI(IG)
  620 CONTINUE
  630 IF(IS) GO TO 660
      DO 640 I =2,NI
      DO 640 KJ=1,NKMAX*NJMAX
      DO 640 IG=1,NG
      SIL(IG,1,KJ,I)=( HIJKL(IG,2,KJ,I+1)*AAI(IG)
     1                -HIJKL(IG,2,KJ,I-1)*DBLE(I-1))*AAL(IG)
  640 CONTINUE
      IF(LS) GO TO 660
      DO 650 I =2,NI
      DO 650 KJ=1,NKMAX*NJMAX
      DO 650 L =2,NL
      DO 650 IG=1,NG
      SIL(IG,L,KJ,I)=( HIJKL(IG,L+1,KJ,I+1)*AAL(IG)
     1                -HIJKL(IG,L-1,KJ,I+1)*DBLE(L-1))*AAI(IG)
     2              +(-HIJKL(IG,L+1,KJ,I-1)*AAL(IG)
     3                +HIJKL(IG,L-1,KJ,I-1)*DBLE(L-1))*DBLE(I-1)
  650 CONTINUE
C
  660 CONTINUE
      IF(SKIPJ) GO TO 860
      IF(SKIPK) GO TO 760
C
C     ----- -SJK- -----
C
      DO 710 I =1,NIMAX
      DO 710 L =1,NLMAX
      DO 710 IG=1,NG
  710 SJK(IG,L,1,1,I)=HIJK(IG,L,2,2,I)*AAK(IG)*AAJ(IG)
      IF(KS) GO TO 730
      DO 720 I =1,NIMAX
      DO 720 K =2,NK
      DO 720 L =1,NLMAX
      DO 720 IG=1,NG
      SJK(IG,L,K,1,I)=( HIJK(IG,L,K+1,2,I)*AAK(IG)
     1                 -HIJK(IG,L,K-1,2,I)*DBLE(K-1))*AAJ(IG)
  720 CONTINUE
  730 IF(JS) GO TO 760
      DO 740 I =1,NIMAX
      DO 740 J =2,NJ
      DO 740 L =1,NLMAX
      DO 740 IG=1,NG
      SJK(IG,L,1,J,I)=( HIJK(IG,L,2,J+1,I)*AAJ(IG)
     1                 -HIJK(IG,L,2,J-1,I)*DBLE(J-1))*AAK(IG)
  740 CONTINUE
      IF(KS) GO TO 760
      DO 750 I =1,NIMAX
      DO 750 J =2,NJ
      DO 750 K =2,NK
      DO 750 L =1,NLMAX
      DO 750 IG=1,NG
      SJK(IG,L,K,J,I)=( HIJK(IG,L,K+1,J+1,I)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J+1,I)*DBLE(K-1))*AAJ(IG)
     2               +(-HIJK(IG,L,K+1,J-1,I)*AAK(IG)
     3                 +HIJK(IG,L,K-1,J-1,I)*DBLE(K-1))*DBLE(J-1)
  750 CONTINUE
C
  760 CONTINUE
      IF(SKIPL) GO TO 860
C
C     ----- -SJL- -----
C
      DO 810 I =1,NIMAX
      DO 810 K =1,NKMAX
      DO 810 IG=1,NG
  810 SJL(IG,1,K,1,I)=FIJKL(IG,2,K,2,I)*AAL(IG)*AAJ(IG)
      IF(LS) GO TO 830
      DO 820 I =1,NIMAX
      DO 820 K =1,NKMAX
      DO 820 L =2,NL
      DO 820 IG=1,NG
      SJL(IG,L,K,1,I)=( FIJKL(IG,L+1,K,2,I)*AAL(IG)
     1                 -FIJKL(IG,L-1,K,2,I)*DBLE(L-1))*AAJ(IG)
  820 CONTINUE
  830 IF(JS) GO TO 860
      DO 840 I =1,NIMAX
      DO 840 J =2,NJ
      DO 840 K =1,NKMAX
      DO 840 IG=1,NG
      SJL(IG,1,K,J,I)=( FIJKL(IG,2,K,J+1,I)*AAJ(IG)
     1                 -FIJKL(IG,2,K,J-1,I)*DBLE(J-1))*AAL(IG)
  840 CONTINUE
      IF(LS) GO TO 860
      DO 850 I =1,NIMAX
      DO 850 J =2,NJ
      DO 850 K =1,NKMAX
      DO 850 L =2,NL
      DO 850 IG=1,NG
      SJL(IG,L,K,J,I)=( FIJKL(IG,L+1,K,J+1,I)*AAL(IG)
     1                 -FIJKL(IG,L-1,K,J+1,I)*DBLE(L-1))*AAJ(IG)
     2               +(-FIJKL(IG,L+1,K,J-1,I)*AAL(IG)
     3                 +FIJKL(IG,L-1,K,J-1,I)*DBLE(L-1))*DBLE(J-1)
  850 CONTINUE
C
  860 CONTINUE
      IF(SKIPK) GO TO 960
      IF(SKIPL) GO TO 960
C
C     ----- -SKL- -----
C
      DO 910 JI=1,NJMAX*NIMAX
      DO 910 IG=1,NG
  910 SKL(IG,1,1,JI)=FIJK(IG,2,2,JI)*AAL(IG)*AAK(IG)
      IF(LS) GO TO 930
      DO 920 JI=1,NJMAX*NIMAX
      DO 920 L =2,NL
      DO 920 IG=1,NG
      SKL(IG,L,1,JI)=( FIJK(IG,L+1,2,JI)*AAL(IG)
     1                -FIJK(IG,L-1,2,JI)*DBLE(L-1))*AAK(IG)
  920 CONTINUE
  930 IF(KS) GO TO 960
      DO 940 JI=1,NJMAX*NIMAX
      DO 940 K =2,NK
      DO 940 IG=1,NG
      SKL(IG,1,K,JI)=( FIJK(IG,2,K+1,JI)*AAK(IG)
     1                -FIJK(IG,2,K-1,JI)*DBLE(K-1))*AAL(IG)
  940 CONTINUE
      IF(LS) GO TO 960
      DO 950 JI=1,NJMAX*NIMAX
      DO 950 K =2,NK
      DO 950 L =2,NL
      DO 950 IG=1,NG
      SKL(IG,L,K,JI)=( FIJK(IG,L+1,K+1,JI)*AAL(IG)
     1                -FIJK(IG,L-1,K+1,JI)*DBLE(L-1))*AAK(IG)
     2              +(-FIJK(IG,L+1,K-1,JI)*AAL(IG)
     3                +FIJK(IG,L-1,K-1,JI)*DBLE(L-1))*DBLE(K-1)
  950 CONTINUE
C
  960 CONTINUE
      RETURN
 1000 CONTINUE
C
C     ----- FIRST DERIVATIVES ONLY -----
C
      IF(SKIPI) GO TO 1030
C
C     ----- -FI- ONLY -----
C
      DO 1010 LKJ=1,NLMAX*NKMAX*NJMAX
      DO 1010 IG =1,NG
       FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
 1010 CONTINUE
      IF(IS) GO TO 1030
      DO 1020 I  =2,NI
      DO 1020 LKJ=1,NLMAX*NKMAX*NJMAX
      DO 1020 IG =1,NG
       FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)
     1              -GI(IG,LKJ,I-1)*DBLE(I-1)
 1020 CONTINUE
C
 1030 CONTINUE
      IF(SKIPJ) GO TO 1130
C
C     ----- -FJ- ONLY -----
C
      DO 1110 I =1,NIMAX
      DO 1110 LK=1,NLMAX*NKMAX
      DO 1110 IG=1,NG
       FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
 1110 CONTINUE
      IF(JS) GO TO 1130
      DO 1120 I =1,NIMAX
      DO 1120 J =2,NJ
      DO 1120 LK=1,NLMAX*NKMAX
      DO 1120 IG=1,NG
       FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG)
     1               -GIJ(IG,LK,J-1,I)*DBLE(J-1)
 1120 CONTINUE
C
 1130 CONTINUE
      IF(SKIPK) GO TO 1230
C
C     ----- -FK- ONLY -----
C
      DO 1210 JI=1,NJMAX*NIMAX
      DO 1210 L =1,NLMAX
      DO 1210 IG=1,NG
       FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
 1210 CONTINUE
      IF(KS) GO TO 1230
      DO 1220 JI=1,NJMAX*NIMAX
      DO 1220 K =2,NK
      DO 1220 L =1,NLMAX
      DO 1220 IG=1,NG
       FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)
     1               -GIJK(IG,L,K-1,JI)*DBLE(K-1)
 1220 CONTINUE
C
 1230 CONTINUE
      IF(SKIPL) GO TO 1330
C
C     ----- -FL- AND -SLL- -----
C
      DO 1310 KJI=1,NKMAX*NJMAX*NIMAX
      DO 1310 IG =1,NG
       FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
 1310 CONTINUE
      IF(LS) GO TO 1330
      DO 1320 KJI=1,NKMAX*NJMAX*NIMAX
      DO 1320 L  =2,NL
      DO 1320 IG =1,NG
       FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)
     1              -GIJKL(IG,L-1,KJI)*DBLE(L-1)
 1320 CONTINUE
C
 1330 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JDXYZV(GI,GIJ,GIJK,GIJKL,HIJK,HIJKL,FIJK,
     1 FIJKL,
     1 NG,NIMAX,NJMAX,NKMAX,NLMAX,NI,NJ,NK,NL,AAI,AAJ,AAK,AAL,
     2 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IS,JS,KS,LS,IP,JP,KP,LP
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,
     1 SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
C
C     ----- VERSION 1 -----
C
      DIMENSION GIJKL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION  GIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   GIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    GI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION HIJKL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION  HIJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION FIJKL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION  FIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   AAI(NG)
      DIMENSION   AAJ(NG)
      DIMENSION   AAK(NG)
      DIMENSION   AAL(NG)
      DIMENSION    FI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION    FJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    FK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION    FL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SII(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION   SJJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SKK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   SLL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SIK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SIL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION   SJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SJL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SKL(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DATA THREE  /3.0D+00/
C
      IS=NI.EQ.1
      JS=NJ.EQ.1
      KS=NK.EQ.1
      LS=NL.EQ.1
      IP=NI.EQ.2
      JP=NJ.EQ.2
      KP=NK.EQ.2
      LP=NL.EQ.2
C
      IF(NDER.EQ.1) GO TO 1000
C
C     ----- FIRST AND SECOND DERIVATIVES -----
C
      IF(SKIPI) GO TO 100
C
C     ----- -FI- AND -SII- -----
C
      DO  10 IG =1,NG
      DO  10 LKJ=1,NLMAX*NKMAX*NJMAX
       FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
      SII(IG,LKJ,1)=( GI(IG,LKJ,3)*AAI(IG)
     1               -GI(IG,LKJ,1)        )*AAI(IG)
   10 CONTINUE
      IF(IS) GO TO 100
      DO  20 IG =1,NG
      DO  20 LKJ=1,NLMAX*NKMAX*NJMAX
       FI(IG,LKJ,2)= GI(IG,LKJ,3)*AAI(IG)
     1              -GI(IG,LKJ,1)
      SII(IG,LKJ,2)=(GI(IG,LKJ,4)*AAI(IG)
     1              -GI(IG,LKJ,2)*THREE  )*AAI(IG)
   20 CONTINUE
      IF(IP) GO TO 100
      DO  30 IG =1,NG
      DO  30 I  =3,NI
      DO  30 LKJ=1,NLMAX*NKMAX*NJMAX
       FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)
     1              -GI(IG,LKJ,I-1)*DBLE(I-1)
      SII(IG,LKJ,I)=(GI(IG,LKJ,I+2)*AAI(IG)
     1              -GI(IG,LKJ,I  )*DBLE(I+I-1)     )*AAI(IG)
     2              +GI(IG,LKJ,I-2)*DBLE((I-1)*(I-2))
   30 CONTINUE
C
  100 CONTINUE
      IF(SKIPJ) GO TO 200
C
C     ----- -FJ- AND -SJJ- -----
C
      DO 110 IG=1,NG
      DO 110 I =1,NIMAX
      DO 110 LK=1,NLMAX*NKMAX
       FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
      SJJ(IG,LK,1,I)=( GIJ(IG,LK,3,I)*AAJ(IG)
     1                -GIJ(IG,LK,1,I)        )*AAJ(IG)
  110 CONTINUE
      IF(JS) GO TO 200
      DO 120 IG=1,NG
      DO 120 I =1,NIMAX
      DO 120 LK=1,NLMAX*NKMAX
       FJ(IG,LK,2,I)= GIJ(IG,LK,3,I)*AAJ(IG)
     1               -GIJ(IG,LK,1,I)
      SJJ(IG,LK,2,I)=(GIJ(IG,LK,4,I)*AAJ(IG)
     1               -GIJ(IG,LK,2,I)*THREE  )*AAJ(IG)
  120 CONTINUE
      IF(JP) GO TO 200
      DO 130 IG=1,NG
      DO 130 I =1,NIMAX
      DO 130 J =3,NJ
      DO 130 LK=1,NLMAX*NKMAX
       FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG)
     1               -GIJ(IG,LK,J-1,I)*DBLE(J-1)
      SJJ(IG,LK,J,I)=(GIJ(IG,LK,J+2,I)*AAJ(IG)
     1               -GIJ(IG,LK,J  ,I)*DBLE(J+J-1)     )*AAJ(IG)
     2               +GIJ(IG,LK,J-2,I)*DBLE((J-1)*(J-2))
  130 CONTINUE
C
  200 CONTINUE
      IF(SKIPK) GO TO 300
C
C     ----- -FK- AND -SKK- -----
C
      DO 210 IG=1,NG
      DO 210 JI=1,NJMAX*NIMAX
      DO 210 L =1,NLMAX
       FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
      SKK(IG,L,1,JI)=( GIJK(IG,L,3,JI)*AAK(IG)
     1                -GIJK(IG,L,1,JI)        )*AAK(IG)
  210 CONTINUE
      IF(KS) GO TO 300
      DO 220 IG=1,NG
      DO 220 JI=1,NJMAX*NIMAX
      DO 220 L =1,NLMAX
       FK(IG,L,2,JI)= GIJK(IG,L,3,JI)*AAK(IG)
     1               -GIJK(IG,L,1,JI)
      SKK(IG,L,2,JI)=(GIJK(IG,L,4,JI)*AAK(IG)
     1               -GIJK(IG,L,2,JI)*THREE  )*AAK(IG)
  220 CONTINUE
      IF(KP) GO TO 300
      DO 230 IG=1,NG
      DO 230 JI=1,NJMAX*NIMAX
      DO 230 K =3,NK
      DO 230 L =1,NLMAX
       FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)
     1               -GIJK(IG,L,K-1,JI)*DBLE(K-1)
      SKK(IG,L,K,JI)=(GIJK(IG,L,K+2,JI)*AAK(IG)
     1               -GIJK(IG,L,K  ,JI)*DBLE(K+K-1)     )*AAK(IG)
     2               +GIJK(IG,L,K-2,JI)*DBLE((K-1)*(K-2))
  230 CONTINUE
C
  300 CONTINUE
      IF(SKIPL) GO TO 400
C
C     ----- -FL- AND -SLL- -----
C
      DO 310 IG =1,NG
      DO 310 KJI=1,NKMAX*NJMAX*NIMAX
       FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
      SLL(IG,1,KJI)=( GIJKL(IG,3,KJI)*AAL(IG)
     1               -GIJKL(IG,1,KJI)        )*AAL(IG)
  310 CONTINUE
      IF(LS) GO TO 400
      DO 320 IG =1,NG
      DO 320 KJI=1,NKMAX*NJMAX*NIMAX
       FL(IG,2,KJI)= GIJKL(IG,3,KJI)*AAL(IG)
     1              -GIJKL(IG,1,KJI)
      SLL(IG,2,KJI)=(GIJKL(IG,4,KJI)*AAL(IG)
     1              -GIJKL(IG,2,KJI)*THREE  )*AAL(IG)
  320 CONTINUE
      IF(LP) GO TO 400
      DO 330 IG =1,NG
      DO 330 KJI=1,NKMAX*NJMAX*NIMAX
      DO 330 L  =3,NL
       FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)
     1              -GIJKL(IG,L-1,KJI)*DBLE(L-1)
      SLL(IG,L,KJI)=(GIJKL(IG,L+2,KJI)*AAL(IG)
     1              -GIJKL(IG,L  ,KJI)*DBLE(L+L-1)     )*AAL(IG)
     2              +GIJKL(IG,L-2,KJI)*DBLE((L-1)*(L-2))
  330 CONTINUE
C
  400 CONTINUE
      IF(SKIPI) GO TO 700
      IF(SKIPJ) GO TO 500
C
C     ----- -SIJ- -----
C
      DO 410 IG=1,NG
      DO 410 LK=1,NLMAX*NKMAX
  410 SIJ(IG,LK,1,1)=GIJ(IG,LK,2,2)*AAJ(IG)*AAI(IG)
      IF(JS) GO TO 430
      DO 420 IG=1,NG
      DO 420 J =2,NJ
      DO 420 LK=1,NLMAX*NKMAX
      SIJ(IG,LK,J,1)=( GIJ(IG,LK,J+1,2)*AAJ(IG)
     1                -GIJ(IG,LK,J-1,2)*DBLE(J-1))*AAI(IG)
  420 CONTINUE
  430 IF(IS) GO TO 500
      DO 440 IG=1,NG
      DO 440 I =2,NI
      DO 440 LK=1,NLMAX*NKMAX
      SIJ(IG,LK,1,I)=( GIJ(IG,LK,2,I+1)*AAI(IG)
     1                -GIJ(IG,LK,2,I-1)*DBLE(I-1))*AAJ(IG)
  440 CONTINUE
      IF(JS) GO TO 500
      DO 450 IG=1,NG
      DO 450 I =2,NI
      DO 450 J =2,NJ
      DO 450 LK=1,NLMAX*NKMAX
      SIJ(IG,LK,J,I)=( GIJ(IG,LK,J+1,I+1)*AAJ(IG)
     1                -GIJ(IG,LK,J-1,I+1)*DBLE(J-1))*AAI(IG)
     2              +(-GIJ(IG,LK,J+1,I-1)*AAJ(IG)
     3                +GIJ(IG,LK,J-1,I-1)*DBLE(J-1))*DBLE(I-1)
  450 CONTINUE
C
  500 CONTINUE
      IF(SKIPK) GO TO 600
C
C     ----- -SIK- -----
C
      DO 510 IG=1,NG
      DO 510 J =1,NJMAX
      DO 510 L =1,NLMAX
  510 SIK(IG,L,1,J,1)=HIJK(IG,L,2,J,2)*AAK(IG)*AAI(IG)
      IF(KS) GO TO 530
      DO 520 IG=1,NG
      DO 520 J =1,NJMAX
      DO 520 K =2,NK
      DO 520 L =1,NLMAX
      SIK(IG,L,K,J,1)=( HIJK(IG,L,K+1,J,2)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J,2)*DBLE(K-1))*AAI(IG)
  520 CONTINUE
  530 IF(IS) GO TO 600
      DO 540 IG=1,NG
      DO 540 I =2,NI
      DO 540 J =1,NJMAX
      DO 540 L =1,NLMAX
      SIK(IG,L,1,J,I)=( HIJK(IG,L,2,J,I+1)*AAI(IG)
     1                 -HIJK(IG,L,2,J,I-1)*DBLE(I-1))*AAK(IG)
  540 CONTINUE
      IF(KS) GO TO 600
      DO 550 IG=1,NG
      DO 550 I =2,NI
      DO 550 J =1,NJMAX
      DO 550 K =2,NK
      DO 550 L =1,NLMAX
      SIK(IG,L,K,J,I)=( HIJK(IG,L,K+1,J,I+1)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J,I+1)*DBLE(K-1))*AAI(IG)
     2               +(-HIJK(IG,L,K+1,J,I-1)*AAK(IG)
     3                 +HIJK(IG,L,K-1,J,I-1)*DBLE(K-1))*DBLE(I-1)
  550 CONTINUE
C
  600 CONTINUE
      IF(SKIPL) GO TO 700
C
C     ----- -SIL- -----
C
      DO 610 IG=1,NG
      DO 610 KJ=1,NKMAX*NJMAX
  610 SIL(IG,1,KJ,1)=HIJKL(IG,2,KJ,2)*AAL(IG)*AAI(IG)
      IF(LS) GO TO 630
      DO 620 IG=1,NG
      DO 620 KJ=1,NKMAX*NJMAX
      DO 620 L =2,NL
      SIL(IG,L,KJ,1)=( HIJKL(IG,L+1,KJ,2)*AAL(IG)
     1                -HIJKL(IG,L-1,KJ,2)*DBLE(L-1))*AAI(IG)
  620 CONTINUE
  630 IF(IS) GO TO 700
      DO 640 IG=1,NG
      DO 640 I =2,NI
      DO 640 KJ=1,NKMAX*NJMAX
      SIL(IG,1,KJ,I)=( HIJKL(IG,2,KJ,I+1)*AAI(IG)
     1                -HIJKL(IG,2,KJ,I-1)*DBLE(I-1))*AAL(IG)
  640 CONTINUE
      IF(LS) GO TO 700
      DO 650 IG=1,NG
      DO 650 I =2,NI
      DO 650 KJ=1,NKMAX*NJMAX
      DO 650 L =2,NL
      SIL(IG,L,KJ,I)=( HIJKL(IG,L+1,KJ,I+1)*AAL(IG)
     1                -HIJKL(IG,L-1,KJ,I+1)*DBLE(L-1))*AAI(IG)
     2              +(-HIJKL(IG,L+1,KJ,I-1)*AAL(IG)
     3                +HIJKL(IG,L-1,KJ,I-1)*DBLE(L-1))*DBLE(I-1)
  650 CONTINUE
C
  700 CONTINUE
      IF(SKIPJ) GO TO 900
      IF(SKIPK) GO TO 800
C
C     ----- -SJK- -----
C
      DO 710 IG=1,NG
      DO 710 I =1,NIMAX
      DO 710 L =1,NLMAX
  710 SJK(IG,L,1,1,I)=HIJK(IG,L,2,2,I)*AAK(IG)*AAJ(IG)
      IF(KS) GO TO 730
      DO 720 IG=1,NG
      DO 720 I =1,NIMAX
      DO 720 K =2,NK
      DO 720 L =1,NLMAX
      SJK(IG,L,K,1,I)=( HIJK(IG,L,K+1,2,I)*AAK(IG)
     1                 -HIJK(IG,L,K-1,2,I)*DBLE(K-1))*AAJ(IG)
  720 CONTINUE
  730 IF(JS) GO TO 800
      DO 740 IG=1,NG
      DO 740 I =1,NIMAX
      DO 740 J =2,NJ
      DO 740 L =1,NLMAX
      SJK(IG,L,1,J,I)=( HIJK(IG,L,2,J+1,I)*AAJ(IG)
     1                 -HIJK(IG,L,2,J-1,I)*DBLE(J-1))*AAK(IG)
  740 CONTINUE
      IF(KS) GO TO 800
      DO 750 IG=1,NG
      DO 750 I =1,NIMAX
      DO 750 J =2,NJ
      DO 750 K =2,NK
      DO 750 L =1,NLMAX
      SJK(IG,L,K,J,I)=( HIJK(IG,L,K+1,J+1,I)*AAK(IG)
     1                 -HIJK(IG,L,K-1,J+1,I)*DBLE(K-1))*AAJ(IG)
     2               +(-HIJK(IG,L,K+1,J-1,I)*AAK(IG)
     3                 +HIJK(IG,L,K-1,J-1,I)*DBLE(K-1))*DBLE(J-1)
  750 CONTINUE
C
  800 CONTINUE
      IF(SKIPL) GO TO 900
C
C     ----- -SJL- -----
C
      DO 810 IG=1,NG
      DO 810 I =1,NIMAX
      DO 810 K =1,NKMAX
  810 SJL(IG,1,K,1,I)=FIJKL(IG,2,K,2,I)*AAL(IG)*AAJ(IG)
      IF(LS) GO TO 830
      DO 820 IG=1,NG
      DO 820 I =1,NIMAX
      DO 820 K =1,NKMAX
      DO 820 L =2,NL
      SJL(IG,L,K,1,I)=( FIJKL(IG,L+1,K,2,I)*AAL(IG)
     1                 -FIJKL(IG,L-1,K,2,I)*DBLE(L-1))*AAJ(IG)
  820 CONTINUE
  830 IF(JS) GO TO 900
      DO 840 IG=1,NG
      DO 840 I =1,NIMAX
      DO 840 J =2,NJ
      DO 840 K =1,NKMAX
      SJL(IG,1,K,J,I)=( FIJKL(IG,2,K,J+1,I)*AAJ(IG)
     1                 -FIJKL(IG,2,K,J-1,I)*DBLE(J-1))*AAL(IG)
  840 CONTINUE
      IF(LS) GO TO 900
      DO 850 IG=1,NG
      DO 850 I =1,NIMAX
      DO 850 J =2,NJ
      DO 850 K =1,NKMAX
      DO 850 L =2,NL
      SJL(IG,L,K,J,I)=( FIJKL(IG,L+1,K,J+1,I)*AAL(IG)
     1                 -FIJKL(IG,L-1,K,J+1,I)*DBLE(L-1))*AAJ(IG)
     2               +(-FIJKL(IG,L+1,K,J-1,I)*AAL(IG)
     3                 +FIJKL(IG,L-1,K,J-1,I)*DBLE(L-1))*DBLE(J-1)
  850 CONTINUE
C
  900 CONTINUE
      IF(SKIPK) GO TO 990
      IF(SKIPL) GO TO 990
C
C     ----- -SKL- -----
C
      DO 910 IG=1,NG
      DO 910 JI=1,NJMAX*NIMAX
  910 SKL(IG,1,1,JI)=FIJK(IG,2,2,JI)*AAL(IG)*AAK(IG)
      IF(LS) GO TO 930
      DO 920 IG=1,NG
      DO 920 JI=1,NJMAX*NIMAX
      DO 920 L =2,NL
      SKL(IG,L,1,JI)=( FIJK(IG,L+1,2,JI)*AAL(IG)
     1                -FIJK(IG,L-1,2,JI)*DBLE(L-1))*AAK(IG)
  920 CONTINUE
  930 IF(KS) GO TO 990
      DO 940 IG=1,NG
      DO 940 JI=1,NJMAX*NIMAX
      DO 940 K =2,NK
      SKL(IG,1,K,JI)=( FIJK(IG,2,K+1,JI)*AAK(IG)
     1                -FIJK(IG,2,K-1,JI)*DBLE(K-1))*AAL(IG)
  940 CONTINUE
      IF(LS) GO TO 990
      DO 950 IG=1,NG
      DO 950 JI=1,NJMAX*NIMAX
      DO 950 K =2,NK
      DO 950 L =2,NL
      SKL(IG,L,K,JI)=( FIJK(IG,L+1,K+1,JI)*AAL(IG)
     1                -FIJK(IG,L-1,K+1,JI)*DBLE(L-1))*AAK(IG)
     2              +(-FIJK(IG,L+1,K-1,JI)*AAL(IG)
     3                +FIJK(IG,L-1,K-1,JI)*DBLE(L-1))*DBLE(K-1)
  950 CONTINUE
C
  990 CONTINUE
      RETURN
 1000 CONTINUE
C
C     ----- FIRST DERIVATIVES ONLY -----
C
      IF(SKIPI) GO TO 1030
C
C     ----- -FI- ONLY -----
C
      DO 1010 IG =1,NG
      DO 1010 LKJ=1,NLMAX*NKMAX*NJMAX
       FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
 1010 CONTINUE
      IF(IS) GO TO 1030
      DO 1020 IG =1,NG
      DO 1020 I  =2,NI
      DO 1020 LKJ=1,NLMAX*NKMAX*NJMAX
       FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)
     1              -GI(IG,LKJ,I-1)*DBLE(I-1)
 1020 CONTINUE
C
 1030 CONTINUE
      IF(SKIPJ) GO TO 1130
C
C     ----- -FJ- ONLY -----
C
      DO 1110 IG=1,NG
      DO 1110 I =1,NIMAX
      DO 1110 LK=1,NLMAX*NKMAX
       FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
 1110 CONTINUE
      IF(JS) GO TO 1130
      DO 1120 IG=1,NG
      DO 1120 I =1,NIMAX
      DO 1120 J =2,NJ
      DO 1120 LK=1,NLMAX*NKMAX
       FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG)
     1               -GIJ(IG,LK,J-1,I)*DBLE(J-1)
 1120 CONTINUE
C
 1130 CONTINUE
      IF(SKIPK) GO TO 1230
C
C     ----- -FK- ONLY -----
C
      DO 1210 IG=1,NG
      DO 1210 JI=1,NJMAX*NIMAX
      DO 1210 L =1,NLMAX
       FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
 1210 CONTINUE
      IF(KS) GO TO 1230
      DO 1220 IG=1,NG
      DO 1220 JI=1,NJMAX*NIMAX
      DO 1220 K =2,NK
      DO 1220 L =1,NLMAX
       FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)
     1               -GIJK(IG,L,K-1,JI)*DBLE(K-1)
 1220 CONTINUE
C
 1230 CONTINUE
      IF(SKIPL) GO TO 1330
C
C     ----- -FL- AND -SLL- -----
C
      DO 1310 IG =1,NG
      DO 1310 KJI=1,NKMAX*NJMAX*NIMAX
       FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
 1310 CONTINUE
      IF(LS) GO TO 1330
      DO 1320 IG =1,NG
      DO 1320 KJI=1,NKMAX*NJMAX*NIMAX
      DO 1320 L  =2,NL
       FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)
     1              -GIJKL(IG,L-1,KJI)*DBLE(L-1)
 1320 CONTINUE
C
 1330 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JKDZER(GIJKL,FIJKL,SIJKL,IJKLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DATA ZERO /0.0D+00/
C
C     ----- ZERO OUT -GIJKL- -----
C
      IF(.NOT.OUT) GO TO 15
      DO 10 N=1,IJKL
      GIJKL(IJKLG(1,N))=ZERO
   10 CONTINUE
   15 CONTINUE
      IF(NDER.EQ.0) RETURN
C
C     ----- ZERO OUT -FIJKL- -----
C
      IF(.NOT.OUT.AND.NDER.NE.2) GO TO 35
      DO 30 N=1,IJKL
      DO 20 M=1,12
      FIJKL(M,IJKLG(1,N))=ZERO
   20 CONTINUE
   30 CONTINUE
C
   35 CONTINUE
      DO 40 I=1,12
   40 FD(I)=ZERO
      IF(NDER.EQ.1) RETURN
C
C     ----- ZERO OUT -SIJKL- -----
C
      IF(.NOT.OUT) GO TO 65
      DO 60 N=1,IJKL
      DO 50 M=1,78
      SIJKL(M,IJKLG(1,N))=ZERO
   50 CONTINUE
   60 CONTINUE
C
   65 CONTINUE
      DO 70 I=1,78
   70 SD(I)=ZERO
      RETURN
      END
C
      SUBROUTINE HND_DSPDFS(NG,NR,IJKLG,GIJKL,FIJKL,SIJKL,XYZ,
     1 FIXYZ,FJXYZ,FKXYZ,FLXYZ,SIIXYZ,SJJXYZ,SKKXYZ,SLLXYZ,
     2 SIJXYZ,SIKXYZ,SILXYZ,SJKXYZ,SJLXYZ,SKLXYZ,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IS,JS,KS,LS
      LOGICAL IJS,IJKS,IJKLS
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1 LGXYZ(4,35)
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION    XYZ(NG*NR,1)
      DIMENSION  FIXYZ(NG*NR,1), FJXYZ(NG*NR,1)
      DIMENSION  FKXYZ(NG*NR,1), FLXYZ(NG*NR,1)
      DIMENSION SIIXYZ(NG*NR,1),SJJXYZ(NG*NR,1)
      DIMENSION SKKXYZ(NG*NR,1),SLLXYZ(NG*NR,1)
      DIMENSION SIJXYZ(NG*NR,1),SIKXYZ(NG*NR,1)
      DIMENSION SILXYZ(NG*NR,1),SJKXYZ(NG*NR,1)
      DIMENSION SJLXYZ(NG*NR,1),SKLXYZ(NG*NR,1)
      DIMENSION XY(NG*NR),XZ(NG*NR),YZ(NG*NR)
      DIMENSION  X(NG*NR), Y(NG*NR), Z(NG*NR)
      DIMENSION SJ(NG*NR),SK(NG*NR),SL(NG*NR)
      DIMENSION DIJSI(NG*NR),DIJSJ(NG*NR),DKLSK(NG*NR),DKLSL(NG*NR)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
C
      IF(SPIJKL) GO TO 1000
C
C     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
C
      IF(NDER.LT.1) GO TO 700
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 640 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 630 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 620 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 610 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 510 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  510 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
C
      IF(SKIPI) GO TO 530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
  520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
  530 IF(SKIPJ) GO TO 550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
  540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
  550 IF(SKIPK) GO TO 570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
  560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
  570 IF(SKIPL) GO TO 600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
  580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
  600 CONTINUE
C
  610 CONTINUE
  620 CONTINUE
  630 CONTINUE
  640 CONTINUE
C
  700 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 940 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 930 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 920 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 910 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 710 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  710 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      DUM=ZERO
      DO 720 IGR=1,NG*NR
  720 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
      IF(SKIPI) GO TO 830
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 820 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
  820 CONTINUE
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
  830 IF(SKIPJ) GO TO 850
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 840 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
  840 CONTINUE
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
  850 IF(SKIPK) GO TO 870
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 860 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
  860 CONTINUE
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
  870 IF(SKIPL) GO TO 900
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 880 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
  880 CONTINUE
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
  900 CONTINUE
C
  910 CONTINUE
  920 CONTINUE
  930 CONTINUE
  940 CONTINUE
C
      RETURN
C
C     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
C
 1000 CONTINUE
      IF(NDER.LT.1) GO TO 3000
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 2640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 2630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 2110 IGR=1,NG*NR
 2110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 2120 IGR=1,NG*NR
 2120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 2130 IGR=1,NG*NR
 2130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 2620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 2210 IGR=1,NG*NR
 2210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 2220 IGR=1,NG*NR
 2220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 2230 IGR=1,NG*NR
 2230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 2610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 2310 IGR=1,NG*NR
 2310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 2320 IGR=1,NG*NR
 2320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 2330 IGR=1,NG*NR
 2330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 2500 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 2500    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 2510 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 2510    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      IF(SKIPI) GO TO 2530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
 2520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
 2530 IF(SKIPJ) GO TO 2550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
 2540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
 2550 IF(SKIPK) GO TO 2570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
 2560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
 2570 IF(SKIPL) GO TO 2600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
 2580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
 2600 CONTINUE
C
 2610 CONTINUE
 2620 CONTINUE
 2630 CONTINUE
 2640 CONTINUE
C
 3000 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 3640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 3630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 3110 IGR=1,NG*NR
 3110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 3120 IGR=1,NG*NR
 3120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 3130 IGR=1,NG*NR
 3130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 3620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 3210 IGR=1,NG*NR
 3210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 3220 IGR=1,NG*NR
 3220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 3230 IGR=1,NG*NR
 3230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 3610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 3310 IGR=1,NG*NR
 3310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 3320 IGR=1,NG*NR
 3320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 3330 IGR=1,NG*NR
 3330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 3400 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 3400    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 3410 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 3410    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      DUM=ZERO
      DO 3420 IGR=1,NG*NR
 3420 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
      IF(SKIPI) GO TO 3530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
 3520 CONTINUE
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
 3530 IF(SKIPJ) GO TO 3550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
 3540 CONTINUE
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
 3550 IF(SKIPK) GO TO 3570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
 3560 CONTINUE
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
 3570 IF(SKIPL) GO TO 3600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
 3580 CONTINUE
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
 3600 CONTINUE
C
 3610 CONTINUE
 3620 CONTINUE
 3630 CONTINUE
 3640 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_DSPDFV(NG,NR,IJKLG,GIJKL,FIJKL,SIJKL,XYZ,
     1 FIXYZ,FJXYZ,FKXYZ,FLXYZ,SIIXYZ,SJJXYZ,SKKXYZ,SLLXYZ,
     2 SIJXYZ,SIKXYZ,SILXYZ,SJKXYZ,SJLXYZ,SKLXYZ,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IS,JS,KS,LS
      LOGICAL IJS,IJKS,IJKLS
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1 LGXYZ(4,35)
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION    XYZ(NG*NR,1)
      DIMENSION  FIXYZ(NG*NR,1), FJXYZ(NG*NR,1)
      DIMENSION  FKXYZ(NG*NR,1), FLXYZ(NG*NR,1)
      DIMENSION SIIXYZ(NG*NR,1),SJJXYZ(NG*NR,1)
      DIMENSION SKKXYZ(NG*NR,1),SLLXYZ(NG*NR,1)
      DIMENSION SIJXYZ(NG*NR,1),SIKXYZ(NG*NR,1)
      DIMENSION SILXYZ(NG*NR,1),SJKXYZ(NG*NR,1)
      DIMENSION SJLXYZ(NG*NR,1),SKLXYZ(NG*NR,1)
      DIMENSION XY(NG*NR),XZ(NG*NR),YZ(NG*NR)
      DIMENSION  X(NG*NR), Y(NG*NR), Z(NG*NR)
      DIMENSION SJ(NG*NR),SK(NG*NR),SL(NG*NR)
      DIMENSION DIJSI(NG*NR),DIJSJ(NG*NR),DKLSK(NG*NR),DKLSL(NG*NR)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
C
      IF(SPIJKL) GO TO 1000
C
C     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
C
      IF(NDER.LT.1) GO TO 700
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 640 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 630 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 620 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 610 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 510 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  510 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
C
      IF(SKIPI) GO TO 530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
  520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
  530 IF(SKIPJ) GO TO 550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
  540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
  550 IF(SKIPK) GO TO 570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
  560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
  570 IF(SKIPL) GO TO 600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
  580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
  600 CONTINUE
C
  610 CONTINUE
  620 CONTINUE
  630 CONTINUE
  640 CONTINUE
C
  700 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 940 I=MINI,MAXI
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 930 J=MINJ,JMAX
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 920 K=MINK,KMAX
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 910 L=MINL,LMAX
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      DO 710 IGR=1,NG*NR
      XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
      XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  710 YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      DUM=ZERO
      DO 720 IGR=1,NG*NR
  720 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
      IF(SKIPI) GO TO 830
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 820 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
  820 CONTINUE
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
  830 IF(SKIPJ) GO TO 850
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 840 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
  840 CONTINUE
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
  850 IF(SKIPK) GO TO 870
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 860 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
  860 CONTINUE
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
  870 IF(SKIPL) GO TO 900
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 880 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
  880 CONTINUE
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
  900 CONTINUE
C
  910 CONTINUE
  920 CONTINUE
  930 CONTINUE
  940 CONTINUE
C
      RETURN
C
C     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
C
 1000 CONTINUE
      IF(NDER.LT.1) GO TO 3000
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 2640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 2630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 2110 IGR=1,NG*NR
 2110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 2120 IGR=1,NG*NR
 2120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 2130 IGR=1,NG*NR
 2130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 2620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 2210 IGR=1,NG*NR
 2210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 2220 IGR=1,NG*NR
 2220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 2230 IGR=1,NG*NR
 2230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 2610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 2310 IGR=1,NG*NR
 2310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 2320 IGR=1,NG*NR
 2320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 2330 IGR=1,NG*NR
 2330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 2500 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 2500    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 2510 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 2510    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      IF(SKIPI) GO TO 2530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
 2520 CONTINUE
      FD( 1)=FD( 1)+DAB(NN)*DUMFX
      FD( 2)=FD( 2)+DAB(NN)*DUMFY
      FD( 3)=FD( 3)+DAB(NN)*DUMFZ
 2530 IF(SKIPJ) GO TO 2550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
 2540 CONTINUE
      FD( 4)=FD( 4)+DAB(NN)*DUMFX
      FD( 5)=FD( 5)+DAB(NN)*DUMFY
      FD( 6)=FD( 6)+DAB(NN)*DUMFZ
 2550 IF(SKIPK) GO TO 2570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
 2560 CONTINUE
      FD( 7)=FD( 7)+DAB(NN)*DUMFX
      FD( 8)=FD( 8)+DAB(NN)*DUMFY
      FD( 9)=FD( 9)+DAB(NN)*DUMFZ
 2570 IF(SKIPL) GO TO 2600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 2580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
 2580 CONTINUE
      FD(10)=FD(10)+DAB(NN)*DUMFX
      FD(11)=FD(11)+DAB(NN)*DUMFY
      FD(12)=FD(12)+DAB(NN)*DUMFZ
 2600 CONTINUE
C
 2610 CONTINUE
 2620 CONTINUE
 2630 CONTINUE
 2640 CONTINUE
C
 3000 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
      IF(.NOT.OUT) RETURN
C
      IJKLN=0
      DO 3640 I=MINI,MAXI
      IS=SPI.AND.I.EQ.1
C
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 3630 J=MINJ,JMAX
      JS=SPJ.AND.J.EQ.1
C
      IF(JS) THEN
         IF(IS) THEN
            DO 3110 IGR=1,NG*NR
 3110       SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
         ELSE
            DO 3120 IGR=1,NG*NR
 3120       SJ(IGR)=DIJSJ(IGR)
         ENDIF
      ELSE
         IF(IS) THEN
            DO 3130 IGR=1,NG*NR
 3130       SJ(IGR)=DIJSI(IGR)
         ENDIF
      ENDIF
      IJS=IS.OR.JS
C
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 3620 K=MINK,KMAX
      KS=SPK.AND.K.EQ.1
C
      IF(KS) THEN
         IF(IJS) THEN
            DO 3210 IGR=1,NG*NR
 3210       SK(IGR)=DKLSK(IGR)*SJ(IGR)
         ELSE
            DO 3220 IGR=1,NG*NR
 3220       SK(IGR)=DKLSK(IGR)
         ENDIF
      ELSE
         IF(IJS) THEN
            DO 3230 IGR=1,NG*NR
 3230       SK(IGR)=SJ(IGR)
         ENDIF
      ENDIF
      IJKS=IJS.OR.KS
C
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 3610 L=MINL,LMAX
      LS=SPL.AND.L.EQ.1
C
      IF(LS) THEN
         IF(IJKS) THEN
            DO 3310 IGR=1,NG*NR
 3310       SL(IGR)=DKLSL(IGR)*SK(IGR)
         ELSE
            DO 3320 IGR=1,NG*NR
 3320       SL(IGR)=DKLSL(IGR)
         ENDIF
      ELSE
         IF(IJKS) THEN
            DO 3330 IGR=1,NG*NR
 3330       SL(IGR)=SK(IGR)
         ENDIF
      ENDIF
      IJKLS=IJKS.OR.LS
C
      IJKLN=IJKLN+1
      NN=IJKLG(1,IJKLN)
      NX=IJKLG(2,IJKLN)
      NY=IJKLG(3,IJKLN)
      NZ=IJKLG(4,IJKLN)
C
      IF(IJKLS) THEN
         DO 3400 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 3400    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
      ELSE
         DO 3410 IGR=1,NG*NR
         XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
         XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 3410    YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
      ENDIF
C
      DUM=ZERO
      DO 3420 IGR=1,NG*NR
 3420 DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
      GIJKL(NN)=GIJKL(NN)+DUM
C
      IF(SKIPI) GO TO 3530
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3520 IGR=1,NG*NR
      DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
 3520 CONTINUE
      FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
      FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
      FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
 3530 IF(SKIPJ) GO TO 3550
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3540 IGR=1,NG*NR
      DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
 3540 CONTINUE
      FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
      FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
      FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
 3550 IF(SKIPK) GO TO 3570
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3560 IGR=1,NG*NR
      DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
 3560 CONTINUE
      FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
      FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
      FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
 3570 IF(SKIPL) GO TO 3600
      DUMFX =ZERO
      DUMFY =ZERO
      DUMFZ =ZERO
      DO 3580 IGR=1,NG*NR
      DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
      DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
      DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
 3580 CONTINUE
      FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
      FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
      FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
 3600 CONTINUE
C
 3610 CONTINUE
 3620 CONTINUE
 3630 CONTINUE
 3640 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKDINV(FIJKL,SIJKL,IJKLG,FC,NC,DF,LDF,NBF,DEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DERINV/INVTYP
      COMMON/HND_DERSHL/FD(3,4),SD(78)
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      DIMENSION FIJKL(3,4,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION FC(NC,1),DF(LDF,1),DEN(1)
C     DIMENSION KAT(4),KC(4)
      DIMENSION KAT(4)
      DIMENSION DD(12,12),DS(3,4,3,4)
      EQUIVALENCE (DD(1,1),DS(1,1,1,1))
      EQUIVALENCE (KAT(1),IIAT)
C
      IF(NDER.EQ.0) RETURN
C
C     ----- TRANSLATIONAL INVARIANCE FOR GRADIENT ELEMENTS -----
C
      GO TO (2000, 200, 300, 400, 400, 600, 700, 800,
     1        900, 900,1100,1200,1300,1400,1500,2000),INVTYP
C
  200 CONTINUE
      DO 210 IXYZ=1,3
  210 FD(IXYZ,1)=- FD(IXYZ,4)
      GO TO 2000
C
  300 CONTINUE
      DO 310 IXYZ=1,3
  310 FD(IXYZ,1)=- FD(IXYZ,3)
      GO TO 2000
C
  400 CONTINUE
      DO 410 IXYZ=1,3
  410 FD(IXYZ,1)=-(FD(IXYZ,3)+FD(IXYZ,4))
      GO TO 2000
C
  600 CONTINUE
      DO 610 IXYZ=1,3
  610 FD(IXYZ,1)=- FD(IXYZ,2)
      GO TO 2000
C
  700 CONTINUE
      DO 710 IXYZ=1,3
  710 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      GO TO 2000
C
  800 CONTINUE
      DO 810 IXYZ=1,3
  810 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      GO TO 2000
C
  900 CONTINUE
      DO 910 IXYZ=1,3
  910 FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,3))
      GO TO 2000
C
 1100 CONTINUE
      DO 1110 IXYZ=1,3
 1110 FD(IXYZ,2)=- FD(IXYZ,1)
      GO TO 2000
C
 1200 CONTINUE
      DO 1210 IXYZ=1,3
 1210 FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,4))
      GO TO 2000
C
 1300 CONTINUE
      DO 1310 IXYZ=1,3
 1310 FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,3))
      GO TO 2000
C
 1400 CONTINUE
      DO 1410 IXYZ=1,3
 1410 FD(IXYZ,3)=-(FD(IXYZ,1)+FD(IXYZ,2))
      GO TO 2000
C
 1500 CONTINUE
      DO 1510 IXYZ=1,3
 1510 FD(IXYZ,4)=-(FD(IXYZ,1)+FD(IXYZ,2)+FD(IXYZ,3))
C
 2000 CONTINUE
C
      DO 2010 I=1,4
      IAT=KAT(I)
      DO 2010 IXYZ=1,3
 2010 DE(IXYZ,IAT)=DE(IXYZ,IAT)+FD(IXYZ,I)
C
      IF(OUT) WRITE(IW,9999)(I,IIAT,DE(I,IIAT),I,JJAT,DE(I,JJAT),
     1                       I,KKAT,DE(I,KKAT),I,LLAT,DE(I,LLAT),I=1,3)
C
      RETURN
 9999 FORMAT(4(' DE(',I2,',',I3,') = ',E18.11))
      END
C
      SUBROUTINE HND_JKDOUT(GIJKL,FIJKL,SIJKL,IJKLG,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXNBF =2048)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL OUT
      LOGICAL BUG1,BUG2,BUG
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/HND_DERSHL/FD(12),SD(78)
      COMMON/HND_SHLPAR/LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      COMMON/HND_SHLNUM/NUMI,NUMJ,NUMK,NUML,IJ,KL,IJKL
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION F(12),S(12,12)
      DATA ZERO /0.0D+00/
      DATA TOL  /1.0D-10/
C
      BUG=.FALSE.
C
      IJKLN=0
      DO 400 I=MINI,MAXI
      PI=PNRM(I)
      JMAX=MAXJ
      IF(IIEQJJ) JMAX=I
      DO 300 J=MINJ,JMAX
      PIJ=PNRM(J)*PI
      KMAX=MAXK
      IF(IJEQKL) KMAX=I
      DO 200 K=MINK,KMAX
      PIJK=PNRM(K)*PIJ
      LMAX=MAXL
      IF(KKEQLL           ) LMAX=K
      IF(IJEQKL.AND.K.EQ.I) LMAX=J
      DO 100 L=MINL,LMAX
      PIJKL=PNRM(L)*PIJK
C
      IJKLN=IJKLN+1
      INUM=LOCI+I
      JNUM=LOCJ+J
      KNUM=LOCK+K
      LNUM=LOCL+L
C
      VAL=GIJKL(IJKLG(1,IJKLN))*PIJKL
      DEN=  DAB(IJKLG(1,IJKLN))
C
      IF(NDER.EQ.0) GO TO 50
C
      DO 10 N=1,12
   10 F(N)=FIJKL( N,IJKLG(1,IJKLN))*PIJKL
      ERR1=ZERO
      DO 20 N=1,12
   20 ERR1=ERR1+F(N)
      BUG1= ABS(ERR1).GT.TOL
      IF(BUG1) THEN
      BUG=.TRUE.
      WRITE(IW,9999) INUM,JNUM,KNUM,LNUM,VAL,DEN,
     1 IJKLG(1,IJKLN),IJKLG(2,IJKLN),IJKLG(3,IJKLN),IJKLG(4,IJKLN)
      WRITE(IW,9998) ERR1
      WRITE(IW,9995) (F(N),N=1,12)
      ENDIF
C
      IF(NDER.GE.1) GO TO 50
C
      DO 30 M=1,12
      DO 30 N=1,M
      S(M,N)=SIJKL(IA(M)+N,IJKLG(1,IJKLN))*PIJKL
   30 S(N,M)=SIJKL(IA(M)+N,IJKLG(1,IJKLN))*PIJKL
      ERR2=ZERO
      DO 40 N=1,12
      DO 40 M=1,12
   40 ERR2=ERR2+S(M,N)
      BUG2= ABS(ERR2).GT.TOL
      IF(BUG2) THEN
      BUG=.TRUE.
      WRITE(IW,9999) INUM,JNUM,KNUM,LNUM,VAL,DEN,
     1 IJKLG(1,IJKLN),IJKLG(2,IJKLN),IJKLG(3,IJKLN),IJKLG(4,IJKLN)
      WRITE(IW,9997) ERR2
      WRITE(IW,9994) ((S(M,N),M=1,12),N=1,12)
      ENDIF
C
   50 CONTINUE
C
  100 CONTINUE
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
C
      IF(.NOT.BUG) WRITE(IW,9996)
C
      RETURN
 9999 FORMAT(' I,J,K,L = ',4I3,2E20.12,4I5)
 9998 FORMAT(' ERROR IN 1-DERIVATIVE INTEGRALS , ERR1 = ',F20.12)
 9997 FORMAT(' ERROR IN 2-DERIVATIVE INTEGRALS , ERR2 = ',F20.12)
 9996 FORMAT(' NO -TRANSLATION INVARIANCE ERROR DETECTED',
     1       ' FOR THIS GROUP OF SHELLS. GOOD GOING')
 9995 FORMAT(3E20.12)
 9994 FORMAT(6E15.8)
      END
C
C     ------- THE SUBROUTINE TAKEN FROM SCF.F IN HONDO ------
C
      SUBROUTINE HND_RFCKR(T,MINK,MAXK,LKT,MINL,MAXL,LLT,NTR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NDIM=35)
      COMMON/HND_SYMSPD/PTR(3,144),DTR(6,288),FTR(10,480),GTR(15,720)
      DIMENSION T(NDIM,NDIM),V(NDIM)
      DATA ZERO /0.0D+00/
C
C     ----- RIGHT MULTIPLY  T  BY  R , RESULT BACK IN  T -----
C
      GO TO (500,400,300,200,100),LLT
C
C     ----- G SHELL -----
C
  100 NG=15*(NTR-1)
      DO 130 K=MINK,MAXK
      DO 120 L=21,35
      DUM=ZERO
      DO 110 N=21,35
  110 DUM=DUM+T(K,N)*GTR(N-20,L-20+NG)
  120 V(L)=DUM
      DO 130 L=21,35
  130 T(K,L)=V(L)
      GO TO 500
C
C     ----- F SHELL -----
C
  200 NF=10*(NTR-1)
      DO 230 K=MINK,MAXK
      DO 220 L=11,20
      DUM=ZERO
      DO 210 N=11,20
  210 DUM=DUM+T(K,N)*FTR(N-10,L-10+NF)
  220 V(L)=DUM
      DO 230 L=11,20
  230 T(K,L)=V(L)
      GO TO 500
C
C     ----- D SHELL -----
C
  300 ND= 6*(NTR-1)
      DO 330 K=MINK,MAXK
      DO 320 L= 5,10
      DUM=ZERO
      DO 310 N= 5,10
  310 DUM=DUM+T(K,N)*DTR(N- 4,L- 4+ND)
  320 V(L)=DUM
      DO 330 L= 5,10
  330 T(K,L)=V(L)
      GO TO 500
C
C     ----- P SHELL -----
C
  400 NP= 3*(NTR-1)
      DO 430 K=MINK,MAXK
      DO 420 L= 2, 4
      DUM=ZERO
      DO 410 N= 2, 4
  410 DUM=DUM+T(K,N)*PTR(N- 1,L- 1+NP)
  420 V(L)=DUM
      DO 430 L= 2, 4
  430 T(K,L)=V(L)
  500 CONTINUE
C
C     ----- LEFT MULTIPLY  T  BY R , RESULT BACK IN  T -----
C
      GO TO (1000,900,800,700,600),LKT
C
C     ----- G SHELL -----
C
  600 NG=15*(NTR-1)
      DO 630 L=MINL,MAXL
      DO 620 K=21,35
      DUM=ZERO
      DO 610 N=21,35
  610 DUM=DUM+GTR(N-20,K-20+NG)*T(N,L)
  620 V(K)=DUM
      DO 630 K=21,35
  630 T(K,L)=V(K)
      GO TO 1000
C
C     ----- F SHELL -----
C
  700 NF=10*(NTR-1)
      DO 730 L=MINL,MAXL
      DO 720 K=11,20
      DUM=ZERO
      DO 710 N=11,20
  710 DUM=DUM+FTR(N-10,K-10+NF)*T(N,L)
  720 V(K)=DUM
      DO 730 K=11,20
  730 T(K,L)=V(K)
      GO TO 1000
C
C     ----- D SHELL -----
C
  800 ND= 6*(NTR-1)
      DO 830 L=MINL,MAXL
      DO 820 K= 5,10
      DUM=ZERO
      DO 810 N= 5,10
  810 DUM=DUM+DTR(N- 4,K- 4+ND)*T(N,L)
  820 V(K)=DUM
      DO 830 K= 5,10
  830 T(K,L)=V(K)
      GO TO 1000
C
C     ----- P SHELL -----
C
  900 NP= 3*(NTR-1)
      DO 930 L=MINL,MAXL
      DO 920 K= 2, 4
      DUM=ZERO
      DO 910 N= 2, 4
  910 DUM=DUM+PTR(N- 1,K- 1+NP)*T(N,L)
  920 V(K)=DUM
      DO 930 K= 2, 4
  930 T(K,L)=V(K)
 1000 CONTINUE
      RETURN
      END
C
C
C     ------ SUBROUTINES TAKE FROM CTL.F IN HONDO
C 
      SUBROUTINE HND_SECOND(SEC,WSEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_CLOCK0/SEC0,WSEC0
      COMMON/HND_CLOCKS/DAT,TIM,CPUSEC,TOTSEC,DELCPU,DELTOT
      DIMENSION DIAL(6)
      EQUIVALENCE (DIAL(1),DAT)
C
C     ----- RETURN ELAPSED - CPU- TIME IN SECONDS IN - SEC- -----
C     ----- RETURN ELAPSED -WALL- TIME IN SECONDS IN -WSEC- -----
C
      CALL HND_SYSCLK(DIAL)
       SEC= CPUSEC
       SEC= SEC- SEC0
      WSEC= TOTSEC
      WSEC=WSEC-WSEC0
      RETURN
      END
C
      SUBROUTINE HND_DAREAD(IDAF,IODA,IX,NX,IDAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_DAFNAV/IXDDAF(2048),NAV10,NAV20
      DIMENSION IODA(1),IX(1)
      DIMENSION IXSDAF(1024)
      EQUIVALENCE (IXSDAF(1),IXDDAF(1))
      DATA IXMAX /1024/
C
      IF(IDAF.EQ.10) NAV10=IODA(IDAR+0)
      IF(IDAF.EQ.20) NAV20=IODA(IDAR+0)
      MAXIX=IXMAX*ISINGL
      LDAR =   NX*ISINGL
C
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) READ(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) READ(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) READ(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) READ(IDAF,REC=NAV20) IXDDAF
      DO 20 I=MIN,MAX
   20 IX(I)=IXDDAF(I-MIN+1)
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 10
      RETURN
      END
C
      SUBROUTINE HND_REWFIL(NFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REWIND NFT
      RETURN
      END
C
C     ------ SUBROUTINE TAKEN FROM CLK.F IN HONDO
C
      SUBROUTINE HND_SYSCLK(DIAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TEN3=1.0D+03)
C
C     ----- -CPU- AND -WALL- TIME -----
C
C     DIAL(1) = DATE
C     DIAL(2) = TIME
C     DIAL(3) = CPUSEC
C     DIAL(4) = TOTSEC
C     DIAL(5) = DELCPU
C     DIAL(6) = DELTOT
C
      DIMENSION DIAL(6)
*     REAL*4 CPUTIM,ETIME_
*     TYPE TB_TYPE
*        SEQUENCE
*        REAL*4 USRTIM
*        REAL*4 SYSTIM
*     END TYPE
c     TYPE (TB_TYPE) ETIME_STRUCT
c     CPUTIM =ETIME_(ETIME_STRUCT)
c     DIAL(3)=CPUTIM
*     TOTTIM =TIMEF()/TEN3    
*     DIAL(4)=TOTTIM          
      RETURN
      END
C
C
C     SUBROUTINE TAKEN FROM MP2.F IN HONDO ------
C
      SUBROUTINE HND_DABMP(II,JJ,KK,LL,DM2,C,PHF,PMP2,DAB,
     1 DABMAX,
     1                 NBF2,NBF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (ZERO=0.0D0, FOUR=4.0D0, PT5=0.5D0)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      LOGICAL OUT
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_NSHEL /
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_SHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),
     1 LGXYZ(4,35)
      COMMON/HND_SHLNRM/PNRM(35)
      COMMON/HND_SHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_MP2DM2/NOC1,NOC2,ISTEP
      COMMON/HND_IJPAIR/IA(MXNBF)
      DIMENSION DM2(NBF2,NBF,*),C(NBF,*),PHF(NBF2),PMP2(NBF2)
      DIMENSION DAB(*)

      MINI = KMIN(II)
      MAXI = KMAX(II)
      MINJ = KMIN(JJ)
      MINK = KMIN(KK)
      MINL = KMIN(LL)
      LOCI = KLOC(II) - MINI
      LOCJ = KLOC(JJ) - MINJ
      LOCK = KLOC(KK) - MINK
      LOCL = KLOC(LL) - MINL

      DABMAX = ZERO

      DO 400 I=MINI,MAXI

      PI = PNRM(I)
      MAXJ = KMAX(JJ)
      IF(IIEQJJ) MAXJ = I

      DO 300 J=MINJ,MAXJ

      PIJ = PI*PNRM(J)
      MAXK = KMAX(KK)
      IF(IJEQKL) MAXK = I

      DO 200 K=MINK,MAXK

      PIJK = PIJ*PNRM(K)
      MAXL = KMAX(LL)
      IF(KKEQLL           ) MAXL = K
      IF(IJEQKL.AND.K.EQ.I) MAXL = J

      DO 100 L=MINL,MAXL

      PIJKL = PIJK*PNRM(L)

         NLA  = LOCI+I
         NSI  = LOCJ+J
         NMU  = LOCK+K
         NNU  = LOCL+L
         LASI = IA(MAX(NLA,NSI)) + MIN(NLA,NSI)
         MUNU = IA(MAX(NMU,NNU)) + MIN(NMU,NNU)

         DFAC = ZERO

C     ----- TRANSFORM NON-SEPARABLE TERM   -----

         DO 10 MI=NOC1,NOC2
   10       DFAC = DFAC + DM2(LASI,NMU,MI-NOC1+1)*C(NNU,MI)
     1                  + DM2(LASI,NNU,MI-NOC1+1)*C(NMU,MI)
         DFAC = DFAC + DFAC

C     ----- ADD SEPARABLE TERM (IF NEEDED) -----

         IF(ISTEP.EQ.1) THEN
         MUSI = IA(MAX(NMU,NSI)) + MIN(NMU,NSI)
         LANU = IA(MAX(NLA,NNU)) + MIN(NLA,NNU)
         MULA = IA(MAX(NMU,NLA)) + MIN(NMU,NLA)
         NUSI = IA(MAX(NNU,NSI)) + MIN(NNU,NSI)

         DFAC = DFAC + FOUR*PHF (MUNU)*PHF(LASI)
     1               -      PHF (MUSI)*PHF(LANU)
     2               -      PHF (MULA)*PHF(NUSI)
     3               + FOUR*PMP2(MUNU)*PHF(LASI)
     4               + FOUR*PMP2(LASI)*PHF(MUNU)
     5               -      PMP2(MUSI)*PHF(LANU)
     6               -      PMP2(MULA)*PHF(NUSI)
     7               -      PMP2(NUSI)*PHF(MULA)
     8               -      PMP2(LANU)*PHF(MUSI)

         ENDIF

         IF(NMU.EQ.NNU  ) DFAC = DFAC*PT5
         IF(NLA.EQ.NSI  ) DFAC = DFAC*PT5
         IF(MUNU.EQ.LASI) DFAC = DFAC*PT5

         IF(ABS(DFAC).GT.DABMAX) DABMAX = ABS(DFAC)
         NN = IGXYZ(1,I) + JGXYZ(1,J) + KGXYZ(1,K) + LGXYZ(1,L)
         DAB(NN) = DFAC*PIJKL
         IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,I,J,K,L,NN,DFAC

  100 CONTINUE
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
      RETURN

 9999 FORMAT(' -DABMP- ',4I4,4I3,I5,E20.12)
      END
C
C     ------ SUBROUTINE TAKEN FROM DER.F IN HONDO ------
C
      SUBROUTINE HND_JKDEND(II)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_INTTYP/NPOPLE,NHONDO
      COMMON/HND_RESTAR/NREST
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_GRAD12/DE(3,MXATOM)
C
C     ----- GET RESTART DATA -----
C
      IF(II.EQ.NSHELL) GO TO 10
      NREST=4
      IST=II+1
      JST=1
      KST=1
      LST=1
      IHONDO=0
      IF(NPOPLE.EQ.0) IHONDO=1
      GO TO 20
   10 NREST=0
      IST=1
      JST=1
      KST=1
      LST=1
   20 CONTINUE
C
C     ----- SAVE GRADIENT + RESTART DATA -----
C
      CALL HND_WRTGRD(DE,NREST,IST,JST,KST,LST)
      IF(NREST.EQ.0) RETURN
C
C     ----- CHECK CPU TIME AND EXIT IF NEEDED -----
C
cd    CALL TEXIT(0,NREST)
cd    IF(TIM.LT.TIMLIM) RETURN
cd    WRITE(IW,9998) NREST,IHONDO
cd    WRITE(IW,9999)
cd    CALL TEXIT(1,NREST)
      RETURN
 9999 FORMAT(' ... WARNING ... THIS JOB MUST BE RESTARTED ... ')
 9998 FORMAT(' ... RESTART PARAMETERS ... ',/,
     1 ' IN -CNTRL- NREST = ',I2,' IN -INTGRL- NHONDO = ',I2)
      END
C
C
C     ------- SUBROUTINE TAKEN FROM DER.F IN HONDO ------
C
      SUBROUTINE HND_OEDWR(DCHRG,NIJ0,ADCHRG,DASI,DASJ,NIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ADCHRG(13,1),DASI(1),DASJ(1)
      DIMENSION  DCHRG(15,1)
C
      DO 20 IJ=1,NIJ
      DO 10 I=1,13
   10 DCHRG( I,IJ+NIJ0)=ADCHRG(I,IJ)
      DCHRG(14,IJ+NIJ0)=DASI(IJ)
      DCHRG(15,IJ+NIJ0)=DASJ(IJ)
   20 CONTINUE
C
      RETURN
      END
C
C
C     ------ SUBROUTINE TAKEN FROM CTL.F IN HONDO ------
C
c*    SUBROUTINE HND_RDGRD(G,N1,N2,N3,N4,N5)
c*    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c*    PARAMETER (MXIODA=255)
c*    PARAMETER (MXATOM=500)
c*    PARAMETER (MXCART=3*MXATOM)
c*    PARAMETER (MXTEMP=MXCART+5)
c*    COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
c*    COMMON/HND_MACHIN/ISINGL,NBITS
c*    DIMENSION G(MXCART),X(MXTEMP),IX(1)
c*    EQUIVALENCE (IX(1),X(1))
C
c*    NX=MXTEMP
c*    CALL HND_DAREAD(IDAF,IODA,X,NX,3)
c*    DO I=1,MXCART
c*       G(I)=X(I)
c*    ENDDO
c*    I0=MXCART*ISINGL

c*    N1=IX(1+I0)
c*    N2=IX(2+I0)
c*    N3=IX(3+I0)
c*    N4=IX(4+I0)
c*    N5=IX(5+I0)
c*    RETURN
c*    END
C
C
C     ------ SUBROUTINE TAKEN FROM CTL.F IN HONDO ------
C
      SUBROUTINE HND_READPK(IS,XX,YY,NXX,NH,IERR,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-READPK-'/
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
      SUBROUTINE HND_TIMIT(IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      CALL HND_SECOND(TIM,WTIM)
       TX= TIM- TI
       TI= TIM
      WTX=WTIM-WTI
      WTI=WTIM
      IF(IFLAG.NE.0.and.ga_nodeid().eq.0)
     &WRITE(IW,9999) TX,TIM,WTX,WTIM
      RETURN
 9999 FORMAT(10X,'ELAPSED - CPU- TIME = ',F10.3,
     1        5X,'TOTAL - CPU- TIME = ',F10.3,/,
     2       10X,'ELAPSED -WALL- TIME = ',F10.3,
     3        5X,'TOTAL -WALL- TIME = ',F10.3)
      END
C

C
C     ------ MODIFIED FROM CTL.F IN HONDO ------
C
      SUBROUTINE HND_HONDO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_MACHIN/ISINGL,NBITS
C
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
C
      DATA MAXBIT /64/
      ISINGL=2
      NBITS =MAXBIT/ISINGL
C
      NSETO = 0
C
      RETURN
      END
C
C     ------ SUBROUTINE TAKEN FROM AIX.F IN HONDO -----
C
      SUBROUTINE HND_DAOPEN(IDAF,IODA,NAV,NDAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- OPEN A DIRECT ACCESS FILE -----
C        IDAF ..... FORTRAN FILE NUMBER
C        IODA ..... INDEX ARRAY OF LENGTH -NDAR-
C        NAV  ..... ASSOCIATED VARIABLE
C        NDAR ..... NUMBER OF -DA- INDEX WORDS
C                   IT MUST BE EQUAL TO 1 + NUMBER
C                   OF LOGICAL RECORDS TO BE
C                   WRITTEN ON -IDAF-.
C
      PARAMETER (MXFIL= 30)
      CHARACTER*80 DSNAME
      CHARACTER*6  STAT
      COMMON/HND_FILDEF/DSNAME(MXFIL),STAT(MXFIL)
      COMMON/HND_DAFNAV/IXDDAF(2048),NAV10,NAV20
      DIMENSION IODA(NDAR)
      LDAR=2048
      DO 10 I=1,NDAR
   10 IODA(I)=0
C
      IF(IDAF.NE.10) GO TO 20
      NAV  =1
      NAV10=1
      OPEN(UNIT=IDAF,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=8192,FORM='UNFORMATTED',FILE=DSNAME(IDAF))
      RETURN
   20 CONTINUE
      NAV  =1
      NAV20=1
      OPEN(UNIT=IDAF,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=8192,FORM='UNFORMATTED',FILE=DSNAME(IDAF))
      RETURN
      END
C
C
      SUBROUTINE HND_SPRTR(D,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      DIMENSION D(1),DD(10)
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
C     II = MAX0( I, J)
C     JJ = MIN0( I, J)
C     IJ = (II*(II-1))/2 + JJ
      IJ = (J-1)*N + I
  140 DD(K) = D(IJ)
      IF(LIST.EQ.0) WRITE (IW,9048) J,(DD(I),I = 1,K)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(DD(I),I = 1,K)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(DD(I),I = 1,K)
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
      END
C---------------------------------------------------------------
      subroutine hnd_cphf_init(rtdb)
      implicit double precision (a-h,o-z)
      integer rtdb
#include "mafdecls.fh"
#include "rtdb.fh"
      parameter (mxatom=500)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_cphf_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
      data zero    /0.0d+00/

c     ----- get nocc and nvir -----

      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen)
c     if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
c    $     nopen = 0
c
        nocc = nclosed 
        nvir = num - nocc 

c     ----- define parameters -----

      nc   = nat*3
      nxyz = nat*3
      norb = num
      nrot = nocc * nvir
c
      nc3  = nc * nc
      num2 =(num*(num+1))/2
      num3 = num* num
      norb2= norb*norb
      nocc2= nocc*nocc
c
c     ----- initialize -fcm- -----
c           -fcm- at x(i00)
      i00=1
      i10=i00+ nc3
      need1=i10-1
c
c     ----- transform -df/dq- to mo basis -----
c           -wax - at x(i10)
c           -aa  - at x(i20)
c           -ta  - at x(i30)
c           -c   - at x(i40)
c           -dfao- at x(i50)
c           -dfmo- at x(i60)
      i00=1
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ num3
c
      i60=i50+ num3
      i70=i60+ num3
      i80=i70+ num3
      need2 =i80-1 
c
c     ----- transform - ta = -1/2 * ds/da - to mo basis -----
c           -wax - at x(i10)
c           -ta  - at x(i30)
c           -c   - at x(i40)
c           -e   - at x(i50)
c           -dsmo- at x(i60)
c           -dsao- at x(i70)
c           -ds  - at x(i90)
      i00=1
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ num3
      i60=i50+ num
      i70=i60+ num3
      i80=i70+ num3
      need3 = i80 -1
c
c     ----- form modified fock matrices -----
c           -wax- at x(i10)
c           -aa - at x(i20)
c           -ta - at x(i30) upon entry
c           -fmo- at x(i40)
c           -xx - at x(i50) 
c           -ix - at x(i60)
      i00=1
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ norb2*nxyz
      i60=i50+ norb2
      i70=i60+ num3
      need4 = i70 -1
c
      need = max(need1,need2,need3,need4)

c     ----- memory block and pointers -----

      if(.not.ma_push_get(mt_dbl,need , 'mem init',i_init,init))
     & call errquit('hnd_cphf_init, malloc of init  failed',911)

c     ----- zero out the elements of the block -----

      do 10 i=1,need
   10 x(i+init-1)=zero
c
      return
      end
c
      subroutine hnd_cphf_popm 
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
      common/hnd_cphf_mem/i_init,init

      if(.not.ma_pop_stack(i_init))
     & call errquit('hnd_cphf_init, ma_pop_stack of init
     & failed',911)

      return
      end
c
      subroutine hnd_cphfx(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
      integer  rtdb,geom,basis
      logical out
      out = .false.
c
      if (out.and.ga_nodeid().eq.0 ) then
          print*,' ..... in subroutine cphfx .....'
      endif

c     ----- cphf integrated with nwchem -----

      call hnd_cphf_init(rtdb)
      call hnd_cphf(rtdb,geom,basis)
      call hnd_cphf_popm
c
      return
      end
c
      subroutine hnd_cphf(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
      parameter (mxatom=500)
      parameter(mxnbf=2048)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "schwarz.fh"
#include "cfock.fh"
c#include "pstat.fh"
      integer  rtdb,geom,basis
c*    integer  g_fc
      integer  g_aodf(3*mxatom), g_modf(3*mxatom)
      integer  g_aods(3*mxatom), g_mods(3*mxatom)
      integer  g_dsta(3*mxatom), g_tmpta(3*mxatom)
      integer  g_modfwa(3*mxatom)
c*    integer  g_vecs,g_eval
c*    integer  g_eval
      integer  g_newaods,g_tmpaods
      integer  g_mafock(1),g_mmfock(1),g_tmpfock(1)
      double precision occ(mxnbf)
      double precision jfac(1),kfac(1)
c*    double precision  occ,evals
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read
      external movecs_read
      character*255 movecs_in,movecs_out
      character*255 movecs
      logical  status
      logical  oskel
      double precision tol2e
      character*256 cphf_rhs
      integer  g_wax,g_wax_tr
      integer  g_rhs
      logical  cphf
      external cphf
      logical  file_write_ga,file_read_ga
      external file_write_ga,file_read_ga
c*    character*(nw_max_path_len) cphf_sol
      character*256 cphf_sol
c*      character*256 fname
      integer  g_ya, g_ya_tr
c
      character*8 errmsg
      logical out
      logical out_grad
      common/hnd_iofile/ir,iw,ip
c*    common/hnd_iofile/ir,iw
      common/hnd_dafile/idaf,nav,ioda(255)
      common/hnd_grad12/de(3,mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_cphf_mem/i_init,init
      dimension errmsg(3)
c
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      data errmsg /'program ','stop in ','-  cphf-'/
      data zero   /0.0d+00/
      data two   /2.0d+00/ 
      data four  /4.0d+00/ 
      data eight /8.0d+00/
      data pt5   /0.5d+00/
      data nft17  /17/
      data nft18  /18/
c
      out      = .false.
      out_grad = .true.

c     ----- get nocc and nvir -----

      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen)

c*    if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
c*   $     nopen = 0
c*
      nocc = nclosed 
      nvir = num - nocc 
      na   = nocc

c     ----- define parameters -----

      nc   = nat*3
      nxyz = nat*3
      norb = num
      nrot = nocc * nvir
c*    if(nvir .eq. 0) return
      nc3  = nc * nc
      num2 =(num*(num+1))/2
      num3 = num* num
      norb2= norb*norb
      nocc2= nocc*nocc

c     ----- initialize -fcm- -----
c           -fcm- at x(i00)

      i00=init
      i10=i00+ nc3
      need=i10-1

c     ----- load fcm from disk to local memory -----

      call hnd_rewfil(nft17)
      do 70 in=1,nc
   70 read(nft17)
      call hnd_hssrd(nft17,x(i00),nc3)

      if(out.and.ga_nodeid().eq.0) write(iw,9993)

c     ----- output force constant matrix -----

      if(out.and.ga_nodeid().eq.0) call hnd_cpfcm(x(i00),nat)

c     ----- transform -df/dq- to mo basis -----
c           -wax - at x(i10)
c           -aa  - at x(i20)
c           -ta  - at x(i30)
c           -c   - at x(i40)
c           -dfao- at x(i50)
c           -dfmo- at x(i60)

      i00=init
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ num3
c
      i60=i50+ num3
      i70=i60+ num3
      i80=i70+ num3
c
      need=i80-1

c     ----- get mo vectors from file -----

       g_vecs=ga_create_atom_blocked (geom, basis, 'mo vectors')

c     ----- figure output input/output MO vectors -----

       if (.not.rtdb_cget(rtdb,'scf:input  vectors',1,movecs_in))
     $      movecs_in = 'atomic'
       if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     $      movecs_out = ' '
       if (movecs.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
            movecs_out = 'movecs'
         else
            movecs_out = movecs_in
         endif
       endif
       movecs = movecs_out

c      ----- read eig-vec and put eng-val in x(i40) -----

       status = movecs_read (movecs, 1, occ, x(i40), g_vecs)

c     ----- rewind nft17 and load into x(i70) -----

      call hnd_rewfil(nft17)

      do in =1,nc
c
        call hnd_hssrd(nft17,x(i70),num2) 
        do 450 j = 1,num
         do 450 i = 1,num

          ii = max0(i,j)
          jj = min0(i,j)
          ij = (ii*(ii-1))/2 + jj

c     ----- expand from triangular matrix to squre matrix ----

  450   x(i50+num*(j-1)+i-1) = x(i70+ij-1)

        if (out.and.ga_nodeid().eq.0) then
            call hnd_sprtr(x(i50),num)
        endif

c     ----- ao-df shell order from hondo to nwchem:x(i50)->x(i70)--

        call hnd_hndnw_tran(x(i50),x(i70),num)

c     ------ create ga block -----

        g_aodf(in) =ga_create_atom_blocked ( geom, basis, 'ao df' )
        g_modf(in) =ga_create_atom_blocked ( geom, basis, 'mo df' )
        call ga_zero(g_aodf(in))
        call ga_zero(g_modf(in))

c     ----- copy from local memory to global array -----

        if (ga_nodeid().eq.0) then
          call ga_put(g_aodf(in),1,num,1,num,x(i70),num)
        endif

c     ----- create MO DF matrix  -----

        if (out) call ga_print(g_aodf(in))
        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs,       1, num,   1, num,
     $                         g_aodf(in),   1, num,   1, num,
     $                         g_modf(in),   1, num,   1, num)
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modf(in),   1, num,  1, num,
     $                         g_vecs,       1, num,  1, num,
     $                         g_aodf(in),   1, num,  1, num )
c
        if (out) call ga_print(g_aodf(in))
        call ga_copy(g_aodf(in),g_modf(in))
        if (out) call ga_print(g_modf(in))
        call ga_symmetrize(g_modf(in))

c     ----- output to standard output (of hondo) ----

        if (out.and.ga_nodeid().eq.0) write(iw,9995) 

        call ga_get(g_modf(in),1,num,1,num, x(i60),num)
        if (out.and.in.eq.nc.and.ga_nodeid().eq.0) then
            print*,'----------     modf   -------------'
            call hnd_sprtr(x(i60),num)
        endif

c     ----- destroy g_aodf(in) -----

        status = ga_destroy (g_aodf(in))

c     ----- modf contribution to -W_x^a- at x(i10)-----

       nbov = nocc*nvir
       g_modfwa(in) =ga_create_atom_blocked 
     &                 (geom, basis, 'modfwa')
       call ga_zero(g_modfwa(in))

       call ga_copy(g_modf(in),g_modfwa(in))
       call ga_dscal(g_modfwa(in),four)
       call ga_get(g_modfwa(in),1,nocc,nocc+1,num,
     &             x(i10+(in-1)*nbov),nocc)

c     ----- modf contribution to -A^a-   at x(i20)-----

       nboo = nocc*nocc
       call ga_get(g_modfwa(in),1,nocc,1,nocc,
     &             x(i20+(in-1)*nboo),nocc)
   
c     ----- release memory pointer of modf + modfwa --

       status = ga_destroy ( g_modf(in))
       status = ga_destroy ( g_modfwa(in))
c
      enddo
c
      if(out.and.ga_nodeid().eq.0) then
         do 110 ixyz=1,nxyz
            write(iw,9977) ixyz
            call hnd_prsq(x(i10+nrot *(ixyz-1)),nvir,nocc,nocc)
            write(iw,9976) ixyz
  110       call hnd_sprtr(x(i20+nocc2*(ixyz-1)),nocc)
      endif

c
c     ----- transform - ta = -1/2 * ds/da - to mo basis -----
c           -wax - at x(i10)
c           -ta  - at x(i30)
c           -c   - at x(i40)
c           -e   - at x(i50)
c           -dsmo- at x(i60)
c           -dsao- at x(i70)
c           -ds  - at x(i90)
c
      i00=init
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ num3
      i60=i50+ num
      i70=i60+ num3
      i80=i70+ num3
      need=i80-1

c      ----- read eig-vec and put eng-val in x(i50) -----

       status = movecs_read (movecs, 1, occ, x(i50), g_vecs)

       if (.not.status) then
         call errquit('hnd_dderjk: could not read mo vectors', 110)
       endif

c     ----- read ds from disk to local memory x(i70) -----

      rewind nft18

      do in=1,nc

c     ----- load ds/dq from disk to local memory x(i70) -----  

       call hnd_hssrd(nft18,x(i70),num2)
       do j = 1,num
        do i = 1,num
         ii = max0(i,j)
         jj = min0(i,j)
         ij = (ii*(ii-1))/2 + jj

c     ----- expand from triangular matrix to squre matrix ----

         x(i60+num*(j-1)+i-1) = x(i70+ij-1)
c*       x(i60+num*(j-1)+i-1) = x(i70+(in-1)*num2+ij-1)
        enddo
       enddo

c     ------ from shell order of hondo to nwchem -----

       call hnd_hndnw_tran(x(i60),x(i70),num)     

c     ------ create ga block -----

       g_aods(in) =ga_create_atom_blocked (geom, basis,'ao ds')
       g_mods(in) =ga_create_atom_blocked (geom, basis,'mo ds' )
       call ga_zero(g_aods(in))
       call ga_zero(g_mods(in))

c     ----- copy from local memory to global array -----

       if (ga_nodeid().eq.0) then
           call ga_put(g_aods(in),1,num,1,num,x(i70),num)
       endif
       
       call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs,       1, num,   1, num,
     $                         g_aods(in),   1, num,   1, num,
     $                         g_mods(in),   1, num,   1, num)
       call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_mods(in),   1, num,  1, num,
     $                         g_vecs,       1, num,  1, num,
     $                         g_aods(in),   1, num,  1, num )
c
       call ga_copy(g_aods(in),g_mods(in))

       call ga_symmetrize(g_mods(in))

c     ----- set -T^a- -----

       call ga_dscal(g_mods(in),-pt5)
       g_dsta(in) =ga_create_atom_blocked (geom, basis, 'ds ta' )
       call ga_copy(g_mods(in),g_dsta(in))

c      ----- -T^a- for -A^a-   ----- 

       call ga_get(g_dsta(in),1,nocc,1,nocc,
     &             x(i30+(in-1)*nboo),nocc)

c     ----- g_dsta contribution to -W_x^a- -----

c%% scale by 8.0 ( not the same as in the paper ?)
       g_tmpta(in) = ga_create_atom_blocked (geom, basis,'tmp ta')
       call ga_zero(g_tmpta(in))       
       call ga_copy(g_dsta(in),g_tmpta(in))
       call ga_dscal(g_tmpta(in),eight)
       call ga_get(g_tmpta(in),1,nocc,nocc+1,num,x(i60),nocc)

c     ----- adding -e*T^a- to -W_x^a- at x(i10)-----

       do ivir = 1,nvir
        do iocc = 1,nocc
           x(i10+(in-1)*nbov+(ivir-1)*nocc+iocc-1) =
     &     x(i10+(in-1)*nbov+(ivir-1)*nocc+iocc-1) +
     &     x(i60+(ivir-1)*nocc+iocc-1)*x(i50+iocc-1)
        enddo
       enddo

c     ----- release memory pointer -----

         status = ga_destroy ( g_aods(in))
         status = ga_destroy ( g_tmpta(in))
         status = ga_destroy ( g_mods(in))

c     ----- dsta (in) for calculations below -----
c          status = ga_destroy ( g_dsta(in))
c
      enddo

      if(out.and.ga_nodeid().eq.0) then
         do 120 ixyz=1,nxyz
            write(iw,9975) ixyz
            call hnd_prsq(x(i10+nrot *(ixyz-1)),nvir,nocc,nocc)
            write(iw,9994) ixyz
  120       call hnd_sprtr(x(i30+nocc2*(ixyz-1)),nocc)
      endif

c     ----- add - ta * tb - contributions to -fcm- -----
c           -ta- at x(i30) upon entry
c           -e - at x(i50) upon entry

      call hnd_cpwab1(x(i00),x(i30),x(i50),ia,nocc,nocc2,nxyz)
      
      if(out.and.ga_nodeid().eq.0) then
             write(iw,9992)
             call hnd_sprtr(x(i00),nc)
      endif
c     ----- form modified fock matrices -----
c           -wax- at x(i10)
c           -aa - at x(i20)
c           -ta - at x(i30) upon entry
c           -fmo- at x(i40)
c           -xx - at x(i50) 
c           -ix - at x(i60) 

      i00=init
      i10=i00+ nc3
      i20=i10+ nrot *nxyz
      i30=i20+ nocc2*nxyz
      i40=i30+ nocc2*nxyz
      i50=i40+ norb2*nxyz
      i60=i50+ norb2
      i70=i60+ num3
      need=i70-1

c     ------ create ga block -----
c     ------ newaods == newao-dsta; tmpaods == newtmp-dsta -----

       g_newaods =ga_create_atom_blocked (geom, basis,'new ao-ds')
       g_tmpaods =ga_create_atom_blocked (geom, basis,'tmp ao-ds' )

       call ga_zero(g_newaods)
       call ga_zero(g_tmpaods)

c     ----- create memory blocks -----

      g_mafock(1)  =ga_create_atom_blocked (geom, basis,'ma  fock')
      g_mmfock(1)  =ga_create_atom_blocked (geom, basis,'mm  fock')
      g_tmpfock(1) =ga_create_atom_blocked (geom, basis,'tmp fock')

      do in = 1,nc

c     ----- form T^a in ao basis -----

        if (out) call ga_print(g_dsta(in))

        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_vecs,       1, num,   1, nocc,
     $                         g_dsta(in),   1, nocc,  1, nocc,
     $                         g_tmpaods ,   1, num,   1, nocc)
       call ga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                         g_tmpaods,    1, num,  1, nocc,
     $                         g_vecs,       1, nocc, 1, num,
     $                         g_newaods,    1, num,  1, num )


c     ----- parameters for modified 2e-fock matrix -----

       nfock   =  1

c     ----- jfac = 8.0; kfac = -4.0 for A^a -----
c     ----- jfac =16.0; kfac = -8.0 for W_x^a -----

       jfac(1) =  8.0
       kfac(1) = -4.0 

c     ----- get tol2e from scf_get_info.F -----
     
       if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
	   tol2e = 1d-9
       endif

c     ----- load symmetry information from scf_get_info.F -----

       If (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
          oskel = sym_number_ops(geom) .gt. 0
       endif

       call ga_copy(g_newaods,g_tmpfock(1))

       if (out) call ga_print(g_tmpfock(1)) 

       call ga_zero(g_mafock(1))

       call int_init(rtdb,1,basis)
       call schwarz_init(geom,basis)
       call scf_get_fock_param(rtdb,tol2e)
       call fock_2e(geom,basis,nfock,jfac,kfac,tol2e,oskel,
     &             g_tmpfock(1),g_mafock(1))
       call fock_2e_tidy(rtdb)
       call int_terminate()
       call schwarz_tidy()

       if (out) call ga_print(g_mafock(1))

c     ----- [occxocc] for A^a and [occxvir] for W_x^a -----

       call ga_zero(g_tmpfock(1))

       call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs,       1, num,    1, num,
     $                         g_mafock(1),  1, num,    1, num,
     $                         g_tmpfock(1), 1, num,    1, num )

       call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_tmpfock(1), 1, num,    1, num,
     $                         g_vecs,       1, num,    1, num,
     $                         g_mmfock(1) , 1, num,    1, num )

c      if (out) call ga_print(g_mmfock(1))

c     ----- download g_mmfock to local memory:x(i60) -----

       call ga_get(g_mmfock(1),1,num,1,num,
     &             x(i60),num)

c     ----- mod-fock contribution to -W_a^x- at x(i10) -----

       do ivir = 1,nvir
        do iocc = 1,nocc
           x(i10+(in-1)*nbov+(ivir-1)*nocc+iocc-1) =
     &     x(i10+(in-1)*nbov+(ivir-1)*nocc+iocc-1) +
     &     two*x(i60+(nocc+ivir-1)*num+iocc-1)
        enddo
       enddo

c     ----- mod-fock contribution to -A^a-   at x(i20) -----

       do jocc = 1,nocc
        do iocc = 1,nocc
           x(i20+(in-1)*nboo+(jocc-1)*nocc+iocc-1) =
     &     x(i20+(in-1)*nboo+(jocc-1)*nocc+iocc-1) +
     &     x(i60+(jocc-1)*num+iocc-1)
        enddo
       enddo
c
       status = ga_destroy ( g_dsta(in))
c
      enddo

c     ----- release memory pointer -----

      status = ga_destroy ( g_newaods )
      status = ga_destroy ( g_tmpaods )
      status = ga_destroy ( g_mafock(1))
      status = ga_destroy ( g_mmfock(1))
      status = ga_destroy ( g_tmpfock(1))
c
      if(out.and.ga_nodeid().eq.0) then
         do 130 ixyz=1,nxyz
            write(iw,9974) ixyz
            call hnd_prsq(x(i10+nrot *(ixyz-1)),nvir,nocc,nocc)
            write(iw,9973) ixyz
  130       call hnd_sprtr(x(i20+nocc2*(ixyz-1)),nocc)
      endif
c
c     ----- add - tb * aa - contributions to -fcm- -----
c           -fcm- at x(i00)
c           -wax- at x(i10)
c           -aa - at x(i20)
c           -tb - at x(i30)
c

      if (out.and.ga_nodeid().eq.0) call hnd_sprtr(x(i00),nc)
      if (out.and.ga_nodeid().eq.0) write(iw,9973) nxyz
      if (out.and.ga_nodeid().eq.0) 
     &    call hnd_sprtr(x(i20+nocc2*(nxyz-1)),nocc)
      if (out.and.ga_nodeid().eq.0) 
     &    call hnd_sprtr(x(i30+nocc2*(nxyz-1)),nocc)

      call hnd_cpwab2(x(i00),x(i20),x(i30),ia,nocc,nocc2,nxyz)

      if(out.and.ga_nodeid().eq.0) then
         write(iw,9991)
         call hnd_sprtr(x(i00),nc)
cd       CALL HND_CPFCM(X(I00),NAT)
      endif

c     ----- cphf solutions -----
c           -fcm    at x(i00)
c           -W_x^a- at x(i10)
c           -A^a-   at x(i20)
c           -Y^a -  at x(i30)

         i00 =init
         i10 =i00+ nc3
         i20 =i10+ nrot *nxyz
         i30 =i20+ nocc2*nxyz
         i40 =i30+ nrot *nxyz

c     ----- create ga blocks -----

      if(.not.ga_create(mt_dbl,nocc,nvir,'g_wax',1,1,g_wax))
     &     call errquit('could not alloc wax array',1)

      if(.not.ga_create(mt_dbl,nvir,nocc,'g_wax_tr',1,1,g_wax_tr))
     &     call errquit('could not alloc tran array',1)

      if(.not.ga_create(mt_dbl,nvir*nocc,1,'cphf_rhs:hnd_cphf',
     $    1,1,g_rhs))call errquit('hnd_cphf: cphf_rhs',0)

      if(.not.ga_create(mt_dbl,nvir,nocc,'g_ya',1,1,g_ya))
     &     call errquit('could not alloc ya array',1)

      if(.not.ga_create(mt_dbl,nocc,nvir,'g_ya_tr',1,1,g_ya_tr))
     &     call errquit('could not alloc ya_tr arrayy',1)

      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)

      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)

c     ----- copy W_x^a from local memory to the global array ----

      do in = 1, nc

       if (ga_nodeid().eq.0) then
          call ga_put(g_wax,1,nocc,1,nvir,x(i10+nbov*(in-1)),nocc)
       endif

       if (out) call ga_print(g_wax)

       call ga_transpose(g_wax,g_wax_tr)

c     ----- reshape W_x^a to g_rhs format:[(nvir*nocc)x1] -----

       call ga_copy_patch('n',g_wax_tr,1,nvir,1,nocc,
     &     g_rhs,1,nvir*nocc,1,1)

       if (out) call ga_print(g_wax_tr)

       if(.not.file_write_ga(cphf_rhs,g_rhs))
     &     call errquit('could not write cphf_rhs',1)

c     ----- call cphf solver ----

       call pstat_on(ps_cphf)
       if(.not.cphf(rtdb))call errquit('could not do cphf',1)
       call pstat_off(ps_cphf)

c     ----- the g_rhs now holds the solution -----

       if(.not.file_read_ga(cphf_sol,g_rhs))
     &    call errquit('could not read cphf_rhs',1)

       if (out) call ga_print(g_rhs)

c     ----- reshape the g_hrs -----

       call ga_copy_patch('n',g_rhs,1,nvir*nocc,1,1,
     &     g_ya,1,nvir,1,nocc)

       if (out) call ga_print(g_ya)

       call ga_transpose(g_ya,g_ya_tr)

       if (out) call ga_print(g_ya_tr)
    
c     ----- download g_ya_tr to local memory:x(i30) -----

       call ga_copy_patch('n',g_ya_tr,1,nocc,1,nvir,
     &     g_rhs,1,nocc*nvir,1,1)

       call ga_get(g_rhs,1,nocc*nvir,1,1,
     &             x(i30+(in-1)*nbov),nbov)

      enddo

c     ----- release memory pointer -----

      if(.not.ga_destroy(g_wax))call errquit(
     &     'could not destroy g_wax',1)

      if(.not.ga_destroy(g_wax_tr))call errquit(
     &     'could not destroy g_wax_tr handle',1)

      if(.not.ga_destroy(g_ya))call errquit(
     &     'could not destroy g_ya',1)

      if(.not.ga_destroy(g_ya_tr))call errquit(
     &     'could not destroy g_ya_tr',1)

      if(.not.ga_destroy(g_rhs))call errquit(
     &     'could not destroy g_ya_tr',1)

c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,9997) mxyz
         call hnd_prsq(x(i10),nxyz,nrot,nrot)
         write(iw,9996) mxyz
         call hnd_prsq(x(i30),nxyz,nrot,nrot)
      endif

c     ----- add - W_x^a * x^b - contributions to -fcm_ab -----

c     ----- create ga block -----

      if(.not.ga_create(mt_dbl,nbov,nc,'g_wax_tot',1,1,g_wax_tot))
     &     call errquit('could not alloc wax_tot array',1)
      if(.not.ga_create(mt_dbl,nbov,nc,'g_xa_tot' ,1,1,g_xa_tot ))
     &     call errquit('could not alloc xa_tot  array',1)

c     ----- copy from local memory to global array -----

      if (ga_nodeid().eq.0) then
       call ga_put(g_wax_tot,1,nbov,1,nc,x(i10),nbov)
       call ga_put(g_xa_tot, 1,nbov,1,nc,x(i30),nbov)
      endif

c     
      if (out) call ga_print(g_wax_tot)
      if (out) call ga_print(g_xa_tot )

c     ----- add -W_x^a*x^b- contribution to -fcm^ab -----

      do iy = 1,nc
         do ix = 1,nc
           x(i00+(iy-1)*nc+ix-1) =
     &     x(i00+(iy-1)*nc+ix-1) -
     & ga_ddot_patch(g_wax_tot,n,1,nbov,ix,ix,
     &               g_xa_tot ,n,1,nbov,iy,iy)

        enddo
      enddo

c     ----- output force constant matrix -----  
      if (out.and.ga_nodeid().eq.0) then
          write(iw,9990)
          call hnd_sprtr(x(i00),nc)
      endif  

c     ----- release memory pointer -----

      if(.not.ga_destroy(g_wax_tot))call errquit(
     &     'could not destroy g_wax_tot',1)

      if(.not.ga_destroy(g_xa_tot ))call errquit(
     &     'could not destroy g_xa_tot',1)

c     ----- for info only, print gradient -----

      if(out_grad.and.ga_nodeid().eq.0)
     &   call hnd_cpgrd(de,nat)

c     ----- symmetrixe -fcm- -----

c     if (ga_nodeid().eq.0) call hnd_cpfcm(x(i00),nat)
      if (ga_nodeid().eq.0) call hnd_hssprt(x(i00),nc)

c$$$c     ------ create ga block for force constant matrix -----
c$$$
c$$$      if(.not.ga_create(mt_dbl,nc,nc,'fcm',nc,nc,
c$$$     $     g_fc))call errquit('cphf: ga_create failed',0)
c$$$      call ga_zero(g_fc)
c$$$
c$$$c     ----- copy fc from local memory to global array -----
c$$$
c$$$      if (ga_nodeid().eq.0) then
c$$$          call ga_put(g_fc,1,nc,1,nc,x(i00),nc)
c$$$      endif
c$$$c
c$$$      call util_file_name('hess',.false.,.false.,fname)
c$$$      if (.not.file_write_ga(fname,g_fc))
c$$$     $   call errquit('rhf_hess: could not write scf_hess',1)
      call rhf_hessian_2file(nc,x(i00))

c     ----- fake dipole derivative -----

c     call hnd_mudset(nc,nc,0)

c
      if (out.and.ga_nodeid().eq.0) write(iw,9989)
c
      return
c
 9998 format(' - W_xy  - matrix ')
 9997 format(' - W_x^a - vector',i5)
 9996 format(' - Y^a   - vector ',i5)
 9995 format(/,' derivative fock matrices in mo basis',/,1x,24(1h-))
 9994 format(' - T^a   - matrix ',i5)
 9993 format('     -W^ab + S^ab - contribution to -fcm- ')
 9992 format(' add - T^a * T^b  - contribution to -fcm- ')
 9991 format(' add - T^a * A^b  - contribution to -fcm- ')
 9990 format(' add -W_x^a * X^b - contribution to -fcm- ')
 9989 format(/,' ...... end of -cphf- contribution',
     1         ' to energy grad. + hess. ...... ')

 9977 format(' after -dfmo- contribution to -W_x^a- ',i3)
 9976 format(' after -dfmo- contribution to -A^a  - ',i3)
 9975 format(' after -T^a - contribution to -W_x^a- ',i3)
 9974 format(' after -fck - contribution to -W_x^a- ',i3)
 9973 format(' after -fck - contribution to -W^a  - ',i3)
      end
c
      SUBROUTINE HND_CPWAB1(FCM,TA,E,IA,NOCC,NOCC2,NXYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FCM(NXYZ,1),TA(NOCC2,1),E(1),IA(1)
      DATA ZERO  /0.0D+00/
      DATA EIGHT /8.0D+00/
C
C     ----- MODIFIED FROM HONDO TO TAKE SQUARE MATRIX -----
C
      DO 20 IXYZ=1,NXYZ
      DO 20 JXYZ=1,NXYZ
      DUM=ZERO
      IJ = 0
      DO 10 IOCC=1,NOCC
      DO 10 JOCC=1,NOCC
      IJ = IJ + 1
      DUM=DUM+(E(IOCC)+E(JOCC))*TA(IJ,IXYZ)*TA(IJ,JXYZ)
   10 CONTINUE
      DUM=DUM*EIGHT
      FCM(IXYZ,JXYZ)=FCM(IXYZ,JXYZ)+DUM
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_CPWAB2(FCM,AA,TA,IA,NOCC,NOCC2,NXYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FCM(NXYZ,1),TA(NOCC2,1),AA(NOCC2,1),IA(1)
      DATA ZERO  /0.0D+00/
C
      DO 20 IXYZ=1,NXYZ
      DO 20 JXYZ=1,NXYZ
      DUM=ZERO
      IJ = 0
      DO 10 IOCC=1,NOCC
      DO 10 JOCC=1,NOCC
      IJ = IJ + 1
      DUM=DUM+TA(IJ,IXYZ)*AA(IJ,JXYZ)+TA(IJ,JXYZ)*AA(IJ,IXYZ)
   10 CONTINUE
      FCM(IXYZ,JXYZ)=FCM(IXYZ,JXYZ)+DUM
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_CPFCM(FCM,NAT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "global.fh"
      CHARACTER*4 DNAM
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      DIMENSION DNAM(3),FC(9),FCM(1)
      DATA DNAM /'d/dx','d/dy','d/dz'/
C
      NC=3*NAT
C
C     ----- OUTPUT FORCE CONSTANT MATRIX -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9998)

      MAX=0
   50 MIN=MAX+1
      MAX=MAX+3
      IF(MAX.GT.NAT) MAX=NAT
      if(ga_nodeid().eq.0) then
         WRITE(IW,9997)
         WRITE(IW,9994) (J,J=MIN,MAX)
         WRITE(IW,9997)
         WRITE(IW,9993) ((DNAM(I),I=1,3),J=MIN,MAX)
         WRITE(IW,9997)
      endif
      JMIN=3*(MIN-1)+1
      JMAX=3* MAX
      I=0
      DO 80 IAT=1,NAT
      DO 70 N=1,3
      I=I+1
      DO 60 J=JMIN,JMAX
      IJ=NC*(J-1)+I
   60 FC(J-JMIN+1)=FCM(IJ)
      KMAX=JMAX-JMIN+1
      if (ga_nodeid().eq.0) then
         WRITE(IW,9995) DNAM(N),(FC(J),J=1,KMAX)
      endif
   70 CONTINUE
   80 CONTINUE
      IF(MAX.LT.NAT) GO TO 50
      RETURN
 9998 FORMAT(/,' ---  Hessian --- ')
 9997 FORMAT(/)
 9996 FORMAT(5X,'Atom',9(6X,I2,4X))
 9995 FORMAT(5X,A4,9F12.7)
 9994 FORMAT(5X,'Atom',3(19X,I2,15X))
 9993 FORMAT(5X,'    ',9(6X,A4,2X))
      END
C
      SUBROUTINE HND_CPGRD(GRD,NAT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "global.fh"
      CHARACTER*4 DNAM
      common/hnd_iofile/ir,iw,ip
C     COMMON/HND_IOFILE/IR,IW
      DIMENSION DNAM(3),GRD(3,1)
      DATA DNAM /'d/dx','d/dy','d/dz'/
C
C     ----- OUTPUT GRADIENT -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+9
      IF(MAX.GT.NAT) MAX=NAT
      if (ga_nodeid().eq.0) then
         WRITE(IW,9997)
         WRITE(IW,9996) (I,I=MIN,MAX)
         WRITE(IW,9997)
         DO 20 N=1,3
   20    WRITE(IW,9995) DNAM(N),(GRD(N,I),I=MIN,MAX)
      endif
      IF(MAX.LT.NAT) GO TO 10
C
      RETURN
 9999 FORMAT(/,' --- Gradient --- ')
 9997 FORMAT(/)
 9996 FORMAT(5X,'Atom',9(6X,I2,4X))
 9995 FORMAT(5X,A4,9F12.7)
      END
C
C     ------ TAKEN FROM CTL.F in HONDO ------
C
      SUBROUTINE HND_PRTRL(D,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A TRIANGULAR MATRIX WITH LABELS
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXBFN=3072)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
C     COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION D(1),DD(10)
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
c%%
C     IF(LIST.EQ.0) WRITE (IW,9028) (BFLAB(I),I = IMIN,IMAX)
C     IF(LIST.EQ.1) WRITE (IW,9128) (BFLAB(I),I = IMIN,IMAX)
C     IF(LIST.EQ.2) WRITE (IW,9228) (BFLAB(I),I = IMIN,IMAX)
c%%
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
  140 DD(K) = D(IJ)
c%%
C     IF(LIST.EQ.0) WRITE (IW,9048) J,BFLAB(J),(DD(I),I = 1,K)
C     IF(LIST.EQ.1) WRITE (IW,9148) J,BFLAB(J),(DD(I),I = 1,K)
C     IF(LIST.EQ.2) WRITE (IW,9248) J,BFLAB(J),(DD(I),I = 1,K)
c%%
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(2X,A8,1X))
 9048 FORMAT(I5,2X,A8,10F11.5)
 9128 FORMAT(15X,7(4X,A8,3X))
 9148 FORMAT(I5,2X,A8,7F15.10)
 9228 FORMAT(15X,7(4X,A8,3X))
 9248 FORMAT(I5,2X,A8,7E15.8)
      END
      subroutine rhf_hessian_2file(n,hess)
      implicit none
#include "global.fh"
c:: writes file in vib_vib format using vib_vib filename default
c:: Note the default filename is set in task_freq
c:: filenames must be made identical.
c
c format of file is ascii lower triangular elements only.
c
c::passed
      integer n ! [input] the rank of the hessian (3*number of atoms)
      double precision hess (n,n) ! [input] the matrix
c::local
      integer h_unit
      parameter (h_unit=47)
      character*255 fname
      integer i,j
      if (ga_nodeid().eq.0) then
c     
c::   -- open default file
         call util_file_name('hess',.false.,.false.,fname)
         open(unit=h_unit,file=fname,form='formatted',status='unknown',
     $        err=99990,access='sequential')
         close(unit=h_unit,status='delete')
         open(unit=h_unit,file=fname,form='formatted',status='new',
     $        err=99991,access='sequential')
c::   -- compiler parinoia      (rak)
         rewind h_unit
c::   -- write information
         do i = 1,n
            do j = 1,i
               write(h_unit,10000) hess(i,j)
            enddo
         enddo
         close(unit=h_unit,status='keep')
      endif
      call ga_sync()
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname,'> as unknown file'
      call errquit('hnd_cphf: fatal error', 911)
99991 write(6,*)' could not open <',fname,'> as new file'
      call errquit('hnd_cphf: fatal error', 911)
      end
      subroutine hnd_get_dens(rtdb,dens,num)
c
c     only RHF(close shell) is considered 
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
      integer           mxatom
      parameter        (mxatom=500)
      integer           num
c
      integer          rtdb
      double precision dens(num*num) 
      integer          g_dens(3) ! max no. of density matrix 
      integer          g_wdens   
      logical          status
c 
c*    double precision c, zan
c*    integer          nat
c*    integer          num
c
      integer          i_init,init,need 
      integer          num3
      integer          i,j,ij
      integer          i00,i10,i20
c
      double precision x
      double precision zero
c
c     common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c     common/hnd_basnum/num
c
      logical debug
      dimension    x(1)
      equivalence (x(1),dbl_mb(1))
c
      data zero   /0.0d+00/
c
      debug=.false.
c
      do i = 1,3
         g_dens(i) = 0
      enddo 
c
      num3=num*num
c     
      need=num3+num3

c
c     allocate memory block
c
      if(.not.ma_push_get(mt_dbl,need , 'mem init',i_init,init))
     & call errquit('hnd_get_dens, malloc of init  failed',911)
c
      i00=init
      i10=i00+ num3
      i20=i10+ num3
c
      do 10 i=1,need
   10 x(i+init-1)=zero
c
c     get density matrices from eigen-v and occ. no. of nwchem
c
      call hnd_dens_path1(rtdb,g_dens,g_wdens)
c
      call ga_get(g_dens(1),1,num,1,num, x(i00),num)
c
c     tranform density matrix from nwchem to hondo shell order
c
      call hnd_nwhnd_tran(x(i00),x(i10),num)
c
c     transform to lower triangle matrix of hondo format
c

      ij = 0
      do i = 1,num
        do j = 1,i
           ij = ij + 1
           dens(1+ij-1) = x(i10 + (j-1)*(num)+(i-1))
        enddo
      enddo
c
      if (debug) call hnd_prtr(dens,num)
c
c     destry ga handles (created in ...) 
c
      do i=1, 3
        status = ga_destroy ( g_dens(i) )
      enddo
        status = ga_destroy ( g_wdens   )
c 
c     release memory block         
c
      if (.not.ma_pop_stack(i_init))
     & call errquit('hnd_get_dens, ma_pop_stack of init
     & failed',911)
c
      return
      end
c
      subroutine hnd_get_wdens(rtdb,wdens,num)
c
c     only RHF(close shell) is considered 
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
      integer           mxatom
      parameter        (mxatom=500)
      integer           num
c
      integer          rtdb
      double precision wdens(num*num) 
      integer          g_dens(3) ! max no. of density matrix 
      integer          g_wdens
      logical          status   
c 
c*    double precision c, zan
c*    integer          nat
c     integer          num
c
      integer          i_init,init,need
      integer          num3
      integer          i,j,ij
      integer          i00,i10,i20 
c
      double precision x
      double precision zero
c    
c     common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c     common/hnd_basnum/num
c
      dimension    x(1)
      equivalence (x(1),dbl_mb(1))
      logical debug
c
      data zero    /0.0d+00/
c
      debug=.false.
c
      do i = 1,3
         g_dens(i) = 0
      enddo 
c
      num3= num*num
c     
      need=num3+num3

c
c     allocate memory block
c
      if(.not.ma_push_get(mt_dbl,need , 'mem init',i_init,init))
     & call errquit('hnd_get_wdens, malloc of init  failed',911)
c
      i00=init
      i10=i00+ num3
      i20=i10+ num3
c
      do 10 i=1,need
   10 x(i+init-1)=zero
c
c     get density matrices from eigen-v and occ. no. of nwchem
c
      call hnd_dens_path1(rtdb,g_dens,g_wdens)
c
      call ga_get(g_wdens  ,1,num,1,num, x(i00),num)
c
c     tranform density matrix from nwchem to hondo shell order
c
      call hnd_nwhnd_tran(x(i00),x(i10),num)
c
c     transform to lower triangle matrix of hondo format
c

      ij = 0
      do i = 1,num
        do j = 1,i
           ij = ij + 1
           wdens(1+ij-1) = x(i10 + (j-1)*(num)+(i-1))
        enddo
      enddo
c
      if (debug) call hnd_prtr(wdens,num)
c
c     destry ga handles (created in ...) 
c
      do i=1, 3
        status = ga_destroy ( g_dens(i) )
      enddo
        status = ga_destroy ( g_wdens   )
c 
c     release memory block         
c
      if (.not.ma_pop_stack(i_init))
     & call errquit('hnd_get_wdens, ma_pop_stack of init
     & failed',911)
c
      return
      end
