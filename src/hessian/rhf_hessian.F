c**********************************************************************
c     subroutines interfaced with and adapted from HONDO
c $Id: rhf_hessian.F,v 1.35 2004-07-27 12:27:38 edo Exp $
c**********************************************************************
      subroutine hnd_der_init(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
c 
#include "hnd_tol.fh"
#include "nwc_const.fh"
      parameter (zero=0.0d+00)
      parameter (mxatom=nw_max_atom)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      parameter (mxnbf =2048)
      parameter (mxtype=5)
      logical  status
      integer  rtdb 
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
     $     ('hnd_der_init: geom_ncent  failed for geom', geom,
     &       GEOM_ERR)
      if(dbug) then
         write(iw,*) ' ... natoms = ', natoms
      endif
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hnd_der_init: bas_numcont failed for basis', basis,
     &       BASIS_ERR)

      if (.not. bas_numbf(basis, num)) call errquit
     $     ('hnd_der_init: bas_numbf failed for basis', basis,
     &       BASIS_ERR)
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
     $       ('hnd_der_init: bas_ce to cnr/bfr failed', iat, BASIS_ERR)
         if(dbug) then
            write(iw,*) ' ... iat, ishlo, ishhi = ',iat,ishlo,ishhi
         endif
c 
c     ----- loop over shells on atoms -----
c
         do ish=ishlo,ishhi
            status=bas_cn2bfr(basis,ish,ishbflo,ishbfhi)
            if(.not.status) call errquit
     $          ('hnd_der_init: bas_cn to bfr failed', ish, BASIS_ERR)
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
                  call errquit('hnd_der_init: illegal -ityp- ', ityp,
     &       BASIS_ERR)
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
     $    (' int_1e_ga_hnd_init: geom_cart_get failed.',911, GEOM_ERR)
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
#include "nwc_const.fh"
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
      parameter (mxatom=nw_max_atom)
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
c
      subroutine hnd_int_1e_stvdd(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "global.fh"
      integer rtdb,geom,basis

      integer ir,iw,ip
      common/hnd_iofile/ir,iw,ip
c
c     ----- routine adapted from HONDO's -STVDD- -----
c
      logical dbug
      logical out 
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... hnd_int_1e_stvdd ... '
      endif
c
c     ----- routine -stvdd- -----
c
      call hnd_stvdd(rtdb,geom,basis)
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
      out =out.or.dbug
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
      integer ir,iw,ip
      common/hnd_iofile/ir,iw,ip
      double precision total_wall,total_cpu
      logical ohondo
      logical dbug
c
      ohondo =.true.
      dbug = .false.
c
      if(ohondo) then

c     ----- fill appropriate arrays for HONDO integrals -----

         call hnd_der_init(rtdb,geom,basis)
         if(dbug) then
            write(iw,*) '-scf_hnd_dder-'
            write(iw,*) 'rtdb  = ',rtdb
            write(iw,*) 'geom  = ',geom
            write(iw,*) 'basis = ',basis
            call util_flush(6)
         endif
   
         call hnd_nw_hnd

c     ----- fill the common block of hnd_molabl

         call hnd_atlabl(rtdb,geom,basis)
       
c     ----- calculate 1e 2nd derivatives -----

         total_wall = -util_wallsec()
         total_cpu  = -util_cpusec()

         call hnd_int_1e_stvdd(rtdb,geom,basis)
       
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
C     ------ SUBROUTINES TAKEN FROM HSS.F IN HONDO ------
C 
      SUBROUTINE HND_STVDD(RTDB,GEOM,BASIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      INTEGER RTDB,GEOM,BASIS
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=nw_max_atom)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DAFILE/IDAF,NAV,IODA(255)
C
C     ----- MEMORY BLOCK   -----
C
      CALL HND_DDER1_INIT1(RTDB,GEOM,BASIS)
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
      CALL HND_DDER1_INIT2(RTDB,GEOM,BASIS)
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
      subroutine hnd_dder1_init1(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
      parameter (mxatom=nw_max_atom)
c
      integer     rtdb,geom,basis 
      character*8 scftyp
      logical     dbug
c 
      common/hnd_iofile/ir,iw,ip
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_dder1_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      if(dbug) then
         write(iw,*) 'in -hnd_dder1_init1-'
      endif
c
      nc=3*nat
      num2=(num*(num+1))/2
      num3= num*num
c     
      need1 = nc*nc+num2*nc
     &        +num2+num2+num3+num3

c     ----- extra memory to hold temporary de, fc, and df -----

      need2 = need1+nc+nc*nc+num2*nc
c
      need  = need2
c
c     ----- memory pointers -----
c
      if(.not.ma_push_get(mt_dbl,need ,
     &    'mem init:rhf_hessian:hnd_dder1_init1:1'
     &    ,i_init,init)) then
         write(iw,*)'The totoal number of basis functions=', num 
         write(iw,*)'The stack memory for dder_1e =', need 
         call errquit('hnd_dder1_init, malloc of init  failed',911,
     &       MA_ERR)
      endif
c
      i00=init
      i01=i00+ nc  *nc
      i10=i01+ num2*nc
      i20=i10+ num2
      i30=i20+ num2
      i40=i30+ num3
      i50=i40+ num3
c
      call hnd_hess_get_dens(rtdb,geom,basis,
     1                       x(i10),x(i20),
     2                       x(i30),x(i40),
     3                       scftyp)
c
      return
      end
      subroutine hnd_hess_get_dens(rtdb,geom,basis,
     1                             d1hnd,d2hnd,   
     2                             d1nw,d2nw,   
     3                             scftyp)
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
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      double precision      d1nw,d2nw
      double precision      d1hnd,d2hnd
c
      integer       rtdb              ! [input] database handle
      integer       geom, basis       ! handles
      character*255 title
      character*255 movecs_in, movecs_out
      character*8   scftyp 
      dimension     d1nw(*),d2nw(*)
      dimension     d1hnd(*),d2hnd(*)
      integer       ir
      integer       iw
      integer       ip
      common /hnd_iofile/ ir,iw,ip
c
      logical       dbug
c
      dbug=.false.
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'hess:input vectors',  1, 
     $     movecs_in )) then
           call util_file_name('movecs', .false.,.false.,movecs_in )
      endif
      call util_file_name_resolve(movecs_in , .false.)
      if (.not. rtdb_cget(rtdb, 'hess:output vectors',  1,
     $     movecs_out)) then
           call util_file_name('movecs', .false.,.false.,movecs_out)
      endif
      call util_file_name_resolve(movecs_out, .false.)
      if(dbug) then
         write(iw,*) '-hnd_hess_get_dens-'
         write(iw,*) 'movecs_in  = ',movecs_in
         write(iw,*) 'movecs_out = ',movecs_out
      endif
c
      call hnd_hess_dens_read(rtdb,geom,basis,movecs_in,
     1                        d1nw,d2nw,
     2                        d1hnd,d2hnd,
     3                        scftyp)
c
      return
      end
      subroutine hnd_hess_dens_read(rtdb,geom,basis,movecs,
     1                              d1nw,d2nw,
     2                              d1hnd,d2hnd,
     3                              scftyp)
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      logical       status
      integer       rtdb, basis, geom
      character*(*) movecs
      character*4   scftype
      character*8   scftyp 
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_molec
  
      integer       max_sh_bf, max_at_bf 
      integer       nat, nbf, nopen , nclosed , nvirt 
      integer                 nopenu, nclosedu, nvirtu  
      integer       ndens, nbf2, i, iset
      character*255 title_molec, basis_name_molec
      character*20  scftype_molec
      integer       nbf_molec, nsets_molec, nmo_molec(2)
      integer                               nmo
c
      integer   l_eval, k_eval, l_occ, k_occ
c
      double precision d1nw(*),d2nw(*)
      double precision d1hnd(*),d2hnd(*)
c
      logical dbug
      dbug=.false.
c
c     get information about basis set
c
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('vecs_read: could not get natoms',0, GEOM_ERR)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('vecs_read: could not get nbf',0, BASIS_ERR)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_sh_bf',0, BASIS_ERR)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_at_bf',0, BASIS_ERR)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2))
     &    call errquit('vecs_read: failed to read header', 911,
     &       DISK_ERR)
      if(dbug) then
         write(6,9999)
     1     movecs
     2   ,title_molec
     3   ,basis_name_molec
     4   ,scftype_molec
     5   ,nbf_molec
     6   ,nsets_molec
     7   ,nmo_molec
      endif
c
      if (nsets_molec.eq.1) then
         scftype = 'RHF'
         scftyp  = 'RHF     '
      else
         scftype = 'UHF'
         scftyp  = 'UHF     '
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2  = 2 * nbf
        ndens = 3
        nmo   = max(nmo_molec(1),nmo_molec(2))
      else if ( scftype .eq. 'RHF' ) then
        nbf2  = nbf
        ndens = 1
        nmo   =     nmo_molec(1)
      endif

C     allocate and initialize global and local memory

C     occupation numbers
      status = MA_push_get(MT_DBL, nbf2,'occ. numbers', l_occ, 
     $                                                  k_occ)
      if (.not.status) 
     $     call errquit('property: could not allocate l_occ',911,
     &       MA_ERR)
C     eigenvalues
      status = MA_push_get(MT_DBL, nbf2,'MO eval', l_eval, k_eval)
      if (.not.status) 
     $     call errquit('property: could not allocate l_eval',911,
     &       MA_ERR)
      
      status = ga_create(MT_DBL, nbf_molec, nmo,' molec MOs',
     $                           nbf_molec, nmo,g_molec)
      if (.not.status)
     $     call errquit('property: GA for molecular MOs',911,
     &       GA_ERR)
c
c     read ( alpha ) orbitals
c
      iset=1
      status = movecs_read(movecs,iset,
     $         dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
      if (.not.status)
     $  call errquit('property: could not read molecular MOs',911,
     &       DISK_ERR)
       
      nclosed = 0
      nopen   = 0
      nvirt   = 0
      do i=1,nmo_molec(1)
         if(dbl_mb(k_occ+i-1).eq.2.0d0) then
            nclosed=nclosed+1
         elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
            nopen=nopen+1
         else
            nvirt=nvirt+1
         endif
      enddo
      if(dbug) then
         write(6,9998) nclosed , nopen , nvirt
      endif
c
c     read ( beta ) orbitals
c
      if( nsets_molec.gt.1 ) then
            iset=2
         status = movecs_read(movecs,iset,
     $            dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
         if (.not.status)
     $     call errquit('property: could not read molecular MOs',911,
     &       DISK_ERR)

         nclosedu= 0
         nopenu  = 0
         nvirtu  = 0
         do i=1,nmo_molec(2)
            if(dbl_mb(k_occ+i-1).eq.2.0d0) then
               nclosedu=nclosedu+1
            elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
               nopenu=nopenu+1
            else
               nvirtu=nvirtu+1
            endif
         enddo
         if(dbug) then
            write(6,9997) nclosedu, nopenu, nvirtu
         endif
      endif

      status = ga_destroy ( g_molec )

      call hnd_hess_dens_make (geom, basis, 
     $     dbl_mb(k_eval), dbl_mb(k_occ), nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c
c     free temporary arrays
c
      status = MA_pop_stack ( l_eval)
      status = MA_pop_stack ( l_occ )

      return
 9999 format(
     1   /,'movecs             = ',a24,
     2   /,'title_molec        = ',a24,
     3   /,'basis_name_molec   = ',a24,
     4   /,'scftype_molec      = ',a24,
     5   /,'nbf_molec          = ',i12,
     6   /,'nsets_molec        = ',i12,
     7   /,'nmo_molec          = ',i12,i12 
     8      )
 9998 format(' nclosed , nopen , nvirt  = ',3i10)
 9997 format(' nclosedu, nopenu, nvirtu = ',3i10)
      end
      subroutine hnd_hess_dens_make (geom, basis,  
     $     evals, occ, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c     read mo vecs 

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis
      integer nbf, nbf2, nclosed, nopen, nclosedu, nopenu
      double precision evals, occ
      double precision d1nw ,d2nw 
      double precision d1hnd,d2hnd
      character*4 scftype
      character*(*) movecs

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read, file_read_ga
      external movecs_read, file_read_ga

      dimension evals(nbf2), occ(nbf2)
      dimension d1nw(*), d2nw(*)
      dimension d1hnd(*), d2hnd(*)

C     local variables
      integer g_vecs, g_vecs2
      integer g_dens, g_dens2
      integer i, j, ij
      integer iopen
      logical status
c
      logical dbug
      dbug=.false.

      if(dbug) then
         write(6,*) 'in -hnd_hess_dens_make'
      endif

C  get MO vectors from file 

      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      g_dens = ga_create_atom_blocked ( geom, basis, 'mo density' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2 = ga_create_atom_blocked (geom, basis, 'uhf mo vectors')
        g_dens2 = ga_create_atom_blocked (geom, basis, 'uhf mo density')
      endif

      status = movecs_read(movecs,1,occ,evals,g_vecs)
      if (.not.status) then
        call errquit('hess_dens_make: could not read mo vectors', 110,
     &       DISK_ERR)
      endif

      if (scftype .eq. 'UHF') then
        status = movecs_read(movecs,2,occ(nbf+1),evals(nbf+1),g_vecs2)
        if (.not.status) then
        call errquit('hess_dens_make: could not read uhf vectors', 110,
     &       DISK_ERR)
        endif
      endif
      
      if(dbug) then
         write(6,*) 'passing by here..... do not stop for now'
      endif
c

      if (scftype .eq. 'RHF') then

        if(nclosed.eq.0) then
          call ga_zero(g_dens)
        else
          call ga_matmul_patch('n','t',2.0d+00,0.0d+00,
     1         g_vecs,1,nbf,1,nclosed,
     2         g_vecs,1,nclosed,1,nbf,
     3         g_dens,1,nbf,1,nbf)
        endif
        if(nopen.gt.0) then
          do iopen=1,nopen
            call ga_matmul_patch('n','t',occ(nclosed+iopen),1.0d+00,
     1           g_vecs,1,nbf,nclosed+iopen,nclosed+iopen,
     2           g_vecs,nclosed+iopen,nclosed+iopen,1,nbf,
     3           g_dens,1,nbf,1,nbf)
          enddo
        endif
        call ga_symmetrize(g_dens)
        call ga_get(g_dens,1,nbf,1,nbf,d1nw,nbf)
 
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d1hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha density'
           call hnd_prtr(d1hnd,nbf) 
        endif

      else  ! UHF

        call ga_zero(g_dens )
        call ga_zero(g_dens2)

        if(nopen.gt.0) then
          do iopen=1,nopen
            call ga_matmul_patch('n','t',occ(iopen),1.0d+00,
     1           g_vecs,1,nbf,iopen,iopen,
     2           g_vecs,iopen,iopen,1,nbf,
     3           g_dens,1,nbf,1,nbf)
          enddo
        endif
        if(nopenu.gt.0) then
          do iopen=1,nopenu
            call ga_matmul_patch('n','t',occ(nbf+iopen),1.0d+00,
     1           g_vecs2,1,nbf,iopen,iopen, 
     2           g_vecs2,iopen,iopen,1,nbf,
     3           g_dens2,1,nbf,1,nbf)
          enddo
        endif

        call ga_symmetrize(g_dens)
        call ga_get(g_dens,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha density'
           call hnd_prtr(d1hnd,nbf) 
        endif

c

        call ga_symmetrize(g_dens2)
        call ga_get(g_dens2,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'beta  density'
           call hnd_prtr(d2hnd,nbf) 
        endif

      endif

c

      status = ga_destroy ( g_vecs )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
      endif

      return
      end
      subroutine hnd_dder1_init2(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
      parameter (mxatom=nw_max_atom)
c
      integer rtdb,geom,basis
      logical dbug
      integer ir,iw,ip
      common/hnd_iofile/ir,iw,ip
c 
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_dder1_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
c
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      if(dbug) then
         write(iw,*) 'in -hnd_dder1_init2-'
      endif
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
      i40=i30+ num3
      i50=i40+ num3
c
      call hnd_hess_get_wdens(rtdb,geom,basis,
     1                        x(i10),x(i20),
     2                        x(i30),x(i40),
     3                        scftyp)
c
      return
      end
      subroutine hnd_hess_get_wdens(rtdb,geom,basis,
     1                             d1hnd,d2hnd,   
     2                             d1nw,d2nw,   
     3                             scftyp)
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
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      double precision      d1nw,d2nw
      double precision      d1hnd,d2hnd
c
      integer       rtdb              ! [input] database handle
      integer       geom, basis       ! handles
      character*255 title
      character*255 movecs_in, movecs_out
      character*8   scftyp 
      dimension     d1nw(*),d2nw(*)
      dimension     d1hnd(*),d2hnd(*)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'prop:input vectors',  1, 
     $     movecs_in )) then
           call util_file_name('movecs', .false.,.false.,movecs_in )
      endif
      call util_file_name_resolve(movecs_in , .false.)
      if (.not. rtdb_cget(rtdb, 'prop:output vectors',  1,
     $     movecs_out)) then
           call util_file_name('movecs', .false.,.false.,movecs_out)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c
      call hnd_hess_wdens_read(rtdb,geom,basis,movecs_in,
     1                        d1nw,d2nw,
     2                        d1hnd,d2hnd,
     3                        scftyp)
c
      return
      end
      subroutine hnd_hess_wdens_read(rtdb,geom,basis,movecs,
     1                              d1nw,d2nw,
     2                              d1hnd,d2hnd,
     3                              scftyp)
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      logical       status
      integer       rtdb, basis, geom
      character*(*) movecs
      character*4   scftype
      character*8   scftyp 
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_molec
  
      integer       max_sh_bf, max_at_bf 
      integer       nat, nbf, nopen , nclosed , nvirt 
      integer                 nopenu, nclosedu, nvirtu  
      integer       ndens, nbf2, i, iset
      character*255 title_molec, basis_name_molec
      character*20  scftype_molec
      integer       nbf_molec, nsets_molec, nmo_molec(2)
      integer                               nmo
c
      integer   l_eval, k_eval, l_occ, k_occ
c
      double precision d1nw(*),d2nw(*)
      double precision d1hnd(*),d2hnd(*)
     
      logical dbug
      integer ir,iw,ip
      common/hnd_iofile/ir,iw,ip
c
      dbug=.false.
      if(dbug) then
         write(iw,*) 'in -hnd_hess_wdens_read-'
      endif
c
c     get information about basis set
c
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('vecs_read: could not get natoms',0, GEOM_ERR)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('vecs_read: could not get nbf',0, BASIS_ERR)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_sh_bf',0, BASIS_ERR)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_at_bf',0, BASIS_ERR)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2))
     &    call errquit('vecs_read: failed to read header', 911,
     &       DISK_ERR)
      if(dbug) then
         write(6,9999)
     1     movecs
     2   ,title_molec
     3   ,basis_name_molec
     4   ,scftype_molec
     5   ,nbf_molec
     6   ,nsets_molec
     7   ,nmo_molec
      endif
c
      if (nsets_molec.eq.1) then
         scftype = 'RHF'
         scftyp  = 'RHF     '
      else
         scftype = 'UHF'
         scftyp  = 'UHF     '
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2  = 2 * nbf
        ndens = 3
        nmo   = max(nmo_molec(1),nmo_molec(2))
      else if ( scftype .eq. 'RHF' ) then
        nbf2  = nbf
        ndens = 1
        nmo   =     nmo_molec(1)
      endif

C     allocate and initialize global and local memory

C     occupation numbers
      status = MA_push_get(MT_DBL, nbf2,'occ. numbers', l_occ, 
     $                                                  k_occ)
      if (.not.status) 
     $     call errquit('property: could not allocate l_occ',911,
     &       MA_ERR)
C     eigenvalues
      status = MA_push_get(MT_DBL, nbf2,'MO eval', l_eval, k_eval)
      if (.not.status) 
     $     call errquit('property: could not allocate l_eval',911,
     &       MA_ERR)
      
      status = ga_create(MT_DBL, nbf_molec, nmo,' molec MOs',
     $                           nbf_molec, nmo,g_molec)
      if (.not.status)
     $     call errquit('property: GA for molecular MOs',911, GA_ERR)
c
c     read ( alpha ) orbitals
c
      iset=1
      status = movecs_read(movecs,iset,
     $         dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
      if (.not.status)
     $  call errquit('property: could not read molecular MOs',911,
     &       DISK_ERR)
       
      nclosed = 0
      nopen   = 0
      nvirt   = 0
      do i=1,nmo_molec(1)
         if(dbl_mb(k_occ+i-1).eq.2.0d0) then
            nclosed=nclosed+1
         elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
            nopen=nopen+1
         else
            nvirt=nvirt+1
         endif
      enddo
      if(dbug) then
         write(6,9998) nclosed , nopen , nvirt
      endif
c
c     read ( beta ) orbitals
c
      if( nsets_molec.gt.1 ) then
            iset=2
         status = movecs_read(movecs,iset,
     $            dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
         if (.not.status)
     $     call errquit('property: could not read molecular MOs',911,
     &       DISK_ERR)

         nclosedu= 0
         nopenu  = 0
         nvirtu  = 0
         do i=1,nmo_molec(2)
            if(dbl_mb(k_occ+i-1).eq.2.0d0) then
               nclosedu=nclosedu+1
            elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
               nopenu=nopenu+1
            else
               nvirtu=nvirtu+1
            endif
         enddo
         if(dbug) then
            write(6,9997) nclosedu, nopenu, nvirtu
         endif
      endif

      status = ga_destroy ( g_molec )

      call hnd_hess_wdens_make (geom, basis, 
     $     dbl_mb(k_eval), dbl_mb(k_occ), nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c
c     free temporary arrays
c
      status = MA_pop_stack ( l_eval)
      status = MA_pop_stack ( l_occ )

      return
 9999 format(
     1   /,'movecs             = ',a24,
     2   /,'title_molec        = ',a24,
     3   /,'basis_name_molec   = ',a24,
     4   /,'scftype_molec      = ',a24,
     5   /,'nbf_molec          = ',i12,
     6   /,'nsets_molec        = ',i12,
     7   /,'nmo_molec          = ',i12,i12 
     8      )
 9998 format(' nclosed , nopen , nvirt  = ',3i10)
 9997 format(' nclosedu, nopenu, nvirtu = ',3i10)
      end
      subroutine hnd_hess_wdens_make (geom, basis,  
     $     evals, occ, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c     read mo vecs 

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis
      integer nbf, nbf2, nclosed, nopen, nclosedu, nopenu
      double precision evals, occ
      double precision d1nw ,d2nw 
      double precision d1hnd,d2hnd
      character*4 scftype
      character*(*) movecs

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read, file_read_ga
      external movecs_read, file_read_ga

      dimension evals(nbf2), occ(nbf2)
      dimension d1nw(*), d2nw(*)
      dimension d1hnd(*), d2hnd(*)

C     local variables
      integer g_vecs, g_vecs2
      integer g_wdens, g_wdens2
      integer i, j, ij
      integer iclosed,iopen
      logical status
c
      logical dbug
      dbug=.false.

      if(dbug) then
         write(6,*) 'in -hnd_hess_wdens_make'
      endif

C  get MO vectors from file 

      g_vecs  = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      g_wdens = ga_create_atom_blocked ( geom, basis, 'mo wensity' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2  = ga_create_atom_blocked (geom, basis,'uhf mo vectors')
        g_wdens2 = ga_create_atom_blocked (geom, basis,'uhf mo wensity')
      endif

      status = movecs_read(movecs,1,occ,evals,g_vecs)
      if (.not.status) then
        call errquit('hess_wdens_make: could not read mo vectors', 110,
     &       DISK_ERR)
      endif

      if (scftype .eq. 'UHF') then
        status = movecs_read(movecs,2,occ(nbf+1),evals(nbf+1),g_vecs2)
        if (.not.status) then
        call errquit('hess_wdens_make: could not read uhf vectors', 110,
     &       DISK_ERR)
        endif
      endif
      
      if(dbug) then
         write(6,*) 'passing by here..... do not stop for now'
      endif
c

      if (scftype .eq. 'RHF') then

        call ga_zero(g_wdens)
        do iclosed=1,nclosed
          occ(iclosed)=occ(iclosed)*evals(iclosed)
          call ga_matmul_patch('n','t',occ(iclosed),1.0d+00,
     1         g_vecs,1,nbf,iclosed,iclosed,
     2         g_vecs,iclosed,iclosed,1,nbf,
     3         g_wdens,1,nbf,1,nbf)
        enddo
        if(nopen.gt.0) then
          do iopen=1,nopen
            occ(nclosed+iopen)=occ(nclosed+iopen)*evals(nclosed+iopen)
            call ga_matmul_patch('n','t',occ(nclosed+iopen),1.0d+00,
     1           g_vecs,1,nbf,nclosed+iopen,nclosed+iopen,
     2           g_vecs,nclosed+iopen,nclosed+iopen,1,nbf,
     3           g_wdens,1,nbf,1,nbf)
          enddo
        endif
        call ga_symmetrize(g_wdens)
        call ga_get(g_wdens,1,nbf,1,nbf,d1nw,nbf)
 
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d1hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha w-density'
           call hnd_prtr(d1hnd,nbf) 
        endif

      else  ! UHF

        call ga_zero(g_wdens )
        call ga_zero(g_wdens2)

        if(nopen.gt.0) then
          do iopen=1,nopen
            occ(iopen)=occ(iopen)*evals(iopen)
            call ga_matmul_patch('n','t',occ(iopen),1.0d+00,
     1           g_vecs,1,nbf,iopen,iopen,
     2           g_vecs,iopen,iopen,1,nbf,
     3           g_wdens,1,nbf,1,nbf)
          enddo
        endif
        if(nopenu.gt.0) then
          do iopen=1,nopenu
            occ(nbf+iopen)=occ(nbf+iopen)*evals(nbf+iopen)
            call ga_matmul_patch('n','t',occ(nbf+iopen),1.0d+00,
     1           g_vecs2,1,nbf,iopen,iopen, 
     2           g_vecs2,iopen,iopen,1,nbf,
     3           g_wdens2,1,nbf,1,nbf)
          enddo
        endif

        call ga_symmetrize(g_wdens)
        call ga_get(g_wdens,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha w-density'
           call hnd_prtr(d1hnd,nbf) 
        endif

c

        call ga_symmetrize(g_wdens2)
        call ga_get(g_wdens2,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'beta  w-density'
           call hnd_prtr(d2hnd,nbf) 
        endif

      endif

c

      status = ga_destroy ( g_vecs )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
      endif

      return
      end
      subroutine hnd_dens_path1(rtdb,g_dens,g_wdens)
      implicit none
#include "errquit.fh"
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
     $     ('gradients: no converged SCF wavefunction available', 0,
     &       RTDB_ERR)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('gradients: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('gradients: no ao basis ', 0, RTDB_ERR)

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
      call hnd_dens_path2(rtdb,geom,basis,movecs_out,g_dens,g_wdens)
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
     &    ('gradients:error destroying geom and basis handles',911,
     &       GEOM_ERR)
c
      return
      end
c
      subroutine hnd_dens_path2 ( rtdb, basis, geom, movecs,
     &                            g_dens,g_wdens)
c
      implicit none
#include "errquit.fh"

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
     & call errquit('path2: could not get natoms',0, GEOM_ERR)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('path2: could not get nbf',0, BASIS_ERR)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('path2: could not get max_sh_bf',0, BASIS_ERR)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('path2: could not get max_at_bf',0, BASIS_ERR)
 
C     get information about type of calculation
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen, RTDB_ERR)

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
*ga:1:0
      status = ga_create( MT_DBL, nocc, nocc,'eigen_diag', 0, nocc, 
     $                       g_eigen_diag )
      call ga_zero ( g_eigen_diag )

C     eigenvalues
      status = MA_alloc_get(MT_DBL, nbf2,'MO evals', l_evals, k_evals)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_evals',1, MA_ERR)
      
C     occupation numbers
      status = MA_alloc_get(MT_DBL, nbf2,'occ. numbers', l_occ, k_occ)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_occ',1, MA_ERR)

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
     $     call errquit('grad: could not destroy g_eigen_diag',1,
     &       GA_ERR)
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
#include "errquit.fh"
#include "mafdecls.fh"
      common/hnd_dder1_mem/i_init,init

      if(.not.ma_pop_stack(i_init))
     & call errquit('hnd_dder1_init, ma_pop_stack of init
     & failed',911, MA_ERR)

      return
      end
C
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
      PARAMETER (MXCART=3*MXATOM)
      LOGICAL OUT
C     LOGICAL SOME
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      DIMENSION DDE(1)
      EQUIVALENCE (DDE(1),DBL_MB(1))
      DIMENSION DNAM(3)
C
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ZERO   /0.0D+00/
      DATA THREE  /3.0D+00/
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION DNAM(3)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL SOME,OUT,DBUG
      DATA ZERO    /0.0D+00/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
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
      IF(DBUG) THEN
         WRITE(IW,*) 'HND-HLFX'
         CALL HND_PRTR(X(I10),NUM)
      ENDIF      
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
         DO I=1,3
            WRITE(IW,9992) DNAM(I),(DE(I,J),J=MIN,MAX)
         ENDDO
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
         DO IC=1,NC
            CALL HND_PRTR(X(I01+(IC-1)*NUM2),NUM)
         ENDDO
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
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
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
#include "global.fh"
C
      OUT =.FALSE.
      SOME=.FALSE.
      SOME=SOME.OR.OUT
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
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
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      DIMENSION DVIJ(3,225),DDVIJ(3,3,225)
      DIMENSION PIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,5,5),  YV(5,5,5),  ZV(5,5,5)
      DIMENSION  DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
      DIMENSION DDXV(5,5,5),DDYV(5,5,5),DDZV(5,5,5)
      DIMENSION W2(5),W4(5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HLFSPD-'/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
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
      TX= SQRT(TT)
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
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
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
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
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
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
      TX= SQRT(TT)
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL OUT
      DATA ZERO   /0.0D+00/
C
      OUT =.FALSE. 
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=nw_max_atom)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
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
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      OUT =.FALSE.
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
      DIMENSION PIJ(225)
      DIMENSION DIDVIJ(3,3,225)
      DIMENSION W2(5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XV(6,5,5),    YV(6,5,5),    ZV(6,5,5)
      DIMENSION    DXV(6,5,5),   DYV(6,5,5),   DZV(6,5,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION DDXVDI(5,5,5),DDYVDI(5,5,5),DDZVDI(5,5,5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- HLFD -'/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
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
      TX= SQRT(TT)
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION W2(5)
      DIMENSION FC(LFC,1),DF(LDF,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION     XV(6,5,5),    YV(6,5,5),    ZV(6,5,5)
      DIMENSION    DXV(6,5,5),   DYV(6,5,5),   DZV(6,5,5)
      DIMENSION  DXVDI(5,5,5), DYVDI(5,5,5), DZVDI(5,5,5)
      DIMENSION DDXVDI(5,5,5),DDYVDI(5,5,5),DDZVDI(5,5,5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- HLFD -'/
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
      TX= SQRT(TT)
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION X(1)
      LOGICAL OUT
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION DNAM(3)
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA ZERO   /0.0D+00/
      DATA NFT17  /17/
C
      OUT =.FALSE.
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
   10    DE(I,1)=ZERO
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
   12    X(I01+I-1)=ZERO     
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
      IF(OUT) THEN          
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
         CALL HND_HSSOUT(X(I00),NC,NC,NC)
         WRITE(IW,9990)
         DO 9250 IC=1,NC
 9250    CALL HND_PRTR(X(I01+(IC-1)*NUM2),NUM)
C
         WRITE(IW,9998)
      ENDIF
C
      ENDIF
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL OUT
      LOGICAL DBUG
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
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
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
      IF(DBUG) THEN       
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
      IF(OUT) THEN        
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
      ENDIF
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                             NI,NJ,CX,CY,CZ
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
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- TVD  -'/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
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
      TX= SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTMOD
      NJ=J
      DO 370 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XS(I,J)=XINT*TX
      YS(I,J)=YINT*TX
      ZS(I,J)=ZINT*TX
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
      TX= SQRT(TT)
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
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                               NI,NJ
     1                              ,CX,CY,CZ
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_GRAD12/DE(3,MXATOM)
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
c     DATA NFT17  /17/
      DATA ERRMSG /'PROGRAM ','STOP IN ','- TVD  -'/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
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
      TX= SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJTMOD
      NJ=J
      DO 370 I=1,LITDER
      NI=I
      CALL HND_DSXYZ
      XS(I,J)=XINT*TX
      YS(I,J)=YINT*TX
      ZS(I,J)=ZINT*TX
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
      TX= SQRT(TT)
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
C     COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_DDER1_MEM/I_INIT,INIT
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      LOGICAL OUT,SOME
      DIMENSION DNAM(3)
      DATA ZERO        /0.0D+00/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA NFT17  /17/
      DATA NFT18  /18/
C
      OUT =.FALSE.
      SOME=.FALSE.
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
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
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      DATA ONE    /1.0D+00/
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
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),
     1CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),
     2KMAX(MXSHEL),NSHELL
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
c     DATA NFT17  /17/
c     DATA NFT18  /18/
      DATA DNAM   /3HE'X,3HE'Y,3HE'Z/
      DATA RLN10  /2.30258D+00/
      DATA ONE    /1.0D+00/
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
C
      TOL=RLN10*ITOL
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
#include "nwc_const.fh"
      CHARACTER*4 CLAB
      character*8 anam,bflab
      character*2 bnam
      PARAMETER (NUMMX=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_IOFILE/IR,IW,IP
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

      SUBROUTINE HND_DVXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                              ,NI,NJ
     1                              ,CX,CY,CZ
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
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                              ,NI,NJ
     1                              ,CX,CY,CZ
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
      SUBROUTINE HND_JKDD(RTDB,GEOM,BASIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      INTEGER RTDB,GEOM,BASIS
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      subroutine hnd_dder2_init(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
      parameter (mxatom=nw_max_atom)
      integer rtdb,geom,basis 
c
      logical dbug
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
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      if(dbug) then
         write(iw,*) 'in -hnd_dder2_init-'
      endif
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

      needa = l2 + l2 + l3 + l3
      need  = needm + needa
   
c     ----- memory block to store 1-e results fc + df -----

      needc = nc*nc +nc*l2
      need  = need + needc
c
      if(.not.ma_push_get(mt_dbl,need,
     &    'mem init:rhf_hessian:hnd_dder2_init:2',
     &    i_init2,init2)) then
         write(iw,*)'The total number of basis functions =',num
         write(iw,*)'The stack memory for dder_2e =',need
         call errquit('hnd_dder2_init, malloc of init2  failed',911,
     &       MA_ERR)
      endif
c
      return
      end
c
      subroutine hnd_dder2_inita(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
      parameter (mxatom=nw_max_atom)
      integer rtdb,geom,basis 
      character*8 scftyp
      logical dbug
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
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      if(dbug) then
         write(6,*) 'in -hnd_dder2_inita-'
      endif
c
      nc= nat*3
      l1= num
      l2=(num*(num+1))/2
      l3= num*num

c     ----- memory block for getting density matrix -----

      needa = l2+l3+l3    

c     ----- total memory block -----

      need  = needm   + needa
c
      i00=iwfn
      i01=i00+ nc* nc
      i02=i01+ nc* l2
      i10=i02
c
      i20   = i10   + l2
      i30   = i20   + l2     
      i40   = i30   + l3
      i50   = i40   + l3
c
      call hnd_hess_get_dens(rtdb,geom,basis,
     1                       x(i10),x(i20),
     2                       x(i30),x(i40),
     3                       scftyp)
      if(dbug) then
         write(6,*) 'in -hnd_dder2_inita-'
         call hnd_prtr(x(i10),l1)
      endif
c
      return
      end
      subroutine hnd_dder2_popm 
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "mafdecls.fh"
      common/hnd_dder2_mem/i_init2,init2
c
      if(.not.ma_pop_stack(i_init2))
     & call errquit('hnd_dder2_popm, ma_pop_stack of init2
     & failed',911, MA_ERR)
c
      return
      end
C
      SUBROUTINE HND_DDERJK(RTDB,GEOM,BASIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
      INTEGER  RTDB,GEOM,BASIS
c
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXNBF =2048)
      CHARACTER*8 ERRMSG
      CHARACTER*8 RHF,UHF,GVB,SCF,TYPMC,TYPCI
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      SOME=.FALSE.
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
      CALL HND_DDER2_INIT(RTDB,GEOM,BASIS)
C
C     ----- RESCALE MEMORY POINTER ----- 
C
      IADDR=INIT2
      CALL HND_JKDMEM(MDER,IADDR,LENGTH)
C
C     ----- GET DENSITY MATRIX ----
C
      CALL HND_DDER2_INITA(RTDB,GEOM,BASIS)
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
      DO I=1,L2
         X(I+IDAB-1)=X(I+I10-1)
      ENDDO
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL OUT
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION FC(NC,NC),DF(LDF,NC)
      DATA NFT17 /17/
C
      ENTRY      HND_JDDCI(FC,NC,DF,LDF,NBF)
C
C     ----- READ DERIVATIVE FOCK MATRICES -----
C     ----- DIAGONAL IS MULTIPLIED BY 2.0 -----
C
      DBUG= .FALSE.
      OUT = .FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9997)
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
       IF(DBUG.AND.GA_NODEID().EQ.0) CALL HND_PRTR (DF(1,IC),NBF)
C
   20 CONTINUE
C
C     ----- READ 1E-HESSIAN -----
C
      CALL HND_HSSRD(NFT17,FC,NC*NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9998)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_HSSOUT(FC,NC,NC,NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9996)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRSQ(DE,NC/3,3,3)
C
      RETURN
 9999 FORMAT(/,' 1E-DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
 9998 FORMAT(/,' 1E-HESSIAN MATRIX',/,1X,14(1H-))
 9997 FORMAT(/,' -JDDHF- ')
 9996 FORMAT(/,' 1E-GRADIENT      ',/,1X,14(1H-))
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
      COMMON/HND_IOFILE/IR,IW,IP
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

c     CALL HND_SECOND(TIM,WTIM)
c     DHND=TIM-TIM0-DG80
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXNBF =2048)
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION FC(NC,NC),DF(LDF,NC)
      DATA NFT17 /17/
      DATA PT5   /0.5D+00/
C 
      DBUG= .FALSE.
      OUT = .FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- WRITE DERIVATIVE FOCK MATRICES  -----
C     ----- OFF-DIAGONAL ARE DIVIDED BY 2.0 -----
C
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9997)
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
      IF(DBUG) CALL HND_PRTR (DF(1,IC),NBF)
C
      CALL HND_HSSWRT(NFT17,DF(1,IC),LDF)
   20 CONTINUE
C
C     ----- WRITE 2E-HESSIAN -----
C
      CALL HND_HSSWRT(NFT17,FC,NC*NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9998)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_HSSOUT(FC,NC,NC,NC)
      IF(OUT.AND.GA_NODEID().EQ.0) WRITE(IW,9996)
      IF(OUT.AND.GA_NODEID().EQ.0) CALL HND_PRSQ(DE,NC/3,3,3)
      RETURN
 9999 FORMAT(/,' 2E-DERIVATIVE FOCK MATRICES',/,1X,24(1H-))
 9998 FORMAT(/,' 2E-HESSIAN MATRIX',/,1X,14(1H-))
 9997 FORMAT(/,' -JDDSAV- ')
 9996 FORMAT(/,' 2E-GRADIENT      ',/,1X,14(1H-))
C
      END
C
C     ------ SUBROUTINES ABOVE FROM HSS.F IN HONDO -------
C
      SUBROUTINE HND_JKDSET(MDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      COMMON/HND_SHLTOL/RTOL,DTOL
      COMMON/HND_SHLORD/MODSHL(512),INVSHL(512)
      COMMON/HND_DSHPRT/OUT
      COMMON/HND_DSHINT/IST,JST,KST,LST
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSHL/MAPSHL(MXSHEL,48),MAPCTR(MXATOM,48)
      COMMON/HND_ISOPAC/INDIN(48),INDOUT(24)
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION MI(48)
      DIMENSION ISO(512,24)
      EQUIVALENCE (ISO(1,1),MAPSHL(1,1))
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
c-    CALL HND_DAREAD(IDAF,IODA,ISO,NISO,5)
c-    DO 130 II=1,NSHELL
c-    DO 110 IT=1,NTWD
c-110 INDOUT(IT)=ISO(II,IT)
c-    CALL HND_ISOOUT(NT)
c-    DO 120 IT=1,NT
c-120 ISO(II,IT)=INDIN(IT)
c-130 CONTINUE
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
#include "hnd_tol.fh"
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL OUT
      LOGICAL SP
      LOGICAL NORM
      LOGICAL VECTOR
      PARAMETER (LENVEC=255)
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MACHNV/VECTOR
      COMMON/HND_MEMORY/MAXCOR,MAXSCM
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
      INEED =ILAST - IWFN
      IFIRST=IWFN
      IADDR =ILAST
      ILAST =ILAST - 1

      IF(OUT) WRITE(IW,9996) IWFN,ICLB,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW,
     6 IFIRST,ILAST,INEED,IADDR
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
     5 ' IABV  ',I8,' ICV   ',I8,' IRW   ',I8,/,
     6 ' IFIRST',I8,' ILAST ',I8,' INEED ',I8,' IADDR ',I8)
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
#include "nwc_const.fh"
C
C     ----- SYMMETRYZE GRADIENT VECTOR -----
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_RESTAR/NREST
      COMMON/HND_DSHINT/IST,JST,KST,LST
      COMMON/HND_IOFILE/IR,IW,IP
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      SOME=.FALSE.
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
#include "nwc_const.fh"
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
      COMMON/HND_IOFILE/IR,IW,IP
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
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
#include "hnd_rys.fh"
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      LOGICAL IEQJ0,KEQL0,SAME0,FIRST
      LOGICAL SAMTYP,SAMSYM
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      DATA ZERO,PT5,FOUR /0.0D+00,0.5D+00,4.0D+00/
      DATA PT25 /0.25D+00/
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
      COMMON/HND_IOFILE/IR,IW,IP
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
#include "hnd_rys.fh"
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
      MAXG=MAXXYZ/MROOTS
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
      YY=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
C
      NG=NG+1
      ABV(1,NG)=AA
      ABV(2,NG)=BB
      ABV(3,NG)=RHO
      ABV(4,NG)=Q4DBDA
      ABV(5,NG)=YY
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
      IF(MROOTS.EQ.1) GO TO 480
      IF(.NOT.SPI) GO TO 420
      DO 410 IROOT=2,MROOTS
      DO 410 IG=1,NUMG
      DIJSI(IG+NUMG*(IROOT-1))=DIJSI(IG)
  410 CONTINUE
  420 IF(.NOT.SPJ) GO TO 440
      DO 430 IROOT=2,MROOTS
      DO 430 IG=1,NUMG
      DIJSJ(IG+NUMG*(IROOT-1))=DIJSJ(IG)
  430 CONTINUE
  440 IF(.NOT.SPK) GO TO 460
      DO 450 IROOT=2,MROOTS
      DO 450 IG=1,NUMG
      DKLSK(IG+NUMG*(IROOT-1))=DKLSK(IG)
  450 CONTINUE
  460 IF(.NOT.SPL) GO TO 480
      DO 470 IROOT=2,MROOTS
      DO 470 IG=1,NUMG
      DKLSL(IG+NUMG*(IROOT-1))=DKLSL(IG)
  470 CONTINUE
  480 CONTINUE
C
      IF(SKIPI) GO TO 520
      DO 510 IRXYZ=2,MROOTS*3
      DO 510 IG=1,NUMG
      AAI(IG+NUMG*(IRXYZ-1))=AAI(IG)
  510 CONTINUE
  520 IF(SKIPJ) GO TO 540
      DO 530 IRXYZ=2,MROOTS*3
      DO 530 IG=1,NUMG
      AAJ(IG+NUMG*(IRXYZ-1))=AAJ(IG)
  530 CONTINUE
  540 IF(SKIPK) GO TO 560
      DO 550 IRXYZ=2,MROOTS*3
      DO 550 IG=1,NUMG
      BBK(IG+NUMG*(IRXYZ-1))=BBK(IG)
  550 CONTINUE
  560 IF(SKIPL) GO TO 580
      DO 570 IRXYZ=2,MROOTS*3
      DO 570 IG=1,NUMG
      BBL(IG+NUMG*(IRXYZ-1))=BBL(IG)
  570 CONTINUE
  580 CONTINUE
C
      LVEC=LVEC+NUMG*MROOTS*3
      NVEC=NVEC+1
C
C     ----- COMPUTE ROOTS AND WEIGHTS FOR QUADRATURE -----
C
      CALL HND_JKWRYS(RWV,ABV,NUMG)
C
C     ----- COMPUTE COEFFICIENTS FOR RECURSION FORMULAE -----
C
      CALL HND_JKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1            ABV,CV,RWV,NUMG,MROOTS)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 2 CENTERS, 2-D ) -----
C
      IF(NUMG*MROOTS*3.LT.MINVEC)
     1 CALL HND_JKGNMS(GNM,NUMG*MROOTS*3,NMAX,MMAX,
     2             B00,B01,B10,C00,D00,F00)
      IF(NUMG*MROOTS*3.GE.MINVEC)
     1 CALL HND_JKGNMV(GNM,NUMG*MROOTS*3,NMAX,MMAX,
     2             B00,B01,B10,C00,D00,F00)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 4 CENTERS, 2-D ) -----
C
      IF(NUMG*MROOTS*3.LT.MINVEC)
     1 CALL HND_JKXYZS(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,GNM,
     2 NUMG*MROOTS*3,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     3 EXPNDI,EXPNDK)
      IF(NUMG*MROOTS*3.GE.MINVEC)
     1 CALL HND_JKXYZV(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,GNM,
     2 NUMG*MROOTS*3,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     3 EXPNDI,EXPNDK)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS FOR DERIVATIVES -----
C
      IF(NUMG*MROOTS*3.LT.MINVEC)
     1 CALL HND_JDXYZS(GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,
     2 NUMG*MROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,LIT,LJT,LKT,LLT,
     3 AAI,AAJ,BBK,BBL,
     4 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
      IF(NUMG*MROOTS*3.GE.MINVEC)
     1 CALL HND_JDXYZV(GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,
     2 NUMG*MROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,LIT,LJT,LKT,LLT,
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
      IF(NUMG*MROOTS.LT.MINVEC)
     1 CALL HND_DSPDFS(NUMG,MROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IF(NUMG*MROOTS.GE.MINVEC)
     1 CALL HND_DSPDFV(NUMG,MROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
C
      ELSE
C
      IF(NUMG*MROOTS.LT.MINVEC)
     1 CALL HND_DDSPDS(NUMG,MROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
      IF(NUMG*MROOTS.GE.MINVEC)
     1 CALL HND_DDSPDV(NUMG,MROOTS,IJKLG,GINT,FINT,SINT,
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_INTTYP/NPOPLE,NHONDO
      COMMON/HND_RESTAR/NREST
      COMMON/HND_TIMEX/TIMLIM,TI,TX,TIM,WTI,WTX,WTIM
      COMMON/HND_IOFILE/IR,IW,IP
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
c*#include "nwc_const.fh"
c*    PARAMETER (MXATOM=NW_MAX_ATOM)
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
      subroutine hnd_cphf_init(rtdb)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
      integer rtdb
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
      parameter (mxatom=nw_max_atom)
      logical dbug
      logical out
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_cphf_mem/i_init,init
      common/hnd_basnum/num
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
      data zero    /0.0d+00/

c
      dbug=.false.
      out =.false.
      out =out.or.dbug

c     ----- get nocc and nvir -----

      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen, RTDB_ERR)
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

      if(.not.ma_push_get(mt_dbl,need ,
     &    'mem init:rhf_hessian:hnd_cphf_init:3',
     &    i_init,init))
     & call errquit('hnd_cphf_init, malloc of init  failed',911, MA_ERR)

c     ----- zero out the elements of the block -----

      do 10 i=1,need
   10 x(i+init-1)=zero
c
      return
      end
c
      subroutine hnd_cphf_popm 
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "mafdecls.fh"
      common/hnd_cphf_mem/i_init,init

      if(.not.ma_pop_stack(i_init))
     & call errquit('hnd_cphf_init, ma_pop_stack of init
     & failed',911, MA_ERR)

      return
      end
c
      subroutine hnd_cphfx(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
      integer  rtdb,geom,basis
      logical dbug
      logical out
      dbug=.false.
      out =.false.
      out =out.or.dbug
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
      subroutine hnd_cphf(rtdb,geom,basis)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "schwarz.fh"
#include "cfock.fh"
c#include "pstat.fh"
      parameter (mxatom=nw_max_atom)
      parameter(mxnbf=2048)
      integer  rtdb,geom,basis
c*    integer  g_fc
      integer  g_aodf(3*mxatom), g_modf(3*mxatom)
      integer  g_aods(3*mxatom), g_mods(3*mxatom)
      integer  g_dsta(3*mxatom), g_tmpta(3*mxatom)
      integer  g_modfwa(3*mxatom)
      integer  g_newaods,g_tmpaods
      integer  g_mafock(1),g_mmfock(1),g_tmpfock(1)
      double precision occ(mxnbf)
      double precision jfac(1),kfac(1)
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
      integer  g_ya, g_ya_tr
c
      logical dbug
      logical out
      logical out_grad
      common/hnd_iofile/ir,iw,ip
      common/hnd_dafile/idaf,nav,ioda(255)
      common/hnd_grad12/de(3,mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_cphf_mem/i_init,init
c
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      data two   /2.0d+00/ 
      data four  /4.0d+00/ 
      data eight /8.0d+00/
      data pt5   /0.5d+00/
      data nft17  /17/
      data nft18  /18/
c
      dbug     = .false.
      out      = .false.
      out      = out.or.dbug
      out_grad = .false.
      if(out) then
         write(iw,*) 'in hnd_cphf'
      endif

c     ----- get nocc and nvir -----

      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen, RTDB_ERR)

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

       if (.not.rtdb_cget(rtdb,'hess:input vectors',1,movecs_in))
     $     then
           call util_file_name('movecs',.false.,.false.,movecs_in)
       endif
       call util_file_name_resolve(movecs_in,.false.)
       if (.not.rtdb_cget(rtdb,'hess:output vectors',1,movecs_out))
     $     then               
           call util_file_name('movecs',.false.,.false.,movecs_out)
       endif
       call util_file_name_resolve(movecs_out,.false.)
       movecs = movecs_in

c      ----- read eig-vec and put eig-val in x(i40) -----

       status = movecs_read (movecs, 1, occ, x(i40), g_vecs)
       if (out.and.ga_nodeid().eq.0) then
           call hnd_prsq(  occ ,1,num,num)
           call hnd_prsq(x(i40),1,num,num)
           call ga_print(g_vecs)
       endif

c     ----- rewind nft17 and load into x(i70) -----

      call hnd_rewfil(nft17)

      do in =1,nc
c
        call hnd_hssrd(nft17,x(i70),num2) 

c     ----- expand from triangular matrix to square matrix ----

        do 450 j = 1,num
         do 450 i = 1,num
          ii = max0(i,j)
          jj = min0(i,j)
          ij = (ii*(ii-1))/2 + jj
  450   x(i50+num*(j-1)+i-1) = x(i70+ij-1)

        if (out.and.ga_nodeid().eq.0) then
            call hnd_prsq(x(i50),num,num,num)
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

c      ----- read eig-vec and put eig-val in x(i50) -----

       status = movecs_read (movecs, 1, occ, x(i50), g_vecs)

       if (.not.status) then
         call errquit('hnd_cphf: could not read mo vectors', 110,
     &       DISK_ERR)
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
     &             g_tmpfock(1),g_mafock(1), .false.)
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

*ga:1:0
      if(.not.ga_create(mt_dbl,nocc,nvir,'g_wax',0,0,g_wax))
     &     call errquit('could not alloc wax array',1, GA_ERR)

*ga:1:0
      if(.not.ga_create(mt_dbl,nvir,nocc,'g_wax_tr',0,0,g_wax_tr))
     &     call errquit('could not alloc tran array',1, GA_ERR)

*ga:1:0
      if(.not.ga_create(mt_dbl,nvir*nocc,1,'cphf_rhs:hnd_cphf',
     $    0,0,g_rhs))call errquit('hnd_cphf: cphf_rhs',0, GA_ERR)

*ga:1:0
      if(.not.ga_create(mt_dbl,nvir,nocc,'g_ya',0,0,g_ya))
     &     call errquit('could not alloc ya array',1, GA_ERR)

*ga:1:0
      if(.not.ga_create(mt_dbl,nocc,nvir,'g_ya_tr',0,0,g_ya_tr))
     &     call errquit('could not alloc ya_tr arrayy',1, GA_ERR)

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
     &     call errquit('could not write cphf_rhs',1, GA_ERR)

c     ----- call cphf solver ----

       call pstat_on(ps_cphf)
       if(.not.cphf(rtdb))call errquit('could not do cphf',1,
     &       CALC_ERR)
       call pstat_off(ps_cphf)

c     ----- the g_rhs now holds the solution -----

       if(.not.file_read_ga(cphf_sol,g_rhs))
     &    call errquit('could not read cphf_rhs',1, DISK_ERR)

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
     &     'could not destroy g_wax',1, GA_ERR)

      if(.not.ga_destroy(g_wax_tr))call errquit(
     &     'could not destroy g_wax_tr handle',1, GA_ERR)

      if(.not.ga_destroy(g_ya))call errquit(
     &     'could not destroy g_ya',1, GA_ERR)

      if(.not.ga_destroy(g_ya_tr))call errquit(
     &     'could not destroy g_ya_tr',1, GA_ERR)

      if(.not.ga_destroy(g_rhs))call errquit(
     &     'could not destroy g_ya_tr',1, GA_ERR)

c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,9997) nxyz
         call hnd_prsq(x(i10),nxyz,nrot,nrot)
         write(iw,9996) nxyz
         call hnd_prsq(x(i30),nxyz,nrot,nrot)
      endif

c     ----- add - W_x^a * x^b - contributions to -fcm_ab -----

c     ----- create ga block -----

*ga:1:0
      if(.not.ga_create(mt_dbl,nbov,nc,'g_wax_tot',0,0,g_wax_tot))
     &     call errquit('could not alloc wax_tot array',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nbov,nc,'g_xa_tot' ,0,0,g_xa_tot ))
     &     call errquit('could not alloc xa_tot  array',1, GA_ERR)

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
     & ga_ddot_patch(g_wax_tot,'N',1,nbov,ix,ix,
     &               g_xa_tot ,'N',1,nbov,iy,iy)

        enddo
      enddo

c     ----- output force constant matrix -----  
      if (out.and.ga_nodeid().eq.0) then
          write(iw,9990)
          call hnd_sprtr(x(i00),nc)
      endif  

c     ----- release memory pointer -----

      if(.not.ga_destroy(g_wax_tot))call errquit(
     &     'could not destroy g_wax_tot',1, GA_ERR)

      if(.not.ga_destroy(g_xa_tot ))call errquit(
     &     'could not destroy g_xa_tot',1, GA_ERR)

c     ----- for info only, print gradient -----

      if(out_grad.and.ga_nodeid().eq.0)
     &   call hnd_cpgrd(de,nat)

c     ----- symmetrixe -fcm- -----

c     if (ga_nodeid().eq.0) call hnd_cpfcm(x(i00),nat)
      if (ga_nodeid().eq.0) call hnd_hssprt(x(i00),nc)

c$$$c     ------ create ga block for force constant matrix -----
c$$$
c$$$      if(.not.ga_create(mt_dbl,nc,nc,'fcm',nc,nc,
c$$$     $     g_fc))call errquit('cphf: ga_create failed',0, GA_ERR)
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
c$$$     $   call errquit('rhf_hess: could not write scf_hess',1,
c&&&     &       GA_ERR)
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      COMMON/HND_IOFILE/IR,IW,IP
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
      subroutine rhf_hessian_2file(n,hess)
      implicit none
#include "errquit.fh"
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
      call errquit('hnd_cphf: fatal error', 911, DISK_ERR)
99991 write(6,*)' could not open <',fname,'> as new file'
      call errquit('hnd_cphf: fatal error', 911, DISK_ERR)
      end
      subroutine hnd_get_dens(rtdb,dens,num)
c
c     only RHF(close shell) is considered 
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
c
      integer           mxatom
      parameter        (mxatom=nw_max_atom)
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
      if(.not.ma_push_get(mt_dbl,need ,
     &    'mem init:rhf_hessian:hnd_get_dens:4',
     &    i_init,init))
     & call errquit('hnd_get_dens, malloc of init  failed',911, MA_ERR)
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
     & failed',911, MA_ERR)
c
      return
      end
c
      subroutine hnd_get_wdens(rtdb,wdens,num)
c
c     only RHF(close shell) is considered 
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
c
      integer           mxatom
      parameter        (mxatom=nw_max_atom)
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
      if(.not.ma_push_get(mt_dbl,need ,
     &    'mem init:rhf_hessian:hnd_get_wdens:5',
     &    i_init,init))
     & call errquit('hnd_get_wdens, malloc of init  failed',911, MA_ERR)
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
     & failed',911, MA_ERR)
c
      return
      end
