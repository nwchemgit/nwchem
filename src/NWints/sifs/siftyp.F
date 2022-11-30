      subroutine siftyp( itypea, itypeb, chrtyp )
c
c  return a character description of the integral type.
c
c  input:
c  itypea, itypeb = generic and specific integral or energy(*) types.
c
c  output:
c  chrtyp = character description. (this should be at least character*8
c           in the calling program.)
c
c  08-oct-90 total energy and convergence types added. -rls
c  01-jun-90 table version. -rls
c  11-aug-89 written by ron shepard.
c
      implicit none
c
c     # dummy...
      integer          itypea,  itypeb
      character*(*)    chrtyp
c
c     # local...
      integer  i
c
      integer  typ1e,  st1e,   end1e
      integer  typ2e,  st2e,   end2e
      integer  typc,   stc,    endc
      integer  typte,  stte,   endte
      integer  typcv,  stcv,   endcv
c
c     # typ1e = number of defined 1-e array types.
c     # typ2e = number of defined 2-e array types.
c     # typc  = number of core energy types.
c     # typte = number of total energy types.
c     # typcv = number of convergence types.
c
      parameter( typ1e = 55, st1e = 1,        end1e = typ1e        )
      parameter( typ2e =  2, st2e = end1e +1, end2e = end1e +typ2e )
      parameter( typc  =  1, stc  = end2e +1, endc  = end2e +typc  )
      parameter( typte = 17, stte = endc  +1, endte = endc  +typte )
      parameter( typcv = 10, stcv = endte +1, endcv = endte +typcv )
c
      integer    ntype
      parameter( ntype = endcv )
c
      integer ltypea(ntype), ltypeb(ntype)
      character*8 lctype(ntype)
c
c     # these array types can be in any order, so new ones can be added
c     # to the end or inserted into the middle as appropriate. -rls
c
c     # warning: case dependent code:  do not change the case of the
c     # following character strings.
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st1e, end1e ) /
     & 0,    0, 'S1(*)',     0,    1, 'T1(*)',
     & 0,    2, 'V1(*)',     0,    3, 'Veff(*)',
     & 0,    4, 'VFC(*)',    0,    5, 'Vref(*)',
     & 0,    6, 'H1(*)',     0,    7, 'D1(*)',
     & 0,    8, 'F(*)',      0,    9, 'Q(*)',
     & 1,    0, 'X(*)',      1,    1, 'Y(*)',
     & 1,    2, 'Z(*)',
     & 1,    3, 'XX(*)',      1,    4, 'XY(*)',
     & 1,    5, 'XZ(*)',      1,    6, 'YY(*)',
     & 1,    7, 'YZ(*)',      1,    8, 'ZZ(*)',
     & 1,    9, 'XXX(*)',     1,   10, 'XXY(*)',
     & 1,   11, 'XXZ(*)',     1,   12, 'XYY(*)',
     & 1,   13, 'XYZ(*)',     1,   14, 'XZZ(*)',
     & 1,   15, 'YYY(*)',     1,   16, 'YYZ(*)',
     & 1,   17, 'YZZ(*)',     1,   18, 'ZZZ(*)',
     &
     & 1,   19, 'XXXX(*)',    1,   20, 'XXXY(*)',
     & 1,   21, 'XXXZ(*)',    1,   22, 'XXYY(*)',
     & 1,   23, 'XXYZ(*)',    1,   24, 'XXZZ(*)',
     & 1,   25, 'XYYY(*)',    1,   26, 'XYYZ(*)',
     & 1,   27, 'XYZZ(*)',    1,   28, 'XZZZ(*)',
     & 1,   29, 'YYYY(*)',    1,   30, 'YYYZ(*)',
     & 1,   31, 'YYZZ(*)',    1,   32, 'YZZZ(*)',
     & 1,   33, 'YZZZ(*)',
     & 1,   34, 'massvel',    1,   35, 'darwin',
     &
     & 2,    0, 'Im(SO:x)',  2,    1, 'Im(SO:y)',
     & 2,    2, 'Im(SO:z)' ,
     & 2,    3, 'Im(px)  ',   2,   4, 'Im(py)  ',
     & 2,    5, 'Im(pz)  ',   2,   6, 'Im(lx)  ',
     & 2,    7, 'Im(ly)  ',   2,   8, 'Im(lz)  '/
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st2e, end2e ) /
     & 3,    0, '1/r12',     3,    1, 'd2(*)'     /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stc, endc ) /
     & 0,   -1, 'Nuc.Rep.' /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stte, endte ) /
     & -1,   0, 'SCF',      -1,   -1, 'MCSCF',
     & -1,  -2, 'MRSDCI',   -1,   -3, 'CPF',
     & -1,  -4, 'ACPF',     -1,   -5, 'LCC-SD',
     & -1,  -6, 'MRPT',     -1,   -7, 'Bk',
     & -1,  -8, 'DV1',      -1,   -9, 'DV2',
     & -1, -10, 'EPOPLE',   -1,  -11, 'S.O. CI',
     & -1, -12, 'SR-SDCI',  -1,  -13, 'UCEPA'  ,
     & -1, -14, 'MR-AQCC',  -1,  -15, 'a4den'  ,
     & -1, -16, 'LRT-AQCC'   /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stcv, endcv ) /
     & -2,   0, 'SCF-D.E.', -2,   -1, 'SCF-D.D1',
     & -2,  -2, 'MC-D.E.',  -2,   -3, 'MC-Wnorm',
     & -2,  -4, 'MC-Knorm', -2,   -5, 'MC-ApxDE',
     & -2,  -6, 'Bk-Resid', -2,   -7, 'CI-Resid',
     & -2,  -8, 'CI-D.E.',  -2,   -9, 'CI-ApxDE' /
c
      do 10 i = 1, ntype
         if ( itypea .eq. ltypea(i) ) then
            if ( itypeb .eq. ltypeb(i) ) then
               chrtyp = lctype(i)
               return
            endif
         endif
10    continue
c     # loop exit means unrecognized type.
      chrtyp = 'Unknown'
c
      return
      end
