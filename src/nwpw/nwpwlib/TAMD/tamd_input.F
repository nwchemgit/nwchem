c $Id$

*     *********************************************
*     *                                           *
*     *              tamd_input                   *
*     *                                           *
*     *********************************************
      subroutine tamd_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "errquit.fh"
c
c

      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 24)
      character*18 dirs(num_dirs)
      data dirs / 'update',
     >            'zfixed',
     >            'print',
     >            'print_shift',
     >            'bond',
     >            'angle',
     >            'dihedral',
     >            'coord_number',
     >            'n-plane',
     >            'x',
     >            'y',
     >            'z',
     >            'tempered',
     >            'boundary',
     >            'potential_filename',
     >            'ld_trace',
     >            'sn2-surface',
     >            '2bonds',
     >            'ld_trace2diff',
     >            'bond_difference2',
     >            'equation',
     >            'group_distance',
     >            'bondings',
     >            'end'/

      integer num_dirs2          ! No. of known directives
      parameter (num_dirs2 = 14)
      character*10 dirs2(num_dirs2)
      data dirs2 / 'index1',
     >             'index2',
     >             'n',
     >             'm',
     >             'r0',
     >             'gamma',
     >             'k',
     >             'mass',
     >             'ztemp',
     >             'zinitial',
     >             'range',
     >             'nrange',
     >             'sprik','n12'/

      integer num_dirs3          ! No. of known directives
      parameter (num_dirs3 = 8)
      character*10 dirs3(num_dirs3)
      data dirs3 / 'gamma',
     >             'k',
     >             'mass',
     >             'ztemp',
     >             'zinitial',
     >             'range',
     >             'nrange',
     >             'nvector'/
c
      logical       value,spline,zfixed
      character*80  rtdb_name
      character*255 test,dir_name,sdir_name
      character*4   elem1,elem2
      integer ind,ind2,ind3           ! Index of matched directive

      integer metaindexmax
      integer nmeta,nxmeta_all,nindxmeta,nparammeta
      integer maxmetacount,metacount,metarayshift,metaraycount
      integer maxmetaprintcount
      integer j,jstart,jlast,jstride,n1,n2,i,d,jj(10)
      integer nka,nkb,nkc,nkd,nka_gs,nkb_gs,nkc_gs,nkd_gs

      integer sindx(50),index_start1(2),index_start2(2)
      integer pmeta(2),ameta(2),bmeta(2),indxmeta(2),sindxmeta(2)
      integer nxmeta(2),parammeta(2),sparammeta(2),ymeta(2)
      real*8  sdist(50),dTtempered

      integer nemach
      character*500 eqnstring

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name
      logical     nwpw_parse_boolean
      external    nwpw_parse_boolean
      integer     nwpw_expression_compile
      external    nwpw_expression_compile
      logical     tamd_bondings_input
      external    tamd_bondings_input

      metaindexmax = 1000
      nmeta = 0
      nxmeta_all = 1
      nindxmeta  = 0
      nparammeta = 0

      if (.not.BA_push_get(mt_int,10,'nxmeta',
     >                    nxmeta(2),nxmeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',0,
     >                MA_ERR)

      if (.not.BA_push_get(mt_int,10,'sindxmeta',
     >                    sindxmeta(2),sindxmeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',0,
     >                MA_ERR)

      if (.not.BA_push_get(mt_int,10,'sparammeta',
     >                    sparammeta(2),sparammeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',0,
     >                MA_ERR)

      if (.not.BA_push_get(mt_int,metaindexmax,'indxmeta',
     >                    indxmeta(2),indxmeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',0,
     >                MA_ERR)

      if (.not.BA_push_get(mt_int,10,'pmmeta',
     >                    pmeta(2),pmeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',2,
     >                MA_ERR)

      if (.not.BA_push_get(mt_dbl,10,'ammeta',
     >                    ameta(2),ameta(1)))
     >   call errquit('cannot allocate stack for tamd_input',2,
     >                MA_ERR)

      if (.not.BA_push_get(mt_dbl,10,'bmmeta',
     >                    bmeta(2),bmeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',3,
     >                MA_ERR)

      if (.not.BA_push_get(mt_dbl,metaindexmax,'parammeta',
     >                    parammeta(2),parammeta(1)))
     >   call errquit('cannot allocate stack for tamd_input',1,
     >                MA_ERR)


 10   if (.not. inp_read()) 
     >  call errquit(
     >  'tamd_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'tamd_input: failed to read keyword',0,INPUT_ERR)

      if (.not. inp_match(num_dirs,.false.,test,dirs,ind))
     >  call errquit(
     >  'tamd_input: unknown_directive',ind,INPUT_ERR)

      goto ( 100, 110, 120,150, 200, 300, 400, 500,600,701,702,703,800,
     >       900,1000,1100,1200,1300,1400,1500,1700,1800,1900,
     >      9999) ind
      call errquit(
     >     'tamd_input: unimplemented directive',ind,INPUT_ERR)

c
c  update
c
 100  if (inp_i(sindx(1))) then
         rtdb_name = 'tamd_update'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,sindx)
         if (.not. value) call errquit(
     >   'tamd_input:failed writing update',0, RTDB_ERR)
      end if

      goto 10

c
c  zfixed [True|False|...]
c
 110  if (inp_a(sdir_name)) then
         zfixed = nwpw_parse_boolean(sdir_name,.true.)
      else
         zfixed = .true.
      end if
      rtdb_name = 'tamd_zfixed'
      value = rtdb_put(rtdb,rtdb_name,mt_log,1,zfixed)
      if (.not. value) call errquit(
     >'tamd_input:failed writing update',0, RTDB_ERR)

      goto 10

c
c  print
c
 120  if (inp_i(sindx(1))) then
         rtdb_name = 'tamd_print'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,sindx)
         if (.not. value) call errquit(
     >   'tamd_input:failed writing update',0, RTDB_ERR)
      end if

      goto 10

c
c  print_shift
c
 150  if (inp_i(sindx(1))) then
         rtdb_name = 'tamd_print_shift'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,sindx)
         if (.not. value) call errquit(
     >   'tamd_input:failed writing update',150, RTDB_ERR)
      end if

      goto 10


c
c  bond
c
 200  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = -1.0
      sdist(6) = -1.0
      sdist(7) = -1.0
      sdist(8) = -1.0
      sindx(3) = 501

 201  if (.not.inp_a(test)) goto 299
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown bond directive',ind3,INPUT_ERR)

      goto ( 210,220,230,240,250,260,270) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 210  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 201
c
c k, default = 0.10 au
c
 220  if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 201
c
c mass, default=300.0 au
c
 230  if (.not.inp_f(sdist(3))) sdist(3)= 300.0d0
      goto 201
c
c ztemp, default = 1000.0K
c
 240  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 201

c
c zinitial
c
 250  if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 201

c
c range
c
 260  if (.not.inp_f(sdist(7))) sdist(7)=-1.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-1.0d0
      goto 201

c
c nrange, default = 501
c
 270  if (.not.inp_i(sindx(3))) sindx(3)=501
      goto 201



 299  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 1
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      nindxmeta  = nindxmeta + 3

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(3)
      int_mb(nxmeta(1)+nmeta) = sindx(3)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if
      end if

      goto 10

c
c  angle
c
 300  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 1.0d0
      sdist(6) = 0.0d0
      sdist(7) = 0.0d0
      sdist(8) = 4.0d0*datan(1.0d0)+0.1d0
      sindx(4) = 501


 301  if (.not.inp_a(test)) goto 399
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown angle directive',ind3,INPUT_ERR)

      goto (310,320,330,340,350,360,370) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)


c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 310  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 301
c
c k, default = 0.10 au
c
 320  if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 301
c
c mass, default=300.0 au
c
 330  if (.not.inp_f(sdist(3))) sdist(3)= 300.0d0
      goto 301
c
c ztemp, default = 1000.0K
c
 340  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 301
c
c zinitial
c
 350  if (.not.inp_f(sdist(5))) sdist(5)= 1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)= 0.0d0
      goto 301
c
c range
c
 360  if (.not.inp_f(sdist(7))) sdist(7)=-1.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-1.0d0
      goto 301
c
c nrange, default = 501
c
 370  if (.not.inp_i(sindx(4))) sindx(4)=501
      goto 301

 399  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 2
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      int_mb(indxmeta(1)+nindxmeta+3) = sindx(3)
      nindxmeta  = nindxmeta + 4

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(4)
      int_mb(nxmeta(1)+nmeta) = sindx(4)

      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing angle',0,RTDB_ERR)

      end if
      end if
      end if

      goto 10

c
c  dihedral
c
 400  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_i(sindx(4))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 1.0d0
      sdist(6) = 0.0d0
      sdist(7) = -4.0d0*datan(1.0d0)
      sdist(8) = 4.0d0*datan(1.0d0)
      sindx(5) = 501

 401  if (.not.inp_a(test)) goto 499
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown dihedral directive',ind3,INPUT_ERR)

      goto (410,420,430,440,450,460,470) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 410  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 401
c
c k, default = 0.10 au
c
 420  if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 401
c
c mass, default=300.0 au
c
 430  if (.not.inp_f(sdist(3))) sdist(3)= 300.0d0
      goto 401
c
c ztemp, default = 1000.0K
c
 440  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 401
c
c zinitial
c
 450  if (.not.inp_f(sdist(5))) sdist(5)= 1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)= 0.0d0
      goto 401
c
c range
c
 460  if (.not.inp_f(sdist(7))) sdist(7)=-4.0d0*datan(1.0d0)
      if (.not.inp_f(sdist(8))) sdist(8)= 4.0d0*datan(1.0d0)
      goto 401
c
c nrange, default = 501
c
 470  if (.not.inp_i(sindx(5))) sindx(5)=501
      goto 401


 499  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 3
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      int_mb(indxmeta(1)+nindxmeta+3) = sindx(3)
      int_mb(indxmeta(1)+nindxmeta+4) = sindx(4)
      nindxmeta  = nindxmeta + 5

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(5)
      int_mb(nxmeta(1)+nmeta) = sindx(5)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 1

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing dihedral',0,RTDB_ERR)

      end if
      end if
      end if
      end if

      goto 10

c
c  coord_number
c
 500  if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'tamd_input:failed allocating index_start1',0,MA_ERR)
      if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start2',
     >                    index_start2(2),index_start2(1)))
     >    call errquit(
     > 'tamd_input:failed allocating index_start2',0,MA_ERR)

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = -1.0
      sdist(6) = -1.0
      sdist(7) = -1.0
      sdist(8) = -1.0
      sdist(9) =  6.0
      sdist(10) = 12.0
      sdist(11) = 3.0
      sdist(12) = 1.0
      sdist(13) = 1.0
      sindx(1) = 501
      n1 = 0
      n2 = 0
 501  if (.not. inp_a(test)) goto 599
      if (.not. inp_match(num_dirs2,.false.,test,dirs2,ind2))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 510,520,530,540,550,560,562,564,566,570,580,590,595,
     >       545 ) ind2
      call errquit(
     >     'tamd_input: unimplemented directive',ind2,INPUT_ERR)


c
c index1
c
 510  n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do

      goto 501
c
c index2 
c
 520  n2 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start2(1)+n2) = j
          n2 = n2+1
         end do
      end do

      goto 501
c
c n
c
 530  if (.not.inp_f(sdist(9))) sdist(9)=6.0d0
      goto 501
c
c m
c
 540  if (.not.inp_f(sdist(10))) sdist(10)=12.0d0
      goto 501
c
c n12
c
 545  if (.not.inp_f(sdist(13))) sdist(13)=1.0d0
      goto 501

c
c r0
c
 550  if (.not.inp_f(sdist(11))) sdist(11)=3.0d0
      goto 501

c
c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 560  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 501
c
c k, default = 0.10 au
c
 562  if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 501

c
c mass, default=300.0 au
c
 564  if (.not.inp_f(sdist(3))) sdist(3)=300.0d0
      goto 501

c
c ztemp, default = 1000.0K
c
 566  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 501

c
c zinitial
c
 570  if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 501

c
c range
c
 580  if (.not.inp_f(sdist(7))) sdist(7)=0.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=10.0d0
      goto 501
c
c nrange, default = 501
c
 590  if (.not.inp_i(sindx(1))) sindx(1)=501
      goto 501

c
c sprik
c
 595  sdist(12)=1.0d0
      goto 501

 599  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 4
      int_mb(indxmeta(1)+nindxmeta+1) = n1
      int_mb(indxmeta(1)+nindxmeta+2) = n2
      do j=1,n1
         int_mb(indxmeta(1)+nindxmeta+2+j) 
     >      = int_mb(index_start1(1)+j-1)
      end do
      do j=1,n2
         int_mb(indxmeta(1)+nindxmeta+2+n1+j) 
     >      = int_mb(index_start2(1)+j-1)
      end do
      nindxmeta  = nindxmeta + 3+n1+n2

      int_mb(sparammeta(1)+nmeta)        = nparammeta
      dbl_mb(parammeta(1)+nparammeta)    = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1)  = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2)  = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3)  = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4)  = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5)  = sdist(6)
      dbl_mb(parammeta(1)+nparammeta+6)  = sdist(7)
      dbl_mb(parammeta(1)+nparammeta+7)  = sdist(8)
      dbl_mb(parammeta(1)+nparammeta+8)  = sdist(9)
      dbl_mb(parammeta(1)+nparammeta+9)  = sdist(10)
      dbl_mb(parammeta(1)+nparammeta+10) = sdist(11)
      dbl_mb(parammeta(1)+nparammeta+11) = sdist(12)
      dbl_mb(parammeta(1)+nparammeta+12) = sdist(13)
      nparammeta  = nparammeta + 13

      nxmeta_all = nxmeta_all*sindx(1)
      int_mb(nxmeta(1)+nmeta) = sindx(1)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing coord_num',0,RTDB_ERR)

      if (.not.BA_pop_stack(index_start2(2))) call errquit(
     > 'tamd_input:failed deallocating index_start2',0,MA_ERR)
      if (.not.BA_pop_stack(index_start1(2))) call errquit(
     > 'tamd_input:failed deallocating index_start1',0,MA_ERR)

      goto 10

c
c  n-plane
c
 600  if (inp_i(sindx(1))) then
      n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = -1.0
      sdist(6) = -1.0
      sdist(7) = -100.0
      sdist(8) =  100.0
      sindx(3) = 501


 601  if (.not.inp_a(test)) goto 699
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 610,620,630,640,650,660,670,680) ind3
      call errquit(
     >   'tamd_inputn: unimplemented directive',ind3,INPUT_ERR)
c
c gamma, default=0.0001 au = 0.06275 kcal/mol
c
 610  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 601
c
c k, default = 0.1 au
c
 620  if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 601
c
c mass
c
 630  if (.not.inp_f(sdist(3))) sdist(3)=300.0d0
      goto 601
c
c ztemp, default = 1000.0K
c
 640  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 601
c
c zinitial
c
 650  if (.not.inp_f(sdist(5))) sdist(5)=0.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=0.0d0
      goto 601
c
c range
c
 660  if (.not.inp_f(sdist(7))) sdist(7)=-100.0d0
      if (.not.inp_f(sdist(8))) sdist(8)= 100.0d0
      goto 601
c
c nrange, default = 501
c
 670  if (.not.inp_i(sindx(3))) sindx(3)=501
      goto 601
c
c normal
c
 680  if (.not.inp_f(sdist(9)))  sdist(3)=0.0d0
      if (.not.inp_f(sdist(10))) sdist(4)=0.0d0
      if (.not.inp_f(sdist(11))) sdist(5)=1.0d0
      goto 601

 699  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 5
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = n1
      do j=1,n1
         int_mb(indxmeta(1)+nindxmeta+2+j)
     >      = int_mb(index_start1(1)+j-1)
      end do
      nindxmeta  = nindxmeta + 3+n1

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      dbl_mb(parammeta(1)+nparammeta+6) = sdist(9)
      dbl_mb(parammeta(1)+nparammeta+7) = sdist(10)
      dbl_mb(parammeta(1)+nparammeta+8) = sdist(11)
      nparammeta  = nparammeta + 9

      nxmeta_all = nxmeta_all*sindx(2)
      int_mb(nxmeta(1)+nmeta) = sindx(3)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing n-plane',0,RTDB_ERR)

      end if
      goto 10

c
c  x,y,z
c
 701  sindx(1) = 1
      goto 704
 702  sindx(1) = 2
      goto 704
 703  sindx(1) = 3
      goto 704

 704  if (inp_i(sindx(2))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = -1.0
      sdist(6) = -1.0
      sdist(7) = -100.0
      sdist(8) =  100.0
      sindx(3) = 501

 705  if (.not.inp_a(test)) goto 799
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 710,720,730,740,750,760,770) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c gamma, default=0.0001 au = 0.06275 kcal/mol
c
 710  if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 705
c
c k, default = 0.1 au
c
 720  if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 705
c
c mass
c
 730  if (.not.inp_f(sdist(3))) sdist(3)=300.0d0
      goto 705
c
c ztemp, default = 1000.0K
c
 740  if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 705

c
c zinitial
c
 750  if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 705
c
c range
c
 760  if (.not.inp_f(sdist(7))) sdist(7)=-100.0d0
      if (.not.inp_f(sdist(8))) sdist(8)= 100.0d0
      goto 705

c
c nrange, default = 501
c
 770  if (.not.inp_i(sindx(3))) sindx(3)=501
      goto 705

 799  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 6
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      nindxmeta  = nindxmeta + 3

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(3)
      int_mb(nxmeta(1)+nmeta) = sindx(3)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing x,y,z',0,RTDB_ERR)

      end if
      goto 10
 


c
c  tempered
c
 800  if (.not.inp_f(sdist(1))) sdist(1)=-1.0d0
      rtdb_name='tamd_tempered'
      if (.not.rtdb_put(rtdb,rtdb_name,mt_dbl,1,sdist))
     >   call errquit('tamd_input:failed writing tempered',
     >                 0,RTDB_ERR)

      goto 10

c
c  boundary
c
 900  if (.not.inp_f(sdist(1))) sdist(1)= 0.0d0
      if (.not.inp_f(sdist(2))) sdist(2)= 1.0d0
      rtdb_name='tamd_boundary'
      if (.not.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist))
     >   call errquit('tamd_input:failed writing boundary',
     >                 0,RTDB_ERR)
      goto 10

c
c  potential_filename
c
 1000  if (inp_a(rtdb_name)) then
          call util_file_name_noprefix(rtdb_name,.false.,
     >                                 .false.,
     >                                 test)
          write(*,'(A,A)') 
     >      "....Reading tamd potential_filename=",rtdb_name
          open(unit=53,file=test,form='formatted')
 
          read(53,*) test
          read(53,*) test,nmeta
          read(53,*) test,(int_mb(nxmeta(1)+i-1),i=1,nmeta)
          read(53,*) test,(int_mb(sindxmeta(1)+i-1),i=1,nmeta)
          read(53,*) test,(int_mb(sparammeta(1)+i-1),i=1,nmeta)
          read(53,*) test,(dbl_mb(ameta(1)+i-1),i=1,nmeta)
          read(53,*) test,(dbl_mb(bmeta(1)+i-1),i=1,nmeta)
          read(53,*) test,(int_mb(pmeta(1)+i-1),i=1,nmeta)
          read(53,*) test,nindxmeta
          read(53,*) test,(int_mb(indxmeta(1)+i-1),i=1,nindxmeta)
          read(53,*) test,nparammeta
          read(53,*) test,(dbl_mb(parammeta(1)+i-1),i=1,nparammeta)
          read(53,*) test,maxmetacount
          read(53,*) test,metacount
          read(53,*) test,metarayshift
          read(53,*) test,metaraycount
          read(53,*) test,maxmetaprintcount
          read(53,*) test,dTtempered
          read(53,*) test,nxmeta_all

          if (.not.BA_push_get(mt_dbl,nxmeta_all,'ymmeta',
     >                         ymeta(2),ymeta(1)))
     >    call errquit('cannot allocate stack, tamd_input',22,
     >                 MA_ERR)

          do d=1,nmeta
             jj(d) = 0
          end do
          do i=1,nxmeta_all
             read(53,'(12F15.6)')
     >            (sdist(d),d=1,nmeta),dbl_mb(ymeta(1)+i-1)

             jj(1) = jj(1) + 1
             if (jj(1).ge.int_mb(nxmeta(1))) then
                read(53,*)
             end if
             do d=1,nmeta-1
                if (jj(d).ge.int_mb(nxmeta(1)+d-1)) then
                   jj(d) = 0
                   jj(d+1) = jj(d+1)+1
                end if
             end do
          end do
          close(53)

          rtdb_name = 'tamd_nindxmeta'
          value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
          rtdb_name = 'tamd_sindxmeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
          rtdb_name = 'tamd_indxmeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                     nindxmeta,int_mb(indxmeta(1)))
          rtdb_name = 'tamd_nparammeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
          rtdb_name = 'tamd_sparammeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                     nmeta,int_mb(sparammeta(1)))
          rtdb_name = 'tamd_parammeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                     nparammeta,dbl_mb(parammeta(1)))
          rtdb_name = 'tamd_ameta'
          value = value.and.
     >            rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
          rtdb_name = 'tamd_bmeta'
          value = value.and.
     >            rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
          rtdb_name = 'tamd_pmeta'
          value = value.and.
     >            rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
          rtdb_name = 'tamd_nxmeta'
          value = value.and.
     >           rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
          rtdb_name = 'tamd_nmeta'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
          rtdb_name = 'tamd_nxmeta_all'
          value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)

          rtdb_name = 'tamd_update'
          value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,maxmetacount)
          rtdb_name = 'tamd_metacount'
          value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,metacount)
          rtdb_name = 'tamd_print_shift'
          value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,metarayshift)
          rtdb_name = 'tamd_metaraycount'
          value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,metaraycount)
          rtdb_name = 'tamd_print'
          value=value.and.
     >          rtdb_put(rtdb,rtdb_name,mt_int,1,maxmetaprintcount)
          rtdb_name = 'tamd_tempered'
          value=value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,dTtempered)
          rtdb_name = 'tamd_ymeta'
          value=value.and.
     >       rtdb_put(rtdb,rtdb_name,mt_dbl,nxmeta_all,dbl_mb(ymeta(1)))

          if (.not. value) call errquit(
     >    'tamd_input:writing potential_filename',0,RTDB_ERR)

          if (.not.BA_pop_stack(ymeta(2)))
     >    call errquit('cannot pop stack for tamd_input',22,
     >                 MA_ERR)
      end if
      goto 10

c
c  ld_trace: local_density_trace ion_index l  ........
c
 1100 if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then

*     *** set defaults for w,sigma,range and nrange ***
      sdist(1) = 0.00005d0
      sdist(2) = 0.10000d0
      sdist(3) = -999.0d0
      sdist(4) = -999.0d0
      sindx(3) = 501

 1101 if (.not.inp_a(test)) goto 1199
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 1110,1120,1130,1140) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c w, default=0.00005 au = 0.03175 kcal/mol
c
 1110 if (.not.inp_f(sdist(1))) sdist(1)=0.00005d0
      goto 1101
c
c sigma, default = 0.1 au
c
 1120  if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 1101
c
c range
c
 1130  if (.not.inp_f(sdist(3))) sdist(3)=-999.0d0
      if (.not.inp_f(sdist(4))) sdist(4)=-999.0d0
      goto 1101
c
c nrange, default = 1101
c
 1140  if (.not.inp_i(sindx(3))) sindx(3)=501
      goto 1101


 1199  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 7
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      nindxmeta  = nindxmeta + 3

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      nparammeta  = nparammeta + 2

      nxmeta_all = nxmeta_all*sindx(3)
      int_mb(nxmeta(1)+nmeta) = sindx(3)
      dbl_mb(ameta(1)+nmeta)  = sdist(3)
      dbl_mb(bmeta(1)+nmeta)  = sdist(4)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if
      end if

      goto 10

c
c  sn2-surface type dr wr ar br dp wp ap bp er ep erp
c
 1200 if (inp_f(sdist(1)).and.
     >    inp_f(sdist(2)).and.
     >    inp_f(sdist(3)).and.
     >    inp_f(sdist(4)).and.
     >    inp_f(sdist(5)).and.
     >    inp_f(sdist(6)).and.
     >    inp_f(sdist(7)).and.
     >    inp_f(sdist(8)).and.
     >    inp_f(sdist(9)).and.
     >    inp_f(sdist(10)).and.
     >    inp_f(sdist(11)).and.
     >    inp_f(sdist(12)) ) then

      rtdb_name = 'tamd_sn2-surface'
      if (.not.rtdb_put(rtdb,rtdb_name,mt_dbl,12,sdist))
     >   call errquit('tamd_input:failed writing sn2-surface',
     >                 0,RTDB_ERR)
     
      end if

      goto 10
c
c  2bonds indx1 indx2 indx3 indx3 d1a d2a d1b d2b
c
 1300 if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_i(sindx(4))) then

      if (inp_f(sdist(3))) then
      if (inp_f(sdist(4))) then
      if (inp_f(sdist(5))) then
      if (inp_f(sdist(6))) then

*     *** set defaults for w,sigma,range and nrange ***
      sdist(1) = 0.00005d0
      sdist(2) = 0.10000d0
      sdist(7) = 0.0d0
      sdist(8) = 1.0d0
      sindx(5) = 501

 1301 if (.not.inp_a(test)) goto 1399
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 1310,1320,1330,1340) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c w, default=0.00005 au = 0.03175 kcal/mol
c
 1310 if (.not.inp_f(sdist(1))) sdist(1)=0.00005d0
      goto 1301
c
c sigma, default = 0.1 au
c
 1320 if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 1301
c
c range
c
 1330 if (.not.inp_f(sdist(7))) sdist(7)=0.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=1.0d0
      goto 1301
c
c nrange, default = 501
c
 1340 if (.not.inp_i(sindx(5))) sindx(5)=501
      goto 1301

 1399  continue

      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 8
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      int_mb(indxmeta(1)+nindxmeta+3) = sindx(3)
      int_mb(indxmeta(1)+nindxmeta+4) = sindx(4)
      nindxmeta  = nindxmeta + 5

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(5)
      int_mb(nxmeta(1)+nmeta) = sindx(5)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if

      goto 10


c
c  ld_trace2diff: local_density_trace ion_index1 ion_index2 l  ........
c
 1400 if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then

*     *** set defaults for w,sigma,range and nrange ***
      sdist(1) = 0.00005d0
      sdist(2) = 0.10000d0
      sdist(3) = -999.0d0
      sdist(4) = -999.0d0
      sindx(4) = 501

 1401 if (.not.inp_a(test)) goto 1499
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 1410,1420,1430,1440) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c w, default=0.00005 au = 0.03175 kcal/mol
c
 1410 if (.not.inp_f(sdist(1))) sdist(1)=0.00005d0
      goto 1401
c
c sigma, default = 0.1 au
c
 1420  if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 1401
c
c range
c
 1430  if (.not.inp_f(sdist(3))) sdist(3)=-999.0d0
      if (.not.inp_f(sdist(4))) sdist(4)=-999.0d0
      goto 1401
c
c nrange, default = 501
c
 1440  if (.not.inp_i(sindx(4))) sindx(4)=501
      goto 1401


 1499  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 9
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      int_mb(indxmeta(1)+nindxmeta+3) = sindx(3)
      nindxmeta  = nindxmeta + 4

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      nparammeta  = nparammeta + 2

      nxmeta_all = nxmeta_all*sindx(4)
      int_mb(nxmeta(1)+nmeta) = sindx(4)
      dbl_mb(ameta(1)+nmeta)  = sdist(3)
      dbl_mb(bmeta(1)+nmeta)  = sdist(4)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if
      end if
      end if

      goto 10

c
c  bond_difference2 indx1 indx2 indx3 
c
 1500 if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 0.0d0
      sdist(6) = 0.0d0
      sdist(7) = -20.0d0
      sdist(8) = 20.0d0
      sindx(5) = 501
 1501 if (.not.inp_a(test)) goto 1599
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto ( 1510,1520,1530,1540,1550,1560,1570,1501) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)
c
c gamma, default=0.0001 
c
 1510 if (.not.inp_f(sdist(1))) sdist(1)=0.00005d0
      goto 1501
c
c k default = 0.1
c
 1520 if (.not.inp_f(sdist(2))) sdist(2)=0.1d0
      goto 1501
c
c mass           
c
 1530 if (.not.inp_f(sdist(2))) sdist(3)=300.0d0
      goto 1501
c
c ztemp          
c
 1540 if (.not.inp_f(sdist(2))) sdist(4)=1000.0d0
      goto 1501
c
c zinitial
c
 1550 if (.not.inp_f(sdist(7))) sdist(5)=0.0d0
      if (.not.inp_f(sdist(8))) sdist(6)=0.0d0
c
c range
c
 1560 if (.not.inp_f(sdist(7))) sdist(7)= -20.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=  20.0d0
      goto 1501
c
c nrange, default = 501
c
 1570 if (.not.inp_i(sindx(5))) sindx(5)=501
      goto 1501

 1599  continue

      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 10
      int_mb(indxmeta(1)+nindxmeta+1) = sindx(1)
      int_mb(indxmeta(1)+nindxmeta+2) = sindx(2)
      int_mb(indxmeta(1)+nindxmeta+3) = sindx(3)
      nindxmeta  = nindxmeta + 4

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(5)
      int_mb(nxmeta(1)+nmeta) = sindx(5)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if
      end if
      end if

      goto 10


c
c  equation
c
 1700 if (inp_a(eqnstring)) then

      nemach = nwpw_expression_compile(rtdb,eqnstring)

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 0.0d0
      sdist(6) = 0.0d0
      sdist(7) = -1.0d0
      sdist(8) = -1.0d0
      sindx(3) = 501

 1701 if (.not.inp_a(test)) goto 1799
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto (1710,1720,1730,1740,1750,1760,1770) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)

c
c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 1710 if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 1701
c
c k, default = 0.10 au
c
 1720 if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 1701
c
c mass, default=300.0 au
c
 1730 if (.not.inp_f(sdist(3))) sdist(3)= 300.0d0
      goto 1701
c
c ztemp, default = 1000.0K
c
 1740 if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 1701

c
c zinitial
c
 1750 if (.not.inp_f(sdist(5))) sdist(5)=0.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=0.0d0
      goto 1701

c
c range
c
 1760 if (.not.inp_f(sdist(7))) sdist(7)=-1.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-1.0d0
      goto 1701

c
c nrange, default = 501
c
 1770 if (.not.inp_i(sindx(3))) sindx(3)=501
      goto 1701


 1799  continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 12
      int_mb(indxmeta(1)+nindxmeta+1) = nemach
      nindxmeta  = nindxmeta + 2

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      nparammeta  = nparammeta + 6

      nxmeta_all = nxmeta_all*sindx(3)
      int_mb(nxmeta(1)+nmeta) = sindx(3)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing bond',0,RTDB_ERR)

      end if

      goto 10


c
c  group_distance
c
 1800 if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'tamd_input:failed allocating index_start1',0,MA_ERR)
      if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start2',
     >                    index_start2(2),index_start2(1)))
     >    call errquit(
     > 'tamd_input:failed allocating index_start2',0,MA_ERR)


*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 0.0d0
      sdist(6) = 0.0d0
      sdist(7) = -1.0d0
      sdist(8) = -1.0d0
      sdist(9) = 6.0d0
      sindx(1) = 501
      n1 = 0
      n2 = 0
 1801 if (.not. inp_a(test)) goto 1899
      if (.not. inp_match(num_dirs2,.false.,test,dirs2,ind2))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto (1810,1820,1830,1801,1801,1860,1862,1864,1866,
     >      1870,1880,1890,1801,1801) ind2
      call errquit(
     >     'tamd_input: unimplemented directive',ind2,INPUT_ERR)
c
c index1
c
 1810 n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do

      goto 1801
c
c index2 
c
 1820 n2 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start2(1)+n2) = j
          n2 = n2+1
         end do
      end do

      goto 1801
c
c n
c
 1830 if (.not.inp_f(sdist(9))) sdist(9)=6.0d0
      goto 1801

c
c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 1860 if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 1801
c
c k, default = 0.10 au
c
 1862 if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 1801

c
c mass, default=300.0 au
c
 1864 if (.not.inp_f(sdist(3))) sdist(3)=300.0d0
      goto 1801

c
c ztemp, default = 1000.0K
c
 1866 if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 1801

c
c zinitial
c
 1870 if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 1801

c
c range
c
 1880 if (.not.inp_f(sdist(7))) sdist(7)=-999.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-999.0d0
      goto 1801
c
c nrange, default = 501
c
 1890 if (.not.inp_i(sindx(1))) sindx(1)=501
      goto 1801



 1899 continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 13
      int_mb(indxmeta(1)+nindxmeta+1) = n1
      int_mb(indxmeta(1)+nindxmeta+2) = n2
      do j=1,n1
         int_mb(indxmeta(1)+nindxmeta+2+j)
     >      = int_mb(index_start1(1)+j-1)
      end do
      do j=1,n2
         int_mb(indxmeta(1)+nindxmeta+2+n1+j)
     >      = int_mb(index_start2(1)+j-1)
      end do
      nindxmeta  = nindxmeta + 3+n1+n2

      int_mb(sparammeta(1)+nmeta)       = nparammeta
      dbl_mb(parammeta(1)+nparammeta)   = sdist(1)
      dbl_mb(parammeta(1)+nparammeta+1) = sdist(2)
      dbl_mb(parammeta(1)+nparammeta+2) = sdist(3)
      dbl_mb(parammeta(1)+nparammeta+3) = sdist(4)
      dbl_mb(parammeta(1)+nparammeta+4) = sdist(5)
      dbl_mb(parammeta(1)+nparammeta+5) = sdist(6)
      dbl_mb(parammeta(1)+nparammeta+6) = sdist(9)
      nparammeta  = nparammeta + 7

      nxmeta_all = nxmeta_all*sindx(1)
      int_mb(nxmeta(1)+nmeta) = sindx(1)
      dbl_mb(ameta(1)+nmeta)  = sdist(7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing coord_num',0,RTDB_ERR)

      if (.not.BA_pop_stack(index_start2(2))) call errquit(
     > 'tamd_input:failed deallocating index_start2',0,MA_ERR)
      if (.not.BA_pop_stack(index_start1(2))) call errquit(
     > 'tamd_input:failed deallocating index_start1',0,MA_ERR)

      goto 10



c
c  bondings
c
 1900 if (tamd_bondings_input(n1,sindx(2),sdist(7))) then

*     **** set defaults for gamma,k,mass,ztemp,zinitial ****
      sdist(1) =    0.00010d0
      sdist(2) =    0.10000d0
      sdist(3) =  300.0d0
      sdist(4) = 1000.0d0
      sdist(5) = 0.0d0
      sdist(6) = 0.0d0
      sdist(n1+7) = -20.0d0
      sdist(n1+8) =  20.0d0
      sindx(1) = 501
 1901 if (.not. inp_a(test)) goto 1999
      if (.not. inp_match(num_dirs3,.false.,test,dirs3,ind3))
     >  call errquit(
     >  'tamd_input: unknown directive',0,INPUT_ERR)

      goto (1910,1920,1930,1940,1950,1960,1970,1901) ind3
      call errquit(
     >     'tamd_input: unimplemented directive',ind3,INPUT_ERR)


c
c gamma, default=0.00010 au = 0.06275 kcal/mol
c
 1910 if (.not.inp_f(sdist(1))) sdist(1)=0.00010d0
      goto 1901
c
c k, default = 0.10 au
c
 1920 if (.not.inp_f(sdist(2))) sdist(2)=0.10d0
      goto 1901

c
c mass, default=300.0 au
c
 1930 if (.not.inp_f(sdist(3))) sdist(3)=300.0d0
      goto 1901

c
c ztemp, default = 1000.0K
c
 1940 if (.not.inp_f(sdist(4))) sdist(4)=1000.0d0
      goto 1901

c
c zinitial
c
 1950 if (.not.inp_f(sdist(5))) sdist(5)=0.0d0
      if (.not.inp_f(sdist(6))) sdist(6)=-0.0d0
      goto 1901

c
c range
c
 1960 if (.not.inp_f(sdist(n1+7))) sdist(n1+7)= -20.0d0
      if (.not.inp_f(sdist(n1+8))) sdist(n1+8)=  20.0d0
      goto 1901
c
c nrange, default = 501
c
 1970 if (.not.inp_i(sindx(1))) sindx(1)=501
      goto 1901



 1999 continue
      int_mb(sindxmeta(1) +nmeta)     = nindxmeta
      int_mb(indxmeta(1)+nindxmeta)   = 14
      int_mb(indxmeta(1)+nindxmeta+1) = n1
      do i=1,n1
         int_mb(indxmeta(1)+nindxmeta + 2*i)   = sindx(2*i)
         int_mb(indxmeta(1)+nindxmeta + 2*i+1) = sindx(2*i+1)
      end do
      nindxmeta  = nindxmeta + 2*n1 + 2

      int_mb(sparammeta(1)+nmeta) = nparammeta
      do i=1,n1+6
         dbl_mb(parammeta(1)+nparammeta+i-1) = sdist(i)
      end do
      nparammeta  = nparammeta + n1 + 6

      nxmeta_all = nxmeta_all*sindx(1)
      int_mb(nxmeta(1)+nmeta) = sindx(1)
      dbl_mb(ameta(1)+nmeta)  = sdist(n1+7)
      dbl_mb(bmeta(1)+nmeta)  = sdist(n1+8)
      int_mb(pmeta(1)+nmeta)  = 0

      nmeta = nmeta + 1

      rtdb_name = 'tamd_nindxmeta'
      value =           rtdb_put(rtdb,rtdb_name,mt_int,1,nindxmeta)
      rtdb_name = 'tamd_sindxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sindxmeta(1)))
      rtdb_name = 'tamd_indxmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nindxmeta,int_mb(indxmeta(1)))
      rtdb_name = 'tamd_nparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nparammeta)
      rtdb_name = 'tamd_sparammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,
     >                 nmeta,int_mb(sparammeta(1)))
      rtdb_name = 'tamd_parammeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                 nparammeta,dbl_mb(parammeta(1)))
      rtdb_name = 'tamd_ameta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(ameta(1)))
      rtdb_name = 'tamd_bmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,nmeta,dbl_mb(bmeta(1)))
      rtdb_name = 'tamd_pmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(pmeta(1)))
      rtdb_name = 'tamd_nxmeta'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,nmeta,int_mb(nxmeta(1)))
      rtdb_name = 'tamd_nmeta'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nmeta)
      rtdb_name = 'tamd_nxmeta_all'
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nxmeta_all)
      if (.not. value) call errquit(
     >'tamd_input:failed writing coord_num',0,RTDB_ERR)

      end if
      goto 10



 9999 continue
      
*     **** pop stack ****
      if (.not.BA_pop_stack(parammeta(2)))
     >   call errquit('cannot pope stack for tamd_input',1,
     >                MA_ERR)
      if (.not.BA_pop_stack(bmeta(2)))
     >   call errquit('cannot pop stack for tamd_input',3,
     >                MA_ERR)
      if (.not.BA_pop_stack(ameta(2)))
     >   call errquit('cannot pop stack for tamd_input',2,
     >                MA_ERR)
      if (.not.BA_pop_stack(pmeta(2)))
     >   call errquit('cannot pop stack for tamd_input',0,
     >                MA_ERR)
      if (.not.BA_pop_stack(indxmeta(2)))
     >   call errquit('cannot pop stack for tamd_input',0,
     >                MA_ERR)
      if (.not.BA_pop_stack(sparammeta(2)))
     >   call errquit('cannot pop stack for tamd_input',0,
     >                MA_ERR)
      if (.not.BA_pop_stack(sindxmeta(2)))
     >   call errquit('cannot pop stack for tamd_input',0,
     >                MA_ERR)
      if (.not.BA_pop_stack(nxmeta(2)))
     >   call errquit('cannot pop stack for tamd_input',0,
     >                MA_ERR)

      return
      end



     
*     ****************************************
*     *                                      *
*     *          tamd_bondings_input         * 
*     *                                      *
*     ****************************************
*
* bondings -1.0 1 2 1.0 2 3
*
      logical function tamd_bondings_input(n,indx,coef)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"  
#include "util.fh"  
c     
      logical ok 
      integer n,indx(*)
      double precision coef(*)
      
      integer itest 
      double precision ftest
 
      ok = .true.
     
      !*** coef and indx ***
      n = 0
      do while (ok.and.inp_f(ftest))
         coef(n+1) = ftest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+1) = itest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+2) = itest
         n = n + 1
      end do
      if (n.eq.0) ok = .false.

      tamd_bondings_input = ok
      return
      end

