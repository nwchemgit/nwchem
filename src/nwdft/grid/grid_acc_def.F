      Subroutine grid_acc_def(rtdb)
c
c$Id: grid_acc_def.F,v 1.17 2000-02-23 02:47:52 edo Exp $
c
      implicit none
c      
#include "inp.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
c
      integer rtdb
      integer num_gridtypes, grid_acc_index,jj,igridin
      parameter (num_gridtypes = 10)
      logical oprint, oprint_grid,lgridin
      character*2 symbol
      character*15 gridtypes(num_gridtypes)
      double precision grid_acc(num_gridtypes)
      double precision job_grid_acc
      character*16 element
      character*16 tag
      double precision eps
      parameter (eps = 1.d-20)
      data gridtypes /'xcoarse', 'coarse', 'medium', 'fine', 'xfine', 
     &                'lebmed','lebfine','lebxfine',
     .                'gausleg','lebedev'/
      data grid_acc /1d-4,1d-5,1d-6,1d-7,1d-8,1d-5,1d-7,1d-8,1d-8,1d-4/
c
      double precision ictr_coord(3), ictr_chg
      integer itype, ictr, iaz, i_atomic_number,iangat
c
      integer mxnrad, mxnang
c
c     Periodic table of elements and suggested number of radial
c     and angular points to achieve desired accuracy:
c     xcoarse = 1d-4
c     coarse  = 1d-5
c     medium  = 1d-6
c     fine    = 1d-7
c     xfine   = 1d-8
c
      integer grid_params
      integer gra
      integer grb,grc
      dimension grid_params(8,2,105)
      dimension gra(8,2,40)
      dimension grb(8,2,40)
      dimension grc(8,2,25)
      integer xc_leb_acc2nrd
      external xc_leb_acc2nrd
      equivalence(grid_params(1,1,1), gra(1,1,1))
      equivalence(grid_params(1,1,41), grb(1,1,1))
      equivalence(grid_params(1,1,81), grc(1,1,1))
* equivalence must come before data statements
c
c      It appears that H (and so likely He as well) can't be reduced as much as
c      previously thought.  Set these to Li - Ne values.
c
c      Data grid_paramS_a/ 24, 28, 36, 48, 58,  5,  7,  9, 13, 15, !  H (1)
c     &                    26, 36, 46, 60, 80,  7, 10, 12, 17, 23, ! He (2)
      Data gra/ 30, 50, 70,100,140, 45, 60,100,
     .     12, 15, 18, 24, 34, 6,  7, 10,!  H (1)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 10,! He (2)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,! Li (3)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,! Be (4)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,!  B (5)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,!  C (6)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,!  N (7)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,!  O (8)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,!  F (9)
     &          30, 50, 70,100,140, 49, 70,100,
     .     12, 15, 18, 24, 34, 6,  7, 11,! Ne (10)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Na (11)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Mg (12)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Al (13)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Si (14)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,!  P (15)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,!  S (16)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Cl (17)
     &          45, 75, 95,125,175, 88,123,125,
     .     12, 18, 24, 30, 44, 7,  8, 11,! Ar (18)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,!  K (19)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Ca (20)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Sc (21)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Ti (22)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,!  V (23)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Cr (24)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Mn (25)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Fe (26)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Co (27)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Ni (28)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Cu (29)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Zn (30)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Ga (31)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Ge (32)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! As (33)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Se (34)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Br (35)
     &          75, 95,110,160,210,112,130,160,
     .     14, 22, 30, 34, 38, 7,  9, 11,! Kr (36)
     &         105,130,155,205,235,123,155,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Rb (37)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Sr (38)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,!  Y (39)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11/! Zr (40)
      Data grb/ 105,130,155,205,235,123,141,205, 
     .     16, 20, 32, 44, 48, 7,  9, 11, ! Nb (41)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Mo (42)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Tc (43)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ru (44)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Rh (45)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pd (46)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ag (47)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Cd (48)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! In (49)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Sn (50)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Sb (51)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Te (52)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,!  I (53)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7, 48, 11,! Xe (54)
     .     105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Cs (55)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ba (56)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! La (57)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ce (58)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pr (59)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Nd (60)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pm (61)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Sm (62)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Eu (63)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Gd (64)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Tb (65)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Dy (66)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ho (67)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Er (68)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Tm (69)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Yb (70) 
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Lu (71)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Hf (72)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ta (73)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,!  W (74)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Re (75)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Os (76)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ir (77)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pt (78)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Au (79)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11/! Hg (80)
      Data grc/105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11, ! Tl (81)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pb (82)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Bi (83)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Po (84)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! At (85)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Rn (86)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Fr (87)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ra (88)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Ac (89)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Th (90)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pa (91)
     &         105,130,155,205,235,123,155,205,
     .     16, 20, 32, 44, 48, 7,  10, 11,!  U (92)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Np (93)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Pu (94)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Am (95)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Cm (96)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Bk (97)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Cf (98)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Es (99)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Fm (100)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Md (101)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! No (102)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Lr (103)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11,! Unq(104)
     &         105,130,155,205,235,123,141,205,
     .     16, 20, 32, 44, 48, 7,  9, 11/! Unp(105)
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
c     initialize iangquad for lebmed and lebfine case
c
      iangquad=0
c
c     Determine desired integration accuracy.
c
      if (.not. rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype))then
c     
c        set default grid to  'MEDIUM'
c     
         gridtype='medium'
      endif
      if (.not. rtdb_get(rtdb, 'dft:ngridin', mt_int, 1, ngridin))
     .     then
         ngridin=0
      else
         if (.not. rtdb_cget(rtdb,'dft:grtagin',ngridin,grtagin))
     &        call errquit('grid_inp: rtdb_put failed', 1000)
         if(.not.rtdb_get(rtdb,'dft:nradin',mt_int,ngridin,nradin))
     &        call errquit('grid_inp: rtdb_get failed', 1000)
         if(.not.rtdb_get(rtdb,'dft:nangin',mt_int,ngridin,nangin))
     &        call errquit('grid_inp: rtdb_get failed', 1000)
         if(.not.rtdb_get(rtdb,'dft:accqrin',mt_int,ngridin,accqrin))
     &        call errquit('grid_inp: rtdb_get failed', 1000)
      endif
      if (.not. inp_match(num_gridtypes, .false., gridtype, gridtypes, 
     &   grid_acc_index))then
         write(6,*)' Unknown gridtype; defaulting to medium grid'
         gridtype='medium'
         grid_acc_index = 3
      endif
      goto (8, 8, 8, 8, 8, 8, 8, 8, 6, 7) grid_acc_index
      call errquit('grid_acc_def: unknown gridtype', grid_acc_index)
    6 continue
c
c     user specified gausleg grid;
c
      if (.not. rtdb_get(rtdb, 'dft:nagrid', mt_int, 1, nagrid)) then
         gridtype='medium'
         grid_acc_index = 3
      endif      
      if (.not. rtdb_get(rtdb,'dft:nradpts', mt_int, 1, nradpts)) then
         gridtype='medium'
         grid_acc_index = 3
      endif      
      goto 8
    7 continue
c
c     lebedev grid;
c
      leb=.true.
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1, iangquad))then
         iangquad = 4
      else
         if (iangquad.lt.1)then
           write(6,*)' IANGQUAD must be in the range: 1-',IAC_MAX,
     *     '; resetting. '
           iangquad = 4
         endif
         if (iangquad.gt.IAC_MAX)then
           write(6,*)' IANGQUAD must be in the range: 1-',IAC_MAX,
     *     '; resetting. '
           iangquad = 4
         endif
      endif
      if (.not. rtdb_get(rtdb,'dft:nradpts', mt_int, 1, nradpts)) then
         gridtype='medium'
         grid_acc_index = 3
      endif      
    8 continue
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     .     job_grid_acc = grid_acc(grid_acc_index)

c
c     shove job_grid_acc into the rtdb
c
      if (.not. rtdb_put(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('grid_acc_def: rtdb_put failed', 911)
c     
c     Define the quadrature size for each atomic center: 
c      - nrad, nang, mxnrad, mxnang
      mxnrad = 0
      mxnang = 0
c     
      do 50 itype = 1, ntypes
c
c        find an atom of this kind in the complete list
c
         do ictr = 1, ncenters
            if (iatype(ictr).eq.itype) then
               iaz = ictr
               if (.not. geom_cent_get(geom, ictr, tag, 
     &            ictr_coord, ictr_chg))call errquit
     &            ('grid_acc_def: geom_cent_get failed', 0)
               goto 40
            endif
         enddo
   40    continue
c        
         if (abs(znuc_atom_type(itype)).lt.EPS) then ! uncharged ghost atom; add grid
c
c           add quadrature grids to uncharged ghost atoms (for BSSE calcs)
c
c           identify atom label following "bq"
c
            if (.not. geom_tag_to_element(tag(3:), symbol, 
     &           element, i_atomic_number)) then
               if (inp_compare(.false.,tag(1:2),'bq')) then 
                  i_atomic_number = 0
               elseif (inp_compare(.false.,tag(1:1),'X')) then 
                  i_atomic_number = 0
               else
                  call errquit
     &              ('grid_acc_def: non-bq center with zero charge', 0)
               endif
            endif
            if (i_atomic_number.eq.0)then
               nrad(itype) = 0
               nang(itype) = 0
            else
               if (grid_acc_index.gt.8)then
                  nrad(itype) = nradpts
                  nang(itype) = nagrid
               else
                     nrad(itype) = grid_params(grid_acc_index,1,
     &                                    i_atomic_number)
                  if(grid_acc_index.ge.6.and.grid_acc_index.le.8) then
                     leb=.true.
                     iangquad=grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
                     nang(itype) = xc_leb_acc2nrd(iangquad)
                  else
                  nang(itype) = grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
                  endif
               endif
c                
               if(nrad(itype).gt.mxnrad)mxnrad = nrad(itype)
               if(nang(itype).gt.mxnang)mxnang = nang(itype)
            endif
         else    ! center is charged
c
c           no quadrature grids on charged ghost atoms 
c
            if (.not. geom_tag_to_element(tag, symbol, 
     &         element, i_atomic_number)) then
               if (symbol .ne. 'bq') call errquit
     &              ('grid_acc_def: center is neither atom nor bq', 0)
            endif
c
            if (i_atomic_number.ne.0)then ! not ghost atom; add grid
c
c
c     check if by atom grid parameters have been redifined by input
c
               
               lgridin=.false.
               if(ngridin.ne.0) then
               do jj=1,ngridin
                  if(grtagin(jj).eq.tag) then
                     lgridin=.true.
                     igridin=jj
                  endif
               enddo
            endif
!               write(0,*) ' jj tag ig ',jj,tag,igridin,lgridin
c     
c
               if(lgridin) then
                 nrad(itype) = nradin(igridin)
                  rad_cutoff(1,itype)=accqrin(igridin)
                  nang(itype) = nangin(igridin)
               if(leb) then
                  iangquad=max(iangquad,nangin(igridin))
                  nang(itype) = xc_leb_acc2nrd(nangin(igridin))
               endif
               elseif (grid_acc_index.gt.8)then
                  nrad(itype) = nradpts
                  if(leb) then
                     nang(itype) = xc_leb_acc2nrd(iangquad)
                  else
                     nang(itype) = nagrid
                  endif
               else
                  nrad(itype) = grid_params(grid_acc_index,1,
     &                 i_atomic_number)
                  if(grid_acc_index.ge.6.and.grid_acc_index.le.8) then
                     leb=.true.
                     iangat=grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
                     nang(itype) = xc_leb_acc2nrd(iangat)
                     iangquad=max(iangquad,iangat)
                  else
                  nang(itype) = grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
                  endif
               endif
!               write(0,*) ' it ',itype,nrad(itype),nang(itype)
c                
               if(nrad(itype).gt.mxnrad)mxnrad = nrad(itype)
               if(nang(itype).gt.mxnang)mxnang = nang(itype)
            else ! atomic number zero; charged ghost atom 
               nrad(itype) = 0
               nang(itype) = 0
            endif 
         endif
   50 continue
c
c     check mxnrad and mxnang
c
      if (mxnrad.gt.MXrad)then
         write(luout,*)'grid_acc_def: Too many radial shells.'
         call errquit(' grid_acc_def: raise MXrad',1)
      end if
      if (.not.leb.and.mxnang.gt.MXagrid)then
         write(LuOut,*)'grid_acc_def: Too many angular grids.'
         call errquit(' grid_acc_def: raise MXagrid',3)
      end if
c
      if (.not. rtdb_put(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_acc_def: rtdb_put failed', 119)
      if (.not. rtdb_put(rtdb, 'dft:mxnang', mt_int, 1, mxnang))
     &   call errquit('grid_acc_def: rtdb_put failed', 119)
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(6,*)' grid_acc_index, job_grid_acc: ',
     &               grid_acc_index, job_grid_acc
         write(6,*)' mxnrad, mxnang: ', 
     &               mxnrad, mxnang
         write(6,*)' nrad(ntypes) ',
     &              (nrad(itype),itype = 1, ntypes)
         write(6,*)' nang(ntypes) ',
     &              (nang(itype),itype = 1, ntypes)
      endif
      return
      end
