      Subroutine grid_acc_def(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
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
      parameter (num_gridtypes = 8)
      logical oprint, oprint_grid,lgridin,qmmm_init
      logical grid_tagdefd
      external grid_tagdefd
      character*2 symbol
      character*15 gridtypes(num_gridtypes)
      double precision grid_acc(num_gridtypes)
      double precision job_grid_acc
      character*16 element
      character*16 tag, tag_mod
      integer lenb
      double precision eps
      parameter (eps = 1.d-20)
      data gridtypes /'xcoarse', 'coarse', 'medium', 'fine', 'xfine',
     H     'huge',
     .                'gausleg','lebedev'/
      data grid_acc /
     .     1d-4,1d-5,1d-6,1d-7,1d-8,1d-10,
     .     1d-8,1d-8/
c
      double precision ictr_coord(3), ictr_chg
      integer itype, ictr, iaz, i_atomic_number,iangat
      integer grptr
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
c     huge    = 1d-10
c
      integer  grid_params(6,2,2,105)
      integer i,j,k,l,iptr
      integer xc_leb_acc2nrd
      external xc_leb_acc2nrd
      integer xc_dau_acc2nrd
      external xc_dau_acc2nrd
* equivalence must come before data statements
c
c      It appears that H (and so likely He as well) can't be reduced as much as
c      previously thought.  Set these to Li - Ne values.
c
      data ((((grid_params(i,j,k,l),i=1,6),j=1,2),k=1,2),l=1,40)/ 
     .     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 45, 60,100,300,  8, 11, 13, 14, 17, 17,!  H (1)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 17, 17,! He (2)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,! Li (3)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,! Be (4)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,!  B (5)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,!  C (6)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,!  N (7)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,!  O (8)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,!  F (9)
     &     30, 50, 70,100,140,300, 12, 15, 18, 24, 34, 34,
     .     21, 35, 49, 70,100,300,  8, 11, 13, 14, 18, 18,! Ne (10)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 13, 15, 18, 18,! Na (11)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,! Mg (12)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,! Al (13)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,! Si (14)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,!  P (15)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,!  S (16)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,! Cl (17)
     &     45, 75, 95,125,175,300, 12, 18, 24, 30, 44, 34,
     .     42, 70, 88,123,125,300,  8, 11, 14, 15, 18, 18,! Ar (18)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,!  K (19)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Ca (20)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Sc (21)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Ti (22)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,!  V (23)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Cr (24)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Mn (25)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Fe (26)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Co (27)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Ni (28)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Cu (29)
     &     75, 95,110,160,210,400, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160,400,  8, 11, 14, 16, 18, 18,! Zn (30)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750, 8, 11, 14, 16, 18, 18,! Ga (31)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750, 8, 11, 14, 16, 18, 19,! Ge (32)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750,  8, 11, 14, 16, 18, 18,! As (33)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750,  8, 11, 14, 16, 18, 18,! Se (34)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750,  8, 11, 14, 16, 18,  18,! Br (35)
     &     75, 95,110,160,210, 750, 14, 22, 30, 34, 38, 38,
     .     75, 95,112,140,160, 750,  8, 11, 14, 16, 18, 18,! Kr (36)
     &    105,130,155,205,235, 400, 16, 20, 32, 44, 48, 48,
     .     84,104,123,155,205, 400, 8, 11, 14, 16, 18, 18,! Rb (37)
     &    105,130,155,205,235, 400,16, 20, 32, 44, 48, 48,
     .     84,104,123,141,205, 400, 8, 11, 14, 16, 18, 18,! Sr (38)
     &    105,130,155,205,235, 400,16, 20, 32, 44, 48, 48,
     .     84,104,123,141,205, 400, 8, 11, 14, 16, 18, 18,!  Y (39)
     &    105,130,155,205,235, 400,16, 20, 32, 44, 48, 48,
     .     84,104,123,141,205, 400, 8, 11, 14, 16, 18, 18/! Zr (40)
      data ((((grid_params(i,j,k,l),i=1,6),j=1,2),k=1,2),l=41,80)/ 
     .     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18,  18,! Nb (41)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Mo (42)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Tc (43)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ru (44)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Rh (45)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Pd (46)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ag (47)
     &     105,130,155,205,235, 400, 16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205, 400, 8, 11, 14, 16, 18, 18,! Cd (48)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! In (49)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Sn (50)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Sb (51)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Te (52)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,!  I (53)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Xe (54)
     .     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Cs (55)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ba (56)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! La (57)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ce (58)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Pr (59)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Nd (60)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Pm (61)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Sm (62)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Eu (63)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Gd (64)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Tb (65)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Dy (66)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ho (67)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Er (68)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Tm (69)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11,1 4, 16, 18, 18,! Yb (70) 
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Lu (71)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Hf (72)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ta (73)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,!  W (74)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Re (75)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Os (76)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ir (77)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Pt (78)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Au (79)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18/! Hg (80)
      data ((((grid_params(i,j,k,l),i=1,6),j=1,2),k=1,2),l=81,105)/ 
     .     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11,1 4, 16, 18, 18, ! Tl (81)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Pb (82)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Bi (83)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Po (84)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! At (85)
     &     105,130,155,205,235,950,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,950, 8, 11, 14, 16, 18, 18,! Rn (86)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Fr (87)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 16, 18, 18,! Ra (88)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Ac (89)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Th (90)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Pa (91)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,155,205,400, 8, 11, 14, 17, 18, 18,!  U (92)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Np (93)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Pu (94)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Am (95)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Cm (96)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Bk (97)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Cf (98)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Es (99)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Fm (100)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Md (101)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! No (102)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Lr (103)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18,! Unq(104)
     &     105,130,155,205,235,400,16, 20, 32, 44, 48, 48,
     .      84,104,123,141,205,400, 8, 11, 14, 17, 18, 18/! Unp(105)
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
c     detect qmmm job
c
      if (.not. rtdb_get(rtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init=.false.
      if (.not. rtdb_get(rtdb,'dft:bqdontcare',mt_log,1,bqdontcare))
     $     bqdontcare=.false.
      bqdontcare=bqdontcare.or.qmmm_init

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
     &        call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dft:nradin',mt_int,ngridin,nradin))
     &        call errquit('grid_inp: rtdb_get failed', 1000, RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dft:nangin',mt_int,ngridin,nangin))
     &        call errquit('grid_inp: rtdb_get failed', 1000, RTDB_ERR)
      endif
      if (.not. inp_match(num_gridtypes, .false., gridtype, gridtypes, 
     &   grid_acc_index))then
         write(LuOut,*)' Unknown gridtype; defaulting to medium grid'
         gridtype='medium'
         grid_acc_index = 8
      endif
c
c     see if we have lebedev
c
      if (.not. rtdb_get(rtdb, 'dft:lebedev', mt_log, 1, leb))
     &   call errquit(' graccdef: grid_inp didnt write leb',0, RTDB_ERR)

      goto (8, 8, 8, 8, 8, 8, 6, 7) grid_acc_index
      call errquit('grid_acc_def: unknown gridtype', grid_acc_index,
     &       INPUT_ERR)
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
      leb=.false.
      goto 8
    7 continue
c
c     lebedev grid;
c
      leb=.true.
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1, iangquad))then
         iangquad = 8
      else
         if (iangquad.lt.1)then
          if(ga_nodeid().eq.0) then
           write(LuOut,*)' IANGQUAD =',iangquad
           write(LuOut,*)' IANGQUAD must be in the range: 1-',IAC_MAX,
     *     '; resetting. '
          end if
          iangquad = 8
         endif
         if (iangquad.gt.IAC_MAX)then
           if(ga_nodeid().eq.0) then
            write(LuOut,*)' IANGQUAD =',iangquad
            write(LuOut,*)' IANGQUAD must be in the range: 1-',IAC_MAX,
     *     '; resetting. '
           end if
           iangquad = 8
         endif
      endif
      if (.not. rtdb_get(rtdb,'dft:nradpts', mt_int, 1, nradpts)) then
         gridtype='medium'
         grid_acc_index = 3
      endif      
      goto 8
    8 continue
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     .     job_grid_acc = grid_acc(grid_acc_index)
      if(leb) then
         grptr=2
      else
         grptr=1
      endif
c
c     shove job_grid_acc into the rtdb
c
      if (.not. rtdb_put(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('grid_acc_def: rtdb_put failed', 911, RTDB_ERR)
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
     &            ('grid_acc_def: geom_cent_get failed', 0, GEOM_ERR)
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
            iptr=3
c hack for nbo
            if(tag(3:4).eq.'gh') iptr=5
            if(tag(1:1).eq.'X'.or.tag(1:1).eq.'x') iptr=2
            lenb=inp_strlen(tag)
            tag_mod=' '
            tag_mod(1:lenb-iptr+1)=tag(iptr:lenb)
            if (.not. geom_tag_to_element(tag_mod, symbol, 
     &           element, i_atomic_number)) then
               if (inp_compare(.false.,tag(1:2),'bq')) then 
                  i_atomic_number = 0
               elseif (inp_compare(.false.,tag(1:1),'X')) then 
                  i_atomic_number = 0
               else
                  call errquit
     &              ('grid_acc_def: non-bq center with zero charge', 0,
     &       INPUT_ERR)
               endif
            endif
            if (i_atomic_number.eq.0)then
c
c     we do not bail out ONLY  when X atom (for Zmatrix)
c     or when we are in a QMMM job
c
               if(iptr.ne.2.and.(.not.bqdontcare)) then
               if(ga_nodeid().eq.0) then
                  write(luout,*) ' *********************************'
                  write(luout,*) ' * INPUT ERROR                   *'
                  write(luout,*) ' * please define bq+element name *'
                  write(luout,*) ' * in the geometry section       *'
                  write(luout,*) ' * e.g. bqSi                     *'
                  write(luout,*) ' *********************************'
                  call util_flush(luout)
               endif
                    call errquit(
     X              ' bq DFT atom with no atom type ',0,INPUT_ERR)
               endif
               if(bqdontcare) then
                  nrad(itype) = grid_params(min(grid_acc_index,6),
     A                 1,grptr, 1)
                  iangat=grid_params(min(grid_acc_index,6),
     A                 2,grptr,       1)
                  nang(itype) = xc_leb_acc2nrd(iangat)
                  iangquad=max(iangquad,iangat)
               else
                  nrad(itype) = 0
                  nang(itype) = 0
               endif
            else
               if (grid_acc_index.eq.7)then
                  nrad(itype) = nradpts
                  nang(itype) = nagrid
               else
c bq handling
               lgridin=grid_tagdefd(tag(3:),ngridin,grtagin,
     I              igridin)
c
               if(lgridin) then
                  nrad(itype) = nradin(igridin)
                  iangat=nangin(igridin)
               else
cfix case for generic bq input
                  if(grid_acc_index.gt.6)then
                     nrad(itype) = nradpts
                     iangat=iangquad
                  else
                     nrad(itype) = grid_params(min(grid_acc_index,6),
     A                    1,grptr, i_atomic_number)
                     iangat=grid_params(min(grid_acc_index,6),
     A                    2,grptr,       i_atomic_number)
                  endif
                  nang(itype) = xc_leb_acc2nrd(iangat)
                  iangquad=max(iangquad,iangat)
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
     &              ('grid_acc_def: center is neither atom nor bq', 0,
     &       INPUT_ERR)
            endif
c
            if (i_atomic_number.ne.0)then ! not ghost atom; add grid
c
c
c     check if by atom grid parameters have been redifined by input
c
               lgridin=grid_tagdefd(tag,ngridin,grtagin,
     I              igridin)
c     
c
               if(lgridin) then
                 nrad(itype) = nradin(igridin)
cedo                  rad_cutoff(1,itype)=accqrin(igridin)
                  nang(itype) = nangin(igridin)
               if(leb) then
                  iangquad=max(iangquad,nangin(igridin))
                  nang(itype) = xc_leb_acc2nrd(nangin(igridin))
               endif
               elseif (grid_acc_index.gt.6)then
                  nrad(itype) = nradpts
                  if(leb) then
                     nang(itype) = xc_leb_acc2nrd(iangquad)
                  else
                     nang(itype) = nagrid
                  endif
               else
                  nrad(itype) = grid_params(grid_acc_index,1,grptr,
     &                 i_atomic_number)
                  if(leb) then
                     iangat=grid_params(grid_acc_index,2,2,
     &                                    i_atomic_number)
                     nang(itype) = xc_leb_acc2nrd(iangat)
                     iangquad=max(iangquad,iangat)
                  else
                     nang(itype) = grid_params(grid_acc_index,2,1,
     &                                    i_atomic_number)
                  endif
               endif
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
         call errquit(' grid_acc_def: raise MXrad',1,
     &       UNKNOWN_ERR)
      end if
      if((.not.leb).and.mxnang.gt.MXagrid)then
         write(LuOut,*)'grid_acc_def: Too many angular grids.'
         call errquit(' grid_acc_def: raise MXagrid',3, UNKNOWN_ERR)
      end if
c
      if (.not. rtdb_put(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_acc_def: rtdb_put failed', 119, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:mxnang', mt_int, 1, mxnang))
     &   call errquit('grid_acc_def: rtdb_put failed', 119, RTDB_ERR)
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(LuOut,*)' grid_acc_index, job_grid_acc: ',
     &               grid_acc_index, job_grid_acc
         write(LuOut,*)' mxnrad, mxnang: ', 
     &               mxnrad, mxnang
         write(LuOut,*)' nrad(ntypes) ',
     &              (nrad(itype),itype = 1, ntypes)
         write(LuOut,*)' nang(ntypes) ',
     &              (nang(itype),itype = 1, ntypes)
      endif
      return
      end
      logical function grid_tagdefd(tag,ngridin,grtagin,
     I     igridin)
      implicit none
#include "inp.fh"
      character*(*) tag,grtagin(*)
      integer ngridin
      integer igridin
c
      integer jj
c      
      grid_tagdefd=.false.
      if(ngridin.ne.0) then
         do jj=1,ngridin
            if(inp_compare(.false.,grtagin(jj),tag)) then
               grid_tagdefd=.true.
               igridin=jj
            endif
         enddo
      endif
      return
      end
