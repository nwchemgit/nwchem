      Subroutine grid_acc_def(rtdb)
c
c$Id: grid_acc_def.F,v 1.7 1999-07-13 21:27:50 edo Exp $
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
      parameter (num_gridtypes = 7)
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
     &                'gausleg', 'lebedev'/
      data grid_acc /1d-4,1d-5,1d-6,1d-7,1d-8,1d-8,1d-8/
c
      character*10 radgrids(5)
      data radgrids /'eumac','murakn','tarad','gausslag','gcinterv'/
      double precision ictr_coord(3), ictr_chg
      integer itype, ictr, iaz, i_atomic_number
      double precision R_MAX
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
      integer grid_params_a
      integer grid_params_b
      dimension grid_params(5,2,105)
      dimension grid_params_a(5,2,54)
      dimension grid_params_b(5,2,51)
      integer xc_leb_acc2nrd
      external xc_leb_acc2nrd
      equivalence(grid_params(1,1,1), grid_params_a(1,1,1))
      equivalence(grid_params(1,1,55), grid_params_b(1,1,1))
* equivalence must come before data statements
c
c      It appears that H (and so likely He as well) can't be reduced as much as
c      previously thought.  Set these to Li - Ne values.
c
c      Data grid_params_a/ 24, 28, 36, 48, 58,  5,  7,  9, 13, 15, !  H (1)
c     &                    26, 36, 46, 60, 80,  7, 10, 12, 17, 23, ! He (2)
      Data grid_params_a/ 30, 50, 70,100,140, 12, 15, 18, 24, 34, !  H (1)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, ! He (2)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, ! Li (3)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, ! Be (4)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, !  B (5)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, !  C (6)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, !  N (7)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, !  O (8)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, !  F (9)
     &                    30, 50, 70,100,140, 12, 15, 18, 24, 34, ! Ne (10)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Na (11)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Mg (12)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Al (13)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Si (14)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, !  P (15)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, !  S (16)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Cl (17)
     &                    45, 75, 95,125,175, 12, 18, 24, 30, 44, ! Ar (18)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, !  K (19)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Ca (20)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Sc (21)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Ti (22)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, !  V (23)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Cr (24)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Mn (25)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Fe (26)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Co (27)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Ni (28)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Cu (29)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Zn (30)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Ga (31)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Ge (32)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! As (33)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Se (34)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Br (35)
     &                    75, 95,110,160,210, 14, 22, 30, 34, 38, ! Kr (36)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Rb (37)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Sr (38)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, !  Y (39)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Zr (40)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Nb (41)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Mo (42)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Tc (43)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ru (44)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Rh (45)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pd (46)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ag (47)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Cd (48)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! In (49)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Sn (50)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Sb (51)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Te (52)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, !  I (53)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48/ ! Xe (54)
      Data grid_params_b/105,130,155,205,235, 16, 20, 32, 44, 48, ! Cs (55)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ba (56)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! La (57)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ce (58)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pr (59)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Nd (60)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pm (61)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Sm (62)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Eu (63)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Gd (64)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Tb (65)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Dy (66)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ho (67)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Er (68)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Tm (69)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Yb (70) 
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Lu (71)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Hf (72)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ta (73)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, !  W (74)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Re (75)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Os (76)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ir (77)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pt (78)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Au (79)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Hg (80)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Tl (81)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pb (82)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Bi (83)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Po (84)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! At (85)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Rn (86)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Fr (87)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ra (88)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Ac (89)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Th (90)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pa (91)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, !  U (92)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Np (93)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Pu (94)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Am (95)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Cm (96)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Bk (97)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Cf (98)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Es (99)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Fm (100)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Md (101)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! No (102)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Lr (103)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48, ! Unq(104)
     &                   105,130,155,205,235, 16, 20, 32, 44, 48/ ! Unp(105)
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
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
      goto (8, 8, 8, 8, 8, 6, 7) grid_acc_index
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
      job_grid_acc = grid_acc(grid_acc_index)
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
               else
                  call errquit
     &              ('grid_acc_def: non-bq center with zero charge', 0)
               endif
            endif
            if (i_atomic_number.eq.0)then
               nrad(itype) = 0
               nang(itype) = 0
            else
               if (grid_acc_index.gt.5)then
                  nrad(itype) = nradpts
                  nang(itype) = nagrid
               else
                  nrad(itype) = grid_params(grid_acc_index,1,
     &                                    i_atomic_number)
                  nang(itype) = grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
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
               elseif (grid_acc_index.gt.5)then
                  nrad(itype) = nradpts
                  if(leb) then
                     nang(itype) = xc_leb_acc2nrd(iangquad)
                  else
                     nang(itype) = nagrid
                  endif
               else
                  nrad(itype) = grid_params(grid_acc_index,1,
     &                                    i_atomic_number)
                  nang(itype) = grid_params(grid_acc_index,2,
     &                                    i_atomic_number)
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
