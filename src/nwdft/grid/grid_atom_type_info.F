      Subroutine grid_atom_type_info
c
c$Id: grid_atom_type_info.F,v 1.12 2003-12-13 02:43:28 marat Exp $
c
      implicit none
#include "errquit.fh"
c      
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      double precision BSrad(105)
      double precision ictr_coord(3), ictr_chg
      double precision autoang, angtoau, eps
      parameter (autoang = 0.529177249d0, angtoau = 1.d0/autoang)
      parameter (eps = 1.d-20)
      integer itype, ictr, iaz, i_atomic_number, icenter, jcenter
c
      integer lcoord, icoord, lcharge, icharge, ltags, itags,iptr
      logical same_atom, same_bq, isbq
c
      character*16 element
      character*16 tag
      character*2 symbol
c
      logical lnewtype
      logical atom_tag_check
      logical oprint, oprint_grid
c
      External atom_tag_check
c
c
c     Table of Bragg-Slater Atomic Radii (Angstroms)
c
c     Bragg-Slater radii: J.C. Slater, Symmetry and Energy Bands in Crystals,
c                         Dover, N.Y. 1972, page 55.
c                         The radii of noble gas atoms are set to be equal 
c                         to the radii of the corresponding halogen atoms.
c                         The radius of At is set to be equal to the radius of
c                         Po. 
c                         The radius of Fr is set to be equal to the radius of
c                         Cs. 
c
c                  H    He   Li   Be    B    C    N    O    F   Ne   
      Data BSrad/0.35,0.35,1.45,1.05,0.85,0.70,0.65,0.60,0.50,0.50,
c                  Na   Mg   Al   Si    P    S   Cl   Ar    K   Ca
     &           1.80,1.50,1.25,1.10,1.00,1.00,1.00,1.00,2.20,1.80,
c                  Sc   Ti    V   Cr   Mn   Fe   Co   Ni   Cu   Zn
     &           1.60,1.40,1.35,1.40,1.40,1.40,1.35,1.35,1.35,1.35,
c                  Ga   Ge   As   Se   Br   Kr   Rb   Sr    Y   Zr
     &           1.30,1.25,1.15,1.15,1.15,1.15,2.35,2.00,1.80,1.55,
c                  Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
     &           1.45,1.45,1.35,1.30,1.35,1.40,1.60,1.55,1.55,1.45,
c                  Sb   Te    I   Xe   Cs   Ba   La   Ce   Pr   Nd
     &           1.45,1.40,1.40,1.40,2.60,2.15,1.95,1.85,1.85,1.85,
c                  Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tm   Yb
     &           1.85,1.85,1.85,1.80,1.75,1.75,1.75,1.75,1.75,1.75,
c                  Lu   Hf   Ta    W   Re   Os   Ir   Pt   Au   Hg
     &           1.75,1.55,1.45,1.35,1.35,1.30,1.35,1.35,1.35,1.50,
c                  Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
     &           1.90,1.80,1.60,1.90,1.90,1.90,2.60,2.15,1.95,1.80,
c                  Pa    U   Np   Pu   Am   Cm   Bk   Cf   Es   Fm
     &           1.80,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,
c                  Md   No   Lr  Unq  Unp
     &           1.75,1.75,1.75,1.55,1.55/
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit(
     .     'grid_atom_type_info: failed to alloc coordinates',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit(
     '     'grid_atom_type_info: failed to alloc charges',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     .     'grid_atom_type_info: failed to alloc center tags',0, MA_ERR)
c
c     Get ncenter tags, coordinates, and charges from the geometry object.
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('gridatom_type_info: geom_cart_get failed',geom,
     &       GEOM_ERR)
c
c     generate number of atom types and atom type array iatype(icenter) 
c
      ntypes = 0
      do icenter = 1, ncenters
c
c        is this a new type of atom?
c
         lnewtype = .true.
         isbq = geom_isbq(geom,icenter)
         do jcenter = 1, icenter - 1
            same_atom = Dbl_MB(icharge + icenter - 1) .eq.
     &         Dbl_MB(icharge + jcenter - 1)
            same_bq = geom_isbq(geom,jcenter) .and. isbq
            same_atom = same_atom .or. same_bq
            if (same_atom .and.
     &          atom_tag_check(Byte_MB(itags + (icenter - 1)*16), 
     &                         Byte_MB(itags + (jcenter - 1)*16))
     &         )then   ! same atom type
               lnewtype = .false.
               iatype(icenter) = iatype(jcenter)
               goto 100
            endif
         enddo 
  100    continue
         if (lnewtype)then
            ntypes = ntypes + 1
            iatype(icenter) = ntypes
         endif 
      enddo 
      if (ntypes.gt.dft_ntags_bsmx)then
        write(LuOut,*) 'grid_atom_type_info:  Too many types of atoms.'
        call errquit(' grid_atom_type_info: raise dft_ntags_bsmx',2,
     &       GEOM_ERR)
      end if
c
c     set up type-indexed znuc array; znuc_atom_type
c
      do itype = 1, ntypes
         do icenter = 1, ncenters
            if (iatype(icenter) .eq. itype)then
c
c              center icenter is of type itype; assign charge
c
               znuc_atom_type(itype) = Dbl_MB(icharge + icenter - 1)
               goto 110 ! next type
            endif 
         enddo 
  110    continue 
      enddo 
c     
c     Define the atomic Bragg-Slater radii for each atom type.
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
     &            ('grid_atom_type_info: geom_cent_get failed', 0,
     &       GEOM_ERR)
               goto 40
            endif
         enddo
   40    continue
c        
         if (abs(znuc_atom_type(itype)).lt.EPS) then ! uncharged ghost atom; 
c
c           identify atom label following "bq"
c
            iptr=3
c hack for nbo
            if(tag(3:4).eq.'gh') iptr=5
            if(tag(1:1).eq.'X'.or.tag(1:1).eq.'x') iptr=2
            if (.not. geom_tag_to_element(tag(iptr:), symbol, 
     &           element, i_atomic_number)) then
               if (inp_compare(.false.,tag(1:2),'bq')) then 
                  i_atomic_number = 0
               elseif (inp_compare(.false.,tag(1:1),'X')) then 
                  i_atomic_number = 0
               else
                  call errquit(
     &            'grid_atom_type_info: non-bq center with zero charge',
     &             0, INPUT_ERR)
     &              
               endif
            endif
c
c
            if (i_atomic_number.eq.0)then
               bsrad_atom_type(itype) = EPS
            else
               bsrad_atom_type(itype) = BSrad(i_atomic_number)*ANGTOAU
            endif
         else    ! center is charged
c
c           no quadrature grids on charged ghost atoms 
c
            if (.not. geom_tag_to_element(tag, symbol, 
     &         element, i_atomic_number)) then
               if (symbol .ne. 'bq') call errquit
     &            ('grid_atom_type_info: center is neither atom nor bq',
     &              0, INPUT_ERR)
            endif
c
            if (i_atomic_number.ne.0)then ! not ghost atom
            ityp2ctr(itype)=ictr
c
               bsrad_atom_type(itype) = BSrad(i_atomic_number)*ANGTOAU
               if (bsrad_atom_type(itype).lt.EPS)then ! no radius found for atom
                  write(LuOut,*)' index ', 
     &               int(abs(znuc_atom_type(itype)) + EPS) 
                  write(LuOut,*)' BSR ',
     &               BSrad(int(abs(znuc_atom_type(itype)) + EPS))
                  write(LuOut,*)' grid_atom_type_info: ',
     &               ' Undefined atomic radius '
                  write(LuOut,*)' for atom type', itype
                  call errquit('Exiting in grid_atom_type_info.',1,
     &       UNKNOWN_ERR)
               endif
            else ! atomic number zero; charged ghost atom 
               bsrad_atom_type(itype) = EPS
            endif 
         endif
   50 continue
c
c     Build logical vector to tag centers as point charge centers or
c     real centers with basis functions
c
      do 55 ictr = 1, ncenters
         itype = iatype(ictr)
         if (bsrad_atom_type(itype).le.EPS)then
            iatype_pt_chg(ictr) = .true.
         else
            iatype_pt_chg(ictr) = .false.
         endif
   55 continue
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('grid_atom_type_info: pop stack failed.',0,
     &       MA_ERR)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('grid_atom_type_info: pop stack failed.',0,
     &       MA_ERR)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('grid_atom_type_info: pop stack failed.',0,
     &       MA_ERR)
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(LuOut,*)' iatype(ncenters) ',
     &              (iatype(ictr),ictr = 1, ncenters)
         write(LuOut,*)' iatype_pt_chg(ncenters) ',
     &              (iatype_pt_chg(ictr),ictr = 1, ncenters)
         write(LuOut,*)' bsrad_atom_type(ntypes) ',
     &              (bsrad_atom_type(itype),itype = 1, ntypes)
         write(LuOut,*)' znuc_atom_type(ntypes) ',
     &              (znuc_atom_type(itype),itype = 1, ntypes)
      endif
      return
      end
      logical function atom_tag_check(atom_tag_i, atom_tag_j)
c
      implicit none
c      
      character*16 atom_tag_i, atom_tag_j
c      write(*,*)' atom_tag_i = ', atom_tag_i
c      write(*,*)' atom_tag_j = ', atom_tag_j
      if (atom_tag_i .eq. atom_tag_j)then
         atom_tag_check = .true.
      else
         atom_tag_check = .false.
      endif      
      return
      end
