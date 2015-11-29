! $Id$
! ===========================================================================
!
!     Goal: create a complete NWChem interface to Molden.
!
!     1. need subroutine which prints MOs in batches (DONE)
!     2. need to dump basis set information (DONE)
!     3. need to dump geometry information (DONE)
!     4. (optional) dump scf and geometry optimization convergence information
!     5. (optional) dump normal modes from frequency calculation
!
!     This code was written on January 18th and 19th, 2009 by Jeff Hammond.
!     Improved 2015 by Vladimir Konjkov.
!
! ===========================================================================
!
#define MOLDEN_UNIT 37
!
! ===========================================================================
!
      subroutine nwmolden_reorder(din, dout, type, sphcart)
      implicit none
      double precision din(*)
      double precision dout(*)
      integer type              ! type (sp/s/p/d/..)
      integer sphcart           ! 0/1 for cartesian/spherical
!
      if (sphcart.eq.0) then
!       In MOLDEN file the following order of
!       cartesian functions is expected:
        if (type.eq.0) then
!        S functions:
          dout(1) = din(1)
        else if (type.eq.1) then
!        P functions:
!         x, y, z
          dout(1) = din(1)
          dout(2) = din(2)
          dout(3) = din(3)
        else if (type.eq.2) then
!        D functions:
!         xx, yy, zz, xy, xz, yz
          dout(1) = din(1)
          dout(2) = din(4)
          dout(3) = din(6)
          dout(4) = din(2)
          dout(5) = din(3)
          dout(6) = din(5)
        else if (type.eq.3) then
!        F functions:
!         xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
          dout(1) = din(1)
          dout(2) = din(7)
          dout(3) = din(10)
          dout(4) = din(4)
          dout(5) = din(2)
          dout(6) = din(2)
          dout(7) = din(6)
          dout(8) = din(9)
          dout(9) = din(8)
          dout(10) = din(5)
        else if (type.eq.4) then
!        G functions:
!         xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
!         xxyy xxzz yyzz xxyz yyxz zzxy
          dout(1) = din(1)
          dout(2) = din(11)
          dout(3) = din(15)
          dout(4) = din(2)
          dout(5) = din(3)
          dout(6) = din(7)
          dout(7) = din(12)
          dout(8) = din(10)
          dout(9) = din(14)
          dout(10) = din(4)
          dout(11) = din(6)
          dout(12) = din(13)
          dout(13) = din(5)
          dout(14) = din(8)
          dout(15) = din(9)
        end if
      else
!       In MOLDEN file the following order of
!       spherical functions is expected:
        if (type.eq.0) then
!        S functions:
          dout(1) = din(1)
        else if (type.eq.1) then
!        P functions:
!         P-1, P_0, P+1
          dout(1) = din(1)
          dout(2) = din(2)
          dout(3) = din(3)
        else if (type.eq.2) then
!        D functions:
!         D_0, D+1, D-1, D+2, D-2
          dout(1) = din(3)
          dout(2) = -din(4)
          dout(3) = din(2)
          dout(4) = din(5)
          dout(5) = din(1)
         else if (type.eq.3) then
!        F functions:
!         F 0, F+1, F-1, F+2, F-2, F+3, F-3
          dout(1) = din(4)
          dout(2) = -din(5)
          dout(3) = din(3)
          dout(4) = din(6)
          dout(5) = din(2)
          dout(6) = -din(7)
          dout(7) = din(1)
        else if (type.eq.4) then
!        G functions:
!         G_0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
          dout(1) = din(5)
          dout(2) = -din(6)
          dout(3) = din(4)
          dout(4) = din(7)
          dout(5) = din(3)
          dout(6) = -din(8)
          dout(7) = din(2)
          dout(8) = din(9)
          dout(9) = din(1)
        end if
      end if
      return
      end
!
! ===========================================================================
!
      integer function nwmolden_nbusfunc(type,sphcart)
      implicit none
      integer type              ! [input] type (sp/s/p/d/..)
      integer sphcart           ! [input] 0/1 for cartesian/spherical
      if (sphcart.eq.1) then
          nwmolden_nbusfunc = 2 * type + 1
      else
          nwmolden_nbusfunc = (type + 1) * (type + 2) / 2
      end if
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_openfile()
      implicit none
#include "errquit.fh"
      character*64 filename
!
      call util_file_prefix('molden', filename )
!
      if (MOLDEN_UNIT.gt.9) then
         open(unit=MOLDEN_UNIT,
     1        file=filename,
     2        status='unknown',
     3        form='formatted',
     4        err=911)
      endif
!
      write(MOLDEN_UNIT,100) '[Molden Format]'
 100  format(a15)
!
      nwmolden_openfile = .true.
      return
!
 911  nwmolden_openfile = .false.
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_closefile()
      implicit none
!
      if (MOLDEN_UNIT.gt.9) then
        close(unit=MOLDEN_UNIT,
     1        status='keep',
     2        err=911)
      endif
!
      nwmolden_closefile = .true.
      return
!
 911  nwmolden_closefile = .false.
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_geom(geom)
      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer geom              ! [input] geometry handler
      integer j
      integer iat               ! index of atom
      integer atn               ! atomic number
      integer natoms            ! no. of atoms
      double precision coord(3)
      double precision qnuc
      character*2 symbol
      character*16 tag
      character*32 element
      character*80 buf
      logical status
!
! this might be wrong if something other than atomic units
! are specified in the NWChem geometry input, etc.
!
      write(MOLDEN_UNIT,100) '[Atoms] AU'
 100  format(a10)
!
      status = geom_check_handle(geom, 'nwmolden_geom')
      if (.not.status)
     1   call errquit('nwmolden_geom: error code = ',1,GEOM_ERR)
!
      status = geom_ncent(geom, natoms)
      if (.not.status)
     1   call errquit('nwmolden_geom: error code = ',2,GEOM_ERR)
!
      do iat=1,natoms
!
         status = geom_cent_get(geom, iat, tag, coord, qnuc)
         if (.not.status)
     1      call errquit('nwmolden_geom: error code = ',3,GEOM_ERR)
!
         status = geom_tag_to_element(tag, symbol, element, atn)
         if (.not.status) symbol = 'X'
!
         write(MOLDEN_UNIT,300) symbol,iat,nint(qnuc),(coord(j),j=1,3)
 300     format(1x,a4,1x,i4,1x,i4,1x,3f15.8)
!
      enddo
!
      nwmolden_geom = .true.
!
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_basis(basis)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "errquit.fh"
      integer basis             ! [input] basis handle
      integer geom              ! geometry handler
      integer iatom             ! index of atom
      integer iprimo            ! index of primitive
      integer icont             ! contraction index
      integer ifirst, ilast     ! first and last contraction index
      integer type              ! type (sp/s/p/d/..)
      integer natoms            ! no. of atoms
      integer nprimo            ! no. of primitives
      integer ngeno             ! no. of contractions
      integer sphcart           ! 0/1 for cartesian/spherical
      logical status
!
      integer max_prim
      parameter (max_prim = 4000)
      double precision exp(max_prim)
      double precision coeff(max_prim)
      character*3 ctype(0:6),cltype(2)
!
      ctype(0)='s'
      ctype(1)='p'
      ctype(2)='d'
      ctype(3)='f'
      ctype(4)='g'
      ctype(5)='h'
      ctype(6)='i'
      cltype(1)='sp'
      cltype(2)='spd'
!
      write(MOLDEN_UNIT,400) '[GTO]'
 400  format(a5)
!
      status = bas_check_handle(basis,'nwmolden_basis')
      if (.not.status)
     1   call errquit('nwmolden_basis: error code = ',1,BASIS_ERR)
!
      status = bas_geom(basis, geom)
      if (.not.status)
     1   call errquit('nwmolden_basis: error code = ',2,BASIS_ERR)
!
      status = geom_ncent(geom, natoms)
      if (.not.status)
     1   call errquit('nwmolden_basis: error code = ',3,BASIS_ERR)
!
      do iatom = 1, natoms
        status = bas_ce2cnr(basis,iatom,ifirst,ilast)
        if (.not.status)
     1      call errquit('nwmolden_basis: error code = ',4,BASIS_ERR)
        write(MOLDEN_UNIT, 450) iatom, 0
450     format(1x,i4,1x,i2)
!
        do icont = ifirst, ilast
!
         status = bas_continfo(basis,icont,type,nprimo,ngeno,sphcart)
         if (.not.status)
     1      call errquit('nwmolden_basis: error code = ',5,BASIS_ERR)
!
         if (nprimo.gt.max_prim)
     1      call errquit('nwmolden_basis: nprimo > max_prim',
     2                   max_prim-nprimo,BASIS_ERR)
!
         if (type.ge.0) then
           write(MOLDEN_UNIT,500) ctype(type),nprimo,0
         elseif (type.lt.0) then
           write(MOLDEN_UNIT,500) cltype(abs(type)),nprimo,0
         endif
 500     format(1x,a4,1x,i4,1x,i2)
!
         status = bas_get_exponent(basis,icont,exp)
         if (.not.status)
     1      call errquit('nwmolden_basis: error code = ',6,BASIS_ERR)

         status = bas_get_coeff(basis,icont,coeff)
         if (.not.status)
     1      call errquit('nwmolden_basis: error code = ',7,BASIS_ERR)
!
         do iprimo = 1, nprimo
            write(MOLDEN_UNIT,700) exp(iprimo),coeff(iprimo)
 700        format(1x,f14.8,1x,f14.8)
         enddo
        enddo
        write(MOLDEN_UNIT,*)
      enddo
!
      nwmolden_basis = .true.
!
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_movecs(basis,nbf,nmo,
     1                                 eval,occ,sym,g_movecs,scftype)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
      integer nbf               ! [input] number of AOs
      integer nmo               ! [input] number of MOs
      double precision eval(*)  ! [input] MO energies
      double precision occ(*)   ! [input] MO occupation numbers
      integer sym(*)            ! [input] MO symmetries
      integer g_movecs          ! [input] GA with vectors
      character*(*) scftype     ! [input] RHF, ROHF, UHF, DFT, or MCSCF
      integer numcont           ! number of mapped contractions
      integer icont             ! contraction index
      integer type              ! type (sp/s/p/d/..)
      integer nprimo            ! no. of primitives
      integer ngeno             ! no. of contractions
      integer ibusfunc          ! basis func with diffrent m
      integer sphcart           ! 0/1 for cartesian/spherical
      integer k_vec, l_vec
      double precision reordered_mo(nbf)
      integer alfabeta, spinref
      integer orb, aof, offset
      integer basis, geom
      integer maxop, maxireps
      parameter (maxop = 120, maxireps=20)
      integer nop, nir,  class_dim(maxireps)
      double precision chars(maxireps*maxireps)
      character*5 cspin, cspins(2)
      character*5 csym
      character*8 zir(maxireps), zclass(maxireps)
      character*20 zname
      logical status
      integer nbusfunc
      integer nwmolden_nbusfunc ! function defined in this file
! external function declarations
      logical sym_char_table
      external sym_char_table
      logical bas_geom
      external bas_geom
!
      status = bas_geom(basis, geom)
      if (.not.status)
     1    call errquit('nwmolden_movecs: error code = ',1,BASIS_ERR)
!
      call sym_group_name(geom, zname)
!
      status = sym_char_table(zname,nop,nir,class_dim,zir,zclass,chars)
      if (.not.status)
     1    call errquit('nwmolden_movecs: error code = ',2,BASIS_ERR)
!
      cspins(1) = 'Alpha'
      cspins(2) = 'Beta '
!
      if (scftype.eq.'RHF')   spinref = 1
      if (scftype.eq.'ROHF')  spinref = 2
      if (scftype.eq.'UHF')   spinref = 2
      if (scftype.eq.'DFT')   spinref = 1 ! ambiguous
      if (scftype.eq.'MCSCF') spinref = 1
!
!     allocate local memory
!
!     TODO: this is not scalable and probably will not work in DDSCF scheme
!
      if (.not.ma_push_get(mt_dbl,nbf,'vector',l_vec,k_vec))
     1    call errquit('nwmolden_movecs: MA alloc. failed',nbf,MA_ERR)
!
      status = bas_numcont(basis, numcont)
      if (.not.status)
     1    call errquit('nwmolden_movecs: error code = ',3,BASIS_ERR)
!
      status = bas_continfo(basis,1,type,nprimo,ngeno,sphcart)
!
      if (sphcart.eq.1) then
        write(MOLDEN_UNIT,50) '[5D]'
        write(MOLDEN_UNIT,50) '[9G]'
      end if
      write(MOLDEN_UNIT,50) '[MO]'
  50  format(a4)
!
      do alfabeta = 1, spinref
         cspin = cspins(alfabeta)
         do orb = 1, nmo
!
!          pull the MOs from a GA
!
           call ga_get(g_movecs, 1, nbf, orb, orb, dbl_mb(k_vec), 1)
!
!          reorder MO
!
           aof = 1
!
           do icont = 1, numcont
              status = bas_continfo(basis,icont,type,
     1                              nprimo,ngeno,sphcart)
              if (.not.status)
     1            call errquit('nwmolden_movecs: error code = ',
     2                         4,BASIS_ERR)
              nbusfunc = nwmolden_nbusfunc(type, sphcart)
              call nwmolden_reorder(dbl_mb(k_vec + aof - 1),
     1                              reordered_mo(aof),
     2                              type, sphcart)
              aof = aof + nbusfunc
           enddo
!
!          print stuff
!
           csym = zir(sym(orb))
           write(MOLDEN_UNIT,100) csym,eval(orb),cspin,occ(orb)
 100       format('Sym= ', a12/
     1            'Ene= ', f12.5/
     2            'Spin= ', a11/
     3            'Occup= ', f10.5)
           do aof = 1, nbf
!              write(MOLDEN_UNIT,200) aof,dbl_mb(k_vec + aof - 1)
              write(MOLDEN_UNIT,200) aof,reordered_mo(aof)
 200          format(i6,f18.12)
           enddo ! aof
        enddo    ! orb
      enddo      ! alfabeta
!
!     de-allocate local memory
!
      if (.not. ma_pop_stack(l_vec))
     1    call errquit('nwmolden_movecs: MA stack is corrupt',5,MA_ERR)
!
      nwmolden_movecs = .true.
!
      return
      end
!
! ===========================================================================
!

