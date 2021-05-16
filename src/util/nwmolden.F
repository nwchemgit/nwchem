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
c           1   2   3   4   5   6
c molden : xx, yy, zz, xy, xz, yz
c nwchem : xx, xy, xz, yy, yz, zz
          dout(1) = din(1)
          dout(2) = din(4)
          dout(3) = din(6)
c conversion factors fron moldens gaussian.f          
          dout(4) = din(2)/sqrt(3d0)
          dout(5) = din(3)/sqrt(3d0)
          dout(6) = din(5)/sqrt(3d0)
        else if (type.eq.3) then
!     F functions:
c            1    2    3    4    5    6    7    8    9   10
c molden : xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
c nwchem : xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz           
          dout(1) = din(1)
          dout(2) = din(7)
          dout(3) = din(10)
c conversion factors fron moldens gaussian.f          
          dout(4) = din(4)/sqrt(5d0)
          dout(5) = din(2)/sqrt(5d0)
          dout(6) = din(3)/sqrt(5d0)
          dout(7) = din(6)/sqrt(5d0)
          dout(8) = din(9)/sqrt(5d0)
          dout(9) = din(8)/sqrt(5d0)
          dout(10) = din(5)/sqrt(5d0)/sqrt(3d0)
        else if (type.eq.4) then
!     G functions:
c           1    2    3    4    5    6    7    8    9
cmolden  xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy
cnwchem  xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz
c          10   11   12   13   14   15
cmolden  xxyy xxzz yyzz xxyz yyxz zzxy
cnwchem  xzzz yyyy yyyz yyzz yzzz zzzz
          dout(1) = din(1)
          dout(2) = din(11)
          dout(3) = din(15)
c conversion factors fron moldens gaussian.f          
          dout(4) = din(2)/sqrt(7d0)
          dout(5) = din(3)/sqrt(7d0)
          dout(6) = din(7)/sqrt(7d0)
          dout(7) = din(12)/sqrt(7d0)
          dout(8) = din(10)/sqrt(7d0)
          dout(9) = din(14)/sqrt(7d0)
          dout(10) = din(4)/sqrt(7d0)/sqrt(5d0)*sqrt(3d0)
          dout(11) = din(6)/sqrt(7d0)/sqrt(5d0)*sqrt(3d0)
          dout(12) = din(13)/sqrt(7d0)/sqrt(5d0)*sqrt(3d0)
          dout(13) = din(5)/sqrt(7d0)/sqrt(5d0)
          dout(14) = din(8)/sqrt(7d0)/sqrt(5d0) 
          dout(15) = din(9)/sqrt(7d0)/sqrt(5d0)
       else if (type.gt.4) then
          call errquit(' molden interface support up to g functions',
     M         type,0)
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
       else if (type.gt.4) then
          call errquit(' molden interface support up to g functions',
     M         type,0)
        end if
      end if
      return
      end
      subroutine moldennw_reorder(din, dout, type, sphcart)
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
c           1   2   3   4   5   6
c molden : xx, yy, zz, xy, xz, yz
c nwchem : xx, xy, xz, yy, yz, zz
          dout(1) = din(1)
          dout(4) = din(2)
          dout(6) = din(3)
c conversion factors fron moldens gaussian.f          
          dout(2) = din(4)*sqrt(3d0)
          dout(3) = din(5)*sqrt(3d0)
          dout(5) = din(6)*sqrt(3d0)
        else if (type.eq.3) then
!     F functions:
c            1    2    3    4    5    6    7    8    9   10
c molden : xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
c nwchem : xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz           
          dout(1)  = din(1)
          dout(7)  = din(2)
          dout(10) = din(3)
c conversion factors fron moldens gaussian.f          
          dout(4)  = din(4 )*sqrt(5d0)
          dout(2)  = din(5 )*sqrt(5d0)
          dout(3)  = din(6 )*sqrt(5d0)
          dout(6)  = din(7 )*sqrt(5d0)
          dout(9)  = din(8 )*sqrt(5d0)
          dout(8)  = din(9 )*sqrt(5d0)
          dout(5)  = din(10)*sqrt(5d0)*sqrt(3d0)
        else if (type.eq.4) then
!     G functions:
c           1    2    3    4    5    6    7    8    9
cmolden  xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy
cnwchem  xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz
c          10   11   12   13   14   15
cmolden  xxyy xxzz yyzz xxyz yyxz zzxy
cnwchem  xzzz yyyy yyyz yyzz yzzz zzzz
          dout(1 ) = din(1 )
          dout(11) = din(2 )
          dout(15) = din(3 )
c conversion factors fron moldens gaussian.f          
          dout(2 ) = din(4 )*sqrt(7d0)
          dout(3 ) = din(5 )*sqrt(7d0)
          dout(7 ) = din(6 )*sqrt(7d0)
          dout(12) = din(7 )*sqrt(7d0)
          dout(10) = din(8 )*sqrt(7d0)
          dout(14) = din(9 )*sqrt(7d0)
          dout(4 ) = din(10)*sqrt(7d0)*sqrt(5d0)/sqrt(3d0)
          dout(6 ) = din(11)*sqrt(7d0)*sqrt(5d0)/sqrt(3d0)
          dout(13) = din(12)*sqrt(7d0)*sqrt(5d0)/sqrt(3d0)
          dout(5 ) = din(13)*sqrt(7d0)*sqrt(5d0)
          dout(8 ) = din(14)*sqrt(7d0)*sqrt(5d0) 
          dout(9 ) = din(15)*sqrt(7d0)*sqrt(5d0)
       else if (type.gt.4) then
          call errquit(' molden interface support up to g functions',
     M         type,0)
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
          dout(3) =  din(1)
          dout(4) = -din(2)
          dout(2) =  din(3)
          dout(5) =  din(4)
          dout(1) =  din(5)
         else if (type.eq.3) then
!        F functions:
!         F 0, F+1, F-1, F+2, F-2, F+3, F-3
          dout(4) =  din(1)
          dout(5) = -din(2)
          dout(3) =  din(3)
          dout(6) =  din(4)
          dout(2) =  din(5)
          dout(7) = -din(6)
          dout(1) =  din(7)
        else if (type.eq.4) then
!        G functions:
!         G_0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
          dout(5) =  din(1)
          dout(6) = -din(2)
          dout(4) =  din(3)
          dout(7) =  din(4)
          dout(3) =  din(5)
          dout(8) = -din(6)
          dout(2) =  din(7)
          dout(9) =  din(8)
          dout(1) =  din(9)
       else if (type.gt.4) then
          call errquit(' molden interface support up to g functions',
     M         type,0)
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
      logical function nwmolden_openfile(rw,filename)
      implicit none
#include "errquit.fh"
#include "inp.fh"
      character*(*) filename
      character*(*) rw
!
!
      if(rw.eq.'write') then
      write(6,*) 'MOLDEN output will be written to the file ',
     C        filename(1:inp_strlen(filename))
      else
      write(6,*) 'MOLDEN output will be read from file ',
     C        filename(1:inp_strlen(filename))
      endif
      if (MOLDEN_UNIT.gt.9) then
         open(unit=MOLDEN_UNIT,
     1        file=filename,
     2        status='unknown',
     3        form='formatted',
     4        err=911)
      endif
!
      if(rw.eq.'write') then
         write(MOLDEN_UNIT,100) '[Molden Format]'
      else
         read(MOLDEN_UNIT,*) 
      endif
         
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
      logical function nwmolden_geom(rw,geom)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "global.fh"
      character*(*) rw          ! [input]
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
      logical status,do_io
!
! this might be wrong if something other than atomic units
! are specified in the NWChem geometry input, etc.
!
      do_io=ga_nodeid().eq.0
      if(do_io) then
      if(rw.eq.'write') then
      write(MOLDEN_UNIT,100) '[Atoms] AU'
 100  format(a10)
      else
         read(MOLDEN_UNIT,*)
      endif
      endif
!
      status = geom_check_handle(geom, 'nwmolden_geom')
      if (.not.status)
     1   call errquit('nwmolden_geom: error code = ',1,GEOM_ERR)
!
      status = geom_ncent(geom, natoms)
      if (.not.status)
     1   call errquit('nwmolden_geom: error code = ',2,GEOM_ERR)
!
      if(do_io) then
      do iat=1,natoms
!
         status = geom_cent_get(geom, iat, tag, coord, qnuc)
         if (.not.status)
     1      call errquit('nwmolden_geom: error code = ',3,GEOM_ERR)
!
         status = geom_tag_to_element(tag, symbol, element, atn)
         if (.not.status) symbol = 'X'
!
         if(rw.eq.'write') then
           write(MOLDEN_UNIT,300) symbol,iat,nint(qnuc),(coord(j),j=1,3)
         else
            read(MOLDEN_UNIT,*)
         endif
 300     format(1x,a4,1x,i4,1x,i4,1x,3f17.10)
!
      enddo
      endif
!
      nwmolden_geom = .true.
!
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_basis(rw,basis)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "errquit.fh"
#include "global.fh"
      character*(*) rw          ! [input]
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
      logical status,do_io
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
      do_io=ga_nodeid().eq.0
      if(do_io) then
      if(rw.eq.'write') then
      write(MOLDEN_UNIT,400) '[GTO]'
 400  format(a5)
      else
         read(MOLDEN_UNIT,*)
      endif
      endif
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
      if(do_io) then
      do iatom = 1, natoms
        status = bas_ce2cnr(basis,iatom,ifirst,ilast)
        if (.not.status)
     1      call errquit('nwmolden_basis: error code = ',4,BASIS_ERR)
      if(rw.eq.'write') then
        write(MOLDEN_UNIT, 450) iatom, 0
 450    format(1x,i4,1x,i2)
      else
         read(MOLDEN_UNIT,*)
      endif
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
         if(rw.eq.'write') then
         if (type.ge.0) then
           write(MOLDEN_UNIT,500) ctype(type),nprimo,0
         elseif (type.lt.0) then
           write(MOLDEN_UNIT,500) cltype(abs(type)),nprimo,0
        endif
      else
         read(MOLDEN_UNIT,*)
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
            if(rw.eq.'write') then
            write(MOLDEN_UNIT,700) exp(iprimo),coeff(iprimo)
 700        format(1x,f20.10,1x,f20.10)
         else
            read(MOLDEN_UNIT,*)
         endif
         enddo
        enddo
      if(rw.eq.'write') then
         write(MOLDEN_UNIT,*)
      else
         read(MOLDEN_UNIT,*)
      endif
      enddo
      endif
!
      nwmolden_basis = .true.
!
      return
      end
!
! ===========================================================================
!
      logical function nwmolden_movecs(rw,basis,nbf,nmo,
     1                                 eval,occ,sym,g_movecs,scftype)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "inp.fh"
      character*(*) rw          ! [input] read or write
      integer nbf               ! [input] number of AOs
      integer nmo               ! [input] number of MOs
      double precision eval(*)  ! [input] MO energies
      double precision occ(*)   ! [input] MO occupation numbers
      integer sym(*)            ! [input] MO symmetries
      integer g_movecs(*)       ! [input] GA with vectors
      character*(*) scftype     ! [input] RHF, ROHF, UHF, DFT, or MCSCF
      integer numcont           ! number of mapped contractions
      integer icont             ! contraction index
      integer type              ! type (sp/s/p/d/..)
      integer nprimo            ! no. of primitives
      integer ngeno             ! no. of contractions
      integer ibusfunc          ! basis func with diffrent m
      integer sphcart           ! 0/1 for cartesian/spherical
      integer k_vec, l_vec,k_reord,l_reord
      double precision reordered_mo(nbf)
      integer alfabeta, spinref
      integer orb, aof, offset, aof_in
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
      integer orb_alfabeta
      integer*4 iostatus
      integer nmo_out
      integer nwmolden_nbusfunc ! function defined in this file
      character*4 cc4
      character*5 cc5
      character*6 cc6
      character*7 cc7
! external function declarations
      logical sym_char_table
      external sym_char_table
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
      spinref=1
      if (scftype.eq.'RHF')   spinref = 1
      if (scftype.eq.'ROHF')  spinref = 1
      if (scftype.eq.'UHF')   spinref = 2
      if (scftype.eq.'DFT')   spinref = 1 ! ambiguous
      if (scftype.eq.'MCSCF') spinref = 1
      if(ga_nodeid().eq.0) write(6,*) ' molden scftype: ',scftype,
     E     ' spinref: ',spinref
!
!     allocate local memory
!
!     TODO: this is not scalable and probably will not work in DDSCF scheme
!
      if (.not.ma_push_get(mt_dbl,nbf,'vector',l_vec,k_vec))
     1    call errquit('nwmolden_movecs: MA alloc. failed',nbf,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf,'vector',l_reord,k_reord))
     1    call errquit('nwmolden_movecs: MA alloc. failed',nbf,MA_ERR)
!
      status = bas_numcont(basis, numcont)
      if (.not.status)
     1    call errquit('nwmolden_movecs: error code = ',3,BASIS_ERR)
!
      status = bas_continfo(basis,1,type,nprimo,ngeno,sphcart)
!
      if(rw.eq.'write')  then
         if (sphcart.eq.1) then
            write(MOLDEN_UNIT,50) '[5D]'
            write(MOLDEN_UNIT,50) '[9G]'
         end if
         write(MOLDEN_UNIT,50) '[MO]'
      elseif (rw.eq.'read')  then
c     rewind since we migh have a different basis
         rewind(MOLDEN_UNIT)
c     search for the [MO] string
         do
            read(MOLDEN_UNIT,'(a4)') cc4
            if(cc4.eq.'[MO]') exit
         enddo
      endif
!
      do alfabeta = 1, spinref
         cspin = cspins(alfabeta)
         do orb = 1, nmo
            if(rw.eq.'write')  then
!
!          pull the MOs from a GA
!
          call ga_get(g_movecs(alfabeta), 1, nbf, orb, orb, 
     D           dbl_mb(k_vec), 1)
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
     1                              dbl_mb(k_reord+aof-1),
     2                              type, sphcart)
              aof = aof + nbusfunc
           enddo
        endif ! for write
!
!          print stuff
!
           orb_alfabeta=(alfabeta-1)*nbf+orb
           if(rw.eq.'write') then
           csym = zir(sym(orb_alfabeta))
           write(MOLDEN_UNIT,100) csym,eval(orb_alfabeta),
     C             cspin,occ(orb_alfabeta)
        else
           read(MOLDEN_UNIT,1000,iostat=iostatus)
     C          cc5,csym,
     c          cc5,eval(orb_alfabeta),
     C          cc6,cspin,
     c          cc7,occ(orb_alfabeta)
           if(iostatus.ne.0) goto 1234
        endif
           do aof = 1, nbf
              if(rw.eq.'write')  then
                 write(MOLDEN_UNIT,200) aof,dbl_mb(k_reord+aof-1)
              else
                 read(MOLDEN_UNIT,200,iostat=iostatus)
     A                aof_in,dbl_mb(k_reord+aof-1)
                 if(iostatus.ne.0) goto 1234
              endif
           enddo                ! aof
           nmo_out=orb
            if(rw.eq.'read') then
!
!          reorder MO
!
           aof = 1
!
           do icont = 1, numcont
              if(.not. bas_continfo(basis,icont,type,
     1                              nprimo,ngeno,sphcart))
     1            call errquit('nwmolden_movecs: error code = ',
     2                         4,BASIS_ERR)
              call moldennw_reorder(dbl_mb(k_reord + aof - 1),
     1                              dbl_mb(k_vec+aof-1),
     2                              type, sphcart)
              aof = aof + nwmolden_nbusfunc(type, sphcart)
           enddo
!
!          put the MOs into a GA
!
          call ga_put(g_movecs(alfabeta), 1, nbf, orb, orb, 
     D           dbl_mb(k_vec), 1)
       endif
        enddo    ! orb
 1234   continue
        if(rw.eq.'read') then
           if(nmo.ne.nmo_out) then
              write(6,*) ' iostat ',iostatus
              write(6,*) ' NWChem: nmo =',nmo
              write(6,*) ' Molden: nmo =',nmo_out
              call errquit('nwmolden_movecs: nmo molden file',0,0)
           endif
           nmo=nmo_out
        endif
      enddo      ! alfabeta
!
!     de-allocate local memory
!
      if (.not. ma_chop_stack(l_vec))
     1    call errquit('nwmolden_movecs: MA stack is corrupt',5,MA_ERR)
!
      nwmolden_movecs = .true.
!
  50  format(a4)
 100       format('Sym= ', a/
     1            'Ene= ', E22.14/
     2            'Spin= ', a/
     3            'Occup= ', f18.10)
 1000       format(a, a/
     1            a, E22.14/
     2            a, a/
     3            a, f18.10)
 200          format(i6,f18.12)
      return
      end
!
! ===========================================================================
!

