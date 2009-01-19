c ===========================================================================
c
c     Goal: create a complete NWChem interface to Molden.
c
c     1. need subroutine which prints MOs in batches (DONE)
c     2. need to dump basis set information (DONE)
c     3. need to dump geometry information (DONE)
c     4. (optional) dump scf and geometry optimization convergence information
c     5. (optional) dump normal modes from frequency calculation
c
c     This code was written on January 18th and 19th, 2009 by Jeff Hammond.
c
c ===========================================================================
c
#define MOLDEN_UNIT 37
c
c ===========================================================================
c
      logical function nwmolden_openfile()
      implicit none
#include "errquit.fh"
      character*64 filename
c
      call util_file_prefix ('molden', filename )
c
      if (MOLDEN_UNIT.gt.9) then
         open(unit=MOLDEN_UNIT,
     1        file=filename,
     2        status='unknown',
     3        form='formatted',
     4        err=911)
      endif
c
      write(MOLDEN_UNIT,100) '[Molden Format]'
 100  format(a15)
c
      nwmolden_openfile = .true.
      return
c
 911  nwmolden_openfile = .false.
      return
      end
c
c ===========================================================================
c
      logical function nwmolden_closefile()
      implicit none
c
      if (MOLDEN_UNIT.gt.9) then
        close(unit=MOLDEN_UNIT,
     1        status='keep',
     2        err=911)
      endif
c
      nwmolden_closefile = .true.
      return
c
 911  nwmolden_closefile = .false.
      return
      end
c
c ===========================================================================
c
      logical function nwmolden_geom(geom)
      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer geom
      integer j
      integer iat
      integer atn
      integer natoms
      double precision coord(3)
      double precision qnuc
      character*2 symbol
      character*16 tag
      character*32 element
      character*80 buf
      logical status
c
      write(MOLDEN_UNIT,100) '[Atoms] AU'
 100  format(a10)
c
      status = geom_check_handle(geom, 'nwmolden_geom')
      if (.not.status) 
     1   call errquit('nwmolden_geom: error code = ',1,GEOM_ERR)
c
      status = geom_ncent(geom, natoms)
      if (.not.status) 
     1   call errquit('nwmolden_geom: error code = ',2,GEOM_ERR)
c
      do iat=1,natoms
c
         status = geom_cent_get(geom, iat, tag, coord, qnuc)
         if (.not.status) 
     1      call errquit('nwmolden_geom: error code = ',3,GEOM_ERR)
c
         status = geom_tag_to_element(tag, symbol, element, atn)
         if (.not.status) symbol = 'X'
c
         write(MOLDEN_UNIT,300) symbol,iat,nint(qnuc),(coord(j),j=1,3)
 300     format(1x,a4,1x,i4,1x,i4,1x,3f15.8)
c
      enddo
c
      nwmolden_geom = .true.
c
      return
      end
c
c ===========================================================================
c
      logical function nwmolden_basis(basis)
      implicit none
#include "bas.fh"
#include "errquit.fh"
      integer basis             ! [input] basis handle
      integer numcont           ! number of mapped contractions
      integer i,j,k
      integer icont             ! contraction index
      integer type              ! type (sp/s/p/d/..)
      integer nprimo            ! no. of primitives
      integer ngeno             ! no. of contractions
      integer sphcart           ! 0/1 for cartesian/shperical
      logical status
c
      integer max_prim
      parameter (max_prim = 100)
      double precision exp(max_prim)
      double precision coeff(max_prim)
      character*3 ctype(0:6),cltype(2)
c
      ctype(0)='s'
      ctype(1)='p'
      ctype(2)='d'
      ctype(3)='f'
      ctype(4)='g'
      ctype(5)='h'
      ctype(6)='i'
      cltype(1)='sp'
      cltype(2)='spd'
c
      write(MOLDEN_UNIT,400) '[GTO]'
 400  format(a5)
c
      status = bas_check_handle(basis,'nwmolden_basis')
      if (.not.status) 
     1   call errquit('nwmolden_basis: error code = ',1,BASIS_ERR)
c
      status = bas_numcont(basis,numcont)
      if (.not.status) 
     1   call errquit('nwmolden_basis: error code = ',2,BASIS_ERR)
c
      do icont = 1, numcont
c
         status = bas_continfo(basis,icont,type,nprimo,ngeno,sphcart)
         if (.not.status) 
     1      call errquit('nwmolden_basis: error code = ',3,BASIS_ERR)
c
         if (nprimo.gt.max_prim) 
     1      call errquit('nwmolden_basis: nprimo > max_prim',
     2                   max_prim-nprimo,BASIS_ERR)
c
         if (type.ge.0) then
           write(MOLDEN_UNIT,500) ctype(type),nprimo,0
         elseif (type.lt.0) then
           write(MOLDEN_UNIT,500) cltype(abs(type)),nprimo,0
         endif
 500     format(1x,a4,1x,i4,1x,i2)
c
         status = bas_get_exponent(basis,icont,exp)
         if (.not.status) 
     1      call errquit('nwmolden_basis: error code = ',4,BASIS_ERR)

         status = bas_get_coeff(basis,icont,coeff)
         if (.not.status) 
     1      call errquit('nwmolden_basis: error code = ',5,BASIS_ERR)
c
         do j = 1, nprimo
            write(MOLDEN_UNIT,700) exp(j),coeff(j)
 700        format(1x,f14.8,1x,f14.8)
         enddo
c
      enddo
c
      nwmolden_basis = .true.
c
      return
      end
c
c ===========================================================================
c
      logical function nwmolden_movecs(basis,nbf,nmo,
     1                                 eval,occ,sym,g_movecs,scftype)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer nbf               ! [input] number of AOs
      integer nmo               ! [input] number of MOs
      double precision eval(*)  ! [input] MO energies
      double precision occ(*)   ! [input] MO occupation numbers
      integer sym(*)            ! [input] MO symmetries
      integer g_movecs          ! [input] GA with vectors
      character*(*) scftype     ! [input] RHF, ROHF, UHF, DFT, or MCSCF
      integer k_vec, l_vec
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
      logical sym_char_table, bas_geom
      external sym_char_table, bas_geom
c
      status = bas_geom(basis, geom)
      if (.not.status) 
     1    call errquit('nwmolden_movecs: error code = ',1,BASIS_ERR)
c
      call sym_group_name(geom, zname)
c
      status = sym_char_table(zname,nop,nir,class_dim,zir,zclass,chars)
      if (.not.status) 
     1    call errquit('nwmolden_movecs: error code = ',2,BASIS_ERR)
c
      cspins(1) = 'Alpha'
      cspins(2) = 'Beta '
c
      if (scftype.eq.'RHF')   spinref = 1
      if (scftype.eq.'ROHF')  spinref = 2
      if (scftype.eq.'UHF')   spinref = 2
      if (scftype.eq.'DFT')   spinref = 1 ! ambiguous
      if (scftype.eq.'MCSCF') spinref = 1
c
      write(MOLDEN_UNIT,50) '[MO]'
  50  format(a4)
c
c     allocate local memory
c
      if (.not.ma_push_get(mt_dbl,nbf,'vector',l_vec,k_vec))
     1    call errquit('nwmolden_movecs: MA alloc. failed',nbf,MA_ERR)
c
      do alfabeta = 1, spinref
         cspin = cspins(spinref)
         do orb = 1, nmo
c
c          pull the MOs from a GA
c
           call ga_get(g_movecs, 1, nbf, orb, orb, dbl_mb(k_vec), 1)
c
c          print stuff
c
           csym = zir(sym(orb))
           write(MOLDEN_UNIT,100) csym,eval(orb),cspin,occ(orb)
 100       format('Sym= ', a12/
     1            'Ene= ', f12.5/
     2            'Spin= ', a11/
     3            'Occup= ', f10.5)
           do aof = 1, nbf
              write(MOLDEN_UNIT,200) aof,dbl_mb(k_vec + aof - 1) 
 200          format(i6,f14.8)
           enddo ! aof
        enddo    ! orb
      enddo      ! alfabeta
c
c     de-allocate local memory
c
      if (.not. ma_pop_stack(l_vec)) 
     1    call errquit('nwmolden_movecs: MA stack is corrupt',0,MA_ERR)
c
      nwmolden_movecs = .true.
c
      return
      end
c
c ===========================================================================
c



