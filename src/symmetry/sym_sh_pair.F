      logical function sym_atom(geom, iat, q1)
      implicit none
      integer geom              ! Geometry handle [input]
      integer iat               ! Atom index [input]
      double precision q1       ! Constituency number [output]
c
c     Return true if (iatom) is the lexically highest atom
c     symmetry equivalent atoms. If true also return the
c     constituency factor q1 (= no. of symmetry equivalent atoms)
c
c     As an optimization inline sym_center_map ... the only reason
c     geomP.fh and mafdecls.fh are included here ... the interface
c     is compatible with that of the actual routine.
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer iat_new
      integer nops, op
      integer n                 ! Counts no. of equivalent pairs
      integer sym_number_ops, sym_center_map_inline
      external sym_number_ops
#include "itri.fh"
      sym_center_map_inline(geom,iat,op) = int_mb(op - 1 + 
     $     sym_center_map_index(geom) + (iat-1)*sym_num_ops(geom))
c
      q1 = 0.0d0
      sym_atom = .false.
c
c     Loop thru operations in the group and map to new pairs
c
      nops = sym_number_ops(geom)
      n = 1                     ! Identity always counts
      do op = 1, nops
c
c     Map centers
c
         iat_new = sym_center_map_inline(geom, iat, op)
c
c     Compare index
c
         if (iat .lt. iat_new) then
            return
         else if (iat .eq. iat_new) then
            n = n + 1
         endif
      end do
c
      q1 = dble(nops + 1) / dble(n)
      sym_atom = .true.
c
      end
      logical function sym_shell(basis, ishell, q1)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"      
#include "global.fh"
      integer basis             ! Basis set handle [input]
      integer ishell            ! Shell index [input]
      double precision q1       ! Constituency number [output]
c
c     Return true if (ishell) is the lexically highest shell
c     of symmetry equivalent shells. If true, also return the
c     constituency factor q1 (= no. of symmetry equivalent shells).
c
      integer ice
      integer geom
      logical sym_atom
      external sym_atom
c
      sym_shell = .false.
c
c     Get geometry handle, centers where shells are located and
c     number of group operations
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_shell: bas_geom?', 0, BASIS_ERR)
      if (.not. bas_cn2ce(basis, ishell, ice)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', ishell, BASIS_ERR)
c
      sym_shell = sym_atom(geom, ice, q1)
c
      end
      logical function sym_shell_pair(basis, ishell, jshell, q2)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"      
#include "global.fh"
      integer basis             ! Basis set handle [input]
      integer ishell, jshell    ! Shell indices [input]
      double precision q2       ! Constituency number [output]
c
c     Return true if (ishell,jshell) is the lexically highest
c     pair of symmetry equivalent shells. If true, also return the
c     constituency factor q2 (= no. of symmetry equivalent pairs).
c
c     This routine uses the exchange symmetry ishell <-> jshell
c     and incorporates a factor of two into q2 to account for
c     this.
c
      integer ice, jce
      integer geom
      logical sym_atom_pair
      external sym_atom_pair
c
      sym_shell_pair = .false.
      if (ishell.lt.jshell) return
c
c     Get geometry handle, centers where shells are located and
c     number of group operations
c

*      write(6,*) ' sym_shell_pair ', ga_nodeid(), ishell, jshell
*      call util_flush(6)
*      call ga_sync()

      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_shell_pair: bas_geom?', 0, BASIS_ERR)
      if (.not. bas_cn2ce(basis, ishell, ice)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', ishell, BASIS_ERR)
      if (.not. bas_cn2ce(basis, jshell, jce)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', jshell, BASIS_ERR)
c
*      write(6,*) ' sym_shell_pair centers', ga_nodeid(), ice, jce
*      call util_flush(6)
*      call ga_sync()
c
      sym_shell_pair = sym_atom_pair(geom, ice, jce, q2)
      if (ishell .ne. jshell) q2 = q2 + q2
c
*      write(6,*) ' sym_shell_pair return ', ga_nodeid(),
*     $     sym_shell_pair, q2
*      call util_flush(6)
*      call ga_sync()
c
      end
      logical function sym_atom_pair(geom, iat, jat, q2)
      implicit none
      integer geom              ! Geometry handle [input]
      integer iat, jat          ! Atom indices [input]
      double precision q2       ! Constituency number [output]
c
c     Return true if (iatom,jatom) is the lexically highest
c     pair of symmetry equivalent atom. If true also return the
c     constituency factor q2 (= no. of symmetry equivalent pairs)
c
c     This routine uses the exchange symmetry iat <-> jat
c     but does not incorporate any factors into q2 to account for
c     this (q2 is point group symmetry only). 
c
c     As an optimization inline sym_center_map ... the only reason
c     geomP.fh and mafdecls.fh are included here ... the interface
c     is compatible with that of the actual routine.
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer iat_new, jat_new, ijat, ijat_new
      integer nops, op
      integer n                 ! Counts no. of equivalent pairs
      integer sym_number_ops, sym_center_map_inline
      external sym_number_ops
#include "itri.fh"
      sym_center_map_inline(geom,iat,op) = int_mb(op - 1 + 
     $     sym_center_map_index(geom) + (iat-1)*sym_num_ops(geom))
c
      q2 = 0.0d0
      sym_atom_pair = .false.
      if (iat .lt. jat) return
      ijat = itri(iat,jat)
c
c     Loop thru operations in the group and map to new pairs
c
      nops = sym_number_ops(geom)
      n = 1                     ! Identity always counts
      do op = 1, nops
c
c     Map centers
c
         iat_new = sym_center_map_inline(geom, iat, op)
         jat_new = sym_center_map_inline(geom, jat, op)
c
c     Compare canonical indices
c
         ijat_new = itri(iat_new, jat_new)
         if (ijat .lt. ijat_new) then
            return
         else if (ijat .eq. ijat_new) then
            n = n + 1
         endif
      end do
c
      q2 = dble(nops + 1) / dble(n)
      sym_atom_pair = .true.
c
      end
      logical function sym_atom_quartet(geom, iat, jat, kat, lat, q4)
      implicit none
      integer geom              ! Geometry handle [input]
      integer iat, jat, kat, lat ! Atom indices [input]
      double precision q4       ! Constituency number [output]
c
c     Return true if (iatom,jatom,katom,latom) is the lexically highest
c     quartet of symmetry equivalent atoms. If true also return the
c     constituency factor q4 (= no. of symmetry equivalent quartets)
c
c     This routine uses the standard three exchange symmetries
c
c     (iat<->jat) <-> (kat<->lat)
c
c     ... it is easy to extend it to not use these
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer iat_new, jat_new, kat_new, lat_new
      integer nops, op, ij, kl, ijkl, ijkl_new, ij_new, kl_new
      integer n                 ! Counts no. of equivalent pairs
      integer sym_number_ops, sym_center_map_inline
      external sym_number_ops
#include "itri.fh"
      sym_center_map_inline(geom,iat,op) = int_mb(op - 1 + 
     $     sym_center_map_index(geom) + (iat-1)*sym_num_ops(geom))
c
      sym_atom_quartet = .false.
c
c     Assume that the code is looping thru a unique set of 
c     indices but maybe not with conventional canonical order
c     ... so don't do the canonical check here
*      if (iat .lt. jat) return  ! Labels must be in canonical order
*      if (iat .lt. kat) return
*      if (kat .lt. lat) return
*      if (iat .eq. kat .and. jat .lt. lat) return
c
      q4 = 0.0d0
      ij   = itri(iat,jat)
      kl   = itri(kat,lat)
      ijkl = itri(ij, kl)
c
c     Loop thru operations in the group and map to new pairs
c
      nops = sym_number_ops(geom)
      n = 1                     ! Identity always counts
      do op = 1, nops
c
c     Map centers
c
         iat_new = sym_center_map_inline(geom, iat, op)
         jat_new = sym_center_map_inline(geom, jat, op)
         kat_new = sym_center_map_inline(geom, kat, op)
         lat_new = sym_center_map_inline(geom, lat, op)
c
c     Compare canonical indices
c
         ij_new   = itri(iat_new,jat_new)
         kl_new   = itri(kat_new,lat_new)
         ijkl_new = itri(ij_new, kl_new)
c
         if (ijkl .lt. ijkl_new) then
            return
         else if (ijkl .eq. ijkl_new) then
            n = n + 1
         endif
      end do
c
      q4 = dble(nops + 1) / dble(n)
      sym_atom_quartet = .true.
c
      end
      logical function sym_shell_quartet(basis,
     $     ishell, jshell, kshell, lshell, q4)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"      
#include "global.fh"
      integer basis             ! Basis set handle [input]
      integer ishell, jshell    ! Shell indices [input]
      integer kshell, lshell    ! Shell indices [input]
      double precision q4       ! Constituency number [output]
c
      integer ice, jce, kce, lce
      integer geom
      logical sym_atom_gen_quartet
      external sym_atom_gen_quartet
c
      q4 = 0.0d0
      sym_shell_quartet = .false.
      if (ishell.lt.jshell) return
      if (kshell.lt.lshell) return
c
c     Get geometry handle, centers where shells are located and
c     number of group operations
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_shell_pair: bas_geom?', 0, BASIS_ERR)
      if (.not. bas_cn2ce(basis, ishell, ice)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', ishell, BASIS_ERR)
      if (.not. bas_cn2ce(basis, jshell, jce)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', jshell, BASIS_ERR)
      if (.not. bas_cn2ce(basis, kshell, kce)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', kshell, BASIS_ERR)
      if (.not. bas_cn2ce(basis, lshell, lce)) call errquit
     $     ('sym_shell_pair: bas_cn2ce', lshell, BASIS_ERR)
c
      sym_shell_quartet = sym_atom_gen_quartet(
     $     geom, ice, jce, kce, lce, q4)
c
      end
      logical function sym_atom_gen_quartet(geom,
     $     iat, jat, kat, lat, q4)
      implicit none
#include "errquit.fh"
      integer geom              ! Geometry handle [input]
      integer iat, jat, kat, lat ! Atom indices [input]
      double precision q4       ! Constituency number [output]
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer iat_new, jat_new, kat_new, lat_new
      integer nops, op, ij, kl, ij_new, kl_new
      integer n                 ! Counts no. of equivalent pairs
      integer sym_number_ops, sym_center_map_inline
      external sym_number_ops
#include "itri.fh"
      sym_center_map_inline(geom,iat,op) = int_mb(op - 1 + 
     $     sym_center_map_index(geom) + (iat-1)*sym_num_ops(geom))
c
      q4 = 0.0d0
      sym_atom_gen_quartet = .false.
c
      if (iat .lt. jat) return  ! Labels must be in canonical order
      if (kat .lt. lat) return
c
      ij   = itri(iat,jat)
      kl   = itri(kat,lat)
c
c     Loop thru operations in the group and map to new pairs
c
      nops = sym_number_ops(geom)
      n = 1                     ! Identity always counts
      do op = 1, nops
c
c     Map centers
c
         iat_new = sym_center_map_inline(geom, iat, op)
         jat_new = sym_center_map_inline(geom, jat, op)
         kat_new = sym_center_map_inline(geom, kat, op)
         lat_new = sym_center_map_inline(geom, lat, op)
c
c     Compare canonical indices
c
         ij_new   = itri(iat_new,jat_new)
         kl_new   = itri(kat_new,lat_new)
c
         if (ij .lt. ij_new) return
         if (ij.eq.ij_new .and. kl.lt.kl_new) return
         if (ij.eq.ij_new .and. kl.eq.kl_new) then
            n = n + 1
         endif
      end do
c
      q4 = dble(nops+1) / dble(n)
      if (abs(q4-nint(q4)).gt.1e-12) call errquit
     $     ('sym_atom_gen_quartet: not divisible', 0, BASIS_ERR)
      sym_atom_gen_quartet = .true.
c
      end
