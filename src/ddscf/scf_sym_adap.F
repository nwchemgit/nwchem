      subroutine scf_movecs_sym_adapt(basis, g_movecs, oprint,
     $     nprint, name, oirs, irs)
      implicit none
#include "bas.fh"
#include "inp.fh"
#include "global.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
      integer basis
      integer g_movecs
      logical oprint
      integer nprint
      character*(*) name
      logical oirs
      integer irs(*)
c
c     if (oirs) return the irreps of each mo in irs()
c
      integer l_irs, k_irs
      integer type, nbf, nvec, i, ilo, ihi, nmixed
      integer geom
c
      integer maxireps
      parameter (maxireps=20)
      integer nir, nop, iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      character*80 buf
c
      if (oprint .and. ga_nodeid() .eq. 0) then
         write(6,*)
         buf = ' '
         write(buf,'(a,a)') 
     $        'Symmetry analysis of molecular orbitals ',
     $        name(1:inp_strlen(name))
         call util_print_centered(6, buf, 29, .true.)
         write(6,*)
         call util_flush(6)
      endif
c
      call ga_inquire(g_movecs, type, nbf, nvec)
      if (.not. ma_push_get(mt_int, nvec, 'irs', l_irs, k_irs))
     $     call errquit('scf_movecs_sym_adapt: memory?', nvec)
c
c     What should the threshold be for reporting contamination?  
c     Probably should be settable from the input or determined
c     from the convergence of the SCF.
c
      if (oscfps) call pstat_on(ps_sym_ad)
      call sym_movecs_adapt(basis, 1d-8, g_movecs, int_mb(k_irs),
     $     nmixed)
      if (oscfps) call pstat_off(ps_sym_ad)
      if (nmixed.ne.0 .and. ga_nodeid().eq.0) then
         write(6,12) nmixed
 12      format(/' !! scf_movecs_sym_adapt:', i5, ' vectors were',
     $        ' symmetry contaminated'/)
         call util_flush(6)
      endif
c
c     Make any near zeros exact zeroes
c
      call ga_screen(g_movecs, 1d-13)
c
c     Must reorthog if mixing was present
c
      if (nmixed .gt. 0) call ga_orthog_mos(basis, g_movecs)
c
      if (oprint .and. ga_nodeid() .eq. 0) then
         if (.not. bas_geom(basis, geom)) 
     $        call errquit('scf_p_v_s:geom',0)
c
         call sym_group_name(geom, group_name)
         if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $        zir, zclass, chars)) 
     $        call errquit('scf_p_v_s: no char table',0)
c
         write(6,*) ' Numbering of irreducible representations: '
         write(6,*)
         do ilo = 1, nir, 5
            ihi = min(nir, ilo+5-1)
            write(6,2) (i,zir(i),i=ilo,ihi)
         enddo
         write(6,*)
         write(6,*) ' Orbital symmetries:'
         write(6,*)
         do ilo = 1, nprint, 5
            ihi = min(nprint, ilo+5-1)
            write(6,2) (i,zir(int_mb(k_irs+i-1)),i=ilo,ihi)
 2          format(1x,5(i5,1x,a8))
         enddo
         write(6,*)
         call util_flush(6)
      endif
c
      if (oirs) then
	 do i = 1, nvec
	    irs(i) = int_mb(k_irs+i-1)
         enddo
      endif
c
      if (.not. ma_pop_stack(l_irs)) call errquit
     $     ('scf_m_s_a: memory corrupt?', 0)
c
      end
