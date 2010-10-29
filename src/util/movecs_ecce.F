      subroutine movecs_ecce_print_on()
      implicit none
#include "ceccemoprint.fh"
      flag = .true.
      end
      subroutine movecs_ecce_print_off()
      implicit none
#include "ceccemoprint.fh"
      flag = .false.
      end
      subroutine movecs_ecce(nbf, nmo, lo, hi, eval, occ, sym,
     $     g_movecs, scftype, key)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "ceccemoprint.fh"
      integer nbf               ! [input] no. of basis functions
      integer nmo               ! [input] no. of MOs
      integer lo, hi            ! [input] range of MOs to print
      double precision eval(*)  ! [input] MO energies
      double precision occ(*)   ! [input] MO occupation numbers
      integer sym(*)            ! [input] MO symmetries
      integer g_movecs          ! [input] GA with vectors
      character*(*) scftype     ! [input] Type of SCF
      character*(*) key         ! [input] Part of keyword for ECCE
c
c     Output all of the evals, occupation numbers and symmetries and a
c     subset of the molecular orbitals to the ECCE print file.
c
c     The ECCE keywords are constructed as 
c
c     <scftype> <key> <molecular orbitals>
c     <scftype> <key> <molecular orbital energies>
c     ...
c
c     scftype is presently either RHF, ROHF, UHF, DFT, or MCSCF.
c
c     Key should be blank for restricted spin methods or alpha
c     or beta for unrestricted methods as appropriate.
c
c     This routine is temporarily being used to print ALL of the
c     MOs to the ECCE file.  This is not suitable for large 
c     calculations which should just print frontier orbitals
c
      integer l_v, k_v, stat1, stat2, range(2), slen, klen
      character*80 ceval, cocc, cvec, cran, csym
c
c     If ECCE is not printing then go home.  Note that only process
c     zero will really know but only it is printing.
c
      stat1 = 0
      call ecce_print_control(stat1, stat2) ! stat2 set to print status
      stat1 = stat2
      call ecce_print_control(stat1, stat2) ! restore previous setting
      if (stat1 .ne. 1) return
c
      if (.not. flag) return    ! Control for geometry optimization
c
c     Get the MOs
c
      if (.not. ma_push_get(mt_dbl, (hi-lo+1)*nbf, 'moecce', l_v, k_v))
     $     call errquit('movecs_ecce: insufficient memory',0, MA_ERR)
      call ga_get(g_movecs, 1, nbf, lo, hi, dbl_mb(k_v), nbf)
c
c     Construct the names
c
      ceval = ' '
      cocc  = ' '
      cvec  = ' '
      cran  = ' '
      csym  = ' '
      slen = inp_strlen(scftype)
      klen = inp_strlen(key)
      if (key .eq. ' ') then
         write(ceval,'(a,1x,a)')
     $      'molecular orbital energies',
     &      scftype(1:slen)
         write(cocc,'(a,1x,a)')
     $       'molecular orbital occupations',
     &       scftype(1:slen) 
         write(cvec,'(a,1x,a)')
     $        'molecular orbital vectors',
     &       scftype(1:slen)
         write(cran,'(a,1x,a)')
     $       'molecular orbital range',
     &       scftype(1:slen)
         write(csym,'(a,1x,a)')
     $       'molecular orbital symmetries',
     &       scftype(1:slen)
      else
         write(ceval,'(a,1x,a,1x,a)')
     &      'molecular orbital energies',
     &      scftype(1:slen), 
     $      key(1:klen)
         write(cocc,'(a,1x,a,1x,a)')
     &       'molecular orbital occupations',
     &       scftype(1:slen), 
     $       key(1:klen)
         write(cvec,'(a,1x,a,1x,a)')
     &       'molecular orbital vectors',
     &       scftype(1:slen),
     $       key(1:klen)
         write(cran,'(a,1x,a,1x,a)')
     &       'molecular orbital range',
     &       scftype(1:slen),
     $       key(1:klen)
         write(csym,'(a,1x,a,1x,a)')
     &       'molecular orbital symmetries',
     &       scftype(1:slen),
     $       key(1:klen)
      endif
      range(1) = lo
      range(2) = hi
c
      call ecce_print1(ceval, mt_dbl, eval, nmo)
      call ecce_print1(cocc,  mt_dbl, occ,  nmo)
      call ecce_print1(csym,  mt_int, sym,  nmo)
      call ecce_print1(cran,  mt_int, range, 2)
      call ecce_print2_dbl_tol(cvec,  dbl_mb(k_v), nbf, nbf, (hi-lo)+1,
     $     1d-3)
c
      if (.not. ma_pop_stack(l_v)) call errquit
     $     ('movecs_ecce: ma stack is corrupt', 0, MA_ERR)
c
      end
