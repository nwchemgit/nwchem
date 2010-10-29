C
C $Id$
C
      subroutine get_wdens(g_dens,ndens,theory,scftype,rtdb,
     *                     g_wdens,nbfs)
c
c  Form the weighted density out of the density and the eigenvalues
c  Allow more than one density to enter this routine because it could
c  be expanded for use in more than one theory type (ex. UHF, MP2)
c  Eigenvalues need to be read from disk.
c  Start with RHF and work our way up.
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer ndens          ! [input]  number of density matrices
      integer g_dens(ndens)  ! [input]  density handles
      integer g_wdens        ! [output] weighted density matrix
      character*3 theory     ! [input]  HF or DFT
      character*4 scftype    ! [input]  RHF, UHF, or ROHF
      integer rtdb           ! [input]  RTDB handle 
      integer nbfs           ! [input]  number of basis functions
c
      integer ilo, ihi, jlo, jhi, index, nrows, ncols, ld
      integer l_evals, k_evals
      character*255 movecs
c
      logical my_diag_patch
      logical moeig_read
      external my_diag_patch
      external moeig_read
c
c Temporarily disable ROHF
c
      if (scftype.eq.'ROHF')
     *  call errquit('get_wdens: ROHF is not supported yet',555,
     &       CAPMIS_ERR)
c
c Get the movecs file name
c
      if (theory.eq.'hf') then
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs))
     *    call errquit('get_wdens: SCF MO vectors not defined',555,
     &       RTDB_ERR)
      else if (theory.eq.'dft') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs))
     *   call errquit('get_wdens: DFT MO vectors not defined',555,
     &       RTDB_ERR)
      else
        call errquit('get_wdens:theory is unrecognized',555,
     &       RTDB_ERR)
      endif
c
c First copy the density (g_dens) to the weighted density (g_wdens)
c
      call ga_copy(g_dens(1), g_wdens)
c
c Get temporary space for the eigenvalues
c
      if (.not. ma_push_get(mt_dbl, nbfs,'MO evals', l_evals, k_evals))
     *  call errquit('get_wdens: could not allocate l_evals',nbfs,
     *  MA_ERR)
c
c Get the eigenvalues
c
      if (.not. moeig_read(movecs, 1, dbl_mb(k_evals)))
     *  call errquit('get_wdens:could not read eigenvectors',555,
     *  DISK_ERR)
c
c Need to form the weighted density matrix which is basically a
c scaling of the diagonal of the density matrix by the eigenvalues.
c
      if (my_diag_patch(g_wdens,ilo,ihi,jlo,jhi)) then
        call ga_access(g_wdens,ilo,ihi,jlo,jhi,index,ld)
        nrows = ihi-ilo+1
        ncols = jhi-jlo+1
        call scale_diag(dbl_mb(index),ld,nrows,dbl_mb(k_evals+ilo-1))
      endif
c
c Free up temporary space for the eigenvalues
c
      if (.not.ma_pop_stack(l_evals))
     *  call errquit('get_wdens:ma free eval',555, MA_ERR)
      call ga_sync()
      end
c
c
      subroutine scale_diag(a, ld, rowcols, vec)
c
c Scales the diagonal of a matrix by the elements of a vector
c
      integer ld, rowcols
      double precision a(ld, rowcols), vec(rowcols)
      integer i
c
      do i = 1, rowcols
        a(i,i) = a(i,i)*vec(i)
      enddo
      end
c
c
      logical function my_diag_patch(g_a,ilo,ihi,jlo,jhi)
c
c Determines if the node has diagonal elements of g_a locally
c
      implicit none
#include "global.fh"
      integer g_a,ilo,ihi,jlo,jhi
      integer ilog, jlog, ihig, jhig
      logical patch_intersect
      external patch_intersect
c
      call ga_distribution(g_a, ga_nodeid(), ilog, ihig, jlog,jhig)
c
c     [ilo:ihi,jlo:jhi] in the first shot is set to be transpose of local patch
c
      ilo = jlog
      jlo = ilog
      ihi = ihig
      jhi = jhig
c
c     now,  compute [ilo:ihi,jlo:jhi] as principal submatrix of g_a
c
      my_diag_patch =  patch_intersect(ilog, ihig, jlog,jhig,
     *                                 ilo,ihi,jlo,jhi)
      end
c
      logical function patch_intersect
     *        (ilo, ihi, jlo, jhi, ilop, ihip, jlop, jhip)
c
      implicit none
c
      integer ilo, ihi, jlo, jhi, ilop, ihip, jlop, jhip
c
c     check if two patches [ilo:ihi,jlo:jhi] and [ilop:ihip,jlop:jhip] intersect
c     if yes, then modify [ilop:ihip,jlop:jhip] to be intersection
c
      patch_intersect = .false.
c
c     check consistency of patch coordinates
      if( ihi .lt. ilo .or. jhi .lt. jlo)     return ! inconsistent
      if( ihip .lt. ilop .or. jhip .lt. jlop) return ! inconsistent
c
c     find the intersection and update (ilop: ihip, jlop: jhip)
      if( ihi .lt. ilop .or. ihip .lt. ilo) return ! do not intersect
      if( jhi .lt. jlop .or. jhip .lt. jlo) return ! do not intersect
      ilop = MAX(ilo,ilop)
      ihip = MIN(ihi,ihip)
      jlop = MAX(jlo,jlop)
      jhip = MIN(jhi,jhip)

      patch_intersect = .true.
      end

