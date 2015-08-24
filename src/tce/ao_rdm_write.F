      logical function ao_rdm_write(filename,k_rdm)
crtdb, basis, filename, 
c     $     scftype, title,
c     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "tce_main.fh"
c
c     Temporary routine
c
c      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
c      integer basis             ! [input] Basis handle(-1 if not accessible)
c      character*(*) filename    ! [input] File to write to
c      character*(*) scftype     ! [input] Type of SCF (dft, scf, mcscf)
c      character*(*) title       ! [input] Title of job that created vectors
c      integer nbf               ! [input] No. of functions in basis
c      integer nsets             ! [input] No. of sets of vectors
c      integer nmo(nsets)        ! [input] No. of vectors in each set
c      integer ldocc             ! [input]
c      integer ldevals           ! [input]
c      double precision occ(ldocc, nsets) ! [input] Occupation numbers 
c      double precision evals(ldevals, nsets) ! [input] Eigenvalues
c      integer g_vecs(nsets)     ! Global array with eigen-vectors
c
      character*(*) filename
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
c      integer lentit
c      integer lenbas
c      integer l_vecs, k_vecs
      integer ok
c, iset, i, j
c      integer geom, ma_type, nelem
c      character*26 date
c      character*32 geomsum, basissum, key
c      character*20 scftype20    ! DECL MUST match movecs_converged routine
c      character*128 basis_name, trans_name
c      double precision energy, enrep
      integer inntsize
      integer k_rdm, i,ind,j
c
c      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
c      if (basis.ne.-1 .and. rtdb.ne.-1) then
c         if (.not. bas_checksum(basis, basissum))
c     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
c         if (.not. bas_geom(basis, geom))
c     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
c         if (.not. geom_checksum(geom, geomsum))
c     $        call errquit('movecs_write: bad geometry',0, GEOM_ERR)
c         if (.not. bas_name(basis, basis_name, trans_name))
c     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
c         if (trans_name .ne. ' ') basis_name = trans_name
c         
c         key = ' '
c         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
c         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
c     $        date = ' '
c      else
c         basissum = ' '
c         geomsum  = ' '
c         date     = ' '
c      endif
c      scftype20 = scftype
c
c
      ok = 0
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it's given a value of zero.
c
c      write(key,'(a,'':energy'')') scftype(1:inp_strlen(scftype))
c      if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy)) then
c        energy=0.0d0
c      endif
c
c get nuclear repulsion energy
c      if (.not. geom_nuc_rep_energy(geom, enrep))
c     $  call errquit('unable to get nuclear rep energy',geom, GEOM_ERR)
c
      if (ga_nodeid() .eq. 0) then
c         if (nsets .gt. 2) then
c            write(6,*) ' movecs_write: nsets > 2 ', nsets
c            goto 10
c         endif
         open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
         do i = 1,nbf
            ind = (i-1)*nbf
            call swrite(unitno, dbl_mb(k_rdm+ind), nbf)
         enddo
c
c     Information about convergence
c
c         write(unitno, err=1001) basissum, geomsum, scftype20, date
c
c     Check that read routines are both consistent with this
c
c         write(unitno, err=1001) scftype20
c         lentit = max(1,inp_strlen(title)) ! 0 length record confuses f2c
c         write(unitno, err=1001) lentit
c         write(unitno, err=1001) title(1:lentit)
c         lenbas = max(1,inp_strlen(basis_name))
c         write(unitno, err=1001) lenbas
c         write(unitno, err=1001) basis_name(1:lenbas)
c         write(unitno, err=1001) nsets
c         write(unitno, err=1001) nbf
c         write(unitno, err=1001) (nmo(i),i=1,nsets)
c         if (.not. ma_push_get(mt_dbl,nbf,'movecs_write',l_vecs,k_vecs))
c     $        call errquit('movecs_write: ma failed', nbf, MA_ERR)
c         do iset = 1, nsets
c            write(unitno, err=1001) (occ(j,iset),j=1,nbf)
c            write(unitno, err=1001) (evals(j,iset),j=1,nbf)
c            do i = 1, nmo(iset)
c               call ga_get(g_vecs(iset), 1, nbf, i, i, dbl_mb(k_vecs),1)
c               call swrite(unitno, dbl_mb(k_vecs), nbf)
c*               write(unitno, err=1001) (dbl_mb(k_vecs+j), j=0,nbf-1)
c            enddo
c         enddo
c         if (.not. ma_pop_stack(l_vecs))
c     $        call errquit('movecs_write: ma pop failed', l_vecs,
c     &       MA_ERR)
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
c         write(unitno, err=1001) energy, enrep
c         
         close(unitno,err=1002)
         ok = 1
      endif
c
   10 call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
c
      ao_rdm_write = ok .eq. 1
c      if (ga_nodeid() .eq. 0 .and.
c     $     util_print('vectors i/o', print_high)) then
c         write(6,22) filename(1:inp_strlen(filename))
c  22     format(/' Wrote molecular orbitals to ',a/)
c         call util_flush(luout)
c      endif
      call ga_sync()
c      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' ao_rdm_write: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' ao_rdm_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
c $Id$
