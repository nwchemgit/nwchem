      logical function et_movecs_read(filename, iset, occ
     &     , evals, g_vecs, energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) filename
      integer iset              ! No. (1,2) of set of vectors to read
      double precision occ(*)   ! Must be at least nbf long (not nmo)
      double precision evals(*) ! Must be at least nbf long (not nmo)
      integer g_vecs
c
      integer nsets             ! No. of sets of vectors
      integer nbf               ! No. of functions in basis
      integer nmo(2)            ! No. of vectors in each set
      integer ok, jset, i, j
      integer l_vecs, k_vecs
      integer unitno
      integer mitob1,mdtob_nbf
      parameter (unitno = 67)
      double precision energy, enrep
c     
      l_vecs = -1               ! An invalid MA handle
c
      if (oscfps) call pstat_on(ps_vecio)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c     Skip over uninteresting bits of the header
c
         read(unitno, err=1001, end=1001) ! convergence info
         read(unitno, err=1001, end=1001) ! scftype
         read(unitno, err=1001, end=1001) ! lentit
         read(unitno, err=1001, end=1001) ! title
         read(unitno, err=1001, end=1001) ! lenbas
         read(unitno, err=1001, end=1001) ! basis_name
         read(unitno, err=1001, end=1001) nsets
         read(unitno, err=1001, end=1001) nbf
         read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)
c
         if (.not. ma_push_get(mt_dbl,nbf,'et_movecs_read'
     &       ,l_vecs,k_vecs))
     $        call errquit('et_movecs_read: ma failed', nbf, MA_ERR)
c
c     Skip over unwanted sets
c
         do jset = 1, iset-1
            read(unitno, err=1001, end=1001)
            read(unitno, err=1001, end=1001)
            do i = 1, nmo(jset)
               read(unitno, err=1001, end=1001)
            enddo
         enddo
c
      if (iset.eq.3) then
        read(unitno, err=1001) energy, enrep
        if (energy.eq.0.0d0.and.ga_nodeid().eq.0) then 
          write(6,*) 
          write(6,*) 'H(RR) or H(PP) not found in file ',filename
          write(6,*) 'The vectors are not converged. If you know'
          write(6,*) 'H(RR) and H(PP), you can still calculate V(RP)'
          write(6,*) 'using this equation : '
          write(6,*) 
          write(6,*) '             H(RP) - S(RP)*0.5*( H(RR) + H(PP) ) '
          write(6,*) 'V(RP) =      ----------------------------------- '
          write(6,*) '                   1 - S(RP)*S(RP)           '
          write(6,*) 
          energy = 0.0d0
        else
          energy = energy - enrep
        endif
      else
c
         read(unitno, err=1001, end=1001) (occ(j),j=1,nbf)
         read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
         do i = 1, nmo(iset)
            call sread(unitno, dbl_mb(k_vecs), nbf)
            call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), 1)
         enddo
      endif
c
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
	 if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('et_movecs_read: pop failed', l_vecs, MA_ERR)
      endif
      call ga_sync()
      mitob1= MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob1, 0) ! Propagate status
      if (ok .eq. 1) then
	 call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, mitob1, 0)
         mdtob_nbf= MA_sizeof(MT_DBL,nbf,MT_BYTE)
	 call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, mdtob_nbf, 0)
	 call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ,   mdtob_nbf, 0)
      endif
c
      et_movecs_read = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. et_movecs_read .and.
     $     util_print('vectors i/o', print_default)) then
c         write(6,22) filename(1:inp_strlen(filename))
c 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(6)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' et_movecs_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' et_movecs_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' et_movecs_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
c $Id$
