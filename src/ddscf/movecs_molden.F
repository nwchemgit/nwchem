      subroutine movecs_molden(rtdb, geom, basis, nset, nmo,
     G     g_vecs, evals, foccs)
C$Id$
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
c
c     Read the input MO's from a molden file
c     no checks are made on the basis set
c     Syntax: vectors input molden filename
c
      integer rtdb, geom, basis, g_vecs(2)
      double precision evals(*),foccs(*)
c
      integer nset ! [in] 1 = rhf, 2 = uhf 
      integer nmo ! [in] might be redefined by molden file??
      integer nbf
      character*80 molden_file
      character*3 scftype
      integer l_ir,k_ir
      logical int_normalize
      external int_normalize
      logical oprint, oprint_info
      character*15 subname
      parameter (subname='movecs_molden')
c
      oprint = util_print('mo guess', print_default)
      oprint_info = (ga_nodeid() .eq. 0 .and. oprint)
c      
      if (.not. rtdb_cget(rtdb, 'vectors:molden', 1, molden_file))
     $     call errquit('movecs_molden: no file',0, 0)
c
      if (oprint_info) then
         write(luout,*)
         call util_print_centered(luout,'Guess from Molden file',15,
     $        .true.)
         write(luout,1)
     $        molden_file(1:inp_strlen(molden_file))
 1       format(/' input vectors filename    = ',a)
         call util_flush(luout)
      endif
c
c     extract MO parameters
c
      if (.not. (bas_numbf(basis, nbf)))
     &     call errquit(subname//': basnumbf failed',1, BASIS_ERR)
c     might need to detect nmo  from reading Molden file
c      nmo=nbf
      
c
      if (.not.MA_Push_Get(mt_int,nset*nbf,'dft:irreps',l_ir,k_ir))
     &   call errquit('dft_scf: cannot allocate irreps',0, MA_ERR)
        call ga_sync()
         if (nset .eq. 1) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
      call dft_nwmolden(rtdb,geom,
     I     nset,evals,foccs,int_mb(k_ir),
     N     nbf, nmo,
     I     g_vecs,scftype,'read')
c         if(.not.nwmolden_movecs('read',bas,nbf,nmo,
c     1        eval,occ,ir,
c     G        g_movecs,scftype))
c     E        call errquit(pname//'nwmolden_movecs failed ',0,0)
         
c
      if (.not.ma_chop_stack(l_ir))
     &   call errquit(subname//': cannot chop stack',10, MA_ERR)
      call ga_sync()
c
c
      return 
      end
c
c

