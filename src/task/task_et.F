      logical function task_et(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "cscf.fh" 
      integer rtdb
c
c  11/26/01 DMAS: this was a copy of task_energy.F 
c  changed 'energy' tasks to 'et', changed 'energy' variable to 'Vab',
c  removed ONIOM stuff, and theories other than SCF and DFT.
c  Also removed 'rtdb_clean' subroutine, and added task_et.o to GNUMakefile.
c
      logical status, ignore
      character*32 theory
      double precision cpu, wall, Vab
      logical task_et_doit
      external task_et_doit
c
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c
      ignore = rtdb_delete(rtdb, 'task:et')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_et: failed to invalidate status',0)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:et: theory not specified',0)
c
      call util_flush(6)
      status = task_et_doit(rtdb,theory,Vab)
      call util_flush(6)
      if (.not. rtdb_put(rtdb, 'task:et', mt_dbl, 1, Vab))
     $     call errquit('task_et: failed storing Vab',0)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_et: failed to store status',0)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_et: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_et: failed storing walltime',0)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('Vab', mt_dbl, Vab, 1)
         call ecce_print_module_exit('task et', 'ok')
      else
         call ecce_print_module_exit('task et', 'failed')
      endif
c         
      task_et = status
c
      end
c===========================================================================
      logical function task_et_doit(rtdb,theory,Vab)
      implicit none
c
c     This function actually calculates Vab
*
* $Id: task_et.F,v 1.4 2001-12-10 16:54:06 d3k958 Exp $
c
*
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
c
      integer rtdb
      character*(*) theory
      double precision Vab
c      double precision c_alpha,c_beta
c      integer nmo(2)            ! No. of vectors in each set
      integer nsets             ! No. of sets of vectors
      integer i, j
      integer k_vecs
      integer l_evals,k_evals
      integer unitno
      integer ga_s,g_vecsa,g_vecsb
      parameter (unitno = 67)
c      character*20 filename
      integer iset              ! No. (1,2) of set of vectors to read
c      integer g_vecs
c
c     Generic NWChem interface to compute Vab. SCF of DFT determinants
c     can be used.
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
c     
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:et (real)   - Vab (electron transfer matrix element)
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      logical movecs_read
      external movecs_read
      logical  scf
      external scf
      logical dft_energy
      external dft_energy
      double precision scr(mem_1e),s(nbf,nbf)
      double precision c_alpha(nbf,nbf),c_beta(nbf,nbf)
      double precision occ_alpha(nbf),occ_beta(nbf)
      double precision evals_alpha(nbf),evals_beta(nbf)
c
      logical status
c
      call ecce_print_module_entry('task et')
c
      task_et_doit = .false.
c
c     Set the theory in the database so that all is consistent for
c     the lower level modules.  Should not really have to do this,
c     but this value is used (appropriately and inappropriately) and
c     things are not yet disentangled.
c
c      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
c     $     call errquit('task_et: setting theory?',0)
c
c
c
c=======================================================
c  1. use the integral API to get the atomic overlap
c     matrix, S
c=======================================================
c
c     load geometry and symmetry info
c     -------------------------------
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
c
c     load the basis set and get info about it
c     ----------------------------------------
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $        call errquit('scf_init: no ao basis set', 0)
c
c     get nbf
c     -------
c
      if (.not.bas_numbf(basis,nbf))
     $     call errquit('bas_numbf failed',0)
c
c     initialize the API 
c     ------------------
c
      call int_init(rtdb, 1, basis)
      call int_acc_std()
c      call int_mem_print()
c
c     normalize the basis
c     -------------------
c
      call int_normalize(rtdb,basis)
c
c     calculate the overlap integrals and store in s(nbf,nbf)
c     ------------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,ga_s))
     $        call errquit('task_et ga_create failed', 0)
c       call int_1e_ga(basis, basis, ga_s, 'overlap', .false.)  
c
       write(6,*) 'overlap integrals'
       write(6,*) '-----------------'
       do i = 1,nbf
       do j = 1,nbf
         call int_1eov(basis,i,basis,j,mem_1e,scr,isz_1e,s(i,j))
       enddo
         write(6,1112) i,(s(i,j),j=1,nbf)
       enddo
c
c      terminate API
c      -------------
       call int_terminate
c
 1112 FORMAT (i2,6x,4(F8.5,5x))
c
c=======================================================
c     2. read MO vectors from files.
c     Read mo vectors for state from file 'movecsa.mo'.
c     Later, read from 'movecsb.mo' the same way.
c======================================================
c
c     allocate occ_alpha & beta, containing 0 or 1 for virtual and occupied orbitals
c     -----------------------------------------------------------------
c
      if (.not.ma_push_get(mt_dbl,nbf,'occ_alpha',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ_alpha',0)
c
      if (.not.ma_push_get(mt_dbl,nbf,'occ_beta',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ_beta',0)
c
c
c     allocate evals_alpha & beta, vector containing eigenvalues
c     ---------------------------------------------
      if (.not.ma_push_get(mt_dbl,nbf,'evals_alpha',l_evals,k_evals))
     $         call errquit('ma_push_get failed for evals_alpha',0)
c
      if (.not.ma_push_get(mt_dbl,nbf,'evals_beta',l_evals,k_evals))
     $         call errquit('ma_push_get failed for evals_beta',0)
c
c     allocate arrays for MO coefficients, c_alpha and c_beta
c     -------------------------------------------------------
      if (.not.ga_create(mt_dbl,nbf,nbf,'c_alpha',0,0,g_vecsa))
     $         call errquit('ga_create failed for c_alpha',0)
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'c_beta',0,0,g_vecsb))
     $         call errquit('ga_create failed for c_beta',0)
c
c     read occ, evals, mo coefficients from files
c     -------------------------------
      if(.not.movecs_read('movecsa.mo',1,occ_alpha,evals_alpha,g_vecsa))
     $   call errquit('movecs_read failed for c_alpha',0)
c
      if(.not.movecs_read('movecsa.mo',2,occ_beta,evals_beta,g_vecsb))
     $         call errquit('movecs_read failed for c_beta',0)
c
      write(6,*) 'occ_alpha',(occ_alpha(i),i=1,nbf)
      write(6,*) 'occ_beta ',(occ_beta(i),i=1,nbf)
      write(6,*) 'evals alpha',(evals_alpha(i),i=1,nbf)
      write(6,*) 'evals beta ',(evals_beta(i),i=1,nbf)
      write(6,*) 
c
      call ga_get(g_vecsa, 1, nbf, 1, nbf, dbl_mb(k_vecs), nbf)
      call ga_get(g_vecsb, 1, nbf, 1, nbf, dbl_mb(k_vecs), nbf)
c
      do i = 1,nbf
      do j = 1,nbf
         call ga_get(g_vecsa, 1, nbf, i, i,dbl_mb(k_vecs), nbf)
         c_alpha(i,j) = dbl_mb(k_vecs+j-1)
      enddo
c         write (6,1111) i,(dbl_mb(k_vecs+j),j=0,nbf-1)
      enddo
c
      do i = 1,nbf
      do j = 1,nbf
         call ga_get(g_vecsb, 1, nbf, i, i,dbl_mb(k_vecs), nbf)
         c_beta(i,j) = dbl_mb(k_vecs+j-1)
      enddo
c         write (6,1111) i,(dbl_mb(k_vecs+j),j=0,nbf-1)
      enddo
c
         write (6,*) 'alpha mo vectors'
         write (6,*) '----------------'
         do i = 1,nbf
           write(6,1111) i,(c_alpha(i,j),j=1,nbf)
         enddo
c
c
         write (6,*) 'beta mo vectors'
         write (6,*) '----------------'
         do i = 1,nbf
           write(6,1111) i,(c_beta(i,j),j=1,nbf)
         enddo
c
 1111    format(i2,1x,1x,4(f13.9,1x))
c
c
c      calculate Vab
c      -------------
c
       Vab=1.0d0
       write(6,*)
       write(6,*)'Vab = ',Vab
c
      status=.true.
      task_et_doit = status
c
      end

