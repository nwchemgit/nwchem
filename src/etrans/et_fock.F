c $Id: et_fock.F,v 1.1 2002-02-19 21:38:23 d3k958 Exp $
       subroutine et_fock(rtdb,nbf,basis,geom,mem2,max2e,
     &                    noc,g_vecs,g_pa,g_pb)
       implicit none
c
c      calculates the two-electron contribution to Vba, omega2
c      using the fock builder.
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer g_pa,g_pb,noc(2)
       integer nbf, geom, basis,max2e, mem2, rtdb
       integer int_theory, nfock
       integer g_p, g_f, g_vecs
c       integer l_gJ, k_gJ, l_gK, k_gK, l_scr, k_scr
c       integer ish, jsh, ksh, lsh, nsh
c       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c       integer i,j

       double precision jfac, kfac, tol2e
c       double precision omega2, omega2_val
c       double precision jfac, kfac, tol2e, val_p,val_f
c
       logical bas_numcon
       external bas_numcon
       logical ga_copy
       external ga_copy
       logical et_movecs_read
       external et_movecs_read
c
       character*3 theory
c
c
c
c  create arrays
c  =============
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'fock',0,0,g_f))
     $        call errquit('task_et ga_create failed', 0)
      call ga_zero(g_f)
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'density',0,0,g_p))
     $        call errquit('task_et ga_create failed', g_p)
      call ga_zero(g_p)
c
c  set fock_2e variables
c  =====================
c
      nfock = 1
      tol2e = 1.0d-7
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('et_calc: couldnt get task:theory',0)
c
      if (theory.eq.'hyb') then
        int_theory=3       
      else if (theory.eq.'dft') then
        int_theory=2       
      else if (theory.eq.'hf') then
        int_theory=1
      else if (theory.eq.'scf') then
        int_theory=1
      else
        call errquit('et_calc: could not understand theory',0)
      endif       
c      
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1,int_theory))
     $     call errquit('et_calc: put of bgj:scf_type failed',0)
c
      call scf_get_fock_param(rtdb,tol2e)
      call fock_force_direct(rtdb)    ! Force Fock build to be direct
c
c ==============================================
c first, try to reproduce the Falpha fock matrix
c from the nwchem output (scf debug printing is 
c needed to print it)
c=================================================
c
      call ga_print(g_vecs)
c
      call ga_dgemm('N','T',nbf,nbf,noc(1),1.d0,g_vecs
     &               ,g_vecs,0.d0,g_p)
c
      call ga_print(g_p)
c
      jfac = 1.0d0
      kfac = 0.0d0
c
      call fock_2e( geom, basis, nfock, jfac, kfac,
     $  tol2e, .false., g_p, g_f, .false. )
c
      call ga_print (g_f)
c
c
c
      return
      end

 
