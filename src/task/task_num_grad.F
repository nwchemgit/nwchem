C> \file task_num_grad.F
C> Calculate gradients numerically
C
C> \ingroup task
C> @{
C>
C> \brief The driver for numerical gradient evaluations
C>
C> \return 
C> - `.TRUE.` if the numerical gradient was calculated
C>    successfully
C> - `.FALSE.` otherwise.
C>
      logical function task_num_grad(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb !< [Input] The RTDB handle
*
* $Id$
*
#include "inp.fh"
      logical analytic
      logical task_energy, num_grad
      external task_energy, num_grad
      character*128 ename, gname
      character*32 theory
      logical ocons,omm
      logical cons_active
      external cons_active
c
c check for analytic user input 
c
      if (rtdb_get(rtdb,'task:analytic', mt_log, 1, analytic)) then
         if (analytic) then
            call errquit(
     $           'task_num_grad:asked for analytic but there is none',0,
     &       RTDB_ERR)
         endif
      endif
      
      if( .not. rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.
      ocons = cons_active(rtdb)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('num_grad: theory missing?',0, RTDB_ERR)
c
c     Nasty little hack ... CCSD module produces energy for multiple theories 
c
      if ((theory.eq.'ccsd(t)') .or. (theory.eq.'ccsd+t(ccsd)'))
     $     theory = 'ccsd'
c
c     TCE module also covers multiple theories
c
      if ((theory.eq.'uccd').or.
     1    (theory.eq.'ulccd').or.
     2    (theory.eq.'uccsd').or.
     3    (theory.eq.'ulccsd').or.
     4    (theory.eq.'uccsdt').or.
     5    (theory.eq.'uccsdtq'))
     6   theory = 'tce'
      if ((theory.eq.'uqcisd').or.
     1    (theory.eq.'ucisd').or.
     2    (theory.eq.'ucisdt').or.
     3    (theory.eq.'ucisdtq'))
     4   theory = 'tce'
      if ((theory.eq.'umbpt2').or.
     1    (theory.eq.'umbpt3').or.
     2    (theory.eq.'umbpt4').or.
     3    (theory.eq.'ump2').or.
     4    (theory.eq.'ump3').or.
     5    (theory.eq.'ump4'))
     6   theory = 'tce'
c
c     Construct the names for the energy and gradient RTDB entries
c
      write(ename,'(a,'':energy'')') theory(1:inp_strlen(theory))
      write(gname,'(a,'':gradient'')') theory(1:inp_strlen(theory))
      if(ocons.or.omm) then
        ename ="task:energy"
      end if
C
      task_num_grad = num_grad(rtdb, task_energy,
     $   ename(1:inp_strlen(ename)), gname(1:inp_strlen(gname)),
     $   theory(1:inp_strlen(theory)), .true.)
c
      end
C> @}
C>
C> \brief Calculate the gradient numerically
C>
C> Calculate the gradient of the function `func(rtdb)` with respect
C> to nuclear displacements using finite differences. The gradient
C> is returned in the runtime database with key `gname`. The 
C> algorithm used takes advantage of the symmetry of the molecule
C> if available. 
C>
C> The expression to be differentiated is passed in as the argument
C> `func` which is a function defined as `logical func(rtdb)`. This
C> function takes the geometry and other necessary information from
C> the runtime database. It puts its value (think of it as the energy)
C> back onto the runtime database with the key `ename`.
C>
C> If `osavestate` is `.TRUE.` then the save/restore state routines
C> can and should be used between finite difference steps and at
C> the end.
C>
C> \return
C> - `.TRUE.` if successful
C> - `.FALSE.` otherwise.
C>
      logical function num_grad(rtdb, func, ename, gname, theory,
     $     osavestate)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "inp.fh"
      integer rtdb          !< [Input] The RTDB handle
      logical func          !< [Input] The function to evaluate 
                            !< declared as `logical func(rtdb)`
      external func
      character *(*) ename  !< [Input] The RTDB key for the function
                            !< values
      character *(*) gname  !< [Input] The RTDB key for the function
                            !< gradient
      character *(*) theory !< [Input] The string specifying the theory
                            !< also used to construct the RTDB key
                            !< for the energy at the reference geometry
      logical osavestate    !< [Input] Use the save state routine
                            !< - .TRUE. use `task_save_state` and 
                            !<   `task_restore_state` routines
                            !< - .FALSE. omit the `task_*_state`
                            !<   routines
c
c     Compute derivatives of function(rtdb) w.r.t. nuclear displacements
c     using numerical finite difference, returning the gradient in rtdb(gname).
c
c     logical func(rtdb)
c     .   Takes the geometry (and other info as necessary) from the
c     .   database and puts its value (think of it as the energy) in
c     .   rtdb(ename).
c
c     If (osavestate) then the save/restore state routines can & should
c     be used between finite difference steps and at the end.
c
c     Uses symmetry and projects out rotations/translations. 
c
      integer geom
      integer nvec, ivec, i, j
      integer l_vec, k_vec, l_tmp, k_tmp, l_coord, k_coord
      integer l_ncoord, k_ncoord, l_grad, k_grad, l_q, k_q, l_t, k_t
      integer natoms, natoms3, ecce_old_print, ecce_junk_print
      logical oprint_mode, oprint_step, oprint_info, osaved
      logical ofivept
      character*128 name, key, tmp
      character*16 thetag
      double precision eref, eplus, eminus, grad, crd(3), q
      double precision eplus2, eminus2, grad3
      double precision delta    ! finite-difference step
      logical use_groups
      logical task_energy, task_save_state, task_restore_state, 
     $     task_delete_state, fd_step_energy
      external task_energy, task_save_state, task_restore_state, 
     $     task_delete_state, fd_step_energy
c
      integer ipcount
      logical rtdb_status
      logical util_sgstart
      external util_sgstart
      integer ncpus
      integer mygroup
      integer num_groups
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
      integer util_sgroup_numgroups
      external util_sgroup_numgroups

      ncpus=ga_nnodes()
      if (.not.rtdb_get(rtdb,'gradient:usegroups',MT_LOG,1,use_groups))
     $   use_groups = .false.
      if (ncpus.ne.1 .and. use_groups) then
       if(.not.util_sgstart(rtdb)) 
     $  call errquit('Error in SUBGROUPS init', 0, RTDB_ERR) 
       mygroup = util_sgroup_mygroup()
       num_groups = util_sgroup_numgroups()
       if (num_groups .eq. 1) call util_sgend(rtdb) ! Not really doing groups
      else ! Do not want to grab any python group info
       mygroup = 1
       num_groups = 1
      end if
      ipcount=mygroup-2
c
      num_grad = .false.
      osaved = .false.
      ofivept = .false.
c
      call ecce_print_control(0, ecce_old_print) ! Disable ECCE printing
c
      call util_print_push()
      call util_print_rtdb_load(rtdb, 'numg')
      oprint_info = util_print('information', print_low) .and. 
     $     ga_nodeid().eq.0
      oprint_mode = util_print('modes', print_high) .and. 
     $     ga_nodeid().eq.0
      oprint_step = util_print('step', print_low) .and. 
     $     ga_nodeid().eq.0
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('num_grad: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('num_grad: no geometry ', 0, RTDB_ERR)
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('num_grad: geom_ncent ?', geom, GEOM_ERR)
      natoms3 = natoms*3
      if (.not. ma_push_get(mt_dbl, natoms3*natoms3, 'tng: vectors',
     $     l_vec, k_vec)) call errquit
     $     ('num_grad: insufficient memory (1) ', natoms3**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, natoms3, 'tng: grad',
     $     l_grad, k_grad)) call errquit
     $     ('num_grad: insufficient memory (2) ', natoms3, MA_ERR)
      if (.not. ma_push_get(mt_dbl, natoms3, 'tng: coord',
     $     l_coord, k_coord)) call errquit
     $     ('num_grad: insufficient memory (3) ', natoms3, MA_ERR)
      if (.not. ma_push_get(mt_dbl, natoms3, 'tng: ncoord',
     $     l_ncoord, k_ncoord)) call errquit
     $     ('num_grad: insufficient memory (4) ', natoms3, MA_ERR)
      if (.not. ma_push_get(mt_dbl, natoms, 'tng: charge',
     $     l_q, k_q)) call errquit
     $     ('num_grad: insufficient memory (5) ', natoms, MA_ERR)
      if (.not. ma_push_get(mt_byte, natoms*16, 'tng: tags',
     $     l_t, k_t)) call errquit
     $     ('num_grad: insufficient memory (6) ', natoms, MA_ERR)
      if (.not. ma_push_get(mt_dbl, natoms3*natoms3, 'tng: temp',
     $     l_tmp, k_tmp)) call errquit
     $     ('num_grad: insufficient memory (7) ', natoms**2, MA_ERR)
c
      if (.not. geom_cart_get(geom, natoms, byte_mb(k_t), 
     $     dbl_mb(k_coord), dbl_mb(k_q))) 
     $     call errquit('num_grad: bad geom',0, GEOM_ERR)
c
c     Save the name, if any, of the reference geometry so we can
c     restore it upon completion.  Also, save the geometry itself
c     so that it can be recovered.
c
      if (.not. rtdb_cget(rtdb,'geometry', 1, name))
     $     name = 'geometry'
      if (.not. geom_rtdb_store(rtdb, geom, 'numg reference')) 
     $     call errquit('num_grad: store of geom failed',0, GEOM_ERR)
c
      if (oprint_info) then
         call util_print_centered(6, 
     $        'NWChem Numerical Gradients Module', 
     $        40, .true.)
         write(6,*)
         write(6,*)
*         write(6,*) ' Reference coordinates (a.u.) ',
*     $        name(1:inp_strlen(name))
*         call output(dbl_mb(k_coord), 1, 3, 1, natoms, 3, natoms, 1)
*         write(6,*)
      endif
c
      call fd_make_vectors(rtdb, geom, natoms, natoms3,
     $     dbl_mb(k_coord), dbl_mb(k_tmp),
     $     dbl_mb(k_vec), nvec, oprint_mode)
c
      if (oprint_info) then
         write(6,987) nvec
 987     format(' No. of totally-symmetric internal modes =', i4/)
      endif
c
      if (.not. ma_pop_stack(l_tmp)) call errquit('tng: ma corrupt',0,
     &       MA_ERR)
c
      if (osavestate) then
c     
c     Compute the energy at the reference geometry, save for restart
c     
         if (.not. task_energy(rtdb)) goto 100 ! Clean up and return
         if (.not. rtdb_get(rtdb, ename, mt_dbl, 1, eref)) then
            if (ga_nodeid() .eq. 0) 
     $           write(6,*) ' num_grad: no energy in database'
            goto 100            ! Clean up and return
         endif
c     
         if (oprint_info) then
            write(6,11) eref
 11         format(' Reference energy ', f20.8/)
            call util_flush(6)
         endif
         if (.not. task_save_state(rtdb, 'numg')) goto 100 
         osaved = .true.
      endif
c
      if (.not.rtdb_get(rtdb,'tng:fivept', mt_log, 1,
     &    ofivept)) ofivept = .false.
c
      if (.not. rtdb_get(rtdb, 'tng:delta', mt_dbl, 1, delta)) 
     $     delta = 0.01d0
      if (oprint_info) then
         write(6,12) delta, ofivept
 12      format(' Finite difference step ',1p,d12.5/
     $          ' Use five point formula ', l1/)
         call util_flush(6)
      endif
c
      call dfill(natoms3, 0.0d0, dbl_mb(k_grad), 1)
c
      do ivec = 1, nvec
       if (num_groups.gt.1) then
         ipcount=ipcount+1
         if (mod(ipcount,num_groups).ne.0) goto 1234
         write (6,*) 'mygroup=',mygroup,' ivec=',ivec
       endif
c
         if (ofivept) then
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t), 2.0d0*delta, eplus2, 
     $           osavestate, ename, func)) 
     $           goto 100
            if (oprint_step) then
               write(6,'('' step '',i3, '' +2 '',f20.8)') ivec, eplus2
               call util_flush(6)
            endif
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t), delta, eplus, 
     $           osavestate, ename, func))
     $           goto 100
            if (oprint_step) then
               write(6,'('' step '',i3, '' +1 '',f20.8)') ivec, eplus
               call util_flush(6)
            endif
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t), -delta, eminus, 
     $           osavestate, ename, func)) 
     $           goto 100
            if (oprint_step) then
               write(6,'('' step '',i3, '' -1 '',f20.8)') ivec, eminus
               call util_flush(6)
            endif
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t), -2.0d0*delta, eminus2, 
     $           osavestate, ename, func)) 
     $           goto 100
            if (oprint_step) then
               write(6,'('' step '',i3, '' -2 '',f20.8)') ivec, eminus2
               call util_flush(6)
            endif
            grad3 = (eplus-eminus) / (2.0d0*delta)
            grad  = (eplus-eminus)*(2.0d0/(3.0d0*delta)) - 
     $           (eplus2-eminus2)/(12.0d0*delta)
            if (oprint_step) then
               write(6,141) grad, grad3
 141           format('            gradient ',7x,f9.6,
     $              ' (3pt was ',f9.6,')')
               call util_flush(6)
            endif
         else
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t), delta, eplus, 
     $           osavestate, ename, func)) 
     $           goto 100
            if (oprint_step) then
               write(6,13) ivec, eplus
 13            format(' step ', i5,' energy ',f20.8)
               call util_flush(6)
            endif
            if (.not.fd_step_energy(rtdb, geom, natoms, dbl_mb(k_coord),
     $           dbl_mb(k_ncoord), dbl_mb(k_vec+(ivec-1)*natoms3), 
     $           dbl_mb(k_grad), dbl_mb(k_q),
     $           byte_mb(k_t),-delta, eminus, 
     $           osavestate, ename, func))
     $           goto 100
            if (oprint_step) then
               write(6,13) -ivec, eminus
               call util_flush(6)
            endif
c     
            grad = (eplus - eminus) / (2.0d0*delta)
            if (oprint_step) then
               write(6,14) grad
 14            format('            gradient ',7x,f9.6)
               call util_flush(6)
            endif
         endif
c     
c     logic should be here to diagnose too large/small energy change
c     
         call daxpy(natoms3, grad, dbl_mb(k_vec+(ivec-1)*natoms3), 1,
     $        dbl_mb(k_grad), 1)
c     
 1234 enddo
c
c     Symmetrize the gradient to get rid of symmetry breaking noise
c
      call sym_grad_symmetrize(geom, dbl_mb(k_grad))
c     
      if (oprint_info) then
         tmp = theory
         call inp_ucase(tmp)
         write(6,1000) tmp(1:inp_strlen(tmp)),
     $        'x','y','z','x','y','z'
         do i = 1, natoms
            if (.not. geom_cent_get(geom, i, thetag, crd, q))
     $           call errquit('task:numg geom_cent_get?',0, GEOM_ERR)
            write(6,2000) i, thetag,
     $           (dbl_mb(k_coord+3*(i-1)+j),j=0,2),
     $           (dbl_mb(k_grad +3*(i-1)+j),j=0,2)
         enddo
         write(6,*)
 1000    format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $        'coordinates',
     $        24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000    format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
         write(6,*)
         call util_flush(6)
      endif
      if (num_groups.gt.1) then
        if (ga_nodeid().gt.0) then
          call dfill (natoms3,0.0d0,dbl_mb(k_grad),1)
        endif
        call util_sgend(rtdb)
        call ga_dgop(1,dbl_mb(k_grad),natoms3,'+')
c
c       Symmetrize the gradient to get rid of symmetry breaking noise
c
        call sym_grad_symmetrize(geom, dbl_mb(k_grad))

        write (6,*) "Print summed up gradient" 
        if (ga_nodeid().eq.0) then
           tmp = theory
           call inp_ucase(tmp)
           write(6,3000) tmp(1:inp_strlen(tmp)),
     $          'x','y','z','x','y','z'
           do i = 1, natoms
              if (.not. geom_cent_get(geom, i, thetag, crd, q))
     $             call errquit('task:numg geom_cent_get?',0, GEOM_ERR)
              write(6,4000) i, thetag,
     $             (dbl_mb(k_coord+3*(i-1)+j),j=0,2),
     $             (dbl_mb(k_grad +3*(i-1)+j),j=0,2)
           enddo
           write(6,*)
 3000      format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $          'coordinates',
     $          24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 4000      format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
        endif
      endif
c
c     store the gradient in gname
c
      if (.not. rtdb_put(rtdb, gname, mt_dbl, natoms3, dbl_mb(k_grad)))
     $     call errquit('num_grad: failed to store gradient',0,
     &       RTDB_ERR)
c
      if (osavestate) then
c     
c     restore the reference energy
c     
         key = ' '
         write(key,'(a,'':energy'')') theory(1:inp_strlen(theory))
         if (.not. rtdb_put(rtdb, key, mt_dbl, 1, eref))
     $        call errquit('num_grad: failed to store energy',0,
     &       RTDB_ERR)
      endif
c     
      num_grad = .true.
c
 100  continue
c
c     Restore state and clean up
c
      if (osavestate .and. osaved) then
         if (.not. task_restore_state(rtdb,'numg'))
     $        call errquit('num_grad: failed to restore state',0,
     &       RTDB_ERR)
         if (.not. task_delete_state(rtdb,'numg'))
     $        call errquit('num_grad: failed to delete state',0,
     &       RTDB_ERR)
      endif
c
c     restore any indirection of the geometry and delete the reference
c
      if (name .ne. 'geometry') then
         if (.not. rtdb_cput(rtdb,'geometry',1,name))
     $        call errquit('num_grad: rtdb corrupt',0, RTDB_ERR)
      endif
      if (.not. geom_rtdb_delete(rtdb, 'numg reference')) 
     $     call errquit('num_grad: delete of ref geom failed',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'scf:converged',mt_log,1,.false.))
     $        call errquit('num_grad: rtdb corrupt',1, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'dft:converged',mt_log,1,.false.))
     $        call errquit('num_grad: rtdb corrupt',2, RTDB_ERR)
c
      if (.not. ma_pop_stack(l_t))call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_q))call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_ncoord))call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_coord)) call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_grad)) call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_vec)) call errquit('tng: ma corrupt',0,
     &       MA_ERR)
      if (.not. geom_destroy(geom)) call errquit('tng: geom corrupt',0,
     &       GEOM_ERR)
      call util_print_pop()
c
      call ecce_print_control(ecce_old_print, ecce_junk_print)
c
      end
      logical function fd_step_energy(rtdb, geom, natoms, coords,
     $     ncoords, vec, grad, q, tags, delta, energy, 
     $     osavestate, ename, func)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb, geom
      integer natoms
      double precision coords(*),ncoords(*),vec(*), grad(*), q(*)
      character*16 tags(*)
      double precision delta
      double precision energy   ! [output]
      logical osavestate
      character *(*) ename
      logical func
      external func
c
      integer i, natoms3
      integer  i_array,l_array,uniquecent 
      logical status
      logical task_restore_state,geom_uniquecent_set
      external task_restore_state,geom_uniquecent_set
c
      fd_step_energy = .false.
      natoms3 = natoms*3
c         
      do i = 1, natoms3
         ncoords(i) = coords(i) + delta*vec(i)
      enddo
      call sym_grad_symmetrize(geom, ncoords) ! Paranoid enforcement
c
c     hack to keep old unique atoms list to avoid getting the C1 list
c     from geom_cart_set
c
        if (.not. geom_ncent_unique(geom,uniquecent))
     &     call errquit('bas_rtdb_do_load: geom_ncent_unique',211,
     &       RTDB_ERR)
      if (.not.MA_push_Get(mt_int,uniquecent, 'iarray',
     &     l_array, i_array))
     &     call errquit('tng: cannot allocate ',0, MA_ERR)
        if (.not. geom_uniquecent_get(geom,uniquecent,
     ,     int_mb(i_array)))
     &     call errquit('tng: geom_uniquecent_get',211,
     &       RTDB_ERR)

      if (.not. geom_cart_set(geom, natoms, tags, ncoords, q))
     $     call errquit('num_grad: bad geom',0, GEOM_ERR)
c
c     restoring old unique atoms list
c
        if (.not. geom_uniquecent_set(geom,uniquecent,
     ,     int_mb(i_array)))
     &     call errquit('bas_rtdb_do_load: geom_uniquecent_get',211,
     &       RTDB_ERR)
      if (.not. ma_pop_stack(l_array)) call errquit(
     $'tng: ma corrupt',0,  MA_ERR)

      if (.not. geom_rtdb_store(rtdb, geom, 'num. grad. geometry')) 
     $     call errquit('num_grad: store of geom failed',0, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'geometry', 1, 'num. grad. geometry'))
     $     call errquit('num_grad: rtdb corrupt',0, RTDB_ERR)
c
      if (osavestate) then
         if (.not. task_restore_state(rtdb,'numg'))
     $        call errquit('num_grad: restore state failed',0, RTDB_ERR)
      endif
      status = func(rtdb)
      fd_step_energy = status
      if ((.not. status) .and. (ga_nodeid().eq.0)) then
         write(6,*) 
         write(6,*) ' num_grad: !! displaced energy failed'
         write(6,*)
         call util_flush(6)
      endif
      if (.not. status) return
c
      if (.not. rtdb_get(rtdb, ename, mt_dbl, 1, energy))
     $     call errquit('num_grad: no energy in database',0, RTDB_ERR)
c
      if (.not. geom_rtdb_delete(rtdb, 'num. grad. geometry'))
     $     call errquit('num_grad: rtdb corrupt (geom) ',0, RTDB_ERR)
      if (.not. rtdb_delete(rtdb,'geometry'))
     $     call errquit('num_grad: rtdb corrupt',0, RTDB_ERR)
c
      end
C>
C> \brief Make a set of orthonormal totally symmetric vectors
C>
C> Make a set of orthonormal vectors that are
C> - totally symmetric
C> - orthogonal to rotations
C> - orthogonal to translations
C>
C> The vectors and the number of vectors is returned. These vectors
C> may be used as displacements of the atomic positions that preserve
C> the symmetry of the molecular geometry.
C>
      subroutine fd_make_vectors(rtdb, geom, natoms, natoms3, 
     $     coords, temp, 
     $     vectors, nvec, oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

      integer rtdb    !< [Input] The RTDB handle
      integer geom    !< [Input] The geometry handle
      integer natoms  !< [Input] The number of atoms
      integer natoms3 !< [Input] The number of nuclear coordinate
c
      double precision coords(3,natoms)          !< [Input] The nuclear
                                                 !< coordinates
      double precision temp(natoms3, natoms3)    !< [Scratch]
      double precision vectors(natoms3, natoms3) !< [Output] The vectors
      !< stored by column, i.e. `vectors(1:natoms3,ivec)` is one vector
      integer nvec    !< [Output] The number of vectors generated
      logical oprint  !< [Input] A print flag
c
c     Make a set of orthonormal vectors that are
c     1) totally symmetric
c     2) orthogonal to rotations
c     3) orthogonal to translations
c
c     Return the vectors and the number of vectors.
c
      integer i3, i, j, k 
cedo  integer info
      integer maxat, nwork, maxat3
      parameter (maxat = 100, nwork=12) ! Nwork must be at least 6
      parameter (maxat3 = 3*maxat)
      logical oactive(natoms)
      integer nactive
      double precision work(natoms3,nwork), eval(natoms3)
      double precision centroid(3), x, y, z, xx, yy, zz, fx
      integer ifirstmax, idamax
      external ifirstmax, idamax
      double precision ddot
      external ddot
      logical omm
c
c     Get list of active atoms.  If there are ANY frozen atoms then
c     we cannot use or translational, or rotational invariance
c     and we will just use the cartesian coordinates of the active atoms
c     plus symmetry (assuming that the frozen atoms preserve symmetry)
c
      call grad_active_atoms(rtdb, natoms, oactive, nactive)
      call dfill(natoms3*natoms3,0.0d0,vectors,1)
      if (nactive .ne. natoms) then
         do i = 1, natoms
            if (oactive(i)) then
               do k = 1,3
                  vectors((i-1)*3+k,(i-1)*3+k) = 1.0d0
               enddo
               nvec = nvec+3
            endif
         enddo
      else
         call dfill(natoms3,1.0d0,vectors,natoms3+1)
      endif
c
c     Project totally symmetric component
c
      do i = 1, natoms3
         call sym_grad_symmetrize(geom, vectors(1,i))
      enddo

      if( .not. rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.
c
      if(.not.omm) then
      if (nactive .eq. natoms) then
c     
c     Construct vectors corresponding to translations and rotations
c
         call dfill(3, 0.0d0, centroid, 1)
         do i = 1, natoms
            do k = 1, 3
               centroid(k) = centroid(k) + coords(k,i)/natoms
            enddo
         enddo
c     
         do k = 1, 3            ! x, y, z translations
            call dfill(natoms3, 0.0d0, work(1,k), 1)
            call dfill(natoms, sqrt(1.0d0/natoms), work(k,k), 3)
         enddo
         do k = 4, 6            ! x, y, z rotations
            do i = 1, natoms
               x = coords(1,i) - centroid(1)
               y = coords(2,i) - centroid(2)
               z = coords(3,i) - centroid(3)
               if (k .eq. 4) then
                  xx = 0.0d0
                  yy = -z
                  zz =  y
               else if (k .eq. 5) then
                  xx =  z
                  yy =  0.0d0
                  zz = -x
               else if (k .eq. 6) then
                  xx = -y
                  yy =  x
                  zz =  0.0d0
               else
c                 HvD: Shutting the compiler up about uninitialized 
c                 variables.
                  xx = 0.0d0
                  yy = 0.0d0
                  zz = 0.0d0
                  call errquit("fd_make_vectors: invalid k",k,UERR)
               endif
               i3 = (i-1)*3
               work(i3+1,k) = xx
               work(i3+2,k) = yy
               work(i3+3,k) = zz
            enddo
            do j = 1, k-1
               fx = ddot(natoms3, work(1,j), 1, work(1,k), 1)
               call daxpy(natoms3, -fx, work(1,j), 1, work(1,k), 1)
            enddo
            fx = sqrt(ddot(natoms3, work(1,k), 1, work(1,k), 1))
            if (fx . gt. 1d-6) then
               call dscal(natoms3, 1.0d0/fx, work(1,k), 1)
            else
               call dfill(natoms3, 0.0d0, work(1,k), 1)
            endif
         enddo
c     
         do k = 1, 6
            call sym_grad_symmetrize(geom, work(1,k))
            fx = sqrt(ddot(natoms3, work(1,k), 1, work(1,k), 1))
            if (fx . gt. 1d-6) then
               call dscal(natoms3, 1.0d0/fx, work(1,k), 1)
            else
               call dfill(natoms3, 0.0d0, work(1,k), 1)
            endif
         enddo
c     
*     if (oprint) then
*     write(6,*) ' symmetrized translation and rotation vectors'
*     call output(work, 1, natoms3, 1, 6, 3*maxat, 6, 1)
*     write(6,*)
*     endif
c     
c     Project trans/rotn from the symmetrized vectors
c     
         do k = 1, 6
            do i = 1, natoms3
               fx = ddot(natoms3, work(1,k), 1, vectors(1,i), 1)
               call daxpy(natoms3, -fx, work(1,k), 1, vectors(1,i), 1)
            enddo
         enddo
      endif
      end if
c
c     Get rid of the zero vectors and orthonormalize the rest
c     by diagonalizing the metric
c     
      call dgemm('t','n',natoms3,natoms3,natoms3,1.0d0, 
     $     vectors, natoms3, vectors, natoms3, 0.0d0, temp, natoms3)
      call util_jacobi(natoms3,temp,natoms3,eval)
cedo      call dsyev('v','u',natoms3,temp,natoms3,eval,work,
cedo     $     3*nwork*maxat,info)
cedo      if (info .ne. 0) call errquit('fd_make_vectors: dsyev failed',0)
*      write(6,*) ' Eigen values of the projected vectors metric '
*      call output(eval, 1, natoms3, 1, 1, natoms3, 1, 1)
*      write(6,*) ' Eigen vectors of the projected vectors metric '
*      call output(temp, 1, natoms3, 1, natoms3, natoms3, natoms3, 1)
c     
c     There should be nvec unit eigen values with all else zero
c     
      nvec = 0
      do i = 1, natoms3
         if (eval(i) .gt. 1d-3) then
            nvec = nvec + 1
            call dcopy(natoms3,temp(1,i),1,vectors(1,nvec),1)
            call dscal(natoms3,1.0d0/sqrt(eval(i)),vectors(1,nvec),1)
         endif
      enddo
      if (oprint) then
        write(6,*) ' Independent symmetric all-internal modes ',  nvec
        call output(vectors, 1, natoms3, 1, nvec, natoms3, nvec, 1)
        call util_flush(6)
      endif
c
c phase align each vector with a positive maximum magnitude value
c
c     do i = 1,nvec
c       k = ifirstmax(natoms3,vectors(1,i),1)
c       if (vectors(k,i).lt.0.0d00) then
c         call dscal(natoms3,-1.0d00,vectors(1,i),1)
c       endif
c     enddo
c
c     Phase align vectors so that their directions are body fixed
c     - Note that the diagonalization of the metric left the vectors
c       in random order. So we need to sort them as well as fixing
c       the phase. We will sort them in order of decreasing dot product
c       with the geometry.
c
      do i = 1, nvec
        eval(i) = ddot(natoms3,coords,1,vectors(1,i),1)
        if (eval(i).lt.0.0d0) then
          call dscal(natoms3,-1.0d0,vectors(1,i),1)
          eval(i) = -eval(i)
        endif
      enddo
      do i = 1, nvec - 1
        k = idamax(nvec-i+1,eval(i),1)
        k = i+k-1
        if (k.ne.i) then
c
c         swap the vectors
c
          xx = eval(i)
          eval(i) = eval(k)
          eval(k) = xx
          do j = 1, natoms3
            xx           = vectors(j,i)
            vectors(j,i) = vectors(j,k)
            vectors(j,k) = xx
          enddo
        endif
      enddo
c
      if (oprint) then
         write(6,*)
     &      ' Phase Corrected ',
     &      'Independent symmetric all-internal modes ',
     &      nvec
         call output(vectors, 1, natoms3, 1, nvec, natoms3, nvec, 1)
         call util_flush(6)
      endif
c
      end
      integer function ifirstmax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     but finds the first element within a predefined threshold
c
      double precision dx(*),dmax, dtmp, diff
      double precision threshold
      parameter (threshold = 1.0d-3)
      integer i,incx,ix,n
c
      ifirstmax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      ifirstmax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do i = 2,n
        dtmp = dabs(dx(ix))
        if(dtmp.lt.dmax) go to 5
        diff = dabs(dmax-dtmp)
        if (diff.gt.threshold) then
          ifirstmax = i
          dmax = dtmp
        endif
 5      ix = ix + incx
      enddo
      return
c
c        code for increment equal to 1
c
 20   continue
      dmax = dabs(dx(1))
      do  i = 2,n
        dtmp = dabs(dx(i))
         if(dtmp.lt.dmax) go to 30
         diff = dabs(dmax-dtmp)
         if (diff.gt.threshold) then
           ifirstmax = i
           dmax = dtmp
         endif
 30      continue
       enddo
      return
      end
C>
C> @}
