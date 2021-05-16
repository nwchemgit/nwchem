      subroutine mcscf_occ2int( nbf, occ, nclosed, nact )
c
c  Convert occupation vector to 
c  number of closed and active orbitals
c
      implicit none
      integer nbf
      double precision occ(nbf)
      integer nclosed
      integer nact
c
      integer i
c
      nclosed = 0
      nact = 0
      do i=1,nbf
        if (occ(i).eq.2.d0) then
          nclosed = nclosed + 1
        else if (occ(i).ne.0.d0) then
          nact = nact + 1
        endif
      enddo
      return
      end


      









c
c  Convert occupation vector
c  to number of electrons (trace 1pdm)
c
      subroutine mcscf_occ2elec( nbf, occ, nelec )
      implicit none
      integer nbf, nelec
      double precision occ(nbf)
      integer i
      double precision xelec

      xelec = 0.d0
      do i=1,nbf
        xelec = xelec + occ(i)
      enddo
      nelec = xelec
      return
      end







c
c Convert the vector representation of the MCSCF
c parameters to a (nbf x nbf) matrix representation
c Accumulates into the matrix (does not copy!)
c
       subroutine mcscf_vec2mat( nbf, nclosed, nact, a, g_v, b, g_m )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision b,a                         ! [input] Additive factors
       integer g_v                                  ! [input] Vector representation
       integer g_m                                  ! [output] Matrix representation
c
c
       integer nvir, vlen
       integer voff, aoff, aend, xoff, xend
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       call ga_dadd_patch( a, g_v, 1, (nvir*nclosed), 1, 1,
     $                     b, g_m, voff, nbf, 1, nclosed,
     $                        g_m, voff, nbf, 1, nclosed )
       if (nact.gt.0) then
         xoff = nclosed*nvir + 1
         xend = nclosed*nvir + nact*nvir
         call ga_dadd_patch( a, g_v, xoff, xend, 1, 1,
     $                       b, g_m, voff, nbf, aoff, aend,
     $                          g_m, voff, nbf, aoff, aend )
         xoff = xend + 1
         xend = xend + nclosed*nact
         call ga_dadd_patch( a, g_v, xoff, xend, 1, 1,
     $                       b, g_m, aoff, aend, 1, nclosed,
     $                          g_m, aoff, aend, 1, nclosed )
       endif
c
       return
       end

       


c
c Convert the matrix representation of the MCSCF
c parameters to vector representation (complement of mcscf_vec2mat)
c Accumulates into the vector (does not copy!)
c
c
       subroutine mcscf_mat2vec( nbf, nclosed, nact, a, g_m, b, g_v )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision b,a                         ! [input] Additive factors
       integer g_m                                  ! [input] Matrix representation
       integer g_v                                  ! [output] Vector representation
c
c
       integer nvir, vlen
       integer voff, aoff, aend, xoff, xend
c
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       xend = nvir*nclosed
       call ga_dadd_patch( a, g_m, voff, nbf, 1, nclosed,
     $                     b, g_v, 1, xend, 1, 1,
     $                        g_v, 1, xend, 1, 1)
       if (nact.gt.0) then
         xoff = xend + 1
         xend = xend + nact*nvir
         call ga_dadd_patch( a, g_m, voff, nbf, aoff, aend,
     $                       b, g_v, xoff, xend, 1, 1,
     $                          g_v, xoff, xend, 1, 1)
         xoff = xend + 1
         xend = vlen
         call ga_dadd_patch( a, g_m, aoff, aend, 1, nclosed,
     $                       b, g_v, xoff, xend, 1, 1,
     $                          g_v, xoff, xend, 1, 1)
       endif
       return
       end









c
c  Create closed and active AO densities from
c  occupation numbers and 1-pdm in active space
c
       subroutine mcscf_occ2dens( nbf, nclosed, nact, dm1,
     $                            g_movecs, g_cdens, g_adens )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
       integer nbf
       integer nclosed
       integer nact
       double precision dm1(nact,nact)
       integer g_movecs
       integer g_cdens
       integer g_adens
c
       integer aoff, aend
       integer g_dm1
c
c
c
       aoff = nclosed + 1
       aend = nclosed + nact
*ga:1:0
       if (.not.ga_create(MT_DBL,nact,nact,'1PDM',nact,0,g_dm1))
     $   call errquit('mcscf_occ2dens: cannot create 1PDM',0, GA_ERR)
       if (ga_nodeid().eq.0)
     $   call ga_put(g_dm1,1,nact,1,nact,dm1,nact)
       call ga_sync()
c
c Active AO density (use closed shell dens for temp)
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_dm1, 1, nact, 1, nact,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_cdens, 1, nact, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_cdens, 1, nact, 1, nbf,
     $                       g_adens, 1, nbf, 1, nbf )
c
c
       if (.not.ga_destroy(g_dm1))
     $   call errquit('mcscf_occ2dens: cannot destroy',0, GA_ERR)
c
c Closed AO density
c
       call ga_dgemm( 'n', 't', nbf, nbf, nclosed, 2.d0,
     $                g_movecs, g_movecs, 0.d0, g_cdens )
       return
       end










      subroutine mcscf_symmetrize_2pdm( nact, dm2, dmsy2 )
      implicit none
      integer nact
      double precision dm2(nact,nact,nact,nact)
      double precision dmsy2(nact,nact,nact,nact)
      integer i,j,k,l
      
      do i=1,nact
        do j=1,nact
          do k=1,nact
            do l=1,nact
              dmsy2(i,j,k,l) = 0.5d0*(dm2(i,j,k,l) + dm2(i,j,l,k))
            enddo
          enddo
        enddo
      enddo
      return
      end

      




















      subroutine mcscf_getmovecs( rtdb, geom, basis, nbf, occ,
     $                            evals, g_movecs )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
c
c  Retrieve MO vectors either from 
c      i)  hcore or atomic guess 
c      ii) projection from another basis
c      ii) restart
c  This is the MCSCF version of "scf_vectors_guess" routine
c  except it does *not* store stuff in common blocks
c
      integer rtdb
      integer geom
      integer basis
      integer nbf
      double precision occ(nbf)
      double precision evals(nbf)
      integer g_movecs
c
c     local variables
c
C      integer nclosed, nact
      integer nelec
      integer g_aodens                         ! Atomic density guess
      character*255 title_vecs, bas_vecs       ! Basis names from movecs_in
      character*255 trans, name                ! Basis names from rtdb
      character*255 movecs_in                  ! Input/output filenames
      character*255 movecs_out
      character*20 scftype_vecs
      integer i, nbf_vecs, nsets, nmo_vecs(2)
c$$$      integer l_oph, k_oph
      integer l_v, k_v
      integer j, j1, j2, nn
      integer nmo
      logical oskel, oprint, status, ovecs, oldvectors
      double precision tol2e, charge, nuclear_charge
c$$$      double precision vphase
c
c     externals
c
      integer ga_create_atom_blocked
      logical movecs_read_header
      logical movecs_read
      external ga_create_atom_blocked
      external movecs_read_header
      external movecs_read
c
      oldvectors = .false.
c
c  Get misc info
c
      oprint = (ga_nodeid().eq.0).and.
     $         (util_print('info',print_default))
      if (.not. rtdb_get(rtdb, 'mcscf:tol2e', MT_DBL, 1, tol2e))
     $    tol2e = 1d-9
      if (.not.rtdb_get(rtdb, 'mcscf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
      if (.not. bas_name(basis, name, trans))
     $     call errquit('mcscf: bas_name?', 0, BASIS_ERR)
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('mcscf: negative no. of electrons ?', nelec, INPUT_ERR)
c
c  Get filenames
c
      if (.not.rtdb_cget(rtdb,'mcscf:input vectors',1,movecs_in)) then
         if (.not.rtdb_cget(rtdb,'scf:input vectors',1,movecs_in)) 
     $        movecs_in = 'atomic'
      endif
c
c  Starting MO's determined from movecs_in string
c
      call ga_zero(g_movecs)
      call scf_lindep(rtdb, geom, basis, g_movecs, nmo)
      if (nmo.ne.nbf)
     $  call errquit('mcscf:linear dependent AO basis',0, BASIS_ERR)
      if (movecs_in .eq. 'hcore' .or. movecs_in .eq. 'atomic') then
c
c  Hcore or atomic guess
c
         g_aodens = ga_create_atom_blocked(geom, basis, 'Guess density')
         call ga_zero(g_aodens)
         if (movecs_in .eq. 'atomic') then
            call guess_dens(rtdb, geom, basis, g_aodens)
         endif
         call rhf_dens_to_mo(rtdb, geom, basis, nelec, nbf, nbf,
     $                       tol2e, g_aodens, g_movecs, oskel, evals )
         if (.not. ga_destroy(g_aodens)) call errquit
     $        ('scf_vectors_guess: destroy of density failed',0,
     &       GA_ERR)
      else if (movecs_in .eq. 'project') then
c
c  Project from a smaller basis set
c
         call movecs_project(rtdb, nmo, geom, basis, g_movecs,1)
      else
c
c  Read vectors from file movecs_in and see if 
c  they make sense
c
         call util_file_name_resolve(movecs_in, .false.)
         status = movecs_read_header( movecs_in, title_vecs, bas_vecs,
     $        scftype_vecs,
     $                                nbf_vecs, nsets, nmo_vecs, 2)
         if (status) then
c     
c  Check that:
c  1) nbf from file = current nbf
c  2) name of basis from file = current basis name
c     
           if (nbf .eq. nbf_vecs .and. ((bas_vecs.eq.trans).or.
     $                                  (bas_vecs.eq.'******'))) then
             if (.not.movecs_read( movecs_in, 1, occ, evals, g_movecs ))
     $            call errquit('mcscf_getmovecs: movecs read',0,
     &                         DISK_ERR)
           else
             if (ga_nodeid().eq.0) write(6,944)
 944         format('Restart movecs mismatch with input')
             status = .false.
           endif
           if (oprint .and. title_vecs.ne.' ') then
             write(6,77) title_vecs(1:inp_strlen(title_vecs))
 77          format(/' Loading old vectors from job with title :'//a/)
             call util_flush(6)
           endif
         endif
         if (.not. status) call errquit('rhf: failed to load vectors',
     &                                    0,DISK_ERR)
         oldvectors = .true.
      endif
c
c     Optional MO vector phase adjustment (for getting exactly the same vectors
c                                          as a previous calculation)
c
c$$$      if (.not.ma_push_get(MT_DBL, nbf,'ophase', l_oph, k_oph))
c$$$     $   call errquit('mcscf_getmovecs: cannot allocate',0)
c$$$       if (rtdb_get(rtdb,'vecphase',MT_DBL,nbf,dbl_mb(k_oph))) then
c$$$         do i=1,nbf
c$$$           vphase = sign(1.d0,dbl_mb(k_oph+i-1))
c$$$           call ga_dscal_patch( g_movecs, 1, nbf, i, i, vphase )
c$$$         enddo
c$$$         if (util_print('input phase',print_high)) then
c$$$           if (ga_nodeid().eq.0) then
c$$$             write(6,831)
c$$$ 831         format('MO vectors phase adjustment')
c$$$             write(6,832) (dbl_mb(k_oph+i-1),i=1,nbf)
c$$$ 832         format(12f6.1)
c$$$           endif
c$$$         endif
c$$$       endif
c$$$       if (.not.ma_pop_stack(l_oph))
c$$$     $   call errquit('mcscf_getmovecs: failed pop stack',0)
c
c  Import MO vectors if file exist
c
      if (ga_nodeid() .eq. 0) then
         inquire(file='vecs',exist=ovecs)
         if (ovecs) then
            if (util_print('information',print_default)) then
               if (ga_nodeid().eq.0) write(6,670)
 670           format(' Importing MO vectors from external file')
            endif
            if (.not.ma_push_get(MT_DBL, nbf,'ophase', l_v, k_v))
     $           call errquit('mcscf_getmovecs: cannot allocate',0,
     &       MA_ERR)
            open(unit=19,file='vecs',status='old',form='formatted')
            read(19,*) nn
            do i=1,nbf
               do j1=1,nbf,6
                  j2 = min(nbf,(j1+5))
                  read(19,913) (dbl_mb(k_v+j-1),j=j1,j2)
 913              format(6e20.11)
               enddo
               call ga_put(g_movecs,1,nbf,i,i,dbl_mb(k_v),nbf)
            enddo
            close(19)
            if (.not.ma_pop_stack(l_v))
     $           call errquit('mcscf_getmovecs: failed pop stack',0,
     &       MA_ERR)
         endif
      endif
      call ga_sync()
c
c
c     Apply any swaps of orbitals
c
      call movecs_swap(rtdb,'mcscf','mcscf',g_movecs,occ,evals)
c
c
c     Orthogonalize all input/guess vectors ... don't mix the occupied
c     with the virtuals.
c
      call ga_orthog_mos(basis, g_movecs)
c
c     Apply phase conventions
c
      call movecs_fix_phase(g_movecs)
c
c
c
      if (util_print('initial vectors', print_high)) then
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,'Initial MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         endif
         call ga_print(g_movecs)
      endif
c
c     Tie the input vectors for the next job to the output vectors
c     of this job
c
       if (.not.rtdb_cget(rtdb,'mcscf:output vectors',1,movecs_out))then
          if (oldvectors) then
             movecs_out = movecs_in
          else
             call util_file_name('movecs', .false., .false., movecs_out)
          endif
       endif

       call util_file_name_resolve(movecs_out, .false.)

       if (.not.rtdb_cput(rtdb,'mcscf:input vectors',1,movecs_out))
     $      call errquit('mcscf: failed resetting input vectors', 0,
     &       RTDB_ERR)
       if (.not.rtdb_cput(rtdb,'mcscf:output vectors',1,movecs_out))
     $      call errquit('mcscf: failed resetting output vectors', 0,
     &       RTDB_ERR)
c
c     Let the user know what is going on
c
       if (ga_nodeid().eq.0 .and. util_print('information',print_low))
     $      then
          write(6,17) movecs_in(1:inp_strlen(movecs_in)),
     $         movecs_out(1:inp_strlen(movecs_out))
 17       format(/
     $         '   input molecular orbitals = ', a/
     $         '  output molecular orbtials = ', a/)
          call util_flush(6)
       endif
c
      end


      subroutine mcscf_movecs_read( rtdb, basis, nbf, occ,
     $                              evals, g_movecs )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
c
c  Read MCSCF MO vectors
c
      integer rtdb
      integer basis
      integer nbf
      double precision occ(nbf)
      double precision evals(nbf)
      integer g_movecs
c
c     local variables
c
      character*255 title_vecs, bas_vecs       ! Names from movecs_in
      character*255 trans, name                ! Basis names from rtdb
      character*20 scftype_vecs
      character*255 v_in                       ! Input filename
      integer nbf_vecs, nsets, nmo, nmo_vecs(2)
      logical status
c
c     externals
c
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
c  Get filenames
c
      if (.not.rtdb_cget(rtdb,'mcscf:input vectors',1,v_in)) then
         call errquit('mcscf: failed to get name of input vectors',0,
     &       RTDB_ERR)
      endif
c
c
c
      if (.not. bas_name(basis, name, trans))
     $     call errquit('scf_init: bas_name?', 0, BASIS_ERR)
      nsets = 1
      nmo = nbf
c
c  Read vectors from file v_in and see if 
c  they make sense
c
      status = movecs_read_header( v_in, title_vecs, bas_vecs,
     $     scftype_vecs,
     $                             nbf_vecs, nsets, nmo_vecs, 2)
      if (status) then
c     
c  Check that:
c  1) nbf from file = current nbf
c  2) name of basis from file = current basis name
c     
        if (nbf .eq. nbf_vecs .and. ((bas_vecs.eq.trans).or.
     $                               (bas_vecs.eq.'******'))) then
          if (.not.movecs_read( v_in, 1, occ, evals, g_movecs ))
     $       call errquit('mcscf_movecs_read: movecs read failed',0,
     &       DISK_ERR)
        else
          if (ga_nodeid().eq.0) write(6,944)
 944      format('Restart movecs mismatch with input')
          status = .false.
        endif
      endif


****      call movecs_swap(rtdb,'mcscf','rohf',g_movecs,occ,evals)

      if (.not. status)
     $  call errquit('mcscf_movecs_read: failed to load vectors',0,
     &       DISK_ERR)
c
c
c
      return
      end



      subroutine mcscf_movecs_write( rtdb, basis, nbf, occ,
     $                               evals, g_movecs )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
c
c  Save MCSCF MO vectors
c
      integer rtdb
      integer basis
      integer nbf
      double precision occ(nbf)
      double precision evals(nbf)
      integer g_movecs
c
c     local variables
c
      character*255 title_vecs                 ! Title
      character*255 v_out                 ! Output filename
      integer nset, nmo
c
c     externals
c
      logical movecs_write
      external movecs_write
c
c  Get filenames
c
      if (.not.rtdb_cget(rtdb,'mcscf:output vectors',1,v_out)) then
         call errquit('mcscf: failed to get name for output vectors',0,
     &       RTDB_ERR)
      endif
c
c
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title_vecs))
     $     title_vecs = ' '
c
c
c
      nset = 1
      nmo = nbf
c
c
c
      if (.not. movecs_write(rtdb, basis, v_out, 'mcscf',
     $     title_vecs, nbf, nset, nmo, occ, nbf, evals, nbf, g_movecs ))
     $   call errquit('mcscf_movecs_write: movecs_write failed', 0,
     &       DISK_ERR)
      return
      end










c
c  $Id$
c
c
c Return the gradient vector from the generalized
c Fock matrix.
c
c                gen     gen
c    g   =  2 ( F    -  F   )                     (1)
c     pq         pq      qp
c
c
c            gen   
c Note that F    = 0      (a=virtual)             (2)
c            ap    
c
c So only the inactive and active rows of the generalized
c Fock matrix are non-zero, so just copy and scale elements, viz,
c
c             gen                       gen
c    g   = 2 F                g   =  2 F
c     ia      ia               ta       ta
c
c But still have to compute inactive-active part of
c the gradient using (1)
c
c
c
       subroutine mcscf_gfock2grad( nbf, nclosed, nact,
     $                              g_gfock, g_grad )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
c
       integer nbf, nclosed, nact               ! [input] Occupation parameters
       integer g_gfock                          ! [input] MCSCF Generalized Fock matrix
       integer g_grad                           ! [output] Gradient
c
c
c
       integer g_tmp
       integer nvir, voff, aoff, aend, xoff, xend
c
c
c
       nvir = nbf - nclosed - nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
c         g   =  F            g  = F
c          ia     ia           ta   ta
c
       call ga_zero(g_grad)
       xend = (nclosed+nact)*nvir
       call ga_copy_patch( 't', g_gfock, 1, aend, voff, nbf,
     $                          g_grad, 1, xend, 1, 1 )
c
c         g   =  F   -  F
c          ti     it     ti
c
*ga:1:0
c
       if (.not.ga_create(MT_DBL,nact,max(1,nclosed),'crap',nact,
     $                    0,g_tmp))
     $      call errquit('mcscf_gfock2grad: cannot allocate',0, GA_ERR)
       call ga_copy_patch('t', g_gfock, 1, nclosed, aoff, aend,
     $                         g_tmp, 1, nact, 1, nclosed )
       call ga_dadd_patch(-1.d0, g_gfock, aoff, aend, 1, nclosed,
     $                     1.d0, g_tmp, 1, nact, 1, nclosed,
     $                           g_tmp, 1, nact, 1, nclosed )
       xoff = xend + 1
       xend = xend + (nact*nclosed)
       call ga_copy_patch( 'n', g_tmp, 1, nact, 1, nclosed,
     $                          g_grad, xoff, xend, 1, 1 )
c
c         g   =  2.g
c          pq       pq
c
       call ga_dscal( g_grad, 2.d0 )
       if (.not.ga_destroy(g_tmp))
     $      call errquit('mcscf_gfock2grad: cannot destroy',0, GA_ERR)
       return
       end








c
c
c   Generate the High-spin ROHF 1-PDM & 2-PDM for 
c   partially occupied orbitals (trivially)
c   from occupation numbers
c   All other active orbital densities are zero
c
c   For debugging purposes
c
c
       subroutine mcscf_rohf_modens( nact, ndbl, nopen, dm1, dm2 )
       implicit none
       integer nact
       integer ndbl
       integer nopen
       double precision dm1(nact,nact)
       double precision dm2(nact,nact,nact,nact)
c
       integer t,u,tt,uu,i,j
c
       t = nact*nact
       call dfill(t,0.d0,dm1,1)
       t = t*t
       call dfill(t,0.d0,dm2,1)
       do t=1,nopen
         tt = ndbl + t
         dm1(tt,tt) = 1.d0
         do u=1,nopen
           if (u.ne.t) then
             uu = ndbl + u
             dm2(tt,tt,uu,uu) = 0.5d0
C             dm2(tt,uu,uu,tt) = 0.5d0
           endif
         enddo
       enddo
       do i=1,ndbl
         dm1(i,i) = 2.d0
         dm2(i,i,i,i) = 1.d0
         do j=1,i-1
           dm2(i,i,j,j) = 2.d0
           dm2(j,j,i,i) = 2.d0
         enddo
       enddo
       do i=1,ndbl
         do t=1,nopen
           tt = ndbl + t
           dm2(i,i,tt,tt) = 1.d0
           dm2(tt,tt,i,i) = 1.d0
         enddo
       enddo
       do i=1,ndbl
         do j=1,ndbl
           if (i.ne.j) dm2(i,j,j,i) = -1.d0
         enddo
       enddo
       do i=1,ndbl
         do t=1,nopen
           tt = ndbl + t
           dm2(i,tt,tt,i) = -0.5d0
           dm2(tt,i,i,tt) = -0.5d0
         enddo
       enddo
c
c
       return
       end



c
c  Approximate the number of *open* shell
c  electrons from diagonal of 1-PDM (totally bogus)
c       
       integer function mcscf_rohf_den2occ( nact, dm1 )
       implicit none
       integer nact
       double precision dm1(nact,nact)
       integer i
       double precision nopen
c
       nopen = 0.d0
       do i=1,nact
         nopen = nopen + dm1(i,i)
       enddo
       mcscf_rohf_den2occ = nint(nopen)
       return
       end








      logical function check_moblocked(g_a,n,jlo,jhi)
      implicit none
#include "global.fh"
#include "msgids.fh"
      integer g_a
      integer n
      integer jlo,jhi
c
      integer gtype, dim1, dim2
      integer ilo, ihi, myid, jstat
      logical status
c
      myid = ga_nodeid()
      call ga_inquire(g_a,gtype,dim1,dim2)
      status = (dim1.eq.n).and.(dim2.eq.n)
      call ga_distribution(g_a,myid,ilo,ihi,jlo,jhi)
      if ((jlo.ne.0).and.(jhi.ne.-1)) then
        status = status.and.((ilo.eq.1).and.(ihi.eq.n))
      endif
      call ga_sync()
      jstat = 0
      if (status) jstat = 1
! igop('*') detects the presence of one or more jstat=0 inputs
#if NWCHEM_USE_IGOP_PROD
      call ga_igop(Msg_ChkMO,jstat,1,'*')
#else
      call ga_igop(Msg_ChkMO,jstat,1,'min')
#endif
      check_moblocked = (jstat.ne.0)
      return
      end

      
      




c
c
c Fill a square matrix with a scalar diagonal patch
c
       subroutine diagfill_patch( g_a, x, ilo, ihi )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer g_a
       integer ilo, ihi
       double precision x
c
       integer ilen, nelt, i, ii
       integer l_p, k_p
c
       ilen = ihi - ilo + 1
       nelt = ilen*ilen
       if (.not.ma_push_get(MT_DBL,nelt,'delta',l_p,k_p))
     $      call errquit('diagfill_patch: cannot allocate',0, MA_ERR)
       call dfill(nelt,0.d0,dbl_mb(k_p),1)
       do i=1,ilen
         ii = (i-1)*ilen + i - 1
         dbl_mb(k_p+ii) = x
       enddo
       if (ga_nodeid().eq.0)
     $    call ga_put(g_a,ilo,ihi,ilo,ihi,dbl_mb(k_p),ilen)
       if (.not.ma_pop_stack(l_p))
     $      call errquit('diagfill_patch: pop stack failed',0, MA_ERR)
       
       return
       end







c
c   Check and zero symmetry forbidden rotations
c
      subroutine mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                            g_v, oscreen, oprint, maxelt )
      implicit none
#include "global.fh"
      integer nbf, nclosed, nact
      integer orbsym(nbf)
      integer g_v
      logical oprint, oscreen
      double precision maxelt
c
      integer nvir, vlen, voff, aoff
      integer ii, i, a, t, aa, tt
      double precision xx, zero
c
      call ga_sync()
c
      nvir = nbf - nclosed - nact
      vlen = (nact+nclosed)*nvir + nact*nclosed
      aoff = nclosed*nvir
      voff = (nact+nclosed)*nvir
      maxelt = 0.d0
      zero = 0.d0
      do i=1+ga_nodeid(),nclosed, ga_nnodes()
        do a=1,nvir
          aa = nclosed + nact + a
          if (orbsym(i).ne.orbsym(aa)) then
            ii = (i-1)*nvir + a
            call ga_get(g_v,ii,ii,1,1,xx,1)
            maxelt = max(maxelt,abs(xx))
            if (oprint) write(6,881) ii,i,aa,abs(xx)
 881        format(5x,i5,5x,2i3,5x,e20.8)
            if (oscreen) call ga_put(g_v,ii,ii,1,1,zero,1)
          endif
        enddo
      enddo
      do t=1+ga_nodeid(),nact, ga_nnodes()
        tt = nclosed + t
        do a=1,nvir
          aa = nclosed + nact + a
          if (orbsym(tt).ne.orbsym(aa)) then
            ii = aoff + (t-1)*nvir + a
            call ga_get(g_v,ii,ii,1,1,xx,1)
            maxelt = max(maxelt,abs(xx))
            if (oprint) write(6,881) ii,tt,aa,abs(xx)
            if (oscreen) call ga_put(g_v,ii,ii,1,1,zero,1)
          endif
        enddo
      enddo
      do i=1+ga_nodeid(),nclosed, ga_nnodes()
        do t=1,nact
          tt = nclosed + t
          if (orbsym(i).ne.orbsym(tt)) then
            ii = voff + (i-1)*nact + t
            call ga_get(g_v,ii,ii,1,1,xx,1)
            maxelt = max(maxelt,abs(xx))
            if (oprint) write(6,881) ii,i,tt,abs(xx)
            if (oscreen) call ga_put(g_v,ii,ii,1,1,zero,1)
          endif
        enddo
      enddo
c
      call ga_sync()
c
      end





c
c  Convert active orbitals to natural orbitals
c  by diagonalizing 1-pdm
c
      subroutine mcscf_naturalorb( nbf, nclosed, nact, dm1, occ,
     $                             evals, g_movecs )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nbf
      integer nclosed
      integer nact
      double precision dm1(nact,nact)
      double precision occ(nbf)
      double precision evals(nbf)
      integer g_movecs
c
      integer l_tmp, k_tmp
      integer tsize, info, i, ii, jj
      integer l_u, k_u, l_x, k_x, l_y, k_y
      integer idamax
      external idamax
c
c  Orbital occupation
c
      call dfill(nbf, 0.0d0, occ, 1)
      do i=1,nclosed                    
        occ(i) = 2.d0                     ! Trivial inactive occupation
      enddo
c
c
c
      tsize = 4*nact
      if (.not.ma_push_get(MT_DBL, (nact*nact), 'u', l_u, k_u))
     $  call errquit('mcscf_naturalorb: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL, (tsize+nact), 'tmp', l_tmp, k_tmp))
     $  call errquit('mcscf_naturalorb: cannot allocate',0, MA_ERR)
      call dsyev( 'V', 'L', nact, dm1, nact, dbl_mb(k_tmp),
     $             dbl_mb(k_tmp+nact), tsize, info )
      if (info.ne.0) then
        if (info.lt.0) then 
          call errquit("mcscf_naturalorb: dsyev has an illegal value "
     $               //"for argument no.",-info,UERR)
        else
          call errquit("mcscf_naturalorb: dsyev failed to converge",
     $                 info,UERR)
        endif
      endif
c
c     Dsyev produced occupation numbers in increasing order, 
c     chemically having them in decreasing order makes more
c     sense. So reorder the vectors and occupation numbers.
c
      do i = 1, nact
        occ(nclosed+nact-i+1) = dbl_mb(k_tmp+i-1)
        jj = (nact-i)*nact
        call dcopy(nact,dm1(1,i),1,dbl_mb(k_u+jj),1)
      enddo
c$$      write(6,'(A)') ' Natural orbital rotation matrix '
c$$      call moints_matprint(nact,nact,dbl_mb(k_u))
      if (.not.ma_pop_stack(l_tmp))
     $   call errquit('mcscf_naturalorb: cannot pop stack',0, MA_ERR)
c
c
c
      if (ga_nodeid().eq.0) then
        if (.not.ma_push_get(MT_DBL, (nbf*nact), 'tmpx', l_x, k_x))
     $       call errquit('mcscf_naturalorb: cannot allocate',0, MA_ERR)
        if (.not.ma_push_get(MT_DBL, (nbf*nact), 'tmpy', l_y, k_y))
     $       call errquit('mcscf_naturalorb: cannot allocate',0, MA_ERR)
        
        call ga_get(g_movecs, 1, nbf, nclosed+1, nclosed+nact,
     $              dbl_mb(k_x), nbf )
        call dgemm( 'n', 'n', nbf, nact, nact, 1.d0, dbl_mb(k_x), nbf,
     $              dbl_mb(k_u), nact, 0.d0, dbl_mb(k_y), nbf )
        call ga_put(g_movecs, 1, nbf, nclosed+1, nclosed+nact,
     $              dbl_mb(k_y), nbf )

        if (.not.ma_pop_stack(l_y))
     $    call errquit('mcscf_naturalorb: cannot pop stack',0, MA_ERR)
        if (.not.ma_pop_stack(l_x))
     $    call errquit('mcscf_naturalorb: cannot pop stack',0, MA_ERR)
      endif
c
c
c
      if (.not.ma_pop_stack(l_u))
     $  call errquit('mcscf_naturalorb: cannot pop stack',0, MA_ERR)
      return
      end

      





c
c  Sort integer array
c
      subroutine util_isort( n, iv )
      implicit none
      integer n
      integer iv(n)
      integer i, j, itmp

      do i=1,n
        do j=1,n-1
          if (iv(j).gt.iv(j+1)) then
            itmp    = iv(j)
            iv(j)   = iv(j+1)
            iv(j+1) = itmp
          endif
        enddo
      enddo
      return
      end









c
c   Individually canonicalize inactive, active and virtual
c   MOs. Use generalized Fock for inactive and active
c   parts. Virtual block is defined as
c
c                  I       A
c          F   =  F    +  F
c           ab     ab      ab
c
c   Maintain order of previous MO columns.
c
c
      subroutine mcscf_canonical( nbf, nclosed, nact, g_ifock, g_afock,
     $                            g_gfock, eps, g_movecs, nflip )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
      integer nbf, nclosed, nact, nclos
      integer g_ifock, g_afock, g_gfock
      double precision eps(*)
      integer g_movecs
      integer nflip
c
      integer g_tmp, g_u, g_v, g_z
      integer i, aoff, aend, nvir
      integer l_cc, k_cc
c
      integer idamax
      external idamax
c
c  Inactive orbitals
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'rotation',nbf,0,g_v))
     $  call errquit('mcscf_canonical: cannot allocate v',0, GA_ERR)
      call ga_zero(g_v)
*ga:1:0
      nclos = max(1,nclosed)
      if (.not.ga_create(MT_DBL,nclos,nclos,'tmp',nclos,0,g_tmp))
     $  call errquit('mcscf_canonical: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nclos,nclos,'u',nclos,0,g_u))
     $  call errquit('mcscf_canonical: cannot allocate u',0, GA_ERR)
      call ga_copy_patch( 'n', g_gfock, 1, nclosed, 1, nclosed,
     $                         g_tmp, 1, nclosed, 1, nclosed )
      call ga_symmetrize( g_tmp )
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
         call ga_pdsyev  (g_tmp, g_u, eps, 0)
#else
      call ga_diag_std( g_tmp, g_u, eps )
#endif
#else
      call ga_diag_std_seq( g_tmp, g_u, eps )
#endif
      call ga_copy_patch( 'n', g_u, 1, nclosed, 1, nclosed,
     $                         g_v, 1, nclosed, 1, nclosed )
      if (.not.ga_destroy(g_u))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
c
c  Active orbitals
c
      aoff = nclosed + 1
      aend = nclosed + nact
*ga:1:0
      if (.not.ga_create(MT_DBL,nact,nact,'tmp',nact,0,g_tmp))
     $  call errquit('mcscf_canonical: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nact,nact,'u',nact,0,g_u))
     $  call errquit('mcscf_canonical: cannot allocate u',0, GA_ERR)
      call ga_copy_patch( 'n', g_gfock, aoff, aend, aoff, aend,
     $                         g_tmp, 1, nact, 1, nact )
      call ga_symmetrize( g_tmp )
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
         call ga_pdsyev  (g_tmp, g_u, eps(aoff), 0)
#else
      call ga_diag_std( g_tmp, g_u, eps(aoff) )
#endif
#else
      call ga_diag_std_seq( g_tmp, g_u, eps(aoff) )
#endif
      call ga_copy_patch( 'n', g_u, 1, nact, 1, nact,
     $                         g_v, aoff, aend, aoff, aend )
      if (.not.ga_destroy(g_u))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
c
c  Virtual orbitals
c
      aoff = nclosed + nact + 1
      nvir = nbf - nclosed - nact
*ga:1:0
      if (.not.ga_create(MT_DBL,nvir,nvir,'tmp',nvir,0,g_tmp))
     $  call errquit('mcscf_canonical: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nvir,nvir,'u',nvir,0,g_u))
     $  call errquit('mcscf_canonical: cannot allocate u',0, GA_ERR)
      call ga_dadd_patch( 2.d0, g_ifock, aoff, nbf, aoff, nbf,
     $                    2.d0, g_afock, aoff, nbf, aoff, nbf,
     $                          g_tmp, 1, nvir, 1, nvir )
      call ga_symmetrize( g_tmp )
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
         call ga_pdsyev  (g_tmp, g_u, eps(aoff), 0)
#else
      call ga_diag_std( g_tmp, g_u, eps(aoff) )
#endif
#else
      call ga_diag_std_seq( g_tmp, g_u, eps(aoff) )
#endif
      call ga_copy_patch( 'n', g_u, 1, nvir, 1, nvir,
     $                         g_v, aoff, nbf, aoff, nbf )
      if (.not.ga_destroy(g_u))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('mcscf_canonical: cannot destroy tmp',0, GA_ERR)
      call dscal( nbf, 0.5d0, eps, 1 )
c
c  Transform MOs
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'rotation',nbf,0,g_z))
     $  call errquit('mcscf_canonical: cannot allocate z',0, GA_ERR)
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_movecs,
     $                g_v, 0.d0, g_z )
      call ga_copy( g_z, g_movecs )
c
c  Probably no need for this
c  since these arrays are soon destroyed
c
      call two_index_transf( g_ifock, g_v, g_v, g_z, g_ifock )
      call two_index_transf( g_afock, g_v, g_v, g_z, g_afock )
      call two_index_transf( g_gfock, g_v, g_v, g_z, g_gfock )
c
c  Check rotation matrix for flipped MOs
c
      if (.not.ma_push_get(MT_DBL, nbf,'col tmp', l_cc, k_cc))
     $   call errquit('mcscf_canonical: cannot allocate cc',0, MA_ERR)
      nflip = 0
      do i=ga_nodeid()+1,nbf,ga_nnodes()
        call ga_get( g_v, 1, nbf, i, i, dbl_mb(k_cc), nbf )
        if (idamax( nbf, dbl_mb(k_cc), 1 ).ne.i) then
          nflip = nflip + 1
        endif
      enddo
      call ga_igop(msg_mcscf_sumflip, nflip, 1, '+')
      if (.not.ma_pop_stack(l_cc))
     $  call errquit('mcscf_canonical: cannot pop stack',0, GA_ERR)
c
c
c
      if (.not.ga_destroy(g_z))
     $  call errquit('mcscf_canonical: cannot destroy v',0, GA_ERR)
      if (.not.ga_destroy(g_v))
     $  call errquit('mcscf_canonical: cannot destroy v',0, GA_ERR)
      return
      end









      subroutine mcscf_save_vec( rtdb, geom, basis, nbf, orbe, occ,
     $                           g_movecs )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb
      integer geom
      integer basis
      integer nbf
      double precision orbe(*), occ(*)
      integer g_movecs
c
c
c
      logical olock
      integer g_tmp
      integer k_oldocc, l_oldocc, k_oldev, l_oldev
      integer i
c
c
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'tmp',nbf,0,g_tmp))
     $   call errquit('mcscf_save_vecs: cannot allocate tmp',0, GA_ERR)
      if (.not.rtdb_get(rtdb, 'mcscf:movecs lock', MT_LOG,
     $     1, olock)) olock = .false.
      if (.not.ma_push_get(MT_DBL, nbf,'tmp evals', l_oldev, k_oldev))
     $  call errquit('mcscf_save_vecs: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL, nbf,'tmp occ', l_oldocc, k_oldocc))
     $  call errquit('mcscf_save_vecs: cannot allocate',0, MA_ERR)
c
c  Read in old vectors and occupations
c  Optionally lock movecs ordering with old vectors
c
      if (olock) then
        call mcscf_movecs_read( rtdb, basis, nbf, dbl_mb(k_oldocc),
     $                          dbl_mb(k_oldev), g_tmp )
        call dcopy(nbf, orbe, 1, dbl_mb(k_oldev), 1 )
        call dcopy(nbf, occ,  1, dbl_mb(k_oldocc), 1)
        call movecs_lock(geom, basis, g_tmp, g_movecs, 
     $       dbl_mb(k_oldev), dbl_mb(k_oldocc))
      else
        call dcopy(nbf, orbe, 1, dbl_mb(k_oldev), 1 )
        call dcopy(nbf, occ,  1, dbl_mb(k_oldocc), 1)
      endif
c
c  Info print
c
      if (util_print('save vectors',print_high).and.
     $   (ga_nodeid().eq.0)) then
        write(6,911) 
 911    format('Saved MO vectors  - orbital energies')
        write(6,881) (dbl_mb(k_oldev+i-1),i=1,nbf)
        write(6,912)
 912    format('Saved MO vectors  - orbital occupations')
        write(6,881) (dbl_mb(k_oldocc+i-1),i=1,nbf)
 881    format(7f12.6)
      endif
c
c  Write out MOs
c
      call mcscf_movecs_write( rtdb, basis, nbf, dbl_mb(k_oldocc),
     $                         dbl_mb(k_oldev), g_movecs )
c
c  Clean up
c
      if (.not.ma_pop_stack(l_oldocc))
     $  call errquit('mcscf_save_vecs: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_oldev))
     $  call errquit('mcscf_save_vecs: cannot pop stack',0, MA_ERR)
      if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_save_vecs: cannot destroy',0, GA_ERR)
      return
      end







      subroutine sym_matrix_check( n, irreps, g_a, isym, tol, nmixed )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"      
#include "msgids.fh"
      integer n
      integer irreps(n)
      integer g_a
      integer isym
      double precision tol
      integer nmixed
c
      integer l_t, k_t
      integer i, j
c
      integer mulx
#include "bitops.fh"
      mulx(i,j) = ieor((i-1),(j-1)) + 1
c
      if (.not.ma_push_get(MT_DBL, n, 'tmp', l_t, k_t))
     $  call errquit('sym_matrix_check: cannot allocate',0, MA_ERR)
      nmixed = 0
      do i=ga_nodeid()+1,n,ga_nnodes()
        call ga_get(g_a, 1, n, i, i, dbl_mb(k_t), n)
        do j=1,n
          if ((mulx(irreps(i),irreps(j)).ne.isym).and.
     $        (abs(dbl_mb(k_t+j-1)).gt.tol)) then
            nmixed = nmixed + 1
          endif
        enddo
      enddo
      call ga_igop(msg_mcscf_symcheck,nmixed,1,'+')

      if (.not.ma_pop_stack(l_t))
     $  call errquit('sym_matrix_check: failed pop stack',0, MA_ERR)
      return
      end







c
c
c   Construct 2x2 unitrary rotation
c   of phi degrees between columns ii & jj
c
      subroutine theta2u( nbf, ii, jj, phi, g_u )
      implicit none
#include "global.fh"
      integer nbf
      integer ii, jj
      double precision phi
      integer g_u
      double precision xx
      integer i

      call ga_zero(g_u)
      do i=ga_nodeid()+1,nbf,ga_nnodes()
        call ga_put(g_u,i,i,i,i,1.d0,1)
      enddo
      call ga_sync()
      if (ga_nodeid().eq.0) then
        xx = cos(phi)
        call ga_put(g_u,ii,ii,ii,ii,xx,1)
        call ga_put(g_u,jj,jj,jj,jj,xx,1)
        xx = sin(phi)
        call ga_put(g_u,ii,ii,jj,jj,xx,1)
        call ga_put(g_u,jj,jj,ii,ii,-xx,1)
      endif
      call ga_sync()
      return
      end

