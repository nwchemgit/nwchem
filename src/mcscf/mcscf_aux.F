c
c  Convert occupation vector to 
c  number of closed and active orbitals
c
      subroutine mcscf_occ2int( nbf, occ, nclosed, nact )
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
      do i=i,nbf
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
c
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
c
c
c
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
c
c
c
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
       if (.not.ga_create(MT_DBL,nact,nact,'1PDM',nact,1,g_dm1))
     $   call errquit('mcscf_occ2dens: cannot create 1PDM',0)
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
     $   call errquit('mcscf_occ2dens: cannot destroy',0)
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
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "bas.fh"
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
      integer i, nbf_vecs, nsets, nmo_vecs(2)
      integer l_oph, k_oph
      integer l_v, k_v
      integer j, j1, j2, nn
      logical oskel, oprint, status, ovecs
      double precision tol2e
      double precision vphase
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
c  Get misc info
c
      if (.not. rtdb_get(rtdb, 'scf:tol2e',
     $     MT_DBL, 1, tol2e)) tol2e = 1.0d-7
      if (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
      if (.not. bas_name(basis, name, trans))
     $     call errquit('scf_init: bas_name?', 0)
c
c  Get filenames
c
       if (.not.rtdb_cget(rtdb,'mcscf:input vectors',1,movecs_in)) then
         if (.not.rtdb_cget(rtdb,'scf:input vectors',1,movecs_in))
     $      movecs_in = 'movecs'
       endif
c
c  Starting MO's determined from movecs_in string
c
      call ga_zero(g_movecs)
      if (movecs_in .eq. 'hcore' .or. movecs_in .eq. 'atomic') then
c
c  Hcore or atomic guess
c
         g_aodens = ga_create_atom_blocked(geom, basis, 'Guess density')
         call ga_zero(g_aodens)
         if (movecs_in .eq. 'atomic') then
            call guess_dens(geom, basis, g_aodens)
         endif
         call rhf_dens_to_mo(rtdb, geom, basis, nelec, nbf,
     $                       tol2e, g_aodens, g_movecs, oskel)
         if (.not. ga_destroy(g_aodens)) call errquit
     $        ('scf_vectors_guess: destroy of density failed',0)
      else if (movecs_in .eq. 'project') then
c
c  Project from a smaller basis set
c
         call movecs_project(rtdb, geom, basis, g_movecs)
      else
c
c  Read vectors from file movecs_in and see if 
c  they make sense
c
         status = movecs_read_header( movecs_in, title_vecs, bas_vecs,
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
     $            call errquit('mcscf_getmovecs: movecs read',0)
            else
               if (ga_nodeid().eq.0) write(6,944)
 944           format('Restart movecs mismatch with input')
               status = .false.
            endif
            if (oprint) then
               write(6,77) title_vecs(1:inp_strlen(title_vecs))
 77            format(/' Loading old vectors from job with title :'//a/)
               call util_flush(6)
            endif
         endif
         call movecs_swap(rtdb,'rohf',g_movecs,occ,evals)
         if (.not. status) call errquit('rhf: failed to load vectors',0)
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
       inquire(file='vecs',exist=ovecs)
       if (ovecs) then
         if (util_print('information',print_default)) then
           if (ga_nodeid().eq.0) write(6,670)
 670       format(' Importing MO vectors from external file')
         endif
         if (.not.ma_push_get(MT_DBL, nbf,'ophase', l_v, k_v))
     $     call errquit('mcscf_getmovecs: cannot allocate',0)
         open(unit=19,file='vecs',status='old',form='formatted')
         read(19,*) nn
         do i=1,nbf
           do j1=1,nbf,6
             j2 = min(nbf,(j1+5))
             read(19,913) (dbl_mb(k_v+j-1),j=j1,j2)
 913         format(6e20.11)
           enddo
           call ga_put(g_movecs,1,nbf,i,i,dbl_mb(k_v),nbf)
         enddo
         if (.not.ma_pop_stack(l_v))
     $     call errquit('mcscf_getmovecs: failed pop stack',0)
       endif
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
      return
      end






c
c  $Id: mcscf_aux.F,v 1.4 1996-02-13 23:00:29 rg240 Exp $
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
       if (.not.ga_create(MT_DBL,nact,nclosed,'crap',nact,1,g_tmp))
     $      call errquit('mcscf_gfock2grad: cannot allocate',0)
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
     $      call errquit('mcscf_gfock2grad: cannot destroy',0)
       return
       end








c
c  $Id: mcscf_aux.F,v 1.4 1996-02-13 23:00:29 rg240 Exp $
c
c Generate the High-spin ROHF 1-PDM & 2-PDM for 
c partially occupied orbitals (trivially)
c from occupation numbers
c All other active orbital densities are zero
c
c For debugging purposes
c
c
       subroutine mcscf_rohf_modens( nopen, nact, dm1, dm2 )
       implicit none
       integer nopen
       integer nact
       double precision dm1(nact,nact)
       double precision dm2(nact,nact,nact,nact)
c
c
       integer t,u
c
c
c
       t = nact*nact
       call dfill(t,0.d0,dm1,1)
       t = t*t
       call dfill(t,0.d0,dm2,1)
       do t=1,nopen
         dm1(t,t) = 1.d0
         do u=1,t-1
           dm2(t,t,u,u) = 0.5d0
           dm2(u,u,t,t) = 0.5d0
           dm2(t,u,t,u) = -0.5d0
           dm2(u,t,u,t) = -0.5d0
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
      call ga_igop(Msg_ChkMO,jstat,1,'*')
      check_moblocked = (jstat.ne.0)
      return
      end

      
      




c
c
c Fill a square matrix with a scalar diagonal patch
c
c
c
       subroutine diagfill_patch( g_a, x, ilo, ihi )
       implicit none
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
     $      call errquit('diagfill_patch: cannot allocate',0)
       call dfill(nelt,0.d0,dbl_mb(k_p),1)
       do i=1,ilen
         ii = (i-1)*ilen + i - 1
         dbl_mb(k_p+ii) = x
       enddo
       if (ga_nodeid().eq.0)
     $    call ga_put(g_a,ilo,ihi,ilo,ihi,dbl_mb(k_p),ilen)
       if (.not.ma_pop_stack(l_p))
     $      call errquit('diagfill_patch: pop stack failed',0)
       
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
      nvir = nbf - nclosed - nact
      vlen = (nact+nclosed)*nvir + nact*nclosed
      aoff = nclosed*nvir
      voff = (nact+nclosed)*nvir
      maxelt = 0.d0
      zero = 0.d0
      do i=1,nclosed
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
      do t=1,nact
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
      do i=1,nclosed
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
      return
      end





c
c  Convert active orbitals to natural orbitals
c  by diagonalizing 1-pdm
c
      subroutine mcscf_naturalorb( nbf, nclosed, nact, dm1, occ,
     $                             g_movecs )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf
      integer nclosed
      integer nact
      double precision dm1(nact,nact)
      double precision occ(nbf)
      integer g_movecs
c
      integer l_tmp, k_tmp
      integer g_v
      integer tsize, info, i, j, ii, jj
c
c
c
      tsize = 4*nact
      if (.not.ma_push_get(MT_DBL, (tsize+nact), 'tmp', l_tmp, k_tmp))
     $  call errquit('mcscf_naturalorb: cannot allocate',0)
      call dsyev( 'V', 'L', nact, dm1, nact, dbl_mb(k_tmp),
     $             dbl_mb(k_tmp+nact), tsize, info )
      do i=1,nclosed                    
        occ(i) = 2.d0                     ! Trivial inactive occupation
      enddo
      do i=1,nact
        occ(nclosed+i) = dbl_mb(k_tmp+nact-i)
      enddo
      if (.not.ma_pop_stack(l_tmp))
     $   call errquit('mcscf_naturalorb: cannot pop stack',0)
c
c
c
      if (.not.ga_create(MT_DBL,nbf,nact,'gtmp',nbf,1,g_v))
     $   call errquit('mcscf_naturalorb: cannot allocate',0)
      call ga_zero(g_v)
      do i=nact,1,-1
        ii = nact - i + 1
        do j=1,nact
          jj = nclosed + j
          call ga_dadd_patch( dm1(j,i), g_movecs, 1, nbf, jj, jj,
     $                        1.d0,     g_v,      1, nbf, ii, ii,
     $                                  g_v,      1, nbf, ii, ii )
        enddo
      enddo
      call ga_copy_patch( 'n', g_v, 1, nbf, 1, nact,
     $             g_movecs, 1, nbf, (nclosed+1), (nclosed+nact) )
      if (.not.ga_destroy(g_v))
     $   call errquit('mcscf_naturalorb: cannot destroy',0)
      return
      end

      




