      subroutine hv_write(lun,ix,len,x)
*
* $Id$
*
      implicit none
      integer lun, ix, len
      double precision x(len)
      integer i, k, khi
      
      write(lun,901) ix,(x(i),i=1,8)
 901  format(i5,8f12.6)
      do k=9,len,8
        khi = min((k+7),len)
        write(lun,902) (x(i),i=k,khi)
 902    format(5x,8f12.6)
      enddo
      return
      end
      



      subroutine hv_reorder( nbf, nclosed, nact, g_x, x, y )
      implicit none
      integer nbf, nclosed, nact
      integer g_x
      double precision x(*), y(*)
c     
      integer nvir, vlen, voff, aoff, aend, xoff, xend, ii
      integer nca, i, j
c     
c     
      nvir = nbf - nclosed - nact
      vlen = (nclosed+nact)*nvir + nclosed*nact
      voff = nclosed + nact + 1
      aoff = nclosed + 1
      aend = nclosed + nact
      nca = nclosed + nact
c     
c     
      xoff = nvir*(nclosed+nact) + 1
      xend = vlen
      call ga_get(g_x,xoff,xend,1,1,x,(nclosed*nact))
      call dcopy((nclosed*nact),x,1,y,1)

      xoff = 1
      xend = nvir*nclosed
      call ga_get(g_x,xoff,xend,1,1,x,(nclosed*nvir))
      ii = nact*nclosed
      do i=1,nvir
        do j=1,nclosed
          y(ii+(i-1)*nca+j) = x((j-1)*nvir+i)
        enddo
      enddo

      xoff = xend + 1
      xend = nvir*(nclosed+nact)
      call ga_get(g_x,xoff,xend,1,1,x,(nact*nvir))
      do i=1,nvir
        do j=1,nact
          y(ii+(i-1)*nca+nclosed+j) = x((j-1)*nvir+i)
        enddo
      enddo

      call dscal(vlen,0.5d0,y,1)
c     
c     
      return
      end



c$$$
c$$$       subroutine mattr(n,m,x,y)
c$$$       implicit none
c$$$       integer n,m
c$$$       double precision x(n,m),y(m,n)
c$$$       integer i,j
c$$$
c$$$       do i=1,n
c$$$         do j=1,m
c$$$           y(j,i) = x(i,j)
c$$$         enddo
c$$$       enddo
c$$$       return
c$$$       end
c$$$




       subroutine hv_writev(lun,len,x)
       implicit none
       integer lun, len
       double precision x(len)
       
       write(lun) x
       return
       end



       subroutine hv_readv(lun,len,x)
       implicit none
       integer lun, len
       double precision x(len)
       
       read(lun) x
       return
       end




       subroutine ga_rowprint( text, g_a )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       character*(*) text
       integer g_a
       integer gtype,vlen,rlen
       integer l_tmp, k_tmp
       integer i
       
       if (ga_nodeid().eq.0) then
         call ga_inquire(g_a,gtype,vlen,rlen)
         if (.not.ma_push_get(MT_DBL,vlen,'tmp',l_tmp,k_tmp))
     $     call errquit('ga_rowprint: cannot allocate',0, MA_ERR)
         call ga_get(g_a,1,vlen,1,1,dbl_mb(k_tmp),1)
         write(6,909) text
 909     format(/,a)
         write(6,910) (dbl_mb(k_tmp+i-1),i=1,vlen)
 910     format(6F12.8)
         if (.not.ma_pop_stack(l_tmp))
     $     call errquit('ga_rowprint: cannot pop stack',0, MA_ERR)
       endif
       call ga_sync()
       return
       end





      subroutine mcscf_twopdm_print(n,dm2)
      implicit none
      integer n
      double precision dm2(n,n,n,n)
      double precision TOL
      integer i,j,k,l,ltop
      data TOL/2.d-1/
      
      do i=1,n
        do j=1,n
          do k=1,n
            ltop = k
            if (k.gt.j) ltop = j
            do l=1,n
              if (abs(dm2(i,j,k,l)).gt.TOL) then
                write(6,771) i,j,k,l,dm2(i,j,k,l)
 771            format(4i5,f18.10)
              endif
            enddo
          enddo
        enddo
      enddo
      return
      end





      subroutine ga_print_x( g_a )
      implicit none
#include "global.fh"
      integer g_a
      integer type, n, m
      integer dim
      parameter(dim=40)
      double precision xtmp(dim*dim)
      integer chi, clo, i, j
      
      call ga_inquire(g_a, type, n, m )
      if ((n*m).gt.(dim*dim)) return
      call ga_get(g_a,1,n,1,n,xtmp,n)
      chi = 0
 11   clo = chi + 1
      chi = min((clo + 7),m)
      do j=1,n
        write(6,771) (xtmp((i-1)*n + j),i=clo,chi)
 771    format(8f12.6)
      enddo
      write(6,*)
      if (chi.lt.n) goto 11
      return
      end



      subroutine mcscf_trprint(n,x)
      implicit none
      integer n
      double precision x(*)
      integer i, ii, j, itop

      do i=1,n
        ii = (i*(i-1))/2
        itop = min(10,i)
        write(6,922) (x(ii+j),j=1,itop)
 922    format(10f12.6)
      enddo
      return
      end




      







#ifdef MCSCF_DEBUGGER
c
c
c   This is the main routine for debugging
c   the MCSCF (this is not distributed)
c
c
      subroutine mcscf_debugger( rtdb, basis, geom, nbf, nclosed, nact,
     $                           g_movecs, evals )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "sym.fh"
#include "pstat.fh"
#include "mcscfprof.fh"
c
       integer rtdb
       integer geom, basis
       integer nbf, nclosed, nact
       integer g_movecs
       double precision evals(*)
c
       integer nvir, noper, nsym
       integer nactel, nela, nelb, multip, orlen
       integer l_occ, k_occ, l_sym, k_sym
       integer l_dm1, k_dm1, l_dm2, k_dm2
       integer l_tmp, k_tmp
       integer g_coul, g_exch
       integer g_grad, g_prod, g_x
       integer mo_lo, mo_hi
       integer i, j, nmixed, clo, chi
       integer blen
       double precision pfac
       double precision eone, etwo, energy, enrep, etrace
       double precision ecore, eci, e0
       double precision citol
       double precision tol2e, gnorm, xx
       logical oskel
       logical ohalf, oblk
c
c---------------------
c  Debugging variables
c
       integer mclosed, mopen                             ! ROHF occupation
       integer g_fcv, g_fpv, g_fcp                        ! ROHF Fock matrices
       integer g_afock, g_ifock, g_gfock
       integer g_u, g_b, g_newgfock, g_tmp, g_tmp2
       integer g_coul2, g_exch2, g_grad2
       integer info

       integer hdim
       parameter(hdim=100)
       double precision crap(1000), hh(hdim*hdim)
       double precision scr(4*hdim), ev(hdim)
       double precision edel4, edel, theta, gg1, gg2, edelF
       double precision phi, xxstep, etgt, ezzz

       double precision cjgtol
       integer iii
       double precision phimin,phimax,phiinc
C       data phimin,phimax,phiinc/-6.3d0,6.3d0,0.1d0/
       data phimin,phimax,phiinc/1.0d0,1.0d0,0.2d0/
      
c---------------------
c
       integer ga_create_atom_blocked, ga_create_JKblocked
       external ga_create_atom_blocked, ga_create_JKblocked
       integer mcscf_rohf_den2occ
       external mcscf_rohf_den2occ
       double precision ga_trace_diag
       external ga_trace_diag
c
       data ohalf/.true./
       data blen/16/
c
c
c
c  Get w.f. parameters
c
       nvir = nbf - nclosed - nact
       orlen = (nclosed*nvir) + (nact*nvir) + (nclosed*nact)
       nsym = sym_number_ops(geom)+1
       if (.not.rtdb_get(rtdb,'scf:skeleton',MT_LOG,1,oskel))
     $   oskel = sym_number_ops(geom).gt.0
       if (.not.geom_nuc_rep_energy( geom, enrep ))
     $   call errquit('mcscf: cannot retrieve nuclear repulsion',0,
     &       GEOM_ERR)
       if (.not.rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e))        
     $      tol2e = 1.d-12                                                   ! Redundant recovered later
c
c  Get electron and spin multiplicity (NB: for info only)
c  Active elec and multiplicity must be set --- no defaults
c
       if (.not.rtdb_get(rtdb,'mcscf:nactelec',MT_INT,1,nactel))
     $   call errquit('number of active electrons not set',0,
     &       RTDB_ERR)
       if (.not.rtdb_get(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $   call errquit('spin multiplicity not set',0, RTDB_ERR)
       nela = (nactel + multip - 1)/2
       nelb = nactel - nela
       if ((mod((nactel + multip - 1),2).ne.0).or.
     $    (nela.lt.0).or.(nelb.lt.0))
     $    call errquit('mcscf: incompatible elec and spin',0,
     &       INPUT_ERR)
       MOPEN = 2
       MCLOSED = (NACTEL - MOPEN)/2
c
c  Print info
c
       write(6,900)
 900   format(///,10x,40('='),
     $          /,17x,'MCSCF Debug Section',
     $          /,10x,40('='),//)
       write(6,901) nbf, nclosed, nact, nactel,
     $              multip, orlen
       write(6,902) (nclosed*nvir),(nact*nvir),
     $              (nclosed*nact)
 901   format(18x,35('-'),/,
     $        20x,'Basis functions:',10x,i5,/,
     $        20x,'Inactive shells:',10x,i5,/,
     $        20x,'Active shells:',12x,i5,/,
     $        20x,'Active electrons:',9x,i5,/,
     $        20x,'Multiplicity:',13x,i5,/,
     $        20x,'Orbital rotations:',8x,i5)
 902   format(25x,'Inact - Virt',9x,i5,/,
     $        25x,'Act - Virt',11x,i5,/,
     $        25x,'Inact - Act',10x,i5)
       write(6,903)
 903   format(18x,35('-'))
c
c  Allocate Fock and gradient matrices
c
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'Act Fock',nbf,0,g_afock))
     $      call errquit('mcscf: cannot allocate active Fock',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'In Fock',nbf,0,g_ifock))
     $      call errquit('mcscf: cannot allocate inactive Fock',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'Gen Fock',nbf,0,g_gfock))
     $      call errquit('mcscf: cannot allocate general Fock',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,orlen,1,'Gradient',0,0,g_grad))
     $      call errquit('rohf_head: cannot allocate',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,orlen,1,'Product',0,0,g_prod))
     $      call errquit('rohf_head: cannot allocate',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,orlen,1,'Arg vec',0,0,g_x))
     $      call errquit('rohf_head: cannot allocate',0, GA_ERR)
c
c  Create occupation vectors
c
       if (.not.ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ))
     $     call errquit('mcscf: cannot allocate',0, MA_ERR)
c
c  Allocate 1- & 2-PDM
c
       if (.not.ma_push_get(MT_DBL, (nact*nact*nact*nact),
     $                      '2P density', l_dm2, k_dm2))
     $     call errquit('mcscf: cannot allocate MO density',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL, (nact*nact),
     $                      '1P density', l_dm1, k_dm1))
     $     call errquit('mcscf: cannot allocate MO density',0, MA_ERR)
c
c  Get orbital symmetries
c 
       if (.not.ma_push_get(MT_INT, nbf, 'MO sym', l_sym, k_sym))
     $     call errquit('mcscf: cannot allocate symmetry',0, MA_ERR)
       call sym_movecs_adapt( basis, 1.d-8, g_movecs,
     $                        int_mb(k_sym), nmixed )
       if (nmixed .ne. 0) call errquit(
     $   'mcscf: symmetry contamination in starting MOs', nmixed,
     &       GEOM_ERR)
c
c  Print orbital info
c
       write(6,550)
 550   format(/,2x,'Starting Orbital Energies')
       write(6,551) (evals(i),i=1,nbf)
 551   format(7f12.6)
       write(6,887) 
 887   format(/,2x,'Orbital Symmetry Irreps')
       write(6,888) (int_mb(k_sym+i),i=0,nbf-1)
 888   format(16i3)
c
c   Allocate operator matrices
c   Memory test required here!
c
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       noper = (nact*(nact+1))/2
       g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
       g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
       if (.not.rtdb_get(rtdb,'mcscf:aoblock',MT_LOG,1,oblk)) then
         if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,1,oblk))
     $      oblk = .false.
       endif
c
c   Initial 4-Index Tranformation
c
       call moints_build_2x( basis, ohalf, oskel,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., blen, oblk )
c
c   Initial core energy
c
       call dfill((nact*nact),0.d0,dbl_mb(k_dm1),1)
       call dfill((nact*nact*nact*nact),0.d0,dbl_mb(k_dm2),1)
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    .false., oskel, tol2e, dbl_mb(k_dm1),
     $                    dbl_mb(k_dm2), g_movecs, g_coul,
     $                    eone, etwo, ecore )
c
c   
c
       cjgtol = 1.d-3
       gnorm = 4.d0
       e0 = ecore + enrep
c
c   Solve CI to set density
c
       citol = 1.d-8
       call mcscf_ifock( geom, basis, nbf, nclosed, nact,
     $                   oskel, tol2e, g_movecs, eone, etwo,
     $                   ecore, g_ifock )
       e0 = ecore + enrep
       call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                     nsym, int_mb(k_sym), e0,
     $                     evals, g_ifock, g_coul, 
     $                     citol, .false., .true., .true.,
     $                     dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c            Natural Orbital Section
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     
c$$$       call mcscf_naturalorb( nbf, nclosed, nact, dbl_mb(k_dm1),
c$$$     $                        dbl_mb(k_occ), evals, g_movecs )
c
c  Print natural orbitals
c
c$$$       write(6,880)
c$$$ 880   format(///,10x,'Natural Orbitals and Occupation')
c$$$       if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
c$$$     $                     l_tmp, k_tmp))
c$$$     $      call errquit('mcscf: cannot allocate local MO',0)
c$$$       call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
c$$$     $                       dbl_mb(k_tmp), nbf)
c$$$       chi = 0
c$$$ 33    clo = chi + 1
c$$$       chi = min((clo + 7),(nact+nclosed))
c$$$       write(6,*)
c$$$       write(6,881) (dbl_mb(k_occ+i-1),i=clo,chi)
c$$$ 881   format(8f12.6)
c$$$       write(6,*)
c$$$       do i=1,nbf
c$$$         write(6,881) (dbl_mb(k_tmp+(j-1)*nbf+i-1),j=clo,chi)
c$$$       enddo
c$$$       if (chi.ne.(nact+nclosed)) goto 33
c$$$       if (.not.ma_pop_stack(l_tmp))
c$$$     $      call errquit('mcscf: cannot pop local MO',0)
c
c  Resolve CI for Natural Orbitals
c
c$$$       write(6,228) 
c$$$ 228   format(/,'Resolving CI for natural orbitals',/)
c$$$       call moints_build_2x( basis, ohalf, oskel,
c$$$     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
c$$$     $                       g_movecs, g_coul, .true.,
c$$$     $                       g_exch, .true., blen, oblk )
c$$$       call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
c$$$     $                   g_coul, g_exch, g_ifock )
c$$$       call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
c$$$     $                     nsym, int_mb(k_sym), e0, evals,
c$$$     $                     g_ifock, g_coul, 
c$$$     $                     citol, .false., .false., .false.,
c$$$     $                     dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
c
c  Print 1-pdm and 2-pdm
c 
       if (util_print('density matrix',print_debug)) then
         if (ga_nodeid().eq.0) then
           write(6,671)
 671       format(/,'<<<<<<< 1pdm density matrix >>>>>>>>>')
           call moints_matprint( nact, nact, dbl_mb(k_dm1) )
           write(6,672)
 672       format(/,'<<<<<<< 2pdm density matrix >>>>>>>>>')
           call mcscf_twopdm_print(nact,dbl_mb(k_dm2))
           write(6,673)
 673       format(/,'<<<<<<< symm. 2pdm density matrix >>>>>>>>>')
           call mcscf_symmetrize_2pdm( nact, dbl_mb(k_dm2), crap )
           call mcscf_twopdm_print(nact,crap)
         endif
       endif
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              MCSCF Energy Trace and Fock Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c   Check energy trace
c     
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    .true., oskel, tol2e, dbl_mb(k_dm1),
     $                    dbl_mb(k_dm2), g_movecs, g_coul,
     $                    eone, etwo, ecore )
       etrace = eone + etwo + enrep
       write(6,674) etrace
 674   format(20x,'Trace Energy: ',f20.14)
c
c   Check Fock energy and gradient
c
       call ga_zero(g_grad)
       call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                  oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2),
     $                  g_movecs, g_coul, eone, etwo, e0,
     $                  g_ifock, g_afock, g_gfock )
       call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
       gnorm = sqrt(ga_ddot(g_grad,g_grad))
       energy = eone + etwo + enrep
       write(6,675) energy
 675   format(20x,'Fock Energy:  ',f20.14)
       write(6,676) gnorm
 676   format(20x,'Gradient norm:',e20.8)
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              MCSCF Hessian vector product Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  Regenerate MO integrals
c
       call moints_build_2x( basis, ohalf, oskel,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., blen, oblk )
c
c
       call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                   g_coul, g_exch, g_afock )
c
c  Finite difference gradient
c
c$$$       call ga_rowprint( ' ==== Analytical Gradient ====', g_grad )
c$$$       call mcscf_fdiff_grad( geom, basis, nbf, nclosed, nact,
c$$$     $                        oskel, tol2e, dbl_mb(k_dm1),
c$$$     $                        dbl_mb(k_dm2), g_movecs, g_coul, g_grad)
c$$$       call ga_rowprint( ' ==== Finite Diff Gradient ====', g_grad )
c$$$       gnorm = sqrt(ga_ddot(g_grad,g_grad))
c$$$       write(6,773) gnorm
c$$$ 773   format(/,10x,'Finite diff gradient norm:',e12.4,/)
c
c  MCSCF Fock matrices and gradient
c
       call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                  oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2),
     $                  g_movecs, g_coul, eone, etwo, e0,
     $                  g_ifock, g_afock, g_gfock )
       call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
c
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c                      B Matrix Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c   Make G matrix on disk 
c   (doesn't work..yet)
c
c$$$       call makeJK( basis, nbf, nclosed, nact, g_movecs )
c$$$       call mcscf_g( nbf, nclosed, nact, dbl_mb(k_dm1), dbl_mb(k_dm2),
c$$$     $               g_ifock, g_afock )
c
c
c
c$$$       if (.not.ga_duplicate(g_coul, g_coul2, 'Coul Copy'))
c$$$     $      call errquit('mcscf: cannot allocate Coulomb copy',0)
c$$$       if (.not.ga_duplicate(g_exch, g_exch2, 'Exch Copy'))
c$$$     $      call errquit('mcscf: cannot allocate Exch copy',0)
c$$$       if (.not.ga_duplicate(g_grad, g_grad2, 'Grad Copy'))
c$$$     $      call errquit('mcscf: cannot allocate Exch copy',0)
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'Unitary',nbf,0,g_u))
c$$$     $      call errquit('mcscf: cannot allocate unitrary',0)
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'Temp',nbf,0,g_tmp))
c$$$     $      call errquit('mcscf: cannot allocate tmp',0) 
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'Temp 2',nbf,0,g_tmp2))
c$$$     $      call errquit('mcscf: cannot allocate tmp',0) 
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'B',nbf,0,g_b))
c$$$     $      call errquit('mcscf: cannot allocate B',0)
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'new G Fock',nbf,0,g_newgfock))
c$$$     $      call errquit('mcscf: cannot allocate new GFock',0)
c
c  Step and generate U = exp(K) (2x2 rotation)
c
c$$$       call ga_zero(g_b)
c$$$       iii = (nclosed+nact-1)*nvir + 1
c$$$       do phi=phimin,phimax,phiinc
c$$$         xxstep = sqrt((1.d0 - cos(phi)))
c$$$         call ga_copy( g_coul, g_coul2 )
c$$$         call ga_copy( g_exch, g_exch2 )
c$$$         call ga_copy( g_grad, g_grad2 )
c$$$         call ga_copy( g_grad2, g_x)
c$$$         call ga_zero(g_u)
c$$$         call theta2u( nbf, nclosed+nact, nclosed+nact+1, phi, g_u )

         
c$$$         PRINT*,'Unitrary transformation'
c$$$         CALL GA_GET(G_U,1,NBF,1,NBF,CRAP,NBF)
c$$$         CALL MOINTS_MATPRINT(NBF,NBF,CRAP)
c
c  Generate B (Gradient response)
c
c$$$         call ga_zero(g_b)
c$$$         call mcscf_b( geom, basis, nbf, nclosed, nact, int_mb(k_sym),
c$$$     $                 dbl_mb(k_dm1), dbl_mb(k_dm2), oskel, tol2e,
c$$$     $                 g_movecs, g_ifock, g_afock, g_coul2, g_exch2,
c$$$     $                 g_u, g_b )
c
c  ~    t
c  A = U B
c
c$$$         call ga_dgemm( 't', 'n', nbf, nbf, nbf, 1.d0, g_b, g_u,
c$$$     $                  0.d0, g_newgfock )
c$$$         call mcscf_gfock2grad(nbf, nclosed, nact, g_newgfock,g_grad2)
C         CALL GA_ROWPRINT('Approximate New Grad',g_grad2)
c$$$         gnorm = sqrt(ga_ddot(g_grad2,g_grad2))
c$$$         call ga_get(g_grad2,iii,iii,1,1,gg2,1)
c
c            t
c  dE = tr( T (A + B) )
c
c$$$         call ga_zero(g_tmp)
c$$$         do i=1+ga_nodeid(),nbf,ga_nnodes()
c$$$           call ga_put(g_tmp,i,i,i,i,1.d0,1)
c$$$         enddo
c$$$         call ga_dadd( 1.d0, g_u, -1.d0, g_tmp, g_tmp )
c$$$         edelF = ga_ddot( g_tmp, g_gfock )
c$$$
c$$$         call ga_transpose( g_b, g_tmp2 )
c$$$         call ga_dadd( 1.d0, g_gfock, 1.d0, g_tmp2, g_newgfock )
c$$$         call ga_transpose( g_newgfock, g_tmp2 )
c$$$         edel4 = ga_ddot( g_tmp, g_tmp2 )
c
c  G operator code is broken!
c
c   (2)           t  ij
c  e    =  sum  (T  G   T )
c             ij           ij
c
c   (1)         t
c  e    =  tr( T A )
c
c$$$         call ga_get(g_tmp,1,nbf,1,nbf,crap,nbf)
c$$$         call mcscf_gt( nbf, nclosed, nact, crap, etgt )
c$$$         ezzz = 2.d0*ga_ddot( g_tmp, g_gfock ) + etgt
c
c  Compute exact new gradient for comparison
c
c$$$         call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_movecs, g_u,
c$$$     $                  0.d0, g_newgfock )
c$$$         call moints_build_2x( basis, ohalf, oskel,
c$$$     $                         mo_lo, mo_lo, mo_hi, 1, nbf,
c$$$     $                         g_newgfock, g_coul2, .true.,
c$$$     $                         g_exch2, .true., blen, oblk )
c$$$         call mcscf_fock( geom, basis, nbf, nclosed, nact,
c$$$     $                    oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2),
c$$$     $                    g_newgfock, g_coul2, eone, etwo, e0,
c$$$     $                    g_u, g_tmp, g_tmp2 )
c$$$         call mcscf_gfock2grad( nbf, nclosed, nact, g_tmp2, g_grad2 )
c$$$         CALL GA_ROWPRINT('Exact New Grad',g_grad2)
c$$$         gnorm = sqrt(ga_ddot(g_grad2,g_grad2))
c$$$         edel = eone + etwo + enrep - energy
c$$$         call ga_get(g_grad2,iii,iii,1,1,gg1,1)
c$$$         write(6,939) ((phi*180.d0)/3.1415629d0),
c$$$     $                 edel,edel4,ezzz,edelF,(edel4-edelF),etgt
c$$$ 939     format('qqqq',f10.2,6e15.6)
c$$$       enddo
c
c
c
c$$$       if (.not.ga_destroy(g_newgfock))
c$$$     $   call errquit('mcscf: cannot destroy unitrary',0)
c$$$       if (.not.ga_destroy(g_u))
c$$$     $   call errquit('mcscf: cannot destroy unitrary',0)
c$$$       if (.not.ga_destroy(g_tmp))
c$$$     $   call errquit('mcscf: cannot destroy unitrary',0)
c$$$       if (.not.ga_destroy(g_tmp2))
c$$$     $   call errquit('mcscf: cannot destroy unitrary',0)
c$$$       if (.not.ga_destroy(g_b))
c$$$     $   call errquit('mcscf: cannot destroy B',0)
c$$$       if (.not.ga_destroy(g_coul2))
c$$$     $   call errquit('mcscf: cannot destroy coul copy',0)
c$$$       if (.not.ga_destroy(g_exch2))
c$$$     $   call errquit('mcscf: cannot destroy exch copy',0)
c$$$       if (.not.ga_destroy(g_grad2))
c$$$     $   call errquit('mcscf: cannot destroy grad copy',0)
c
c
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c      Explicit Hessian contruction and Eigenvalue Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c$$$       if (orlen.le.hdim) then
c$$$         call makeJK( basis, nbf, nclosed, nact, g_movecs )
c$$$         call dfill((orlen*orlen),0.d0,hh,1)
c$$$         call hmat( nbf, nclosed, nact, orlen, dbl_mb(k_dm1),
c$$$     $              dbl_mb(k_dm2), g_ifock, g_afock, g_gfock,
c$$$     $              g_coul, g_exch, hh)
c$$$
c$$$         write(6,951) 
c$$$ 951     format('Hessian diagonal')
c$$$         write(6,222) (hh((i-1)*orlen+i),i=1,orlen)
c$$$ 222     format(8f16.11)
c$$$         call dsyev('V','L',orlen,hh,orlen,ev,scr,(4*hdim),info)
c$$$         write(6,953)
c$$$ 953     format('Hessian eigenvalues')
c$$$         write(6,222) (ev(i),i=1,orlen)
c$$$       endif
c
c
c  Make explicit Hessian matrix from vector products
c
c$$$       call dfill((orlen*orlen),0.d0,hh,1)
c$$$       call mcscf_hessmake( geom, basis, nbf, nclosed, nact,
c$$$     $                      oskel, orlen, g_movecs, dbl_mb(k_dm1),
c$$$     $                      dbl_mb(k_dm2), g_ifock, g_afock,
c$$$     $                      g_gfock, g_coul, g_exch, g_x, g_prod,
c$$$     $                      hh )
c
c
c  Finite difference Hessian
c
c$$$       call mcscf_fdiff_hess( geom, basis, nbf, nclosed, nact,
c$$$     $                        oskel, tol2e, dbl_mb(k_dm1),
c$$$     $                        dbl_mb(k_dm2), g_movecs, g_coul,
c$$$     $                        g_grad )
c
c  Free temporaries for debugging section
c
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              ROHF Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  ROHF occupation and trivial density
c
c$$$       call mcscf_occ2int( nbf, dbl_mb(k_occ), mclosed, mopen )
c$$$       write(6,903) mclosed, mopen
c$$$ 903   format(' ROHF Occupation',6x,'(closed):',i3,5x,'(open):',i3)
c$$$       call mcscf_rohf_modens(mopen,nact,dbl_mb(k_dm1),dbl_mb(k_dm2))
c
c
c  Fock build and 1e-Hessian vector product
c
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'CV Fock',nbf,0,g_fcv))
c$$$     $      call errquit('mcscf: cannot allocate active Fock',0)
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'PV Fock',nbf,0,g_fpv))
c$$$     $      call errquit('mcscf: cannot allocate active Fock',0)
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nbf,nbf,'CP Fock',nbf,0,g_fcp))
c$$$     $      call errquit('mcscf: cannot allocate active Fock',0)
c
c  ROHF Gradient
c
c$$$       call rohf_fock( geom, basis, nclosed, nopen, tol2e, g_movecs,
c$$$     $                 eone, etwo, g_fcv, g_fpv, g_fcp, oskel )
c$$$       call rohf_fock2grad( nbf, nclosed, nact,
c$$$     $                      g_fcv, g_fpv, g_fcp, g_grad )
c$$$       call ga_rowprint( '==== ROHF Gradient ====', g_grad )
c
c  ROHF Hessian vector product
c
c$$$       call ga_copy(g_grad, g_x)
c$$$       call ga_zero(g_prod)
c$$$       pflg = 2
c$$$       lshift = 0.d0
c$$$       call rohf_hessv_xx( basis, geom, nbf, nclosed, nact, pflg,
c$$$     $                     g_movecs, oskel, g_fcv, g_fpv, g_fcp,
c$$$     $                     tol2e, lshift, g_x, g_prod )
c$$$
c$$$       call ga_rowprint( '==== ROHF Product', g_prod )
c
c  Make ROHF Hessian
c
c$$$       call rohf_hessmake( basis, geom, nbf, nclosed, nact, 
c$$$     $                 g_movecs, oskel, g_fcv, g_fpv, g_fcp,
c$$$     $                 tol2e, g_x, g_prod )
c$$$
c
c  Create ROHF Hessian diagonal to compare with
c
c$$$       call rohf_hxxx( nbf, nclosed, nact, 0.d0, g_fcv, g_fpv,
c$$$     $                 g_fcp, g_prod )
c
c  Deallocate stuff
c
c$$$       if (.not.ga_destroy(g_fcv))
c$$$     $      call errquit('mcscf: cannot destroy MO vectors',0)
c$$$       if (.not.ga_destroy(g_fpv))
c$$$     $      call errquit('mcscf: cannot destroy MO vectors',0)
c$$$       if (.not.ga_destroy(g_fcp))
c$$$     $      call errquit('mcscf: cannot destroy MO vectors',0)
c$$$       if (.not.ga_destroy(g_grad))
c$$$     $      call errquit('mcscf: cannot destroy gradient',0)
c$$$       if (.not.ga_destroy(g_prod))
c$$$     $      call errquit('mcscf: cannot destroy product',0)
c$$$       if (.not.ga_destroy(g_x))
c$$$     $      call errquit('mcscf: cannot destroy product',0)
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              Final Cleanup
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
 401   continue
c
       if (.not.ma_pop_stack(l_sym))
     $      call errquit('mcscf: cannot pop stack?',0, MA_ERR)
       if (.not.ma_pop_stack(l_dm1))
     $      call errquit('mcscf: cannot pop stack?',0, MA_ERR)
       if (.not.ma_pop_stack(l_dm2))
     $      call errquit('mcscf: cannot pop stack?',0, MA_ERR)
       if (.not.ma_pop_stack(l_occ))
     $      call errquit('mcscf: cannot pop stack?',0, MA_ERR)
       if (.not.ga_destroy(g_exch))
     $      call errquit('mcscf: cannot destroy exchange',0, GA_ERR)
       if (.not.ga_destroy(g_coul))
     $      call errquit('mcscf: cannot destroy Coulomb',0, GA_ERR)

       if (.not.ga_destroy(g_grad))
     $      call errquit('rohf_head: cannot destroy gradient',0, GA_ERR)
       if (.not.ga_destroy(g_prod))
     $      call errquit('rohf_head: cannot destroy product',0, GA_ERR)
       if (.not.ga_destroy(g_x))
     $      call errquit('rohf_head: cannot destroy product',0, GA_ERR)
       if (.not.ga_destroy(g_afock))
     $      call errquit('rohf_head: cannot destroy afock',0, GA_ERR)
       if (.not.ga_destroy(g_ifock))
     $      call errquit('rohf_head: cannot destroy ifock',0, GA_ERR)
       if (.not.ga_destroy(g_gfock))
     $      call errquit('rohf_head: cannot destroy gfock',0, GA_ERR)

       write(6,771) 
 771   format(///,10x,40('='),
     $          /,17x,'END DEBUG SECTION',
     $          /,10x,40('='),//)
c
c
c
       return
       end









       subroutine mcscf_hessmake( geom, basis, nbf, nclosed, nact,
     $                            oskel, orlen, g_movecs, dm1, dm2,
     $                            g_ifock, g_afock, g_gfock,
     $                            g_coul, g_exch,
     $                            g_x, g_prod, hh )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
       integer geom, basis
       integer nbf, nclosed, nact
       logical oskel
       integer orlen
       double precision dm1(*), dm2(*)
       integer g_movecs, g_ifock, g_afock, g_gfock
       integer g_coul, g_exch
       integer g_x, g_prod
       double precision hh(orlen,orlen)
c
c
C       INTEGER L_HD, K_HD
       integer l_hv, k_hv, l_hy, k_hy
       integer nvir, incr, xoff, i, j, ii, pflg
       double precision xx, lshift, tol2e
       data pflg/2/
       data lshift/0.d0/
       data tol2e/1.d-12/
c
c
c
       if (.not.ma_push_get(MT_DBL, orlen, 'H', l_hv, k_hv))
     $     call errquit('mcscf: cannot allocate',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL, orlen, 'Hy', l_hy, k_hy))
     $     call errquit('mcscf: cannot allocate',0, MA_ERR)
       nvir = nbf - nclosed - nact
C       goto 10000
       open(unit=88,file='hess',form='unformatted',
     $      status='unknown')
       write(88) nclosed,nact,nvir
       write(88) (nclosed*nvir),((nclosed+nact)*nvir),orlen
       xx = 1.d0
       incr = 0
c
c
c
       do j=1,nclosed+nact
         do i=1,nvir
           ii = (j-1)*nvir + i
           call ga_zero(g_x)
           call ga_put(g_x,ii,ii,1,1,xx,1)
           call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                       tol2e, .false., pflg, lshift, dm1, dm2,
     $                       g_movecs, g_ifock, g_afock, g_gfock,
     $                       g_coul, g_exch, g_x, g_prod )
           call ga_get(g_prod,1,orlen,1,1,dbl_mb(k_hv),orlen)
           call hv_writev(88,orlen,dbl_mb(k_hv))
           call dcopy(orlen,dbl_mb(k_hv),1,hh(1,ii),1)
         enddo
       enddo
c
       xoff = nvir*(nclosed+nact)
       do j=1,nclosed
         do i=1,nact
           ii = xoff + (j-1)*nact + i
           call ga_zero(g_x)
           call ga_put(g_x,ii,ii,1,1,xx,1)
           call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                       tol2e, .false., pflg, lshift, dm1, dm2,
     $                       g_movecs, g_ifock, g_afock, g_gfock,
     $                       g_coul, g_exch, g_x, g_prod )
           call ga_get(g_prod,1,orlen,1,1,dbl_mb(k_hv),orlen)
           call hv_writev(88,orlen,dbl_mb(k_hv))
           call dcopy(orlen,dbl_mb(k_hv),1,hh(1,ii),1)
         enddo
       enddo
c
c
c
c$$$       if (.not.ma_push_get(MT_DBL, orlen, 'Hy', l_hd, k_hd))
c$$$     $     call errquit('mcscf: cannot allocate',0)
c$$$c
c$$$c
c$$$c
c$$$       open(unit=11,file='hessian.ascii',form='formatted',
c$$$     $      status='unknown')
c$$$           call hv_write(11,incr,orlen,dbl_mb(k_hv))
c$$$           call hv_writev(12,orlen,dbl_mb(k_hv))
c$$$       open(unit=12,file='hessian',form='unformatted',
c$$$     $      status='unknown')
c$$$       close(11)
c$$$       close(12)

c$$$10000  continue
c$$$       xx = 1.d0
c$$$       do i=1,orlen
c$$$         call ga_zero(g_x)
c$$$         call ga_put(g_x,i,i,1,1,xx,1)
c$$$         call mcscf_hessv( geom, basis, nbf, nclosed, nact,
c$$$     $                     tol2e, oskel, dm1, dm2, g_movecs,
c$$$     $                     g_ifock, g_afock, g_gfock,
c$$$     $                     g_coul, g_exch, g_x, g_prod )
c$$$         call ga_get(g_prod,i,i,1,1,dbl_mb(k_hd+i-1),1)
c$$$       enddo
c$$$       write(6,900) 
c$$$ 900   format('Exact Hessian diagonal')
c$$$       write(6,901) (dbl_mb(k_hd+i-1),i=1,orlen)
c$$$ 901   format(10f12.6)

c$$$       if (.not.ma_pop_stack(l_hd))
c$$$     $      call errquit('mcscf: damn',0)
c
c  Clean up
c
       if (.not.ma_pop_stack(l_hy))
     $      call errquit('mcscf: damn',0, MA_ERR)
       if (.not.ma_pop_stack(l_hv))
     $      call errquit('mcscf: damn',0, MA_ERR)

       return
       end



      
c
c  Make J and K integrals for closed + active operators
c  and dump to disk
c
       subroutine makeJK( basis, nbf, nclosed, nact, g_movecs )
       implicit none
#include "errquit.fh"
#include "global.fh"
       integer basis, nbf, nclosed, nact, g_movecs
       integer mo_lo, mo_hi, noper, nn, i
       integer g_coul, g_exch
       double precision tmp(1000)
       logical oskel
       integer ga_create_JKblocked
       external ga_create_JKblocked
       data oskel/.false./

       mo_lo = 1
       mo_hi = (nclosed+nact)
       noper = (mo_hi*(mo_hi+1))/2
       nn = nbf*nbf
       g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
       g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
       call moints_build_6x( basis, oskel,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., 16, .false., .false. )
       open(unit=88,file='JKints',status='unknown',form='unformatted')
       do i=1,noper
         call ga_get(g_coul,1,nn,i,i,tmp,1)
         call hv_writev(88,nn,tmp)
       enddo
       do i=1,noper
         call ga_get(g_exch,1,nn,i,i,tmp,1)
         call hv_writev(88,nn,tmp)
       enddo
       close(88)
       if (.not. ga_destroy(g_coul)) call errquit('mcscf: ga?',0,
     &       GA_ERR)
       if (.not. ga_destroy(g_exch)) call errquit('mcscf: ga?',0,
     &       GA_ERR)
       return
       end










#endif

