      subroutine int_2e_sifs_a( geom, basis, ablklen, tol2e, oskel,
     $     tmp, blocks, nblock, aoint2, ninfo, info, ifmt, ibvtyp, 
     &     ibitv, clabs)
      implicit none
#include "errquit.fh"
c     
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "geom.fh"
#include "cfock.fh"
#include "sym.fh"
#include "util.fh"
#include "cscfps.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
c     
      integer geom, basis
      integer ablklen
      double precision tol2e
      logical oskel
*     double precision dij(nfock*ablklen),dik(nfock*ablklen)
*     double precision dli(nfock*ablklen),djk(nfock*ablklen)
*     double precision dlj(nfock*ablklen),dlk(nfock*ablklen)
*     double precision fij(nfock*ablklen),fik(nfock*ablklen)
*     double precision fli(nfock*ablklen),fjk(nfock*ablklen)
*     double precision flj(nfock*ablklen),flk(nfock*ablklen)
      double precision tmp(ablklen)
*     integer vg_dens(nfock)
*     integer vg_fock(nfock)
      integer blocks(2,*)
      integer nblock
c
      integer ncomputed, nwritten, num
      common/sif2e/ncomputed, nwritten, num
c
      integer natoms, nproc, ijkl, next
      integer ib, jb, kb, lb, lbhi, i
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)     ! (i/j/k, lo/hi)
      integer blklen
      double precision sij, smax, skl ! schwarz info
      double precision denmax, dentol, temp, integ_acc
      logical otest
      double precision schwarz_atom_block
      external schwarz_atom_block
      logical sym_atom_block_quartet
      external sym_atom_block_quartet
c     
      integer nxtask
      external nxtask
      integer l_sifval, l_sifbuf, l2rec, n2max
      integer k_sifval, k_sifbuf 
      integer*4 ninfo, ifmt, ibvtyp, ibitv, aoint2, nipv, itypea, 
     &          itypeb, wait, nrec, reqnum, ierr, iwait
      integer*4 info(ninfo)
      integer*4 clabs(4,*)
      integer*4 msame,nmsame,nomore
      parameter(msame=0, nmsame=1, nomore= 2)

      integer l_i, l_j, l_k, l_l, l_g
      integer k_i, k_j, k_k, k_l, k_g
      logical status
      integer max2e, mem2
c
c     Figure out the threshold on which to screen just integrals
c     by dividing by maximum density matrix element.  10x for safety
c     since really wanted max over atomic blocks
c
*     denmax = 1d-10            ! To avoid a hard zero
*     do i = 1, nfock
*        call ga_maxelt(vg_dens(i), temp)
*        denmax = max(temp, denmax)
*     enddo
*     dentol = min(dentolmax,tol2e/(10.0d0*denmax)) ! For screen integs only
c
c     SEE AO_REPLICATED.F to make the next parameters consistent
c     
c     NOTE THAT ONLY USER SPECIFIED INTEGRAL ACCURACY IS NOW USED
c
*     if (intacc .eq. 0.0d0) then
*        integ_acc = min(1d-10,max(0.01d0*dentol,1d-20)) ! Variable
*         call int_acc_set(integ_acc)
*     else
cgk debug
*     write(*,*)'gk: entered int_2e_sifs_a'
cgk end
      integ_acc = intacc     ! User controlled
cgk debug
*     write(*,*)'gk: intacc = ', intacc
cgk end
*     call int_acc_set(integ_acc)
      integ_acc = 1d-10
cgk debug
*     write(*,*)'gk: integ_acc = ', integ_acc
cgk end
*     endif
c
c     How many atoms are there?
c     
      ncomputed = 0
      nwritten=0
      nrec=0
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('fock_2e: geom_ncent failed', 0, GEOM_ERR)
cgk debug
*     write(*,*)'gk: natoms = ', natoms
cgk end
      if (natoms .gt. nw_max_atom) call errquit
     $     ('fock_2e: hard dimension failed', natoms, UNKNOWN_ERR)
cgk debug
*     write(*,*)'gk: nw_max_atom = ', nw_max_atom
cgk end
c     Get sifs memory here.
c     values
      n2max=info(5)
      if (.not. MA_push_get(MT_DBL, n2max, 'int_2e_sifs_a:sifbuf',
     &     l_sifval, k_sifval))
     $     call errquit('int_2e_sifs_a: ma failed getting sifbuf', 
     &        n2max, MA_ERR)
c     buffer
      l2rec=info(4)
      if (.not. MA_push_get(MT_DBL, l2rec, 'int_2e_sifs_a:sifbuf',
     &     l_sifbuf, k_sifbuf))
     $     call errquit('int_2e_sifs_a: ma failed getting sifbuf', 
     &        l2rec, MA_ERR)

      call int_mem_2e4c(max2e, mem2)
      mem2 = max(mem2,max2e)

      status = .true.
      status = status .and.
     $     ma_push_get(MT_DBL, max2e, 'int_2e_sifs_b: buf', l_g, k_g)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'int_2e_sifs_b: buf', l_i, k_i)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'int_2e_sifs_b: buf', l_j, k_j)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'int_2e_sifs_b: buf', l_k, k_k)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'int_2e_sifs_b: buf', l_l, k_l)
      if (.not. status) call errquit('fock_2e_b: not enuf memory',0,
     &       MA_ERR)
      
c     
c     How many processors are there?
c     
      nproc = ga_nnodes()
cgk debug
*     write(*,*)'gk: nproc = ', nproc
cgk end
c     
      ijk_prev(1,1) = -1
      ijk_prev(2,1) = -1
      ijk_prev(3,1) = -1
      ijk_prev(1,2) = -1
      ijk_prev(2,2) = -1
      ijk_prev(3,2) = -1
c
      blklen = ablklen
c     
      ijkl = 0 ! this should force work to node zero
      next = nxtask(nproc, 1)
c     
c     Loop thru blocked atomic quartets
c
      smax = schwarz_max()
cgk debug
*     write(*,*)'gk: starting loop through atomic quartets'
*     write(*,*)'gk:  smax = ', smax
*     write(*,*)'gk:  next = ', next
*     write(*,*)'gk:  nblock = ', nblock
*     write(*,*)'gk:  dentol = ', dentol
cgk end
      do ib = nblock, 1, -1
         iatlo = blocks(1,ib)
         iathi = blocks(2,ib)
         do jb = 1, ib
            jatlo = blocks(1,jb)
            jathi = blocks(2,jb)
*           sij = schwarz_atom_block(iatlo,iathi,jatlo,jathi)
*
*     Cannot use symmetry of pairs since the blocking of atoms means
*     that we may have (kl) > (ij).
*
*           if (sij*smax .ge. dentol) then
               do kb = ib, 1, -1
                  katlo = blocks(1,kb)
                  kathi = blocks(2,kb)
                  lbhi = kb
                  if (ib .eq. kb) lbhi = jb
                  do lb = 1, lbhi
                     latlo = blocks(1,lb)
                     lathi = blocks(2,lb)
*                    skl = schwarz_atom_block(katlo,kathi,latlo,lathi)
*                    if (oskel) otest=sym_atom_block_quartet
*    $                    (geom, iatlo, iathi, jatlo, jathi, 
*    $                       katlo, kathi, latlo, lathi)
                     otest=.true.
c     
c     Load balance over non-zero interactions
c     
cgk but we need to need to send blocks of integrals to node zero
cgk to write out, so is it necessary to do this in parallel?
cgk debug
*     write(*,*)'gk: otest = ', otest
*     write(*,*)'gk: ijkl = ', ijkl
cgk end
 
                     if (otest .and. (ijkl .eq. next)) then
cgk debug
*     write(*,*)'gk: calling sifs_2e_task'
cgk end
                        call sifs_2e_task(
     $                       geom, basis, oskel, iatlo, jatlo, katlo, 
     &                       latlo, iathi, jathi, kathi, lathi, 
     &                       ijk_prev, tmp, tol2e, aoint2, ninfo, info, 
     &                       ifmt, ibvtyp, dbl_mb(k_sifval), ibitv, 
     &                       dbl_mb(k_sifbuf), clabs, int_mb(k_i),
     &                       int_mb(k_j), int_mb(k_k), int_mb(k_l),
     &                       dbl_mb(k_g))
cgk debug
*     write(*,*)'gk: back from sifs_2e_task'
cgk end
c
                        next = nxtask(nproc, 1)
                     end if
                     if (otest) ijkl = ijkl + 1
                  end do
               enddo
*           endif
         end do
      end do

      nipv=4
      itypea=3
      itypeb=0
      iwait=1

      if (0 .eq. ga_nodeid()) then
cgk debug
*     write(*,*)'gk: before sifew2'
*     write(*,*)'gk: nwritten = ', nwritten
*     write(*,*)'gk: num = ', num
cgk end
         nwritten=nwritten+num
         call sifew2( aoint2, info, nipv, num,  
     &        nomore, itypea, itypeb, ifmt,
     &        ibvtyp, dbl_mb(k_sifval), clabs, ibitv,
     &        dbl_mb(k_sifbuf), iwait, nrec, reqnum,
     &        ierr )
cgk debug
*     write(*,*)'gk: after sifew2'
*     write(*,*)'gk: nwritten = ', nwritten
*     write(*,*)'gk: num = ', num
cgk end
         nwritten=nwritten+num
      endif
      
      next = nxtask(-nproc, 1)
      call ga_sync()

      status = status .and. ma_pop_stack(l_l)
      status = status .and. ma_pop_stack(l_k)
      status = status .and. ma_pop_stack(l_j)
      status = status .and. ma_pop_stack(l_i)
      status = status .and. ma_pop_stack(l_g)
      if (.not. status) call errquit
     $     ('int_2e_sifs_b: pop stack failed?', 0, MA_ERR)

      if (.not. MA_pop_stack(l_sifbuf)) call errquit
     $     ('int_2e_sifs_a: pop failed at l_bfnlab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifval)) call errquit
     $     ('int_2e_sifs_a: pop failed at l_bfnlab', 0, GA_ERR)
cgk debug
*     write(*,*)'gk: leaving int_2e_sifs_a'
cgk end
c     
      end
