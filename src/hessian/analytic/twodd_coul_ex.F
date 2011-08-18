      subroutine twodd_coul_ex(
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $     d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $     d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6,
     $     d_ik7, d_jl7, d_il7, d_jk7,
     $     d_ik8, d_jl8, d_il8, d_jk8,
     $     Hxy, g_dens, g_force, blen, geom, basis, 
     $     nproc, nat, lscr, scr, leri, eri, labels, 
     $     maxq, list, q4, tol2e,
     $     nsh, oactive, oskel, scftype, omp2,
     $     nopen,nbf,pdm2,pdm2a,pdm2b,pdm2c,pdm2d,coeff,
     $     nshblocks, shmap, shglo, shghi, bfglo, bfghi, bfmap, rbfmap,
     $     bftoat, shbflo, shbfhi, jfac, kfac)
c     $Id$
c
c  This routine was essentially stolen from gradients/grad2.F.  I will 
c  work to make these routines one after I get things working.
c
c     This routine evaluates the 2-electron integral contribution to the
c     Hessian. The quantities jfac and kfac are used only to scale
c     contributions as appropriate for Coulomb fitting and DFT 
c     exchange-correlation contributions. Neither the Coulomb fitting
c     or the DFT terms are evaluated here.
c
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
#include "stdio.fh"
C-----------------------parameters---------------------------------
      integer blen              ! [input] max dimension of density block
      double precision          ! [scratch] 6 blocks per possible density
     $     d_ij(blen,blen),  d_kl(blen,blen),  d_ik(blen,blen),
     $     d_jl(blen,blen),  d_il(blen,blen),  d_jk(blen,blen),
     $     d_ij2(blen,blen), d_kl2(blen,blen), d_ik2(blen,blen), 
     $     d_jl2(blen,blen), d_il2(blen,blen), d_jk2(blen,blen),
     $     d_ij3(blen,blen), d_kl3(blen,blen), d_ik3(blen,blen), 
     $     d_jl3(blen,blen), d_il3(blen,blen), d_jk3(blen,blen),
     $     d_ij4(blen,blen), d_kl4(blen,blen), d_ik4(blen,blen), 
     $     d_jl4(blen,blen), d_il4(blen,blen), d_jk4(blen,blen)
      double precision    ! dummies needed for the so calls to grad_get_dens_all
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6,
     $     d_ik7, d_jl7, d_il7, d_jk7,
     $     d_ik8, d_jl8, d_il8, d_jk8
c     
      integer nat               ! [input] no. of atoms
      double precision Hxy(3, nat, 3, nat) ! [input/output] hessian
      integer g_dens(*)         ! [input] array of global density handles
      integer g_force           ! [input] handle to total gradients (global)
      integer basis, geom       ! [input] familiar handles
      integer nproc             ! [input] No. of processors
      integer lscr              ! [input] Size of scratch for integrals
      double precision scr(lscr) ![scratch] Scratch space for integrals
      integer leri              ! [input] Size of eri buffer
      double precision eri(78*leri) ! [scratch] Derivative integral buffer
      integer labels(leri,4)    ! [scratch] Labels for integrals
      integer maxq              ! [input] Max. no. of shell quartets
      integer list(maxq,4)      ! [scratch] Buffer for shell labels
      double precision q4(maxq) ! [scratch] Buffer for quartet prefactors
      double precision tol2e    ! [input] Integral screening threshold
      integer nsh               ! [input] Number of shells
      logical oactive(*)        ! [input] T if atom is active
      logical oskel             ! [input] T is using petite list symmetry
      character*(*) scftype     ! [input] SCF wavefunction uhf/rohf/rhf/mcscf
      logical omp2              ! [input] T if doing MP2
      integer nopen             ! [input] No. of MCSCF active orbitals
      integer nbf               ! [input] No. of basis functions (MCSCF)
      double precision pdm2(nopen,nopen,nopen,nopen) ! [scratch] MCSCF buffer
      double precision pdm2a(nopen,nopen,nopen,blen) ! [scratch] MCSCF buffer
      double precision pdm2b(nopen,nopen,blen,blen) ! [scratch] MCSCF buffer
      double precision pdm2c(nopen,blen,blen,blen) ! [scratch] MCSCF buffer
      double precision coeff(nbf,nopen) ! [input] MCSCF active space MOs
      double precision pdm2d(blen,blen,blen,blen) ! [scratch] AO 2-pdm
      integer nshblocks         ! [input] No. of shell blocks
      integer shmap(nsh)        ! [input] Map from new to old shell order
      integer shglo(nshblocks)  ! [input] First new shell in block
      integer shghi(nshblocks)  ! [input] Last new shell in block
      integer bfglo(nshblocks)  ! [input] First new basis func in block
      integer bfghi(nshblocks)  ! [input] Last new basis func in block
      integer bfmap(nbf)        ! [input] Map from new to old basis func order
      integer rbfmap(nbf)       ! [input] Map from old to new basis func order
      integer bftoat(nbf)       ! [input] Map from old basis func to atom
      integer shbflo(nsh)       ! [input] First new bf in new shell
      integer shbfhi(nsh)       ! [input] Last new bf in new shell
      double precision jfac(*), kfac(*) ! [input] DFT scaling for Coulomb/Exchange
C-------------------------local variables-----------------------------
      integer  next, nint, ijklblock
      double precision scale, q4max

      integer nxtask, task_size
      external nxtask

      double precision psum
      integer i, j, k, l
      integer iilo, jjlo, kklo, lllo
      integer iihi, jjhi, kkhi, llhi
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer lblockhi, jshtop, kshtop, lshtop
      integer ishblock, jshblock, kshblock, lshblock
      integer iish, jjsh, kksh, llsh, iiat1, iiat2, iat1, iat2
      integer nq, integ, atoms(4), ibuf, icart1, icart2, itemp
*      integer numq
      double precision smax, sij, sijkl, block_eff, p, val
*     double precision integ_acc
      double precision half
      data half /0.5d+00/
      
      integer ibflo_prev, jbflo_prev, kbflo_prev ! for caching
      
      logical oij, okl, oikjl, ouhf, orohf, omcscf, omore
      logical ogso
      logical intbdd_init4c, intbdd_2e4c
      external intbdd_init4c, intbdd_2e4c
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
      logical oprint
c
      oprint = util_print('twodd_cont',print_debug)
      smax  = schwarz_max()
      q4max = 8.0d0*(sym_number_ops(geom) + 1)*10000.0d0 ! 10000 for max(psum)
c
      call int_acc_std()
*      integ_acc = min(1d-9,max(0.01d0*tol2e/q4max,1d-17))
*      call int_acc_set(integ_acc)
c
      ibflo_prev = -1
      jbflo_prev = -1
      kbflo_prev = -1

c     MCSCF is the same as RHF except for the 2pdm piece

*      numq = 0

      omcscf= (scftype .eq. 'MCSCF')
      ouhf  = (scftype .eq. 'UHF')
      orohf = (scftype .eq. 'ROHF')
      ogso = .false.    ! this is for SO contribution for gradient, not hessian
c     
c Silent scaling to handle the "standard" values for RHF
c
      if (.not.ouhf) kfac(1) = -2.0d0*kfac(1)
c
      task_size = nshblocks*(nshblocks+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
      if(oprint) write(luout,*) 
     .     'nshblocks and task_size ',nshblocks,task_size
c     
c     Parallel loop over quartets of reordered shell blocks
c     
      ijklblock = 0
      next = nxtask(nproc,task_size)
      do ishblock = nshblocks, 1, -1
         ishlo = shglo(ishblock)
         ishhi = shghi(ishblock)
         ibflo = bfglo(ishblock)
         ibfhi = bfghi(ishblock)
         idim  = ibfhi - ibflo + 1
         if (omcscf) call grad_mcscf_a
     $        (nbf,nopen,blen,ibflo,ibfhi,coeff,pdm2,pdm2a)
         do jshblock = 1, ishblock
            jshlo = shglo(jshblock)
            jshhi = shghi(jshblock)
            jbflo = bfglo(jshblock)
            jbfhi = bfghi(jshblock)
            jdim  = jbfhi - jbflo + 1
            if (omcscf) call grad_mcscf_b
     $           (nbf,nopen,blen,idim,jbflo,jbfhi,coeff,pdm2a,pdm2b)
            do kshblock = ishblock, 1, -1
               kshlo = shglo(kshblock)
               kshhi = shghi(kshblock)
               kbflo = bfglo(kshblock)
               kbfhi = bfghi(kshblock)
               kdim  = kbfhi - kbflo + 1
               if (omcscf) call grad_mcscf_c
     $              (nbf,nopen,blen,idim,jdim,kbflo,kbfhi,
     $              coeff,pdm2b,pdm2c)
               lblockhi = kshblock
               if (ishblock .eq. kshblock) lblockhi = jshblock
               do lshblock = 1, lblockhi
                  lshlo = shglo(lshblock)
                  lshhi = shghi(lshblock)
                  lbflo = bfglo(lshblock)
                  lbfhi = bfghi(lshblock)
                  ldim  = lbfhi - lbflo + 1
                  if (next .eq. ijklblock) then
                     if (omcscf) call grad_mcscf_d
     $                    (nbf,nopen,blen,idim,jdim,kdim,ldim,
     $                    lbflo,lbfhi,coeff,pdm2c,pdm2d)
c     
c     Get blocks of the one-particle densities
c     
                     call grad_get_dens_all(omp2,ouhf,orohf,
     $                    d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk,
     $                    d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $                    d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $                    d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
     $                    d_ik5, d_jl5, d_il5, d_jk5,
     $                    d_ik6, d_jl6, d_il6, d_jk6,
     $                    d_ik7, d_jl7, d_il7, d_jk7,
     $                    d_ik8, d_jl8, d_il8, d_jk8,
     $                    blen, 
     $                    ibflo, ibfhi, jbflo, jbfhi, 
     $                    kbflo, kbfhi, lbflo, lbfhi, 
     $                    ibflo_prev, jbflo_prev, kbflo_prev,
     $                    g_dens, ogso)
c     
c     Build the list of integral shell quartets in NWChem shell labelling
c     
                     oij = ishblock.eq.jshblock
                     okl = kshblock.eq.lshblock
                     oikjl = ishblock.eq.kshblock .and. 
     $                    jshblock.eq.lshblock
c     
                     nq = 0
c     
                     do iish = ishlo, ishhi
                        atoms(1)=bftoat(bfmap(shbflo(iish)))
                        iilo = shbflo(iish)-ibflo+1
                        iihi = shbfhi(iish)-ibflo+1
                        ish = shmap(iish)
                        jshtop = jshhi
                        if (oij) jshtop = iish
                        do jjsh = jshlo, jshtop
                           atoms(2)=bftoat(bfmap(shbflo(jjsh)))
                           jjlo = shbflo(jjsh)-jbflo+1
                           jjhi = shbfhi(jjsh)-jbflo+1
                           jsh = shmap(jjsh)
                           sij = schwarz_shell(ish,jsh)
                           if (sij*smax*q4max .lt. tol2e) goto 100 ! Next jjsh
                           kshtop = kshhi
                           if (oikjl) kshtop = iish
                           do kksh = kshlo, kshtop
                              atoms(3)=bftoat(bfmap(shbflo(kksh)))
                              kklo = shbflo(kksh)-kbflo+1
                              kkhi = shbfhi(kksh)-kbflo+1
                              ksh = shmap(kksh)
                              lshtop = lshhi
                              if (okl) lshtop = kksh
                              if (oikjl .and. iish.eq.kksh) 
     $                             lshtop = jjsh
                              do llsh = lshlo, lshtop
                                 atoms(4)=bftoat(bfmap(shbflo(llsh)))
                                 lllo = shbflo(llsh)-lbflo+1
                                 llhi = shbfhi(llsh)-lbflo+1
                                 lsh = shmap(llsh)
c     
c     Check on sparsity, symmetry, active list, single center
c     
                                 sijkl = sij*schwarz_shell(ksh,lsh)
                                 if(oprint) then
               write(luout,*) 'sijkl, q4max and tol2e',sijkl,q4max,tol2e
               call util_flush(luout)
                                 endif
C     
                                 if (sijkl*q4max .lt. tol2e) goto 200 ! next llsh
                                 scale = 1.0d0
                                 if(oprint) then
               write(luout,*) 'before the oskel test'
               call util_flush(luout)
               endif
                                 if (oskel) then
                                    if (.not. sym_atom_quartet(geom,
     $                                   atoms(1), atoms(2), atoms(3), 
     $                                   atoms(4), scale)) goto 200
                                 endif
                                 if(oprint) then
               write(luout,*) 'before the oactive test'
               call util_flush(luout)
               endif
                                 if (.not. (
     $                                oactive(atoms(1)) .or.
     $                                oactive(atoms(2)) .or.
     $                                oactive(atoms(3)) .or.
     $                                oactive(atoms(4)))) goto 200
                                 if (oprint) then
               write(luout,*) 'before the atoms test ',atoms(1),
     $           atoms(2),atoms(3),atoms(4)
               call util_flush(luout)
               endif
                                 if (
     $                                atoms(1).eq.atoms(2) .and.
     $                                atoms(2).eq.atoms(3) .and.
     $                                atoms(3).eq.atoms(4)) goto 200
                                 if (oprint) then
               write(luout,*) 'after the atoms test'
               call util_flush(luout)
               endif
c     
c     Scale according to permutation symmetry of the shell labels
c     
                                 if (.not.oij .or. ish.ne.jsh) 
     $                                scale = scale*2.0d0
                                 if (.not.okl .or. ksh.ne.lsh) 
     $                                scale = scale*2.0d0
                                 if (.not.oikjl .or. ish.ne.ksh 
     $                                .or. jsh.ne.lsh) 
     $                                scale = scale*2.0d0
c     
c     Make the twopdm just for the current shell block.  Note that
c     pdm2d will have holes with junk in it for quartets we did not
c     compute ... cannot zero this since it will be too expensive.
c     
                                 call grad_make_twopdm(
     $                                omp2,ouhf,orohf,omcscf,
     $                                d_ij,  d_kl,  d_ik,  
     $                                d_jl,  d_il,  d_jk,
     $                                d_ij2, d_kl2, d_ik2,
     $                                d_jl2, d_il2, d_jk2, 
     $                                d_ij3, d_kl3, d_ik3, 
     $                                d_jl3, d_il3, d_jk3, 
     $                                d_ij4, d_kl4, d_ik4,
     $                                d_jl4, d_il4, d_jk4, 
     $                                d_ik5, d_jl5, d_il5, d_jk5,
     $                                d_ik6, d_jl6, d_il6, d_jk6,
     $                                d_ik7, d_jl7, d_il7, d_jk7,
     $                                d_ik8, d_jl8, d_il8, d_jk8,
     $                                blen, 
     $                                iilo, jjlo, kklo, lllo,
     $                                iihi, jjhi, kkhi, llhi,
     $                                pdm2d, psum, jfac, kfac, ogso)
c     
c     Final screening including magnitude of the density
c     
                                 sijkl = sijkl*psum*scale
                                 if (oprint) then
        write(luout,*) 'sijkl,psum,scale ',sijkl,psum,scale
        call util_flush(luout)
        endif
                                 if (sijkl .gt. tol2e) then
                                    nq = nq + 1
                                    q4(nq) = scale
                                    list(nq,1) = ish
                                    list(nq,2) = jsh
                                    list(nq,3) = ksh
                                    list(nq,4) = lsh
                                 endif
c     
 200                             continue
                              end do ! next llsh
                           end do
 100                       continue
                        end do  ! next jjsh
                     end do


                     if (.not. intbdd_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit
     $                    ('twodd_coul_ex:txs init?',nq, UNKNOWN_ERR)
c     
                     if (oprint) then
       write(luout,*) 'basis, lists,nq,q4,tol2e,leri,lscr ',
     $  basis, list(1,1),list(1,2),list(1,3),list(1,4),
     $  nq,q4(1),tol2e,leri,lscr 
       endif
 10                  omore = intbdd_2e4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1),labels(1,2), 
     $                    labels(1,3), labels(1,4), 
     $                    eri, leri, nint, lscr, scr)
c     
c
                     if (oprint) then
      write(luout,*) ' nint ', nint, omore
      endif
c 
                     ibuf = 1
                     if (nint .gt. 0) then
                        do integ = 1, nint
                           i = labels(integ,1)
                           j = labels(integ,2)
                           k = labels(integ,3)
                           l = labels(integ,4)
                           if (oprint) then
      write(luout,*) ' b ' ,i,j,k,l
      endif
                           atoms(1) = bftoat(i)
                           atoms(2) = bftoat(j)
                           atoms(3) = bftoat(k)
                           atoms(4) = bftoat(l)
                           if (oprint) then
      write(luout,*) ' a ' ,atoms
      endif
                           i = rbfmap(i) - ibflo + 1
                           j = rbfmap(j) - jbflo + 1
                           k = rbfmap(k) - kbflo + 1
                           l = rbfmap(l) - lbflo + 1
                           if (oprint) then
      write(luout,*) ' c ', i,j,k,l
      endif
                           p = pdm2d(l,k,j,i)
                           if (oprint) then
      write(luout,*) ' p ', p
      endif
                           do iiat1 = 1, 4
                           iat1 = atoms(iiat1)
                           do iiat2 = iiat1, 4
                              iat2 = atoms(iiat2)
                           if ((iiat1.eq.iiat2).or.
     $                       (iat1.ne.iat2)) then
                           do icart1 = 1, 3
                              itemp = 1
                              if (iiat2.eq.iiat1) itemp = icart1
c                             if (iat2.eq.iat1) itemp = icart2
                              do icart2 = itemp, 3
                              val = p*eri(ibuf)
                              if ((iat1.eq.iat2).and.
     $                          (icart1.eq.icart2)) val = val * half
                              if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      write(luout,*) 'icart1,iat1,icart2,iat2 ',icart1,iat1,icart2,iat2
      call util_flush(luout)
      endif
      endif
                              Hxy(icart1,iat1,icart2,iat2) =  
     $                             Hxy(icart1,iat1,icart2,iat2) + 
     $                             val
                              Hxy(icart2,iat2,icart1,iat1) =  
     $                            Hxy(icart2,iat2,icart1,iat1) + 
     $                            val
                              ibuf = ibuf + 1
                              end do
                           end do
                           else
                             val = p*eri(ibuf)*2.0d+00
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(1,iat1,1,iat1) = Hxy(1,iat1,1,iat1)
     $                         + val
                             val = p*(eri(ibuf+1) + eri(ibuf+3))
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(1,iat1,2,iat1) = Hxy(1,iat1,2,iat1)
     $                         + val
                             Hxy(2,iat1,1,iat1) = Hxy(2,iat1,1,iat1)
     $                         + val
                             val = p*(eri(ibuf+2) + eri(ibuf+6))
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(1,iat1,3,iat1) = Hxy(1,iat1,3,iat1)
     $                         + val
                             Hxy(3,iat1,1,iat1) = Hxy(3,iat1,1,iat1)
     $                         + val
                             val = p*eri(ibuf+4)*2.0d+00
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(2,iat1,2,iat1) = Hxy(2,iat1,2,iat1)
     $                         + val
                             val = p*(eri(ibuf+5) + eri(ibuf+7))
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(2,iat1,3,iat1) = Hxy(2,iat1,3,iat1)
     $                         + val
                             Hxy(3,iat1,2,iat1) = Hxy(3,iat1,2,iat1)
     $                         + val
                             val = p *eri(ibuf+8)*2.0d+00
                             if (oprint) then
      if (abs(val).gt. 1.0d+00) then
      write(luout,*) 'p and eri is ',p,eri(ibuf),ibuf
      write(luout,*) 'val is ',val
      call util_flush(luout)
      endif
      endif
                             Hxy(3,iat1,3,iat1) = Hxy(3,iat1,3,iat1)
     $                         + val
                             ibuf = ibuf + 9
                           endif
                           end do
                           end do
                        end do    ! 1 to nint
                     end if
                     if (oprint) then
                 if (omore) write(luout,*) 'going back for more ints'
                     call util_flush(luout)
                     endif
                     if (omore) goto 10
c     
                     next = nxtask(nproc,task_size)
                  end if
c     
                  ijklblock = ijklblock + 1
c     
               end do
            end do
         end do
      end do
c
      call int_acc_std()
c
      next = nxtask(-nproc,task_size)
      call ga_sync()
c     
c Silent scaling back to the "standard" values for RHF
c
      if (.not.ouhf) kfac(1) = kfac(1)/(-2.0d0)
c
      return
      end
