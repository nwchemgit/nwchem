c
c     Two-electron part (modified /src/gradients/grad2.F)
c
      subroutine tddft_grad_grad_2e(
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $     d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $     d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4,
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6,
     $     d_ik7, d_jl7, d_il7, d_jk7,
     $     d_ik8, d_jl8, d_il8, d_jk8,
     $     d_ik9, d_jl9, d_il9, d_jk9,
     $     frc, g_dens, g_force, blen, geom, basis, 
     $     nproc, nat, lscr, scr, leri, eri, labels,    
     $     maxq, list, q4, tol2e,
     $     nsh, oactive, oskel, scftype, omp2,
     $     nopen,nbf,pdm2,pdm2a,pdm2b,pdm2c,pdm2d,coeff,
     $     nshblocks, shmap, shglo, shghi, bfglo, bfghi, bfmap, rbfmap,
     $     bftoat, shbflo, shbfhi, jfac, kfac, exso)
c
      implicit none
c
#include "errquit.fh"
#include "sym.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
c
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
     $     d_jl4(blen,blen), d_il4(blen,blen), d_jk4(blen,blen),
     $     d_ik5(blen,blen), 
     $     d_jl5(blen,blen), d_il5(blen,blen), d_jk5(blen,blen),
     $     d_ik6(blen,blen), 
     $     d_jl6(blen,blen), d_il6(blen,blen), d_jk6(blen,blen),
     $     d_ik7(blen,blen), 
     $     d_jl7(blen,blen), d_il7(blen,blen), d_jk7(blen,blen),
     $     d_ik8(blen,blen), 
     $     d_jl8(blen,blen), d_il8(blen,blen), d_jk8(blen,blen),
     $     d_ik9(blen,blen), 
     $     d_jl9(blen,blen), d_il9(blen,blen), d_jk9(blen,blen)
c     
      integer nat               ! [input] no. of atoms
      double precision frc(3, nat) ! [input/output] 2-e piece of gradient
      integer g_dens(*)         ! [input] array of global density handles
      integer g_force           ! [input] handle to total gradients (global)
      integer basis, geom       ! [input] familiar handles
      integer nproc             ! [input] No. of processors
      integer lscr              ! [input] Size of scratch for integrals
      double precision scr(lscr) ![scratch] Scratch space for integrals
      integer leri              ! [input] Size of eri buffer
      double precision eri(3,4,leri) ! [scratch] Derivative integral buffer
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
      double precision jfac, kfac ! [input] DFT scaling for Coulomb/Exchange
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
      integer iish, jjsh, kksh, llsh, iiat, iat
      integer nq, integ, atoms(4)
*      integer numq
      double precision smax, sij, sijkl, block_eff, p
*     double precision integ_acc
      
      integer ibflo_prev, jbflo_prev, kbflo_prev ! for caching
      
      logical oij, okl, oikjl, ouhf, orohf, omcscf, omore, exso
      logical otddft, otda
      logical otriplet
      integer ipol
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
cDEBUG
      logical oroot
      oroot = ga_nodeid().eq.0
cDEBUG
c     
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

c     MCSCF is the same as RHF except for the 2pdm piece

*      numq = 0
c
      omcscf = (scftype .eq. 'MCSCF')
      ouhf   = (scftype .eq. 'UHF')
      orohf  = (scftype .eq. 'ROHF')
      otda   = (scftype .eq. 'STDADFT') .or. (scftype .eq. 'TTDADFT')
     +    .or. (scftype .eq. 'UTDADFT')
      otddft = (scftype .eq. 'STDDFT') .or. (scftype .eq. 'TTDDFT')
     +    .or. (scftype .eq. 'UTDDFT') .or. otda
      if ((scftype .eq. 'STDADFT').or.(scftype .eq. 'STDDFT').or.
     +    (scftype .eq. 'TTDADFT').or.(scftype .eq. 'TTDDFT')) then
        ipol = 1
      else
        ipol = 2
      endif
      if ((scftype .eq. 'TTDADFT').or.(scftype .eq. 'TTDDFT')) then
        otriplet = .true.
      else
        otriplet = .false.
      endif
cDEBUG
c     write(*,*)'*** grad2: scftype   = ',scftype
c     write(*,*)'*** grad2: otda      = ',otda
c     write(*,*)'*** grad2: otddft    = ',otddft
c     write(*,*)'*** grad2: ipol      = ',ipol
c     write(*,*)'*** grad2: nsh       = ',nsh
c     write(*,*)'*** grad2: nshblocks = ',nshblocks
cDEBUG
c
cDEBUG
c     write(*,*)'*** grad2: scftype   = ',scftype
c     write(*,*)'*** grad2: otda      = ',otda
c     write(*,*)'*** grad2: otddft    = ',otddft
c     write(*,*)'*** grad2: ipol      = ',ipol
c     write(*,*)'*** grad2: nsh       = ',nsh
c     write(*,*)'*** grad2: nshblocks = ',nshblocks
cDEBUG
c
      task_size = nshblocks*(nshblocks+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
cDEBUG
c     if (oroot) write(*,*)'*** task_size truncation disabled'
c     if (oroot) write(*,*)'*** task_size truncation enabled'
cDEBUG
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
cDEBUG
c     write(*,*)"*** 2dm_lkji(",
c    +                     lbflo,":",lbfhi,",",kbflo,":",kbfhi,",",
c    +                     jbflo,":",jbfhi,",",ibflo,":",ibfhi,")"
cDEBUG
                     call tddft_grad_get_dens(otda,ipol,kfac,
     $                      d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk,
     $                      d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $                      d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $                      d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
     $                      d_ik5, d_jl5, d_il5, d_jk5,
     $                      d_ik6, d_jl6, d_il6, d_jk6, 
     $                      d_ik7, d_jl7, d_il7, d_jk7, 
     $                      d_ik8, d_jl8, d_il8, d_jk8, 
     $                      d_ik9, d_jl9, d_il9, d_jk9, 
     $                      blen, 
     $                      ibflo, ibfhi, jbflo, jbfhi, 
     $                      kbflo, kbfhi, lbflo, lbfhi, 
     $                      g_dens)
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
cDEBUG
                           if (iihi-iilo+1.le.0 .or. 
     +                         jjhi-jjlo+1.le.0) goto 100
cDEBUG
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
c     
                                 if (sijkl*q4max .lt. tol2e) goto 200 ! next llsh
cDEBUG
                                 if (kkhi-kklo+1.le.0 .or. 
     +                               llhi-lllo+1.le.0) goto 200
cDEBUG
cDEBUG
c     write(*,*)'*** passed schwarz'
cDEBUG
                                 scale = 1.0d0
                                 if (oskel) then
                                    if (.not. sym_atom_quartet(geom,
     $                                   atoms(1), atoms(2), atoms(3), 
     $                                   atoms(4), scale)) goto 200
                                 endif
cDEBUG
c     write(*,*)'*** passed symmetry'
cDEBUG
                                 if (.not. (
     $                                oactive(atoms(1)) .or.
     $                                oactive(atoms(2)) .or.
     $                                oactive(atoms(3)) .or.
     $                                oactive(atoms(4)))) goto 200
cDEBUG
c     write(*,*)'*** passed oactive'
cDEBUG
                                 if (
     $                                atoms(1).eq.atoms(2) .and.
     $                                atoms(2).eq.atoms(3) .and.
     $                                atoms(3).eq.atoms(4)) goto 200
cDEBUG
c     write(*,*)'*** passed atoms'
cDEBUG
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
cDEBUG
                                 call tddft_grad_make_2pdm(otda,
     $                                  otriplet,ipol,jfac,kfac,
     $                                  d_ij,  d_kl,  d_ik,  
     $                                  d_jl,  d_il,  d_jk,
     $                                  d_ij2, d_kl2, d_ik2,
     $                                  d_jl2, d_il2, d_jk2, 
     $                                  d_ij3, d_kl3, d_ik3, 
     $                                  d_jl3, d_il3, d_jk3, 
     $                                  d_ij4, d_kl4, d_ik4,
     $                                  d_jl4, d_il4, d_jk4, 
     $                                  d_ik5, d_jl5, d_il5, d_jk5,
     $                                  d_ik6, d_jl6, d_il6, d_jk6,
     $                                  d_ik7, d_jl7, d_il7, d_jk7,
     $                                  d_ik8, d_jl8, d_il8, d_jk8,
     $                                  d_ik9, d_jl9, d_il9, d_jk9,
     $                                  blen, 
     $                                  iilo, jjlo, kklo, lllo,
     $                                  iihi, jjhi, kkhi, llhi,
     $                                  pdm2d, psum)
c     
c     Final screening including magnitude of the density
c     
                                 sijkl = sijkl*psum*scale
                                 if (sijkl .gt. tol2e) then
                                    nq = nq + 1
                                    q4(nq) = scale
                                    list(nq,1) = ish
                                    list(nq,2) = jsh
                                    list(nq,3) = ksh
                                    list(nq,4) = lsh
                                 endif
c     
 200                          end do ! next llsh
                           end do
 100                    end do  ! next jjsh
                     end do

*                     numq = numq + nq

c     
                     if (.not. intbd_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit('grad2:txs init?',nq,
     &       INT_ERR)
c     
 10                  omore = intbd_2e4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1),labels(1,2), 
     $                    labels(1,3), labels(1,4), 
     $                    eri, leri, nint, lscr, scr)
c     
*     write(6,*) ' nint ', nint, omore
c     
                     if (nint .gt. 0) then
                        do integ = 1, nint
                           i = labels(integ,1)
                           j = labels(integ,2)
                           k = labels(integ,3)
                           l = labels(integ,4)
*     write(6,*) ' b ' ,i,j,k,l
                           atoms(1) = bftoat(i)
                           atoms(2) = bftoat(j)
                           atoms(3) = bftoat(k)
                           atoms(4) = bftoat(l)
*     write(6,*) ' a ' ,atoms
                           i = rbfmap(i) - ibflo + 1
                           j = rbfmap(j) - jbflo + 1
                           k = rbfmap(k) - kbflo + 1
                           l = rbfmap(l) - lbflo + 1
*     write(6,*) ' c ', i,j,k,l
                           p = pdm2d(l,k,j,i)
*     write(6,*) ' p ', p
                           do iiat = 1, 4
                              iat = atoms(iiat)
                              frc(1,iat) =  frc(1,iat) + 
     $                             p*eri(1,iiat,integ)
                              frc(2,iat) =  frc(2,iat) + 
     $                             p*eri(2,iiat,integ)
                              frc(3,iat) =  frc(3,iat) + 
     $                             p*eri(3,iiat,integ)
                           end do
                        end do
                     end if
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
*      write(6,*) ' numq ', numq
c     
      end
c $Id$
