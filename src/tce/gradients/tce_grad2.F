      subroutine tce_grad2(d_2pdm_ao,k_2pdm_ao_offset,
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $     d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $     d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4,
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6,
     $     d_ik7, d_jl7, d_il7, d_jk7,
     $     d_ik8, d_jl8, d_il8, d_jk8,
     $     frc, g_dens, g_force, blen, geom, basis, 
     $     nproc, nat, lscr, scr, leri, eri, labels,    
     $     maxq, list, q4, tol2e,
     $     nsh, oactive, oskel, scftype, omp2,
     $     nopen,nbf,pdm2,pdm2a,pdm2b,pdm2c,pdm2d,coeff,
     $     nshblocks, shmap, shglo, shghi, bfglo, bfghi, bfmap, rbfmap,
     $     bftoat, shbflo, shbfhi, jfac, kfac, exso)
c     $Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "apiP.fh"
#include "inp.fh"
C-----------------------parameters---------------------------------
c->pengdong.fan
      integer d_2pdm_ao,k_2pdm_ao_offset,k_a,l_a,dima
c<-pengdong.fan
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
     $     d_jl8(blen,blen), d_il8(blen,blen), d_jk8(blen,blen)
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
      integer i, j, k, l,i1,j1,k1,l1
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
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
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

      omcscf= (scftype .eq. 'MCSCF')
      ouhf  = (scftype .eq. 'UHF')
      orohf = (scftype .eq. 'ROHF')

      task_size = nshblocks*(nshblocks+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
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
c         if (omcscf) call grad_mcscf_a
c     $        (nbf,nopen,blen,ibflo,ibfhi,coeff,pdm2,pdm2a)
         do jshblock = 1, ishblock
            jshlo = shglo(jshblock)
            jshhi = shghi(jshblock)
            jbflo = bfglo(jshblock)
            jbfhi = bfghi(jshblock)
            jdim  = jbfhi - jbflo + 1
c            if (omcscf) call grad_mcscf_b
c     $           (nbf,nopen,blen,idim,jbflo,jbfhi,coeff,pdm2a,pdm2b)
            do kshblock = ishblock, 1, -1
               kshlo = shglo(kshblock)
               kshhi = shghi(kshblock)
               kbflo = bfglo(kshblock)
               kbfhi = bfghi(kshblock)
               kdim  = kbfhi - kbflo + 1
c               if (omcscf) call grad_mcscf_c
c     $              (nbf,nopen,blen,idim,jdim,kbflo,kbfhi,
c     $              coeff,pdm2b,pdm2c)
               lblockhi = kshblock
               if (ishblock .eq. kshblock) lblockhi = jshblock
               do lshblock = 1, lblockhi
                  lshlo = shglo(lshblock)
                  lshhi = shghi(lshblock)
                  lbflo = bfglo(lshblock)
                  lbfhi = bfghi(lshblock)
                  ldim  = lbfhi - lbflo + 1
                  if (next .eq. ijklblock) then
c                     if (omcscf) call grad_mcscf_d
c     $                    (nbf,nopen,blen,idim,jdim,kdim,ldim,
c     $                    lbflo,lbfhi,coeff,pdm2c,pdm2d)
c     
c     Get blocks of the one-particle densities
c     
c                     call grad_get_dens_all(omp2,ouhf,orohf,
c     $                    d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk,
c     $                    d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
c     $                    d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
c     $                    d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
c     $                    d_ik5, d_jl5, d_il5, d_jk5,
c     $                    d_ik6, d_jl6, d_il6, d_jk6, 
c     $                    d_ik7, d_jl7, d_il7, d_jk7, 
c     $                    d_ik8, d_jl8, d_il8, d_jk8, 
c     $                    blen, 
c     $                    ibflo, ibfhi, jbflo, jbfhi, 
c     $                    kbflo, kbfhi, lbflo, lbfhi, 
c     $                    ibflo_prev, jbflo_prev, kbflo_prev,
c     $                    g_dens, exso)
c     
c     Get blocks of the two-particle densities
                     dima=idim*jdim*kdim*ldim
                     if (.not.ma_push_get(mt_dbl,dima,'noname',
     1               l_a,k_a))
     1               call errquit('tce_grad2: ma problem',3,ma_err)
c
                     call get_hash_block(d_2pdm_ao,dbl_mb(k_a),
     1                    dima,
     1                    int_mb(k_2pdm_ao_offset),
     1                    (ishblock-1)+(kshblock-1)*nshblocks+
     1                    (jshblock-1)*nshblocks*nshblocks+
     1                    (lshblock-1)*nshblocks*nshblocks*nshblocks)
c
c                      print *,'PDM2D FOR I,J,K,L SHBLOCK'
c                      print *, ishblock,jshblock,kshblock,lshblock
c
                      do i=1,idim
                      do j=1,jdim
                      do k=1,kdim
                      do l=1,ldim
                         pdm2d(l,k,j,i) = dbl_mb(k_a+
     1                   (i-1)+(k-1)*idim+(j-1)*idim*kdim
     1                  +(l-1)*idim*kdim*jdim)
c                         write(6,'(a,i5,i5,i5,i5,a,f20.16)') 'PDM2D',
c     1                   l,k,j,i,'=',pdm2d(l,k,j,i)
                      enddo
                      enddo
                      enddo
                      enddo
                      if(.not.ma_pop_stack(l_a))
     1                 call errquit('tce_grad2: ma problem',4,ma_err)
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
c     
                                 if (sijkl*q4max .lt. tol2e) goto 200 ! next llsh
                                 scale = 1.0d0
                                 if (oskel) then
                                    if (.not. sym_atom_quartet(geom,
     $                                   atoms(1), atoms(2), atoms(3), 
     $                                   atoms(4), scale)) goto 200
                                 endif
                                 if (.not. (
     $                                oactive(atoms(1)) .or.
     $                                oactive(atoms(2)) .or.
     $                                oactive(atoms(3)) .or.
     $                                oactive(atoms(4)))) goto 200
                                 if (
     $                                atoms(1).eq.atoms(2) .and.
     $                                atoms(2).eq.atoms(3) .and.
     $                                atoms(3).eq.atoms(4)) goto 200
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
c->pengdong.fan
c                                 call tce_grad_make_twopdm(
c     $                                omp2,ouhf,orohf,omcscf,
c     $                                d_ij,  d_kl,  d_ik,
c     $                                d_jl,  d_il,  d_jk,
c     $                                d_ij2, d_kl2, d_ik2,
c     $                                d_jl2, d_il2, d_jk2,
c     $                                d_ij3, d_kl3, d_ik3,
c     $                                d_jl3, d_il3, d_jk3,
c     $                                d_ij4, d_kl4, d_ik4,
c     $                                d_jl4, d_il4, d_jk4,
c     $                                d_ik5, d_jl5, d_il5, d_jk5,
c     $                                d_ik6, d_jl6, d_il6, d_jk6,
c     $                                d_ik7, d_jl7, d_il7, d_jk7,
c     $                                d_ik8, d_jl8, d_il8, d_jk8,
c     $                                blen,
c     $                                iilo, jjlo, kklo, lllo,
c     $                                iihi, jjhi, kkhi, llhi,
c     $                                pdm2d, psum, jfac, kfac, exso)
c<-pengdong.fan
c                                 do i=iilo,iihi
c                                    i1=bfmap(ibflo+i-1)
c                                    do j=jjlo,jjhi
c                                       j1=bfmap(jbflo+j-1)
c                                       do k=kklo,kkhi
c                                          k1=bfmap(kbflo+k-1)
c                                          do l=lllo,llhi
c                                            l1=bfmap(lbflo+l-1)
cc            pdm2d(l,k,j,i)=dbl_mb(k_2pdm_ao+
cc     1      (i1-1)+(j1-1)*nbf+(k1-1)*nbf*nbf+(l1-1)*nbf*nbf*nbf)
c            pdm2d(l,k,j,i)=dbl_mb(k_2pdm_ao+
c     1      (ibflo+i-1-1)+(jbflo+j-1-1)*nbf+(kbflo+k-1-1)*nbf*nbf
c     1     +(lbflo+l-1-1)*nbf*nbf*nbf)
c
c                                          enddo
c                                       enddo
c                                    enddo
c                                 enddo
                                 psum=0.d0
                                 do i=iilo,iihi
                                    do j=jjlo,jjhi
                                       do k=kklo,kkhi
                                          do l=lllo,llhi
                                             psum=psum+pdm2d(l,k,j,i)**2
                                          enddo
                                       enddo
                                    enddo
                                 enddo
c     
c     Final screening including magnitude of the density
c     
                                 sijkl = sijkl*psum*scale
                                 if (sijkl .gt. tol2e) then
c                                    write(6,*) 'scale=',scale
                                    nq = nq + 1
                                    q4(nq) = scale
                                    list(nq,1) = ish
                                    list(nq,2) = jsh
                                    list(nq,3) = ksh
                                    list(nq,4) = lsh
                                 endif
c     
  200                         end do ! next llsh
                           end do
  100                   end do  ! next jjsh
                     end do

*                     numq = numq + nq

c                     write(6,*) 'Is this line in any loops?!!!'

c     
                     if (.not. intbd_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit('grad2:txs init?',nq,
     &       INT_ERR)
c     
   10                omore = intbd_2e4c(
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
c                           write(6,*)'bbbbbbbbb'
c                           write(6,'(a,4i4)') 'l,k,j,i',l,k,j,i
                           atoms(1) = bftoat(i)
                           atoms(2) = bftoat(j)
                           atoms(3) = bftoat(k)
                           atoms(4) = bftoat(l)
*     write(6,*) ' a ' ,atoms
c->pdfan
c                           p = dbl_mb(k_2pdm_ao+(i-1)+(j-1)*nbf
c     1                         +(k-1)*nbf*nbf+(l-1)*nbf*nbf*nbf)
c<-pdfan
                           i = rbfmap(i) - ibflo + 1
                           j = rbfmap(j) - jbflo + 1
                           k = rbfmap(k) - kbflo + 1
                           l = rbfmap(l) - lbflo + 1
*     write(6,*) ' c ', i,j,k,l
                           p = pdm2d(l,k,j,i)
c                           write(6,*) 'ccccccccccccc'
c                           write(6,'(a,i4,i4,i4,i4,f20.16)') '2pdm_ao',
c     1                     l,k,j,i,p
*     write(6,*) ' p ', p
                           do iiat = 1, 4
                              iat = atoms(iiat)
                              frc(1,iat) =  frc(1,iat) + 
     $                             p*eri(1,iiat,integ)
                              frc(2,iat) =  frc(2,iat) + 
     $                             p*eri(2,iiat,integ)
                              frc(3,iat) =  frc(3,iat) + 
     $                             p*eri(3,iiat,integ)
c->debug
c                              write(6,'(a,2i4,f20.16)')'iiat,integ,p',
c     1                         iiat,integ,p
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iiat,integ, eri(1)',
c     1                         iiat,integ,eri(1,iiat,integ)
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iat,integ, frc(1)',
c     1                         iat,integ,frc(1,iat)
c
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iiat,integ, eri(2)',
c     1                         iiat,integ,eri(2,iiat,integ)
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iat,integ, frc(2)',
c     1                         iat,integ,frc(2,iat)
c
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iiat,integ, eri(3)',
c     1                         iiat,integ,eri(3,iiat,integ)
c                              write(6,'(a,2i4,f20.16)') 
c     1                         'iat,integ, frc(3)',
c     1                         iat,integ,frc(3,iat)
c<-debug
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
c      subroutine grad_mcscf_a
c     $     (nbf,nopen,blen,plo,phi,coeff,pdm2,pdm2a)
c      implicit none
c      integer nbf, nopen, blen
c      integer plo, phi
c      double precision coeff(nbf,nopen)
c      double precision pdm2(nopen,nopen,nopen,nopen)
c      double precision pdm2a(nopen,nopen,nopen,blen)
cc     
cc     pdm2a(l,k,j,p) = sum(i) C(p,i) pdm2(l,k,j,i)
cc     
cc     for p in given shell(plo<=p<=phi)
cc     
c      integer i, j, k, l        ! Active MOs
c      integer p                 ! AO
c      double precision cpi
cc     
c      do p = 1,(phi-plo+1)
c         do j = 1, nopen
c            do k = 1, nopen
c               do l = 1, nopen
c                  pdm2a(l,k,j,p) = 0.0d0
c               end do
c            end do
c         end do
c      end do
cc     
c      do i = 1, nopen
c         do p =  1,(phi-plo+1)
c            cpi = coeff(p+plo-1,i)
c            if (abs(cpi) .gt. 1d-10) then
c               do j = 1, nopen
c                  do k = 1, nopen
c                     do l = 1, nopen
c                        pdm2a(l,k,j,p) = pdm2a(l,k,j,p) + 
c     $                       pdm2(l,k,j,i)*cpi
c                     end do
c                  end do
c               end do
c            end if
c         end do
c      end do
cc     
c      end
c      subroutine  grad_mcscf_b
c     $     (nbf,nopen,blen,pdim,qlo,qhi,coeff,pdm2a,pdm2b)
c      implicit none
c      integer nbf, nopen, blen
c      integer pdim,qlo, qhi
c      double precision coeff(nbf,nopen)
c      double precision pdm2a(nopen,nopen,nopen,blen)
c      double precision pdm2b(nopen,nopen,blen,blen)
cc     
cc     pdm2b(l,k,q,p) = sum(j) C(q,j) pdm2a(l,k,j,p)
cc     
cc     for q in given shell(qlo<=q<=qhi)
cc     
c      integer j, k, l           ! Active MOs
c      integer p, q              ! AO
c      double precision cqj
cc     
c      do p = 1, pdim
c         do q = 1,(qhi-qlo+1)
c            do k = 1, nopen
c               do l = 1, nopen
c                  pdm2b(l,k,q,p) = 0.0d0
c               end do
c            end do
c         end do
c      end do
cc     
c      do p =  1, pdim
c         do q = 1,(qhi-qlo+1)
c            do j = 1, nopen
c               cqj = coeff(q+qlo-1,j)
c               if (abs(cqj) .gt. 1d-10) then
c                  do k = 1, nopen
c                     do l = 1, nopen
c                        pdm2b(l,k,q,p) = pdm2b(l,k,q,p) + 
c     $                       pdm2a(l,k,j,p)*cqj
c                     end do
c                  end do
c               end if
c            end do
c         end do
c      end do
cc     
c      end
c      subroutine  grad_mcscf_c
c     $     (nbf,nopen,blen,pdim,qdim,rlo,rhi,coeff,pdm2b,pdm2c)
c      implicit none
c      integer nbf, nopen, blen
c      integer pdim,qdim,rlo,rhi
c      double precision coeff(nbf,nopen)
c      double precision pdm2b(nopen,nopen,blen,blen)
c      double precision pdm2c(nopen,blen,blen,blen)
cc     
cc     pdm2c(l,r,q,p) = sum(k) C(r,k) pdm2b(l,k,q,p)
cc     
cc     for r in given shell(rlo<=r<=rhi)
cc     
c      integer k, l              ! Active MOs
c      integer p, q, r           ! AO
c      double precision crk
cc     
c      do p = 1, pdim
c         do q = 1, qdim
c            do r = 1,(rhi-rlo+1)
c               do l = 1, nopen
c                  pdm2c(l,r,q,p) = 0.0d0
c               end do
c            end do
c         end do
c      end do
cc     
c      do p =  1, pdim
c         do q = 1, qdim
c            do r = 1,(rhi-rlo+1)
c               do k = 1, nopen
c                  crk = coeff(r+rlo-1,k)
c                  if (abs(crk) .gt. 1d-10) then
c                     do l = 1, nopen
c                        pdm2c(l,r,q,p) = pdm2c(l,r,q,p) + 
c     $                       pdm2b(l,k,q,p)*crk
c                     end do
c                  end if
c               end do
c            end do
c         end do
c      end do
cc     
c      end
c      subroutine  grad_mcscf_d
c     $     (nbf,nopen,blen,pdim,qdim,rdim,sdim,slo,shi,
c     $     coeff,pdm2c,pdm2d)
c      implicit none
c      integer nbf, nopen, blen
c      integer pdim,qdim,rdim,sdim,slo,shi
c      double precision coeff(nbf,nopen)
c      double precision pdm2c(nopen,blen,blen,blen)
c      double precision pdm2d(blen,blen,blen,blen)
cc     
cc     pdm2d(s,r,q,p) = sum(l) C(s,l) pdm2c(l,r,q,p)
cc     
cc     for s in given shell(slo<=s<=shi)
cc     
c      integer l                 ! Active MOs
c      integer p, q, r, s        ! AO
c      double precision pdm2clrqp
cc     
c      do p =  1, pdim
c         do q = 1, qdim
c            do r = 1, rdim
c               do s = 1, sdim
c                  pdm2d(s,r,q,p) = 0.0d0
c               end do
c               do l = 1, nopen
c                  pdm2clrqp = pdm2c(l,r,q,p)
c                  if (abs(pdm2clrqp) .gt. 1d-10) then
c                     do s = 1, sdim
c                        pdm2d(s,r,q,p) = pdm2d(s,r,q,p) + 
c     $                       pdm2clrqp*coeff(s+slo-1,l)
c                     end do
c                  end if
c               end do
c            end do
c         end do
c      end do
cc     
c      end
c      subroutine grad_make_twopdm(omp2,ouhf,orohf,omcscf,
c     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk,
c     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
c     $     d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
c     $     d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4,
c     &     d_ik5, d_jl5, d_il5, d_jk5, 
c     &     d_ik6, d_jl6, d_il6, d_jk6, 
c     &     d_ik7, d_jl7, d_il7, d_jk7, 
c     &     d_ik8, d_jl8, d_il8, d_jk8, 
c     $     blen, ilo, jlo, klo, llo, ihi, jhi, khi, lhi, pdm2d, psum,
c     $     jfac, kfac, exso)
c      implicit none
c#include "errquit.fh"
c      logical omp2,ouhf,orohf,omcscf,exso   ! [input] Wavefunction info
c      integer blen
c      double precision          ! [input] AO density blocks
c     $     d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
c     $     d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen),
c     $     d_ij2(blen,blen), d_kl2(blen,blen),
c     $     d_ik2(blen,blen), d_jl2(blen,blen),
c     $     d_il2(blen,blen), d_jk2(blen,blen),
c     $     d_ij3(blen,blen), d_kl3(blen,blen),
c     $     d_ik3(blen,blen), d_jl3(blen,blen),
c     $     d_il3(blen,blen), d_jk3(blen,blen),
c     $     d_ij4(blen,blen), d_kl4(blen,blen),
c     $     d_ik4(blen,blen), d_jl4(blen,blen),
c     $     d_il4(blen,blen), d_jk4(blen,blen),
c     $     d_ik5(blen,blen), d_jl5(blen,blen),
c     $     d_il5(blen,blen), d_jk5(blen,blen),
c     $     d_ik6(blen,blen), d_jl6(blen,blen),
c     $     d_il6(blen,blen), d_jk6(blen,blen),
c     $     d_ik7(blen,blen), d_jl7(blen,blen),
c     $     d_il7(blen,blen), d_jk7(blen,blen),
c     $     d_ik8(blen,blen), d_jl8(blen,blen),
c     $     d_il8(blen,blen), d_jk8(blen,blen)
c
c      double precision pdm2d(blen,blen,blen,blen) ! [output] 2pdm
c      double precision psum     ! [output] Norm of this block of twopdm
c      integer ilo, jlo, klo, llo, ihi, jhi, khi, lhi
cc
c      double precision jfac, kfac ! [input] Factors for DFT
cc
c      double precision den2
c      integer i, j, k, l
c      double precision half, quarter, eigth, sixteenth
c      parameter(half = 0.5D0, quarter= 0.25D0, eigth = 0.125D0,
c     $     sixteenth = 0.0625D0)
cc     
cc     Compute the shell block of the 2-pdm
cc     
c      if (omp2) then  
c         if (ouhf) then         ! UHF + MP2
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = quarter *(
c     $                       (d_ij(i,j)  + d_ij2(i,j))
c     $                       *(d_kl3(k,l) + d_kl4(k,l))
c     $                       +(d_ij3(i,j) + d_ij4(i,j))
c     $                       *(d_kl(k,l) + d_kl2(k,l)))
c     $                       - eigth *(d_il(i,l) * d_jk3(j,k) 
c     $                       + d_jk(j,k) * d_il3(i,l)
c     $                       + d_jl(j,l) * d_ik3(i,k)
c     $                       + d_ik(i,k) * d_jl3(j,l)
c     $                       + d_il2(i,l) * d_jk4(j,k) 
c     $                       + d_jk2(j,k) * d_il4(i,l)
c     $                       + d_jl2(j,l) * d_ik4(i,k)
c     $                       + d_ik2(i,k) * d_jl4(j,l))
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         else if (orohf) then
c            call errquit('only RHF or UMP2 grads available',110,
c     &       CAPMIS_ERR)
c         else                   ! RHF + MP2
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = quarter*(d_ij(i,j)*d_kl2(k,l)
c     $                       +  d_ij2(i,j) * d_kl(k,l))
c     $                       - sixteenth *(d_il(i,l) * d_jk2(j,k) 
c     $                       + d_jk(j,k) * d_il2(i,l)
c     $                       + d_jl(j,l) * d_ik2(i,k)
c     $                       + d_ik(i,k) * d_jl2(j,l))
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         end if
c      else
c         if (ouhf) then         ! UHF
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = half*jfac*(d_ij(i,j) * d_kl(k,l))
c     $                       - eigth*kfac*(d_ik(i,k) * d_jl(j,l)
c     $                       +  d_il(i,l) * d_jk(j,k)
c     $                       +  d_ik2(i,k) * d_jl2(j,l)
c     $                       +  d_il2(i,l) * d_jk2(j,k))
cc
c                        if(exso)then  !d_lj=-d_jl, d_kj=-d_jk, i*i=-1
c                           den2 = den2 + eigth*kfac*(
c     &                          -  d_ik3(i,k) * d_jl3(j,l) ! Im(aa) and Im(bb)
c     $                          -  d_il3(i,l) * d_jk3(j,k)
c     $                          -  d_ik4(i,k) * d_jl4(j,l)
c     $                          -  d_il4(i,l) * d_jk4(j,k))
cc     
c                           den2 = den2 - eigth*kfac*(
c     &                          -  d_ik5(i,k) * d_jl5(j,l) ! Re(ab) and Re(ba)
c     $                          -  d_il5(i,l) * d_jk5(j,k)
c     $                          +  d_ik6(i,k) * d_jl6(j,l)
c     $                          +  d_il6(i,l) * d_jk6(j,k))
cc     
c                           den2 = den2 + eigth*kfac*(
c     &                          -  d_ik7(i,k) * d_jl7(j,l) ! Im(ab) and Im(ba)
c     $                          -  d_il7(i,l) * d_jk7(j,k)
c     $                          +  d_ik8(i,k) * d_jl8(j,l)
c     $                          +  d_il8(i,l) * d_jk8(j,k))
c                        endif
cc     
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         else if (orohf) then   ! ROHF
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = half *(d_ij(i,j) * d_kl(k,l))
c     $                       + quarter *(d_ij(i,j) * d_kl2(k,l)
c     $                       +  d_ij2(i,j) * d_kl(k,l))
c     $                       + eigth *(d_ij2(i,j) * d_kl2(k,l))
c     $                       - eigth *(d_ik(i,k) * d_jl(j,l)
c     $                       +  d_il(i,l) * d_jk(j,k))
c     $                       - sixteenth *(d_ik(i,k) * d_jl2(j,l)
c     $                       +  d_il2(i,l) * d_jk(j,k)
c     $                       +  d_ik2(i,k) * d_jl(j,l)
c     $                       +  d_il(i,l) * d_jk2(j,k)
c     $                       +  d_ik2(i,k) * d_jl2(j,l)
c     $                       +  d_il2(i,l) * d_jk2(j,k))
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         else if (omcscf) then  ! MCSCF
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = half *(d_ij(i,j)*d_kl(k,l))
c     $                       - eigth * (d_ik(i,k) * d_jl(j,l)
c     $                       +  d_il(i,l) * d_jk(j,k))
c     $                       + 0.5d0*pdm2d(l,k,j,i)
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         else                   ! RHF
c            do i = ilo, ihi
c               do j = jlo, jhi
c                  do k = klo, khi
c                     do l = llo, lhi
c                        den2 = half*jfac*(d_ij(i,j)*d_kl(k,l))
c     $                       - eigth*kfac*(d_ik(i,k) * d_jl(j,l)
c     $                       +  d_il(i,l) * d_jk(j,k))
c                        pdm2d(l,k,j,i) = den2
c                     end do
c                  end do
c               end do
c            end do
c         end if
c      end if
cc
c      psum = 0.0d0
c      do i = ilo, ihi
c         do j = jlo, jhi
c            do k = klo, khi
c               do l = llo, lhi
c                  psum = psum + pdm2d(l,k,j,i)**2
c               end do
c            end do
c         end do
c      end do
c      psum = sqrt(psum)
cc
c      end
