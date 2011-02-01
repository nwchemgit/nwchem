      subroutine detci_onepdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                         nexa, nexb, osym, exa, exb, map,
     $                         g_civec, g_civect, onepdm )
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer nexa                            ! [input] Alpha excitations
      integer nexb                            ! [input] Beta excitations
      integer osym(norb)                      ! [input] Orbital irreps     
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer map(norb)                       ! [input] Orbital index map
      integer g_civec                         ! [input] CI-vector
      integer g_civect                        ! [input] CI-vector transpose
      double precision onepdm(norb,norb)      ! [output] 1-pdm 
c
      integer ia, ja, ib, jb
      integer myid, rlo, rhi, cilo, cihi, dilo, dihi
      integer l_t, k_t, k_s, lds
      integer g_dtmp
      integer ii, jj, iex, ph
      double precision xx
c
c
c
      myid = ga_nodeid()
*ga:1:0
      if (.not.(ga_create(MT_DBL, norb, norb, 'd', norb, 0, g_dtmp )))
     $   call errquit('detci_twopdm: cannot create global',0, GA_ERR)
      call ga_zero(g_dtmp)
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_onepdm: wrong distrib for CI vector',0,
     &       INPUT_ERR)
      call ga_distribution( g_civect, myid, rlo, rhi, dilo, dihi )
      if (((dilo.ne.0).and.(dihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstra)))
     $  call errquit('detci_onepdm: wrong distrib for trp(CI) vec',0,
     &       INPUT_ERR)
c
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      do ia=1,nstra
        call ga_get( g_civec, 1, nstrb, ia, ia, dbl_mb(k_t), 1 )
        do iex=1,nexa
          ja = exa(1,iex,ia)
          if ((ja.ge.cilo).and.(ja.le.cihi)) then
            ph = exa(4,iex,ia)
            ii = map(exa(5,iex,ia))
            jj = map(exa(6,iex,ia))
            call ga_access(g_civec, 1, nstrb, ja, ja, k_s, lds )
            xx = ddot(nstrb, dbl_mb(k_t), 1, dbl_mb(k_s), 1) * ph
            call ga_release( g_civec, 1, nstrb, ja, ja )
            call ga_acc( g_dtmp, ii, ii, jj, jj, xx, 1, 1.d0 )
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
      call ga_sync()

c
c
c
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      do ib=1,nstrb
        call ga_get( g_civect, 1, nstra, ib, ib, dbl_mb(k_t), 1 )
        do iex=1,nexb
          jb = exb(1,iex,ib)
          if ((jb.ge.dilo).and.(jb.le.dihi)) then
            ph = exb(4,iex,ib)
            ii = map(exb(5,iex,ib))
            jj = map(exb(6,iex,ib))
            call ga_access( g_civect, 1, nstra, jb, jb, k_s, lds )
            xx = ddot( nstra, dbl_mb(k_t), 1, dbl_mb(k_s), 1 ) * ph
            call ga_release( g_civect, 1, nstra, jb, jb )
            call ga_acc( g_dtmp, ii, ii, jj, jj, xx, 1, 1.d0 )
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
      call ga_sync()
      call ga_get(g_dtmp, 1, norb, 1, norb, onepdm, norb )
c
c
c      
      if (.not.ga_destroy(g_dtmp))
     $  call errquit('detci_onepdm:cannot destroy tmp',0, GA_ERR)
      return
      end














c
c
c  Two-particle density matrix
c  Some contributions in following routine
c
c
c              1
c  Gamma     = - <0| E  E   - delta  E   |0>
c       ijkl   2      ij kl        jk il
c
c
      subroutine detci_twopdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                       nexa, nexb, osym, exa, exb, map,
     $                       g_civec, g_civect, twopdm )
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "detci.fh"
      integer norb                                    ! [input] Orbitals
      integer nsym                                    ! [input] Irreps
      integer nela                                    ! [input] Alpha electrons
      integer nelb                                    ! [input] Beta electrons
      integer nstra                                   ! [input] Alpha strings
      integer nstrb                                   ! [input] Beta strings
      integer nexa                                    ! [input] Alpha excitations
      integer nexb                                    ! [input] Beta excitations
      integer osym(norb)                              ! [input] Orbital irreps     
      integer exa(6,nexa,nstra)                       ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)                       ! [input] Beta excitation lookup table
      integer map(norb)                               ! [input] Orbital index map
      integer g_civec                                 ! [input] CI-vector
      integer g_civect                                ! [input] CI-vector transposed
      double precision twopdm(norb,norb,norb,norb)    ! [output] 2-pdm 
c
c     
      integer ia, ib, ja, jb, ka, kb
      integer iex, jex, ij, kl, jl
      integer i, j, l, k
      double precision xx, yy, p1, p2
      integer l_t, k_t, l_s, k_s, lds
**      integer g_dentmp
      integer myid, nn, rlo, rhi, cilo, cihi, dilo, dihi
c
c
      myid = ga_nodeid()
      nn = norb*norb
*ga:1:0
*      if (.not.(ga_create(MT_DBL, nn, nn, 'g dens', nn, 0, g_dentmp )))
*     $   call errquit('detci_twopdm: cannot create global',0)
*      call ga_zero(g_dentmp)
      call ga_distribution( g_civect, myid, rlo, rhi, dilo, dihi )
      if (((dilo.ne.0).and.(dihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstra)))
     $  call errquit('detci_sigmaaa: wrong distrib for trp(CI)',0,
     &       INPUT_ERR)
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_sigmaaa: wrong distrib for CI vector',0,
     &       INPUT_ERR)
c
c
c
c   a  a
c  E  E      contribution
c   ij kl
c
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0, MA_ERR)
c
      if (cilo.gt.0 .and. cihi.gt.0) 
     $     call ga_access(g_civec, 1, nstrb, cilo, cihi, k_s, lds )
      call ga_sync()
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          p1 = exa(4,iex,ia)
          j  = map(exa(5,iex,ia))
          i  = map(exa(6,iex,ia))
          ij = (j-1)*norb + i
          do jex=1,nexa
            ka = exa(1,jex,ja)
            if ((ka.ge.cilo).and.(ka.le.cihi)) then
              p2 = exa(4,jex,ja) * p1 * 0.5d0
              l  = map(exa(5,jex,ja))
              k  = map(exa(6,jex,ja))
              kl = (l-1)*norb + k
*              call ga_access(g_civec, 1, nstrb, ka, ka, k_s, lds )
              xx = ddot(nstrb, dbl_mb(k_s+(ka-cilo)*lds), 1, 
     $             dbl_mb(k_t), 1 ) * p2
*              call ga_release( g_civec, 1, nstrb, ka, ka )
*              call ga_acc( g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
            endif
          enddo
        enddo
      enddo
      if (cilo.gt.0 .and. cihi.gt.0) 
     $     call ga_release(g_civec, 1, nstrb, cilo, cihi)
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
c
c   b  b
c  E  E      contribution
c   ij kl
c
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      if (dilo.gt.0 .and. dihi.gt.0)
     $     call ga_access(g_civect, 1, nstra, dilo, dihi, k_s, lds )
      call ga_sync()
      do ib=1,nstrb                                                  ! Parallel
        call ga_get(g_civect,1,nstra,ib,ib,dbl_mb(k_t),1)
        do iex=1,nexb
          jb = exb(1,iex,ib)
          p1 = exb(4,iex,ib)
          j  = map(exb(5,iex,ib))
          i  = map(exb(6,iex,ib))
          ij = (j-1)*norb + i
          do jex=1,nexb
            kb = exb(1,jex,jb)
            if ((kb.ge.dilo).and.(kb.le.dihi)) then
              p2 = exb(4,jex,jb) * p1 * 0.5d0
              l  = map(exb(5,jex,jb))
              k  = map(exb(6,jex,jb))
              kl = (l-1)*norb + k
*              call ga_access(g_civect, 1, nstra, kb, kb, k_s, lds )
              xx = ddot(nstra, dbl_mb(k_s+(kb-dilo)*lds), 1,
     $             dbl_mb(k_t), 1 ) * p2
*              call ga_release( g_civect, 1, nstra, kb, kb )
*              call ga_acc( g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
            endif
          enddo
        enddo
      enddo
      if (dilo.gt.0 .and. dihi.gt.0)
     $     call ga_release(g_civect, 1, nstra, dilo, dihi)
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
c
c   a  
c  E  delta      contribution
c   il     jk
c
c
      call ga_sync()
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          if ((ja.ge.cilo).and.(ja.le.cihi)) then
            p1 = exa(4,iex,ia)
            i  = map(exa(5,iex,ia))
            l  = map(exa(6,iex,ia))
            call ga_access(g_civec, 1, nstrb, ja, ja, k_s, lds )
            xx = ddot( nstrb, dbl_mb(k_t), 1, dbl_mb(k_s), 1 )
            call ga_release(g_civec, 1, nstrb, ja, ja)
            yy = -0.5d0*xx*p1
            do j=1,norb
              ij = (map(j)-1)*norb + i
              jl = (l-1)*norb + map(j)
*              call ga_acc(g_dentmp, ij, ij, jl, jl, yy, 1, 1.d0 )
              twopdm(i,map(j),map(j),l) = twopdm(i,map(j),map(j),l) + yy
            enddo
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
c
c   b
c  E  delta      contribution
c   il     jk
c
c
      call ga_sync()
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0, MA_ERR)
      if  (dilo.gt.0 .and. dihi.gt.0) 
     $     call ga_access(g_civect, 1, nstra, dilo, dihi, k_s, lds )
      do ib=1,nstrb
        call ga_get(g_civect, 1, nstra, ib, ib, dbl_mb(k_t), 1 )
        do iex=1,nexb
          jb = exb(1,iex,ib)
          if ((jb.ge.dilo).and.(jb.le.dihi)) then
            p1 = exb(4,iex,ib)
            i  = map(exb(5,iex,ib))
            l  = map(exb(6,iex,ib))
*            call ga_access(g_civect, 1, nstra, jb, jb, k_s, lds )
            xx = ddot(nstra, dbl_mb(k_t), 1,
     $           dbl_mb(k_s+(jb-dilo)*lds), 1 )
*            call ga_release( g_civect, 1, nstra, jb, jb )
            yy = -0.5d0*xx*p1
            do j=1,norb
              ij = (map(j)-1)*norb + i
              jl = (l-1)*norb + map(j)
*              call ga_acc(g_dentmp, ij, ij, jl, jl, yy, 1, 1.d0 )
              twopdm(i,map(j),map(j),l) = twopdm(i,map(j),map(j),l) + yy
            enddo
          endif
        enddo
      enddo
      if  (dilo.gt.0 .and. dihi.gt.0) 
     $     call ga_release(g_civect, 1, nstra, dilo, dihi)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0, MA_ERR)
      call ga_sync()
*      call ga_get(g_dentmp,1,nn,1,nn,twopdm,nn)
c
c  Clean up
c
*      if (.not.ga_destroy(g_dentmp))
*     $  call errquit('detci_twopdm: cannot destroy global dens',0)
          
      return
      end

















c                 a  b         b  a                 
c  Special for   E  E    and  E  E   contributions
c                 ij kl        ij kl              
c
c  Analogous to sigma alpha-beta routine
c  
c
      subroutine detci_twopdm_ab( norb, nsym, nela, nelb, nstra, nstrb,
     $                            nexa, nexb, nekla, neklb,
     $                            osym, map, exa, exb, ataba, atabb, 
     $                            vrhs, vlhs, vphase,
     $                            cprime, sprime, g_civec, twopdm )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb                                     ! [input] Orbitals                             
      integer nsym                                     ! [input] Irreps                               
      integer nela                                     ! [input] Alpha electrons                      
      integer nelb                                     ! [input] Beta electrons                       
      integer nstra                                    ! [input] Alpha strings                        
      integer nstrb                                    ! [input] Beta strings                         
      integer nexa                                     ! [input] Alpha excitations                    
      integer nexb                                     ! [input] Beta excitations                     
      integer nekla                                    ! [input] Max non-zero alpha strings for E_kl  
      integer neklb                                    ! [input] Max non-zero beta strings for E_kl   
      integer map(norb)                                ! [input] External remapping of indices        
      integer osym(norb)                               ! [input] Orbital irreps                       
      integer ataba(0:norb,0:nela,nsym)                    ! [input] Alpha arc weights                    
      integer atabb(0:norb,0:nelb,nsym)                    ! [input] Beta arc weights                     
      integer exa(6,nexa,nstra)                        ! [input] Alpha excitation lookup table        
      integer exb(6,nexb,nstrb)                        ! [input] Beta excitation lookup table         
      integer vrhs(nekla)                              ! [scratch] Array of RHS Strings for E_kl      
      integer vlhs(nekla)                              ! [scratch] Array of LHS Strings for E_kl      
      integer vphase(nekla)                            ! [scratch] Array of Phases for E_kl           
      double precision cprime(nstrb,nekla)             ! [scratch] Gathered CI-vector                 
      double precision sprime(nstrb,nekla)             ! [scratch] Gathered sigma vector              
      integer g_civec                                  ! [input] CI-vector                            
      double precision twopdm(norb,norb,norb,norb)     ! [output] Two-particle density 
c
      integer i, k, l
      integer ii, jj, kk, ll
      integer ij, ji, kl, lk
      integer ai, aj
      integer iph
      integer jstr, kstr
      integer ne_kl, nidx, iex
      integer relv(detci_maxelec)
      integer lelv(detci_maxelec)
      integer oidx(detci_maxorb)
      integer ip(detci_maxelec)
*      integer g_dentmp
      integer ploop, numnodes, next, myid, nn
      double precision phase, xx, dummy(1)
      double precision tx
      integer nxtask
      external nxtask
c      
c
c   Initialize parallel stuff
c
      tx = util_cpusec()
      nn = norb*norb
*ga:1:0
*      if (.not.(ga_create(MT_DBL, nn, nn, 'g dens', nn, 0, g_dentmp )))
*     $   call errquit('detci_twopdm_ab: cannot create global',0)
*      call ga_zero(g_dentmp)
*      if (ga_nodeid().eq.0)
*     $   call ga_put( g_dentmp, 1, nn, 1, nn, twopdm, nn )
      call ga_sync()
      ploop = -1
      numnodes = ga_nnodes()
      myid = ga_nodeid()
      next = nxtask(numnodes, 1)
c
c   Parallel loop over all alpha excitation operators
c
c            t
c      E  = a a
c       ij   i j
c
      do ai=1,norb
        do aj=1,norb
          ploop = ploop + 1
          if (ploop.eq.next) then
            ij = (map(aj)-1)*norb + map(ai)
            ji = (map(ai)-1)*norb + map(aj)
            ii = map(ai)
            jj = map(aj)
            ne_kl = 0
c
c   Vector of orbital indices except create/annih indices
c   Initialize pointer vector
c
            nidx = 0
            do i=1,norb
              if ((i.ne.ai).and.(i.ne.aj)) then
                nidx = nidx + 1
                oidx(nidx) = i
              endif
            enddo
            do i=1,nela-1
              ip(i) = i
            enddo
c
c   Loop through all strings for nidx and (nela-1)
c   Insert orbital index k and l to create 
c   LHS and RHS strings where
c
c        |LHS> = E  |RHS>
c                 kl
c
c   Put indices into gather/scatter arrays
c   Note: special case when nela = norb then
c
c       E  !RHS> = 0   for k != l  
c        kl
c
c   thus ne_kl = 0
c     
            ne_kl = 0
            if (nela.eq.1) then
              ne_kl = 1
              vrhs(i) = aj
              vlhs(i) = ai
              vphase(i) = 1
            else if ((norb.ne.nela).or.(ai.eq.aj)) then
 101          continue
              ne_kl = ne_kl + 1
              iph = 1
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           aj,relv,iph)
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           ai,lelv,iph)
              vphase(ne_kl) = iph
              vrhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,
     $                                    ataba,relv)
              vlhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,
     $                                    ataba,lelv)
              if (detci_getnextelv(nidx,(nela-1),ip)) goto 101
            endif
c
c   Gather in left and right CI columns
c
            call detci_cigather( nstrb, nstra, ne_kl, g_civec, vlhs,
     $                           vphase, cprime , dummy, .false.)
            call detci_cigather0( nstrb, nstra, ne_kl, g_civec, vrhs,
     $                            sprime )
c
c   Loop over all beta strings
c
            do jstr=1,nstrb
              do iex=1,nexb
                kstr  = exb(1,iex,jstr)
                phase = exb(4,iex,jstr)
                l     = map(exb(5,iex,jstr))
                k     = map(exb(6,iex,jstr))
                kk    = k
                ll    = l
                kl    = (l-1)*norb + k
                lk    = (k-1)*norb + l
                xx = ddot(ne_kl, cprime(jstr,1), nstrb,
     $                    sprime(kstr,1), nstrb ) * phase * 0.5d0
*                call ga_acc(g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
*                call ga_acc(g_dentmp, lk, lk, ji, ji, xx, 1, 1.d0 )
                 twopdm(ii,jj,kk,ll) = twopdm(ii,jj,kk,ll) + xx
		 twopdm(ll,kk,jj,ii) = twopdm(ll,kk,jj,ii) + xx
              enddo
            enddo
c
c   End parallel task
c
            next = nxtask(numnodes, 1)
          endif
        enddo
      enddo
      next = nxtask(-numnodes, 1)
c
c   Cleanup
c
*      call ga_get(g_dentmp, 1, nn, 1, nn, twopdm, nn )
*      if (.not.ga_destroy(g_dentmp))
*     $  call errquit('detci_twopdm_ab: cannot destroy global dens',0)
*      return
      end
