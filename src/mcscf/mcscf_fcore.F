c
c        
c  Form
c
c       F    = h                                core Hamiltonian
c        ij     ij
c
c                              rs    rs
c       F    =  h   +  sum ( 2J   - K   )       inactive Fock element
c        rs      rs       i    ii    ii
c
c
       subroutine mcscf_fcore( basis, nbf, nclosed, nact, osym,
     $                         g_movecs, g_coul, g_exch, g_fcore )
       implicit none
#include "mafdecls.fh"       
#include "global.fh"
       integer basis
       integer nbf
       integer nclosed
       integer nact
       integer osym(nbf)
       integer g_movecs
       integer g_coul
       integer g_exch
       integer g_fcore
       integer nn, t, u, tu, jlo, jhi, k_j, k_k, ld1
       logical ga_check_JKblocked
       external ga_check_JKblocked
       DOUBLE PRECISION CRAP(200)
c
c  1e Hamiltononian
c
       call moints_1e( nbf, basis, g_movecs, g_fcore)
       call moints_1esym( nbf, g_fcore, osym )
c
c$$$       CALL GA_GET(G_FCORE,1,NBF,1,NBF,CRAP,NBF)
c$$$       PRINT*,'<<<<<<< 1e Hamiltonian >>>>>>>>'
c$$$       CALL MAT_PRINT(NBF,NBF,CRAP)
c
c  Inactive-active Coulomb interation
c
       nn = nbf*nbf
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_fcore: wrong distrib. for Coulomb',0)
       do t=1,nact
         do u=1,t
           tu = ((t-1)*t)/2 + u
           if ((tu.ge.jlo).and.(tu.le.jhi)) then
             call ga_access(g_coul,1,nn,tu,tu,k_j,ld1)
             call mcscf_fcore01( nbf, nclosed, nact, t, u, 2.d0,
     $                           dbl_mb(k_j), g_fcore)
             call ga_release(g_coul,1,nn,tu,tu)
           endif
         enddo
       enddo
c
c  Inactive-active Exchange interaction
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_fcore: wrong distrib. for exchange',0)
       do t=1,nact
         do u=1,t
           tu = ((t-1)*t)/2 + u
           if ((tu.ge.jlo).and.(tu.le.jhi)) then
             call ga_access(g_exch,1,nn,tu,tu,k_k,ld1)
             call mcscf_fcore01( nbf, nclosed, nact, t, u, -1.d0,
     $                           dbl_mb(k_k), g_fcore)
             call ga_release(g_exch,1,nn,tu,tu)
           endif
         enddo
       enddo
c
c
c
       return
       end






       subroutine mcscf_fcore01( nbf, nclosed, nact, t, u, fact,
     $                           xmo, g_fcore )
       implicit none
#include "global.fh"
       integer nbf
       integer nclosed
       integer nact
       integer t, u
       double precision fact
       double precision xmo(nbf,nbf)
       integer g_fcore
       double precision xx
       integer i, tt, uu

       xx = 0.d0
       do i=1,nclosed
         xx = xx + xmo(i,i)
       enddo
       tt = nclosed + t
       uu = nclosed + u
       call ga_acc(g_fcore, tt, tt, uu, uu, xx, 1, fact )
       if (t.ne.u) call ga_acc(g_fcore, uu, uu, tt, tt, xx, 1, fact )
       return
       end

         
           
           



       subroutine moints_1esym( nbf, g_a, osym )
       implicit none
#include "global.fh"
       integer nbf
       integer g_a
       integer osym(nbf)
       integer i,j,jj,isym
       double precision xx,yy
c
       isym = 0
       do i=1,nbf
         yy = 0.d0
         do j=1,i-1
           call ga_get(g_a,i,i,j,j,xx,1)
           if (abs(xx).gt.yy) then
             yy = abs(xx)
             jj = j
           endif
         enddo
         if (yy.gt.1.d-8) then
           osym(i) = osym(jj)
         else
           isym = isym + 1
           osym(i) = isym
         endif
       enddo
       return
       end



           
