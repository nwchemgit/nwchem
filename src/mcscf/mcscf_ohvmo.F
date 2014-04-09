       subroutine mcscf_hessv_2e_mo( nbf, nclosed, nact,
     $                               dm2, g_coul, g_exch,
     $                               g_x, g_ax )
*
* $Id$
*
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "mcscfprof.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_x                                  ! [input] Argument matrix
       integer g_ax                                 ! [output] Hessian product (in matrix format)
c
c
       integer nvir, vlen, voff, aoff, aend
       integer l_xit, k_xit, l_xua, k_xua, l_xju, k_xju
c
c
c
       if (omcscfprof) call pstat_on(ps_hv2emo)
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact

c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'xit',l_xit,k_xit))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0, MA_ERR)
       if (nclosed.gt.0) then
         call ga_get(g_x,aoff,aend,1,nclosed,dbl_mb(k_xit),nact)
       endif
       call mcscf_hessv_2e_mo_avia( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xit), g_ax )
       if (.not.ma_pop_stack(l_xit))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0,
     &       MA_ERR)
c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'xua',l_xua,k_xua))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0,
     &       MA_ERR)
       call ga_get(g_x,voff,nbf,aoff,aend,dbl_mb(k_xua),nvir)
       call mcscf_hessv_2e_mo_avav( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xua), g_ax )
       call mcscf_hessv_2e_mo_iaav( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xua), g_ax)
       if (.not.ma_pop_stack(l_xua))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0,
     &       MA_ERR)
c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'xju',l_xju,k_xju))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0, MA_ERR)
       if (nclosed.gt.0) then
         call ga_get(g_x,aoff,aend,1,nclosed,dbl_mb(k_xju),nact)
       endif
       call mcscf_hessv_2e_mo_iaia( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xju), g_ax )
       if (.not.ma_pop_stack(l_xju))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0,
     &       MA_ERR)
c
c
c
       if (omcscfprof) call pstat_off(ps_hv2emo)
       return
       end













       subroutine mcscf_hessv_2e_mo_avia( nbf, nclosed, nact, dm2,
     $                                    g_coul, g_exch, xit, g_ax )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       double precision xit(nact,nclosed)           ! [input] Local argument (Inactive-Active)
       integer g_ax                                 ! [output] Hessian product (in matrix format)
c
       integer nvir, vlen, voff, aoff, aend
       integer l_jj, k_jj, l_z, k_z, l_dt, k_dt, k_j, ld1
       integer v, x, vx, jlo, jhi, nn
       double precision fac
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       nn = nbf*nbf
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'crapj',l_jj,k_jj))
     $      call errquit('mcscf_hessv_2e_mo_avia: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'crapz',l_z,k_z))
     $      call errquit('mcscf_hessv_2e_mo_avia: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nact),'crapz',l_dt,k_dt))
     $      call errquit('mcscf_hessv_2e_mo_avia: cannot allocate',0,
     &       MA_ERR)
c
c
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avia: wrong distrib.',0,
     &       GA_ERR)
       call dfill((nact*nvir),0.d0,dbl_mb(k_z),1)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = -4.d0
           if (v.eq.x) fac = -2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call dcopy((nact*nact),dm2(1,1,v,x),1,dbl_mb(k_dt),1)
             call daxpy((nact*nact),1.d0,dm2(1,1,x,v),1,dbl_mb(k_dt),1)       ! Form symmetric 2-PDM for indices (v,x)
             call ga_access(g_coul,1,nn,vx,vx,k_j,ld1)
             call mcscf_avia_01( nbf, nclosed, nact, nvir, fac, xit,
     $                           dbl_mb(k_dt), dbl_mb(k_j),
     $                           dbl_mb(k_jj), dbl_mb(k_z) )
             call ga_release(g_coul,1,nn,vx,vx)
           endif
         enddo
       enddo
c
c
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avia: wrong distrib.',0,
     &       GA_ERR)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = -4.d0
           if (v.eq.x) fac = -2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call ga_access(g_exch,1,nn,vx,vx,k_j,ld1)
             call mcscf_avia_02( nbf, nclosed, nact, nvir, v, x, fac,
     $                           xit, dm2, dbl_mb(k_j), dbl_mb(k_jj),
     $                           dbl_mb(k_z))
             call ga_release(g_exch,1,nn,vx,vx)
           endif
         enddo
       enddo
       call ga_acc(g_ax,voff,nbf,aoff,aend,dbl_mb(k_z),nvir,1.d0)
c
c
       if (.not.ma_pop_stack(l_dt))
     $      call errquit('mcscf_hessv_2e_mo_avia: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_z))
     $      call errquit('mcscf_hessv_2e_mo_avia: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_jj))
     $      call errquit('mcscf_hessv_2e_mo_avia: pop stack failed',0,
     &       MA_ERR)
       return
       end









       subroutine mcscf_avia_01( nbf, nclosed, nact, nvir, fac,
     $                           xit, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       double precision fac
       double precision xit(nact,nclosed)
       double precision dm2(nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nvir,nact)
       double precision z(nvir,nact)
       integer voff

       voff = nclosed + nact + 1


       call dgemm( 'n', 't', nvir, nact, nclosed,
     $             1.d0, eri(voff,1), nbf, xit, nact,
     $             0.d0, tmp, nvir )
       call dgemm( 'n', 'n', nvir, nact, nact,
     $             fac, tmp, nvir, dm2, nact,
     $             1.d0, z, nvir )
       end









       subroutine mcscf_avia_02( nbf, nclosed, nact, nvir, v, x,
     $                           fac, xit, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       integer v, x
       double precision fac
       double precision xit(nact,nclosed)
       double precision dm2(nact,nact,nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nvir,nact)
       double precision z(nvir,nact)
       integer voff, u, a, t
       double precision xx

       voff = nclosed + nact + 1

       call dgemm( 'n', 't', nvir, nact, nclosed,
     $             1.d0, eri(voff,1), nbf, xit, nact,
     $             0.d0, tmp, nvir )

       do u=1,nact
         do a=1,nvir
           do t=1,nact
             xx = dm2(t,v,u,x) + dm2(t,v,x,u)
             z(a,u) = z(a,u) + tmp(a,t)*xx*fac
           enddo
         enddo
       enddo

       call dgemm( 't', 't', nvir, nact, nclosed,
     $             1.d0, eri(1,voff), nbf, xit, nact,
     $             0.d0, tmp, nvir )

       do u=1,nact
         do a=1,nvir
           do t=1,nact
             xx = dm2(t,x,u,v) + dm2(t,x,v,u)
             z(a,u) = z(a,u) + tmp(a,t)*xx*fac
           enddo
         enddo
       enddo
       return
       end









       subroutine mcscf_hessv_2e_mo_iaav( nbf, nclosed, nact, dm2,
     $                                    g_coul, g_exch, xua, g_ax )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                     ! [input] Basis functions
       integer nclosed                                 ! [input] Closed shells
       integer nact                                    ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)       ! [input] 2PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       double precision xua((nbf-nclosed-nact),nact)   ! [input] Local argument (Active-Virtual)
       integer g_ax                                    ! [output] Hessian product (in matrix format)
c
       integer nvir, vlen, voff, aoff, aend
       integer l_jj, k_jj, l_z, k_z, l_dt, k_dt, k_j, ld1
       integer v, x, vx, jlo, jhi, nn
       double precision fac
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       nn = nbf*nbf
c
c
       if (.not.ma_push_get(MT_DBL,(nclosed*nact),'crapj',l_jj,k_jj))
     $      call errquit('mcscf_hessv_2e_mo_iaav: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nclosed*nact),'crapz',l_z,k_z))
     $      call errquit('mcscf_hessv_2e_mo_iaav: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nact),'crapz',l_dt,k_dt))
     $      call errquit('mcscf_hessv_2e_mo_iaav: cannot allocate',0,
     &       MA_ERR)
c
c
c
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avia: wrong distrib.',0,
     &       GA_ERR)
       call dfill((nact*nclosed),0.d0,dbl_mb(k_z),1)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = -4.d0
           if (v.eq.x) fac = -2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call dcopy((nact*nact),dm2(1,1,v,x),1,dbl_mb(k_dt),1)
             call daxpy((nact*nact),1.d0,dm2(1,1,x,v),1,dbl_mb(k_dt),1)       ! Form symmetric 2-PDM for indices (v,x)
             call ga_access(g_coul,1,nn,vx,vx,k_j,ld1)
             call mcscf_iaav_01( nbf, nclosed, nact, nvir, fac, xua,
     $                           dbl_mb(k_dt), dbl_mb(k_j),
     $                           dbl_mb(k_jj), dbl_mb(k_z) )
             call ga_release(g_coul,1,nn,vx,vx)
           endif
         enddo
       enddo
c
c
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avia: wrong distrib.',0,
     &       GA_ERR)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = -4.d0
           if (v.eq.x) fac = -2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call ga_access(g_exch,1,nn,vx,vx,k_j,ld1)
             call mcscf_iaav_02( nbf, nclosed, nact, nvir, v, x, fac,
     $                           xua, dm2, dbl_mb(k_j), dbl_mb(k_jj),
     $                           dbl_mb(k_z))
             call ga_release(g_exch,1,nn,vx,vx)
           endif
         enddo
       enddo
       if (nclosed.gt.0) then
         call ga_acc(g_ax,aoff,aend,1,nclosed,dbl_mb(k_z),nact,1.d0)
       endif
c
c
       if (.not.ma_pop_stack(l_dt))
     $      call errquit('mcscf_hessv_2e_mo_iaav: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_z))
     $      call errquit('mcscf_hessv_2e_mo_iaav: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_jj))
     $      call errquit('mcscf_hessv_2e_mo_iaav: pop stack failed',0,
     &       MA_ERR)

       return
       end








       subroutine mcscf_iaav_01( nbf, nclosed, nact, nvir, fac,
     $                           xua, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       double precision fac
       double precision xua(nvir,nact)
       double precision dm2(nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nact,nclosed)
       double precision z(nact,nclosed)
       integer voff

       voff = nclosed + nact + 1

       call dgemm( 't', 'n', nact, nclosed, nvir,
     $             1.d0, xua, nvir, eri(voff,1), nbf,
     $             0.d0, tmp, nact )
       call dgemm( 'n', 'n', nact, nclosed, nact,
     $             fac, dm2, nact, tmp, nact,
     $             1.d0, z, nact )
       end







       subroutine mcscf_iaav_02( nbf, nclosed, nact, nvir, v, x,
     $                           fac, xua, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       integer v, x
       double precision fac
       double precision xua(nvir,nact)
       double precision dm2(nact,nact,nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nact,nclosed)
       double precision z(nact,nclosed)
       integer voff, u, i, t
       double precision xx

       voff = nclosed + nact + 1

       call dgemm( 't', 'n', nact, nclosed, nvir,
     $             1.d0, xua, nvir, eri(voff,1), nbf,
     $             0.d0, tmp, nact )

       do i=1,nclosed
         do t=1,nact
           do u=1,nact
             xx = dm2(u,x,t,v) + dm2(x,u,t,v)
             z(t,i) = z(t,i) + tmp(u,i)*xx*fac
           enddo
         enddo
       enddo

       call dgemm( 't', 't', nact, nclosed, nvir,
     $             1.d0, xua, nvir, eri(1,voff), nbf,
     $             0.d0, tmp, nact )

       do i=1,nclosed
         do t=1,nact
           do u=1,nact
             xx = dm2(u,v,t,x) + dm2(v,u,t,x)
             z(t,i) = z(t,i) + tmp(u,i)*xx*fac
           enddo
         enddo
       enddo
       return
       end











       subroutine mcscf_hessv_2e_mo_avav( nbf, nclosed, nact, dm2,
     $                                    g_coul, g_exch, xta, g_ax )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                     ! [input] Basis functions
       integer nclosed                                 ! [input] Closed shells
       integer nact                                    ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)       ! [input] 2PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       double precision xta((nbf-nclosed-nact),nact)   ! [input] Local argument (Active-Virtual)
       integer g_ax                                    ! [output] Hessian product (in matrix format)
c
       integer nvir, vlen, voff, aoff, aend
       integer l_jj, k_jj, l_z, k_z, l_dt, k_dt, k_j, ld1
       integer v, x, vx, jlo, jhi, nn
       double precision fac
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       nn = nbf*nbf
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'crapj',l_jj,k_jj))
     $      call errquit('mcscf_hessv_2e_mo_avav: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'crapz',l_z,k_z))
     $      call errquit('mcscf_hessv_2e_mo_avav: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nact),'crapz',l_dt,k_dt))
     $      call errquit('mcscf_hessv_2e_mo_avav: cannot allocate',0,
     &       MA_ERR)
c
c
c
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avav: wrong distrib.',0,
     &       GA_ERR)
       call dfill((nact*nvir),0.d0,dbl_mb(k_z),1)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = 4.d0
           if (v.eq.x) fac = 2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call dcopy((nact*nact),dm2(1,1,v,x),1,dbl_mb(k_dt),1)
             call daxpy((nact*nact),1.d0,dm2(1,1,x,v),1,dbl_mb(k_dt),1)       ! Form symmetric 2-PDM for indices (v,x)
             call ga_access(g_coul,1,nn,vx,vx,k_j,ld1)
             call mcscf_avav_01( nbf, nclosed, nact, nvir, fac, xta,
     $                           dbl_mb(k_dt), dbl_mb(k_j),
     $                           dbl_mb(k_jj), dbl_mb(k_z) )
             call ga_release(g_coul,1,nn,vx,vx)
           endif
         enddo
       enddo
c
c
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_avia: wrong distrib.',0,
     &       GA_ERR)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = 4.d0
           if (v.eq.x) fac = 2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call ga_access(g_exch,1,nn,vx,vx,k_j,ld1)
             call mcscf_avav_02( nbf, nclosed, nact, nvir, v, x, fac,
     $                           xta, dm2, dbl_mb(k_j), dbl_mb(k_jj),
     $                           dbl_mb(k_z))
             call ga_release(g_exch,1,nn,vx,vx)
           endif
         enddo
       enddo
       call ga_acc(g_ax,voff,nbf,aoff,aend,dbl_mb(k_z),nvir,1.d0)
c
c
       if (.not.ma_pop_stack(l_dt))
     $      call errquit('mcscf_hessv_2e_mo_avav: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_z))
     $      call errquit('mcscf_hessv_2e_mo_avav: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_jj))
     $      call errquit('mcscf_hessv_2e_mo_avav: pop stack failed',0,
     &       MA_ERR)

       return
       end






       subroutine mcscf_avav_01( nbf, nclosed, nact, nvir, fac,
     $                           xta, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       double precision fac
       double precision xta(nvir,nact)
       double precision dm2(nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nvir,nact)
       double precision z(nvir,nact)
       integer voff

       voff = nclosed + nact + 1

       call dgemm( 'n', 'n', nvir, nact, nvir,
     $             1.d0, eri(voff,voff), nbf, xta, nvir,
     $             0.d0, tmp, nvir )
       call dgemm( 'n', 'n', nvir, nact, nact,
     $             fac, tmp, nvir, dm2, nact,
     $             1.d0, z, nvir )
       end





       subroutine mcscf_avav_02( nbf, nclosed, nact, nvir, v, x,
     $                           fac, xta, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       integer v, x
       double precision fac
       double precision xta(nvir,nact)
       double precision dm2(nact,nact,nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nvir,nact)
       double precision z(nvir,nact)
       integer voff, u, a, t
       double precision xx

       voff = nclosed + nact + 1

       call dgemm( 'n', 'n', nvir, nact, nvir,
     $             1.d0, eri(voff,voff), nbf, xta, nvir,
     $             0.d0, tmp, nvir )

       do t=1,nact
         do a=1,nvir
           do u=1,nact
             xx = dm2(v,u,t,x) + dm2(u,v,t,x)
             z(a,t) = z(a,t) + tmp(a,u)*xx*fac
           enddo
         enddo
       enddo

       call dgemm( 't', 'n', nvir, nact, nvir,
     $             1.d0, eri(voff,voff), nbf, xta, nvir,
     $             0.d0, tmp, nvir )

       do t=1,nact
         do a=1,nvir
           do u=1,nact
             xx = dm2(x,u,t,v) + dm2(u,x,t,v)
             z(a,t) = z(a,t) + tmp(a,u)*xx*fac
           enddo
         enddo
       enddo
       return
       end











       subroutine mcscf_hessv_2e_mo_iaia( nbf, nclosed, nact, dm2,
     $                                    g_coul, g_exch, xju, g_ax )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                     ! [input] Basis functions
       integer nclosed                                 ! [input] Closed shells
       integer nact                                    ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)       ! [input] 2PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       double precision xju(nact,nclosed)              ! [input] Local argument (Active-Virtual)
       integer g_ax                                    ! [output] Hessian product (in matrix format)
c
       integer nvir, vlen, voff, aoff, aend
       integer l_jj, k_jj, l_z, k_z, l_dt, k_dt, k_j, ld1
       integer v, x, vx, jlo, jhi, nn
       double precision fac
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       nn = nbf*nbf
c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'crapj',l_jj,k_jj))
     $      call errquit('mcscf_hessv_2e_mo_iaia: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'crapz',l_z,k_z))
     $      call errquit('mcscf_hessv_2e_mo_iaia: cannot allocate',0,
     &       MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nact*nact),'crapz',l_dt,k_dt))
     $      call errquit('mcscf_hessv_2e_mo_iaia: cannot allocate',0,
     &       MA_ERR)
c
c
c
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_iaia: wrong distrib.',0,
     &       GA_ERR)
       call dfill((nact*nclosed),0.d0,dbl_mb(k_z),1)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = 4.d0
           if (v.eq.x) fac = 2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call dcopy((nact*nact),dm2(1,1,v,x),1,dbl_mb(k_dt),1)
             call daxpy((nact*nact),1.d0,dm2(1,1,x,v),1,dbl_mb(k_dt),1)       ! Form symmetric 2-PDM for indices (v,x)
             call ga_access(g_coul,1,nn,vx,vx,k_j,ld1)
             call mcscf_iaia_01( nbf, nclosed, nact, nvir, fac, xju,
     $                           dbl_mb(k_dt), dbl_mb(k_j),
     $                           dbl_mb(k_jj), dbl_mb(k_z) )
             call ga_release(g_coul,1,nn,vx,vx)
           endif
         enddo
       enddo
c
c
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_hessv_2e_mo_iaia: wrong distrib.',0,
     &       GA_ERR)
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           fac = 4.d0
           if (v.eq.x) fac = 2.d0
           if ((vx.ge.jlo).and.(vx.le.jhi)) then
             call ga_access(g_exch,1,nn,vx,vx,k_j,ld1)
             call mcscf_iaia_02( nbf, nclosed, nact, nvir, v, x, fac,
     $                           xju, dm2, dbl_mb(k_j), dbl_mb(k_jj),
     $                           dbl_mb(k_z))
             call ga_release(g_exch,1,nn,vx,vx)
           endif
         enddo
       enddo
       if (nclosed.gt.0) then
         call ga_acc(g_ax,aoff,aend,1,nclosed,dbl_mb(k_z),nact,1.d0)
       endif
c
c
c
       if (.not.ma_pop_stack(l_dt))
     $      call errquit('mcscf_hessv_2e_mo_iaia: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_z))
     $      call errquit('mcscf_hessv_2e_mo_iaia: pop stack failed',0,
     &       MA_ERR)
       if (.not.ma_pop_stack(l_jj))
     $      call errquit('mcscf_hessv_2e_mo_iaia: pop stack failed',0,
     &       MA_ERR)

       return
       end











       subroutine mcscf_iaia_01( nbf, nclosed, nact, nvir, fac,
     $                           xju, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       double precision fac
       double precision xju(nact,nclosed)
       double precision dm2(nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nact,nclosed)
       double precision z(nact,nclosed)
       integer voff

       voff = nclosed + nact + 1

       call dgemm( 'n', 'n', nact, nclosed, nclosed,
     $             1.d0, xju, nact, eri(1,1), nbf,
     $             0.d0, tmp, nact )
       call dgemm( 'n', 'n', nact, nclosed, nact,
     $             fac, dm2, nact, tmp, nact,
     $             1.d0, z, nact )
       end











       subroutine mcscf_iaia_02( nbf, nclosed, nact, nvir, v, x,
     $                           fac, xju, dm2, eri, tmp, z )
       implicit none
       integer nbf, nclosed, nact, nvir
       integer v, x
       double precision fac
       double precision xju(nact,nclosed)
       double precision dm2(nact,nact,nact,nact)
       double precision eri(nbf,nbf)
       double precision tmp(nact,nclosed)
       double precision z(nact,nclosed)
       integer voff, i, u, t
       double precision xx

       voff = nclosed + nact + 1

       call dgemm( 'n', 'n', nact, nclosed, nclosed,
     $             1.d0, xju, nact, eri(1,1), nbf,
     $             0.d0, tmp, nact )

       do i=1,nclosed
         do t=1,nact
           do u=1,nact
             xx = dm2(v,t,u,x) + dm2(t,v,u,x)
             z(t,i) = z(t,i) + tmp(u,i)*xx*fac
           enddo
         enddo
       enddo

       call dgemm( 'n', 't', nact, nclosed, nclosed,
     $             1.d0, xju, nact, eri(1,1), nbf,
     $             0.d0, tmp, nact )

       do i=1,nclosed
         do t=1,nact
           do u=1,nact
             xx = dm2(x,t,u,v) + dm2(t,x,u,v)
             z(t,i) = z(t,i) + tmp(u,i)*xx*fac
           enddo
         enddo
       enddo

       return
       end
