      subroutine grid_loop(dogradient,
     ,     nops,
     ,     oprint,grad,g_dens,
     ,     dosympts,sub4,largenode,grid_written,
     ,     ncoef_max,natoms,
     ,     nqtot,iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,
     ,     do_nmrks_r, tmat,
     ,     l3d,triplet,
     &     tau, ttau, kske, Mmat, Mmat2,
     ,     iqat,
     ,     qsh,aqbuf,rij,rho_int,
     D     func,d_qwght,rq0,hh,ppp,ff,zzz,
     D     drhonuc,ddelrhonuc,d2scr,force,
     D     acc_ao_gauss,
     L     do_force,do_hess,ldew,gradchi,hesschi,der3chi,oactive,
     I     iandex,max_at_bf,curatoms,g_rhs,
     I     nbf_ao_mxnbf_ce2,nactive)
     
      implicit none
#include "global.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
      logical dogradient
      integer nops
      logical dosympts,sub4,oprint,largenode
      logical grid_written
      integer iqat(*) ! [scr]
      integer g_dens(*),ncoef_max,natoms
      integer qsh(3,*)
      double precision aqbuf(4,*)
      double precision rij,rho_int
      logical grad,kske
      integer iniz(*),ifin(*),nscr,iscr(*),ibf_ao(*),
     ,     g_vxc(*),nmat
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r,
     ,      l3d,triplet
      double precision rho(*),delrho(*),hess(*),
     ,     qwght(*),qxyz(3,*),xyz(3,*),
     ,     expo(*), ccoef(*), Bmat(*), Fmat(*)
      double precision Pmat(*), rq(*),  
     R     rdens_atom(*), cetobfr(*), 
     ,     rchi_atom(*),rdelchi_atom(*),
     ,     Ex, Ec, amat(*),cmat(*),
     ,     Amat2(*), Cmat2(*),bp(*),
     ,     tmat(*), tau(*), Mmat(*), 
     &     ttau(*), Mmat2(*)
CA
      double precision func(*),d_qwght(*),
     D     rq0(*),hh(*),ppp(*),ff(*),zzz(*),force(*)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*),
     ,     acc_ao_gauss
      integer iandex(*),max_at_bf,curatoms(*),g_rhs(*)
      logical do_force,do_hess,ldew,gradchi,hesschi,der3chi,
     ,     oactive(*)
      integer nbf_ao_mxnbf_ce2,nactive
c
      integer k_symops(2),l_symops(2),
     k     k_gridscr,l_gridscr,k_iscr,l_iscr
      integer nxtask,n1,n2,icount
      external nxtask
      integer me,nproc
      integer i_a,iang_done,nwrbuf
      double precision symfact
      integer nqsh_at,iqsh,iaia,nxyz
      double precision rad
      integer nsub4,ncu,nq_out,fnxyz,nq_cube
      double precision rmax,rmin,rad_sh
      integer iii,irsh,nqtotc,ictr,nq
      logical cube_done
      integer nqtot,ncubenz,ntot_cube
      integer grid_nxyz,grid_nsub4
      external grid_nxyz,grid_nsub4
c
      nproc=ga_nnodes()
      me=ga_nodeid()
c
      ntot_cube=0
      nqtotc=0
      ncubenz=0
c
c     parallel loop init
c
      icount = 1
      n1 = 0
      n2 = nxtask(nproc,icount)

      if(dosympts) call grid_sympts0(geom,nops,k_symops,l_symops,
     ,     k_gridscr,l_gridscr,k_iscr,l_iscr,nqmax)

      do 111 i_a =1,ncenters
         iang_done=-1
         nwrbuf=1
         if(oskel) then
            if (.not.sym_atom(geom, i_a, symfact))  goto 111
         endif
c
c     get shells centers on i_a
c
         nqsh_at=0
         do iqsh=1,nqshells
            if(qsh(2,iqsh).eq.i_a) then
               nqsh_at=nqsh_at+1
               iqat(nqsh_at)=iqsh
            endif
         enddo
c
c     get no of segment
c
         iaia=iatype(i_a)
         rad=dble(rad_cutoff(1,iaia))
         nxyz=grid_nxyz(rad,lcube(iaia))
         if(sub4) then
            nsub4=grid_nsub4(nxyz)
         endif
c
c     loop over cubes centers on i_a
c
         do ncu=1,ncubes_type(iatype(i_a))
            nq_out=0
            if (n1.eq.n2 )then
               fnxyz=int_mb(k_nxyz+(iatype(i_a)-1)*ncubesmx+ncu-1)
               call grid_rmaxmin(fnxyz,rmax,rmin,lcube(iaia),nxyz)
               nq_cube=0
               do iii = 1, nqsh_at
                  iqsh=iqat(iii)
                  irsh = qsh(1,iqsh)         
                  rad_sh=rpts(irsh,iatype(i_a))
c     
c     check if radius between rmin and rmax
c     
                  if(rad_sh.ge.rmin.and.rad_sh.le.rmax) then
                     if(dogradient) then
                        call dftg_loop(
     ,                    nq_out,iqsh,irsh,rad_sh,
     ,                    dosympts,
     ,                    fnxyz,nxyz,nops,grad,
     ,                    nq_cube,nqtotc,symfact,
     ,                    iang_done,
     ,                    dbl_mb(k_gridscr),
     .                    int_mb(k_iscr),dbl_mb(k_symops(1)),
     .                    dbl_mb(k_symops(2)),
     ,                    qsh, aqbuf,rij, 
     ,                    g_dens, ncoef_max, natoms, 
     &                    rho_int,
     &                    iniz, ifin, rho, delrho, 
     &                    hess, qwght, qxyz, xyz,iscr,
     &                    expo, Bmat, Fmat, 
     &                    Pmat, rq,  rdens_atom, cetobfr, 
     ,                    ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,                    Ex, Ec, amat,cmat,g_vxc,
     ,                    Amat2, Cmat2,
     ,                    nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     &                    tau, ttau, kske, Mmat, Mmat2,
     D                    func,d_qwght,rq0,hh,ppp,ff,zzz,
     D                    drhonuc,ddelrhonuc,d2scr,force,
     D                    acc_ao_gauss,
     L                    do_force,do_hess,ldew,gradchi,hesschi,
     L                    der3chi,oactive,
     I                    iandex,max_at_bf,curatoms,g_rhs,
     I                    nbf_ao_mxnbf_ce2,nactive)
c
                        nqtot=nqtotc
                     else
c
                        call grid_quadvw1(ncu,nq_out,iqsh,irsh,rad_sh,
     ,                    dosympts,
     ,                    fnxyz,nxyz,nops,grad,
     ,                    nq_cube,nqtotc,rmax,symfact,
     ,                    iang_done,nwrbuf,
     ,                    dbl_mb(k_gridscr),
     .                    int_mb(k_iscr),dbl_mb(k_symops(1)),
     .                    dbl_mb(k_symops(2)),
     ,                    qsh, aqbuf,rij, 
     ,                    g_dens, ncoef_max, natoms, 
     &                    rho_int, nqtot, 
     &                    iniz, ifin, rho, delrho, 
     &                    hess, qwght, qxyz, xyz,iscr,nscr,
     &                    expo, ccoef, Bmat, Fmat, 
     &                    Pmat,  rq,  rdens_atom, cetobfr, 
     ,                    ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,                    Ex, Ec, amat,cmat,g_vxc,
     ,                    Amat2, Cmat2,
     ,                    nmat, do_scf, do_cpks_l, do_cpks_r,
     ,                    do_nmrks_r, tmat,
     ,                    l3d,triplet,
     &                    tau, ttau, kske, Mmat, Mmat2)
                     endif
                  endif
               enddo            ! end loop over angular shells

               if(nq_cube.ne.0) then
                  cube_done=.true.
                  ictr = qsh(2,iqsh)
                   if(dogradient) then
                   else
                      nq=nqmax
                   call grid_quadvc(cube_done,nwrbuf,
     ,                 aqbuf,nq,nq_out,
     ,                 ntot_cube,ictr,qxyz, qwght, 
     ,                 rmax,grad,
     .                 g_dens, ncoef_max, natoms, 
     &                 rho_int, nqtot, 
     &                 iniz, ifin, rho, delrho, 
     &                 hess,  xyz,iscr,nscr,
     &                 expo, ccoef, Bmat, Fmat, 
     &                 Pmat, rq, rdens_atom, cetobfr, 
     ,                 ibf_ao,rchi_atom,rdelchi_atom,
     ,                 Ex, Ec, amat,cmat,g_vxc, amat2,cmat2,
     ,                 nmat, do_scf, do_cpks_l, do_cpks_r,
     ,                 do_nmrks_r, tmat,l3d,
     ,                 triplet,
     &                 tau, ttau, kske,Mmat, Mmat2 )
                   endif
                  ncubenz=ncubenz+1
               endif
               n1 = n1 + 1
               n2 = nxtask(nproc,icount)
            else
               n1 = n1 + 1
            endif
         enddo                  ! end loop over cubes
 111  continue                  ! end loop over atoms
      n1 = nxtask(-nproc,icount)
c
      if(oprint.and.ncubenz.ne.0) write(LuOut,16) ' me ',me,
     ,     ' NQTOTC NQ/cube',nqtotc,nqtotc/ncubenz
      if(dosympts)call grid_sympts2(l_symops,l_gridscr)
 16   format(A,I3,A,2I8)
      call grid_buckend()
c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      if(largenode) return
      call grid_file_write(nq, n_tot_tasks, rad, aqbuf,0)
      grid_written=.true.
c
      return
      end
