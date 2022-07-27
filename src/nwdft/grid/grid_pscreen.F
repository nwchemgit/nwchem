C> \ingroup nwdft_xc
C> @{
C>
C> \file grid_pscreen.F
C> Screen the grid points
C>
C> @}
C>
C> \ingroup nwdft_x
C> @{
C>
C> \brief Screen the grid points and associated quantities
C>
C> Estimate the maximum contribution that grid points may make
C> based on the amplitudes of the basis functions across all atoms,
C> the maximum density matrix element on an atom, and the weight of
C> each grid point. If a grid point contributes too little it is 
C> eliminated. 
C>
C> Returned are:
C> * The number of remaining grid points
C> * The compressed list of grid point coordinates and weights
C> * The compressed list of basis function values
C> * The compressed list of basis function gradients
C>
      subroutine grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,chi,delchi,
     ,     heschi,
     ,     mbf_ao,grad,dolap,ipol,iscratch,adft)
c
c$Id$
c
      implicit none
#include "mafdecls.fh"
      logical adft
      integer natoms !< [Input] The number of atoms
      integer nqpts !< [In/output] The number of grid points
      double precision qwght(*) !< [In/output] The grid point weights
      double precision qxyz(3,*)!< [In/output] The grid point coordinates
      double precision rchi_atom(*) !< [Input] The overall basis function
                                    !< amplitude for each atom
      double precision rdelchi_atom(*) !< [Input] The overall basis
                                       !< function gradient amplitude
                                       !< for each atom
      double precision chi(*)   !< [Input] The basis function values on
                                !< the grid
      double precision delchi(*) !< [Input] The basis function gradient
                                 !< values on the grid
      double precision heschi(*) !< [Input] The basis function hessian
                                 !< values on the grid
      double precision rdens_atom(natoms,natoms,*) !< [Input] The 
      !< density matrix size for each pair of atoms and spin
      integer iniz(*) !< [Input] The initial basis function for each
                      !< atom
      integer mbf_ao  !< [Input] The number of basis functions
      integer iscratch(*) !< [Scratch]
      double precision tol_rho !< [Input] The tolerance on the density
      logical grad !< [Input] If `.true.` the basis function gradients
                   !< have to be handled in addition to the basis
                   !< functions
      logical dolap !< [Input] If `.true.` the basis function hessian
                    !< have to be handled
      integer ipol !< [Input] The number of spin channels present
c
      
      double precision dabsmax
      external dabsmax
c
      integer nzero,i,j,ispin
      double precision pmax,funcmax,ff,delmax
c
      nzero=0
c
c     get pmax ang funcmax
c
      pmax=0d0
      funcmax=dabsmax(natoms,rchi_atom)
      delmax=funcmax
      if(grad) delmax=max(delmax,dabsmax(natoms,rdelchi_atom))
      do ispin=1,ipol
         do i=1,natoms
            if(iniz(i).ne.0) 
     .           pmax=max(pmax,dabsmax(natoms,rdens_atom(1,i,ispin)))
         enddo
      enddo
      ff=abs(funcmax*delmax*pmax)
      if (dolap) ff=max(ff,abs(2d0*pmax*delmax**2))
      if (adft) ff=delmax*pmax
      do i=1,nqpts
         if(abs(ff*qwght(i)).lt.tol_rho) then
            nzero=nzero+1
            qwght(i)=0d0
         endif
      enddo
      if(nzero.ne.0) then
         call grid_pcompress(nqpts,qwght,chi,qxyz,delchi,
     H        heschi,
     ,        mbf_ao,1d-19,grad,dolap,iscratch)
      endif
      return
      end
      subroutine grid_pcompress(nq,qwght,chi,qxyz,delchi,
     H     heschi,
     .     mbf_ao,toll,grad,dolap,iscr)
      implicit none
      integer nq,mbf_ao
      double precision qwght(*),chi(*),qxyz(3,*),delchi(*),
     H heschi(*)
      logical grad,dolap
      double precision toll ![input]
c
      integer q,nzero,mq0,mq1,j,iscr(*),
     ,     nq_orig,llq
#if defined(DEBUG)
      integer nqxmbf
      nqxmbf=nq*mbf_ao
#endif
c        
c       compress grid
c
      if(nq.eq.0 )return
      nq_orig=nq
c
c     BFs
c
      nzero=0
      do q=1,nq
         if(abs(qwght(q)).ge.toll) then
            nzero=nzero+1
            iscr(nzero)=q
         endif
      enddo

      if(nzero.eq.nq) return
      nq=nzero
      mq0=0
      mq1=0
      do j=1,mbf_ao
         do q=1,nq
            chi(mq1+q)=chi(mq0+iscr(q))
         enddo
         mq0=mq0+nq_orig
         mq1=mq1+nq
      enddo
      if(grad) then
         mq0=0
         mq1=0
         do j=1,mbf_ao*3
            do q=1,nq
               delchi(mq1+q)=delchi(mq0+iscr(q))
            enddo
            mq0=mq0+nq_orig
            mq1=mq1+nq
         enddo
      endif
      if(dolap) then
         mq0=0
         mq1=0
         do j=1,mbf_ao*6
            do q=1,nq
               heschi(mq1+q)=heschi(mq0+iscr(q))
            enddo
            mq0=mq0+nq_orig
            mq1=mq1+nq
         enddo
      endif
c
c     grid pts
c
      do q=1,nq
         qwght(q)=qwght(iscr(q))
         qxyz(1,q)=qxyz(1,iscr(q))
         qxyz(2,q)=qxyz(2,iscr(q))
         qxyz(3,q)=qxyz(3,iscr(q))
      enddo
#if defined(DEBUG)
      if(nqxmbf/mbf_ao.ne.nq) then
      write(0,*) ' NONO: bf ',nqxmbf/mbf_ao,' ww ',nq
      endif
#endif
      return
      end
C>
C> @}
