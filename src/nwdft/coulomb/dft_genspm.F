      subroutine dft_genspm(lmax,iga_dens,basis,
     &     natoms,xyz,Pmat,max_at_bf,
     &     charge,dipole)
C$Id: dft_genspm.F,v 1.1 1995-05-03 19:27:50 og845 Exp $
      implicit none
#include <mafdecls.h>
#include <tcgmsg.fh>
#include <global.fh>
#include <bas.fh>
#include "stdio.fh"
      integer int_nint,nxtask
      external int_nint,nxtask
      double precision ddot
      external ddot

      integer natoms
      integer iga_dens,basis,lmax
      double precision xyz(3,natoms)
      double precision Pmat(max_at_bf*max_at_bf)
      integer max_at_bf

      double precision charge(natoms)
      double precision dipole(3,natoms)
      integer lscr,iscr,iMP,lMP
      integer lmpint              

      integer ifirst,ilast,iat,idim,ish,ifirst_sh,ilast_sh,nbf_i,jfunc
      integer jfirst,jlast,jat,jdim,jsh,jfirst_sh,jlast_sh,nbf_j
      logical status
      integer nshells,max_scratch,maxg
      integer num_int,lval,ij_int
      integer istart,ipmat,Pmat_point,Pmat_point0,istartd
      integer n1,n2,icount,nproc

c     common /cmm_block4/ xcrd, ycrd, zcrd
c     common /cmm_block6/ acharge
c     double precision xcrd(NMAX), ycrd(NMAX), zcrd(NMAX), acharge(NMAX)
c     integer NMAX, MMAX, NCELLMAX
c     parameter (NMAX=100000,MMAX=2**MAXFIN, NCELLMAX=MMAX*MMAX*MMAX)
c     common / block4 / xcrd, ycrd, zcrd
c     common / block6 / xp, yp, zp
c     double precision xcrd(NMAX), ycrd(NMAX), zcrd(NMAX)
c     double precision xp(NMAX), yp(NMAX), zp(NMAX) 

c     
c     evaluation of # of multipole ints
c     
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit('Exiting in dft_gensmp.',3)

      lmpint=0
      do ish=1,nshells
        do jsh=1,ish
          ij_int = int_nint(basis, ish, basis, jsh, 0,0, 0,0)
          num_int=0
          do lval = 0,lmax
            num_int = num_int+ 1 + ij_int*(lval+1)*(lval+2)/2
          enddo
          lmpint=max(lmpint,num_int)
        enddo
      enddo

      nproc  = ga_nnodes()

      if(.not.MA_Push_Get(mt_dbl,lmpint,' multipole ints',lMP,iMP))
     &     call errquit(' exiting in dft_gensmp: insuff stack',1)

      call int_mem_2e3c(maxg, max_scratch)
      if(.not.MA_Push_Get(mt_dbl,max_scratch,' scratch',
     &     lscr,iscr))
     &     call errquit(' exiting in dft_gensmp: insuff stack',2)


      call dfill(natoms,0.d0,charge,1)
      if(lmax.gt.0) then
      call dfill(natoms*3,0.d0,dipole,1)
      endif

      icount=1
      n1 = nxtask(nproc,icount)
      n2 = 0

      do iat=1,natoms
        if (.not.bas_ce2cnr(basis,iat,ifirst,ilast))
     &       call errquit('Exiting in dft_gensmp',1)

        do jat=1,natoms
          if(n1.eq.n2) then
          if (.not.bas_ce2cnr(basis,jat,jfirst,jlast))
     &         call errquit('Exiting in dft_gensmp',2)

          call get_atom_block(iga_dens, basis,
     $         iat, jat, Pmat, idim, jdim)
          Pmat_point0=1
          do ish=ifirst,ilast
            Pmat_point=Pmat_point0
            
            if (.not. bas_cn2bfr( basis,ish,
     &           ifirst_sh,ilast_sh))
     &           call errquit('Exiting in dft_gensmp.',5)
            nbf_i=ilast_sh-ifirst_sh+1

            do jsh=jfirst,jlast
              
              if (.not. bas_cn2bfr( basis,jsh,
     &             jfirst_sh,jlast_sh))
     &             call errquit('Exiting in dft_gensmp.',5)
              nbf_j=jlast_sh-jfirst_sh+1

              call int_mpole(basis, jsh, basis, ish, lmax, xyz(1,iat),
     &             max_scratch, DBL_MB(iscr), lmpint, DBL_MB(iMP))
              
c             
c             charges
c             

              istart=iMP
              istartd=iMP+nbf_i*nbf_j
              do jfunc=1,nbf_j

               charge(iat)=charge(iat) + 
     &               ddot(nbf_i,Pmat(Pmat_point),1,DBL_MB(istart),1)
                istart=istart+nbf_i

                if(lmax.gt.0)then
c                 
c                 dipoles
c                 

                  dipole(1,iat) = dipole(1,iat) +
     +                 ddot(nbf_i,Pmat(Pmat_point),1,DBL_MB(istartd),1)      
                  istartd=istartd+nbf_i
                  dipole(2,iat) = dipole(2,iat) +
     +                 ddot(nbf_i,Pmat(Pmat_point),1,DBL_MB(istartd),1)      
                  istartd=istartd+nbf_i
                  dipole(3,iat) = dipole(3,iat) +
     +                 ddot(nbf_i,Pmat(Pmat_point),1,DBL_MB(istartd),1)      
                  istartd=istartd+nbf_i

                endif

                Pmat_point=Pmat_point+idim

              enddo

            enddo
            Pmat_point0=Pmat_point0+nbf_i

          enddo
          n1 = nxtask(nproc,icount)
          end if
          
          n2 = n2 + 1

        enddo
      enddo

      n1 = nxtask(-nproc,icount)

      status = ma_pop_stack(lscr)
      status = ma_pop_stack(lMP)

      write(LuOut,*) 
      write(LuOut,*) ' SPM '
      write(LuOut,*) 
      if(lmax.gt.0) then
        write(LuOut,*) ' atom  charge        dipole'
        do iat=1,natoms
          
          write(LuOut,20) iat,charge(iat),(dipole(ish,iat),ish=1,3)
   20     format(i5,f12.6,2x,3f12.6)
        enddo
      else
        write(LuOut,*) ' atom  charge '
        do iat=1,natoms
          
          write(LuOut,10) iat,charge(iat)
   10     format(i5,f12.6)
        enddo
      endif

      return
      end
