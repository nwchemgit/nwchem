      subroutine dft_genspm(lmax,iga_dens,basis,
     &     natoms,centerl,Pmat,max_at_bf,
     &     charge,dipole,coord)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "stdio.fh"
      integer int_nint,nxtask
      external int_nint,nxtask
      double precision ddot
      external ddot

      integer natoms
      integer iga_dens,basis,lmax
      integer max_at_bf
      double precision Pmat(max_at_bf*max_at_bf)

      double precision charge(natoms)
      double precision dipole(3,natoms)
      double precision coord(3,natoms)
      integer lscr,iscr,iMP,lMP
      integer lmpint              

      integer ifirst,ilast,iat,idim,ish,ifirst_sh,ilast_sh,nbf_i,jfunc
      integer jfirst,jlast,jat,jdim,jsh,jfirst_sh,jlast_sh,nbf_j
      integer max_scratch
*rak:      integer num_int,lval,ij_int,maxg,nshells
      integer istart,Pmat_point,Pmat_point0,istartd
      integer n1,n2,icount,nproc
      double precision centerl(3)

*rak:c     
*rak:c     evaluation of # of multipole ints
*rak:c     
*rak:      if( .not. bas_numcont(basis,nshells) )
*rak:     &     call errquit('Exiting in dft_gensmp.',3)
*rak:
*rak:      lmpint=0
*rak:      do ish=1,nshells
*rak:        do jsh=1,ish
*rak:          ij_int = int_nint(basis, ish, basis, jsh, 0,0, 0,0)
*rak:          num_int=0
*rak:          do lval = 0,lmax
*rak:            num_int = num_int+ 1 + ij_int*(lval+1)*(lval+2)/2
*rak:          enddo
*rak:          lmpint=max(lmpint,num_int)
*rak:        enddo
*rak:      enddo

      nproc  = ga_nnodes()

      call int_mem_dipole(lmpint, max_scratch, basis, basis, lmax)

      if(.not.MA_Push_Get(mt_dbl,lmpint,' multipole ints',lMP,iMP))
     &  call errquit('dft_genspm: cannot allocate multipole ints',0,
     &       MA_ERR)

      if(.not.MA_Push_Get(mt_dbl,max_scratch,' scratch',lscr,iscr))
     &  call errquit('dft_genspm: cannot allocate scratch',0,
     &       MA_ERR)


      call dfill(natoms,0.d0,charge,1)
      if(lmax.gt.0) then
      call dfill(natoms*3,0.d0,dipole,1)
      endif

      icount=1
      n1 = nxtask(nproc,icount)
      n2 = 0

      do iat=1,natoms
        if (.not.bas_ce2cnr(basis,iat,ifirst,ilast))
     &       call errquit('Exiting in dft_gensmp',1, BASIS_ERR)

        do jat=1,natoms
          if(n1.eq.n2) then
          if (.not.bas_ce2cnr(basis,jat,jfirst,jlast))
     &         call errquit('Exiting in dft_gensmp',2, BASIS_ERR)

          call get_atom_block(iga_dens, basis,
     $         iat, jat, Pmat, idim, jdim)
          Pmat_point0=1
          do ish=ifirst,ilast
            Pmat_point=Pmat_point0
            
            if (.not. bas_cn2bfr( basis,ish,
     &           ifirst_sh,ilast_sh))
     &           call errquit('Exiting in dft_gensmp.',5, BASIS_ERR)
            nbf_i=ilast_sh-ifirst_sh+1

            do jsh=jfirst,jlast
              
              if (.not. bas_cn2bfr( basis,jsh,
     &             jfirst_sh,jlast_sh))
     &             call errquit('Exiting in dft_gensmp.',5, BASIS_ERR)
              nbf_j=jlast_sh-jfirst_sh+1

              call int_mpole(basis, jsh, basis, ish, lmax, centerl,
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

      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_genspm: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lMP))
     &  call errquit('dft_genspm: cannot pop stack',0, MA_ERR)

      write(LuOut,*) 
      write(LuOut,*) ' SPM '
      write(LuOut,*) 
      if(lmax.gt.0) then
        write(LuOut,*) ' atom    x         y          z',
     * '         charge    dipole'
        do iat=1,natoms
          
          write(LuOut,20) iat,(coord(ish,iat),ish=1,3),
     *         charge(iat),(dipole(ish,iat),ish=1,3)
   20     format(i5,3f10.4,1x,f10.4,1x,3f10.4)
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
