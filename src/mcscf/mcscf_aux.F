c
c  Convert occupation vector to 
c  number of closed and active orbitals
c
      subroutine mcscf_occ2int( nbf, occ, nclosed, nact )
      implicit none
      integer nbf
      double precision occ(nbf)
      integer nclosed
      integer nact
c
      integer i
c
      nclosed = 0
      nact = 0
      do i=1,nbf
        if (occ(i).eq.2.d0) then
          nclosed = nclosed + 1
        else if (occ(i).ne.0.d0) then
          nact = nact + 1
        endif
      enddo
      return
      end


      









c
c  Convert occupation vector
c  to number of electrons (trace 1pdm)
c
      subroutine mcscf_occ2elec( nbf, occ, nelec )
      implicit none
      integer nbf, nelec
      double precision occ(nbf)
      integer i
      double precision xelec

      xelec = 0.d0
      do i=i,nbf
        xelec = xelec + occ(i)
      enddo
      nelec = xelec
      return
      end







c
c Convert the vector representation of the MCSCF
c parameters to a (nbf x nbf) matrix representation
c Accumulates into the matrix (does not copy!)
c
c
c

       subroutine mcscf_vec2mat( nbf, nclosed, nact, a, g_v, b, g_m )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision b,a                         ! [input] Additive factors
       integer g_v                                  ! [input] Vector representation
       integer g_m                                  ! [output] Matrix representation
c
c
       integer nvir, vlen
       integer voff, aoff, aend, xoff, xend
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
c
       call ga_dadd_patch( a, g_v, 1, (nvir*nclosed), 1, 1,
     $                     b, g_m, voff, nbf, 1, nclosed,
     $                        g_m, voff, nbf, 1, nclosed )
       xoff = nclosed*nvir + 1
       xend = nclosed*nvir + nact*nvir
       call ga_dadd_patch( a, g_v, xoff, xend, 1, 1,
     $                     b, g_m, voff, nbf, aoff, aend,
     $                        g_m, voff, nbf, aoff, aend )
       xoff = xend + 1
       xend = xend + nclosed*nact
       call ga_dadd_patch( a, g_v, xoff, xend, 1, 1,
     $                     b, g_m, aoff, aend, 1, nclosed,
     $                        g_m, aoff, aend, 1, nclosed )
c
c
c
       return
       end

       


c
c Convert the matrix representation of the MCSCF
c parameters to vector representation (complement of mcscf_vec2mat)
c Accumulates into the vector (does not copy!)
c
c
c

       subroutine mcscf_mat2vec( nbf, nclosed, nact, a, g_m, b, g_v )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision b,a                         ! [input] Additive factors
       integer g_m                                  ! [input] Matrix representation
       integer g_v                                  ! [output] Vector representation
c
c
       integer nvir, vlen
       integer voff, aoff, aend, xoff, xend
c
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
c
       xend = nvir*nclosed
       call ga_dadd_patch( a, g_m, voff, nbf, 1, nclosed,
     $                     b, g_v, 1, xend, 1, 1,
     $                        g_v, 1, xend, 1, 1)
       xoff = xend + 1
       xend = xend + nact*nvir
       call ga_dadd_patch( a, g_m, voff, nbf, aoff, aend,
     $                     b, g_v, xoff, xend, 1, 1,
     $                        g_v, xoff, xend, 1, 1)
       xoff = xend + 1
       xend = vlen
       call ga_dadd_patch( a, g_m, aoff, aend, 1, nclosed,
     $                     b, g_v, xoff, xend, 1, 1,
     $                        g_v, xoff, xend, 1, 1)
       return
       end









c
c  Create closed and active AO densities from
c  occupation numbers and 1-pdm in active space
c

       subroutine mcscf_occ2dens( nbf, nclosed, nact, dm1,
     $                            g_movecs, g_cdens, g_adens )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
c
       integer nbf
       integer nclosed
       integer nact
       double precision dm1(nact,nact)
       integer g_movecs
       integer g_cdens
       integer g_adens
c
       integer aoff, aend
       integer g_dm1
c
c
c
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ga_create(MT_DBL,nact,nact,'1PDM',nact,1,g_dm1))
     $   call errquit('mcscf_occ2dens: cannot create 1PDM',0)
       if (ga_nodeid().eq.0)
     $   call ga_put(g_dm1,1,nact,1,nact,dm1,nact)
       call ga_sync()
c
c Active AO density (use closed shell dens for temp)
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_dm1, 1, nact, 1, nact,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_cdens, 1, nact, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_cdens, 1, nact, 1, nbf,
     $                       g_adens, 1, nbf, 1, nbf )
c
c
       if (.not.ga_destroy(g_dm1))
     $   call errquit('mcscf_occ2dens: cannot destroy',0)
c
c Closed AO density
c
       call ga_dgemm( 'n', 't', nbf, nbf, nclosed, 2.d0,
     $                g_movecs, g_movecs, 0.d0, g_cdens )
       return
       end










      subroutine mcscf_symmetrize_2pdm( nact, dm2, dmsy2 )
      implicit none
      integer nact
      double precision dm2(nact,nact,nact,nact)
      double precision dmsy2(nact,nact,nact,nact)
      integer i,j,k,l
      
      do i=1,nact
        do j=1,nact
          do k=1,nact
            do l=1,nact
              dmsy2(i,j,k,l) = 0.5d0*(dm2(i,j,k,l) + dm2(i,j,l,k))
            enddo
          enddo
        enddo
      enddo
      return
      end

      





















