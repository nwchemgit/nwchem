c
      subroutine print_dens_bit(d_ij, ndens, ibflo, ibfhi,
     $                    jbflo, jbfhi, blen)
      implicit none
c
      integer ndens, ibflo, ibfhi, jbflo, jbfhi, blen
      double precision d_ij(ndens,blen,blen)
      integer ndns, id, jd
c
      write(6,*) 'density bit'
      do ndns = 1, ndens
        write(6,*) 'ndns is ',ndns
        do id = 1, 1+ ibfhi -ibflo
          do jd = 1, 1+ jbfhi -jbflo
            write(6,*) id,jd,d_ij(ndns,id,jd)
          enddo
        enddo
      enddo
c
      return
      end
c $Id$
