c***********************************************************************
c
c-->  Construct the density matrix.
c
c***********************************************************************
c
      subroutine dft_densm(idens,ivecs,nbf,ipol,noc)
C$Id: dft_densm.F,v 1.1 1995-04-27 03:16:18 og845 Exp $
      implicit real*8 (a-h,o-z)

c
      fact=2.d0/ipol
      zero=0.d0
         call ga_dgemm('n', 't', nbf, nbf, noc, fact, ivecs,
     $        ivecs, zero, idens)
      return
      end
