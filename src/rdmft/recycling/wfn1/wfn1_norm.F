C> \ingroup wfn1
C> @{
C>
C> \brief Compute the norm of the wfn1 gradient
C>
C> The gradient of the energy wrt the wavefunction parameters is 
C> defined by a variety of coefficients. This routine computes a 
C> 2-norm of this gradient by summing the square of all numbers 
C> and finally taken the square root of that sum.
C>
      double precision function wfn1_norm(nbf,r_dnoa,r_dnob,
     +       r_dofa,r_dofb,r_dl1a,r_dl1b,r_dl2a,r_dl2b,r_dl3a,r_dl3b)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      double precision r_dnoa(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< alpha natural orbitals
      double precision r_dnob(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< beta natural orbitals
      double precision r_dofa(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< alpha occupation functions
      double precision r_dofb(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< beta occupation functions
      double precision r_dl1a(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< alpha first lagrangian
      double precision r_dl1b(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< beta first lagrangian
      double precision r_dl2a(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< alpha second lagrangian
      double precision r_dl2b(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< beta second lagrangian
      double precision r_dl3a(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< alpha third lagrangian
      double precision r_dl3b(nbf,nbf) !< [Input] Derivatives wrt the
                                       !< beta third lagrangian
c
c     Local variables
c
      integer ix, iy !< Counters
c
      double precision sum
      double precision sum1
c
c     Code
c
      sum = 0.0d0
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dnoa(ix,iy)*r_dnoa(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dnoa = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dnob(ix,iy)*r_dnob(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dnob = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dofa(ix,iy)*r_dofa(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dofa = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dofb(ix,iy)*r_dofb(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dofb = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl1a(ix,iy)*r_dl1a(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl1a = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl1b(ix,iy)*r_dl1b(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl1b = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl2a(ix,iy)*r_dl2a(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl2a = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl2b(ix,iy)*r_dl2b(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl2b = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl3a(ix,iy)*r_dl3a(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl3a = ',sum1
cDEBUG
      sum1 = 0.0d0
      do iy = 1, nbf
        do ix = 1, nbf
          sum1 = sum1 + r_dl3b(ix,iy)*r_dl3b(ix,iy)
        enddo
      enddo
      sum = sum + sum1
cDEBUG
      write(*,*)'r_dl3b = ',sum1
cDEBUG
      wfn1_norm = sqrt(sum)
c
      return
      end
C> @}

