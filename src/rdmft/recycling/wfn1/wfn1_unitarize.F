C> \ingroup wfn1
C> @{
C>
C> \brief Make a nearly unitary matrix properly unitary
C>
C> Use an iterative algorithm to make a nearly unitary matrix truly
C> unitary. The equation
C> \f{eqnarray*}{
C>   C' &=& C*\left(\frac{3}{2}I-\frac{1}{2}C^T*S*C\right)
C> \f}
C> produces a matrix \f$C'\f$ that is closer to unitarity than \f$C\f$.
C> This equation is iteratively applied until \f$C\f$ is unitary.
C>
C> In practice I have found that occasionally the above approach does
C> not seem to reach machine precision. Worse, when it cannot reach
C> machine precision the method starts to diverge and eventually
C> the calculation explodes with a floating point exception. As this
C> code is a proof of principle implementation it is more important
C> to get robust results than performance. So it is desirable to use
C> something else. 
C>
C> The simplest alternative is to use a Modified Gramm-Schmidt 
C> orthogonalization (repeated should not be needed as the matrix
C> should be nearly unitary, nevertheless we will check and bomb if
C> suspicious results are found).
C>
      subroutine wfn1_unitarize(nbf,c,ov,s,cp,ouse_ov)
      implicit none
      integer nbf !< [Input] The dimension of the matrices
      double precision c(nbf,nbf)  !< [In/Output] The matrix being 
                                   !< unitarized
      double precision ov(nbf,nbf) !< [Input] The overlap matrix
      double precision s(nbf,nbf)  !< [Scratch] mainly \f$C*C^T\f$
      double precision cp(nbf,nbf) !< [scratch]
      logical ouse_ov              !< [Input] If `.TRUE.` use the 
                                   !< overlap matrix as metric, 
                                   !< otherwise use the unit matrix
                                   !< as metric.
c
      double precision dnrm2 !< Another norm
      double precision dnrm0 !< Old norm of the error
      double precision dnrm !< Norm of the error
      integer iter          !< Iteration counter
      integer i,j           !< counters
c
      double precision ddot
      external         ddot
c
cDEBUG
c     write(50,*)'wfn1_unitarize: ',nbf
c     do i=1,nbf
c       write(50,*)i,(c(i,j),j=1,nbf)
c     enddo
cDEBUG
c     iter = 0
c     dnrm0 = 10.0d0
c     dnrm = 1.0d0
c     do while (dnrm.lt.dnrm0.and.dnrm.gt.1.0d-15)
c        iter = iter + 1
c        if (ouse_ov) then
c          call dgemm('n','n',nbf,nbf,nbf,1.0d0,ov,nbf,c,nbf,0.0d0,
c    &                cp,nbf)
c          call dgemm('t','n',nbf,nbf,nbf,1.0d0,c,nbf,cp,nbf,0.0d0,
c    &                s,nbf)
c        else
c          call dgemm('t','n',nbf,nbf,nbf,1.0d0,c,nbf,c,nbf,0.0d0,s,nbf)
c        endif
c        do i = 1, nbf
c          s(i,i) = s(i,i) - 1.0d0
c        enddo
c        dnrm0 = dnrm
c        dnrm = 0.0d0
c        do j = 1, nbf
c          do i = 1, nbf
c            dnrm = max(dnrm,abs(s(i,j)))
c          enddo
c        enddo
cDEBUG
c        write(50,*)'*** wfn1_unitarize: it, norm = ',iter,dnrm
cDEBUG
c        if (.not.(dnrm.lt.dnrm0.and.dnrm.gt.1.0d-15)) then
c
c          My result is not going to improve so get out
c
c          goto 10
c        endif
c        call dscal(nbf*nbf,-0.5d0,s,1)
c        do i = 1, nbf
c          s(i,i) = s(i,i) + 1.0d0
c        enddo
c        call dgemm('n','n',nbf,nbf,nbf,1.0d0,c,nbf,s,nbf,0.0d0,cp,nbf)
c        call dcopy(nbf*nbf,cp,1,c,1)
c     enddo
c10   continue
c
cDEBUG
c     close(50)
c     write(*,*)'*** wfn1_unitarize: input'
c     call hess_hssout(c,nbf,nbf,nbf)
cDEBUG
      dnrm0 = 1.0d0
      if (dnrm0.gt.2.0d-15) then
c
c       Did not achieve machine precision, do modified Gramm-Schmidt
c
        dnrm0 = 0.0d0
        do i = 1, nbf
          if (ouse_ov) then
            call dgemv('n',nbf,nbf,1.0d0,ov,nbf,c(1,i),1,0.0d0,cp,1)
            dnrm2 = ddot(nbf,c(1,i),1,cp,1)
          else
            dnrm2 = ddot(nbf,c(1,i),1,c(1,i),1)
          endif
          dnrm2 = 1.0d0/dsqrt(dnrm2)
          call dscal(nbf,dnrm2,c(1,i),1)
          do j = 1, i-1
            if (ouse_ov) then
              call dgemv('n',nbf,nbf,1.0d0,ov,nbf,c(1,i),1,0.0d0,cp,1)
              dnrm2 = ddot(nbf,c(1,j),1,cp,1)
            else
              dnrm2 = ddot(nbf,c(1,j),1,c(1,i),1)
            endif
            call daxpy(nbf,-dnrm2,c(1,j),1,c(1,i),1)
          enddo
          if (ouse_ov) then
            call dgemv('n',nbf,nbf,1.0d0,ov,nbf,c(1,i),1,0.0d0,cp,1)
            dnrm2 = ddot(nbf,c(1,i),1,cp,1)
          else
            dnrm2 = ddot(nbf,c(1,i),1,c(1,i),1)
          endif
          dnrm0 = max(dnrm0,abs(dnrm2-1.0d0))
          dnrm2 = 1.0d0/dsqrt(dnrm2)
          call dscal(nbf,dnrm2,c(1,i),1)
        enddo
      endif
      if (dnrm0.gt.1.0d-13) then
        write(0,*)'WARNING: wfn1_unitarize: Gramm-Schmidt inaccurate'//
     &            ': ',dnrm0
      endif
c
c     Check unitarity
c
      if (ouse_ov) then
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,ov,nbf,c,nbf,0.0d0,
     &             cp,nbf)
        call dgemm('t','n',nbf,nbf,nbf,1.0d0,c,nbf,cp,nbf,0.0d0,
     &             s,nbf)
      else
        call dgemm('t','n',nbf,nbf,nbf,1.0d0,c,nbf,c,nbf,0.0d0,s,nbf)
      endif
      do i = 1, nbf
        s(i,i) = s(i,i) - 1.0d0
      enddo
      dnrm = 0.0d0
      do j = 1, nbf
        do i = 1, nbf
          dnrm = max(dnrm,abs(s(i,j)))
        enddo
      enddo
      if (dnrm.gt.1.0d-14) then
        write(0,*)'WARNING: wfn1_unitarize: unitarity test failed: ',
     &            dnrm
cDEBUG
c       call hess_hssout(s,nbf,nbf,nbf)
c       call hess_hssout(c,nbf,nbf,nbf)
cDEBUG
      endif
      return
      end
C> @}
