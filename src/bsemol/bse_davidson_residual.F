      subroutine bse_davidson_residual(nroots,ntrials,omega,x,ym,r,dia,
     $                                 w,npoles,ndim,norm)
      implicit none

      logical main
      integer nroots,npoles,ntrials,ndim
      double precision omega(nroots)
      double precision ym(npoles,nroots),r(npoles)
      double precision x(npoles,nroots),dia(npoles),w(npoles,nroots)

      integer iroot
      double precision norm(nroots)

      ndim = 0
      do iroot=1,nroots
        r(:) = ym(:,iroot) - omega(iroot)*x(:,iroot)
        norm(iroot) = maxval(abs(r))
        call ga_dgop((/1/),norm(iroot),1,"maxabs")
        if (norm(iroot).gt.1d-5) then
          ndim = ndim + 1
          w(:,ndim) = r(:)/(dia(:) - omega(iroot))
        endif
      enddo

      end subroutine




