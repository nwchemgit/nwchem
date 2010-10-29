* $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t2_init3 (n,ldQ,m_min,m_max,j,h,tol,ai,bi,
     &    alpha,beta,gamma,prefactor,temp,ind,Q)
*                                                                      *
*   Calculate values of Q^3_{mm}, Q^4_{m+1m} and Q^4_{mm+1}            *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of functions to evaluate                          *
*   ldQ (inp) - leading dimension of array Q                           *
*   m_min (inp) - minimum value of m                                   *
*   m_max (inp) - maximum value of m                                   *
*   l_a (inp) - (maximum) angular momentum of functions on centre a    *
*   l_b (inp) - (maximum) angular momentum of functions on centre b    *
*   tol (inp) - maximum relative error in bessel functions             *
*   ai (inp) - 1/a                                                     *
*   bi (inp) - 1/b                                                     *
*   alpha (inp) -  a/2sqrt(c)                                          *
*   beta (inp) - b/2sqrt(c)                                            *
*   gamma (inp) - 1/sqrt(c)                                            *
*   prefactor (inp) - exponential prefactor (see calling routine)      *
*   temp - work array                                                  *
*   ind (scr) - index array for bessel driver                          *
*   Q (out) - uncontracted Q integrals                                 *
*                                                                      *
*   Written by W. A. de Jong                                           *
*                                                                      *
************************************************************************
      implicit none
      integer h,i,j,ja,jb,jl,js,k,ldQ,m,m_min,m_max,mm,mp,n,na,nb,nl,
     &    ns,nq,ind(n)
      double precision ai(n),bi(n),alpha(n),beta(n),gamma(n),
     &    prefactor(n),temp(n,18),Q(ldQ,m_min:m_max,2),
     &    tol,one,two,cutoff,big,afac,bfac
      parameter (one = 1.0d00, two = 2.0d00,
     &    cutoff = 24.5d00, big = 75.0d00)
*
*   Gather arguments for quadrature and single power series
*   Quadrature starts at alpha*beta=cutoff
*
      if (n .eq. 0) return
      nq = 0
      ns = 0
      do i = 1,n
        if (alpha(i)*beta(i) .ge. cutoff) then
          nq = nq+1
          ind(i) = nq
        else 
          ns = ns-1
          ind(i) = ns
        end if
      end do
*
*   Subdivide arguments for different ranges in the quadrature
*   and into alpha > beta and alpha < beta for the power series
*
      na = nq
      nb = n+1
      ns = 0
      nl = nq+1
      do i = 1,n
        if (ind(i) .gt. 0) then
          if (alpha(i)+beta(i) .le. big) then
            ns = ns+1
            ind(i) = ns
            temp(ns,1) = alpha(i)
            temp(ns,2) = beta(i)
            temp(ns,3) = gamma(i)
          else
            nl = nl-1
            ind(i) = nl
            temp(nl,1) = alpha(i)
            temp(nl,2) = beta(i)
            temp(nl,3) = gamma(i)
          end if
        else
          if (alpha(i) .le. beta(i)) then
            na = na+1
            ind(i) = na
            temp(na,1) = alpha(i)
            temp(na,2) = beta(i)
            temp(na,3) = gamma(i)
          else
            nb = nb-1
            ind(i) = nb
            temp(nb,2) = alpha(i)
            temp(nb,1) = beta(i)
            temp(nb,3) = gamma(i)
          end if
        end if
      end do
      js = 1
      jl = nl
      ja = nq+1
      jb = nb
      nl = nq-ns
      nb = n-na
      na = na-nq
*
*   Evaluate Q^3_{mm} functions. Note that the function for m_max is
*   not required for the recursion.
*
      do m = m_min,max(m_max-1,m_min)
        if (ns .gt. 0) call ecp_t2_ghq (3,m,m,ns,12,temp(js,1),
     &      temp(js,2),temp(js,3),temp(js,6),temp(js,4),temp(js,5),tol)
        if (nl .gt. 0) call ecp_t2_ghq (3,m,m,nl,6,temp(jl,1),
     &      temp(jl,2),temp(jl,3),temp(jl,6),temp(jl,4),temp(jl,5),tol)
        if (na .gt. 0) call ecp_t2_p3pow (na,m,0,temp(ja,1),temp(ja,2),
     &      temp(ja,3),temp(ja,4),temp(ja,6),temp(ja,7),temp(ja,5),tol)
        if (nb .gt. 0) call ecp_t2_p3pow (nb,m,0,temp(jb,1),temp(jb,2),
     &      temp(jb,3),temp(jb,4),temp(jb,6),temp(jb,7),temp(jb,5),tol)
        do i = 1,n
          Q(i,m,1) = prefactor(i)*temp(ind(i),5)
        end do
      end do
*
      if (m_min .eq. m_max) return
*
*   Evaluate Q^3_{mm-1} or Q^3_{m-1m} functions
*
      k = j-h
      do m = m_max-1,max(m_max-2,m_min),-1
        if (ns .gt. 0) call ecp_t2_ghq (4,m+j,m+h,ns,12,temp(js,1),
     &      temp(js,2),temp(js,3),temp(js,6),temp(js,4),temp(js,5),tol)
        if (nl .gt. 0) call ecp_t2_ghq (4,m+j,m+h,nl,6,temp(jl,1),
     &      temp(jl,2),temp(jl,3),temp(jl,6),temp(jl,4),temp(jl,5),tol)
        if (na .gt. 0) call ecp_t2_p3pow (na,m+h,k,temp(ja,1),
     &      temp(ja,2),temp(ja,3),temp(ja,4),temp(ja,6),temp(ja,7),
     &      temp(ja,5),tol)
        if (nb .gt. 0) call ecp_t2_p3pow (nb,m+j,-k,temp(jb,1),
     &      temp(jb,2),temp(jb,3),temp(jb,4),temp(jb,6),temp(jb,7),
     &      temp(jb,5),tol)
        do i = 1,n
          Q(i,m,2) = prefactor(i)*temp(ind(i),5)
        end do
      end do
*
      do m = m_max-2,m_min+1,-1
        mm = m-1
        mp = m+1
        afac = 2*mp+k
        bfac = 2*mp-k
        do i = 1,n
          Q(i,mm,2) = Q(i,mp,2) 
     &        + afac*Q(i,m+j,1)*ai(i)
     &        + bfac*Q(i,m+h,1)*bi(i)
        end do
      end do
*
      return
      end
