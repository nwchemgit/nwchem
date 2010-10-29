C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_3j_prod (l1,l2,l3,m1,m2,result)
*                                                                      *
*   Calculate product of 3j symbols which appears in the angular       * 
*   integral over three spherical tensors. Product is                  *
*   (l1 l2 l3 / m1 m2 m3) (l1 l2 l3 / 0 0 0).                          *
*                                                                      *
*   l1,l2,l3 - l values                                                *
*   m1,m2 - m values; m3 = -(m1+m2)                                    *
*   result - value of product of 3j symbols                            *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer l1,l2,l3,m1,m2,m3
      integer l_min,l_mid,l_max,m_min,m_mid,m_max,
     &    ll,ll_min,ll_mid,ll_max,l_sum
      integer i,j,k,l,m,n,phase
      double precision result,wa,wb
*
      result = 0.0d00
      m3 = -(m1+m2)
      if (abs(m1) .gt. l1) return
      if (abs(m2) .gt. l2) return
      if (abs(m3) .gt. l3) return
*
*   Order angular momenta by size.
*
      if (l2 .gt. l1) then
        l_max = l2
        m_max = m2
        l_min = l1
        m_min = m1
      else
        l_max = l1
        m_max = m1
        l_min = l2
        m_min = m2
      end if
      if (l3 .gt. l_max) then
        l_mid = l_max
        m_mid = m_max
        l_max = l3
        m_max = m3
      else if (l3. lt. l_min) then
        l_mid = l_min
        m_mid = m_min
        l_min = l3
        m_min = m3
      else
        l_mid = l3
        m_mid = m3
      end if
      phase = (abs(m1)+abs(m2)+abs(m3))/2
      if (m_min .lt. 0) then
        m_min = -m_min
        m_mid = -m_mid
        m_max = -m_max
      end if
*
*   Special code for l_min = 0
*
      if (l_min .eq. 0) then
        wb = 1
        wa = 2*l_max+1
        result = wb/wa
*
*   Special code for l_min = 1
*
      else if (l_min .eq. 1) then
        wa = (2*l_max+1)*(2*l_max-1)
        if (m_min .eq. 0) then
          wb = (l_max+m_mid)*(l_max-m_mid)
        else
          wb = (l_max+m_mid)*(l_max+m_mid+1)/2
        end if
        result = sqrt(wb)/wa
*
*   General code for l_min > 1
*
      else
*
*     (a+b-c-1)!!(b+c-a-1)!!(c+a-b-1)!!/(a+b+c+1)!!
*
        l_sum = l1+l2+l3
        ll_max = 2*l_max
        ll_min = 2*l_min
        ll_mid = 2*l_mid
        m = 1
        do l = l_sum-ll_min+1,l_sum+1,2
          m = m*l
        end do
        wa = m
        m = 1
        do l = 3,l_sum-ll_mid-1,2
          m = m*l
        end do
        do l = 3,l_sum-ll_max-1,2
          m = m*l
        end do
        wb = m
        result = wb/wa
*
*     [(c+gamma)!(c-gamma)!/(a+alpha)!(a-alpha)!(b+beta)!(b-beta)!]^1/2
*
        m = 1
        n = 1
        ll_max = l_max+abs(m_max)
        ll_mid = l_mid+abs(m_mid)
        if (ll_max .ge. ll_mid) then
          do i = ll_mid+1,ll_max
            m = m*i
          end do
        else
          do i = ll_max+1,ll_mid
            n = n*i
          end do
        end if
        ll_max = l_max-abs(m_max)
        ll_mid = l_mid-abs(m_mid)
        if (ll_max .ge. ll_mid) then
          do i = ll_mid+1,ll_max
            m = m*i
          end do
        else
          do i = ll_max+1,ll_mid
            n = n*i
          end do
        end if
        do i = 1,l_min-m_min
          n = n*i
        end do
        do i = 1,l_min+m_min
          n = n*i
        end do
        wa = m
        wb = n
        result = result*sqrt(wa/wb)
*
*     Remaining factor
*
        ll = l_min+l_mid-l_max
        if (ll .gt. 0) then
          i = l_min+m_min
          j = l_min-m_min
          m = l_mid+m_mid
          n = l_mid-m_mid
          ll_min = max(0,ll-i,ll-n)
          ll_max = min(ll,j,m)
          k = 1
          do l = 1,ll_min
            k = -k*j*m/l
            j = j-1
            m = m-1
          end do
          do l = 1,ll-ll_min
            k = k*i*n/l
            i = i-1
            n = n-1
          end do
          l_sum = k
          do l = ll_min+1,ll_max
            i = i+1
            n = n+1
            k = -k*j*m*(ll-l+1)/(i*n*l)
            j = j-1
            m = m-1
            l_sum = l_sum+k
          end do
          wa = l_sum
          result = result*wa
        end if
	phase = phase+ll/2
      end if
*
      if (mod(abs(phase+m_max),2) .eq. 1) result = -result
*
      return
      end
