C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a matrix that can be used to enumerate vectors
C>
C> In order to create a delta-SCF approach to RDMFT a way is needed to
C> order a set of orbitals. The ordering should be such that two sets
C> of vectors will be ordered according to similarity. 
C>
C> Consider a diagonal matrix \f$ \diag(A_{11}, ..., A_{NN}) \f$ where
C> \f$ A_{ii} = i N^3 \f$. Also there is a set of vectors \f$ v_i \f$
C> such that \f$ (v_i|v_j) = \delta_{ij} \f$. Then the expectation 
C> value is \f$ R = \sum_k A_{kk}v_{ki}v_{ki} = \sum_k A_{kk}w_k \f$
C> where \f$w_k = v_{ki}v_{ki} \f$. Obviously all numbers \f$w_k\f$ are
C> non-negative and they sum to \f$1\f$. If we create \f$v_2\f$ from
C> \f$v_1\f$ by changing coefficient \f$v_{11}\f$ then some other
C> coefficients must change as well to maintain the normalization. 
C> As the change of \f$\partial R/\partial v_{ki}\f$ is different for
C> every \f$k\f$ we have that two vectors with different values of the
C> coefficients generate different values for \f$R\f$.
C>
C> The diagonal matrix discussed above does not solve the problem when
C> two vectors only differ in the signs of some of the coefficients.
C> However the expectation values of two vectors that differ only in
C> the signs of some coefficients will differ in the sign of the cross
C> terms. Hence if we give every off-diagonal element a unique value
C> then a change of sign in the cross terms must lead to a change in
C> \f$ R \f$. 
C>
C> Hence we construct the matrix
C> \f{eqnarray*}{
C>   A_{ii} &=& i N^3 \\\\
C>   A_{ij} &=& A_{ji} \\\\
C>   A_{ij} &=& (i-2)(i-1)/2+j \; i \gt j
C> \f}
C> that we can use to compute the expectation values that we may use
C> to order the vectors. We also need to transform the matrix so it
C> can be used on the natural orbitals.
C> 
      subroutine wfn1_nxt2_enumerate(nbf,a_no,a_ao,s12,w)
      implicit none
      integer nbf !< [Input] the number of basis functions
      double precision a_no(nbf,nbf) !< [Output] matrix A in an 
                                     !< orthogonal basis
      double precision a_ao(nbf,nbf) !< [Output] matrix A in the
                                     !< AO basis
      double precision s12(nbf,nbf)  !< [Input] Sqrt of the overlap 
                                     !< matrxi
      double precision w(nbf,nbf)    !< [Scratch]
c
      integer ii !< Counter
      integer jj !< Counter
c
      double precision one  !< the number 1.0
      double precision zero !< the number 0.0
      parameter(zero = 0.0d0)
      parameter(one  = 1.0d0)
c
      do ii = 1, nbf
c       a_no(ii,ii) = ii*(dble(nbf)**3)
        a_no(ii,ii) = ii
      enddo
      do jj = 1, nbf
        do ii = jj+1, nbf
c         a_no(ii,jj) = (ii-2)*(ii-1)/2+jj
          a_no(ii,jj) = nbf*(dble(ii-2)*dble(ii-1)/2+jj+1)
          a_no(jj,ii) = a_no(ii,jj)
        enddo
      enddo
      call dgemm('n','n',nbf,nbf,nbf,one,s12,nbf,a_no,nbf,zero,w,nbf)
      call dgemm('n','n',nbf,nbf,nbf,one,w,nbf,s12,nbf,zero,a_ao,nbf)
c
      end
C>
C> @}
