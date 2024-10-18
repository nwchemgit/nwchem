c
C> \brief Calculate the screened COSMO charges
C>
C> ## Introduction ##
C>
C> In the COSMO model [1] the charges to represent the solvation effects
C> are obtained from solving a linear system of equations \f$ Ax=b \f$.
C> In this system \f$ b \f$ is the electrostatic potential at the 
C> point charge positions. The matrix \f$ A \f$ consists of the 
C> electrostatic interaction between two unit charges at the point
C> charge position, i.e.
C> \f{eqnarray*}{
C>   A_{\mu\nu} &=& ||t_\mu - t_\nu||^{-1} \\\\
C>   A_{\mu\mu} &=& 3.8 |S_\mu|^{-1/2} \\\\
C>              &=& a_{\mathrm{diag}} |S_\mu|^{-1/2}
C> \f}
C> where \f$ S_\mu \f$ is the surface area associated with the COSMO
C> charge (see [1] Eqs. 7a and 7b). The quantity \f$a_{\mathrm{diag}}\f$
C> is evaluated using Eq.(B1). This is done in two steps in that the
C> this quantity is evaluated for the unit sphere in `hnd_cossph`, any
C> remaining scale factors are applied in the evaluation of \f$A\f$.
C>
C> The original COSMO approach has problems when charges approach each
C> other and the interaction becomes singular. Therefore we have 
C> added a model where the interactions are smooth [4]
C> which is outlined below.
C>
C> ## Surface charges satisfying Gauss's theorem ##
C>
C> Solving these equations gives the "raw" COSMO charges \f$ x \f$.
C> These charges should sum up to the total charge contained within 
C> the Solvent Accessible Surface. For 2 reasons that will be in
C> practice not be exactly true:
C>
C> - The discretization of the SAS is not perfect
C>
C> - The electron distribution being represented with Gaussian functions
C>   extends beyond the SAS
C>
C> Therefore the raw COSMO charges are corrected by adding corrections
C> based on a Lagrange multiplier technique [4]. The corresponding
C> equations can be derived by starting from an energy expression
C> in terms of the solute charge distribution and the surface charges.
C> This energy expression including the Lagrange term is
C> \f{eqnarray*}{
C>   E(Q,q,\lambda) &=& \frac{1}{2}Q^TCQ + Q^TBq
C>                   +  \frac{1}{2f(\epsilon)}q^TAq
C>                   +  \lambda\left(f(\epsilon)Q_{in}+\sum_i q_i\right)
C> \f}
C> where \f$Q\f$ is the charge distribution of the solute (including
C> both nucleii and electrons), \f$q\f$ is the surface charge
C> distribution, \f$f(\epsilon)\f$ is the dielectric screening constant
C> as discussed below, \f$Q_{in}\f$ is the charge within the cavity
C> or equivalently the solute charge, and finally \f$A\f$, \f$B\f$,
C> and \f$C\f$ are Coulomb matrices.
C>
C> From this equation the surface charges can be derived by minimizing
C> \f$E\f$ wrt \f$q\f$ and \f$\lambda\f$. This yields
C> \f{eqnarray*}{
C>    \frac{\partial E}{\partial q} 
C>       &=& BQ + \frac{1}{f(\epsilon)}Aq + \Lambda = 0\\\\
C>    \frac{\partial E}{\partial \lambda}
C>       &=& f(\epsilon)Q_{in} + \sum_i q_i = 0
C> \f}
C> where \f$\Lambda\f$ is a vector of which each element is
C> \f$\lambda\f$, i.e. \f$\forall_i, \Lambda_i = \lambda\f$. 
C> Next we get
C> \f{eqnarray*}{
C>   q &=& -f(\epsilon)A^{-1}\left(BQ+\Lambda\right) \\\\
C>   \lambda &=& \frac{Q_{in}-\sum_i\left[A^{-1}BQ\right]_i}{
C>                     \sum_{ij}A^{-1}_{ij}}
C> \f}
C> Because \f$E\f$ is variationally optimized wrt \f$Q\f$, \f$q\f$,
C> and \f$\lambda\f$ the gradient expression only involves derivatives
C> of \f$A\f$, \f$B\f$ and \f$C\f$ just like the original COSMO
C> gradients [1].
C>
C> Previously the surface charge correction was implemented by scaling
C> the raw COSMO charges.
C> However, this led to complications with neutral molecules
C> where the correct integrated surface charge is 0. Hence the 
C> correction factor would be 0 as well, eliminating the COSMO charges
C> and hence all solvation effects. This problem had been patched by
C> calculating the COSMO charges for the nucleii and the electrons 
C> separately. This led to a cumbersome and expensive algorithm. So
C> the use of a Lagrange constraint is a solution that is
C> preferred over scaling the charges. The difference between the two
C> approaches should be small provided the discretization is fine
C> enough.
C>
C> In the COSMO model [1] it is realized that dielectric screening
C> scales as
C> \f{eqnarray*}{
C>   f(\epsilon) &= \frac{\epsilon-1}{\epsilon+a}, & 0\le a \le 2
C> \f}
C> Klamt and Sch&uuml;&uuml;rmann suggested to use \f$ a = 1/2 \f$, 
C> essentially based on an argument that the true value should not 
C> exceed 1 (see appendix A). Stevanovich and Truong [3] realized that
C> the screened charges should be such that the Gauss theorem holds,
C> which requires \f$ a = 0 \f$ (see Eq. (5)). Based on this physical
C> motivation the latter value is used by default.
C>
C> The linear system of equations to be solved may be tackled in 2
C> different ways. For small systems a direct solver is appriopriate,
C> whereas for large systems an iterative solver is used. These solvers
C> have different implications. If \f$ N \f$ is the number of COSMO
C> charges then
C>
C> - the direct solver uses \f$ O(N^2) \f$ memory and \f$ O(N^3) \f$
C>   instructions
C>
C> - the iterative solver uses \f$ O(N) \f$ memory and \f$ k*O(N^2) \f$
C>   instructions (\f$ k \f$ is the number of iterations to convergence)
C>
C> The costs are based on the assumptions that for the direct solver
C> the matrix is stored explicitly and the inverse is not stored but
C> recalculated every time. For the iterative solver only a 
C> matrix-vector multiplication is implemented that regenerates the
C> matrix elements every time, and the number of iterations is roughly
C> independent from the values of the matrix. 
C>
C> Based on these assumptions the iterative solver is optimal for large
C> systems both with respect to memory requirements as well as compute
C> requirements. The case for the iterative solver can be improved 
C> further by parallelizing the matrix-vector multiplication which
C> reduces the compute cost per processor to \f$ k*O(N^2)/N_{proc} \f$.
C> Furthermore if we start the iterative solver in each SCF cycle with
C> the solution from the previous iteration rather than \f$ x=0 \f$ then
C> \f$ k \f$ may be reduced as well. In practice \f$ k \f$ does not
C> depend strongly on the initial value of \f$ x \f$, reductions by at
C> most a factor 2 are seen when the SCF is nearly converged.
C>
C> ## Matching terms up ##
C>
C> When the COSMO charges have been calculated as described above on 
C> would like to evaluate the \f$E(Q,q,\lambda)\f$ expression to
C> obtain the energy of a solvated system. To do this successfully
C> we need to account for every contribution and see where that is
C> calculated. The \f$E(Q,q,\lambda)\f$ lists a number of terms:
C>
C> - The \f$\frac{1}{2}Q^TCQ\f$ term includes all terms that involve
C>   the quantum system in vacuum. In particular \f$Q\f$ combines the
C>   nuclear charges as well as the electron density. 
C>
C> - The \f$Q^TBq\f$ term includes all interactions between the charges
C>   of the QM system and the COSMO point charges.
C>
C> - The \f$\frac{1}{2f(\epsilon)}q^TAq\f$ term includes the self
C>   interaction between the COSMO charges.
C>
C> In reality the QM code must optimize the electronic structure under
C> the influence of the COSMO charges. So while the total energy of
C> the solvated system is calculated as 
C> \f{eqnarray*}{
C>    E_{\mathrm{tot}} &=&
C>    E_1 + E_2 + E_{\mathrm{nuc}} + E_{\mathrm{COSMO}}
C> \f}
C> these terms involve non-trivial combinations of terms in 
C> \f$E(Q,q,\lambda)\f$. The contributions included are
C>
C> - The term \f$E_1\f$ includes the electron-nuclear attraction 
C>   as well as the electron-COSMO charge interaction.
C>
C> - The term \f$E_2\f$ is just the 2-electron term.
C>
C> - The term \f$E_{\mathrm{nuc}}\f$ just contains the nuclear-nuclear
C>   repulsion.
C>
C> - The \f$E_{\mathrm{COSMO}}\f$ is a catch all entity that makes
C>   \f$E_{\mathrm{tot}}\f$ match \f$E(Q,q,\lambda)\f$.
C>
C> In particular we need to establish what \f$E_{\mathrm{COSMO}}\f$
C> actually is. Because the various terms partition in different ways
C> we need to separate \f$Q\f$ into the charge of the nucleii and the
C> electrons, i.e. \f$Q = Q_n + Q_e\f$. Within this context the terms 
C> breakdown in the following way:
C>
C> - \f$\frac{1}{2}Q^TCQ=\frac{1}{2}(\sum_i\frac{p_i^2}{m_e}+Q_n^TCQ_n
C>   +2Q_n^TCQ_e+Q_e^TCQ_e)\f$
C>
C> - \f$Q^TBq=Q_n^TBq+Q_e^TBq\f$
C>
C> - \f$\frac{1}{2f(\epsilon)}q^TAq\f$
C>
C> In the code we have
C>
C> - \f$E_1=\frac{1}{2}\sum_i\frac{p_i^2}{m_e}+Q_e^TCQ_n+Q_e^TBq\f$
C>
C> - \f$E_2=\frac{1}{2}Q_e^TCQ_e\f$
C>
C> - \f$E_{\mathrm{nuc}}=\frac{1}{2}Q_n^TCQ_n\f$
C>
C> Equating \f$E(Q,q,\lambda)\f$ to \f$E_{\mathrm{tot}}\f$ we find
C> that
C>
C> - \f$E_{\mathrm{COSMO}}=Q_n^TBq+\frac{1}{2f(\epsilon)}q^TAq\f$
C>
C> Essentially this expression is evaluated at "COSMO contribution
C> Alternative 1". Alternatively one could start from Ref.[1] Eq.(11)
C> where \f$E(Q,q)=\frac{1}{2}Q^T(C-BA^{-1}B)Q=\frac{1}{2}(Q^TCQ+Q^TBq)\f$.
C> Equating \f$E(Q,q)\f$ to \f$E_{\mathrm{tot}}\f$ from which we find
C> that
C>
C> - \f$E_{\mathrm{COSMO}}=\frac{1}{2}(Q_n^TBq-Q_e^TBq)\f$
C>
C> This expression is evaluated at "COSMO contribution Alternative 2".
C> Both these equations must produce the same result if no charge 
C> corrections are imposed.
C>
C> In the case Lagrangian multipliers are used to enforce the Gauss
C> theorem the solution for the COSMO charges have an extra term 
C> involving the Lagrangian multiplier. Substituting this solution into
C> \f$E(Q,q,\lambda)\f$ we get
C> \f{eqnarray*}{
C>   E(Q,q,\lambda) &=& \frac{1}{2}Q^TCQ+Q^TBq
C>                   -  \frac{1}{2}qAA^{-1}(BQ+\Lambda) \\\\
C>    &=& \frac{1}{2}Q^TCQ+\frac{1}{2}Q^TBq-\frac{1}{2}\lambda\sum_iq_i \\\\
C>    &=& \frac{1}{2}Q^TCQ+\frac{1}{2}Q^TBq-\frac{1}{2}\lambda Q_{in}
C> \f}
C> Following the same arguments as above we find that the expression
C> starting from Ref.[1] Eq.(11) gets an extra term in the COSMO
C> energy. I.e. at "COSMO contribution Alternative 2" we should have
C> \f{eqnarray*}{
C>   E_{\mathrm{COSMO}}&=&\frac{1}{2}(Q_n^TBq-Q_e^TBq)
C>                      +\frac{1}{2}\lambda Q_{\mathrm{in}}
C> \f}
C> I.e. adding a Lagrangian constraint to Ref.[1] Eq.(8) introduces a
C> shift in Ref.[1] Eq.(11) that equals the product of the Lagrange
C> multiplier and the target value of the charge. This also follows
C> from the interpretation of the Lagrange multiplier as the rate of
C> change of the energy as a function of the constraint value:
C> \f{eqnarray*}{
C>   \lambda &=& \frac{\partial E}{\partial Q_{\mathrm{in}}}
C> \f}
C> This is explained in [8] and it is also commonly used in economics
C> where the constraint values can be chosen delibirately to shift
C> the Lagrangian in the desired direction (the constraints are referred
C> to as "choice sets", see e.g. [9]). This feature also has 
C> implications for methods such a "charge equilibration", "density
C> matrix functional theory", and "constrained DFT".
C>
C> In any case, here the constraint was introduced simply to obtain
C> the correct surface charge. The associated shift in the energy is
C> unintentional and hence when using \f$E(Q,q,\lambda)\f$ to calculate
C> the energy we must correct for this shift ot obtain physically
C> meaningful results (see the subtraction of `elambda` at "Alternative
C> 1").
C>
C> ## Singularity free surface charge self-interaction ##
C>
C> The cavity in continuum solvation models is constructed by creating
C> spherical cavities around all atoms and merging these volumes. The
C> Solvent Accessible Surface (SAS) is created representing the
C> spherical surface around every atom with points and eliminating the 
C> points that fall inside the sphere around a neighboring atom.
C> Klamt et al. [1] suggested representing the surface charge by point
C> charges at the surface discretization points. This leads to
C> singularities in the solvation energy when some of these points come
C> together, typically in the vicinity of the boundary between spheres.
C>
C> The singularities in the surface self-interaction energy need to be
C> addressed to ensure that sensible geometry optimizations are 
C> possible. York and Karplus [4] suggested formally replacing the 
C> surface point charges by Gaussian charge distributions. This leads
C> to an interaction of the form
C> \f{eqnarray*}{
C>   A(r_i,r_j) &=& \frac{\mathrm{erf}(\zeta_{ij}r_{ij})}{r_{ij}}
C> \f}
C> This interaction is relatively easily implemented in force fields
C> where all charges are point charges. Implementing this for the
C> interaction between a point charge and an electron distribution
C> is conceptually harder to do.
C>
C> A consistent implementation would need to use this expression for
C> both the \f$A\f$ and \f$B\f$ matrix. Scalmani and Frisch [5] have
C> sought to do this exactly by explicitly representing the surface
C> charges as Gaussians. This could be done of course using the 
C> charge density fitting integrals but it would require some 
C> engineering as the surface charges would need to be stored in a
C> geometry object. Lange and Herbert [6,7] have followed York and
C> Karplus more closely and applied this approximation only to
C> matrix \f$A\f$ and not to the surface charge-electron interaction.
C> However, even applying this interaction just for \f$A\f$ is 
C> problematic as the expression still contains singularities at
C> \f$r_{ij} = 0\f$ that have to be handled explicitly to avoid floating
C> point exceptions.
C>
C> Considering the use of the potential above in more detail we have
C> \f{eqnarray*}{
C>   A(r_i,r_j) &=& \frac{\mathrm{erf}(\zeta_{ij}r_{ij})}{r_{ij}} \\\\
C>   A(r_i,r_i) &=& \lim_{r_{ij}\to 0}\frac{\mathrm{erf}(\zeta_{ii}r_{ij})}{r_{ij}}(F(r_i))^{-1} \\\\
C>              &=& \zeta_i\sqrt{2/\pi}(F(r_i))^{-1} \\\\
C>   \zeta_{ij} &=& \frac{\zeta_i\zeta_j}{
C>                        \left(\zeta_i^2+\zeta_j^2\right)^{1/2}} \\\\
C>   \zeta_i    &=& \frac{\zeta}{R_I\sqrt{w_i}} \\\\
C>              &=& \frac{\zeta}{\sqrt{|S_i|}} \\\\
C>              &=& \frac{\zeta\sqrt{M}}{R_I\sqrt{2\pi}}
C> \f}
C> where \f$R_I\f$ is the radius of the cavity around atom \f$I\f$,
C> \f$w_i\f$ is the weight of point \f$i\f$ on the unit sphere, hence 
C> \f$R_I\sqrt{w_i}\f$ is equivalent to the surface \f$|S_i|\f$ of the
C> point, \f$\zeta\f$ is a width parameter for the Gaussian distribution
C> that has been optimized to reproduce the Born solvation energy,
C> \f$M\f$ is the number of discretization points on the sphere.
C> In [4] Table 1 it is shown that \f$\zeta\f$ is essentially
C> \f$4.90\f$ for Lebedev grids (for the Boundary Element Mesh we
C> use it is \f$1.00\f$). Finally \f$F_i\f$ is the switching function
C> defined below.
C> From this the limit of two point charges approaching eachother can
C> be established as
C> \f{eqnarray*}{
C>   \lim_{r_j \to r_i} A(r_i,r_j) &=& \zeta_{ij}\frac{2}{\sqrt{\pi}}
C> \f}
C> When two point charges come so close together that \f$r_{ij} < C\f$
C> then this last expression has to be used, otherwise the regular
C> expression for \f$A(r_i,r_j)\f$ should be used.
C>
C> The gradient of this expression is given by
C> \f{eqnarray*}{
C>   \nabla_M A(r_i,r_j)
C>   &=& -\left(\mathrm{erf}(\zeta_{ij}r_{ij})-
C>        \frac{2\zeta_{ij}}{\sqrt{\pi}}
C>        e^{-\zeta_{ij}^2r_{ij}^2}\right)
C>        \frac{\nabla_M r_{ij}}{r_{ij}^2} \\\\
C>   \nabla_M A(r_i,r_i)
C>   &=& -A(r_i,r_i)\sum_B\frac{\partial F}{\partial R_B}\nabla_M r_i
C> \f}
C>
C> ## Continuous switching functions ##
C>
C> In order to obtain a smooth function wrt the nuclear coordinates
C> it is necessary that when the atoms move the surface areas associated
C> with a point charge change smoothly. The approach suggested by 
C> York and Karplus [4] proposes to use multiple radii around an atom. 
C> Each atom has an inner radius \f$R_{in}\f$ and an outer radius
C> \f$R_{out}\f$. The difference between them is the switching radius 
C> \f$R_{sw} = R_{out} - R_{in}\f$. The areas of cavity surface points
C> are multiplied with a weighting factor based on their relative
C> position. In practice the diagonal elements of \f$A\f$ are scaled 
C> by the inverse of the surface areas generating corresponding weights
C> as
C> \f{eqnarray*}{
C>   r_{A,i} - R_B \ge R_{out},&& W_{Ai,B} = 1 \\\\
C>   r_{A,i} - R_B \le R_{in}, && W_{Ai,B} = \infty \\\\
C>   R_{in} < r_{A,i}-R_B < R_{out},
C>      && W_{Ai,B} = 1/f\left(\frac{r_{A,i}-R_B-R_{in}}{R_{SW}}\right)
C> \f}
C> where the function \f$f(r)\f$ is given by
C> \f{eqnarray*}{
C>   f(r) &=& r^3\left(10-15r+6r^2\right) \\\\
C>   \frac{\partial f(r)}{\partial r} &=& 30r^2(r-1)^2
C> \f}
C> The weighting function overall is given by
C> \f{eqnarray*}{
C>    F_{Ai} = \prod_{B\neq A} W_{Ai,B}
C> \f}
C> It can be shown that if \f$F_{Ai}\f$ is close to \f$0\f$ then
C> the corresponding point will not contribute to the energy
C> expression and can be eliminated. This is most easily shown by
C> considering the energy expression and substituting the surface
C> charge expression. If \f$F_{Ai}\f$ goes to 0 then the diagonal
C> element \f$A(r_i,r_i)\f$ approaches infinity. This means that the
C> corresponding row and column in \f$A^{-1}\f$ goes to zero, and
C> the energy expression depends only on \f$A^{-1}\f$ hence the 
C> corresponding point will not contribute to the energy.
C>
C> One interesting observation is that matrix \f$A\f$ is used in
C> linear system of equations \f$Aq = BQ\f$. The condition number of
C> the matrix \f$A\f$ determines the ratio of the relative error in
C> \f$x\f$ and the relative error in \f$BQ\f$. In particular, the 
C> larger the condition number the larger the relative error in \f$x\f$
C> for a given relative error in \f$BQ\f$. The condition number
C> is given by
C> \f{eqnarray*}{
C>    K(A) = \left|\frac{\lambda_{\mathrm{max}}(A)}{
C>                       \lambda_{\mathrm{min}}(A)}\right|
C> \f}
C> where \f$\lambda_{\mathrm{max}}(A)\f$ is the maximum eigenvalue, 
C> and \f$\lambda_{\mathrm{min}}(A)\f$ is the minimal eigenvalue of 
C> \f$A\f$. The approach by York and Karplus for eliminating surface
C> charges is based on raising the condition number to infinity, at
C> which point the vector \f$q\f$ must become very inaccurate. Hence
C> the cutoff for eliminating charges must be chosen carefully to
C> limit the condition number of the remaining part of \f$A\f$.
C>
C> ### References ###
C>
C>   [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>   "COSMO: a new approach to dielectric screening in solvents with
C>    explicit expressions for the screening energy and its gradient",
C>   <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>   <a href="https://doi.org/10.1039/P29930000799">
C>   10.1039/P29930000799</a>.
C>
C>   [2] M.A. Aguilar, F.J. Olivares del Valle, J. Tomasi,
C>   "Nonequilibrium solvation: An ab initio quantummechanical method
C>    in the continuum cavity model approximation",
C>   <i>J. Chem. Phys.</i> (1993) <b>98</b>, pp 7375-7384, DOI:
C>   <a href="https://doi.org/10.1063/1.464728">
C>   10.1063/1.464728</a>.
C>
C>   [3] E.V. Stefanovich, T.N. Truong,
C>   "Optimized atomic radii for quantum dielectric continuum solvation
C>    models", <i>Chem. Phys. Lett.</i> (1995) <b>244</b>, pp 65-74,
C>   DOI:
C>   <a href="https://doi.org/10.1016/0009-2614(95)00898-E">
C>   10.1016/0009-2614(95)00898-E</a>.
C>
C>   [4] D.M. York, M. Karplus,
C>   "A smooth solvation potential based on the conductor-like 
C>    screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>   pp 11060-11079, DOI:
C>   <a href="https://doi.org/10.1021/jp992097l">
C>   10.1021/jp992097l</a>.
C>
C>   [5] G. Scalmani, M.J. Frisch,
C>   "Continuous surface charge polarizable continuum models of
C>    solvation. I. General formalism", <i>J. Chem. Phys.</i> (2010)
C>   <b>132</b>, 114110, DOI:
C>   <a href="https://doi.org/10.1063/1.3359469">
C>   10.1063/1.3359469</a>.
C>
C>   [6] A.W. Lange, J.M. Herbert, 
C>   "Polarizable continuum reaction-field solvation models affording
C>    smooth potential energy surfaces", <i>J. Phys. Chem. Lett.</i>
C>   (2010) <b>1</b>, pp 556-561, DOI:
C>   <a href="https://doi.org/10.1021/jz900282c">
C>   10.1021/jz900282c</a>.
C>
C>   [7] A.W. Lange, J.M. Herbert,
C>   "A smooth, nonsingular, and faithful discretization scheme for
C>    polarizable continuum models: The switching/Gaussian approach",
C>   <i>J. Chem. Phys.</i> (2010) <b>133</b>, 244111, DOI:
C>   <a href="https://doi.org/10.1063/1.3511297">
C>   10.1063/1.3511297</a>.
C>
C>   [8] D. Klein,
C>   "Lagrange multipliers without permanent scarring",
C>   <a href="http://www.cs.berkeley.edu/~klein/papers/lagrange-multipliers.pdf">
C>   http://www.cs.berkeley.edu/~klein/papers/lagrange-multipliers.pdf</a>
C>   [accessed Oct 14, 2014].
C>
C>   [9] P. Milgrom, I. Segal,
C>   "Envelope theorems for arbitrary choice sets",
C>   <i>Econometrica</i> (2002) <b>70</b>, 583-601, DOI:
C>   <a href="https://doi.org/10.1111/1468-0262.00296">
C>   10.1111/1468-0262.00296</a>.
C>
      subroutine hnd_coschg(g_dens,ndens,rtdb,geom,basis,nat,nefc,
     &                      efcc,efcs,efcz,efczz,efciat,ratm,ecos,
     &                      cosmo_file_in)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "bq.fh"
#include "stdio.fh"
#include "prop.fh"
#include "util.fh"
#include "inp.fh"
#include "cosmo_params.fh"
#include "cosmoP.fh"
#include "util_params.fh"
c
      integer ndens         !< [Input] the number of density matrices
      integer g_dens(ndens) !< [Input] the handle for the density
                            !< matrices
      integer rtdb          !< [Input] the RTDB handle
      integer geom          !< [Input] the molecular geometry handle
      integer basis         !< [Input] the molecular basis set handle
      integer nat           !< [Input] the number of atoms
      integer nefc          !< [Input] the number of COSMO charges
c
      double precision efcc(3,nefc) !< [Input] the COSMO charge
                                    !< coordinates
      double precision efcs(nefc)   !< [Input] the COSMO charge
                                    !< surface area
c
      double precision efcz(nefc)   !< [Output] the COSMO charges
      double precision efczz(nefc)  !< [Input[ the COSMO zeta value for
                                    !< each surface charge
      integer          efciat(nefc) !< [Input] the atom associated
                                    !< with each surface charge
      double precision ratm(nat)    !< [Input] the atom radii
      double precision ecos !< [Output] the energy contribution due to
                            !< the COSMO charges
      logical  status
c
      logical  dbug,more,out,direct,noall,all,elec,nucl,iefc_done
      character*16 at_tag
      integer istrlen
      character*255 cosmo_file_in
      character*(nw_max_path_len) cosmo_file
      integer fn
c      integer lineq ! 0: fast direct solver, 1: slow iterative solver
c      integer minbem 
c      integer maxbem 
c      integer do_cosmo_model 
c      integer ificos ! 0 use octahedron, 1 use icosahedron tesselation
c      integer cosmo_sccor ! do correction?
      integer iat ! counter over atoms
      integer jef ! counter over COSMO charges
      integer l_i10, i10
      integer l_i11, i11
      integer l_i12, i12
      integer l_i20, i20
      integer l_i21, i21
      integer l_i22, i22
      integer l_i30, i30
      integer l_i40, i40
      integer l_i50, i50
      integer l_i60, i60
      integer l_i70, i70
      integer l_i80, i80
      integer l_i90, i90
      integer i,ipt,ief ! counters
      integer i_init ! number of ints in memory requirement vector
      integer init ! memory requirement vector
      integer ierr ! error flag
      integer iep ! memory offset of b from Ax=b
      integer ieq ! memory offset of x from Ax=b
      integer nodcmp ! flag specifying how to handle errors
      integer need ! the amount of memory needed
c
      integer l_epot, l_xyzpt, l_zanpt ! memory handles
      integer k_epot, k_xyzpt, k_zanpt ! memory offsets
c
      double precision charge ! the total QM region charge
      double precision chgnuc ! the total nuclear charge
      integer          nelec  ! the total number of electrons
      double precision chgfac ! scale factor for COSMO charges
      double precision chgcos ! the total COSMO surface charge
      double precision chgcvg ! the convergence of the COSMO charges
      double precision chgina ! the inv(A) COSMO charge
      double precision corcos ! the COSMO charge correction
      double precision errcos ! the COSMO charge error
      double precision delchg ! charge difference
      double precision aii,aij,bij,chgief,dij,deta,dum,oldief
      double precision atmefc ! atom - COSMO charge interaction
      double precision efcefc ! COSMO charge - COSMO charge interaction
      double precision elcefc ! electron - COSMO charge interaction
      double precision allefc ! total QM - COSMO charge interaction
      double precision zan ! charges
      double precision xi, xj, xn, xp ! X-coordinates
      double precision yi, yj, yn, yp ! Y-coordinates
      double precision zi, zj, zn, zp ! Z-coordinates
      double precision qi, qj ! charges
      double precision rr ! distance
      double precision solnrg ! solvation energy
      double precision dlambda ! lambda (surface charge correction)
      double precision elambda ! lambda dependent energy term
      double precision ecos2
      double precision pi
      double precision zero, pt5, one, two ! constants
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      double precision zetai, zetaj, zetaij
      logical stat
      logical oprint_energies
c
      logical util_io_unit
      external util_io_unit
c
c MN solvation models --> 
c
      double precision gspol, gstote, estote, gspoldyn, gspolneq
      double precision disp_cosmo_vem
      double precision dvem1, dvem2, dvem3, dvem4, dvem5
      double precision espolgsrf, espol, espolin, espoldyn, delwvem
      double precision wgsrf_cosmo_vem, wstar_cosmo_vem
      double precision wcgsrf_cosmo_vem, w_cosmo_vem, wold_cosmo_vem
      double precision tolw
      logical do_cosmo_smd
      integer do_cosmo_vem, istep_cosmo_vem
c
c     possible values of istep_cosmo_vem: 
c     =0 (initial GS SCF to get GSRF+GS MOs), 
c     =1 (TDDFT + ES density to get ESRF),
c     =2 (calculate expect. values of the new Fock operator containing ESRF and repeat =1)  
c     =3 (when the VEM excitation calculation is converged either exit
c         or proceed to GS emission calculation)
c
      integer iter_cosmo_vem
      integer l_efczfx, k_efczfx  ! memory handles for fixed cosmo-vem charges
      integer g_vem(3) ! ga handles for cosmo-vem GS potential, GS charges, and ES noneq charges
      integer nefc_rtdb
      integer type1, ndim, dim1

      double precision, external :: ddot
c
      oprint_energies = util_print("cosmo_energies",print_never)
      oprint_energies = oprint_energies.and.(ga_nodeid().eq.0)
      
c
      wgsrf_cosmo_vem = zero
      wstar_cosmo_vem = zero
      tolw = 0.00001d0
      iter_cosmo_vem = 0
      if(dielecinf.ne.0d0) then
         dvem1=(dielecinf-1d0)*dielec/dielecinf/(dielec-1d0)
         dvem2=(dielec-dielecinf)/dielecinf/(dielec-1d0)
      else
         dvem1=0d0
         dvem2=0d0
      endif
c
      if (.not. rtdb_get
     $ (rtdb,'cosmo:istep_cosmo_vem',mt_int,1,istep_cosmo_vem))
     $  call errquit('hnd_coschg: cannot get istep_cosmo_vem',
     $  0,rtdb_err)
      if (.not. rtdb_get
     $ (rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $  call errquit('hnd_coschg: cannot get do_cosmo_vem',
     $  0,rtdb_err)
      if (.not. rtdb_get
     $ (rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $  call errquit('hnd_coschg: cannot get do_cosmo_smd',
     $  0,rtdb_err)
c
      if (do_cosmo_vem.ne.0) then
c
c     --- an array for cosmo-vem charges
c
        if (.not.ma_push_get
     $ (mt_dbl,nefc,"hnd_coschg:efczfx",l_efczfx,k_efczfx))
     &   call errquit("hnd_coschg: malloc efczfx failed",0,MA_ERR)
      endif
c
      if (istep_cosmo_vem.eq.2.or.(do_cosmo_vem.eq.2.and.
     $ istep_cosmo_vem.eq.3)) then
c restore charges from GA g_vem(3), pass them through this subroutine
c without change, and save them under "cosmo_bq_efc" along with
c tesserae's geometries for later use in nwdft/scf_dft/dft_fockbld.F 
c
        call ga_sync()
c
        if(.not.rtdb_get(rtdb,'cosmo:g_vem',mt_int,3,g_vem)) then
           if(ga_nodeid().eq.0) write(luout,*)
     W          ' W: failed to read g_vem3: resetting istep_cosmo_vem'
           istep_cosmo_vem=0 
        else
           call nga_get(g_vem(3), 1, nefc, dbl_mb(k_efczfx), nefc)       
c
           call ga_sync()
        endif
      endif
c
c <-- MN solvation models
c
c     ----- number of electrons and charge from rtdb ... -----
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, chgnuc))
     $     call errquit('hnd_coschg: geom_nuc_charge failed',
     $                  0, GEOM_ERR)
      nelec = nint(chgnuc - charge)
      if (nelec .le. 0)
     $     call errquit('hnd_coschg: negative no. of electrons ?',
     $                  nelec, INPUT_ERR)
      if (abs(chgnuc - charge - dble(nelec)) .gt. 1d-8)
     $     call errquit('hnd_coschg: non-integral no. of electrons ?',
     $                  0, INPUT_ERR)
c
      more=.false.
      dbug=.false.
      dbug=dbug.or.more
      out =.false.
      out =out.or.dbug
c
      pi = acos(-1.0d0)
c
      dbug=dbug.and.ga_nodeid().eq.0
      more=more.and.ga_nodeid().eq.0
      out =out .and.ga_nodeid().eq.0
c
      iefc_done=.false.
      all=.false.
      elec=.false.
      nucl=.false.
      ecos=zero
      elambda=zero
      istrlen = 0
c
c     ----- get electrostatic potential at surface points -----
c
c     --- total field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i10",l_i10,i10))
     &   call errquit("hnd_coschg: malloc i10 failed",913,MA_ERR)
c     --- nuclear field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i11",l_i11,i11))
     &   call errquit("hnd_coschg: malloc i11 failed",914,MA_ERR)
c     --- electron field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i12",l_i12,i12))
     &   call errquit("hnd_coschg: malloc i12 failed",914,MA_ERR)
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon_cosmo(basis,geom,g_dens(ndens),efcc,nefc,
     &                      dbl_mb(i12),0)
      call dscal(nefc,-1d0,dbl_mb(i12),1)
c
c     ----- nuclear contribution -----
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_coschg: geom_ncent',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      do iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_coschg: geom_cent_get',911,GEOM_ERR)
      enddo ! iat
c
      do ipt=1,nefc
         xp = efcc(1,ipt)
         yp = efcc(2,ipt)
         zp = efcc(3,ipt)
         dbl_mb(i11+ipt-1) = 0.0d0
         do i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            dbl_mb(i11+ipt-1) = dbl_mb(i11+ipt-1) + zan/rr
         enddo ! i
         dbl_mb(i10+ipt-1) = dbl_mb(i11+ipt-1) + dbl_mb(i12+ipt-1)
      enddo ! ipt
c
c     ----- get surface charges -----
c
c     ----- set up the memory based on the solver -----
c     lineq = 0: fast direct solver, lineq = 1: slow iterative solver
c
      if (lineq.eq.0) then
c
        stat = .true.
        stat = stat .and. ma_push_get(mt_dbl,2*nefc,"hnd_coschg i20",
     &                                l_i20,i20)
        stat = stat .and. ma_push_get(mt_dbl,2*nefc,"hnd_coschg i21",
     &                                l_i21,i21)
        stat = stat .and. ma_push_get(mt_dbl,2*nefc,"hnd_coschg i22",
     &                                l_i22,i22)
        stat = stat .and. ma_push_get(mt_dbl,nefc*nefc,"hnd_coschg i30",
     &                                l_i30,i30)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i40",
     &                                l_i40,i40)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i50",
     &                                l_i50,i50)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i80",
     &                                l_i80,i80)
c
c      check memory
c
       if (.not.stat) then
        call errquit("hnd_coschg: out of memory: lineq = 0 ",950,MA_ERR)
       endif
c
      else if (lineq.eq.1) then
c
        stat = .true.
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i20",
     &                                l_i20,i20)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i21",
     &                                l_i21,i21)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i22",
     &                                l_i22,i22)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i30",
     &                                l_i30,i30)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i40",
     &                                l_i40,i40)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i50",
     &                                l_i50,i50)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i60",
     &                                l_i60,i60)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i70",
     &                                l_i70,i70)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i80",
     &                                l_i80,i80)
c
c      check memory
c
       if (.not.stat) then
        call errquit("hnd_coschg: out of memory: lineq = 1 ",950,MA_ERR)
       endif
c
      else 
       call errquit("hnd_coschg: unknown solver",911,INPUT_ERR)
      end if ! lineq 
c
      if(out) then
         if(lineq.eq.1) then
            write(luout,*) 'use iterative method for -lineq-'
         else
            write(luout,*) 'use in-memory method for -lineq-'
         endif
      endif
c
 10   if (cosmo_sccor.eq.COSMO_SCCOR_NO) then
        all  = .true.
        nucl = .false.
        elec = .false.
      else if (cosmo_sccor.eq.COSMO_SCCOR_LAGRA) then
        all  = .true.
        nucl = .false.
        elec = .false.
      else if (cosmo_sccor.eq.COSMO_SCCOR_SCALE) then
        all  = .false.
        if (.not.nucl.and..not.elec) then
          nucl = .true.
          elec = .false.
        else if (nucl.and..not.elec) then
          nucl = .false.
          elec = .true.
        endif
      endif

c
      if(lineq.eq.0) then
c
c        ===== in memory =====
c
         if(out) then
            write(luout,*) '-lineq- in memory'
         endif
c
c        ----- calculate q* = A^{-1}BQ -----
c
         iep=-99
         ieq=-99
         if (all) then
           iep=i10
           ieq=i20
         else if (nucl) then 
           iep=i11
           ieq=i21
         else if (elec) then
           iep=i12
           ieq=i22
         endif

         call dcopy(nefc,dbl_mb(iep),1,dbl_mb(ieq),1)
         call dfill(nefc,1d0,dbl_mb(ieq+nefc),1)

         call hnd_cosmata(nat,nefc,efcc,efcs,efczz,efciat,ratm,
     $                    dbl_mb(i30))
c
         nodcmp=0
         call hnd_linequ(dbl_mb(i30),nefc,dbl_mb(ieq),nefc,
     1                   dbl_mb(i40),dbl_mb(i50),ierr,nodcmp)
         if (ierr.ne.0) call errquit("hnd_coschg: hnd_linequ failed",
     &                               ierr,UERR)
c
c        ----- calculate t = A^{-1}1 -----
c
         if (cosmo_sccor.eq.COSMO_SCCOR_LAGRA) then
c
           call dcopy(nefc,dbl_mb(ieq+nefc),1,dbl_mb(i80),1)
c
         endif
c
      else
c
c        ===== iterative =====
c
         if(out) then
            write(luout,*) 'iterative -lineq-'
         endif
c
c        ----- calculate qraw = q* = A^{-1}BQ -----
c
         if (.not.rtdb_get(rtdb,'cosmo:n_raw',mt_int,1,nefc_rtdb))
     c        nefc_rtdb=0

         iep=-99
         ieq=-99
         if (all) then
           iep=i10
           ieq=i20
         else if (nucl) then
           iep=i11
           ieq=i21
         else if (elec) then
           iep=i12
           ieq=i22
         endif

         if((nefc_rtdb.eq.0).or.(nefc_rtdb.ne.nefc)) then
            call dfill(nefc,zero,dbl_mb(ieq),1)
            status = rtdb_delete(rtdb,'cosmo:qraw')
            status = rtdb_delete(rtdb,'cosmo:qrawn')
            status = rtdb_delete(rtdb,'cosmo:qrawe')
         else
         if (all) then
           if (.not.rtdb_get(rtdb,'cosmo:qraw',mt_dbl,nefc,
     d           dbl_mb(i20))) then
c             If no raw total COSMO charges were found initialize x=0
              call dfill(nefc,zero,dbl_mb(i20),1)
           endif
         else if (nucl) then
           if (.not.rtdb_get(rtdb,'cosmo:qrawn',mt_dbl,nefc,
     d           dbl_mb(i21)))then
c             If no raw nuclear COSMO charges were found initialize x=0
              call dfill(nefc,zero,dbl_mb(i21),1)
           endif
        else if (elec) then
           if (.not.rtdb_get(rtdb,'cosmo:qrawe',mt_dbl,nefc,
     d           dbl_mb(i22)))then
c             If no raw electron COSMO charges were found initialize x=0
              call dfill(nefc,zero,dbl_mb(i22),1)
           endif
           endif
        endif
c
         direct=.true.
c
c        ----- solve ... -----
c

         call hnd_cg(nat,dbl_mb(iep),dbl_mb(ieq),nefc,
     1        dbl_mb(i40),dbl_mb(i50),dbl_mb(i60),
     d        dbl_mb(i70),
     2                   efcc,efcs,efczz,efciat,ratm)
c
         if (.not.rtdb_put(rtdb,'cosmo:n_raw',mt_int,1,nefc))
     c        call errquit('hnd_coschg: could not store n_raw COSMO ',
     1                  nefc,RTDB_ERR)

         if (all) then
            if (.not.rtdb_put(rtdb,'cosmo:qraw',mt_dbl,nefc,
     d           dbl_mb(i20))) then
             call errquit('hnd_coschg: could not store raw COSMO '
     1                  //'charge',nefc,RTDB_ERR)
           endif
         else if (nucl) then
            if (.not.rtdb_put(rtdb,'cosmo:qrawn',mt_dbl,nefc,
     d           dbl_mb(i21)))then
             call errquit('hnd_coschg: could not store raw nuclear '
     &                  //'COSMO charge',nefc,RTDB_ERR)
           endif
         else if (elec) then
            if (.not.rtdb_put(rtdb,'cosmo:qrawe',mt_dbl,nefc,
     d           dbl_mb(i22)))then
             call errquit('hnd_coschg: could not store raw electron '
     &                  //'COSMO charge',nefc,RTDB_ERR)
           endif
         endif
c
c        ----- calculate t = A^{-1}1 -----
c
         if (cosmo_sccor.eq.COSMO_SCCOR_LAGRA) then
            if (.not.rtdb_get(rtdb,'cosmo:rawt',mt_dbl,nefc,
     d           dbl_mb(i80))) then
c             If no raw A^{-1}1 charges were found initialize x=0
              call dfill(nefc,zero,dbl_mb(i80),1)
           endif
c
           direct=.true.
c
c          ----- solve ... -----
c
           call dfill(nefc,one,dbl_mb(i30),1)
           call hnd_cg(nat,dbl_mb(i30),dbl_mb(i80),nefc,
     1          dbl_mb(i40),dbl_mb(i50),dbl_mb(i60),
     d          dbl_mb(i70),
     2                     efcc,efcs,efczz,efciat,ratm)
c
           if (.not.rtdb_put(rtdb,'cosmo:rawt',mt_dbl,nefc,dbl_mb(i80)))
     1        call errquit('hnd_coschg: could not store raw A^{-1}1',
     2                     nefc,RTDB_ERR)
         endif
      endif
c
c     ----- correct the COSMO charges ... -----
c
      if (cosmo_sccor.eq.COSMO_SCCOR_SCALE) then
c
c       ----- correct the COSMO surface charge by scaling -----
c
        chgcos=zero
        do ief=1,nefc
          chgief=dbl_mb(ief+ieq-1)
          chgcos=chgcos+chgief
        enddo
        if (all) then
c
c         ----- should not get here: it does not make sense for neutrals
c
          errcos=charge-chgcos
          chgfac=charge/chgcos
        else if (nucl) then
          errcos=chgnuc-chgcos
          chgfac=chgnuc/chgcos
        else if (elec) then
          errcos=-dble(nelec)-chgcos
          chgfac=-dble(nelec)/chgcos
        endif
        call dscal(nefc,chgfac,dbl_mb(ieq),1)
      else if (cosmo_sccor.eq.COSMO_SCCOR_LAGRA) then
c
c       ----- correct the COSMO surface charge using Lagrangian -----
c
        chgcos=zero
        chgina=zero ! inverse of A
        do ief=1,nefc
           chgief=dbl_mb(ief+i20-1)
           chgcos=chgcos+chgief
           chgina=chgina+dbl_mb(i80+ief-1)
        enddo
        errcos=charge-chgcos
        dlambda=errcos/chgina
        chgcos=zero
        call daxpy(nefc,dlambda,dbl_mb(i80),1,dbl_mb(i20),1)
        do ief=1,nefc
           chgcos=chgcos+dbl_mb(ief+i20-1)
        enddo
        elambda = 0.5d0*screen*dlambda*chgcos
      endif
      if (cosmo_sccor.eq.COSMO_SCCOR_SCALE) then
        if (.not.all.and..not.elec) goto 10
      endif
      if (cosmo_sccor.eq.COSMO_SCCOR_SCALE) then
c
c       ----- compute combined COSMO charges (nuclear + elec) -----
c
        do ief=1,nefc
           dbl_mb(ief+i20-1) = dbl_mb(ief+i21-1) + dbl_mb(ief+i22-1)
        enddo
      endif
c
c     ----- charge screening due to the dielectric medium -----
c
c     ----- set screening factor -----
c
      chgfac=screen
c
c     ----- apply screening factor -----
c
      chgcos=zero
      call dscal(nefc,-chgfac,dbl_mb(i20),1)
      do ief=1,nefc
         chgcos=chgcos+dbl_mb(ief+i20-1)
      enddo
c
c     ----- store effective charges in -efcz- ... -----
c           check convergence ...
c
      chgcvg=zero
      do ief=1,nefc
         oldief=efcz(ief)
         chgief=   dbl_mb(ief+i20-1)
         delchg=abs(chgief-oldief)
         if(delchg.gt.chgcvg) then
            chgcvg=delchg
         endif
         efcz(ief)=dbl_mb(ief+i20-1)
      enddo
c
c MN solvation models -->
c
      if (istep_cosmo_vem.eq.2 .and. do_cosmo_vem.ne.0) then
c fixed ES nonequilibrium cosmo-vem charges to be added to the Fock matrix
c later on (excitation VEM)
        call dcopy(nefc,dbl_mb(k_efczfx),1,efcz,1)
      endif
      if (istep_cosmo_vem.eq.3 .and. do_cosmo_vem.eq.2) then
c GS RF containing ES cosmo-vem charges (slow portion) to calculate
c the VEM energy in the case of nonequilibrium emission to the GS
        call dscal(nefc,dvem1,efcz,1)
        call daxpy(nefc,dvem2,dbl_mb(k_efczfx),1,efcz,1)
      endif
c
c <-- MN solvation models
c
      if(all) then
         if(out) then
            write(luout,9987) chgcvg
         endif
      endif

      if(dbug) then
         write(luout,9998)
         do ief=1,nefc
            write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
         enddo
      endif
c
c     ----- calculate energy terms from -cosmo- charges -----
c
c     Below are two formulations of the COSMO energy ecos. When no
c     surface charge corrections are used the results should be 
c     identical. If surface charge corrections are employed using
c     scaling factors or Lagrangian multipliers then additional 
c     correction terms are needed to ensure the expressions agree.
c
      call hnd_cos_energy(nat,nefc,chgfac,efcc,efcs,efcz,efczz,efciat,
     &                    ratm,dbl_mb(k_xyzpt),dbl_mb(k_zanpt),
     &                    dbl_mb(i10),
     &                    allefc,atmefc,elcefc,efcefc)
      solnrg= allefc+efcefc
c     COSMO contribution Alternative 1 (Ref.[1] Eq.(8))
c     ecos  = atmefc+efcefc-elambda
      ecos  = atmefc+efcefc
      if (oprint_energies) then
         write(luout,*)'Alternative 1'
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9995) efcefc
         write(luout,9983) elambda
         write(luout,9988) solnrg
         write(luout,9989) allefc,(-two*efcefc)
         write(luout,9994) ecos  
      endif
c
c     ----- other form of the solvation energy ... -----
c
      allefc = ddot(nefc,efcz,1,dbl_mb(i10),1)
      atmefc = ddot(nefc,efcz,1,dbl_mb(i11),1)
      elcefc = ddot(nefc,efcz,1,dbl_mb(i12),1)
      solnrg= pt5* allefc + elambda
c     COSMO contribution Alternative 2 (Ref.[1] Eq.(11))
c     ecos = pt5*(atmefc-elcefc)+elambda
      if (do_cosmo_smd) then
        ecos = pt5*(atmefc-elcefc)
      else
c        ecos = efcefc + atmefc
        ecos = pt5*(atmefc-elcefc) + elambda
      endif

      if (oprint_energies) then
         write(luout,*)'Alternative 2'
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9989) allefc
         write(luout,9988) solnrg
         write(luout,9994) ecos 
      endif

      if(dbug) then
         write(luout,9998)
         do ief=1,nefc
            write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
         enddo
      endif
c
      if(out) then
         write(luout,9993)
      endif
c
c     ----- save -cosmo- charges and energy to -rtdb- -----
c
      if (.not. rtdb_put(rtdb,'cosmo:energy',mt_dbl,1,ecos))
     &   call errquit('hnd_coschg: rtdb put failed for ecos',911,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     $   call errquit('hnd_coschg: rtdb put failed for efcz',912,
     &       rtdb_err)
c
c     ----- for the time being, save in 'geometry' object -----
c
      if(out) then
         write(luout,*) 'in -hnd_coschg ... for -efc- geom = ',geom
      endif
      status=bq_set(cosmo_bq_efc,nefc,efcz,efcc)
      if (.not.status) then
         call errquit('hnd_coschg: bq_set failed !', 0,
     &       geom_err)
      endif
      status=bq_rtdb_store(rtdb,cosmo_bq_efc)
      if (.not.status) then
         call errquit('hnd_coschg: bq_rtdb_store failed !', 0,
     &       geom_err)
      endif
     
c
c     ----- printing cosmo charges for bq module -----
c
      istrlen = inp_strlen(cosmo_file_in)
      if (istrlen.le.0) cosmo_file_in = "cosmo.xyz" 
      call util_file_name(cosmo_file_in,.false.,.false.,cosmo_file)
      call util_file_name_resolve(cosmo_file,.false.)
      if(ga_nodeid().eq.0) then
        if(.not.util_io_unit(80,90,fn))
     +     call errquit('cannot get free unit', 0,
     +       0)
c
        open(unit=fn,form="formatted",file=cosmo_file)
        if (dbug) then
          write(*,*) "printing cosmo charges for bq module",
     +     cosmo_file
        end if
c
        write(fn,*) nefc
        write(fn,*) "cosmo charges (= -Bq-charge) (coordinates in Angstr
     +om) format: Bq x y z q"
        do ief=1,nefc
           write(fn,*) 
     +      "Bq",
     +      efcc(1,ief)*cau2ang,
     +      efcc(2,ief)*cau2ang,
     +      efcc(3,ief)*cau2ang,
     +     -efcz(  ief)
        end do
        close(fn)
      end if
c
c MN solvation models -->
c
      if (do_cosmo_vem.ne.0.or.do_cosmo_smd) then
      ecos = pt5*(atmefc-elcefc)
      endif   
c
c save V_GS and Q_GS for VEM calculation
c
      call ga_sync()
c
      if (do_cosmo_vem.ne.0.and.istep_cosmo_vem.eq.0) then
         if(.not.rtdb_get(rtdb,'cosmo:g_vem',mt_int,3,g_vem)) 
     c        call hnd_create_gvems(rtdb,nefc,g_vem)
         call ga_check_handle(g_vem(1), 'gvem1')
         call nga_inquire(g_vem(1), type1, ndim, dim1)
      if(dim1.ne.nefc) then
         if(ga_nodeid().eq.0) write(luout,*)
     W        ' W: dims1.ne.nefc. creating g_vems'
         call hnd_destroy_gvems(rtdb,g_vem)
         call hnd_create_gvems(rtdb,nefc,g_vem)
      endif

       call ga_zero(g_vem)
       call nga_put(g_vem(1), 1, nefc, dbl_mb(i10),nefc)
       call ga_sync()
       call nga_put(g_vem(2), 1, nefc, efcz, nefc)
       call ga_sync()
      end if
c
c      save 1/2*V_gs*Q_gs for VEM calculation or for SMD output     
c
      if ((do_cosmo_vem.eq.0.and.do_cosmo_smd).or.
     $ (do_cosmo_vem.ne.0.and.istep_cosmo_vem.eq.0).or.
     $ (do_cosmo_vem.eq.2.and.istep_cosmo_vem.eq.3)) then
        gspol = zero
        do ief=1,nefc
         gspol = gspol + dbl_mb(ief+i10-1) * efcz(ief) * pt5
        end do 
        if(.not.rtdb_put(rtdb,'cosmo:gspol',mt_dbl,1,gspol))
     $  call errquit('hnd_coschg: cannot put gspol in rtdb',
     $  0,rtdb_err)
      endif
c
c      save Gp_gs_neq and Gp_gs_dyn for VEM emission calculation  
c
      if (do_cosmo_vem.eq.2.and.istep_cosmo_vem.eq.3) then
        if(.not.rtdb_get(rtdb,
     $ 'cosmo:espol',mt_dbl,1,espol)) call errquit
     $ ('hnd_coschg: cannot get espol from rtdb',
     $ 0,rtdb_err)
        gspolneq = gspol - espol * dvem2       
        do ief=1,nefc
         gspolneq = gspolneq + dbl_mb(ief+i10-1) * pt5 * dvem2 *
     $    dbl_mb(k_efczfx+ief-1)
        end do
        if(.not.rtdb_put(rtdb,'cosmo:gspolneq',mt_dbl,1,gspolneq))
     $  call errquit('hnd_coschg: rtdb put failed for gspolneq',
     $  0,rtdb_err)
c
        gspoldyn = zero
        do ief=1,nefc
         gspoldyn = gspoldyn + pt5 * dbl_mb(ief+i10-1) * 
     $ (efcz(ief) - dvem2 * dbl_mb(k_efczfx+ief-1))
        end do
        if(.not.rtdb_put(rtdb,'cosmo:gspoldyn',mt_dbl,1,gspoldyn))
     $  call errquit('hnd_coschg: rtdb put failed for gspoldyn',
     $  0,rtdb_err)
      endif
c
      call ga_sync()
c
      if (do_cosmo_vem.ne.0.and.istep_cosmo_vem.eq.1) then
c
c      at this stage the excited-state cosmo charges should have been
c      calculated in tddft_grad_compute_g.F.
c
c      calculate noneq charges from current ES charges and saved GS charges
c
         if (do_cosmo_vem.eq.2) then
            dvem1=1d0
            dvem2=0d0
         endif
        espolgsrf = zero
        espolin = zero
        espoldyn = zero
        delwvem = zero
       if(.not.rtdb_get(rtdb,'cosmo:g_vem',mt_int,3,g_vem))
     $  call errquit('hnd_coschg: cannot get g_vem from rtdb',
     $  0,rtdb_err)
        do ief=1,nefc
         call nga_get(g_vem(1), ief, ief, dvem3, 1)
         call nga_get(g_vem(2), ief, ief, dvem4, 1)
         espolgsrf=espolgsrf+(dbl_mb(ief+i10-1)-pt5*dvem3)*dvem4
         espolin=espolin+(dbl_mb(ief+i10-1)-pt5*dvem3)*dvem2*dvem4
         espoldyn=espoldyn+pt5*dbl_mb(ief+i10-1)*dvem1*efcz(ief)
         delwvem=delwvem+pt5*(dbl_mb(ief+i10-1)-dvem3)*
     $ (efcz(ief)-dvem4)*dvem1
         dvem5=dvem1*efcz(ief)+dvem2*dvem4
         dbl_mb(k_efczfx+ief-1)=dvem5
        end do
c
        call ga_sync()
c
        if (.not. rtdb_get(rtdb, 'dft:gstote', mt_dbl, 1, gstote))
     $ call errquit('hnd_coschg: cannot get gstote from rtdb',
     $ 0,rtdb_err)
        if(.not.rtdb_get(rtdb,'cosmo:gspol',mt_dbl,1,gspol))
     $  call errquit('hnd_coschg: cannot get gspol from rtdb',
     $ 0,rtdb_err)
        if (.not.rtdb_get
     $ (rtdb,'tddft:wgsrf_cosmo_vem',mt_dbl,1,wgsrf_cosmo_vem))
     $  call errquit
     $('hnd_coschg: cannot get wgsrf_cosmo_vem from rtdb',
     $ 0,rtdb_err)
        if(.not.rtdb_get(rtdb,
     $ 'cosmo:wold_cosmo_vem',mt_dbl,1,wold_cosmo_vem)) 
     $ wold_cosmo_vem = zero
        status = rtdb_get
     $ (rtdb,'tddft:wstar_cosmo_vem',mt_dbl,1,wstar_cosmo_vem)
c
        espol=espolin+espoldyn
        wcgsrf_cosmo_vem=wgsrf_cosmo_vem+delwvem
        w_cosmo_vem=wstar_cosmo_vem-delwvem
c
        if (status) then
         estote=gstote+w_cosmo_vem
         if(.not.rtdb_put(rtdb,
     $ 'cosmo:estote',mt_dbl,1,estote)) call errquit
     $ ('hnd_coschg: cannot put estote in rtdb',
     $ 0,rtdb_err)
         if(.not.rtdb_put(rtdb,
     $ 'cosmo:espol',mt_dbl,1,espol)) call errquit
     $ ('hnd_coschg: cannot put espol in rtdb',
     $ 0,rtdb_err)
         if (dabs(w_cosmo_vem-wold_cosmo_vem).le.tolw) then
          istep_cosmo_vem = 3
         endif
        endif
c
        if(.not.rtdb_get(rtdb,
     $ 'cosmo:iter_cosmo_vem',mt_int,1,iter_cosmo_vem)) 
     $  iter_cosmo_vem = 0
        iter_cosmo_vem = iter_cosmo_vem + 1
        if(.not.rtdb_put(rtdb,
     $ 'cosmo:iter_cosmo_vem',mt_int,1,iter_cosmo_vem)) 
     $  call errquit
     $ ('hnd_coschg: cannot put iter_cosmo_vem in rtdb',
     $ 0,rtdb_err)
c
        wold_cosmo_vem = w_cosmo_vem
        if(.not.rtdb_put(rtdb,
     $ 'cosmo:wold_cosmo_vem',mt_dbl,1,wold_cosmo_vem)) call errquit
     $ ('hnd_coschg: cannot put wold_cosmo_vem in rtdb',
     $ 0,rtdb_err)
c
        if(ga_nodeid().eq.0) then
         if (.not.status) then
          write (luout,9960)
          write (luout,9961) iter_cosmo_vem
          write (luout,9962) gstote
          write (luout,9963) gspol,gspol*27.211399d0
          write (luout,9964) (gstote+wgsrf_cosmo_vem)
          write (luout,9965) espolgsrf,espolgsrf*27.211399d0
          write (luout,9966) wgsrf_cosmo_vem,wgsrf_cosmo_vem*27.211399d0
          write (luout,9967) (gstote+wcgsrf_cosmo_vem)
          write (luout,9968) espol,espol*27.211399d0
          write (luout,9969) espoldyn,espoldyn*27.211399d0
          write (luout,9970) delwvem,delwvem*27.211399d0
          write (luout,9971) 
     $ wcgsrf_cosmo_vem,wcgsrf_cosmo_vem*27.211399d0
         else 
          write (luout,9960)
          write (luout,9961) iter_cosmo_vem
          write (luout,9962) gstote
          write (luout,9963) gspol,gspol*27.211399d0
          write (luout,9964) (gstote+wgsrf_cosmo_vem)
          write (luout,9965) espolgsrf,espolgsrf*27.211399d0
          write (luout,9966) wgsrf_cosmo_vem,wgsrf_cosmo_vem*27.211399d0
          write (luout,9972) estote
          write (luout,9973) espol,espol*27.211399d0
          write (luout,9969) espoldyn,espoldyn*27.211399d0
          write (luout,9970) delwvem,delwvem*27.211399d0
          write (luout,9974)
     $ w_cosmo_vem,w_cosmo_vem*27.211399d0
         endif
        end if
c
c if do_cosmo_vem=1 g_vem(3)[Q_neq] contains nonequilibrium charges composed of
c slow GS charges and fast ES charges; if do_cosmo_vem=2 Q_neq contains
c full (equilibrium) ES charges 
c
       call ga_zero(g_vem(3))
       call nga_put(g_vem(3), 1, nefc, dbl_mb(k_efczfx),nefc)
c
       call ga_sync()
c
        if (
     $ rtdb_get(rtdb,'cosmo:disp_cosmo_vem',mt_dbl,1,disp_cosmo_vem))
     $    then
          if(ga_nodeid().eq.0) write (luout,9975)
     $ disp_cosmo_vem,disp_cosmo_vem*27.211399d0
        endif
c
        if (iter_cosmo_vem.gt.8) then
          istep_cosmo_vem = 3
          if(ga_nodeid().eq.0) write (luout,9976)
        else
          if (istep_cosmo_vem.eq.3.and.
     $    ga_nodeid().eq.0) write (luout,9977)
        endif
c
      end if
c     end of print Vi and Qi for VEM calculation
c
      if (do_cosmo_vem.ne.0.and.(istep_cosmo_vem.eq.1.or.
     $ istep_cosmo_vem.eq.3)) then
c change istep_cosmo_vem from 1 (TDDFT) to 2 (SCF on fixed MOs)
       if (istep_cosmo_vem.eq.1) istep_cosmo_vem = 2
        if (.not.
     $ rtdb_put(rtdb,'cosmo:istep_cosmo_vem',mt_int,1,istep_cosmo_vem))
     $  call errquit('hnd_coschg: cannot put istep_cosmo_vem in rtdb',
     $ 0,rtdb_err)
      endif
c
c <-- MN solvation models
c
c     ----- release memory block -----
c
      if(.not.ma_chop_stack(l_i10)) call
     &  errquit('hnd_coschg, ma_pop_stack of init failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- charges',
     2       /,10x,15(1h-))
 9998 format(4x,'iefc',6x,'x',5x,6x,'y',5x,6x,'z',5x,5x,'q',4x,
     1     /,1x,53(1h-))
 9997 format(1x,i7,3f12.6,f10.6)
 9995 format(' -efcefc- energy = ',f20.12)
 9994 format(' -ecos  - energy = ',f20.12)
 9993 format(' ...... end of -coschg- ......')
 9992 format(' cosmo potential at -ief = ',i6,f16.10)
 9991 format(' -atmefc- energy = ',f20.12)
 9990 format(' -elcefc- energy = ',f20.12)
 9989 format(' -allefc- energy = ',f20.12,f20.12)
 9988 format(' -solnrg- energy = ',f20.12)
 9987 format(' -cosmo- charges convergence = ',f10.6)
 9986 format(' -wfntyp- = ',a8)
 9985 format(' -scftyp- = ',a8)
 9984 format(' applied screening factor -chgfac- = ',f10.6)
 9983 format(' -lambda- energy = ',f20.12)
c
c MN solvation models -->
c
 9960 format(/,
     $'                          COSMO-VEM solvation results',/,
     $'                          ---------------------------',/
     $' Reference for the VEM model:',/,
     $' Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.;',
     $' Guido, C. A.; Mennucci, B.;',/,' Scalmani, G.; Frisch, M. J.',
     $' Chem. Sci. 2011, 2, 2143',/)
 9961 format(1x,
     $'excitation spectrum data:  GS = initial state, ES = final state',
     $ /,1x,'iteration #',I1)
 9962 format(1x,
     $'(1)          GS equilibrium total free energy = ',
     $ f20.10)
 9963 format(1x,
     $'(2)               GS polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 9964 format(1x,
     $'(3)                 GSRF ES total free energy = ',
     $ f20.10)
 9965 format(1x,
     $'(4)          GSRF ES polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 9966 format(1x,
     $'(5)          GSRF excitation energy (3) - (1) = ',
     $ f20.10,' (',f10.4,' eV)')
 9967 format(1x,
     $'(6)                cGSRF ES total free energy = ',
     $ f20.10)
 9968 format(1x,
     $'(7)         cGSRF ES polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 9969 format(1x,
     $'(8)        fast polarization component of (7) = ',
     $ f20.10,' (',f10.4,' eV)')
 9970 format(1x,
     $'(9)                 1/2 * delV * delQdyn term = ',
     $ f20.10,' (',f10.4,' eV)')
 9971 format(1x,
     $'(10)        cGSRF excitation energy (6) - (1) = ',
     $ f20.10,' (',f10.4,' eV)')
 9972 format(1x,
     $'(6)                  VEM ES total free energy = ',
     $ f20.10)
 9973 format(1x,
     $'(7)           VEM ES polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 9974 format(1x,
     $'(10) VEM vertical excitation energy (6) - (1) = ',
     $ f20.10,' (',f10.4,' eV)')
 9975 format(1x,
     $'(11)    SMSSP solute-solvent dispersion shift = ',
     $ f20.10,' (',f10.4,' eV)')
 9976 format(/,1x,
     $ 'Number of VEM iterations has reached maximum of 9')
 9977 format(/,1x,
     $ 'VEM vertical excitation energy converged')
c
c <-- MN solvation models
c
      end
      subroutine hnd_create_gvems(rtdb,nefc,g_vem)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb
      integer nefc
      integer g_vem(*)
c      
      if(.not.nga_create(mt_dbl, 1, nefc, 'V_GS', nefc, g_vem(1)))
     c     call errquit('creategvems: failed ga_create',0,ga_err)
      if(.not.nga_create(mt_dbl, 1, nefc, 'Q_GS', nefc, g_vem(2)))
     c     call errquit('creategvems: failed ga_create',0,ga_err)
      if(.not.nga_create(mt_dbl, 1, nefc, 'Q_neq', nefc, g_vem(3)))
     c     call errquit('creategvems: failed ga_create',0,ga_err)
      if(.not.rtdb_put(rtdb,'cosmo:g_vem',mt_int,3,g_vem))
     c     call errquit('creategvems: failed rtdbput',0,rtdb_err)
        return
        end
      subroutine hnd_destroy_gvems(rtdb,g_vem)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb
      integer g_vem(*)
c      
          if (.not. ga_destroy(g_vem(1))) call errquit
     $        ('hnd_destr: could not destroy g_vem(1)', 0, ga_err)
          if (.not. ga_destroy(g_vem(2))) call errquit
     $         ('hnd_destr: could not destroy g_vem(2)', 0, ga_err)
          if (.not. ga_destroy(g_vem(3))) call errquit
     $         ('hnd_destr: could not destroy g_vem(3)', 0, ga_err)
          if (.not. rtdb_delete(rtdb,'cosmo:g_vem')) call
     C         errquit(' hnd_destr: rtdb_del failed',0,rtdb_err)
          return
          end
