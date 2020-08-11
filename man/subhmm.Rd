
\name{subhmm}
\alias{subhmm}
\title{The SubHMM Model for Identifying Tumor Subclones}
\description{
  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.                        
}
\usage{
subhmm(logR, logOR, chrLoc, purity=0.8, ploidy=2, clonal.prop=0.8, 
       logR.var=0.5, logOR.var=0.5, df=3,
       genoStates=c("0", "A", "AA", "AB","AAA","AAB","AAAA",
                    "AABB", "AAAB","AAAAA","AAABB","AAAAB"), 
                    prob0=NULL,mainclone.trans=NULL, subclone.trans=NULL,
                    subclone.prob=c(0.5, 0.5), subclone.prob.first=NULL,
                    options=NULL)
}
\arguments{
  \item{logR}{ Vector of (non-missing) log(Ratio) values. No default }
  \item{logOR}{ Vector of log(Odds Ratio) values. Note that this vector can have
               missing values and will be squared inside the function. No default. }
  \item{chrLoc}{ Matrix or data frame with two columns, where column 1 is the chromosome
                 and column 2 is the location for each genetic marker with
                 log(Ratio) values given in \code{logR}. No default. }
  \item{purity}{ Initial value for the tumor purity. The default is 0.8.}
  \item{ploidy}{ Initial value for the ploidy. The default is 2.}
  \item{clonal.prop}{ Initial value for the constant clonal proportion parameter
                     at the first stage (see details). The default is 0.5.}
  \item{logR.var}{ Initial value for the variance component of logR. The default is 0.5.}
  \item{logOR.var}{Initial value for the variance component of logOR. The default is 0.5.}
  \item{df}{Initial value for the degrees of freedom. The default is 3.}
  \item{genoStates}{Character vector of hidden genotype states. The default is \cr
     c("0", "A", "AA", "AB", "AAA", "AAB", "AAAA", \cr
       "AABB", "AAAB", "AAAAA", "AAABB", "AAAAB").
    The state "0" denotes homozygous deletion, and will always be included.}
  \item{prob0}{NULL or a vector of initial probabilities for \code{genoStates}. 
               If NULL, then it will be set to \code{rep(1/length(genoStates), length(genoStates))}. 
               The default is NULL.}
  \item{mainclone.trans}{NULL or estimated transition matrix of mainclone genotypes. 
       If \code{length(genoStates) = J}, then this matrix has the form: \cr
          \tabular{lccc}{
                            \tab  P(genoStates[1]) \tab ... \tab P(genoStates[J]) \cr
           P(genoStates[1]) \tab  p11              \tab ... \tab p1J              \cr
              ...           \tab  ...              \tab ... \tab ...              \cr     
           P(genoStates[J]) \tab  pJ1              \tab ... \tab pJJ          
          }
               If NULL, then it will be set to \cr
               \code{matrix(c(rep(c(1-(J-1)/5000,rep(1/5000,J)),(J-1)),1-(J-1)/5000),J,J)},
               where \code{J = length(genoStates)}.
               The default is NULL.}
  \item{subclone.trans}{NULL or estimated transition matrix of subclone existence. 
               This matrix has the form: \cr
          \tabular{lcc}{
                          \tab  P(no subclone) \tab P(subclone) \cr
           P(no subclone) \tab  p11            \tab p12         \cr
           P(subclone)    \tab  p21            \tab p22          
          }
               If NULL, then it will be set to 
               \code{matrix(c(0.999, 0.001, 0.001, 0.999), nrow=2)}.
               The default is NULL.}
  \item{subclone.prob}{Vector of initial probability of subclone existence, ie 
                    \code{c(P(no subclone), P(subclone))}.
                  The default is \code{c(0.5, 0.5)}}.

  \item{subclone.prob.first}{Vector of initial probabilities of a subclone genotype at the first
                 locus where the subclone region begins. The default is NULL.}
  \item{options}{List of options (see details). The default is NULL.}
}
\details{
 This model consists of two stages, where stage 1 utilizes an EM algorithm for 
  parameter estimation, and stage 2 identifies the subclone regions. 
 See the reference for details of the SubHMM algorithm. \cr
 Missing values are allowed for \code{logOR} values as long as \code{logR} is observed. 
  To prevent numerical errors in the \code{\link[base]{gamma}} function, \code{df.min} should be 
  a positive value, and \code{df.max} should not be too large. Similarly, a positive
  lower bound (\code{logOR.var.min}) needs to be set for \code{logOR.var} to prevent 
  overflow errors and to prevent \code{\link[stats]{dchisq}} from taking a very long
  time to compute when the non-centrality parameter and argument are large.
  Internally, the algorithm transforms the parameters, so the value of the 
  loglikelihood is based on the transformed parameters. However, parameter estimates 
  that are printed to the console (\code{print} > 1) will be their untransformed 
  values.
  The algorithm terminates when either a likelihood condition defined by
  \code{loglike.eps} or a parameter condition defined by \code{parm.eps} is satisfied. Either 
  condition can be turned off by setting that stopping tolerance to a value less than or
  equal to zero.

  \bold{options}
  \itemize{
    \item \code{maxiter} Maximum number of iterations for the algorithm. The default is 1000.
    \item \code{logOR2.min} Minimum value for \code{logOR^2} to prevent numerical difficulties in
                  the algorithm.
                The default is 1e-6.
    \item \code{logOR.var.min} Minimum value for \code{logOR.var} to prevent numerical difficulties in
                  the algorithm.
                The default is 1e-5.
    \item \code{df.min} Lower bound for the \code{df} parameter in the optimization. 
                The default is 1e-5.
    \item \code{df.max} Upper bound for the \code{df} parameter in the optimization. 
                The default is 100.
    \item \code{loglike.eps} A stopping tolerance for the algorithm (see details). The algorithm will stop
                when two succesive log-likelihood values differ by less than \code{loglike.eps}.
                The default is 0.01.
    %\item \code{parm.eps} A stopping tolerance for the algorithm (see details). The algorithm will stop
    %            when the maximum difference between two succesive sets of estimated 
    %            parameters differ by less than \code{parm.eps}.
    %            The default is 1e-7.
    \item \code{genLen.min} Minimum genetic length of the subclone region to be identified. 
                The default is 2.5e7.
    \item \code{print} Integer to print information to the console. More information is printed with
                a larger value for this option (see details).
             The default is 1.
    \item \code{plot} 0 or 1 to create plots when the analysis has finished.
             The default is 1.
    \item \code{stage2} 0 or 1 to conduct the stage 2 analysis (see details).
             The default is 1.
    \item \code{optim.method} "BFGS" or "L-BFGS-B".
             The default is "BFGS".


  } % END: itemize
}
\value{
  A list with the following names and descriptions: 

  \itemize{
    \item \code{converged} Convergence status
    \item \code{parms} Vector of estimated parameters \code{ploidy}, 
         \code{purity}, \code{logR.var}, \code{logOR.var}, \code{df}.
    \item \code{cov.parms} Estimated covariance matrix of \code{parms}.
    \item \code{se.parms} Estimated standard errors of \code{parms}.
    \item \code{asym.se.varlogR} Estimated asymptotic standard error of 
             the variance of \code{logR}.
    \item \code{logR.est} The expected value of \code{logR}.
    \item \code{logOR.est} The expected value of \code{logOR}.
    \item \code{prob.stage1} Matrix of posterior probabilities for each genotype state 
                              at stage 1. The column names have the form 
        \code{P(M=G|no S) = } the probability of mainclone genotype G without
        subclone existence, and \code{P(S=G1|M=G2) = } the probability of
        subclone genotype G1 given mainclone genotype G2 and subclone existence.
    \item \code{mainclone.genotype} Mainclone genotype index in \code{genoStates}.
    \item \code{mainclone.segments} Mainclone segments.
    \item \code{subclone.prob} Probabilities of subclone genotypes corresponding to 
        each subclone region.
    \item \code{subclone.ind} Subclone region-locus-based indicator.
    \item \code{subclone.genotype} Subclone genotype corresponding to subclone regions by 
        choosing the largest probability of subclone genotypes for each region according
        to \code{subclone.prob}.
    \item \code{subclone.regions} List of subclone region locations.
    \item \code{clonal.prop.region} Region-specific clonal proportion estimates.
    \item \code{clonal.prop.est}  Clonal proportion estimate.
    \item \code{loglike} Final log-likelihood value.
    \item \code{mainclone.trans.est} Estimate of the \code{mainclone.trans} matrix.
    \item \code{subclone.trans.est} Estimate of the \code{subclone.trans} matrix.
    \item \code{subclone.prob.first} Probabilities of a subclone genotype at the first
                 locus where the subclone region begins.
    \item \code{chrLoc} The chromosome and location for each locus.
  } % END: itemize
  Other items in the list are objects needed for the stage 2 analysis.

  \bold{NOTE: The value for \code{logR.var} in the object \code{se.parms} is
    the estimate of the standard error for the variance component of 
    \code{logR} (\code{logR.var}). This is different from the estimated
    standard error of the variance of \code{logR}, which is \code{asym.se.varlogR}.
  }

 
}
\references{
  Choo-Wosoba H., Albert P.S., Zhu B.
  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
}


\author{Hyoyoung Choo-Wosoba, Paul S. Albert, and Bin Zhu <bin.zhu@nih.gov>}
\examples{
  data(data, package="subHMM")

  # Toy example so that it runs faster
  op  <- list(plot=0)
  ret <- subhmm(logR, logOR, chrLoc, genoStates=c("", "A", "AA", "AB"), options=op)
  
}
\keyword{ subHMM, ASCNA } 
\seealso{\code{\link{subhmm_stage2}}}
