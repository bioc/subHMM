
\name{subhmm_stage2}
\alias{subhmm_stage2}
\title{Stage 2 of the SubHMM Model for Identifying Tumor Subclones}
\description{
  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.                        
}
\usage{
subhmm_stage2(obj, genLen.min=2.5e7, plot=1)
}
\arguments{
  \item{obj}{Th object returned from \code{\link{subhmm}}}
  \item{genLen.min}{Minimum genetic length of the subclone region to be identified. 
                    The default is 2.5e7.}
  \item{plot}{0 or 1 to create plots when the analysis has finished.
             The default is 1.}
}
\details{
  Identifies the subclone regions and produces plots.
}
\value{
  A list with the following names and descriptions: 

  \itemize{
    \item \code{subclone.regions} List of subclone region locations.
    \item \code{subclone.genotype} Subclone genotype corresponding to subclone regions by 
        choosing the largest probability of subclone genotypes for each region according
        to \code{subclone.prob}.
    \item \code{subclone.ind} Subclone region-locus-based indicator.
    \item \code{subclone.prob} Probabilities of subclone genotypes corresponding to 
        each subclone region.
    \item \code{clonal.prop.region} Region-specific clonal proportion estimates.
    \item \code{logR.est} The expected value of \code{logR}
    \item \code{logOR.est} The expected value of \code{logOR}
  } % END: itemize
  Other items in the list are objects needed for creating plots.

}
\references{
  Choo-Wosoba H., Albert P.S., Zhu B.
  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
}
\examples{
  data(data, package="subHMM")

  # Toy example that it runs faster
  op   <- list(stage2=0)
  ret1 <- subhmm(logR, logOR, chrLoc, genoStates=c("", "A", "AA", "AB"), options=op)
  ret2 <- subhmm_stage2(ret1)

}


\author{Hyoyoung Choo-Wosoba, Paul S. Albert, and Bin Zhu <bin.zhu@nih.gov>}
\keyword{ subHMM, ASCNA } 
\seealso{\code{\link{subhmm}}}
