test_subhmm <- function() {

  set.seed(123)
  logR   <- runif(100)
  logOR  <- runif(100)
  chrLoc <- cbind(1, 1:100)
  op     <- list(plot=0, stage2=0, print=0)
  parms  <- c(1.41560274,  0.72360610,  0.07978076,  0.33319020, 98.86430255)

  ret    <- subhmm(logR, logOR, chrLoc, genoStates=c("0", "A", "AA", "AB"), 
                   options=op)

  checkEqualsNumeric(parms, ret$parms, tolerance=1e-4)
}
