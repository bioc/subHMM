test_compute_hatp0 <- function() {

  J        <- 4
  mainJ    <- 2
  tauJ     <- 1
  hatp0_M  <- matrix(c(rep(c(1-tauJ/10,rep(1/10, J/mainJ)),tauJ),1-tauJ/10), J/mainJ, J/mainJ)
  hatp0_m  <- c(0.9, 0.1, 0.1, 0.9)
  hatp0_t0 <- c(rep(1/tauJ,tauJ), rep(1/tauJ,tauJ*tauJ))
  vec      <- c(0.81, 0.09, 0.09, 0.01, 
              0.09, 0.81, 0.01, 0.09, 
              0.09, 0.01, 0.81, 0.09,  
              0.01, 0.09, 0.09, 0.81)
  mat      <- matrix(data=vec, nrow=4, ncol=4, byrow=TRUE)

  ret <- subHMM:::compute_hatp0(hatp0_M, hatp0_m, hatp0_t0, mainJ, tauJ, J)
  checkEqualsNumeric(mat, ret, tolerance=1e-8)

}
