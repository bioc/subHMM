test_callGenoStates <- function() {

  gs  <- c("A", "AA", "AB","AAA","AAB")
  zk  <- c("0", "A", "AA", "AB","AAA","AAB")
  ct  <- c(0, 1, 2, 2, 3, 3)
  mz  <- c(0, 1, 2, 1, 3, 2)
  pz  <- c(0, 0, 0, 1, 0, 1)
  unq <- c("A", "B")
  nul <- "0"

  ret <- subHMM:::callGenoStates(gs)
  
  
  checkEquals(ret$zk,          zk)
  checkEquals(ret$ctz0,        ct)
  checkEquals(ret$mz,          mz)
  checkEquals(ret$pz,          pz)
  checkEquals(ret$unq_alleles, unq)
  checkEquals(ret$gs.null,     nul)
  
}
