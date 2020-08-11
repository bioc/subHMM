subhmm <- function(logR, logOR, chrLoc, purity=0.8, ploidy=2, clonal.prop=0.8, 
                    logR.var=0.5, logOR.var=0.5, df=3,
                    genoStates=c("0", "A", "AA", "AB","AAA","AAB","AAAA","AABB",
                                 "AAAB","AAAAA","AAABB","AAAAB"), prob0=NULL, 
                    mainclone.trans=NULL, subclone.trans=NULL,
                    subclone.prob=c(0.5, 0.5), subclone.prob.first=NULL,
                    options=NULL) {

  options       <- checkOptions(options)
  maxiter       <- options$maxiter
  df.min        <- options$df.min
  df.max        <- options$df.max
  loglike.eps   <- options$loglike.eps
  parm.eps      <- options$parm.eps
  print         <- options$print

  if (length(logR) != length(logOR)) stop("ERROR: logR and logOR must have the same length")
  if ((purity >= 1) || (purity <= 0)) stop("ERROR: purity must be betwen 0 and 1")
  if (ploidy <= 0) stop("ERROR: ploidy must be positive")
  if ((clonal.prop <= 0) || (clonal.prop >= 1)) stop("ERROR: clonal.prop must be positive")
  if (logR.var <= 0) stop("ERROR: logR.var must be positive")
  if (logOR.var <= 0) stop("ERROR: logOR.var must be positive")
  if (df <= 0) stop("ERROR: df must be positive")

  main_zk                     <- genoStates
  tmp                         <- callGenoStates(main_zk)
  main_zk                     <- tmp$zk
  alleles                     <- tmp$unq_alleles
  gs.null                     <- tmp$gs.null
  names(main_zk)              <- NULL
  mainJ                       <- length(main_zk)
  if (mainJ < 2) stop("ERROR: too few genotype states")
  tauJ                        <- mainJ - 1
  allzk_main                  <- c(main_zk, rep(main_zk,each=tauJ))
  sub_zk                      <- c(sapply(1:mainJ, function(x) main_zk[-x]))
  mz                          <- str_count(allzk_main, alleles[1]) # A allele frequency corresponding to z
  mz_sub                      <- str_count(sub_zk, alleles[1])     # A allele of subclone genotype
  pz                          <- str_count(allzk_main, alleles[2]) # B allele frequency corresponding to z
  pz_sub                      <- str_count(sub_zk, alleles[2])     # B allele of subclone genotype
  ctzs                        <- nchar(allzk_main) # main genotype copy number
  ctzs[which(allzk_main==gs.null)] <- 0
  ctms                         <- nchar(sub_zk)  # subclone genotype copy number
  ctms[which(sub_zk==gs.null)] <- 0

  J                            <- length(allzk_main)
  lr                           <- logR
  logor2                       <- logOR^2
  pr0                          <- purity
  lpr0                         <- log((1-pr0)/pr0)
  sigma20_lr                   <- logR.var
  lsigma20_lr                  <- log(sigma20_lr) 
  lpsi0                        <- log(ploidy, base=2)
  lsigma20_lor                 <- log(logOR.var) 
  r0                           <- prob0
  lsz0                         <- log((1-clonal.prop)/clonal.prop)
  P5_LB                        <- log(options$logOR.var.min) 
  v0                           <- df
  V0_LB                        <- df.min 
  if (V0_LB <= 0) stop("ERROR: with df.min")
  if (v0 <= V0_LB) stop("ERROR: with df or df.min")
  V0_UB                       <- df.max 
  if (V0_UB < V0_LB) stop("ERROR: with df.max")
  if (v0 >= V0_UB) stop("ERROR: with df or df.max")
  if (!is.finite(gamma(df.min/2))) warning("df.min may be too small")
  if (!is.finite(gamma(df.max/2))) warning("df.max may be too large")
  if (!is.finite(lpr0)) stop("ERROR with purity")
  if (!is.finite(lpsi0)) stop("ERROR with ploidy")
  if (!is.finite(lsigma20_lr)) stop("ERROR with logR.var")
  if (!is.finite(lsigma20_lor)) stop("ERROR with logOR.var")
  if (!is.finite(lsz0)) stop("ERROR with clonal.prop")
  if (!is.finite(P5_LB)) stop("ERROR: logOR.var.min is too small")
  REPARM <- 1
  if (options$optim.method != "BFGS") REPARM <- 0

  if (is.null(r0)) r0 <- rep(1/mainJ, mainJ)
  if (length(r0) != mainJ) stop("ERROR with prob0")
  

  # Prevent log(0)
  LOG0ARG <- 1e-300
  LOG0    <- log(LOG0ARG)
  DEBUG   <- print > 9
  
  main_geneidx <- c(1:mainJ,rep(1:mainJ, each=tauJ))
  u0           <- subclone.prob
  hatp0_M      <- mainclone.trans
  if (is.null(hatp0_M)) hatp0_M <- matrix(c(rep(c(1-tauJ/5000,rep(1/5000, J/mainJ)),tauJ),1-tauJ/5000), J/mainJ, J/mainJ)
  if (nrow(hatp0_M) != mainJ)  stop("ERROR: with mainclone.trans")
  if (ncol(hatp0_M) != mainJ)  stop("ERROR: with mainclone.trans")

  hatp0_m      <- as.numeric(subclone.trans)
  if (!length(hatp0_m)) hatp0_m <- c(0.999, 0.001, 0.001, 0.999)
  if (length(hatp0_m) != 4) stop("ERROR: with subclone.trans")

  hatp0_t0 <- subclone.prob.first
  if (!length(hatp0_t0)) hatp0_t0 <- c(rep(1/tauJ,tauJ), rep(1/tauJ,tauJ*tauJ))
  if (length(hatp0_t0) != mainJ*tauJ) stop("ERROR: hatp0_t0 has the wrong length")

  # hatp0
  hatp0 <- compute_hatp0(hatp0_M, hatp0_m, hatp0_t0, mainJ, tauJ, J)

  # Check chrLoc
  check_chrLoc(chrLoc, length(lr), options$stage2)

  # Remove missing. logOR can have missing values
  tmp             <- is.finite(lr)
  tmp[is.na(tmp)] <- FALSE
  if (!all(tmp)) {
    lr     <- lr[tmp]
    logor2 <- logor2[tmp] 
    logOR  <- logOR[tmp] 
    chrLoc <- chrLoc[tmp, , drop=FALSE]
  }  
  n <- length(lr)
  if (!n) stop("ERROR: no observations")

  # Replace small values of logor2 with minLogOR2
  tmp <- logor2 < options$logOR2.min
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) logor2[tmp] <- options$logOR2.min
  

  notMissing  <- is.finite(logor2)
  noMissing   <- all(notMissing)
  if (!noMissing) logor2[!notMissing] <- -999999999.0

  
  ###################
  # Local functions #
  ###################
  logL_ <- function(parm, pzks1){

    #if (DEBUG) print("BEGIN: logL_")
    p6 <- parm[6]
    p1 <- parm[1]
    p2 <- parm[2]
    p3 <- parm[3]
    p4 <- parm[4]
    p5 <- parm[5]

    # Reparameterize v0
    if (REPARM) {
      p62 <- p6*p6
      p6  <- V0_LB + (V0_UB - V0_LB)*p62/(1 + p62) 
    }

    mu0_lrs    <- mu0_lr(J, p1, p2, p3)
    mu0_lors   <- mu0_lor(J, p2, p3)
    sigma20_lr <- exp(p4)
    ep5        <- exp(-p5)
    ncp        <- mu0_lors^2*ep5

    # Problem when ep5 is large (Inf), for now set dmat to 0.
    # Perhaps set bounds for parm[5]
    if (is.finite(ep5)) {    
      #dmat       <- dchisq(rep(logor2*ep5, each=J), df=1, ncp=rep(ncp, times=n), log=T)
      vec1             <- rep(NA, n)
      vec1[notMissing] <- logor2[notMissing]
      vec1             <- rep(vec1*ep5, each=J)
      vec2             <- rep(ncp, times=n)
      dmat             <- dchisq(vec1, df=1, ncp=vec2, log=T)
      dmat             <- matrix(dmat, nrow=n, ncol=J, byrow=TRUE)
    } else {
      stop("ERROR: set bounds on parm[5]")
      dmat             <- matrix(0, nrow=n, ncol=J, byrow=TRUE)
    }
    mat  <- (-0.5*log(2*pi*sigma20_lr) + lgamma((p6+1)/2) - lgamma(p6/2) + 0.5*p6*log(p6/2) -
             0.5*(p6+1)*log(0.5*(p6+(lr-matrix(mu0_lrs, nrow=n, ncol=J, byrow=TRUE))^2/sigma20_lr)) )*pzks1
  
    if (noMissing) {
      mat <- mat - pzks1*(p5 - dmat)
    } else {
      mat[notMissing, ] <- mat[notMissing, , drop=FALSE] - pzks1[notMissing, , drop=FALSE]*(p5 - dmat[notMissing, , drop=FALSE])
    }
    logf <- sum(mat)

    #if (DEBUG) print("END: logL_")

    return(-logf)   

  } # END: logL_

  

  akh_all.ca_all <- function(a1d, v0, sigma20_lr, lsigma20_lor, mu0_lors, mu0_lrs, hatp0) {
    if (DEBUG) print("BEGIN: akh_all.ca_all")

    c1              <- 1/sum(a1d) 
    a1h             <- c1*a1d
    #akh_all.new     <- matrix(NA, n,J)
    #akh_all.new[1,] <- a1h
    #ca_all.new      <- numeric(n)
    #ca_all.new[1]   <- c1
    a1h             <- matrix(a1h, nrow=1)

    pt5v0                  <- 0.5*v0
    pt5v0POWpt5v0          <- pt5v0^pt5v0
    v0Plus1Over2           <- (v0+1)/2
    gamma_v0Over2          <- gamma(v0/2)

    tmp1  <- pt5v0POWpt5v0*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2  <- 2*sigma20_lr
    tmp3  <- exp(-lsigma20_lor)
    tmp4  <- mu0_lors^2*tmp3
    ret_code    <- as.integer(-1)
    ret_akh_all <- as.numeric(rep(-9999, n*J))
    ret_ca_all  <- as.numeric(rep(-9999, n))


    tmp <- .C("akh_ca", as.integer(n), as.integer(J), as.integer(mainJ), as.integer(tauJ), 
                   as.numeric(pt5v0), as.numeric(pt5v0POWpt5v0), as.numeric(v0Plus1Over2),
                   as.numeric(gamma_v0Over2), as.numeric(tmp1), as.numeric(tmp2), 
                   as.numeric(tmp3), as.numeric(tmp4), as.numeric(logor2),
                   as.numeric(lr), as.integer(notMissing), as.numeric(mu0_lrs), 
                   as.numeric(a1h), as.numeric(hatp0), as.numeric(c1), 
                   ret_code=ret_code, ret_akh_all=ret_akh_all, 
                   ret_ca_all=ret_ca_all, PACKAGE="subHMM");
    if (tmp$ret_code) stop("ERROR in akh_ca")
    ret_ca_all  <- tmp$ret_ca_all
    ret_akh_all <- matrix(tmp$ret_akh_all, byrow=TRUE, nrow=n, ncol=J)
    if (DEBUG) print("END: akh_all.ca_all")

    return(list(akh_all=ret_akh_all, ca_all=ret_ca_all))
  
  } # END: akh_all.ca_all

  bkh_all.cb_all <- function() {
    if (DEBUG) print("BEGIN: bkh_all.cb_all")
    tmp0 <- 0.5*v0
    tmp1 <- (tmp0)^(tmp0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2 <- gamma(v0/2)
    tmp3 <- 2*sigma20_lr
    tmp4 <- (v0+1)/2
    tmp5 <- exp(-lsigma20_lor)
    tmp6 <- mu0_lors^2*tmp5
  
    ret_code    <- as.integer(-1)
    bkh_all.new <- as.numeric(rep(-9999, n*J))
    cb_all.new  <- as.numeric(rep(-9999, n))

    tmp <- .C("bkh_cb", as.integer(n), as.integer(J), as.numeric(tmp0), as.numeric(tmp1), 
              as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5),
              as.numeric(tmp6), as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
              as.numeric(mu0_lrs), as.numeric(hatp0),  
              ret_code=ret_code, ret_bkh_all=bkh_all.new, 
              ret_cb_all=cb_all.new, PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR: bad return code in bkh_cb")
    if (DEBUG) print("END: bkh_all.cb_all")
    list(cb_all=tmp$ret_cb_all, bkh_all=matrix(tmp$ret_bkh_all, nrow=n, ncol=J, byrow=TRUE))

  } # END: bkh_all.cb_all



  locFunc_hatp0 <- function(akh_all, lca_all, lcb_all, ca_all, bkh_all, pzks1, 
                  sigma20_lr, lsigma20_lor, mu0_lors, v0, mu0_lrs,
                  hatp0, hatp0_M, hatp0_m, hatp0_t0) {
    if (DEBUG) print("BEGIN: locFunc_hatp0")  
    t1        <- 2*sigma20_lr 
    t2        <- exp(-lsigma20_lor)
    t3        <- sum(akh_all[n,]) 
    tveca     <- reverseCumSum(lca_all)
    tvecb     <- reverseCumSum(lcb_all)

    tvec      <- exp(tveca[3:n]-tvecb[2:(n-1)])*ca_all[2:(n-1)]/t3
    ncpVecJ   <- mu0_lors^2*t2
    tmp0      <- 0.5*v0
    tmp1      <- tmp0^tmp0*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2      <- gamma(v0/2)
    tmp3      <- (v0+1)/2
    
    if (notMissing[n]) {
      vec <- ca_all[n]*tmp1/(tmp2*(tmp0+((lr[n]-mu0_lrs)^2/t1))^tmp3)*t2*dchisq(logor2[n]*t2,df=1, ncp=ncpVecJ)/sum(akh_all[n,])
    } else {
      vec <- ca_all[n]*tmp1/(tmp2*(tmp0+((lr[n]-mu0_lrs)^2/t1))^tmp3)/sum(akh_all[n,])
    }   

    ret_code <- -1
    ret_r0s  <- rep(0, mainJ*tauJ + mainJ)
    tveca    <- NULL
    tvecb    <- NULL
    hatp0_M  <- rep(0, mainJ*mainJ)
    hatp0_m  <- rep(0, 4)
    hatp0_t0 <- rep(0, mainJ*tauJ)
  
    if (DEBUG) print("BEGIN: C code")
    ret <- .C("locFunc_hatp0", as.integer(n), as.integer(J), as.integer(mainJ), 
                 as.integer(tauJ), as.numeric(t1), as.numeric(t2), 
                 as.numeric(tvec), as.numeric(ncpVecJ),
                 as.numeric(tmp0), as.numeric(tmp1), as.numeric(tmp2), as.numeric(tmp3), 
                 as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
                 as.numeric(akh_all), as.numeric(bkh_all), as.numeric(mu0_lrs), 
                 as.numeric(vec), as.numeric(pzks1),
                 ret_code=as.integer(ret_code), ret_r0s=as.numeric(ret_r0s), 
                 hatp0=as.numeric(hatp0), hatp0_M=as.numeric(hatp0_M),
                 hatp0_m=as.numeric(hatp0_m), 
                 hatp0_t0=as.numeric(hatp0_t0), PACKAGE="subHMM")
    if (DEBUG) print("END: C code") 
    if (ret$ret_code) stop("ERROR in locFunc_hatp0")
    
    ret <- list(r0s=ret$ret_r0s, hatp0=matrix(ret$hatp0, nrow=J, ncol=J, byrow=FALSE),
                hatp0_m=ret$hatp0_m, hatp0_t0=ret$hatp0_t0,
                hatp0_M=matrix(ret$hatp0_M, nrow=mainJ, ncol=mainJ, byrow=FALSE))
    if (DEBUG) print("END: locFunc_hatp0")
    return(ret)

  } # END: locFunc_hatp0

  get_a1d <- function(v0, sigma20_lr, mu0_lrs, r0s, lsigma20_lor, mu0_lors) {
    if (DEBUG) print("BEGIN: get_a1d")
    pt5v0                  <- 0.5*v0
    pt5v0POWpt5v0          <- pt5v0^pt5v0
    v0Plus1Over2           <- (v0+1)/2
    gamma_v0Over2          <- gamma(v0/2)
    gamma_v0Plus1Over2     <- gamma(v0Plus1Over2)
    twoPi                  <- 2*pi
    twoPiSimgma20lrPOWmpt5 <- (twoPi*sigma20_lr)^-0.5

    tmp1    <- pt5v0POWpt5v0*gamma_v0Plus1Over2*twoPiSimgma20lrPOWmpt5
    tmp2    <- gamma_v0Over2*(pt5v0+((lr[1]-mu0_lrs)^2/(2*sigma20_lr)))^v0Plus1Over2
    a1d     <- (tmp1/tmp2)*r0s

    if (notMissing[1]) {
      tmp3     <- dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=(mu0_lors)^2*exp(-lsigma20_lor))
      tmp4     <- exp(-lsigma20_lor)*tmp3
      a1d      <- a1d*tmp4
    }

    if (DEBUG) print("END: get_a1d") 
    a1d

  } # END: get_a1d

  get_cb_bkh <- function(v0, sigma20_lr, lsigma20_lor, mu0_lors, mu0_lrs, hatp0) {
    if (DEBUG) print("BEGIN: get_cb_bkh")
    # Scaling HMM with forward-backward alg. in E-step
    tmp0 <- 0.5*v0
    tmp1 <- (tmp0)^(tmp0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2 <- gamma(v0/2)
    tmp3 <- 2*sigma20_lr
    tmp4 <- (v0+1)/2
    tmp5 <- exp(-lsigma20_lor)
    tmp6 <- mu0_lors^2*tmp5

    ret_code    <- as.integer(-1)
    bkh_all.new <- as.numeric(rep(-9999, n*J))
    cb_all.new  <- as.numeric(rep(-9999, n))

    tmp <- .C("bkh_cb", as.integer(n), as.integer(J), as.numeric(tmp0), as.numeric(tmp1), 
              as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5),
              as.numeric(tmp6), as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
              as.numeric(mu0_lrs), as.numeric(hatp0),  
              ret_code=ret_code, ret_bkh_all=bkh_all.new, 
              ret_cb_all=cb_all.new, PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR in bkh_cb")
    if (DEBUG) print("END: get_cb_bkh")
    list(cb_all=tmp$ret_cb_all, bkh_all=matrix(tmp$ret_bkh_all, nrow=n, ncol=J, byrow=TRUE))
    
  } # END: get_cb_bkh 

  

  # These functions are needed to evaluate logL_ function
  mu0_lr <- function(nj, lpsi0, lpr0, lsz0) {
    #if (DEBUG) print("BEGIN: mu0_lr")
    ret      <- rep(NA, nj)
    pr0      <- 1/(1+exp(lpr0))
    sz0      <- 1/(1+exp(lsz0))
    jv1      <- 1:mainJ
    jv2      <- (mainJ+1):nj
    ret[jv1] <- log((2+pr0*(ctzs[jv1]-2)), base=2)-lpsi0
    tmp1     <- ctms[jv2 - mainJ]
    ret[jv2] <- log(2+(tmp1-2)*pr0+(ctzs[jv2]-tmp1)*sz0*pr0, base=2)-lpsi0
    #if (DEBUG) print("END: mu0_lr")
    ret
  
  } # END: mu0_lr

  mu0_lor <- function(nj, lpr0, lsz0) {
    #if (DEBUG) print("BEGIN: mu0_lor")
    ret      <- rep(NA, nj)
    pr0      <- 1/(1+exp(lpr0))
    sz0      <- 1/(1+exp(lsz0))
    jv1      <- 1:mainJ
    jv2      <- (mainJ+1):nj
    arg1     <- mz[jv1]*pr0+(1-pr0)
    arg2     <- pz[jv1]*pr0+(1-pr0)
    tmp      <- arg1 < LOG0ARG
    if (any(tmp)) arg1[tmp] <- LOG0ARG
    tmp      <- arg2 < LOG0ARG
    if (any(tmp)) arg2[tmp] <- LOG0ARG
    ret[jv1] <- log(arg1) - log(arg2)
    tmpj     <- jv2 - mainJ
    tmp1     <- mz_sub[tmpj]
    tmp2     <- pz_sub[tmpj]
    arg1     <- 1+(mz[jv2]-tmp1)*sz0*pr0+(tmp1-1)*pr0   
    arg2     <- 1+(pz[jv2]-tmp2)*sz0*pr0+(tmp2-1)*pr0
    tmp      <- arg1 < LOG0ARG
    if (any(tmp)) arg1[tmp] <- LOG0ARG
    tmp      <- arg2 < LOG0ARG
    if (any(tmp)) arg2[tmp] <- LOG0ARG
    ret[jv2] <- log(arg1) - log(arg2)
    #if (DEBUG) print("END: mu0_lor") 
    ret
  
  } # END: mu0_lor

  hes_logL1 <- function(parm){
    #if (DEBUG) print("BEGIN: hes_logL1")
    mu0_lrs     <- mu0_lr(J, parm[1], parm[2], lsz0)
    mu0_lors    <- mu0_lor(J, parm[2], lsz0)
    ret_code    <- as.integer(-1)
    ret_loglike <- as.numeric(0)

  

    tmp <- .C("C_hes_logL1", as.integer(n), as.integer(J), as.numeric(parm), 
              as.numeric(mu0_lrs), as.numeric(mu0_lors), as.numeric(pi), 
              as.numeric(logor2), as.numeric(lr), as.integer(notMissing), 
              as.numeric(r0s), as.numeric(hatp0),  
              ret_code=ret_code,
              ret_loglike=ret_loglike, PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR in C_hes_logL1")
    #if (DEBUG) print("END: hes_logL1")

    return(tmp$ret_loglike)

  } # END: hes_logL1

  myoptimC <- function(parm, pzks1) {
    
    # From reparameterizing v0 in likelihood
    p6      <- parm[6]
    parm[6] <- sqrt((p6 - V0_LB)/(V0_UB - p6)) 

    # From reparameterizing logOR.var in likelihood
    parm[5] <- sqrt(parm[5] - P5_LB) 

    ret_code <- as.integer(-1)
    nparm    <- as.integer(length(parm))
    reltol   <- sqrt(.Machine$double.eps)
    tmp <- .C("C_myoptimC", as.numeric(parm), nparm, as.numeric(P5_LB),
              as.numeric(V0_LB), as.numeric(V0_UB), 
              as.integer(n), as.integer(J), as.integer(mainJ), as.numeric(lr), as.numeric(logor2),
              as.integer(notMissing), as.numeric(pzks1), as.numeric(ctzs), as.numeric(ctms), 
              as.numeric(mz), as.numeric(pz), as.numeric(mz_sub), as.numeric(pz_sub), 
              as.numeric(reltol), ret_code=ret_code, PACKAGE="subHMM")

    if (tmp$ret_code) stop("ERROR in C_myoptimC")
    ret    <- tmp[[1]]
    if (!all(is.finite(ret))) stop("ERROR with optim: non-finite estimated parameters")
    p6     <- ret[6]
    p62    <- p6*p6
    ret[6] <- V0_LB + (V0_UB - V0_LB)*(p62/(1 + p62)) # From rep
    p5     <- ret[5]
    ret[5] <- P5_LB + p5*p5 
    

    return(ret)

  } # END: myoptimC

  myoptimR <- function(parm, pzks1) {

    tmp <- optim(parm, logL_, gr=NULL, pzks1, method="L-BFGS-B", control=list(trace=0), 
                 lower=c(rep(-Inf, 4), P5_LB, V0_LB), upper=c(rep(Inf,5), V0_UB))
    if (tmp$convergence) stop("ERROR: optim did not converge")
    ret <- tmp$par   

    if (!all(is.finite(ret))) stop("ERROR with optim: non-finite estimated parameters")  

    ret

  } # END: myoptimR

  myoptim <- function(parm, pzks1) {

    if (DEBUG) print("BEGIN: optim")
    if (options$optim.method == "BFGS") {
      ret <- myoptimC(parm, pzks1)
    } else {
      ret <- myoptimR(parm, pzks1)
    }

    

    if (DEBUG) print("END: optim")
    ret

  } # END: myoptim

  compute_hes <- function(theta1) {

    if (DEBUG) print("BEGIN: hessian")

    return(hessian(hes_logL1, theta1))

    hes_logL_theta0 <- function(x) { hes_logL1(x + theta1) }
    h          <- 0.1*abs(theta1) + 1e-4*(abs(theta1) < 1.78e-5)
    hes_theta1 <- hessian(hes_logL1, theta1)
    var_theta1 <- try(solve(-hes_theta1), silent=TRUE)
    vec        <- diag(var_theta1)
    if (any(vec <= 0)) {
      margs      <- list(eps=1e-4, r=2)
      hes_theta1 <- hessian(hes_logL_theta0, rep(0, length(theta1)), method.args=margs)
    }

    if (DEBUG) print("END: hessian")

    hes_theta1

  } # END: compute_hes

  compute_cov <- function() {

    prm1          <- getParmNames(1)
    prm0          <- getParmNames(0)
    theta1        <- c(lpsi0,lpr0,lsigma20_lr,lsigma20_lor,v0)
    names(theta1) <- prm1
    se_theta      <- NULL
    asymse_varlr  <- NULL
    hes_theta1    <- compute_hes(theta1)
    var_theta1    <- try(solve(-hes_theta1), silent=TRUE)
    if ("try-error" %in% class(var_theta1)) {
      warning("Hessian matrix could not be inverted, covariances will not be computed")
      var_theta1  <- NULL
      var_theta   <- NULL
      asymvar1_lr <- NULL
    } else {
      theta1_dev   <- diag(c(2^lpsi0*log(2), -(1+exp(lpr0))^-2*exp(lpr0),
                           exp(lsigma20_lr),exp(lsigma20_lor),1))
      var_theta    <- t(theta1_dev)%*%var_theta1%*%theta1_dev
      se_theta     <- sqrt(diag(var_theta)) # asymtotic standard errors of the global estimators given data and lsz0
      var1_s2lr    <- diag(var_theta)[3]
      var1_s2v     <- diag(var_theta)[5]
  
      if (v0 < 2) {
        warning("df was estimated to be less than two, asymptotic variance of logR will not be computed")
        asymvar1_lr <- NULL 
      } else {
        asymvar1_lr  <- 1/(v0-2)^2*((v0^2*var1_s2lr)-4*v0/(v0-2)*exp(lsigma20_lr)*var_theta[3,5]+4/((v0-2))^2*(exp(lsigma20_lr))^2*var1_s2v)
        asymse_varlr <- sqrt(asymvar1_lr)
      }
      rownames(var_theta) <- prm0
      colnames(var_theta) <- prm0
      names(se_theta)     <- prm0
    }

    list(cov.parms=var_theta, se.parms=se_theta, asym.se.varlogR=asymse_varlr)

  } # END: compute_cov

  EM_alg_C <- function() {

    iargs     <- c(mainJ, tauJ, J, n, maxiter, print, DEBUG)
    reltol    <- sqrt(.Machine$double.eps)
    dargs     <- c(pr0, lpr0, sigma20_lr, lsigma20_lr, lpsi0, logOR.var, lsigma20_lor, lsz0, 
                   P5_LB, V0_LB, V0_UB, v0, loglike.eps, parm.eps, reltol)
    niter     <- -1
    conv      <- 0
    parms     <- c(lpsi0,lpr0,lsz0,lsigma20_lr,lsigma20_lor, v0)
    ret1      <- mu0_lrs
    ret2      <- mu0_lors
    retLL     <- 0
    ret_pzks1 <- rep(-9999, n*J)
    ret_hatp0 <- rep(-9999, J*J)

    tmp    <- .C("C_EM_alg", as.integer(iargs), as.numeric(dargs), as.numeric(lr), as.numeric(logor2), 
                 as.numeric(mz), as.numeric(mz_sub), as.numeric(pz), as.numeric(pz_sub),
                 as.numeric(ctzs), as.numeric(ctms), as.numeric(r0), as.numeric(u0), 
                 hatp0_M=as.numeric(hatp0_M), hatp0_m=as.numeric(hatp0_m), hatp0_t0=as.numeric(hatp0_t0),
                 r0s=as.numeric(r0s), ret_parms=as.numeric(parms), ret_conv=as.integer(conv), 
                 ret_niter=as.integer(niter), ret_mu0_lrs=as.numeric(ret1), 
                 ret_mu0_lors=as.numeric(ret2), ret_loglike=as.numeric(retLL),
                 ret_pzks1=as.numeric(ret_pzks1), 
                 ret_hatp0=as.numeric(ret_hatp0), PACKAGE="subHMM")
    
    list(converged=tmp$ret_conv, parms=tmp$ret_parms, niter=tmp$ret_niter,
         mu0_lrs=tmp$ret_mu0_lrs, mu0_lors=tmp$ret_mu0_lors, r0s=tmp$r0s, 
         hatp0_m=tmp$hatp0_m, hatp0_t0=tmp$hatp0_t0, loglike=tmp$ret_loglike,
         hatp0_M=matrix(tmp$hatp0_M, nrow=mainJ, ncol=mainJ, byrow=FALSE),
         pzks1=matrix(tmp$ret_pzks1, nrow=n, ncol=J, byrow=TRUE),
         hatp0=matrix(tmp$ret_hatp0, nrow=J, ncol=J, byrow=TRUE))

  } # END: EM_alg_C

  EM_alg_R <- function() {

    a1d       <- get_a1d(v0, sigma20_lr, mu0_lrs, r0s, lsigma20_lor, mu0_lors)
    tmp       <- akh_all.ca_all(a1d, v0, sigma20_lr, lsigma20_lor, mu0_lors, mu0_lrs, hatp0)
    ca_all    <- tmp$ca_all
    akh_all   <- tmp$akh_all

    logL_Y    <- rep(NA, maxiter+1)
    logL0     <- sum(-log(ca_all[1:n])) + log(sum(akh_all[n,]))
    logL_Y[1] <- logL0

    if (print) {
      cat(paste("Initial loglikelihood = ", formatValue(logL0), "\n", sep=""))
      if (print > 1) {
        PARM0 <- formatValue(c(ploidy,purity,clonal.prop,logR.var,logOR.var, df))
        names(PARM0) <- getParmNames(0, all=1)
        cat("Initial estimates:\n")
        print(PARM0)
      }
    }

    converged  <- FALSE
    for (mm in 1:maxiter) {    

      # Scaling HMM with forward-backward alg. in E-step
      tmp              <- get_cb_bkh(v0, sigma20_lr, lsigma20_lor, mu0_lors, mu0_lrs, hatp0)
      cb_all           <- tmp$cb_all
      bkh_all          <- tmp$bkh_all

      pzks             <- akh_all*bkh_all/sum(akh_all[n,])
      lca_all          <- log(ca_all)
      lcb_all          <- log(cb_all)
      tmp1             <- reverseCumSum(lca_all)
      tmp2             <- reverseCumSum(lcb_all)

      pzks1            <- matrix(data=NA, nrow=n, ncol=J)
      pzks1[n ,]       <- pzks[n,]/cb_all[n]
      pzks1[1:(n-1), ] <- exp(tmp1[2:n] - tmp2[1:(n-1)])*pzks[1:(n-1), ]

      tmp              <- locFunc_hatp0(akh_all, lca_all, lcb_all, ca_all, bkh_all, pzks1,
                                        sigma20_lr, lsigma20_lor, mu0_lors, v0, mu0_lrs,
                                        hatp0, hatp0_M, hatp0_m, hatp0_t0) 
      r0s              <- tmp$r0s
      hatp0            <- tmp$hatp0
      hatp0_M          <- tmp$hatp0_M
      hatp0_m          <- tmp$hatp0_m
      hatp0_t0         <- tmp$hatp0_t0
      lca_all          <- NULL
      lcb_all          <- NULL
      pzks             <- NULL
      tmp              <- NULL
      tmp1             <- NULL
      tmp2             <- NULL
  
      PARM0            <- c(lpsi0,lpr0,lsz0,lsigma20_lr,lsigma20_lor, v0)
      PARM1            <- myoptim(PARM0, pzks1) # uses pzks
      lpsi0            <- PARM1[1]
      lpr0             <- PARM1[2]
      lsz0             <- PARM1[3]
      lsigma20_lr      <- PARM1[4]
      lsigma20_lor     <- PARM1[5]
      v0               <- PARM1[6]
      sigma20_lr       <- exp(lsigma20_lr)
      pr0              <- 1/(1+exp(lpr0))
      mu0_lrs          <- mu0_lr(J,lpsi0, lpr0, lsz0)
      mu0_lors         <- mu0_lor(J, lpr0, lsz0)

      a1d              <- get_a1d(v0, sigma20_lr, mu0_lrs, r0s, lsigma20_lor, mu0_lors)
  
      c1               <- 1/sum(a1d)
      a1h              <- c1*a1d
      akh_all          <- matrix(NA, n,J)
      akh_all[1,]      <- a1h
      ca_all           <- numeric(n)
      ca_all[1]        <- c1
      tmp              <- akh_all.ca_all(a1d, v0, sigma20_lr, lsigma20_lor, mu0_lors, mu0_lrs, hatp0) 
      ca_all           <- tmp$ca_all
      akh_all          <- tmp$akh_all

      logL1            <- sum(-log(ca_all)) + log(sum(akh_all[n, ]))
      logL_Y[mm+1]     <- logL1
      converged        <- checkStop(logL0, logL1, PARM0, PARM1, loglike.eps, parm.eps, mm, print) 
      if (converged) break
      logL0        <- logL1

    } # END: for (mm in 1:M)

    list(converged=converged, parms=PARM1, niter=mm,
         mu0_lrs=mu0_lrs, mu0_lors=mu0_lors, r0s=r0s, 
         hatp0_m=hatp0_m, hatp0_t0=hatp0_t0, loglike=logL1,
         hatp0_M=hatp0_M, hatp0=hatp0, pzks1=pzks1)


  } # END: EM_alg_R

  checkParms <- function(parms, tol=1e-4) {

    if (abs(parms[5] - P5_LB) < tol) {
      warning("logOR.var has reached its lower bound logOR.var.min")
    }
    p6 <- parms[6]
    if (abs(p6 - df.min) < tol) {
      warning("df has reached its lower bound df.min")
    }
    if (abs(p6 - df.max) < tol) {
      warning("df has reached its upper bound df.max")
    }

    NULL

  } # END: checkParms

  if (print) cat("\nBegin stage 1\n")

  # Compute r0s
  r0s <- c(r0*u0[1], c(sapply(1:mainJ, function(j) sapply(1:tauJ, function(t) r0[j]*u0[2]*hatp0_t0[(j-1)*tauJ+t]))))

  mu0_lrs   <- mu0_lr(J,lpsi0, lpr0, lsz0)
  mu0_lors  <- mu0_lor(J, lpr0, lsz0)


  if (options$codeIndex == 2) {
    tmp <- EM_alg_C()
  } else {
    tmp <- EM_alg_R()
  } 
  akh_all    <- NULL
  ca_all     <- NULL
  a1d        <- NULL  
  converged  <- tmp$converged
  niter      <- tmp$niter
  logL_Y     <- tmp$loglike
  PARM1      <- tmp$parms
  mu0_lrs    <- tmp$mu0_lrs
  mu0_lors   <- tmp$mu0_lors
  r0s        <- tmp$r0s 
  hatp0_m    <- tmp$hatp0_m
  hatp0_t0   <- tmp$hatp0_t0
  hatp0_M    <- tmp$hatp0_M  
  hatp0      <- tmp$hatp0
  pzks1      <- tmp$pzks1
  tmp        <- NULL
  gc()

  if (!converged) {
    warning("Maximum number of iterations reached without convergence")
  } else if (print) {
    cat("Algorithm converged in ", niter, " iterations\n", sep="")
  }
  
  # Check parms
  checkParms(PARM1, tol=1e-4)

  # Obtain estimates from the first step - global estimates
  pr0         <- 1/(1+exp(PARM1[2]))
  psi0        <- 2^PARM1[1]
  sz0         <- 1/(1+exp(PARM1[3]))
  sigma20_lr  <- exp(PARM1[4])
  sigma20_lor <- exp(PARM1[5])
  v0          <- PARM1[6] 

  idx_hgtype        <- max.col(pzks1, ties.method="first")
  subc_hidx         <- ifelse(idx_hgtype >= (mainJ+1), 1, 0)
  main_hgtype       <- numeric(n)
  tmp1              <- subc_hidx == 0
  if (any(tmp1)) main_hgtype[tmp1] <- max.col(pzks1[tmp1, 1:mainJ], ties.method="first")
  tmp1              <- max.col(pzks1[,(mainJ+1):J], ties.method="first")
  tmp2              <- main_geneidx[-c(1:mainJ)]
  tmp3              <- subc_hidx == 1
  if (any(tmp3)) main_hgtype[tmp3] <- tmp2[tmp1[tmp3]]

  hat_logr       <- numeric(n)
  hat_logor      <- numeric(n)
  tmp            <- which(subc_hidx==0)
  hat_logr[tmp]  <- mu0_lrs[main_hgtype[tmp]]
  hat_logor[tmp] <- mu0_lors[main_hgtype[tmp]]

  # Compute variances, asymptotic standard errors
  covlist <- compute_cov()

  theta.orig             <- c(psi0, pr0, sigma20_lr, sigma20_lor, v0)
  names(theta.orig)      <- getParmNames(0)
  hatp0_m                <- matrix(hatp0_m, nrow=2)
  rownames(hatp0_m)      <- c("P(no subclone)", "P(subclone)")
  colnames(hatp0_m)      <- rownames(hatp0_m)
  rownames(hatp0_M)      <- main_zk
  colnames(hatp0_M)      <- main_zk
  tmp <- paste("P(M=", main_zk, "|no S)", sep="")
  for (i in 1:mainJ) tmp <- c(tmp, paste("P(S=", genoStates[-1], "|M=", genoStates[i], ")", sep=""))
  colnames(pzks1)        <- tmp

  # Get clonal segments
  cs  <- try(getCloneSegments(main_hgtype, main_zk, mz, pz, chrLoc), silent=TRUE)
  if ("try-error" %in% class(cs)) cs <- NULL 

  ret <- list(converged=converged, parms=theta.orig, 
              cov.parms=covlist$cov.parms, se.parms=covlist$se.parms, 
              asym.se.varlogR=covlist$asym.se.varlogR,
              logR=lr, logOR=logOR, 
              logR.est=hat_logr, logOR.est=hat_logor,
              prob.stage1=pzks1, mainclone.genotype=main_hgtype, 
              mainclone.segments=cs, 
              clonal.prop.est=sz0, subclone.prob.first=hatp0_t0,
              loglike=logL_Y,
              mainclone.trans.est=hatp0_M, subclone.trans.est=hatp0_m,
              subclone.ind=subc_hidx,
              chrLoc=chrLoc, mainJ=mainJ, genoStates=main_zk,
              lsigma20_lr=lsigma20_lr, v0=v0, notMissing=notMissing,
              lsigma20_lor=lsigma20_lor, ctzs=ctzs, mz=mz, pz=pz,
              LOG0ARG=LOG0ARG, mz_sub=mz_sub, pz_sub=pz_sub, 
              logOR2=logor2, mu0_lrs=mu0_lrs, mu0_lors=mu0_lors,
              lsz0=lsz0, J=J, lpr0=lpr0, lpsi0=lpsi0, ctms=ctms,
              sub_zk=sub_zk, tauJ=tauJ) 

  if (options$stage2) {
    if (print) cat("Begin stage 2\n")  
    ret2 <- subhmm_stage2(ret, genLen.min=options$genLen.min, plot=options$plot)
    tmp  <- names(ret2)
    for (nm in tmp) ret[[nm]] <- ret2[[nm]]
  }

  ret

} # END: subhmm

getParms.origScale <- function(PARM1) {

  ret    <- PARM1
  ret[2] <- 1/(1+exp(PARM1[2]))
  ret[1] <- 2^PARM1[1]
  ret[3] <- 1/(1+exp(PARM1[3]))
  ret[4] <- exp(PARM1[4])
  ret[5] <- exp(PARM1[5])

  ret
 
} # END: getParms.origScale

getParmNames <- function(which, all=0) {

  if (all) {
    if (!which) {
      # Orignal scale
      ret <- c("ploidy", "purity", "clonal.prop", "logR.var", "logOR.var", "df")
    } else {
      ret <- c("logBase2Ploidy", "logOneMinusPurityOverPurity", 
              "logOneMinusClonalpOverClonalp",
              "logLogR.var", "logLogOR.var", "df")
    }
  } else {
    if (!which) {
      # Orignal scale
      ret <- c("ploidy", "purity", "logR.var", "logOR.var", "df")
    } else {
      ret <- c("logBase2Ploidy", "logOneMinusPurityOverPurity", 
              "logLogR.var", "logLogOR.var", "df")
    }
  }

  ret

} # END: getParmNames

checkStop <- function(logL0, logL1, parm0, parm1, loglikeTol, parmTol, iter, print) {

  v1  <- abs(logL0 - logL1)
  v2  <- max(abs(parm0 - parm1))
  if (print) {
    str <- paste("Iter ", iter, " loglike = ", formatValue(logL1, digits=4),
                                " loglike diff = ", formatValue(v1),
                                " max parm diff = ", formatValue(v2), "\n", sep="")
    cat(str)
    if (print > 1) {
      cat("Estimates:\n")
      parm1 <- getParms.origScale(parm1) 
      tmp   <- formatValue(parm1)
      names(tmp) <- getParmNames(0, all=1)
      print(tmp)
    }
  }
  ret <- (v1 < loglikeTol) || (v2 < parmTol)

  ret

} # END: checkStop

formatValue <- function(val, digits=4) {
  as.numeric(formatC(val, format="f", digits=digits))
} # END: formatValue

# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function to compute a cummaulative sum in "reverse order"
reverseCumSum <- function(vec) {

  n   <- length(vec)
  ret <- vec[n:1]
  ret <- cumsum(ret)
  ret <- ret[n:1]

  ret
}


# Function to check the genotype states
checkGenoStates <- function(gs, unq) {

  temp <- nchar(gs) > 0
  gs   <- gs[temp]
  gs2  <- gsub(unq[1], "@", gs, fixed=TRUE)
  gs2  <- gsub(unq[2], unq[1], gs2, fixed=TRUE)
  gs2  <- gsub("@", unq[2], gs2, fixed=TRUE)
    
  temp <- gs2 %in% gs
  if (any(temp)) {
    err <- gs2[temp]
    str <- paste(err, collapse=",", sep="")
    str <- paste("ERROR: with genoStates ", str, sep="")
    stop(str)
  }

  NULL

} # END: checkGenoStates

# Function to compute copy number, etc
callGenoStates <- function(zk) {

  # Initially set "0" to "" and make it be first
  zk    <- unique(removeWhiteSpace(zk))
  zk    <- gsub("0", "", zk, fixed=TRUE)
  tmp   <- nchar(zk) > 0
  zk    <- zk[tmp]
  zk    <- c("", zk)
  J     <- length(zk)
  ctz0  <- nchar(zk)
  tlist <- strsplit(zk, "", fixed=TRUE)
  vec   <- unlist(tlist)
  unq   <- unique(vec)
  temp  <- nchar(unq) > 0
  unq   <- unq[temp]
  nunq  <- length(unq)
  if ((!nunq) || (nunq > 2)) stop("ERROR with genoStates")
  if (nunq == 1) unq <- c(unq, unq)
  mz    <- rep(NA, J)
  pz    <- rep(0, J)

  # Determine the major/minor allele. We want the major to be first
  n1 <- sum(vec %in% unq[1])
  n2 <- sum(vec %in% unq[2])
  if (n2 > n1) unq <- unq[2:1]

  # Compute mz, pz and also normalize the genotype states
  for (j in 1:J) {
    vec   <- tlist[[j]]
    mz[j] <- sum(vec %in% unq[1]) 
    if (nunq > 1) pz[j] <- sum(vec %in% unq[2])
    zk[j] <- paste(sort(vec), collapse="", sep="") 
  } 
  if (length(unique(unq)) > 1) checkGenoStates(zk, unq) 

  # Reset "" to "0"
  zk[1] <- "0"

  list(zk=zk, ctz0=ctz0, mz=mz, pz=pz, unq_alleles=unq, gs.null="0")

} # END: callGenoStates

compute_hatp0 <- function(hatp0_M, hatp0_m, hatp0_t0, mainJ, tauJ, J) {

  hatp01 <- NULL
  for(h in 1:mainJ){
    hatp_0s <- numeric(mainJ*tauJ)
    for(hh in 1:mainJ){
      hatp_0s[((hh-1)*tauJ+1):(hh*tauJ)] <- hatp0_M[h,hh]*hatp0_m[2]*hatp0_t0[((hh-1)*tauJ+1):(hh*tauJ)]
    }
    hatp01 <- rbind(hatp01, hatp_0s)
  }

  hatp0 <- matrix(0,J,J)
  for(i in 1:mainJ){
    hatp0[i,] <- c(hatp0_M[i,]*hatp0_m[1], hatp01[i,] )
  } 

  hatp01 <- NULL
  hatp11 <- NULL
  pt     <- diag(tauJ)
  for(h in 1:mainJ){
    hatp_11s <- matrix(0,tauJ,mainJ*tauJ)
    hatp_11s[,((h-1)*tauJ+1):(h*tauJ)] <- hatp0_M[h,h]*hatp0_m[4]*pt
    for(hh in c(1:mainJ)[-h]){
      hatp_11s[,((hh-1)*tauJ+1):(hh*tauJ)] <- hatp0_M[h,hh]*hatp0_m[4]*matrix(c(hatp0_t0[((hh-1)*tauJ+1):(hh*tauJ)]), byrow=T, nrow=tauJ, ncol=tauJ)
    }
    hatp11 <- rbind(hatp11, hatp_11s)
  }
  pt <- NULL
  hatp0[(mainJ+1):J,(mainJ+1):J] <- hatp11
  
  hatp11 <- NULL
  hatp10 <- NULL
  for(i in 1:mainJ){
    temp <- matrix(c(hatp0_M[i,]*hatp0_m[3]), nrow=1)
    temp1 <- matrix(rep(temp,tauJ), byrow=T, ncol=mainJ)
    hatp10 <- rbind(hatp10, temp1)
  }
  hatp0[(mainJ+1):J,1:mainJ] <- hatp10
  
  hatp0

} # END: compute_hatp0

check_chrLoc <- function(chrLoc, n, stage2) {

  if (length(chrLoc)) {
    nr <- nrow(chrLoc)
    nc <- ncol(chrLoc)
    if (!length(nr)) nr <- 0
    if (!length(nc)) nc <- 0
    if (nc != 2) stop("ERROR: chrLoc must be a matrix of data frame with two columns")
    if (nr != n) stop(paste("ERROR: chrLoc must have ", n, " rows", sep=""))
  } else {
    stop("ERROR: chrLoc must be specified")
  }

  NULL

} # END: check_chrLoc

printAbsSum <- function(obj) {
  val <- sum(abs(obj))
  val2 <- formatC(val, format="g", digits=20)
  print(c(val, val2))
}

checkOptions <- function(options) {

  valid <- c("maxiter", "logOR2.min", "logOR.var.min",
             "df.min", "df.max", "loglike.eps", "parm.eps", 
             "genLen.min", "print", "plot", "optim.method", "stage2",
             "codeIndex")
  def   <- list(1000, 1e-6, 1e-5,
             1e-5, 100, 0.01, 0, 
             2.5e7, 1, 1, "BFGS", 1,
             2)

  options    <- default.list(options, valid, def)

  nm    <- names(options)
  tmp   <- !(nm %in% valid)
  err   <- nm[tmp]
  if (length(err)) {
    tmp <- paste(err, collapse=", ", sep="")
    stop(paste("ERROR: ", tmp, " are not valid options"))
  }
  
  if (options$maxiter < 0) stop("ERROR: with options$maxiter")
  if (options$genLen.min < 1) stop("ERROR: with options$genLen.min")
  if (options$logOR2.min < 0) stop("ERROR: with options$logOR2.min")
  if (options$logOR.var.min < 0) stop("ERROR: with options$logOR.var.min")
  if ((options$loglike.eps < 0) && (options$parm.eps < 0)) stop("ERROR: with loglike.eps and parm.eps")
  if (options$df.min < 0) stop("ERROR: with options$df.min")
  if (options$df.max < 0) stop("ERROR: with options$df.max")
  if (options$df.max < options$df.min) stop("ERROR: with options$df.max")
  if (!(options$optim.method %in% c("BFGS", "L-BFGS-B"))) stop("ERROR: with options$optim.method")
  if (options$optim.method != "BFGS") options$codeIndex <- 0

  options
  
} # END: checkOptions

default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

