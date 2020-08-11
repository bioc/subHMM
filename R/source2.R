subhmm_stage2 <- function(obj, genLen.min=2.5e7, plot=1) {

  lsigma20_lr         <- obj$lsigma20_lr
  v0                  <- obj$v0
  notMissing          <- obj$notMissing
  pzks1               <- obj$prob.stage1
  lr                  <- obj$logR
  lsigma20_lor        <- obj$lsigma20_lor
  logor2              <- obj$logOR2
  mainJ               <- obj$mainJ
  ctzs                <- obj$ctzs
  mz                  <- obj$mz
  pz                  <- obj$pz
  LOG0ARG             <- obj$LOG0ARG
  mz_sub              <- obj$mz_sub
  pz_sub              <- obj$pz_sub
  subc_hidx           <- obj$subclone.ind
  main_hgtype         <- obj$mainclone.genotype
  main_zk             <- obj$genoStates
  mu0_lrs             <- obj$mu0_lrs
  mu0_lors            <- obj$mu0_lors
  lsz0                <- obj$lsz0
  J                   <- obj$J
  lpr0                <- obj$lpr0
  lpsi0               <- obj$lpsi0
  ctms                <- obj$ctms
  sub_zk              <- obj$sub_zk
  tauJ                <- obj$tauJ
  n                   <- length(lr)

  subclone.genotype   <- NA
  col_subidx          <- NULL
  plot.flag           <- 1
  sub_prob            <- NULL
  sub_hgtype_dec      <- NULL
  sub_hgtype_fin      <- NULL 
  sz_ests             <- NULL
  ok_subzoneh_spl     <- NULL
  hat_logr            <- NULL
  hat_logor           <- NULL 

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

  logL_reg <- function(parm, sigma20_lr_reg, lr_reg, pzks1_reg, tmp1_reg, missFlag, noMiss, tmp2_reg, vec1){
    #if (DEBUG) print("BEGIN: logL_reg")
    mu0_lrs    <- mu0_lr(J,lpsi0, lpr0, parm)
    mu0_lors   <- mu0_lor(J, lpr0, parm)
    sigma20_lr <- sigma20_lr_reg
    n          <- length(lr_reg)
   
    logf <- 0
    tmp1 <- tmp1_reg
    mat  <- matrix(lr_reg, nrow=n, ncol=length(mu0_lrs), byrow=FALSE) - matrix(mu0_lrs, nrow=n, ncol=length(mu0_lrs), byrow=TRUE)
    mat  <- 0.5*(v0+1)*log(0.5*(v0+mat^2/sigma20_lr))
    mat  <- (tmp1 - mat)*pzks1_reg
    if (missFlag) {
      m    <- length(mu0_lors)
      mat2 <- dchisq(matrix(vec1, nrow=n, ncol=m, byrow=FALSE), df=1, 
                     ncp=matrix(mu0_lors^2*tmp2_reg, nrow=n, ncol=m, byrow=TRUE), log=T)
      mat[noMiss, ] <- mat[noMiss, , drop=FALSE] - pzks1_reg[noMiss,, drop=FALSE]*lsigma20_lor + 
                           pzks1_reg[noMiss, , drop=FALSE]*mat2[noMiss,, drop=FALSE]
    }
    logf <- sum(mat)
    #if (DEBUG) print("END: logL_reg")

    return(-logf)  
 
  } # END: logL_reg

  lsz_est <- function(sub_reg) {
    #if (DEBUG) print("BEGIN: lsz_est")

    # Non-changing objects
    n_reg          <- length(sub_reg)
    if (!n_reg) return(NA)
    sigma20_lr_reg <- exp(lsigma20_lr) 
    tmp1_reg       <- -0.5*log(2*pi*sigma20_lr_reg)+lgamma((v0+1)/2)-lgamma(v0/2)+0.5*v0*log(v0/2)
    noMiss         <- notMissing[sub_reg]
    missFlag       <- any(!noMiss) 
    pzks1_reg      <- pzks1[sub_reg, , drop=FALSE]
    lr_reg         <- lr[sub_reg]
    tmp2_reg       <- exp(-lsigma20_lor)
    
    vec1           <- logor2[sub_reg]*tmp2_reg
    vec1[!noMiss]  <- NA

    m_est_s <- optim(lsz0, logL_reg, gr=NULL, 
                     sigma20_lr_reg, lr_reg, pzks1_reg, tmp1_reg, missFlag, noMiss, tmp2_reg, vec1,
                     method="BFGS", control=list(trace=0))
    lsz0    <- m_est_s$par
    #if (DEBUG) print("END: lsz_est")

    return(lsz0)

  } # END: lsz_est

  get_sub_prob <- function() {
    #if (DEBUG) print("BEGIN: get_sub_prob")
    nn    <- length(ok_subzoneh_spl)
    if (!nn) {
      #if (DEBUG) print("END: get_sub_prob")
      return(NULL)
    }
    ret   <- matrix(0, nn, mainJ+1) 
    colnames(ret) <- c("SubRegion", main_zk)
    wlist <- list()
    nzk   <- length(main_zk)
    for (i in 1:nzk) wlist[[i]] <- which(sub_zk %in% main_zk[i])+mainJ
    jvec  <- (mainJ+1):J
    for(ss in 1:length(ok_subzoneh_spl)){
      sub_reg <- ok_subzoneh_spl[[ss]]
      mat     <- pzks1[sub_reg, , drop=FALSE]
      denom   <- rowSums(mat[, jvec, drop=FALSE])
      vec     <- rep(NA, nzk)
      for (i in 1:nzk) {
        cols   <- wlist[[i]]
        vec[i] <- mean(rowSums(mat[, cols, drop=FALSE])/denom)
      }
      ret[ss,] <- c(sz_ests[ss,2], vec)
    }
    #if (DEBUG) print("END: get_sub_prob")
    ret

  } # END: get_sub_prob

  if (!all(subc_hidx==0)) {
    # identify subclone regions based on the model
    tmp             <- get_subzone(obj$chrLoc, subc_hidx, main_hgtype, genLen.min)
    ok_subzoneh_spl <- tmp[["ok_subzoneh_spl", exact=TRUE]]
    subc_hidx       <- tmp[["subc_hidx", exact=TRUE]]
    plot.flag       <- 2

    if (!all(subc_hidx==0)) {
      if (length(ok_subzoneh_spl)) {
        sz_ests <- data.frame(Proportion=NA, SubRegion=1:length(ok_subzoneh_spl), 
                              Chr="", Location="", stringsAsFactors=FALSE)
        for(ss in 1:length(ok_subzoneh_spl)){
          lsz1             <- lsz_est(ok_subzoneh_spl[[ss]])
          tmp              <- getChrLocStrs(obj$chrLoc, ok_subzoneh_spl[[ss]])
          sz_ests[ss, 1:2] <- c(1/(1+exp(lsz1)),ss) 
          sz_ests[ss, 3:4] <- c(tmp$chr, tmp$loc) 
        }    
      }

      subclone.ind             <- subc_hidx
      subclone.genotype        <- NA

      ##################
      # BEGIN sub_prob #
      ##################
      sub_prob       <- get_sub_prob()

      hat_logr       <- numeric(n)
      hat_logor      <- numeric(n)
      tmp            <- which(subc_hidx==0)
      hat_logr[tmp]  <- mu0_lrs[main_hgtype[tmp]]
      hat_logor[tmp] <- mu0_lors[main_hgtype[tmp]]

      # final copy number profile considering subclone 
      if (!is.null(sub_prob)) {
        #if (DEBUG) print("BEGIN: final copy number")
        col_subidx <- rep("dark blue",n)
        plot.flag  <- 1 
        sub_hgtype_dec <- apply(matrix(sub_prob[,2:(mainJ+1)],ncol=mainJ),1, function(x) which.max(x))
        sub_hgtype_fin <- rep(NA, n)
        if (length(ok_subzoneh_spl)) {
          for(ss in 1:length(ok_subzoneh_spl)){
            sub_reg <- ok_subzoneh_spl[[ss]]
            if (!length(sub_reg)) next
            lszi    <- log((1-sz_ests[ss,1])/sz_ests[ss,1])
            temp0   <- (tauJ*(unique(main_hgtype[sub_reg])-1)+1):(tauJ*(unique(main_hgtype[sub_reg])-1)+tauJ)
            state_i <- temp0[which(sub_zk[temp0]==main_zk[sub_hgtype_dec[ss]])]+mainJ
            if (length(state_i) && is.finite(state_i)) {
              hat_logr[sub_reg]       <- mu0_lr(J, lpsi0, lpr0, lszi)[state_i]
              hat_logor[sub_reg]      <- mu0_lor(J, lpr0, lszi)[state_i]
            }
            sub_hgtype_fin[sub_reg] <- sub_hgtype_dec[ss]
            col_subidx[sub_reg]     <- "green"
          }
          plot.flag <- 3
        }
        subclone.genotype <- sub_hgtype_fin   
        #if (DEBUG) print("END: final copy number")
      } # END: if (!is.null(sub_prob))

    } # END: if (!all(subc_hidx==0))
  } # END: if (!all(subc_hidx==0))

  if (all(subc_hidx==0)) warning("Not all stage 2 objects can be computed")

  ret <- list(subclone.regions=ok_subzoneh_spl,
              subclone.genotype=subclone.genotype,
              subclone.ind=subc_hidx,
              subclone.prob=sub_prob,
              clonal.prop.region=sz_ests,
              logR.est=hat_logr, logOR.est=hat_logor, 
              plot.flag=plot.flag, col_subidx=col_subidx)

  if (plot) try(subhmm_plot(obj, ret), silent=FALSE)

  ret

} # END: subhmm_stage2

get_subzone <- function(chrloc, subc_hidx, main_hgtype, lowerBound) {

  ok_subzoneh_spl <- NULL
  n               <- nrow(chrloc)
  ret_subc_hidx   <- numeric(n)

  subzoneh     <- which(subc_hidx==1)
  if (!length(subzoneh)) return(list(ok_subzoneh_spl=ok_subzoneh_spl, subc_hidx=ret_subc_hidx))
  subzoneh_spl <- split(subzoneh, cumsum(c(1,diff(subzoneh)!=1 | diff(main_hgtype[subzoneh])!=0)))
  if (!length(subzoneh_spl)) return(list(ok_subzoneh_spl=ok_subzoneh_spl, subc_hidx=ret_subc_hidx))

  gen_length <- NULL
  for(i in 1:length(subzoneh_spl)){
    nsubi <- range(subzoneh_spl[[i]])
    if (all(is.finite(nsubi))) {  
      temp0 <- unlist(tapply(chrloc[nsubi[1]:nsubi[2],2], chrloc[nsubi[1]:nsubi[2],1], function(x) range(x)[2]-range(x)[1]))
      temp1 <- matrix(c(i, sum(temp0)), ncol=2)
      gen_length <- rbind(gen_length, temp1)
    }
  }

  if (length(gen_length)) { 
    idx_oksubzone <- gen_length[which(gen_length[,2]>=lowerBound),1]
    if (length(idx_oksubzone)) {
      ok_subzoneh_spl <- sapply(idx_oksubzone, function(x) list(subzoneh_spl[[x]]))
      ret_subc_hidx[unlist(ok_subzoneh_spl)] <- 1
    }
  }

  list(ok_subzoneh_spl=ok_subzoneh_spl, subc_hidx=ret_subc_hidx)

} # END: get_subzone 

getChrLabels <- function(chrloc) {

  chrs  <- unfactor(chrloc[, 1])
  chrs  <- unique(chrs)
  tmp   <- is.finite(as.numeric(chrs))
  chrs1 <- sort(as.numeric(chrs[tmp]))
  chrs2 <- sort(chrs[!tmp])
  ret   <- c(as.character(chrs1), chrs2)

  ret

} # END: getChrLabels

subhmm_plot12 <- function(obj1, obj2) {

  vec          <- obj1$logR
  n            <- length(vec)
  mainJ        <- obj1$mainJ
  chrloc_label <- as.vector(tapply(1:n, obj1$chrLoc[, 1], function(x) round(median(x))))
  clabels      <- getChrLabels(obj1$chrLoc)

  logR.est <- obj2[["logR.est", exact=TRUE]]
  if (!length(logR.est)) logR.est <- obj1$logR.est
  logOR.est <- obj2[["logOR.est", exact=TRUE]]
  if (!length(logOR.est)) logOR.est <- obj1$logOR.est

  par(mfrow=c(2,1))
  plot(1:n, vec, pch=".", col="dark blue", ylab="logR",xlab="",cex=0.5,xaxt='n')
  if (length(logR.est)) points(1:n, logR.est, col="red",pch=".")
  axis(side=1, cex.axis=0.7,at=chrloc_label, labels=clabels)
  
  vec <- obj1$logOR
  plot(1:n, vec, pch=".",cex=0.5,xaxt='n', col="dark blue", xlab="",ylab="logOR", ylim=c(min(na.omit(vec)), max(na.omit(vec))))
  if (length(logOR.est)) {
    points(1:n,  logOR.est, col="red",pch=".")
    points(1:n, -logOR.est, col="red",pch=".")
  }
  axis(side=1,cex.axis=0.7,at=chrloc_label, labels=clabels)
  
  par(mfrow=c(2,1))
  plot(1:n, obj1$mainclone.genotype, mgp = c(5, 1, 0), pch=".", xaxt='n',ylab="main genotype",
       ylim=c(1,mainJ), xlab="", yaxt="n",col="red")
  axis(2, at=1:mainJ, labels=obj1$genoStates, las=2)
  axis(side=1,cex.axis=0.8, at=chrloc_label, labels=clabels)
  
  plot(1:n, obj2$subclone.ind, pch=".", xaxt='n',ylab="subclone existence",xlab="", yaxt="n",col="red")
  axis(2, at=c(0,1), labels=c("no subclone","subclone"))
  axis(side=1, cex.axis=0.8,at=chrloc_label, labels=clabels)

  NULL

} # END: subhmm_plot12

subhmm_plot3 <- function(obj1, obj2) {

  vec          <- obj1$logR
  n            <- length(vec)
  mainJ        <- obj1$mainJ
  chrloc_label <- as.vector(tapply(1:n, obj1$chrLoc[, 1], function(x) round(median(x))))
  clabels      <- getChrLabels(obj1$chrLoc)

  logR.est <- obj2[["logR.est", exact=TRUE]]
  if (!length(logR.est)) logR.est <- obj1$logR.est
  logOR.est <- obj2[["logOR.est", exact=TRUE]]
  if (!length(logOR.est)) logOR.est <- obj1$logOR.est
 
  par(mfrow=c(2,1))
  plot(1:n, vec, pch=".", col=obj2$col_subidx, ylab="logR",xlab="",cex=0.5,xaxt='n')
  if (length(logR.est)) points(1:n, logR.est, col="red",pch=".")
  axis(side=1, cex.axis=0.7,at=chrloc_label, labels=clabels)

  vec <- obj1$logOR
  plot(1:n, vec, pch=".",cex=0.5,xaxt='n', col=obj2$col_subidx, xlab="",ylab="logOR", ylim=c(min(na.omit(vec)), max(na.omit(vec))))
  if (length(logOR.est)) {
    points(1:n,  logOR.est, col="red",pch=".")
    points(1:n, -logOR.est, col="red",pch=".")
  }
  axis(side=1,cex.axis=0.7,at=chrloc_label, labels=clabels)

  par(mfrow=c(3,1))

  plot(1:n, obj1$mainclone.genotype, mgp = c(5, 1, 0), pch=".", xaxt='n',ylab="main genotype", 
       ylim=c(1,mainJ), xlab="", yaxt="n",col="red")
  axis(2, at=1:mainJ, labels=obj1$genoStates, las=2)
  axis(side=1,cex.axis=0.8, at=chrloc_label, labels=clabels)

  plot(1:n, obj2$subclone.ind, pch=".", xaxt='n',ylab="subclone existence",xlab="", yaxt="n",col="red")
  axis(2, at=c(0,1), labels=c("no subclone","subclone"))
  axis(side=1, cex.axis=0.8,at=chrloc_label, labels=clabels)

  plot(1:n, obj2$subclone.genotype, mgp = c(5, 1, 0),xaxt='n',pch=".", ylab="subclone genotype", 
       yaxt="n", xlab="", ylim=c(1,mainJ),col="red" )
  axis(2, at=1:mainJ, labels=obj1$genoStates, las=2)
  axis(side=1,cex.axis=0.8, at=chrloc_label, labels=clabels)

  NULL

} # END: subhmm_plot3

subhmm_plot <- function(obj1, obj2) {

  if (obj2$plot.flag == 3) {
    subhmm_plot3(obj1, obj2)
  } else {
    subhmm_plot12(obj1, obj2)
  }

  NULL

} # END: subhmm_plot

unfactor <- function(fac, fun=NULL) {

  # fac   Factor
  # fun   Function like as.character or as.numeric, etc

  if (is.factor(fac)) {
    ret <- levels(fac)[fac]
  } else {
    ret <- fac
  }

  if (!is.null(fun)) ret <- fun(ret)

  ret

} # END: unfactor

getCloneSegments <- function(clone.genotype, genoStates, mz, pz, chrloc_thin) {

  ngs                          <- length(genoStates)
  genoStates.mapping           <- data.frame(state=1:ngs, genoState=genoStates, 
                                             nMajor=mz[1:ngs], nMinor=pz[1:ngs], 
                                             stringsAsFactors=FALSE)
  idx = is.na(clone.genotype)
  clone.genotype[idx] = -9
  uchrs = unique(chrloc_thin[, 1])
  nchrs = length(uchrs)

  clone.segs = NULL
  for(j in 1:nchrs) {
    chrj            = uchrs[j]    
    idx             = chrloc_thin[, 1] == chrj
    cur.chrloc_thin = chrloc_thin[idx, , drop=FALSE]

    cur.clone.genotype = clone.genotype[idx]
    cur.clone.genotype.length = length(cur.clone.genotype)
    idx = which(cur.clone.genotype[-1] != cur.clone.genotype[-cur.clone.genotype.length])
    idx.start  = c(1, idx+1)
    idx.end    = c(idx, cur.clone.genotype.length)

    start      = cur.chrloc_thin[idx.start, 2]
    end        = cur.chrloc_thin[idx.end, 2]
    genotype   = cur.clone.genotype[idx.start]
    idx        = match(genotype, genoStates.mapping[, 1])

    segs       = data.frame(chr=chrj, startpos=as.integer(start), 
                    endpos=as.integer(end), nMajor=genoStates.mapping[idx, 'nMajor'], 
                    nMinor = genoStates.mapping[idx, 'nMinor'], check.names=FALSE, 
                    stringsAsFactors=FALSE)
    idx        = is.na(segs[, 'nMajor']) | is.na(segs[, 'nMinor'])
    segs       = segs[!idx, , drop = F]
    clone.segs = rbind(clone.segs, segs)
  }

  clone.segs

}

getChrLocStrs <- function(chrLoc, rows) {

  if (!length(rows)) return(chr="", loc="")
  chrs   <- chrLoc[rows, 1]
  locs   <- as.numeric(chrLoc[rows, 2])
  uchrs  <- sort(unique(chrs))
  nchrs  <- length(uchrs)

  for (i in 1:nchrs) {
    chr  <- uchrs[i]
    tmp  <- chrs %in% chr
    loc2 <- locs[tmp]
    a    <- min(loc2, na.rm=TRUE)
    b    <- max(loc2, na.rm=TRUE)
    lstr <- paste(a, b, sep="-")
    if (i == 1) {
      chr.str <- chr
      loc.str <- lstr
    } else {
      chr.str <- paste(chr.str, chr, sep=",")
      loc.str <- paste(loc.str, lstr, sep=",")
    }
  }

  list(chr=chr.str, loc=loc.str)  

} # END: getChrLocStrs
