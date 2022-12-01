
#' TRACE: fasT and Robust Ancestry Coordinate Estimation
#' 
#' Adapted from http://csg.sph.umich.edu/chaolong/LASER
#' 
#' @param refGeno A matrix of genotypes of reference individuals. Provide probes as rows and samples as columns.
#' @param studyGeno A matrix of genotypes of study samples. Provide probes as rows and samples as columns.
#' @param MIN_LOCI Minimum number of non-missing loci required
#' @param DIM Number of PCs in the reference to match
#' @param DIM_HIGH Number of PCs for sample-specific PCA
#' @param MAX_ITER Maximum iterations for the projection Procrustes analysis
#' @param THRESHOLD Convergence criterion for the projection Procrustes analysis
#' @param cpu Number of CPU.
#' @return  A list containing
#' \item{refPC}{Top PCs in the reference}
#' \item{studyPC}{Top PCs in study samples}
#' @export
trace <- function(refGeno, studyGeno, MIN_LOCI=100, DIM=4, DIM_HIGH=20, MAX_ITER=10000, THRESHOLD=0.000001, cpu=1){
  ## Get intersected probes
  probes <- intersect(rownames(refGeno), rownames(studyGeno))
  data(refGeno_1KGP3_SNP_failQC)
  probes <- probes[!(probes %in% refGeno_1KGP3_SNP_failQC)] # rm probes with HWE<1e-20 or F_MISSING>=0.05 in 1KGP
  refGeno <- refGeno[probes,]
  studyGeno <- studyGeno[probes,]
  if(length(probes) < MIN_LOCI){
    print(paste0("ERROR: Number of loci less than ", MIN_LOCI, ". Exit."))
    return(NA)
  }
  
  ## PCA on reference individuals
  refMeans <- rowMeans(refGeno)
  refSDs <- apply(refGeno, 1, sd)
  refGeno_scale <- (refGeno - matrix(rep(refMeans, ncol(refGeno)), ncol=ncol(refGeno))) / 
    matrix(rep(refSDs, ncol(refGeno)), ncol=ncol(refGeno))
  svd <- svd(t(refGeno_scale)) # "d" "u" "v"
  refPC <- svd$u * matrix(rep(svd$d, each=nrow(svd$u)), nrow=nrow(svd$u))
  rownames(refPC) <- colnames(refGeno)
  colnames(refPC) <- paste0("PC", 1:ncol(refPC))
  refPC <- refPC[, 1:DIM]
  
  ## PCA and Procrustes analysis for each study individuals (one-by-one)
  cl <- makeCluster(cpu)
  registerDoParallel(cl)
  studyPC_list <- foreach(sp=1:ncol(studyGeno), .export=c("procrustes", "pprocrustes")) %dopar% {
    ## PCA on reference and study individuals
    X_scale <- (studyGeno[,sp] - refMeans) / refSDs
    svd <- svd(t(cbind(refGeno_scale, Study=X_scale))) # "d" "u" "v"
    pc <- svd$u * matrix(rep(svd$d, each=nrow(svd$u)), nrow=nrow(svd$u))
    rownames(pc) <- c(colnames(refGeno), colnames(studyGeno)[sp])
    colnames(pc) <- paste0("PC", 1:ncol(pc))
    refPC_new <- pc[1:ncol(refGeno), 1:DIM_HIGH]
    PC_one <- pc[nrow(pc), 1:DIM_HIGH, drop=F]
    
    ## Procrustes analysis
    proc <- pprocrustes(refPC_new, refPC, MAX_ITER=MAX_ITER, THRESHOLD=THRESHOLD)
    if(proc$epsilon > THRESHOLD){
      print("Warning: Projection Procrustes analysis doesn't converge")
    }
    rotPC_one = proc$rho * PC_one %*% proc$A + proc$b
    if(DIM_HIGH > DIM){
      rotPC_one <- rotPC_one[, 1:DIM, drop=F]
    }
    rotPC_one
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  studyPC <- do.call(rbind, studyPC_list)
  colnames(studyPC) <- paste0("PC", 1:DIM)
  list(refPC=refPC, studyPC=studyPC)
}

#' Standard Procrustes Analysis
#' 
#' Adapted from http://csg.sph.umich.edu/chaolong/LASER
#' 
#' @param refPC_new Top PCs in the combination of reference samples and one study sample
#' @param refPC Top PCs in the reference samples
#' @param PROCRUSTES_SCALE Fit the scaling parameter to maximize similarity
#' @return Procrustes Analysis results
#' @export
procrustes <- function(refPC_new, refPC, PROCRUSTES_SCALE=0){
  NUM <- nrow(refPC_new)
  # Center to mean
  Xm <- colMeans(refPC_new)
  Ym <- colMeans(refPC)
  Xc <- refPC_new - matrix(rep(Xm, each=NUM), ncol=length(Xm))
  Yc <- refPC - matrix(rep(Ym, each=NUM), ncol=length(Ym))
  # SVD
  C <- t(Yc) %*% Xc
  svd <- svd(C)
  # Transformation
  trXX <- sum(diag(t(Xc) %*% Xc)) # trace of a matrix: the sum of diagonal of a matrix
  trYY <- sum(diag(t(Yc) %*% Yc))
  trS <- sum(svd$d)
  A <- svd$v %*% t(svd$u)
  if(PROCRUSTES_SCALE==1){ # Orthogonal Procrustes analysis, match variance between refPC_new and refPC
    rho = sqrt(trYY/trXX)
  }else{
    rho = trS/trXX
  }
  b <- Ym - matrix(rho * Xm, nrow=1) %*% A
  # New coordinates and similarity score
  refPC_rot <- rho * refPC_new %*% A + matrix(rep(b, each=NUM), nrow=NUM)
  # Z <- refPC - refPC_rot
  # d <- sum(diag(t(Z) %*% Z))
  # D <- d / trYY
  # t <- sqrt(1-D)
  list(refPC_rot=refPC_rot, rho=rho, A=A, b=b)
}

#' Projection Procrustes Analysis
#' 
#' Adapted from http://csg.sph.umich.edu/chaolong/LASER
#' 
#' @param refPC_new Top PCs in the combination of reference samples and one study sample
#' @param refPC Top PCs in the reference samples
#' @param MAX_ITER Maximum iterations for the projection Procrustes analysis
#' @param THRESHOLD Convergence criterion for the projection Procrustes analysis
#' @param PROCRUSTES_SCALE Fit the scaling parameter to maximize similarity
#' @return Projection Procrustes Analysis results
#' @export
pprocrustes <- function(refPC_new, refPC, MAX_ITER=10000, THRESHOLD=0.000001, PROCRUSTES_SCALE=0){
  NUM = nrow(refPC_new)
  DimX = ncol(refPC_new)
  DimY = ncol(refPC)
  if(DimX < DimY){
    print("Error: dimension of refPC cannot be higher than dimension of refPC_new.")
    return(NA)
  }else if(DimX==DimY){
    procrustes(refPC_new, refPC)
  }else{
    Z <- matrix(0, nrow=NUM, ncol=DimX-DimY)
    for(i in 1:MAX_ITER){
      W <- cbind(refPC, Z)
      proc <- procrustes(refPC_new, W)
      Znew <- proc$refPC_rot[, -(1:DimY)]
      Zm <- colMeans(Znew)
      Zc <- Znew - matrix(rep(Zm, each=NUM), nrow=NUM)
      Zd <- Znew - Z
      epsilon <- sum(diag(t(Zd) %*% Zd)) / sum(diag(t(Zc) %*% Zc))
      if(epsilon<=THRESHOLD){
        break
      }else{
        Z = Znew
      }
    }
    # X2 <- refPC_rot[, -((DimY+1):DimX)]
    # proc2 <- procrustes(X2, refPC)
    proc$epsilon <- epsilon
    proc
  }
}

