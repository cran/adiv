evodiss_family <-
function(phyl, comm, method = NULL, abundance = TRUE, squareroot = TRUE, diag = FALSE, upper = FALSE, tol = 1e-8){

  m <- comm
  METHODS <- c("Phylo - JACCARD", "Phylo - SOCKAL & MICHENER", "Phylo - SOCKAL & SNEATH",
        "Phylo - ROGERS & TANIMOTO", "Phylo - CZEKANOWSKI", "Phylo - GOWER & LEGENDRE S9",
        "Phylo - OCHIAI", "Phylo - SOKAL & SNEATH", "Phylo - Phi of PEARSON",
        "Phylo - GOWER & LEGENDRE S2")
  if (is.null(method)) {
        cat("1 = JACCARD index (1901) S3 coefficient of GOWER & LEGENDRE\n")
        cat("s1 = a/(a+b+c) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("2 = SOCKAL & MICHENER index (1958) S4 coefficient of GOWER & LEGENDRE \n")
        cat("s2 = (a+d)/(a+b+c+d) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("3 = SOCKAL & SNEATH(1963) S5 coefficient of GOWER & LEGENDRE\n")
        cat("s3 = a/(a+2(b+c)) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("4 = ROGERS & TANIMOTO (1960) S6 coefficient of GOWER & LEGENDRE\n")
        cat("s4 = (a+d)/(a+2(b+c)+d) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("5 = CZEKANOWSKI (1913) or SORENSEN (1948) S7 coefficient of GOWER & LEGENDRE\n")
        cat("s5 = 2*a/(2*a+b+c) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("6 = S9 index of GOWER & LEGENDRE (1986)\n")
        cat("s6 = (a-(b+c)+d)/(a+b+c+d) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("7 = OCHIAI (1957) S12 coefficient of GOWER & LEGENDRE\n")
        cat("s7 = a/sqrt((a+b)(a+c)) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("8 = SOKAL & SNEATH (1963) S13 coefficient of GOWER & LEGENDRE\n")
        cat("s8 = ad/sqrt((a+b)(a+c)(d+b)(d+c)) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("9 = Phi of PEARSON = S14 coefficient of GOWER & LEGENDRE\n")
        cat("s9 = (ad-bc)/sqrt((a+b)(a+c)(b+d)(d+c)) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter\n")
        cat("10 = S2 coefficient of GOWER & LEGENDRE\n")
        cat("s10 =  a/(a+b+c+d) --> d = 1 - s or d = sqrt(1 - s) depending on squareroot parameter and unit self-similarity\n")
        cat("Select an integer (1-10): ")
        method <- as.integer(readLines(n = 1))
  }
  nsp <- ncol(m)
  ncom <- nrow(m)
  if(ncom < 2) stop("At least two rows for comm are required")
  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  if(is.null(colnames(m))) stop("comm must have names for column")
  if(any(!colnames(m) %in%tipLabels(tre4))) stop("comm contains tip names that are not available in phyl")
  if(any(m<0)) stop("comm should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
  if(any(!colnames(m) %in%tipLabels(tre4))) stop("comm contains tip names that are not available in phyl")
  if(!hasEdgeLength(tre4)){
      treape <- as(tre4, "phylo")
      tre4 <- as(compute.brlen(treape, 1), "phylo4")
  }
  if(!isRooted(tre4)){
        treape <- as(tre4, "phylo")
        treape$root.edge <- 0
        tre4 <- as(treape, "phylo4")
  }

  if(!hasNodeLabels(tre4)) nodeLabels(tre4) <- names(nodeLabels(tre4))
  else{
     e <- nodeLabels(tre4)
     e[is.na(e)] <- names(e[is.na(e)])
     nodeLabels(tre4) <- e
  }

  a <- edgeLength(tre4)
  b <- a[getEdge(tre4, rootNode(tre4))]
  if(is.na(b)){
  ab <- a
  ab[getEdge(tre4, rootNode(tre4))] <- 0
  edgeLength(tre4) <- ab
  }
  tre4 <- subset(tre4, tips.exclude=tipLabels(tre4)[!tipLabels(tre4)%in%colnames(m)])

  des <- lapply(as.vector(nodeLabels(tre4)), function(x) names(descendants(tre4, x, type="tips")))
  des <- lapply(des, function(x) x[x%in%colnames(m)])
  fun <- function(namestips){
      return(rowSums(m[, namestips]))
  }
  abundancesnodes <- cbind.data.frame(lapply(des, fun))
  mBabtot <- cbind(abundancesnodes, m)
  colnames(mBabtot) <- c(nodeLabels(tre4), colnames(m))

  branchlengths <- getEdge(tre4, colnames(mBabtot), missing = "OK")

  branchlengths <- edgeLength(tre4)[branchlengths]

  if(any(is.na(branchlengths))) stop("the lengths of some branches are missing in the phylogenetic tree; note that lengths of zero are allowed")

  if(!abundance) {
        df <- as.data.frame(mBabtot)
        df[df > 0] <- 1
        df <- as.matrix(df)
  	a <- df %*%diag(branchlengths)%*% t(df)
  	b <- df %*%diag(branchlengths)%*% (1 - t(df))
  	c <- (1 - df) %*% diag(branchlengths) %*% t(df)
  	d <- sum(branchlengths) - a - b - c
   }
   else{
        mBabtot <- as.data.frame(mBabtot)
        combi1 <- rep(1:(ncom-1), (ncom-1):1)
        combi2 <- unlist(sapply(2:ncom, function(i) i:ncom))
  	    a <- sapply(1:length(combi1), function(i) sum(branchlengths*sapply(mBabtot[c(combi1[i], combi2[i]), ], min)))
        A <- matrix(0, ncom, ncom)
        A[col(A)<row(A)] <- a
        a <- A+t(A)
        a <- a + diag(sapply(1:ncom, function(i) sum(mBabtot[i, ]*branchlengths)))
  	    b <- sapply(1:length(combi1), function(i) sum(branchlengths*sapply(mBabtot[c(combi1[i], combi2[i]), ], max)) -
             sum(branchlengths*mBabtot[combi2[i], ]))
  	    c <- sapply(1:length(combi1), function(i) sum(branchlengths*sapply(mBabtot[c(combi1[i], combi2[i]), ], max)) -
             sum(branchlengths*mBabtot[combi1[i], ]))
        B <- matrix(0, ncom, ncom)
        B[col(B)<row(B)] <- b
        C <- matrix(0, ncom, ncom)
        C[col(C)<row(C)] <- c
        b <- B+t(C)
        c <- C+t(B)
        b <- b + diag(sapply(1:ncom, function(i) sum(mBabtot[i, ]*branchlengths)))
        c <- c + diag(sapply(1:ncom, function(i) sum(mBabtot[i, ]*branchlengths)))
  	    d <- sapply(1:length(combi1), function(i) sum(branchlengths*sapply(mBabtot, max))-sum(branchlengths*sapply(mBabtot[c(combi1[i], combi2[i]), ], max)))
        D <- matrix(0, ncom, ncom)
        D[col(D)<row(D)] <- d
        d <- D+t(D)
        d <- d + diag(sapply(1:ncom, function(i) sum(branchlengths*sapply(mBabtot, max))-sum(mBabtot[i, ]*branchlengths)))
   }
    method <- method[1]
    if (method == 1) {
        d <- a/(a + b + c)
    }
    else if (method == 2) {
        d <- (a + d)/(a + b + c + d)
    }
    else if (method == 3) {
        d <- a/(a + 2 * (b + c))
    }
    else if (method == 4) {
        d <- (a + d)/(a + 2 * (b + c) + d)
    }
    else if (method == 5) {
        d <- 2 * a/(2 * a + b + c)
    }
    else if (method == 6) {
        d <- (a - (b + c) + d)/(a + b + c + d)
    }
    else if (method == 7) {
        d <- a/sqrt((a + b) * (a + c))
    }
    else if (method == 8) {
        d <- a * d/sqrt((a + b) * (a + c) * (d + b) * (d + c))
        d[is.na(d)] <- 0
    }
    else if (method == 9) {
        dis <- (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *
            (d + c))
        dis[(a+b) < tol | (a+c) < tol | (b+d) < tol | (d+c) < tol] <- 0
        d <- dis
    }
    else if (method == 10) {
        d <- a/(a + b + c + d)
        diag(d) <- 1
    }
    else stop("Non convenient method")
    if(squareroot) d <- sqrt(1 - d)
    else d <- 1 - d
    d <- as.dist(d)
    attr(d, "Size") <- ncom
    attr(d, "Labels") <- rownames(m)
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)

}
