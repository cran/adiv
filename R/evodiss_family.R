evodiss_family <-
function(phyl, comm, method = NULL, abundance = TRUE, squareroot = TRUE, diag = FALSE, upper = FALSE, tol = 1e-8){

  m <- comm
  METHODS <- c("Phylo - Jaccard", "Phylo - Sockal & Michener", "Phylo - Gower & Legendre S5",
        "Phylo - Rogers & Tanimoto", "Phylo - Czekanowski", "Phylo - Gower & Legendre S9",
        "Phylo - Ochiai", "Phylo - Gower & Legendre S13", "Phylo - Phi of PEARSON",
        "Phylo - Gower & Legendre S2", "Phylo - Kulczynski", "Phylo - Gower & Legendre S11", "Phylo - Gower & Legendre S8")
  if (is.null(method)) {
        cat("Below are the possible methods for similarities (s) between communities\n")
        cat("dissimilarities between communities are defined as d = 1 - s or d = sqrt(1 - s) depending on the squareroot parameter\n")
        cat("1 = Jaccard index (1901) S3 coefficient of Gower and Legendre (1986) = a/(a+b+c)\n")
        cat("2 = Sokal & Michener index (1958) S4 coefficient of Gower and Legendre (1986) = (a+d)/(a+b+c+d)\n")
        cat("3 = Sokal & Sneath (1963) S5 coefficient of Gower and Legendre (1986) = a/(a+2(b+c))\n")
        cat("4 = Rogers & Tanimoto (1960) S6 coefficient of Gower and Legendre (1986) = (a+d)/(a+2(b+c)+d)\n")
        cat("5 = Czekanowski (1913) or SORENSEN (1948) S7 coefficient of Gower and Legendre (1986) = 2*a/(2*a+b+c)\n")
        cat("6 = S9 index of Gower and Legendre (1986) = (a-(b+c)+d)/(a+b+c+d)\n")
        cat("7 = Ochiai (1957) S12 coefficient of Gower and Legendre (1986) = a/sqrt((a+b)(a+c))\n")
        cat("8 = Sokal & Sneath (1963) S13 coefficient of Gower and Legendre (1986) = ad/sqrt((a+b)(a+c)(d+b)(d+c))\n")
        cat("9 = Phi of Pearson = S14 coefficient of Gower and Legendre (1986) = (ad-bc)/sqrt((a+b)(a+c)(b+d)(d+c))\n")
        cat("10 = S2 coefficient of Gower and Legendre (1986) =  a/(a+b+c+d)\n")
        cat("11 = Kulczynski index; S10 coefficient of Gower and Legendre (1986) = 0.5 * (a/(a+b) + a/(a+c)) \n")
        cat("12 = S11 coefficient of Gower and Legendre (1986) = 0.25 * (a/(a+b) + a/(a+c) + d/(b+d) + d/(c+d)) \n")
        cat("13 = S8 coefficient of Gower and Legendre (1986)  = (a+d)/(a+0.5*(b+c)+d)\n")
        cat("14 = Simpson coefficient = a/min(b,c) \n")
        cat("Select an integer (1-14): ")
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
        dis <- (a * d - b * c)/sqrt((a + b) * (a + c) * 
            (b + d) * (d + c))
        dis[(a+b) < tol | (a+c) < tol | (b+d) < tol | (d+c) < tol] <- 0
        d <- dis
    }
    else if (method == 10) {
        d <- a/(a + b + c + d)
        diag(d) <- 1
    }
    else if (method == 11) {
        dis <- 0.5 * (a/(a+b) + a/(a+c))
        dis[a < tol] <- 0
        d <- dis
    }
    else if (method == 12) {
        dis <- 0.25 * (a/(a+b) + a/(a+c) + d/(b+d) + d/(c+d))
        dis[a < tol & d < tol] <- 0
        disA <- 0.25 * (a/(a+b) + a/(a+c))
        dis[a < tol & d >= tol] <- disA[a < tol & d >= tol]
        disD <- 0.25 * (d/(b+d) + d/(c+d))
        dis[a >= tol & d < tol] <- disD[a >= tol & d < tol]
        d <- dis
    }
    else if (method == 13) {
        d <- (a+d)/(a+0.5*(b+c)+d)
    }
    else if (method == 14) {
        minbc <- (b<=c)*b+(b>c)*c
        d <- a/(a+minbc)
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
