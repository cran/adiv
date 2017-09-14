evodiss_ternaryplot <-
function(phyl, comm, abundance = TRUE, tol = 1e-8, ...){
  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  m <- comm
  nsp <- ncol(m)
  ncom <- nrow(m)
  if(ncom < 2) stop("At least two rows for m are required")
  if(is.null(colnames(m))) stop("m must have names for column")
  if(any(!colnames(m) %in%tipLabels(tre4))) stop("m contains tip names that are not available in phyl")
  if(any(m<0)) stop("m should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
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
   if(ncom==2){
       a <- as.vector(a)[2]
       b <- as.vector(b)[2]
       c <- as.vector(c)[2]
   }
   else{
       a <- as.vector(a)[rep(c(0, (1:(ncom-2))*ncom), (ncom-1):1)+unlist(sapply(2:ncom, function(x) x:ncom))]
       b <- as.vector(b)[rep(c(0, (1:(ncom-2))*ncom), (ncom-1):1)+unlist(sapply(2:ncom, function(x) x:ncom))]
       c <- as.vector(c)[rep(c(0, (1:(ncom-2))*ncom), (ncom-1):1)+unlist(sapply(2:ncom, function(x) x:ncom))]
   }
   tabc <- cbind.data.frame(a,b,c)
   if(abundance){
       colnames(tabc) <- LETTERS[c(1,3,2)]
       tabc <- tabc[, LETTERS[1:3]]
   }
   else
       colnames(tabc) <- letters[1:3]
   if(is.null(rownames(m))) rownames(m) <- paste("com", 1:ncom, sep="")
   part2 <- rep(rownames(m), (ncom-1):0)
   part1 <- unlist(sapply(2:ncom, function(x) rownames(m)[x:ncom]))
   rownames(tabc) <- paste(part1, part2, sep=":")
   triangle.label(tabc, labels=rownames(tabc), ...)

}
