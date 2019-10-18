evodiss <-
function(phyl, comm, method = NULL, q = NULL, w = c("evoab", "even", "speciesab"), diag = FALSE, upper = FALSE, tol = 1e-8){

  m <- comm
  nsp <- ncol(m)
  ncom <- nrow(m)
  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl

  if (is.null(method)) {
        cat("enter the name of the method to be used\n")
        cat("for calculating evodissimilarity (see help file):")
        method <- readLines(n = 1)
  }

  method <- method[1]
  if(!method%in%c("Minkowski","Euclidean","Manhattan","Chord","ScaledCanberra","Divergence","BC","MH","LG","Hellinger","chi2","Hill","Renyi", "C", "U","S"))
    stop("Unavailable method")

  r <- q
  if(is.null(colnames(m))) stop("comm must have names for column")
  if(ncom < 2) stop("At least two rows for comm are required")
  if(is.null(colnames(m))) stop("comm must have names for column")
  if(any(!colnames(m) %in%tipLabels(tre4))) stop("comm contains tip names that are not available in phyl")
  if(any(m<0)) stop("comm should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded in comm")
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

  if(method%in%c("Hill", "Renyi", "C", "U", "S")){
  tab2 <- (t(t(mBabtot)*branchlengths))
  if(is.numeric(w) & length(w)==ncom & all(w>0)) w <- w
  else if(w[1] == "evoab") w <- rowSums(tab2)
  else if(w[1] == "even") w <- rep(1/ncom, ncom)
  else if(w[1] == "speciesab") w <- rowSums(m)
  else stop("Incorrect definition of w")
  names(w) <- rownames(m)
  }
  hillgamma <- function(x, branch, q, wcom){
        funhillgamma <- function(y, q) {
            b <- branch[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- exp(-sum(b*y*log(y)))
            }
            else{
            	resi <- (sum(b*y^q))^(1/(1-q))
            }
            return(resi)
        }
        xmean <- sapply(x, function(u) sum(u*wcom))
        res <- funhillgamma(xmean, q)
        return(res)

    }
    hillalpha <- function(x, branch, q, wcom){

            if(abs(q-1) < tol){
                tab <- as.data.frame(x*wcom)
                tab2 <- tab
                tab2[tab2 < tol] <- 1
                tab[tab < tol] <- 0
                tab3 <- tab * log(tab2)
                resi <- exp(-(sum(branch * sapply(tab3, sum))))/2
            }
            else{
                xq <- as.data.frame(x^q)
                xq[x < tol] <- 0
                resi <- (sum(branch * sapply(as.data.frame(xq*wcom^q), sum)))^(1/(1-q))/2
            }
            return(resi)

    }

  method <- method[1]

  if (method == "Euclidean") {
        tab <- t(t(mBabtot)*sqrt(branchlengths))
        d <- dist(tab)
  }
  else if (method == "Minkowski") {
    if(r < tol) stop("q must be positive with Minkowski index")
    tab <- t(t(mBabtot)*(branchlengths)^(1/r))
    d <- matrix(0, ncom, ncom)
    fun1 <- function(x) {
        (sum((abs(tab[x[1], ] - tab[x[2], ]))^r))^(1/r)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun1))
  }
  else if (method == "Manhattan") {
    tab <- t(t(mBabtot)*(branchlengths))
    d <- matrix(0, ncom, ncom)
    fun1 <- function(x) {
        sum((abs(tab[x[1], ] - tab[x[2], ])))
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun1))
  }
  else if (method == "Chord") {
    tab <- t(t(mBabtot)*(branchlengths)^(1/2))
    d <- matrix(0, ncom, ncom)
    fun2 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        w0 <- 2*(1 - sum(p * q)/sqrt(sum(p * p))/sqrt(sum(q * q)))
        if(abs(w0)<tol) return(0)
        w <- sqrt(w0)
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun2))
  }
  else if (method == "ScaledCanberra") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun3 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        b <- branchlengths[(p + q) > 0]
        ps <- p[(p + q) > 0]
        qs <- q[(p + q) > 0]
        w <- sum(b * abs(ps - qs)/(ps + qs))/sum(branchlengths)
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun3))
  }

   else if (method == "Divergence") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun4 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        b <- branchlengths[(p + q) > 0]
        ps <- p[(p + q) > 0]
        qs <- q[(p + q) > 0]
        w <- sqrt(sum(b * ((ps - qs)/(ps + qs))^2)/sum(branchlengths))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun4))
  }
  else if (method == "BC") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun5 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        b <- branchlengths[(p + q) > 0]
        ps <- p[(p + q) > 0]
        qs <- q[(p + q) > 0]
        w <- sum(abs(b * ps - b * qs))/sum(b * ps + b * qs)
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun5))
  }
    else if (method == "MH") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun5 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        b <- branchlengths[(p + q) > 0]
        ps <- p[(p + q) > 0]
        qs <- q[(p + q) > 0]
        w <- sum(b*abs(ps - qs)^2)/sum(b * ps^2 + b * qs^2)
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun5))
  }
  else if (method == "LG") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun6 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        w <- sqrt(sum(branchlengths *( p/sum(branchlengths *p) - q/sum(branchlengths * q))^2))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun6))
  }
 else if (method == "Hellinger") {
    tab <- mBabtot
    d <- matrix(0, ncom, ncom)
    fun6 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        w <- sqrt(sum(branchlengths *(sqrt( p/sum(branchlengths *p)) - sqrt(q/sum(branchlengths * q)))^2))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun6))
  }
  else if (method == "chi2") {
    tab <- mBabtot[, colSums(mBabtot) > 0]
    d <- matrix(0, ncom, ncom)
    fun7 <- function(x) {
        p <- tab[x[1], ]
        q <- tab[x[2], ]
        b <- branchlengths[colSums(mBabtot) > 0]
        w <- sqrt(sum(b * (sum(b * colSums(tab)) / colSums(tab)) *
              (p/sum(b * p) - q/sum(b * q))^2))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun7))
  }
  else if (method == "Hill"){
    if(r < 0) stop("q should be nonnegative")
    tab1 <- mBabtot[, colSums(mBabtot) > 0]
    branchlengths <- branchlengths[colSums(mBabtot) > 0]
    tab2 <- (t(t(tab1)*branchlengths))
    composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))
    d <- matrix(0, ncom, ncom)
    fun9 <- function(x) {
        vresgamma <- hillgamma(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        vresalpha <- hillalpha(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        beta <- vresgamma / vresalpha
        w <- beta - 1
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun9))
  }
  else if (method == "Renyi"){
    if(r < 0) stop("q should be nonnegative")
    tab1 <- mBabtot[, colSums(mBabtot) > 0]
    branchlengths <- branchlengths[colSums(mBabtot) > 0]
    tab2 <- (t(t(tab1)*branchlengths))
    composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))
    d <- matrix(0, ncom, ncom)
    fun9 <- function(x) {
        vresgamma <- hillgamma(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        vresalpha <- hillalpha(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        beta <- vresgamma / vresalpha
        w <- log(beta, 2)
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun9))
  }
  else if (method == "C"){
    if(r < 0) stop("q should be nonnegative")
    tab1 <- mBabtot[, colSums(mBabtot) > 0]
    branchlengths <- branchlengths[colSums(mBabtot) > 0]
    tab2 <- (t(t(tab1)*branchlengths))
    composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))
    d <- matrix(0, ncom, ncom)
    fun9 <- function(x) {
        vresgamma <- hillgamma(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        vresalpha <- hillalpha(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        beta <- vresgamma / vresalpha
        if(abs(r-1)<tol){
            w <- log(beta, 2)
        }
        else
            w <- 1-((1/beta)^(r-1)-(1/2)^(r-1))/(1-(1/2)^(r-1))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun9))
  }
  else if (method == "U"){
    if(r < 0) stop("q should be nonnegative")
    tab1 <- mBabtot[, colSums(mBabtot) > 0]
    branchlengths <- branchlengths[colSums(mBabtot) > 0]
    tab2 <- (t(t(tab1)*branchlengths))
    composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))
    d <- matrix(0, ncom, ncom)
    fun9 <- function(x) {
        vresgamma <- hillgamma(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        vresalpha <- hillalpha(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        beta <- vresgamma / vresalpha
        if(abs(r-1)<tol){
            w <- log(beta, 2)
        }
        else
            w <- 1-((1/beta)^(1-r)-(1/2)^(1-r))/(1-(1/2)^(1-r))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun9))
  }
  else if (method == "S"){
    if(r < 0) stop("q should be nonnegative")
    tab1 <- mBabtot[, colSums(mBabtot) > 0]
    branchlengths <- branchlengths[colSums(mBabtot) > 0]
    tab2 <- (t(t(tab1)*branchlengths))
    composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))
    d <- matrix(0, ncom, ncom)
    fun9 <- function(x) {
        vresgamma <- hillgamma(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        vresalpha <- hillalpha(composition[c(x[1], x[2]), ], branchlengths, r, w[c(x[1], x[2])]/sum(w[c(x[1], x[2])]))
        beta <- vresgamma / vresalpha
        w <- 1-2*((1/beta)-(1/2))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, fun9))
  }
 else stop("Non convenient method")
    attr(d, "Size") <- ncom
    attr(d, "Labels") <- rownames(m)
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- method
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
