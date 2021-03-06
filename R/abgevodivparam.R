abgevodivparam <-
function(phyl, comm, w = c("evoab","even","speciesab"), method = c("hillCJC", "hillR", "tsallis"), 
    q = 2, option = c("multiplicative", "additive", "proportional", "C", "U", "V", "S", "Renyi"), tol = 1e-8){
 
  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  m <- comm
  method <- method[1]
  option <- option[1]
  if(!method%in%c("tsallis", "hillR", "hillCJC")) stop("Unavailable method")
  if(!(is.numeric(q) | is.integer(q))) stop("Incorrect definition for q")
  if(is.null(colnames(m))) stop("m must have names for column")
  if(any(colSums(m)==0)){
      nsp <- ncol(m)
      nspreal <-length((1:nsp)[colSums(m) > 0])
      if(nspreal>1)
      	m <- m[, colSums(m) > 0, drop = FALSE] 
      # warning("Species with zero summed abundances over all sites have been discarded")
  }
  ncom <- nrow(m)
  if(any(m<0)) stop("m should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
  if(any(!colnames(m)%in%tipLabels(tre4))) stop("m contains tip names that are not available in phyl")
  
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

  tab1 <- mBabtot
  tab2 <- (t(t(mBabtot)*branchlengths))
 
  composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))

  if(is.numeric(w) & length(w)==ncom & all(w>0)) w <- w/sum(w)
  else if(w[1] == "evoab") w <- rowSums(tab2)/sum(tab2)
  else if(w[1] == "even") w <- rep(1/ncom, ncom)
  else if(w[1] == "speciesab") w <- rowSums(m)/sum(m)
  else stop("Incorrect definition of w")

    tsallis <- function(x, q){
        funtsallis <- function(y, q) {
            b <- branchlengths[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- -sum(b*y*log(y))
            }
            else{
            	resi <- (1-sum(b*y^q))/(q-1)
            }
            return(resi)
        }
        res <- sapply(as.data.frame(t(x)), funtsallis, q)
        return(res)
    }
    hillgamma <- function(x, q){
        funhillgamma <- function(y, q) {
            b <- branchlengths[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- exp(-sum(b*y*log(y)))
            }
            else{
            	resi <- (sum(b*y^q))^(1/(1-q))
            }
            return(resi)
        }
        xmean <- sapply(x, function(u) sum(u*w))
        res <- funhillgamma(xmean, q)
        return(res)

    }
    hillalpha <- function(x, q){

            if(abs(q-1) < tol){
                tab <- as.data.frame(x*w)
                tab2 <- tab
                tab2[tab2 < tol] <- 1
                tab[tab < tol] <- 0
                tab3 <- tab * log(tab2)
                resi <- exp(-(sum(branchlengths * sapply(tab3, sum))))/ncom
            }
            else{
                xq <- as.data.frame(x^q)
                xq[x < tol] <- 0
                resi <- (sum(branchlengths * sapply(as.data.frame(xq*w^q), sum)))^(1/(1-q))/ncom
            }
            return(resi)

    }
    hillalphaR <- function(x, q){

            if(abs(q-1) < tol){
                tab <- as.data.frame(x)
                tab2 <- tab
                tab2[tab2 < tol] <- 1
                tab[tab < tol] <- 0
                tab3 <- tab * log(tab2)
                tab3 <- sweep(tab3, 2, branchlengths, "*")
                tab3 <- as.data.frame(tab3*w)
                resi <- exp(-(sum(tab3)))
            }
            else{
                xq <- as.data.frame(x^q)
                xq[x < tol] <- 0
                tab <- sweep(xq, 2, branchlengths, "*")
                tab <- as.data.frame(tab*w)
                resi <- (sum(tab))^(1/(1-q))
            }
            return(resi)

    }
    funq <- function(q){
        if(method == "tsallis"){
            vres <- tsallis(composition, q)
            alpha <- sum(w*vres)
            compositiongamma <- as.data.frame(t(sapply(as.data.frame(composition), function(u) sum(u*w))))
            gamma <- tsallis(compositiongamma, q)
            if(option == "additive"){
                beta <- gamma - alpha
            }
            else if(option %in% c("multiplicative", "C", "U", "V", "S", "Renyi")){
                beta <- gamma / alpha
            }
            else if(option == "proportional"){
                beta <- (gamma - alpha)/gamma
            }
            else stop("The option chosen is not available")
            v <- c(alpha, beta, gamma)
            names(v) <- c("Alpha", "Beta", "Gamma")    
            return(v)
        }
        if(method == "hillCJC"){
            vresgamma <- hillgamma(composition, q) 
            vresalpha <- hillalpha(composition, q)
            if(option == "additive"){
                beta <- vresgamma - vresalpha
            }
            else if(option == "multiplicative"){
                beta <- vresgamma / vresalpha
            }
            else if(option == "proportional"){
                beta <- (vresgamma - vresalpha)/vresgamma
            }
            else if(option == "C"){
                betabrut <- vresgamma / vresalpha
                if(abs(q-1) < tol){
                beta <- log(betabrut)/log(ncom)
                }
                else{ 
                beta <- (betabrut^(1 - q) - 1) / (ncom^(1 - q) - 1)
                }
            }
            else if(option == "U"){
                betabrut <- vresgamma / vresalpha
                if(abs(q-1) < tol){
                beta <- log(betabrut)/log(ncom)
                }
                else{ 
                beta <- (1 - (1 / betabrut)^(1 - q)) / (1 - (1 / ncom)^(1 - q))
                }
            }
            else if(option == "V"){
                betabrut <- vresgamma / vresalpha
                beta <- (betabrut - 1) / (ncom - 1)
            }
            else if(option == "S"){
                betabrut <- vresgamma / vresalpha
                beta <- (1 - (1 / betabrut)) / (1 - (1 / ncom))
            }
            else if(option == "Renyi"){
                betabrut <- vresgamma / vresalpha
                beta <- log(betabrut, ncom)
            }
            else stop("The option chosen is not available")
            if(option == "C" | option == "U" | option == "S" | option == "V" | option=="Renyi"){
                v <- c(vresalpha, betabrut, beta, vresgamma)
                names(v) <- c("Alpha", "Beta", "Transformed.beta", "Gamma")
            }
            else{
                v <- c(vresalpha, beta, vresgamma)
                names(v) <- c("Alpha", "Beta", "Gamma")            
            }            
            return(v)
        }
            if(method == "hillR"){
            vresgamma <- hillgamma(composition, q) 
            vresalpha <- hillalphaR(composition, q)
            if(option == "additive"){
                beta <- vresgamma - vresalpha
            }
            else if(option %in% c("multiplicative", "C", "U", "V", "S", "Renyi")){
                beta <- vresgamma / vresalpha
            }
            else if(option == "proportional"){
                beta <- (vresgamma - vresalpha)/vresgamma
            }
            else stop("The option chosen is not available")
            v <- c(vresalpha, beta, vresgamma)
            names(v) <- c("Alpha", "Beta", "Gamma")            
            return(v)
        }
    }
    if(length(q)==1){
        v <- funq(q)
        class(v) <- "abgevodivparam"
        return(v)
    }
    if ( length(q) > 1){
           calcul1 <- sapply(q, funq)
           tab1 <- cbind.data.frame(calcul1)
           listtotale <- list()
           listtotale$q <- q
           listtotale$div <- tab1
           class(listtotale) <- "abgevodivparam"
           return(listtotale)
    }

}
