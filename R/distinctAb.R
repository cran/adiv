distinctAb <- function(comm, disORtree, method = c("Q", "KY", "KstarI"), palpha = 2, option = c("asymmetric", "symmetric"), tol = 1e-10){

    option <- option[1]
    method <- method[1]
    palpha <- palpha[1]
    if (!method %in% c("Q", "KY", "KstarI")) 
        stop("unconvenient method")
    if (any(comm < 0)) 
        stop("comm should contain nonnegative values")
    commgardees <- (1:nrow(comm))[rowSums(comm) >= tol]
    commdiv <- comm[commgardees, ]
    ncommdiv <- nrow(commdiv)
    nsp <- ncol(comm)
    if(inherits(disORtree, "phylo") | inherits(disORtree, "phylo4") | inherits(disORtree, "hclust")){
        arg.phyl <- .checkphyloarg(disORtree)
        phyl.phylo <- arg.phyl$phyl.phylo
        tre4 <- arg.phyl$phyl
        if (!hasEdgeLength(tre4))
            phyl.phylo <- compute.brlen(phyl.phylo, 1)
        if(is.ultrametric(phyl.phylo) | option == c("symmetric"))
            dis <- cophenetic.phylo(phyl.phylo)/2
        else
            dis <- matrix(rep(diag(vcv.phylo(phyl.phylo))), nsp, nsp)-vcv.phylo(phyl.phylo)
        commdiv <- commdiv[, phyl.phylo$tip.label]
    }
    else if(inherits(disORtree, "matrix") | inherits(disORtree, "dist") )
         dis <- as.matrix(disORtree)
    else stop("uncovenient definition of disORtree")
    if (is.null(colnames(comm))) 
        stop("comm must have names for columns")
    if (any(!colnames(comm) %in% colnames(dis))) 
        stop("comm contains species names that are not available in disORtree")
    commdiv <- commdiv[, colnames(dis)]
    Freq <- sweep(commdiv, 1, rowSums(commdiv), "/")
    funi <- function(i){
        if(method == "Q"){
            commi <- as.vector(unlist(Freq[i,]))
            Matdis <- as.matrix(dis)
            ooi <- colSums(Matdis*commi)
        }
        if(method == "KY"){
            dis01 <- dis/max(dis)
            commi <- as.vector(unlist(Freq[i,]))
            Matdis <- as.matrix(dis01)
            ordii <- 1-colSums(Matdis*commi)
            if(abs(palpha-1) > tol){
                ooi <- (1-ordii^(palpha-1))/(palpha-1)
                ooi[abs(ordii)<tol] <- 1/(palpha-1)   
            }
            else{
                ordii[abs(ordii) < tol] <- 1
                ooi <- log(1/ordii)
            }
        }
        if(method == "KstarI"){
            dissort <- apply(dis, 1, sort)
            dissort <- t(apply(dissort, 2, diff))
            Freqord <- sapply(1:nsp, function(x) Freq[i, order(dis[x,])])
            Freqord <- apply(Freqord, 2, cumsum)
            Freqord <- Freqord[-nrow(Freqord), ]
            if(abs(palpha-1)<1e-10){
                Freqord[abs(Freqord) < tol] <- 1
                ooi <- colSums(t(dissort)*(log(1/Freqord)))
            }
            else{
                ooi <- colSums( (t(dissort)*(1 - (Freqord)^(palpha-1) ))/(palpha-1))
            }
         }
        return(ooi)
    }
    RESoo <- sapply(1:ncommdiv, funi)
    colnames(RESoo) <- rownames(comm)
    PA <- commdiv
    PA[PA>0] <- 1   
    RESoop <- t(RESoo)
    RESoop[PA<1] <- NA 
    Contr <- t(RESoo)*Freq 
    Contr[PA<1] <- 0
    if(method == "Q")
        AbRare <- 1-Freq 
    else{
        if(abs(palpha-1)<1e-10)       
            AbRare <- log(1/Freq)
        else{
           AbRare <- (1-Freq^(palpha-1))/(palpha-1)
           AbRare[Freq<tol] <- 1/(palpha-1)
        }
    }
    FunPhyRare <- t(RESoo) / AbRare
    FunPhyRarep <- FunPhyRare
    FunPhyRarep[PA<1] <- NA
    RESoop[PA<1] <- NA

    if(method == "Q" | palpha > 1 + tol) RES <- list(TotContr = Contr, EffOriPres = RESoop, EffOriAll = t(RESoo), DistinctPres = FunPhyRarep, DistinctAll = FunPhyRare, Rarity = AbRare) 
    else {
        AbRare[PA<1] <- NA
        RES <- list(TotContr = Contr, EffOriPres = RESoop, DistinctPres = FunPhyRarep, Rarity = AbRare)
    } 
    return(RES)
}
