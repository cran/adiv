EqRSintra <-
function(comm, dis = NULL, structures = NULL,
    option=c("eq", "normed1", "normed2"), formula = c("QE", "EDI"), tol = 1e-8, metmean = c("harmonic", "arithmetic")){

    df <- comm
    metmean <- metmean[1]
    if(!option[1]%in%c("eq", "normed1", "normed2")) stop("unavailable option, please modify; option can be eq, normed1, or normed2")
    if(!(is.data.frame(df) | is.matrix(df))) stop("df must be a matrix or a data frame")
    if(!(is.data.frame(df))) df <- as.data.frame(df)
    if(is.null(colnames(df)) | is.null(rownames(df))) stop("df must have both row names and column names")
    if(length(colnames(df)[colSums(df)>0])==1) stop("df must contain at least two species with positive occurrence in at least one site")
    if (is.null(dis)){
        dis <- as.dist((matrix(1, ncol(df), ncol(df)) - diag(rep(1, ncol(df)))))
        attributes(dis)$Labels <- colnames(df)
        formula <- "QE"
    }
    if(!inherits(dis, "dist")) stop("dis must be of class dist")
    if(!formula[1]%in%c("QE","EDI")) stop("formula can be either QE or EDI")
    if(any(!colnames(df)%in%attributes(dis)$Labels)) stop("column names in df are missing in dis")
    else{
        d <- as.matrix(dis)[colnames(df), colnames(df)]
        if(formula[1]=="EDI"){
            op <- options()$warn
            options(warn=-1)
            if(!is.euclid(as.dist(d))) stop("dis should be Euclidean")
            options(warn=op)
            d <- d^2/2 # Euclidean Diversity Index
        }
        else{
            op <- options()$warn
            options(warn=-1)
            if(!is.euclid(sqrt(as.dist(d))))  stop("dis should be squared Euclidean")
            options(warn=op)
        }
        if(any(d>1)) d <- d/max(d)
    }
    nsp <- ncol(df)
    if(any(rowSums(df)<tol)) warnings("empty sites have been removed")
    if(sum(df)<tol) stop("df must contain nonnegative values and the sum of its values must be positive")
    dfold <- df
    df <- df[rowSums(df)>tol, ]
    ncomm <- nrow(df)
    if(nrow(df)==1) stop("df must contain at least two non-empty sites")
    if(!is.null(structures)){
        if(!inherits(structures, "data.frame")) stop("structures should be a data frame or NULL")
        if(!nrow(structures)==nrow(dfold)) stop("incorrect number of rows in structures")
        structures <- structures[rowSums(dfold)>0, , drop=FALSE]
        for(i in 1:ncol(structures)) structures[,i] <-factor(structures[,i])
        if(!is.null(rownames(structures)) & !is.null(rownames(df))){
            e <- sum(abs(match(rownames(df), rownames(structures))-(1:ncomm)))
            if(e>1e-8) warning("be careful that rownames in df should be in the same order as rownames in structures")
        }
        checknested <- function(forstru){
            n <- ncol(forstru)
            for (i in 1:(n - 1)) {
                tf <- table(forstru[, c(i, i + 1)])
                niv <- apply(tf, 1, function(x) sum(x != 0))
                if (any(niv != 1)) {
                    stop(paste("non hierarchical design for structures, column", i, "is not nested in column", i + 1))
                }
            }
        }
        if(ncol(structures)> 1) checknested(structures)
    }
    P <- as.data.frame(sweep(df, 1, rowSums(df), "/"))
        if(is.null(structures)) w <- rep(1/nrow(df), nrow(df))
        else{
            nc <- ncol(structures)
            fun <- function(i){
                x <- table(structures[, i], structures[, i-1])
                x[x>0] <- 1
                x <- rowSums(x)
                v <- x[structures[, i]]
                v <- 1/v
                return(v)
            }
            if(ncol(structures)==1){
                firstw <- table(structures[, 1])
                w <- 1/firstw[structures[, 1]]/length(unique(structures[, 1]))
            }
            else {
                listw <- lapply(2:nc, fun)
                firstw <- table(structures[, 1])
                firstw <- 1/firstw[structures[, 1]]
                finalw <- 1/length(unique(structures[, ncol(structures)]))
                forw <- cbind.data.frame(as.vector(firstw), as.vector(listw), as.vector(rep(finalw, nrow(structures))))
                w <- apply(forw, 1, prod)
            }
        }
        df <- P * w

    if(!is.null(structures)){
        if(length(levels(factor(structures[, 1])))==1) stop("All sites belong to a unique level in the first column of structures, remove this first column in structures")
        if(length(levels(factor(structures[, 1])))==nrow(df)) stop("Each site belongs to a distinct level in the first column of structures, this first column is useless, remove it and re-run")
    }
    diversity <- function(x){
        if(sum(x)==0) return(0)
        else return(t(x)%*%d%*%x)
    }

    op <- options()$warn
    options(warn=-1)
    a <- apqe(as.data.frame(t(df)), sqrt(2*as.dist(d)), structures)
    options(warn=op)
    dfprop <- sweep(df, 1, rowSums(df), "/")
    divsites <- sapply(as.data.frame(t(dfprop)), diversity)
    if(metmean=="arithmetic")
        alphaEQ <- sum(w*(1/(1-divsites)))
    else
    alphaEQ <- 1/(1-sum(w*divsites))
    nlev <- nrow(a$results)
    beta1 <- (1-sum(a$results[-c(1, nlev), 1]))/(1-a$results[nlev ,1])
    if(nlev==3){
        gamma <- beta1*alphaEQ
        if(option[1]=="eq"){
            res <- cbind.data.frame(c(beta1, alphaEQ, gamma))
            rownames(res) <- c("Inter-sites", "Intra-sites", "Gamma")
            colnames(res) <- "Equivalent numbers"
        }
        else if(option[1]=="normed1"){
            beta1 <- (1-1/beta1)/(1-1/ncomm)
            res <- cbind.data.frame(beta1)
            rownames(res) <- c("Inter-sites")
            colnames(res) <- "Normed inter-site diversity"
	    }
	    else if(option[1]=="normed2"){
            beta1 <- (beta1-1)/(ncomm-1)
            res <- cbind.data.frame(beta1)
            rownames(res) <- c("Inter-sites")
            colnames(res) <- "Normed inter-site diversity"
	    }
	    return(res)
    }
    else{
        nc <- ncol(structures)
        funnc <- function(i){
            dfilist <- split(df, as.factor(structures[, i]))
            if(i > 1){
                strulist <- split(structures, as.factor(structures[, i]))
                strulist <- lapply(strulist, function(x) x[, -(i:nc), drop=FALSE])
            }
            wi <- tapply(w, as.factor(structures[, i]), sum)
            fun <- function(j){
                tab <- dfilist[[j]]
                if(i > 1)
                    thestru <- strulist[[j]]
                else
                    thestru <- NULL
                op <- options()$warn
                options(warn=-1)
                ai <- apqe(as.data.frame(t(tab)), sqrt(2*as.dist(d)), thestru)
                nlev <- nrow(ai$results)
                numi <- ai$results[nlev, 1] - ai$results[1, 1]
                gammai <- ai$results[nlev, 1]
                options(warn=op)
                ratioi <- (1-numi)/(1-gammai)
                if(option[1]=="eq"){
                    return(ratioi)
                }
                else if(option[1]=="normed1"){
                    if(!is.null(thestru)) len <- length(levels(factor(thestru[, ncol(thestru)])))
                    else len <- nrow(tab)
                    if((ratioi-1)<1e-10) res <-0
                    else
                    res <- (1 - 1/ratioi)/(1 - 1/len)
                    return(res)
                }
                else if(option[1]=="normed2"){
                    if(!is.null(thestru)) len <- length(levels(factor(thestru[, ncol(thestru)])))
                    else len <- nrow(tab)
                    if((ratioi-1)<1e-10) res <-0
                    else
                    res <- (ratioi - 1)/(len - 1)
                    return(res)
                }
            }
            resi <- unlist(lapply(1:length(dfilist), fun))
            if(metmean=="arithmetic")
                res <- sum(wi*resi)
            else
                res <- 1/sum(wi*(1/resi))
            return(res)
        }
        resbeta <- sapply(nc:1, funnc)

        if(option[1]=="eq"){
            gamma <- beta1*prod(resbeta)*alphaEQ
            res <- cbind.data.frame(c(beta1, resbeta, alphaEQ, gamma))
            inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep=""))
            intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep=""))
            intrainter <- paste(inter[-(nc+1)], intra[-1])
            rownames(res) <- c(inter[nc+1], rev(intrainter), intra[1], "Gamma")
            colnames(res) <- "Equivalent numbers"

        }
        else if(option[1]=="normed1"){
            beta1N <- (1 - 1/beta1) / (1 - 1/length(unique(structures[, nc])))
            res <- cbind.data.frame(c(beta1N, resbeta))
            inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep=""))
            intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep=""))
            intrainter <- paste(inter[-(nc+1)], intra[-1])
            rownames(res) <- c(inter[nc+1], rev(intrainter))
            colnames(res) <- "Normed contributions to beta diversity"
        }
        else if(option[1]=="normed2"){
            beta1N <- (beta1 - 1) / (length(unique(structures[, nc])) - 1)
            res <- cbind.data.frame(c(beta1N, resbeta))
            inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep=""))
            intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep=""))
            intrainter <- paste(inter[-(nc+1)], intra[-1])
            rownames(res) <- c(inter[nc+1], rev(intrainter))
            colnames(res) <- "Normed contributions to beta diversity"
        }
    	return(res)
    }

}
