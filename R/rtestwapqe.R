rtestwapqe <-
function(comm, dis = NULL, structures = NULL,
formula = c("QE", "EDI"), wopt = c("even", "speciesab"), popt = c("aggregated", "independent"), level = 1, nrep = 99, alter = c("greater", "less", "two-sided"), tol = 1e-8){

    df <- comm
    dfold <- df
    df <- df[rowSums(df)>0, ]
    ncomm <- nrow(df)
    if(!wopt[1]%in%c("even", "speciesab")) stop("Incorrect definition of wopt")
    if(!popt[1]%in%c("aggregated", "independent")) stop("Incorrect definition of popt")
    popt <- popt[1]
    if(popt == "independent")
        df <- round(df)
    if(!is.null(structures)){
        if(!inherits(structures, "data.frame")) stop("structures should be a data frame or NULL")
        if(!nrow(structures)==nrow(dfold)) stop("incorrect number of rows in structures")
        structures <- structures[rowSums(dfold)>0, , drop=FALSE]
        structures <- as.data.frame(apply(structures, 2, factor))
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
    dfbrut <- df
    P <- as.data.frame(sweep(df, 1, rowSums(df), "/"))
    if(wopt[1] == "speciesab"){
        w <- rowSums(df)/sum(df)
    }
    else if(wopt[1] == "even"){
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
                w <- 1/firstw[structures[, 1]]/length(levels(structures[, 1]))
            }
            else {
                listw <- lapply(2:nc, fun)
                firstw <- table(structures[, 1])
                firstw <- 1/firstw[structures[, 1]]
                finalw <- 1/length(levels(structures[, ncol(structures)]))
                forw <- cbind.data.frame(as.vector(firstw), as.vector(listw), as.vector(rep(finalw, nrow(structures))))
                w <- apply(forw, 1, prod)
            }
        }
        df <- P * w
    }
    else if(is.numeric(wopt) & length(wopt) == nrow(df) & sum(wopt) > tol){
        if(!is.null(names(w)) & all(rownames(df)%in%w)) w <- w[rownames(df)]
        w <- w/sum(w)
        if(any(w<=tol)) {
            warnings("sites with weights of zero in w have been removed")
            df <- df[w>tol, ]
            structures <- structures[w>tol, ]
            w <- w[w>tol]
            w <- w/sum(w)
        }
        df <- P * w
    }
    else stop("incorrect definition of wopt")

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
    }
    d <- as.dist(d)
    alter <- alter[1]

    if(is.null(structures)){
        fun <- function(i){
            if(popt == "aggregated"){
            dfperm <- as.data.frame(sapply(df, sample))
            if(any(rowSums(dfperm)<tol)) return(NA)
            }
            else{
            dfperm <- as.data.frame(r2dtable(1,rowSums(df),colSums(df))[[1]])
            }
            if(wopt[1]!="speciesab"){
                dfperm <- as.data.frame(sweep(dfperm, 1, rowSums(dfperm), "/"))
                dfperm <- dfperm*w
            }
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(dfperm)), dis = sqrt(2*d), NULL)$results
            options(warn=op)
            res <- a[1, ]/a[2, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), NULL)$results
        options(warn=op)
        obs <- a[1, ]/a[2, ]
        sim <- ressim[!is.na(ressim)]
        res <- as.randtest(obs=obs, sim=sim, alter=alter)
        res$call <- match.call()
    }
    else if(level == 1){
        aggr.permut <- function(x){
            listval <- split(1:ncomm, as.factor(structures[, 1]))
            posiori <- as.vector(unlist(listval))
            fun0 <- function(v){
            if(length(v)==1) return(v)
                else return(sample(v))
            }
            listval2 <- lapply(listval, fun0)
            posiori2 <- as.vector(unlist(listval2))
            x[posiori] <- x[posiori2]
            return(x)
        }
        inde.permut <- function(x){
            Tperm <- as.data.frame(r2dtable(1,rowSums(x),colSums(x))[[1]])
            rownames(Tperm) <- rownames(x)
            colnames(Tperm) <- colnames(x)
            return(Tperm)
        }
        fun <- function(i){
            if(popt=="aggregated") {
            dfperm <- sapply(df, aggr.permut)
            rownames(dfperm) <- rownames(df)
            dfperm <- as.data.frame(dfperm)
            if(any(rowSums(dfperm)<tol)) return(NA)
            }
            else{
               dfsplit <- split(df, as.factor(structures[, 1]))
               dfsplitR <- lapply(dfsplit, inde.permut)
               dfperm <- as.data.frame(matrix(0, nrow(df), ncol(df)))
               rownames(dfperm) <- rownames(df)
               colnames(dfperm) <- colnames(df)
               for(i in 1:length(dfsplitR)) {
               dfperm[rownames(dfsplitR[[i]]), ] <- dfsplitR[[i]]
            }}
            if(wopt[1]!="speciesab"){
                dfperm <- as.data.frame(sweep(dfperm, 1, rowSums(dfperm), "/"))
                dfperm <- dfperm*w
            }
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(dfperm)), dis = sqrt(2*d), structures)$results
            options(warn=op)
            res <- a[nrow(a) - 2, ] / a[nrow(a) - 1, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures)$results
        options(warn=op)
        obs <- a[nrow(a) - 2, ] / a[nrow(a) - 1, ]
        sim <- ressim[!is.na(ressim)]
        res <- as.randtest(obs=obs, sim=sim, alter=alter)
        res$call <- match.call()
    }
    else if((level-1)==ncol(structures) & level==2){
        fun <- function(i){
	        e <- sample(ncomm)
	        strusim <- structures[e, , drop=FALSE]
            rownames(strusim) <- rownames(structures)
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures = strusim)$results
            options(warn=op)
            res <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures)$results
        options(warn=op)
        obs <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
        res <- as.randtest(obs=obs, sim=ressim, alter=alter)
        res$call <- match.call()
    }
    else if((level-1)==ncol(structures)){
        strulev <- structures[!duplicated(structures[level-2]), level-1]
	    names(strulev) <- unique(structures[level-2])
        fun <- function(i){
            strusim <- structures
	        strulevperm <- sample(strulev)
	        names(strulevperm) <- names(strulev)
	        strusim[, level-1] <- strulevperm[strusim[, level-2]]
            rownames(strusim) <- rownames(df)
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures = strusim)$results
            options(warn=op)
            res <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures)$results
        options(warn=op)
        obs <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
        res <- as.randtest(obs=obs, sim=ressim, alter=alter)
        res$call <- match.call()
    }
    else if(level==2){
        strulev <- as.character(structures[, level-1])
	    names(strulev) <- paste("c", 1:ncomm)
	    strulevsup <- structures[, level]
	    listbase <- split(strulev, strulevsup)
        fun <- function(i){
	        fun0 <- function(x){
	            strulevperm <- sample(x)
	            names(strulevperm) <- names(x)
                return(strulevperm)
	        }
	        listbase2 <- lapply(listbase, fun0)
	        names(listbase2) <- NULL
            listbase2 <- unlist(listbase2)
	        strusim <- structures
	        strusim[, level-1] <- as.factor(listbase2[names(strulev)])
            rownames(strusim) <- rownames(df)
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures = strusim)$results
            options(warn=op)
            res <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
            return(res)	
	    }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures)$results
        options(warn=op)
        obs <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
        res <- as.randtest(obs=obs, sim=ressim, alter=alter)
        res$call <- match.call()
    }
    else{
        if((level-1)>ncol(structures)) stop("level should be between 1 and ", ncol(structures)+1)
        strulev <- as.character(structures[!duplicated(structures[level-2]), level-1])
	    names(strulev) <- unique(structures[level-2])
	    strulevsup <- structures[!duplicated(structures[level-2]), level]
	    listbase <- split(strulev, strulevsup)
        fun <- function(i){
	        fun0 <- function(x){
	            strulevperm <- sample(x)
	            names(strulevperm) <- names(x)
                return(strulevperm)
	        }
	        listbase2 <- lapply(listbase, fun0)
	        names(listbase2) <- NULL
            listbase2 <- unlist(listbase2)
	        strusim <- structures
	        strusim[, level-1] <- as.factor(listbase2[strusim[, level-2]])
            rownames(strusim) <- rownames(df)
            op <- options()$warn
            options(warn=-1)
            a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures = strusim)$results
            options(warn=op)
            res <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
            return(res)	
	    }
        ressim <- sapply(1:nrep, fun)
        op <- options()$warn
        options(warn=-1)
        a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures)$results
        options(warn=op)
        obs <- a[nrow(a) - 1 - level, ] / a[nrow(a) - level, ]
        res <- as.randtest(obs=obs, sim=ressim, alter=alter)
        res$call <- match.call()
    }
    return(res)
}
