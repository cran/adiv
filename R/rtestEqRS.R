rtestEqRS <-
function(comm, dis = NULL, structures = NULL,
formula = c("QE", "EDI"), option=c("normed1", "normed2", "eq"), popt = c("aggregated", "independent"), level = 1, nrep = 99,
alter = c("greater", "less", "two-sided"), tol = 1e-8){

    df <- comm
    if(!option[1]%in%c("normed1", "normed2", "eq")) stop("Incorrect definition of option")
    if(!popt[1]%in%c("aggregated", "independent")) stop("Incorrect definition of popt")
    popt <- popt[1]
    if(popt == "independent")
        df <- round(df)
    if(option[1]=="eq")
        indexk <- 0
    else
        indexk <- 2
    dfold <- df
    df <- df[rowSums(df)>0, ]
    ncomm <- nrow(df)
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
            res <- EqRS(dfperm, dis = dis, NULL, formula = formula, option = option, tol = tol)[1, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        obs <- EqRS(df, dis = dis, NULL, formula = formula, option = option, tol = tol)[1, ]
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
            res <- EqRS(dfperm, dis = dis, structures, formula = formula, option = option, tol = tol)
            res <- res[nrow(res)- 2 + indexk, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
        obs <- obs[nrow(obs)-2 + indexk, ]
        sim <- ressim[!is.na(ressim)]
        res <- as.randtest(obs=obs, sim=sim, alter=alter)
        res$call <- match.call()
    }
    else if((level-1)==ncol(structures) & level==2){
        fun <- function(i){
	        e <- sample(ncomm)
	        strusim <- structures[e, , drop=FALSE]
            rownames(strusim) <- rownames(df)
            res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
            res <- res[1, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
	    obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
        obs <- obs[1, ]
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
            res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
            res <- res[1, ]
            return(res)
        }
        ressim <- sapply(1:nrep, fun)
        obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
        obs <- obs[1, ]
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
            res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
            res <- res[nrow(res) - 3 + indexk, ]
            return(res)	
	    }
        ressim <- sapply(1:nrep, fun)
        obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
        obs <- obs[nrow(obs) - 3 + indexk, ]
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
            res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
            res <- res[nrow(res) - 2 + indexk - level + 1, ]
            return(res)	
	    }
        ressim <- sapply(1:nrep, fun)
        obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
        obs <- obs[nrow(obs) - 2 + indexk - level + 1, ]
        res <- as.randtest(obs=obs, sim=ressim, alter=alter)
        res$call <- match.call()
    }
    return(res)
}
