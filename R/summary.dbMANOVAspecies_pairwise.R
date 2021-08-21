summary.dbMANOVAspecies_pairwise <- function(object, DIGITS = 3, ...){

   dbspair <- object
   if(!inherits(dbspair, "dbMANOVAspecies_pairwise")) stop("Argument object must be of class dbMANOVAspecies_pairwise")
   DIGITS <- DIGITS[1]

   funsum <- function(LList) {
        TT <- LList$test 
        if(length(TT$pvalue)>1){
        if(!is.null(TT$adj.pvalue) && !TT$adj.method%in%"none")
        TAB <- as.data.frame(t(cbind.data.frame(SES = TT$expvar[, 1], pvalue = TT$pvalue, adj.pvalue = TT$adj.pvalue)))
        else 
        TAB <- as.data.frame(t(cbind.data.frame(SES = TT$expvar[, 1], pvalue = TT$pvalue)))
        colnames(TAB) <- TT$names
        return(TAB)
        }
        else{
        VEC <- c(TT$expvar[1], TT$pvalue)
        names(VEC) <- c("SES", "pvalue") 
        return(VEC)
        }
    }
    funsumglobal <- function(LList) {
        TT <- LList$test 
        VEC <- c(SES = TT$expvar[1], pvalue = TT$pvalue)
        return(VEC)
    }

    if(names(dbspair)[1] == "global_test"){
        Npairs <- length(dbspair[[2]])
        LLspecies <- dbspair[[2]]
        RESspecies <- lapply(LLspecies, funsum)
        LLglobal <- dbspair[[1]]
        RESglobal <- lapply(LLglobal, funsumglobal)
        if(nrow(RESspecies[[1]]) == 3)
        RESglobal <- lapply(RESglobal, function(x) c(x, "NA"))
        RES <- list()
        for(i in 1:length(RESglobal)){
            res <- cbind.data.frame(GLOBAL = RESglobal[[i]], RESspecies[[i]])
            rownames(res) <- rownames(RESspecies[[i]])
            RES[[i]] <- res
        }
        names(RES) <- names(RESspecies)
    }
    else {
        Npairs <- length(dbspair)
        RES <- lapply(dbspair, funsum)
    }

    if(inherits(RES[[1]], "data.frame")){
        if(colnames(RES[[1]])[1] == "GLOBAL")
            spenames <- c("GLOBAL", attributes(dbspair)$species.names)
        else 
            spenames <- attributes(dbspair)$species.names
        FUNttab <- function(TTAB){
            ttab <- as.data.frame(matrix(NA, length(spenames), nrow(RES[[1]])))
            rownames(ttab) <- spenames
            colnames(ttab) <- rownames(RES[[1]])
            ttab[colnames(TTAB), rownames(TTAB)] <- t(TTAB)
            return(ttab)
        }
        RES <- lapply(RES, FUNttab)
        nniv <- length(RES)
        ncol <- ncol(RES[[1]])
        nam1 <- names(RES) 
        nam2 <- colnames(RES[[1]])
        nam1 <- rep(nam1, rep(ncol,nniv))
        nam2 <- rep(nam2, nniv)
        RES <- cbind.data.frame(RES)
        colnames(RES) <- paste(nam1, nam2, sep=".")
    }
    if(!inherits(RES[[1]], "data.frame")){
        RES <- cbind.data.frame(RES)
    }
    if(!is.null(DIGITS)){
        RES[!is.na(RES)] <- round(RES[!is.na(RES)], digits = DIGITS)
    }
    return(RES)
}
