rlqESLTP <-
function(dudiE, dudiS, dudiL, dudiT, dudiP, ...){
    if(!is.null(dudiE))
	tabE <- dudiE$li/sqrt(dudiE$eig[1])
    if(!is.null(dudiS))
	tabS <- dudiS$li/sqrt(dudiS$eig[1])
    if(!is.null(dudiP))
	tabP <- dudiP$li/sqrt(dudiP$eig[1])
    if(!is.null(dudiT))
	tabT <- dudiT$li/sqrt(dudiT$eig[1])
    if(!is.null(dudiE)&!is.null(dudiS)){
	tabES <- cbind.data.frame(tabE, tabS)
	names(tabES) <- c(paste("E", 1:ncol(tabE), sep = ""),
		paste("S", 1:ncol(tabS), sep = ""))
	pcaES <- dudi.pca(tabES, scale = FALSE, row.w = dudiL$lw, scannf = FALSE,
		nf = (length(dudiE$eig) + length(dudiS$eig)))
    }
    else{
        if(!is.null(dudiE))
            pcaES <- dudiE
        else{
            if(!is.null(dudiS))
                pcaES <- dudiS
            else 
                stop("one table giving site attributes is missing")
        }
    }
    if(!is.null(dudiT)&!is.null(dudiP)){    
	tabTP <- cbind.data.frame(tabT, tabP)
	names(tabTP) <- c(paste("T", 1:ncol(tabT), sep = ""),
		paste("P", 1:ncol(tabP), sep = ""))
	pcaTP <- dudi.pca(tabTP, scale = FALSE, row.w = dudiL$cw, scannf = FALSE,
		nf = (length(dudiT$eig) + length(dudiP$eig)))
    }
    else{
        if(!is.null(dudiT))
            pcaTP <- dudiT
        else{
            if(!is.null(dudiP))
                pcaTP <- dudiP
            else
                stop("one table giving species attributes is missing")
        }
    }

	X <- rlq(pcaES, dudiL, pcaTP, ...)

    if(!is.null(dudiE)){
	U <- as.matrix(X$l1) * unlist(X$lw)
	U <- data.frame(as.matrix(pcaES$tab[, 1:ncol(tabE)]) %*% U[1:ncol(tabE), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
	names(U) <- names(X$lR)
	X$lR_givenE <- U
	}
    if(!is.null(dudiS)){
	U <- as.matrix(X$l1) * unlist(X$lw)
    if(!is.null(dudiE)){
	U <- data.frame(as.matrix(pcaES$tab[, -(1:ncol(tabE))]) %*% U[-(1:ncol(tabE)), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
    }
    else{
	U <- data.frame(as.matrix(pcaES$tab) %*% U[, 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
    }    
	names(U) <- names(X$lR)
	X$lR_givenS <- U
    }
    if(!is.null(dudiT)){
	U <- as.matrix(X$c1) * unlist(X$cw)
	U <- data.frame(as.matrix(pcaTP$tab[, 1:ncol(tabT)]) %*% U[1:ncol(tabT), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
	names(U) <- names(X$lQ)
	X$lQ_givenT <- U
    }
    if(!is.null(dudiP)){
	U <- as.matrix(X$c1) * unlist(X$cw)
    if(!is.null(dudiT)){
	U <- data.frame(as.matrix(pcaTP$tab[, -(1:ncol(tabT))]) %*% U[-(1:ncol(tabT)), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
    }
    else{
	U <- data.frame(as.matrix(pcaTP$tab) %*% U[, 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)    
    }
	names(U) <- names(X$lQ)
	X$lQ_givenP <- U
    }
	X$row.w <- dudiL$lw

	X$col.w <- dudiL$cw

    X$dudiL <- dudiL
    
    X$dudiR <- pcaES
    X$dudiQ <- pcaTP
	
	class(X) <- c("rlqESLTP", "rlq", "dudi")

	return(X)
}
