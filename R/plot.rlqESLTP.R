plot.rlqESLTP <-
function(x, which = NULL, phyl = NULL, xy = NULL, traits = NULL,
	env = NULL, type = NULL, ax = 1, ...){

    if(is.null(which)) stop("Specify which graph you would like to plot")
	
	if(which == "S"){
	if(is.null(xy)) stop("xy required")
    if(!is.null(x$lR_givenE)){
	par(mfrow = c(1, 3))
	s.value(xy, x$lR_givenE[, ax], zmax = max(x$lR[, ax]), addaxes = F, clegend = 2,
		sub = "environment-based", csub = 2, include.origin = FALSE)
	s.value(xy, x$lR_givenS[, ax], zmax = max(x$lR[, ax]), addaxes = F, clegend = 2,
		sub = "space-based", csub = 2, include.origin = FALSE)
	s.value(xy, x$lR[,ax], zmax = max(x$lR[, ax]), addaxes = F, clegend = 2,
		sub = "global", csub = 2, include.origin = FALSE)
    }
    else
	s.value(xy, x$lR[,ax], zmax = max(x$lR[, ax]), addaxes = F, clegend = 2,
		sub = "global", csub = 2, include.origin = FALSE)        
	}

	if(which == "P"){
	if(is.null(phyl)) stop("phyl required")
     arg.phyl <- .checkphyloarg(phyl)
     phyl <- arg.phyl$phyl.phylo

    if(!is.null(x$lQ_givenT)){
    parmar <- par()$mar
    par(mar=rep(.1,4))
     CB <- cbind.data.frame(x$lQ_givenT[phyl$tip.label, ax], x$lQ_givenP[phyl$tip.label, ax], x$lQ[phyl$tip.label, ax])
     colnames(CB) <- c("trait-based", "phylogeny-based","global")  
     X.4d <- phylo4d(phyl, as.matrix(CB))
table.phylo4d(X.4d, show.node.label=FALSE, symbol="squares", center=FALSE, scale=FALSE)
    par(mar=parmar)
    }
    else{
    parmar <- par()$mar
    par(mar=rep(.1,4))
	CB <- as.data.frame(x$lQ[phyl$tip.label, ax])
      colnames(CB) <- "global"
     X.4d <- phylo4d(phyl, as.matrix(CB))
table.phylo4d(X.4d, show.node.label=FALSE, symbol="squares", center=FALSE, scale=FALSE)
    par(mar=parmar)
     }   
	}

	if(which == "T" | which == "E"){
		
        mfrow = n2mfrow(length(type))
        par(mfrow = mfrow)

		if(which == "T"){
			ltab <- traits
			w <- x$col.w
            sco1 <- x$lQ
            sco2 <- x$mQ
		}
		else{
			ltab <- env
			w <- x$row.w
            sco1 <- x$lR
            sco2 <- x$mR
		}
		if (is.data.frame(ltab)) ltab <- list(ltab)
		for(i in 1:length(ltab)){
			if(type[i] == "Q"){
				thetab <- ltab[[i]]
                if(!any(is.na(thetab))){
                    thetabS <- scalewt(thetab, w)
                    corS <- (t(thetabS)%*%diag(w)%*%sco2[, ax])[, 1]
                }
                else{
                    funcorS <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        sco2sna <- sco2[!is.na(x), ax]
                        wsna <- w[!is.na(x)]
                        thetabSsna <- scalewt(xsna, wsna)
                        corSsna <- t(thetabSsna)%*%diag(wsna)%*%sco2sna
                        return(corSsna)
                    }
                    corS <- sapply(1:ncol(thetab), funcorS)
                    names(corS) <- names(thetab)
                }
				dotchart(sort(corS), labels = rownames(corS)[order(corS)],
					main = "Pearson correlation")
				abline(v = 0)

			}
			if(type[i] == "O"){
				thetab <- ltab[[i]]
				thetab <- as.data.frame(apply(thetab, 2, rank))

                if(!any(is.na(thetab))){
                    thetabS <- scalewt(thetab, w)
                    corS <- t(thetabS)%*%diag(w)%*%scalewt(rank(sco2[, ax]), w)
                }
                else{
                    funcorS <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        wsna <- w[!is.na(x)]
                        sco2sna <- scalewt(rank(sco2[!is.na(x), ax]), wsna)
                        thetabSsna <- scalewt(xsna, wsna)
                        corSsna <- t(thetabSsna)%*%diag(wsna)%*%sco2sna
                        return(corSsna)
                    }
                    corS <- sapply(1:ncol(thetab), funcorS)
                    names(corS) <- names(thetab)
                }
				dotchart(sort(corS), labels = rownames(corS)[order(corS)],
					main = "Spearman correlation")
				abline(v = 0)
				
			}
			if(type[i] == "N"){
                thetab <- ltab[[i]]
                funmod <- function(unx){

                    if(!any(is.na(unx))){
                    mod <- model.matrix(~-1+factor(unx))
                    colnames(mod) <- levels(factor(unx))
                    rownames(mod) <- rownames(thetab)
                    return(as.data.frame(mod))
                    }
                    else{
                        mod <- model.matrix(~-1+factor(unx))
                        correctedtab <- matrix(NA, nrow(thetab), ncol(mod))
                        correctedtab[as.numeric(rownames(mod)), ] <- mod
                        colnames(correctedtab) <- levels(factor(unx))
                        rownames(correctedtab) <- rownames(thetab)
                        return(as.data.frame(correctedtab))
                    }
                    }
                    res <- cbind.data.frame(apply(thetab, 2, funmod))
                    res[is.na(res)] <- 0
                    sco.distri(sco1[, ax], res)
			}
			if(type[i] == "F" | type[i] == "B" | type[i] == "D"){
                thetab <- ltab[[i]]
                thetab[is.na(thetab)] <- 0
				sco.distri(sco1[, ax], thetab)
			}
            if(type[i] == "C"){
                thetab <- ltab[[i]]
                if(!any(is.na(thetab))){
                    alphat <- t(t(thetab * 2 * pi)/attributes(thetab)$max)
                    alphatcos <- scalewt(cos(alphat), w)
                    alphatsin <- scalewt(sin(alphat), w)
                    rxc <- t(alphatcos)%*%diag(w)%*%sco2[, ax]
                    rxs <- t(alphatsin)%*%diag(w)%*%sco2[, ax]
                    rcs <- diag(t(alphatsin)%*%diag(w)%*%alphatcos)
                    corC <- (sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 -
                    rcs^2)))[, 1]
                }
                else{
                    funcorC <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        sco2sna <- sco2[!is.na(x), ax]
                        wsna <- w[!is.na(x)]
                        alphat <- xsna * 2 * pi/attributes(thetab)$max[j]
                        alphatcos <- scalewt(cos(alphat), wsna)
                        alphatsin <- scalewt(sin(alphat), wsna)
                        rxc <- t(alphatcos)%*%diag(wsna)%*%sco2sna
                        rxs <- t(alphatsin)%*%diag(wsna)%*%sco2sna
                        rcs <- diag(t(alphatsin)%*%diag(wsna)%*%alphatcos)
                        corCsna <- sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 - rcs^2))
                        return(corCsna)
                    }
                    corC <- sapply(1:ncol(thetab), funcorC)
                    names(corC) <- names(thetab)
                }
           		dotchart(sort(corC), labels = rownames(corC)[order(corC)],
					main = "Circular correlation")
				abline(v = 0)

            }
		}
	}

    par(mfrow=c(1, 1))

}
