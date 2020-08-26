wapqe <-
function(comm, dis = NULL, structures = NULL,
    formula = c("QE", "EDI"), wopt = c("even", "speciesab"), tol = 1e-8){

    df <- comm
    dfold <- df
    df <- df[rowSums(df)>0, ]
    ncomm <- nrow(df)
    if(!is.null(structures)){
        if(!inherits(structures, "data.frame")) stop("structures should be a data frame or NULL")
        if(!nrow(structures)==nrow(dfold)) stop("incorrect number of rows in structures")
        structures <- structures[rowSums(dfold)>0, , drop=FALSE]
        for(i in 1:ncol(structures)) structures[,i] <-factor(structures[,i]) 
        if(!is.null(rownames(structures))  & !is.null(rownames(df))){
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
        if(!is.null(names(wopt)) & all(rownames(df)%in%wopt)) w <- wopt[rownames(df)]
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
    ncomm <- nrow(df)

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
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2*d), structures = structures)
    return(a$results)

}
