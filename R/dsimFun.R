dsimFun <-
function(df, vartype=c("Q","N","M","P"), method=1:5, type=c("similarity", "dissimilarity")){

    type <- type[1]
    if(!type%in% c("dissimilarity", "similarity")) stop("type must be either dissimilarity or similarity")
    meantype <- method[1]
    if(!meantype%in%(1:5)) stop("Incorrect definition of method")
    fun0 <- function(i){
    df0 <- as.matrix(df[[i]])
    type <- type[1]
    vartype0 <- vartype[i]
    if(vartype0=="Q" | vartype0=="N"){
        if(type=="dissimilarity")
            return(daisy(df0, metric = "gower")*ncol(df0))
        else
            return((1-as.matrix(daisy(df0, metric = "gower")))*ncol(df0))
    }
    if(vartype0=="P"){
        df0 <- sweep(df0, 1, rowSums(df0), "/")
    }
    if(vartype0=="P" | vartype0=="M"){
        A <- df0%*%t(df0)
        B <- diag(A)%*%t(rep(1, nrow(df0)))
        C <- rep(1, nrow(df0))%*%t(diag(A))
        if(meantype==4) S <- A/sqrt(B)/sqrt(C)
        else if(meantype==3){
            S <- 2*A/(B+C)
        }
        else if(meantype==1){
            S <- A/(2*B+2*C-3*A)
        }
        else if(meantype==2){
            S <- A/(B+C-A)
        }
        else S <- 4*A/(2*A+B+C)

        rownames(S)<-colnames(S)<-rownames(df0)
        if(type=="dissimilarity")
            return(as.dist(1-S))
        else
            return(S)
    }
    }
    if(inherits(df, "ktab")){
        listdsim <- lapply(1:length(df$blo), fun0)
        res <- listdsim[[1]]
        if(length(listdsim)>1){
            for(i in 2:length(listdsim))
                res <- res + listdsim[[i]]
        }
        nk <- length(vartype[vartype!="Q" & vartype!="N"])
        nk <- nk + sum(df$blo[vartype=="Q" | vartype=="N"])
        return(res/nk)
    }
    else{
    df <- as.matrix(df)
    type <- type[1]
    vartype <- vartype[1]
    if(vartype=="Q" | vartype=="N"){
        if(type=="dissimilarity")
            return(daisy(df, metric = "gower"))
        else
            return(1-as.matrix(daisy(df, metric = "gower")))
    }
    if(vartype=="P"){
        df <- sweep(df, 1, rowSums(df), "/")
    }
    if(vartype=="P" | vartype=="M"){
        A <- df%*%t(df)
        B <- diag(A)%*%t(rep(1, nrow(df)))
        C <- rep(1, nrow(df))%*%t(diag(A))
        if(meantype==4) S <- A/sqrt(B)/sqrt(C)
        else if(meantype==3){
            S <- 2*A/(B+C)
        }
        else if(meantype==1){
            S <- A/(2*B+2*C-3*A)
        }
        else if(meantype==2){
            S <- A/(B+C-A)
        }
        else S <- 4*A/(2*A+B+C)

        rownames(S)<-colnames(S)<-rownames(df)
        if(type=="dissimilarity")
            return(sqrt(1-S))
        else
            return(S)
    }
    }
}
