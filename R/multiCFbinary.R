multiCFbinary <-
function (Ktab, w.attributes=lapply(Ktab, function(x) rep(1,ncol(x))), w.traits=rep(1/length(Ktab),length(Ktab)), labels=rownames(Ktab[[1]]), solution=c(2,1))
{
    if(!inherits(Ktab, "ktab"))
        K <- ktab.list.df(Ktab)
    else K <- Ktab
    nK <- length(K$blo)
    for(i in 1:nK){
    df <- K[[i]]
    if (any(df < 0))
        stop("non negative value expected in df")
    dfs <- apply(df, 1, sum)
    if (any(dfs == 0))
        stop("row with all zero value")
    }
    wA <- w.attributes
    wT <- w.traits
    solution <- solution[1]
    if(solution == 1){
        listC <- lapply(1:nK, function(k) CFbinary(K[[k]], wA[[k]]))
        matC <- wT[1]*listC[[1]]
        if(nK>1){
        for (i in 2:nK){
            matC <- matC + wT[i] * listC[[i]]
        }
        }
        colnames(matC)<-rownames(matC)<-labels
        return(matC)
    }
    if(solution == 2){
        listC <- lapply(1:nK, function(k) CFbinary(K[[k]], wA[[k]]))
        listC <- lapply(listC, function(x) diag(1/diag(x))%*%x%*%diag(1/diag(x)))
        matC <- wT[1]*listC[[1]]
        if(nK>1){
        for (i in 2:nK){
            matC <- matC + wT[i] * listC[[i]]
        }
        }
            W <- diag(1/diag(matC))
            matC <- W%*%matC%*%W
            colnames(matC)<-rownames(matC)<-labels
            return(matC)
    }
}
