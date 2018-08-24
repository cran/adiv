dsimcom <-
function(comm, Sigma = NULL, method = 1:5, option=c("relative", "absolute"), type=c("similarity", "dissimilarity")){

    type <- type[1]
    if(!type%in% c("dissimilarity", "similarity")) stop("type must be either dissimilarity or similarity")
# In df: communities as rows.
    df <- comm
    S <- Sigma
    if(!inherits(df, "data.frame"))
        df <- as.data.frame(df)
    df <- as.data.frame(t(df))
    if(is.null(S))
        S <- diag(rep(1, nrow(df)))
    S[S<0] <- 0
    S[S>1] <- 1
    method <- method[1]
    option <- option[1]
    if(option=="relative")
        dfp <- t(t(df)/colSums(df))
    else
        dfp <- as.matrix(df)
    fun2 <- function(dfp, S){
        A <- t(dfp)%*%S%*%dfp
        B <- matrix(1, ncol(dfp), ncol(dfp))%*%diag(diag(A))
        C <- diag(diag(A))%*%matrix(1, ncol(dfp), ncol(dfp))
        s <- A / (B+C-A)
        rownames(s)<-colnames(s)<-colnames(df)
        s[s<0]<-0
        s[s>1]<-1
        return(s)
    }
    fun1 <- function(dfp, S){
        A <- t(dfp)%*%S%*%dfp
        B <- matrix(1, ncol(dfp), ncol(dfp))%*%diag(diag(A))
        C <- diag(diag(A))%*%matrix(1, ncol(dfp), ncol(dfp))
        s <- A / (2*B+2*C-3*A)
        rownames(s)<-colnames(s)<-colnames(df)
        s[s<0]<-0
        s[s>1]<-1
        return(s)
    }
    fun3 <- function(dfp, S){
        A <- t(dfp)%*%S%*%dfp
        B <- matrix(1, ncol(dfp), ncol(dfp))%*%diag(diag(A))
        C <- diag(diag(A))%*%matrix(1, ncol(dfp), ncol(dfp))
        s <- 2*A / (B+C)
        rownames(s)<-colnames(s)<-colnames(df)
        s[s<0]<-0
        s[s>1]<-1
        return(s)
    }
    fun4 <- function(dfp, S){
        C <- t(dfp)%*%S%*%dfp
        W <- diag(1/sqrt(diag(C)))
        Scom <- W%*%C%*%W
        rownames(Scom)<-colnames(Scom)<-colnames(df)
        Scom[Scom<0]<-0
        Scom[Scom>1]<-1
        return(Scom)
    }
    fun5 <- function(dfp, S){
        A <- t(dfp)%*%S%*%dfp
        B <- matrix(1, ncol(dfp), ncol(dfp))%*%diag(diag(A))
        C <- diag(diag(A))%*%matrix(1, ncol(dfp), ncol(dfp))
        s <- 4*A / (2*A+B+C)
        rownames(s)<-colnames(s)<-colnames(df)
        s[s<0]<-0
        s[s>1]<-1
        return(s)
    }
    if(type=="similarity"){
    if(method == 1)
        return(fun1(dfp, S))
    if(method == 2)
        return(fun2(dfp, S))
    if(method == 3)
        return(fun3(dfp, S))
    if(method == 4)
        return(fun4(dfp, S))
    if(method == 5)
        return(fun5(dfp, S))
    }
    else{
    if(method == 1)
        return(as.dist(1-fun1(dfp, S)))
    if(method == 2)
        return(as.dist(1-fun2(dfp, S)))
    if(method == 3)
        return(as.dist(1-fun3(dfp, S)))
    if(method == 4)
        return(as.dist(1-fun4(dfp, S)))
    if(method == 5)
        return(as.dist(1-fun5(dfp, S)))
    }
}
